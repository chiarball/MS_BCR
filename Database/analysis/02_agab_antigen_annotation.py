#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Annotate Ab-Ag DB from 'metadata' column:
- Parse metadata (Python-literal-like dict in a CSV/TSV string)
- Extract target_name, target_pdb, target_uniprot
- Build *only* antigen:
    * If all three targets are empty or 'origin' -> antigen = NA
    * Else prefer UniProt (from target_uniprot or PDB→UniProt), using:
        - protein_description as antigen
        - fallback to target_name / accession / pdb id if needed
- Save CSV incrementally with checkpoints per chunk
- Includes persistent cache, polite throttling, progress with elapsed time and ETA

NOTE: taxonomy/family are NO LONGER stored or exported.
"""

import ast
import csv
import json
import math
import re
import sys
import time
import threading
from pathlib import Path
from typing import Dict, Optional, Tuple, List

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# =========================
# CONFIG (adjust as needed)
# =========================

# New clean input (output of Step 1)
INPUT_CSV = Path("/doctorai/chiarba/AbAg_database/clean/agab_cdr3_extracted1.tsv")

# New clean output, to avoid overwriting the old one
OUTPUT_CSV = INPUT_CSV.with_name("agab_cdr3_annotated1.csv")

# Set to an integer (e.g., 300) to run a small sample first; set to None for full DB.
SAMPLE_N = None

# Chunk size for incremental processing/saving (tune as you like)
CHUNK_SIZE: int = 20_000

# HTTP settings
TIMEOUT = 15                # seconds per request
REQS_PER_SEC = 5            # polite cap ~5 req/s
BACKOFF_FACTOR = 0.5        # for exponential backoff on retries

# Persistent caches (on disk) to avoid re-querying the same accessions/IDs
CACHE_DIR = INPUT_CSV.parent
UNIPROT_CACHE_PATH = CACHE_DIR / "uniprot_cache.json"
PDBMAP_CACHE_PATH  = CACHE_DIR / "pdb2uniprot_cache.json"

# Optional local mapping to normalize family names (kept but unused)
FAMILY_MAPPING_CSV = CACHE_DIR / "family_mapping.csv"

# =========================
# Globals (sessions, cache)
# =========================

_session = None
_session_lock = threading.Lock()
_LAST_CALL_TS = 0.0

def get_session() -> requests.Session:
    """Create a single shared requests.Session with retries/backoff."""
    global _session
    with _session_lock:
        if _session is None:
            s = requests.Session()
            retries = Retry(
                total=5,
                backoff_factor=BACKOFF_FACTOR,
                status_forcelist=[429, 500, 502, 503, 504],
                allowed_methods=["GET"],
                raise_on_status=False
            )
            adapter = HTTPAdapter(max_retries=retries, pool_connections=20, pool_maxsize=50)
            s.mount("http://", adapter)
            s.mount("https://", adapter)
            _session = s
        return _session

def throttle():
    """Ensure we do not exceed REQS_PER_SEC."""
    global _LAST_CALL_TS
    min_interval = 1.0 / max(1, REQS_PER_SEC)
    now = time.time()
    sleep_s = min_interval - (now - _LAST_CALL_TS)
    if sleep_s > 0:
        time.sleep(sleep_s)
    _LAST_CALL_TS = time.time()

def load_json_cache(path: Path) -> dict:
    if path.exists():
        try:
            with open(path, "r") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}

def save_json_cache(path: Path, obj: dict):
    try:
        with open(path, "w") as f:
            json.dump(obj, f)
    except Exception:
        pass

UNIPROT_CACHE = load_json_cache(UNIPROT_CACHE_PATH)
PDBMAP_CACHE  = load_json_cache(PDBMAP_CACHE_PATH)

# Optional family normalization map (kept but not used downstream)
FAMILY_MAP: Optional[Dict[str, str]] = None
if FAMILY_MAPPING_CSV.exists():
    try:
        _dfmap = pd.read_csv(FAMILY_MAPPING_CSV).dropna()
        if {"raw", "value"}.issubset(set(_dfmap.columns)):
            FAMILY_MAP = dict(zip(_dfmap["raw"].astype(str), _dfmap["value"].astype(str)))
    except Exception:
        FAMILY_MAP = None

# =========================
# Utility functions
# =========================

def safe_literal_eval(s) -> dict:
    """Safely parse a Python-literal-like string dict (single quotes allowed)."""
    if isinstance(s, dict):
        return s
    if not isinstance(s, str):
        return {}
    s = s.strip()
    if not s:
        return {}
    try:
        return ast.literal_eval(s)
    except Exception:
        try:
            js = re.sub(r"'", '"', s)
            return json.loads(js)
        except Exception:
            return {}

def norm_field(x: Optional[str]) -> str:
    """Normalize to stripped string."""
    if x is None:
        return ""
    if not isinstance(x, str):
        x = str(x)
    return x.strip()

def is_empty_or_origin(x: str) -> bool:
    x = norm_field(x)
    return (x == "") or (x.lower() == "origin")

def looks_like_uniprot_accession(x: str) -> bool:
    """Loose check for UniProt accession formats (case-insensitive)."""
    x = norm_field(x)
    if not x:
        return False
    xu = x.upper()  # important: make it case-insensitive
    # Very loose patterns for UniProt accessions
    return bool(re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[A-Z0-9]{6,10}", xu))

def looks_like_pdb_id(x: str) -> bool:
    """PDB IDs are usually 4 chars alphanumeric."""
    x = norm_field(x)
    return bool(re.fullmatch(r"[0-9a-zA-Z]{4}", x))

def normalize_family(name: Optional[str]) -> Optional[str]:
    """Apply optional local mapping to standardize family names."""
    if not name:
        return None
    if FAMILY_MAP and name in FAMILY_MAP:
        return FAMILY_MAP[name]
    return name

# =========================
# Remote lookups
# =========================

def heuristic_tax_family_from_text(s: str):
    """
    Infer taxonomy/family from free-text antigen/target_name.
    Kept for compatibility, but taxonomy/family are not exported anymore.
    """
    if not s:
        return (None, None)
    txt = str(s).lower().strip().replace("_", " ").replace("-", " ")

    RULES = [
        # Mammals
        (r"\b(human|h\.?\s*sapiens|homo sapiens)\b", "Homo sapiens", None),
        (r"\b(mouse|m\.?\s*musculus|mus musculus)\b", "Mus musculus", None),
        # Viruses
        (r"\b(sars\s*cov\s*2|sars2|2019[- ]?ncov|ncov)\b", "Severe acute respiratory syndrome coronavirus 2", "Betacoronavirus"),
        (r"\b(hiv[- ]?1|hiv1|hiv)\b", "Human immunodeficiency virus 1", "Retroviridae"),
        (r"\b(influenza|ha[0-9])\b|\bneuraminidase\b", "Influenza A virus", "Orthomyxoviridae"),
        (r"\baav\s*2\b|\badeno[- ]?associated virus 2\b|\baav\b", "Adeno-associated virus 2", "Parvoviridae"),
        (r"\bhbv\b|\bhepatitis b\b", "Hepatitis B virus", "Hepadnaviridae"),
        (r"\bhcv\b|\bhepatitis c\b", "Hepatitis C virus", "Flaviviridae"),
        (r"\bdengue\b|\bzika\b|\byellow fever\b", "Flavivirus", "Flaviviridae"),
        # Bacteria (organisms)
        (r"\bstaphylococcus aureus\b|\bs\.?\s*aureus\b", "Staphylococcus aureus", "Staphylococcaceae"),
        (r"\bescherichia coli\b|\be\.?\s*coli\b", "Escherichia coli", "Enterobacteriaceae"),
        (r"\bshigella\b", "Shigella", "Enterobacteriaceae"),
        (r"\bvibrio cholerae\b|\bv\.?\s*cholerae\b", "Vibrio cholerae", "Vibrionaceae"),
        (r"\bbordetella pertussis\b", "Bordetella pertussis", "Alcaligenaceae"),
        (r"\bcorynebacterium diphtheriae\b", "Corynebacterium diphtheriae", "Corynebacteriaceae"),
        (r"\bclostridium botulinum\b|\bc\.?\s*botulinum\b", "Clostridium botulinum", "Clostridiaceae"),
        (r"\bclostridium tetani\b|\bc\.?\s*tetani\b", "Clostridium tetani", "Clostridiaceae"),
        (r"\bpseudomonas aeruginosa\b|\bp\.?\s*aeruginosa\b", "Pseudomonas aeruginosa", "Pseudomonadaceae"),
        (r"\bbacillus anthracis\b", "Bacillus anthracis", "Bacillaceae"),
        (r"\bstreptococcus pyogenes\b", "Streptococcus pyogenes", "Streptococcaceae"),
        # Bacterial toxins
        (r"\bbotulinum toxin\b|\bbo?t?ox\b", "Clostridium botulinum", "AB toxin (zinc endopeptidase)"),
        (r"\btetanus toxin\b", "Clostridium tetani", "AB toxin (zinc endopeptidase)"),
        (r"\bdiphtheria toxin\b", "Corynebacterium diphtheriae", "AB toxin (ADP-ribosyltransferase)"),
        (r"\bcholera toxin\b", "Vibrio cholerae", "AB5 toxin (ADP-ribosyltransferase)"),
        (r"\bshiga(?:-like)? toxin\b|\bstx\b", "Shigella dysenteriae (or E. coli Stx)", "AB5 toxin (N-glycosidase)"),
        (r"\bpertussis toxin\b", "Bordetella pertussis", "AB5 toxin (ADP-ribosyltransferase)"),
        (r"\bexotoxin a\b", "Pseudomonas aeruginosa", "ADP-ribosyltransferase toxin"),
        (r"\balpha[- ]?toxin\b|\bhaemo?lysin\b", "Staphylococcus aureus (or others)", "Pore-forming toxin"),
        # Plant/ricin-like
        (r"\bricin\b", "Ricinus communis", "Ribosome-inactivating protein"),
        (r"\bsaporin\b", "Saponaria officinalis", "Ribosome-inactivating protein"),
        # Parasites
        (r"\bplasmodium falciparum\b|\bp\.?\s*falciparum\b", "Plasmodium falciparum", "Plasmodiidae"),
        (r"\btrypanosoma brucei\b|\bt\.?\s*brucei\b", "Trypanosoma brucei", "Trypanosomatidae"),
        (r"\bleishmania\b", "Leishmania", "Trypanosomatidae"),
        (r"\btoxoplasma gondii\b", "Toxoplasma gondii", "Sarcocystidae"),
        (r"\bschistosoma\b", "Schistosoma", "Schistosomatidae"),
        # Fungi
        (r"\bcandida albicans\b", "Candida albicans", "Saccharomycetaceae"),
        (r"\baspergillus fumigatus\b", "Aspergillus fumigatus", "Trichocomaceae"),
        (r"\bcryptococcus neoformans\b", "Cryptococcus neoformans", "Tremellaceae"),
        # Human proteins (examples)
        (r"\b(pd-?1|pdcd1)\b", "Homo sapiens", "Immunoglobulin superfamily"),
        (r"\bher2\b|\berbb2\b", "Homo sapiens", "Receptor tyrosine kinases"),
        (r"\bleptin\b|\blep\b", "Homo sapiens", "Cytokine hormones"),
        (r"\bmyc\b", "Homo sapiens", "bHLH transcription factors"),
        (r"\bvegf[a-z]?\b", "Homo sapiens", "Cystine-knot growth factors"),
        (r"\btgf-?b(?:eta)?\b", "Homo sapiens", "TGF-beta family"),
    ]

    for pat, tax, fam in RULES:
        try:
            if re.search(pat, txt, flags=re.IGNORECASE):
                return (tax, fam)
        except re.error:
            continue

    if " human" in txt:
        return ("Homo sapiens", None)
    if " mouse" in txt:
        return ("Mus musculus", None)
    if " sars" in txt:
        return ("Severe acute respiratory syndrome coronavirus 2", "Betacoronavirus")
    if " hiv" in txt:
        return ("Human immunodeficiency virus 1", "Retroviridae")
    if "toxin" in txt:
        return (None, "Bacterial toxin (unspecified)")

    return (None, None)

def uniprot_fetch(accession: str) -> dict:
    """
    Robust UniProt fetch:
      1) Try entry endpoint /uniprotkb/{ACC} (most reliable).
      2) Fallback to search endpoint query=accession:ACC with fields.
    Returns dict with: protein_description, organism, lineage, keywords, pfam[], interpro[].
    Uses persistent cache. Case-insensitive accession.
    """
    acc = accession.strip().upper()
    if acc in UNIPROT_CACHE:
        return UNIPROT_CACHE[acc]

    sess = get_session()

    def _parse_uniprot_entry_json(rec: dict) -> dict:
        # protein description
        pd = rec.get("proteinDescription", {}) or {}
        rn = pd.get("recommendedName", {}) or {}
        protein_desc = (
            (rn.get("fullName") or {}).get("value")
            or ((pd.get("submissionNames") or [{}])[0].get("fullName") or {}).get("value")
            or ((pd.get("alternativeNames") or [{}])[0].get("fullName") or {}).get("value")
        )
        organism = (rec.get("organism") or {}).get("scientificName")
        lineage = (rec.get("organism") or {}).get("lineage") or []
        keywords = [k.get("value") for k in (rec.get("keywords") or []) if isinstance(k, dict) and "value" in k]

        pfam, interpro = [], []
        for xr in rec.get("uniProtKBCrossReferences") or []:
            db = xr.get("database")
            if db == "Pfam":
                pfam.append({"id": xr.get("id"), "name": ((xr.get("properties") or [{}])[0].get("value") if xr.get("properties") else None)})
            elif db == "InterPro":
                interpro.append({"id": xr.get("id"), "name": ((xr.get("properties") or [{}])[0].get("value") if xr.get("properties") else None)})

        return {
            "protein_description": protein_desc,
            "organism": organism,
            "lineage": lineage,
            "keywords": keywords,
            "pfam": pfam,
            "interpro": interpro,
        }

    # (1) Entry endpoint
    try:
        throttle()
        r = sess.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=TIMEOUT)
        if r.status_code == 200:
            rec = r.json()
            out = _parse_uniprot_entry_json(rec)
            UNIPROT_CACHE[acc] = out
            save_json_cache(UNIPROT_CACHE_PATH, UNIPROT_CACHE)
            return out
    except Exception:
        pass

    # (2) Search endpoint as fallback
    try:
        throttle()
        url = (
            "https://rest.uniprot.org/uniprotkb/search"
            f"?query=accession:{acc}"
            "&fields=accession,protein_description,organism_name,organism_id,organism_lineage,keywords,"
            "database(Pfam),database(InterPro)"
            "&format=json&size=1"
        )
        r = sess.get(url, timeout=TIMEOUT)
        if r.status_code == 200:
            data = r.json()
            if data.get("results"):
                rec = data["results"][0]
                out = _parse_uniprot_entry_json(rec)
                UNIPROT_CACHE[acc] = out
                save_json_cache(UNIPROT_CACHE_PATH, UNIPROT_CACHE)
                return out
    except Exception:
        pass

    UNIPROT_CACHE[acc] = {}
    save_json_cache(UNIPROT_CACHE_PATH, UNIPROT_CACHE)
    return {}

def pdb_to_uniprot(pdb_id: str) -> Optional[str]:
    """
    Map PDB ID to a UniProt accession.
    Order:
      - PDBe/SIFTS
      - RCSB fallback
    Cached on disk.
    """
    key = pdb_id.strip().lower()
    if key in PDBMAP_CACHE:
        return PDBMAP_CACHE[key]

    sess = get_session()

    # PDBe/SIFTS
    try:
        throttle()
        r = sess.get(f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{key}", timeout=TIMEOUT)
        if r.status_code == 200:
            data = r.json()
            mappings = (data.get(key) or {}).get("UniProt") or {}
            for acc in mappings.keys():
                PDBMAP_CACHE[key] = acc.upper()
                save_json_cache(PDBMAP_CACHE_PATH, PDBMAP_CACHE)
                return PDBMAP_CACHE[key]
        elif r.status_code == 404:
            pass
    except Exception:
        pass

    # RCSB fallback
    try:
        throttle()
        r = sess.get(f"https://data.rcsb.org/rest/v1/core/entry/{key}", timeout=TIMEOUT)
        if r.status_code == 200:
            data = r.json()
            # try to dig uniprot IDs (varies by entry); use polymer entities endpoint if needed
            throttle()
            r2 = sess.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{key}/1", timeout=TIMEOUT)
            if r2.status_code == 200:
                e = r2.json()
                dbrefs = (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("reference_sequence_identifiers") or []
                for ref in dbrefs:
                    if ref.get("database_name", "").lower() == "uniprot":
                        PDBMAP_CACHE[key] = ref.get("database_accession", "").upper() or None
                        save_json_cache(PDBMAP_CACHE_PATH, PDBMAP_CACHE)
                        return PDBMAP_CACHE[key]
    except Exception:
        pass

    PDBMAP_CACHE[key] = None
    save_json_cache(PDBMAP_CACHE_PATH, PDBMAP_CACHE)
    return None

def uniprot_to_best_pdbs(acc: str, max_n: int = 3) -> List[str]:
    """
    Return up to max_n PDB IDs mapped to a UniProt accession using PDBe 'best_structures'.
    """
    acc = acc.strip().upper()
    sess = get_session()
    try:
        throttle()
        r = sess.get(f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{acc}", timeout=TIMEOUT)
        if r.status_code != 200:
            return []
        data = r.json().get(acc, [])
        pdbs = [d.get("structure_id") for d in data if d.get("structure_id")]
        out = []
        for p in pdbs:
            if p not in out:
                out.append(p)
        return out[:max_n]
    except Exception:
        return []

def choose_family_from_pfam_interpro(pfam_list, interpro_list) -> Optional[str]:
    """Prefer a clean family/domain name from Pfam or InterPro (kept for compatibility)."""
    if pfam_list:
        for item in pfam_list:
            name = (item or {}).get("name")
            if name and len(name) >= 3:
                return name
        if pfam_list[0].get("id"):
            return f"Pfam:{pfam_list[0]['id']}"
    if interpro_list:
        for item in interpro_list:
            name = (item or {}).get("name")
            if name and len(name) >= 3:
                return name
        if interpro_list[0].get("id"):
            return f"InterPro:{interpro_list[0]['id']}"
    return None

def family_from_keywords(keywords: list) -> Optional[str]:
    """Fallback if Pfam/InterPro did not yield something descriptive."""
    if not keywords:
        return None
    for k in keywords:
        if "family" in (k or "").lower():
            return k
    candidates = ["toxin", "capsid", "envelope protein", "glycoprotein", "nucleoprotein"]
    for k in keywords:
        if any(c in (k or "").lower() for c in candidates):
            return k
    return None

def pdb_fetch_taxonomy(pdb_id: str) -> Optional[str]:
    """
    Attempt to fetch source organism scientific name from PDB via RCSB.
    Taxonomy is not exported anymore, but this is kept for full compatibility.
    """
    key = pdb_id.strip().lower()
    sess = get_session()
    try:
        for ent in range(1, 9):
            throttle()
            r = sess.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{key}/{ent}", timeout=TIMEOUT)
            if r.status_code != 200:
                continue
            e = r.json()
            src = e.get("rcsb_entity_source_organism") or e.get("rcsb_polymer_entity", {}).get("rcsb_entity_source_organism") or []
            if isinstance(src, dict):
                src = [src]
            for s in src:
                sci = s.get("scientific_name") or (s.get("ncbi_taxonomy") or {}).get("scientific_name")
                if sci:
                    return sci
            tax = e.get("entity_src_gen", []) or e.get("entity_src_nat", [])
            if isinstance(tax, dict):
                tax = [tax]
            for s in tax:
                sci = s.get("organism_scientific")
                if sci:
                    return sci
    except Exception:
        pass
    return None

# =========================
# Resolution logic per row
# =========================

class Resolved:
    def __init__(self, antigen: Optional[str], taxonomy: Optional[str], family: Optional[str]):
        self.antigen = antigen
        self.taxonomy = taxonomy
        self.family = family

def resolve_antigen_tax_family(target_name: str, target_pdb: str, target_uniprot: str) -> Resolved:
    """
    Resolution priority:
    1) If all three empty/origin -> NA, NA, NA (handled by caller)
    2) If UniProt accession available -> fetch UniProt
    3) Else if PDB id available -> map to UniProt via PDBe/RCSB, then fetch UniProt
       If mapping/fetch fails, try taxonomy directly from PDB (RCSB)
    4) Else if target_name provided -> antigen=target_name, taxonomy via heuristics if possible

    NOTE: caller now uses ONLY res.antigen. taxonomy/family are ignored.
    """
    tn = norm_field(target_name)
    tp = norm_field(target_pdb)
    tu = norm_field(target_uniprot)

    # UniProt first
    if tu and not is_empty_or_origin(tu) and looks_like_uniprot_accession(tu):
        u = uniprot_fetch(tu.upper())
        if u:
            antigen = u.get("protein_description") or (tn if tn else tu)
            return Resolved(antigen or None, u.get("organism"), None)

    # PDB -> UniProt -> UniProt fetch
    if tp and not is_empty_or_origin(tp) and looks_like_pdb_id(tp):
        acc = pdb_to_uniprot(tp)
        if acc:
            u = uniprot_fetch(acc.upper())
            if u:
                antigen = u.get("protein_description") or (tn if tn else acc)
                return Resolved(antigen or None, u.get("organism"), None)

        # Fallback: taxonomy from PDB directly (kept, but not exported)
        tax_from_pdb = pdb_fetch_taxonomy(tp)
        if tax_from_pdb:
            antigen = tn if tn else tp
            return Resolved(antigen, tax_from_pdb, None)

    # Fallback: target_name only + heuristic taxonomy (taxonomy ignored downstream)
    if tn and not is_empty_or_origin(tn):
        tax_h, _ = heuristic_tax_family_from_text(tn)
        return Resolved(tn, tax_h, None)

    # Nothing resolved
    return Resolved(None, None, None)

# =========================
# Progress helpers
# =========================

def human_time(seconds: float) -> str:
    """Format seconds into H:MM:SS."""
    seconds = max(0, int(seconds))
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60
    if h > 0:
        return f"{h:d}:{m:02d}:{s:02d}"
    return f"{m:d}:{s:02d}"

def print_progress(done_rows: int, total_rows: int, start_ts: float):
    """Print progress line with percentage, elapsed, ETA (rough)."""
    elapsed = time.time() - start_ts
    pct = (done_rows / total_rows) * 100 if total_rows else 100.0
    if done_rows > 0:
        rate = done_rows / max(elapsed, 1e-6)
        remaining = total_rows - done_rows
        eta = remaining / max(rate, 1e-6)
    else:
        eta = float("inf")
    sys.stdout.write(
        f"\rProgress: {done_rows}/{total_rows} ({pct:5.1f}%) | Elapsed: {human_time(elapsed)} | ETA: {('??:??' if math.isinf(eta) else human_time(eta))}"
    )
    sys.stdout.flush()

# --- Robust CSV reader (same as before) ---

import io

def read_flexible_csv(path: Path,
                      prefer_python_engine: bool = True,
                      report: bool = True) -> pd.DataFrame:
    """
    Read a messy CSV/TSV robustly.
    """
    def _clean_columns(cols):
        cleaned = []
        for c in cols:
            s = str(c).replace("\ufeff", "").strip().strip('"').strip("'")
            canon = re.sub(r"[\s\-]+", "_", s).lower()
            if canon in {"metadata", "meta_data"}:
                cleaned.append("metadata")
            else:
                cleaned.append(s)
        return cleaned

    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        header_line = ""
        while True:
            line = fh.readline()
            if not line:
                break
            if line.strip():
                header_line = line.rstrip("\n\r")
                break
        if not header_line:
            raise ValueError("Empty file or no non-empty header line found.")

    sniffed_sep = None
    try:
        dialect = csv.Sniffer().sniff(header_line, delimiters=[",", ";", "\t", "|"])
        sniffed_sep = dialect.delimiter
    except Exception:
        pass
    seps = [sniffed_sep] if sniffed_sep else []
    for s in [",", ";", "\t", "|"]:
        if s not in seps:
            seps.append(s)

    def parse_header_cols(sep: str) -> list:
        reader = csv.reader([header_line], delimiter=sep, quotechar='"')
        cols = next(reader)
        return _clean_columns(cols)

    engines = ["python", "c"] if prefer_python_engine else ["c", "python"]
    last_err = None

    for sep in seps:
        try:
            header_cols = parse_header_cols(sep)
            if report:
                print(f"[READ] Candidate sep='{sep}' -> {len(header_cols)} columns.")
            if len(header_cols) < 2:
                continue
        except Exception as e:
            last_err = e
            continue

        for eng in engines:
            try:
                if report:
                    print(f"[READ] Trying sep='{sep}' engine='{eng}' with explicit header (strict)...")
                df = pd.read_csv(
                    path,
                    sep=sep,
                    engine=eng,
                    low_memory=False,
                    header=None,
                    names=header_cols,
                    skiprows=1,
                    encoding="utf-8-sig",
                )
                if report:
                    print(f"[READ] Success with sep='{sep}' engine='{eng}' (strict)")
                return df
            except Exception as e:
                last_err = e

            try:
                if report:
                    print(f"[READ] Retry sep='{sep}' engine='{eng}' with on_bad_lines='skip'...")
                df = pd.read_csv(
                    path,
                    sep=sep,
                    engine=eng,
                    low_memory=False,
                    header=None,
                    names=header_cols,
                    skiprows=1,
                    encoding="utf-8-sig",
                    on_bad_lines="skip",
                )
                if report:
                    print(f"[READ] Success (skipping bad lines) with sep='{sep}' engine='{eng}'")
                return df
            except Exception as e2:
                last_err = e2

    raise last_err

# =========================
# Main
# =========================

def main():
    start_ts = time.time()

    if not INPUT_CSV.exists():
        print(f"ERROR: Input file not found: {INPUT_CSV}", file=sys.stderr)
        sys.exit(1)

    df = read_flexible_csv(INPUT_CSV)

    if "metadata" not in df.columns:
        print("ERROR: 'metadata' column not found.", file=sys.stderr)
        sys.exit(1)

    if SAMPLE_N is not None and len(df) > SAMPLE_N:
        df = df.sample(SAMPLE_N, random_state=42).copy()
        print(f"[INFO] Sampling {SAMPLE_N} rows for a quick run.")

    total_rows = len(df)
    print(f"[INFO] Loaded {total_rows} rows.")

    # Parse metadata and extract targets
    target_name_list, target_pdb_list, target_uniprot_list = [], [], []
    for val in df["metadata"].tolist():
        md = safe_literal_eval(val)
        target_name_list.append(norm_field(md.get("target_name", "")))
        target_pdb_list.append(norm_field(md.get("target_pdb", "")))
        target_uniprot_list.append(norm_field(md.get("target_uniprot", "")))

    df["target_name"] = target_name_list
    df["target_pdb"] = target_pdb_list
    df["target_uniprot"] = target_uniprot_list

    # Prepare output column (antigen only)
    df["antigen"] = pd.NA

    rows = len(df)
    start = 0
    processed = 0

    # Write header immediately (empty file with columns)
    df.iloc[0:0].to_csv(OUTPUT_CSV, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"[INFO] Writing to: {OUTPUT_CSV}")

    while start < rows:
        end = min(start + CHUNK_SIZE, rows)
        sub = df.iloc[start:end].copy()

        antigens = []
        for tn, tp, tu in zip(sub["target_name"], sub["target_pdb"], sub["target_uniprot"]):
            if all(is_empty_or_origin(x) for x in (tn, tp, tu)):
                antigens.append(pd.NA)
            else:
                res = resolve_antigen_tax_family(tn, tp, tu)
                antigens.append(res.antigen if res.antigen is not None else pd.NA)

            processed += 1
            if processed % 200 == 0 or processed == rows:
                print_progress(processed, rows, start_ts)

        df.loc[sub.index, "antigen"] = antigens

        # Append this chunk to CSV
        df.iloc[:end].to_csv(OUTPUT_CSV, index=False, quoting=csv.QUOTE_MINIMAL)

        start = end

    # Summary report (only antigen now)
    name_empty    = df["target_name"].apply(is_empty_or_origin)
    pdb_empty     = df["target_pdb"].apply(is_empty_or_origin)
    uniprot_empty = df["target_uniprot"].apply(is_empty_or_origin)
    all_empty     = name_empty & pdb_empty & uniprot_empty
    any_present   = ~all_empty

    antigen_na_with_targets = df["antigen"].isna() & any_present

    print(f"[REPORT] empty targets — name:{int(name_empty.sum())} "
          f"pdb:{int(pdb_empty.sum())} uniprot:{int(uniprot_empty.sum())} "
          f"| all_three_empty:{int(all_empty.sum())}")

    print(f"[REPORT] with ≥1 target present — antigen_NA:{int(antigen_na_with_targets.sum())}")

    # Unresolved rows with at least one target present (antigen NA)
    mask_any_present = ~(df["target_name"].apply(is_empty_or_origin) &
                         df["target_pdb"].apply(is_empty_or_origin) &
                         df["target_uniprot"].apply(is_empty_or_origin))

    mask_antigen_na = df["antigen"].isna() & mask_any_present
    unresolved = df[mask_antigen_na][
        ["target_name", "target_pdb", "target_uniprot", "antigen"]
    ]
    print("\n[DEBUG] Unresolved rows (with ≥1 target and antigen NA):")
    print(unresolved.to_string(index=False))

    # Export unresolved-with-targets (antigen NA only)
    export_cols = [c for c in [
        "dataset", "heavy_sequence", "light_sequence", "scfv",
        "affinity_type", "affinity", "processed_measurement",
        "cdr3_aa", "cdr3_aa_len", "cdr3_filtered",
        "metadata",
        "target_name", "target_pdb", "target_uniprot",
        "antigen"
    ] if c in df.columns]

    unresolved_df = df.loc[mask_antigen_na, export_cols].copy()
    unresolved_path = OUTPUT_CSV.with_name(OUTPUT_CSV.stem + "_targets_present_antigen_NA.csv")
    unresolved_df.to_csv(unresolved_path, index=False)
    print(f"[EXPORT] Unresolved with targets → {len(unresolved_df)} rows written to: {unresolved_path}")

    print()
    total_elapsed = time.time() - start_ts
    print(f"[DONE] Wrote: {OUTPUT_CSV}")
    print(f"[STATS] Elapsed: {human_time(total_elapsed)} | Rows: {rows} | Sample mode: {SAMPLE_N is not None}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Partial results saved if any.", file=sys.stderr)
        sys.exit(130)
