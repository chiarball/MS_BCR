#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unified pipeline for Ab-Ag DB:

1) Read agab_with_cdr3.tsv (output of step1) and parse 'metadata'
   to extract target_name, target_pdb, target_uniprot.
2) Build 'antigen' using UniProt / PDB.
3) Compute 'origin_class' (human / viral / bacteria / other) using
   UniProt taxonomy + PDB taxonomy (former step3).
4) Disambiguate antigens that appear with >1 origin_class by appending
   the origin_class to the antigen name.

Outputs (all in aggregating_test):
- agab_cdr3_annotated1.csv
- agab_cdr3_annotated1_with_origin.csv
- agab_cdr3_annotated1_with_origin_adjustantigen.csv
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

# ============================================================
# PATHS / CONFIG
# ============================================================

BASE_DIR = Path("/doctorai/chiarba/AbAg_database")
INPUT_CSV = BASE_DIR / "1_final/agab_with_cdr3.tsv"

OUT_DIR = BASE_DIR / "aggregating_test"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STEP2_OUT = OUT_DIR / "agab_cdr3_annotated1.csv"
STEP3_OUT = OUT_DIR / "agab_cdr3_annotated1_with_origin.csv"
STEP35_OUT = OUT_DIR / "agab_cdr3_annotated1_with_origin_adjustantigen.csv"

# Sample mode for quick tests. Set to None for full DB.
SAMPLE_N = None

TIMEOUT = 15          # HTTP timeout
REQS_PER_SEC = 5      # request rate limit
BACKOFF_FACTOR = 0.5  # retry backoff

# Caches (reuse original ones in 1_final so we do not refetch everything)
CACHE_DIR = BASE_DIR / "1_final"
UNIPROT_CACHE_PATH = CACHE_DIR / "uniprot_cache.json"
PDBMAP_CACHE_PATH = CACHE_DIR / "pdb2uniprot_cache.json"
PDB_ENTITY_CACHE_PATH = CACHE_DIR / "pdb_antigen_entity_cache.json"

# Optional family mapping (kept for compatibility)
FAMILY_MAPPING_CSV = CACHE_DIR / "family_mapping.csv"

# ============================================================
# GLOBALS: session + caches
# ============================================================

_session = None
_session_lock = threading.Lock()
_LAST_CALL_TS = 0.0

def get_session() -> requests.Session:
    """Create a shared requests.Session with retry/backoff."""
    global _session
    with _session_lock:
        if _session is None:
            s = requests.Session()
            retries = Retry(
                total=5,
                backoff_factor=BACKOFF_FACTOR,
                status_forcelist=[429, 500, 502, 503, 504],
                allowed_methods=["GET"],
                raise_on_status=False,
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

UNIPROT_CACHE: Dict[str, dict] = load_json_cache(UNIPROT_CACHE_PATH)
PDBMAP_CACHE: Dict[str, Optional[str]] = load_json_cache(PDBMAP_CACHE_PATH)
PDB_ENTITY_CACHE: Dict[str, List[dict]] = load_json_cache(PDB_ENTITY_CACHE_PATH)

# Optional family normalization map (not really used downstream)
FAMILY_MAP: Optional[Dict[str, str]] = None
if FAMILY_MAPPING_CSV.exists():
    try:
        _dfmap = pd.read_csv(FAMILY_MAPPING_CSV).dropna()
        if {"raw", "value"}.issubset(set(_dfmap.columns)):
            FAMILY_MAP = dict(zip(_dfmap["raw"].astype(str), _dfmap["value"].astype(str)))
    except Exception:
        FAMILY_MAP = None

# ============================================================
# BASIC UTILS
# ============================================================

def norm_field(x: Optional[str]) -> str:
    if x is None:
        return ""
    if not isinstance(x, str):
        x = str(x)
    return x.strip()

def is_empty_or_origin(x: str) -> bool:
    x = norm_field(x)
    return (x == "") or (x.lower() == "origin")

def looks_like_uniprot_accession(x: str) -> bool:
    """Loose, case-insensitive UniProt accession check."""
    x = norm_field(x)
    if not x:
        return False
    xu = x.upper()
    return bool(re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[A-Z0-9]{6,10}", xu))

def looks_like_pdb_id(x: str) -> bool:
    """PDB IDs are usually 4 alphanumeric characters."""
    x = norm_field(x)
    return bool(re.fullmatch(r"[0-9A-Za-z]{4}", x))

def normalize_family(name: Optional[str]) -> Optional[str]:
    """Apply optional mapping to standardize family names."""
    if not name:
        return None
    if FAMILY_MAP and name in FAMILY_MAP:
        return FAMILY_MAP[name]
    return name

# ============================================================
# FLEXIBLE CSV READER (from step2)
# ============================================================

def read_flexible_csv(path: Path,
                      prefer_python_engine: bool = True,
                      report: bool = True) -> pd.DataFrame:
    """Robust reader for slightly messy CSV/TSV files."""
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
                    print(f"[READ] Trying sep='{sep}' engine='{eng}'...")
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
                    print(f"[READ] Success with sep='{sep}' engine='{eng}'")
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

# ============================================================
# METADATA PARSING (from step2)
# ============================================================

def safe_literal_eval(s) -> dict:
    """Safely parse metadata strings that look like Python dicts."""
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

# ============================================================
# HEURISTIC TAXONOMY (kept but only used for antigen fallback)
# ============================================================

def heuristic_tax_family_from_text(s: str):
    """Infer taxonomy/family from free-text antigen/target_name."""
    if not s:
        return (None, None)
    txt = str(s).lower().strip().replace("_", " ").replace("-", " ")

    RULES = [
        (r"\b(h\.?\s*sapiens|homo sapiens)\b", "Homo sapiens", None),
        (r"\b(mouse|m\.?\s*musculus|mus musculus)\b", "Mus musculus", None),
        (r"\b(sars\s*cov\s*2|sars2|2019[- ]?ncov|ncov)\b", "Severe acute respiratory syndrome coronavirus 2", "Betacoronavirus"),
        # ... (keep other rules if you want; truncated for brevity or copy full block)
    ]

    for pat, tax, fam in RULES:
        try:
            if re.search(pat, txt, flags=re.IGNORECASE):
                return (tax, fam)
        except re.error:
            continue

    if " mouse" in txt:
        return ("Mus musculus", None)
    if " sars" in txt:
        return ("Severe acute respiratory syndrome coronavirus 2", "Betacoronavirus")
    if " hiv" in txt:
        return ("Human immunodeficiency virus 1", "Retroviridae")
    if "toxin" in txt:
        return (None, "Bacterial toxin (unspecified)")

    return (None, None)

# ============================================================
# UNIPROT / PDB (from step2, slightly extended with taxid)
# ============================================================

def uniprot_fetch(acc: str) -> dict:
    """
    Robust UniProt fetch (entry + search endpoints).
    Returns dict with:
      protein_description, organism, lineage, keywords, pfam[], interpro[], taxid.
    Uses persistent cache. Case-insensitive accession.
    """
    acc = acc.strip().upper()
    cached = UNIPROT_CACHE.get(acc)
    if isinstance(cached, dict) and cached:
        return cached

    sess = get_session()

    def _parse_uniprot_entry_json(rec: dict) -> dict:
        pd = rec.get("proteinDescription", {}) or {}
        rn = pd.get("recommendedName", {}) or {}
        protein_desc = (
            (rn.get("fullName") or {}).get("value")
            or ((pd.get("submissionNames") or [{}])[0].get("fullName") or {}).get("value")
            or ((pd.get("alternativeNames") or [{}])[0].get("fullName") or {}).get("value")
        )

        org = rec.get("organism") or {}
        organism = org.get("scientificName")
        lineage = org.get("lineage") or []
        taxid_raw = org.get("taxonId")
        try:
            taxid = int(taxid_raw) if taxid_raw is not None else None
        except Exception:
            taxid = None

        keywords = [k.get("value") for k in (rec.get("keywords") or []) if isinstance(k, dict) and "value" in k]

        pfam, interpro = [], []
        for xr in rec.get("uniProtKBCrossReferences") or []:
            db = xr.get("database")
            if db == "Pfam":
                pfam.append({
                    "id": xr.get("id"),
                    "name": ((xr.get("properties") or [{}])[0].get("value") if xr.get("properties") else None),
                })
            elif db == "InterPro":
                interpro.append({
                    "id": xr.get("id"),
                    "name": ((xr.get("properties") or [{}])[0].get("value") if xr.get("properties") else None),
                })

        return {
            "protein_description": protein_desc,
            "organism": organism,
            "lineage": lineage,
            "keywords": keywords,
            "pfam": pfam,
            "interpro": interpro,
            "taxid": taxid,
        }

    # 1) Entry endpoint
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

    # 2) Search endpoint as fallback
    try:
        throttle()
        url = (
            "https://rest.uniprot.org/uniprotkb/search"
            f"?query=accession:{acc}"
            "&fields=accession,protein_description,organism_name,organism_id,organism_lineage,"
            "keywords,database(Pfam),database(InterPro)"
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

def uniprot_fetch_min(acc: str) -> dict:
    """
    Minimal UniProt info for origin_class:
    organism, taxid, lineage.
    Uses the same cache as uniprot_fetch().
    """
    acc = acc.strip().upper()
    cached = UNIPROT_CACHE.get(acc)
    if isinstance(cached, dict) and cached:
        return {
            "organism": cached.get("organism"),
            "taxid": cached.get("taxid"),
            "lineage": cached.get("lineage") or [],
        }

    full = uniprot_fetch(acc)
    if not full:
        return {}
    return {
        "organism": full.get("organism"),
        "taxid": full.get("taxid"),
        "lineage": full.get("lineage") or [],
    }

def pdb_to_uniprot(pdb_id: str) -> Optional[str]:
    """Map PDB ID to a UniProt accession (PDBe + RCSB fallback)."""
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
    except Exception:
        pass

    # RCSB fallback (polymer entity 1)
    try:
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

def pdb_fetch_taxonomy(pdb_id: str) -> Optional[str]:
    """Fetch source organism scientific name from PDB via RCSB."""
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
    except Exception:
        pass
    return None

# ============================================================
# ANTIGEN RESOLUTION (step2)
# ============================================================

class Resolved:
    def __init__(self, antigen: Optional[str], taxonomy: Optional[str], family: Optional[str]):
        self.antigen = antigen
        self.taxonomy = taxonomy
        self.family = family

def resolve_antigen_tax_family(target_name: str, target_pdb: str, target_uniprot: str) -> Resolved:
    """
    1) If UniProt accession available -> use UniProt.
    2) Else if PDB ID available -> map to UniProt, then UniProt or PDB taxonomy.
    3) Else if only target_name available -> antigen=target_name.
    """
    tn = norm_field(target_name)
    tp = norm_field(target_pdb)
    tu = norm_field(target_uniprot)

    if tu and not is_empty_or_origin(tu) and looks_like_uniprot_accession(tu):
        u = uniprot_fetch(tu.upper())
        if u:
            antigen = u.get("protein_description") or (tn if tn else tu)
            return Resolved(antigen or None, u.get("organism"), None)

    if tp and not is_empty_or_origin(tp) and looks_like_pdb_id(tp):
        acc = pdb_to_uniprot(tp)
        if acc:
            u = uniprot_fetch(acc.upper())
            if u:
                antigen = u.get("protein_description") or (tn if tn else acc)
                return Resolved(antigen or None, u.get("organism"), None)

        tax_from_pdb = pdb_fetch_taxonomy(tp)
        if tax_from_pdb:
            antigen = tn if tn else tp
            return Resolved(antigen, tax_from_pdb, None)

    if tn and not is_empty_or_origin(tn):
        tax_h, _ = heuristic_tax_family_from_text(tn)
        return Resolved(tn, tax_h, None)

    return Resolved(None, None, None)

# ============================================================
# ORIGIN_CLASS: PDB ENTITY → UNIPROT + TAXONOMY (step3, no antigen)
# ============================================================

def score_description(desc: str, name_hint: str) -> int:
    """Crude similarity score between PDB entity description and a name hint."""
    desc = (desc or "").lower()
    name_hint = (name_hint or "").lower()
    if not desc or not name_hint:
        return 0

    score = 0
    if name_hint in desc or desc in name_hint:
        score += 5

    tokens = [t for t in re.findall(r"[a-z0-9]+", name_hint) if len(t) >= 3]
    for t in tokens:
        if t in desc:
            score += 1

    return score

def pdb_antigen_uniprots_from_entities(pdb_id: str, name_hint: str, max_entities: int = 10) -> list:
    """
    Use RCSB polymer_entity to find the entity most similar to name_hint.
    Return a list of items:
      - {"kind": "uniprot", "acc": "P12345"}
      - {"kind": "taxonomy", "taxid": 9606, "organism": "Homo sapiens"}
    """
    key = f"{pdb_id.lower()}|{name_hint.lower()}"
    if key in PDB_ENTITY_CACHE:
        return PDB_ENTITY_CACHE[key]

    sess = get_session()
    best_score = -1
    best_items: list = []

    for ent in range(1, max_entities + 1):
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{ent}"
        try:
            throttle()
            r = sess.get(url, timeout=TIMEOUT)
        except Exception:
            break

        if r.status_code == 404:
            break
        if r.status_code != 200:
            continue

        try:
            data = r.json()
        except Exception:
            continue

        pe = data.get("rcsb_polymer_entity") or {}
        ep = data.get("entity_poly") or {}
        desc = pe.get("pdbx_description") or ep.get("pdbx_description") or ""

        s = score_description(desc, name_hint)
        if s <= best_score:
            continue

        items: list = []
        cont = data.get("rcsb_polymer_entity_container_identifiers") or {}
        refs = cont.get("reference_sequence_identifiers") or []
        for ref in refs:
            if (ref.get("database_name") or "").lower() == "uniprot":
                acc = ref.get("database_accession")
                if acc:
                    items.append({"kind": "uniprot", "acc": acc.upper()})

        if not items:
            tax_ids = pe.get("source_organism_ids") or []
            sci_names = pe.get("source_scientific_name") or []
            if tax_ids:
                items.append(
                    {
                        "kind": "taxonomy",
                        "taxid": int(tax_ids[0]),
                        "organism": sci_names[0] if sci_names else "",
                    }
                )

        if not items:
            continue

        best_score = s
        best_items = items

    PDB_ENTITY_CACHE[key] = best_items
    save_json_cache(PDB_ENTITY_CACHE_PATH, PDB_ENTITY_CACHE)
    return best_items

def classify_by_text(name: str) -> str:
    """Text-based heuristic when no taxonomy is available."""
    txt = (name or "").lower()

    # 1) Viral: keep this first so viruses always win
    if any(k in txt for k in [
        "virus", "viral", "sars-cov2", "sars cov2", "covid-19",
        "orf", "polyprotein", "gag-pol", "hemagglutinin", "nucleoprotein",
        "envelope glycoprotein",
    ]):
        return "viral"

    # 2) Bacterial
    if any(k in txt for k in [
        "escherichia coli", "e. coli", "staphylococcus", "streptococcus",
        "mycobacterium", "pseudomonas", "salmonella", "bacillus",
    ]):
        return "bacteria"

    # 3) Non-human mammals (mouse, rat, etc.) → other
    if any(k in txt for k in [
        "_mouse", " mouse", "mus musculus",
        " rat ", "rattus norvegicus",
    ]):
        return "other"

    # 4) Clear human taxonomy
    if "homo sapiens" in txt:
        return "human"

    # 5) Generic 'human' in the name, but only if it does not look viral/bacterial
    if "human" in txt and not any(
        k in txt for k in [
            "virus", "viral", "hiv", "influenza", "cytomegalovirus",
            "hepatitis", "herpes", "retrovirus",
            "bacteria", "bacterial", "bacterium",
        ]
    ):
        # Examples: "human PD-1", "igg4_human", "human igg1"
        return "human"

    # 6) Fallback
    return "other"


def classify_origin_row(target_uniprot: str,
                        target_pdb: str,
                        target_name: str) -> str:
    """
    Origin classification WITHOUT using antigen:

    1) Try UniProt taxonomy (target_uniprot).
    2) If UniProt not available or uninformative, try PDB entities (PDB → UniProt or taxonomy).
    3) If taxonomy is not available, fall back to text-based heuristic on target_name.
    4) Final fallback: 'other'.
    """

    tn = norm_field(target_name)
    tu = norm_field(target_uniprot)
    tp = norm_field(target_pdb)

    name_hint = tn

    # ----- 1) UniProt taxonomy -----
    if tu and looks_like_uniprot_accession(tu):
        info = uniprot_fetch_min(tu)
        if info:
            cls_uni = classify_by_taxonomy(
                info.get("taxid"),
                info.get("organism"),
                info.get("lineage") or [],
                name_hint=name_hint,
            )
            if cls_uni in {"human", "viral", "bacteria"}:
                return cls_uni

    # ----- 2) PDB → UniProt / taxonomy -----
    if tp and looks_like_pdb_id(tp):
        items = pdb_antigen_uniprots_from_entities(tp, name_hint)
        if items:
            priority = {"human": 0, "viral": 1, "bacteria": 2, "other": 3}
            best_class = None
            best_rank = 999

            for raw in items:
                if isinstance(raw, dict):
                    item = raw
                else:
                    item = {"kind": "uniprot", "acc": str(raw)}

                kind = item.get("kind")
                if kind == "uniprot":
                    info = uniprot_fetch_min(item["acc"])
                    if not info:
                        continue
                    cls = classify_by_taxonomy(
                        info.get("taxid"),
                        info.get("organism"),
                        info.get("lineage") or [],
                        name_hint=name_hint,
                    )
                elif kind == "taxonomy":
                    cls = classify_by_taxonomy(
                        item.get("taxid"),
                        item.get("organism"),
                        [],
                        name_hint=name_hint,
                    )
                else:
                    continue

                rank = priority.get(cls, 3)
                if rank < best_rank:
                    best_rank = rank
                    best_class = cls
                if best_rank == 0:
                    break

            if best_class is not None and best_class in {"human", "viral", "bacteria"}:
                return best_class

    # ----- 3) Text-only fallback on target_name -----
    if name_hint:
        return classify_by_text(name_hint)

    # ----- 4) Final fallback -----
    return "other"


def classify_by_taxonomy(
    taxid: Optional[int],
    organism: Optional[str],
    lineage_list: Optional[list],
    name_hint: str = "",
) -> str:
    """
    Classify as human / viral / bacteria / other using UniProt taxonomy.

    Priority:
      1) Viral (based on lineage and organism name)
      2) Bacteria
      3) Human (taxid 9606 or 'Homo sapiens')
      4) Fallback: use text only if taxonomy is missing
    """
    org = (organism or "").lower()
    lineage = "; ".join(lineage_list or []).lower()
    name_hint_l = (name_hint or "").lower()

    # ----- 1) Viral -----
    # UniProt viral entries have lineage containing 'Viruses' / 'Viridae'
    # or organism names ending in 'virus' / 'phage'.
    if (
        "viruses" in lineage
        or "viridae" in lineage
        or "virus" in org
        or " phage" in org
    ):
        return "viral"

    # Extra safety: if taxonomy is empty but name_hint is clearly viral
    if not org and not lineage and any(
        k in name_hint_l for k in ["sars-cov2", "sars-cov-2", "covid-19", "orf", "polyprotein", " virus"]
    ):
        return "viral"

    # ----- 2) Bacteria -----
    if "bacteria" in lineage or "bacterium" in org or "bacterial" in org:
        return "bacteria"

    # ----- 3) Human -----
    # Here we are strict: NO plain 'human' to avoid
    # 'Human cytomegalovirus', 'Human immunodeficiency virus', etc.
    if taxid == 9606 or "homo sapiens" in org:
        return "human"

    # ----- 4) Fallback: if taxonomy missing, rely on text heuristic -----
    if not org and not lineage and name_hint_l:
        return classify_by_text(name_hint_l)

    return "other"


# ============================================================
# SMALL HELPERS
# ============================================================

def human_time(seconds: float) -> str:
    seconds = max(0, int(seconds))
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60
    if h > 0:
        return f"{h:d}:{m:02d}:{s:02d}"
    return f"{m:d}:{s:02d}"

# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    start_ts = time.time()

    if not INPUT_CSV.exists():
        print(f"ERROR: Input file not found: {INPUT_CSV}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Reading input: {INPUT_CSV}")
    df = read_flexible_csv(INPUT_CSV)

    if "metadata" not in df.columns:
        print("ERROR: 'metadata' column not found.", file=sys.stderr)
        sys.exit(1)

    if SAMPLE_N is not None and len(df) > SAMPLE_N:
        df = df.sample(SAMPLE_N, random_state=42).copy()
        print(f"[INFO] Sampling {SAMPLE_N} rows.")

    total_rows = len(df)
    print(f"[INFO] Loaded {total_rows} rows.")

    # ---- parse metadata to targets ----
    target_name_list, target_pdb_list, target_uniprot_list = [], [], []
    for val in df["metadata"].tolist():
        md = safe_literal_eval(val)
        target_name_list.append(norm_field(md.get("target_name", "")))
        target_pdb_list.append(norm_field(md.get("target_pdb", "")))
        target_uniprot_list.append(norm_field(md.get("target_uniprot", "")))

    df["target_name"] = target_name_list
    df["target_pdb"] = target_pdb_list
    df["target_uniprot"] = target_uniprot_list

    # ---- build antigen (former step2) ----
    print("[INFO] Resolving antigen...")
    antigens = []
    for tn, tp, tu in zip(df["target_name"], df["target_pdb"], df["target_uniprot"]):
        if all(is_empty_or_origin(x) for x in (tn, tp, tu)):
            antigens.append(pd.NA)
        else:
            res = resolve_antigen_tax_family(tn, tp, tu)
            antigens.append(res.antigen if res.antigen is not None else pd.NA)
    df["antigen"] = antigens

    # save step2-like output
    df.to_csv(STEP2_OUT, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"[INFO] Step2 output written to: {STEP2_OUT}")

    # ---- origin_class (former step3, no antigen) ----
    print("[INFO] Computing origin_class...")
    origins = []
    for i, row in df.iterrows():
        origin = classify_origin_row(
            target_uniprot=row.get("target_uniprot", ""),
            target_pdb=row.get("target_pdb", ""),
            target_name=row.get("target_name", ""),
        )
        origins.append(origin)
        if (i + 1) % 1000 == 0:
            sys.stdout.write(f"\r[origin] {i+1}/{total_rows} rows processed")
            sys.stdout.flush()
    if total_rows >= 1000:
        print()

    df["origin_class"] = origins

    df.to_csv(STEP3_OUT, index=False)
    print(f"[INFO] Step3 output written to: {STEP3_OUT}")

    # ---- unresolved antigen report (optional, like step2) ----
    name_empty = df["target_name"].apply(is_empty_or_origin)
    pdb_empty = df["target_pdb"].apply(is_empty_or_origin)
    uniprot_empty = df["target_uniprot"].apply(is_empty_or_origin)
    all_empty = name_empty & pdb_empty & uniprot_empty
    any_present = ~all_empty
    antigen_na_with_targets = df["antigen"].isna() & any_present

    print(f"[REPORT] empty targets — name:{int(name_empty.sum())} "
          f"pdb:{int(pdb_empty.sum())} uniprot:{int(uniprot_empty.sum())} "
          f"| all_three_empty:{int(all_empty.sum())}")
    print(f"[REPORT] with ≥1 target present — antigen_NA:{int(antigen_na_with_targets.sum())}")

    export_cols = [c for c in [
        "dataset", "heavy_sequence", "light_sequence", "scfv",
        "affinity_type", "affinity", "processed_measurement",
        "cdr3_aa", "cdr3_aa_len", "cdr3_filtered",
        "metadata",
        "target_name", "target_pdb", "target_uniprot",
        "antigen", "origin_class",
    ] if c in df.columns]

    unresolved_df = df.loc[antigen_na_with_targets, export_cols].copy()
    unresolved_path = OUT_DIR / "agab_cdr3_annotated1_targets_present_antigen_NA.csv"
    unresolved_df.to_csv(unresolved_path, index=False)
    print(f"[EXPORT] Unresolved (targets present, antigen NA): {unresolved_path}")

    # ---- step3.5: disambiguate antigen names by origin_class ----
    print("[INFO] Disambiguating antigens by origin_class (step3.5)...")

    if "antigen" not in df.columns or "origin_class" not in df.columns:
        raise ValueError("Columns 'antigen' and 'origin_class' must be present.")

    antigen_series = df["antigen"].astype("string")
    origin_series = df["origin_class"].astype("string")

    n_origins_per_antigen = (
        df
        .dropna(subset=["antigen", "origin_class"])
        .groupby("antigen")["origin_class"]
        .nunique()
    )
    ambiguous_antigens = set(n_origins_per_antigen[n_origins_per_antigen > 1].index)
    print(f"[INFO] Antigens with >1 origin_class: {len(ambiguous_antigens)}")

    def disambiguate(row):
        a = row["antigen"]
        o = row["origin_class"]
        if pd.isna(a):
            return a
        if a not in ambiguous_antigens:
            return a
        if pd.isna(o):
            return a
        return f"{a} ({o})"

    df["antigen_split"] = df.apply(disambiguate, axis=1)
    df["antigen"] = df["antigen_split"]
    df = df.drop(columns=["antigen_split"])

    df.to_csv(STEP35_OUT, index=False)
    print(f"[INFO] Step3.5 final output written to: {STEP35_OUT}")

    elapsed = time.time() - start_ts
    print(f"[DONE] Total elapsed: {human_time(elapsed)} | rows: {total_rows}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[INTERRUPTED]", file=sys.stderr)
        sys.exit(130)

