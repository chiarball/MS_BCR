#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 3: MS-only CDR3, AGAB–MS matching (lev 0–2), and Fisher test on antigen.

Pipeline:
1) Build MS_only1 from MS_HC_naive_h0h1.tsv:
   - normalize CDR3
   - length >= MIN_LEN
   - exclude all cluster_h1 that appear in subjects starting with 'H' (HC)
   - drop duplicate CDR3
   -> save to /doctorai/chiarba/AbAg_database/clean/MS_only1_cdr3.csv

2) Match AGAB vs MS_only1:
   - AGAB: /doctorai/chiarba/AbAg_database/clean/agab_cdr3_annotated1.csv
   - MS-only1: above file
   - progressive Levenshtein with k = 1..2
   - 4-mer anchors prefilter, min_hits = 2
   - keep matches at the smallest k giving any hit for each AGAB
   -> save matches to /doctorai/chiarba/AbAg_database/clean/agab_vs_ms_lev_matches1.tsv

3) Fisher test (antigen only):
   - use full MS/HC table for metadata and cluster filtering:
       /doctorai/chiarba/analysis/MS_HC_naive_h0h1.tsv
   - filter length >= MIN_LEN and remove H-clusters (same rule as step 1)
   - use matches with lev <= 2
   - compute Fisher enrichment per antigen, separately for lev 0, 1, 2
   - enrich matches with MS and AGAB metadata
   - write summary tables

Code comments are in English.
"""

import os
import sys
import time
import argparse
import tempfile
import shutil
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Set

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact

# Be nice on shared servers
# os.environ.setdefault("OMP_NUM_THREADS", "1")
# os.environ.setdefault("MKL_NUM_THREADS", "1")
# os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
# os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
# try:
#     os.nice(10)
# except Exception:
#     pass

# =========================
# CONFIG
# =========================

AGAB_PATH      = "/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/agab_cdr3_annotated3_disambiguated.csv"
MS_FULL_PATH   = "/doctorai/chiarba/analysis/MS_HC_naive_h0h1.tsv" #OLD
# MS_FULL_PATH   = "/doctorai/chiarba/dataset/FINAL/MS_BCR_DB_Final_Version.tsv"
MS_ONLY_PATH   = "/doctorai/chiarba/AbAg_database/clean/MS_only1_cdr3.csv"
MATCHES_PATH   = "/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/agab_vs_ms_lev_matches1.tsv"
OUT_PREFIX     = "/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/fisher_results_annot2"

MIN_LEN        = 5       # AA length filter
MAX_LEV        = 2       # stop at lev = 2
MIN_ANCHOR_HITS = 2      # minimum shared 4-mers for candidate admission

# Column name expectations
COLS_AGAB = dict(cdr3="cdr3_aa", antigen="antigen", origin="origin_class")
COLS_MS   = dict(cdr3="cdr3_aa", subject="subject_number", tissue="tissue", cluster="cluster_h1")
COLS_MATCH = dict(cdr3_agab="cdr3_agab", cdr3_ms="cdr3_ms", lev="lev")

# =========================
# Levenshtein backend
# =========================

BACKEND = None
def setup_backend():
    """Pick the fastest available edit-distance backend."""
    global BACKEND, edlib, rfuzz, Lev
    edlib = rfuzz = Lev = None
    try:
        import edlib  # type: ignore
        BACKEND = "edlib"
        return
    except Exception:
        pass
    try:
        from rapidfuzz.distance import Levenshtein as rfuzz  # type: ignore
        BACKEND = "rapidfuzz"
        return
    except Exception:
        pass
    try:
        import Levenshtein as Lev  # type: ignore
        BACKEND = "python-Levenshtein"
        return
    except Exception:
        pass
    BACKEND = "builtin"

def _lev_banded(a: str, b: str, k: int) -> int:
    """Banded Levenshtein (Ukkonen). Returns k+1 if distance > k."""
    n, m = len(a), len(b)
    if n == 0:
        return m
    if m == 0:
        return n
    if abs(n - m) > k:
        return k + 1
    if n > m:
        a, b, n, m = b, a, m, n
    band = 2 * k + 1
    prev = [k + 1] * band
    curr = [k + 1] * band

    def bi(i, j):
        return j - i + k

    for j in range(0, min(m, k) + 1):
        prev[bi(0, j)] = j

    for i in range(1, n + 1):
        for t in range(band):
            curr[t] = k + 1
        j_min, j_max = max(0, i - k), min(m, i + k)
        for j in range(j_min, j_max + 1):
            bidx = bi(i, j)
            sub = prev[bi(i - 1, j - 1)] + (a[i - 1] != b[j - 1]) if 0 <= bi(i - 1, j - 1) < band else k + 1
            ins = curr[bi(i, j - 1)] + 1 if 0 <= bi(i, j - 1) < band else k + 1
            dele = prev[bi(i - 1, j)] + 1 if 0 <= bi(i - 1, j) < band else k + 1
            curr[bidx] = min(sub, ins, dele)
        if min(curr) > k:
            return k + 1
        prev, curr = curr, prev
    d = prev[bi(n, m)]
    return d if d <= k else k + 1

def edit_distance(a: str, b: str, k: int) -> int:
    """Return edit distance with early-exit cutoff k when available (k+1 if >k)."""
    if abs(len(a) - len(b)) > k:
        return k + 1
    if BACKEND == "edlib":
        res = edlib.align(a, b, task="distance", mode="NW", k=k)
        d = res["editDistance"]
        return k + 1 if d < 0 else d
    elif BACKEND == "rapidfuzz":
        d = rfuzz.distance(a, b, score_cutoff=k)
        return d if d is not None else k + 1
    elif BACKEND == "python-Levenshtein":
        return Lev.distance(a, b)
    else:
        return _lev_banded(a, b, k)

# =========================
# Normalization & helpers
# =========================

def normalize_cdr3_series(s: pd.Series) -> pd.Series:
    """Uppercase and strip leading C / trailing W."""
    return (
        s.astype("string")
         .str.upper()
         .str.replace(r"^C", "", regex=True)
         .str.replace(r"W$", "", regex=True)
    )

def read_cdr3(path: str, col: str = "cdr3_aa", vgene_col=None) -> pd.DataFrame:
    """Read CSV/TSV; keep & normalize a single CDR3 column; drop duplicates."""

    # Choose separator based on file extension
    if isinstance(path, pd.DataFrame):
        df = path
    elif path.endswith(".csv"):
        df = pd.read_csv(path, sep=",")
    elif path.endswith(".tsv"):
        df = pd.read_csv(path, sep="\t")
    else:
        # Fallback to automatic detection if extension is unknown
        df = pd.read_csv(path, sep=None, engine="python")

    # If expected column is missing, try to recover from typical issues
    if col not in df.columns:
        # Some files may have header truncated as 'cdr3_'
        if "cdr3_" in df.columns:
            # Rename truncated header to the expected name
            df = df.rename(columns={"cdr3_": col})
        else:
            raise ValueError(
                f"Column '{col}' not found in {path}. "
                f"Available columns: {list(df.columns)}"
            )

    out = df[[col]].rename(columns={col: "cdr3_aa"}).copy()
    out["cdr3_aa"] = normalize_cdr3_series(out["cdr3_aa"])
    out = out[out["cdr3_aa"].notna() & (out["cdr3_aa"] != "")]
    out["len"] = out["cdr3_aa"].str.len().astype("int32")
    out = out[out["len"] >= MIN_LEN]
    if vgene_col is not None:
        print(f"colomns are {df.columns}")
        out['v_gene']= df[ vgene_col]
    out = out.drop_duplicates(subset=["cdr3_aa"]).reset_index(drop=True)
    return out


# =========================
# 4-mer anchors and index
# =========================

def tokens_pms4(s: str) -> tuple:
    """Return up to four 4-aa anchors (prefix, Q1-centered, middle, suffix)."""
    n = len(s)
    if n < 4:
        return tuple()
    anchors = []
    anchors.append(s[:4])
    anchors.append(s[max(0, n // 4 - 2): max(0, n // 4 - 2) + 4])
    anchors.append(s[(n - 4) // 2: (n - 4) // 2 + 4])
    anchors.append(s[-4:])
    out, seen = [], set()
    for a in anchors:
        if len(a) == 4 and a not in seen:
            seen.add(a)
            out.append(a)
    return tuple(out)

class MSIndex:
    """Length- and 4-mer anchor-based inverted index over MS with min_hits threshold."""
    def __init__(self, ms_df: pd.DataFrame):
        self.by_len: Dict[int, List[str]] = {int(L): sub["cdr3_aa"].tolist()
                                             for L, sub in ms_df.groupby("len")}
        self.idx: Dict[int, Dict[str, List[int]]] = {}
        for L, seqs in self.by_len.items():
            tok2idx: Dict[str, List[int]] = defaultdict(list)
            for j, s in enumerate(seqs):
                for t in tokens_pms4(s):
                    tok2idx[t].append(j)
            self.idx[L] = tok2idx

    def candidates(self, Ls: List[int], anchors: tuple, min_hits: int = 2) -> Dict[int, List[int]]:
        """Return candidate indices per length with at least 'min_hits' shared anchors."""
        res: Dict[int, List[int]] = {}
        if not anchors:
            return res
        for L in Ls:
            tokmap = self.idx.get(L)
            if not tokmap:
                continue
            counts = defaultdict(int)
            for t in anchors:
                for j in tokmap.get(t, ()):
                    counts[j] += 1
            posts = [j for j, c in counts.items() if c >= min_hits]
            if posts:
                res[L] = sorted(posts)
        return res

# =========================
# File lock (for MS_only1 update, if used)
# =========================

class FileLock:
    """Advisory file lock using a .lock file (POSIX flock if available)."""
    def __init__(self, target_path: str):
        self.lock_path = target_path + ".lock"
        self.fh = None
    def __enter__(self):
        self.fh = open(self.lock_path, "w")
        try:
            import fcntl
            fcntl.flock(self.fh.fileno(), fcntl.LOCK_EX)
        except Exception:
            pass
        return self
    def __exit__(self, exc_type, exc, tb):
        try:
            if self.fh:
                try:
                    import fcntl
                    fcntl.flock(self.fh.fileno(), fcntl.LOCK_UN)
                except Exception:
                    pass
                self.fh.close()
        finally:
            pass

# =========================
# Build MS_only1 from MS/HC table
# =========================

def filter_ms_exclude_H_clusters(ms: pd.DataFrame, return_HC_df: bool = False) -> pd.DataFrame:
    """Remove every row whose cluster_h1 appears in ANY subject_number starting with 'H'."""
    subj_col = COLS_MS["subject"]
    clus_col = COLS_MS["cluster"]
    clusters_with_H = set(
        ms.loc[ms[subj_col].astype(str).str.startswith("H", na=False), clus_col]
          .dropna()
          .unique()
          .tolist()
    )
    if len(clusters_with_H) == 0:
        return ms.copy()
    if return_HC_df:
        return (ms[ms[clus_col].isin(clusters_with_H) | ms[clus_col].isna()].copy())
    else:   
        return (ms[~ms[clus_col].isin(clusters_with_H) | ms[clus_col].isna()].copy())

def build_ms_only1(ms_full_path: str, ms_only_path: str):
    """Create MS_only1_cdr3.csv from full MS/HC table, excluding H-clusters and duplicates."""
    ms = pd.read_csv(ms_full_path, sep="\t")
    for col in [COLS_MS["cdr3"], COLS_MS["subject"], COLS_MS["tissue"], COLS_MS["cluster"]]:
        if col not in ms.columns:
            raise ValueError(f"Column '{col}' not found in {ms_full_path}. Columns: {list(ms.columns)}")

    ms[COLS_MS["cdr3"]] = normalize_cdr3_series(ms[COLS_MS["cdr3"]])
    ms = ms[(ms[COLS_MS["cdr3"]].notna()) & (ms[COLS_MS["cdr3"]] != "")]
    ms["len"] = ms[COLS_MS["cdr3"]].str.len().astype("int32")
    ms = ms[ms["len"] >= MIN_LEN].copy()

    ms_filt = filter_ms_exclude_H_clusters(ms)

    ms_cdr3_only = (
        ms_filt[[COLS_MS["cdr3"]]]
        .drop_duplicates(subset=[COLS_MS["cdr3"]])
        .rename(columns={COLS_MS["cdr3"]: "cdr3_aa"})
        .reset_index(drop=True)
    )
    ms_cdr3_only.to_csv(ms_only_path, index=False)
    print(f"[STEP 3.1] MS_only1 CDR3 written to: {ms_only_path} "
          f"(N unique CDR3: {len(ms_cdr3_only)})")

# =========================
# Core matching logic (progressive k up to MAX_LEV)
# =========================

def match_all(agab_df: pd.DataFrame,
              ms_df: pd.DataFrame,
              max_lev: int,
              throttle_every: int,
              sleep_sec: float,
              min_hits: int = 2):
    """
    Progressive deepening over k: try k=1, then 2, ... up to max_lev.
    For each AGAB sequence:
      * Evaluate all candidates at k=1; if any matches are found, keep ONLY those and stop for this AGAB.
      * Otherwise try k=2; if found, keep ONLY k=2 matches; otherwise continue.
    Returns:
      rows: list of (cdr3_agab, cdr3_ms, lev) with matches retained.
      best_for_ms: dict ms_seq -> (best_lev, agab_seq) computed over retained rows only.
      checked: total number of Levenshtein calls executed.
    """
    ms_index = MSIndex(ms_df)
    rows: List[Tuple[str, str, int]] = []
    best_for_ms: Dict[str, Tuple[int, str]] = {}
    checked = 0
    t0 = time.time()

    agab_list = agab_df["cdr3_aa"].tolist()
    n_agab = len(agab_list)
    hit_by_k = defaultdict(int)

    for i, a in enumerate(agab_list, 1):
        La = len(a)
        anchors = tokens_pms4(a)

        for k in range(1, max_lev + 1):
            Lms = [L for L in range(La - k, La + k + 1) if L in ms_index.by_len]
            if not Lms:
                continue

            cand_map = ms_index.candidates(Lms, anchors, min_hits=min_hits)
            if not cand_map:
                continue

            local_rows: List[Tuple[str, str, int]] = []
            for Lm, idxs in cand_map.items():
                B = ms_index.by_len[Lm]
                for j in idxs:
                    b = B[j]
                    d = edit_distance(a, b, k)
                    checked += 1
                    if d <= k:
                        local_rows.append((a, b, d))

            if local_rows:
                rows.extend(local_rows)
                hit_by_k[k] += 1
                for (_, b, d) in local_rows:
                    if (b not in best_for_ms) or (d < best_for_ms[b][0]) or (d == best_for_ms[b][0] and a < best_for_ms[b][1]):
                        best_for_ms[b] = (d, a)
                break

        if i % 1000 == 0:
            elapsed = time.time() - t0
            rate = i / max(elapsed, 1e-9)
            eta = (n_agab - i) / max(rate, 1e-9)
            k1, k2 = hit_by_k.get(1, 0), hit_by_k.get(2, 0)
            print(f"[INFO] processed {i}/{n_agab} AGAB | "
                  f"matches: {len(rows)} | checked: {checked} | "
                  f"elapsed: {elapsed:.1f}s | ETA: {eta/3600:.1f}h | "
                  f"k1={k1} k2={k2}",
                  flush=True)

    return rows, best_for_ms, checked


def match_all_ighv(agab_df: pd.DataFrame,
              ms_df: pd.DataFrame,
              max_lev: int,
              throttle_every: int,
              sleep_sec: float,
              min_hits: int = 2,
              vgene_col='v_gene'):
    ms_index = MSIndex(ms_df)

    # Build maps cdr3 -> set(v_genes) for quick lookup
    def _build_vmap(df: pd.DataFrame, cdr3_col: str = "cdr3_aa", vcol: str = vgene_col) -> Dict[str, Set[str]]:
        vm: Dict[str, Set[str]] = {}
        if vcol not in df.columns:
            # no vgene column present -> empty sets
            for s in df[cdr3_col].astype(str).tolist():
                vm[s] = set()
            return vm
        for _, r in df[[cdr3_col, vcol]].drop_duplicates().iterrows():
            s = str(r[cdr3_col])
            v = r[vcol]
            if pd.isna(v):
                continue
            # allow multiple genes separated by '|' or ';' if present
            if isinstance(v, str) and (("|" in v) or (";" in v)):
                parts = [p.strip() for delim in ("|", ";") for p in v.split(delim)]
                vals = {p for p in parts if p}
            else:
                vals = {str(v).strip()} if v != "" else set()
            vm.setdefault(s, set()).update(vals)
        # ensure every sequence has an entry
        for s in ms_index.by_len.values():
            for seq in s:
                vm.setdefault(seq, set())
        return vm

    ms_vmap = _build_vmap(ms_df, cdr3_col="cdr3_aa", vcol=vgene_col)
    agab_vmap = _build_vmap(agab_df, cdr3_col="cdr3_aa", vcol=vgene_col)

    rows: List[Tuple[str, str, int, bool]] = []
    best_for_ms: Dict[str, Tuple[int, str]] = {}
    checked = 0
    t0 = time.time()

    agab_list = agab_df["cdr3_aa"].tolist()
    n_agab = len(agab_list)
    hit_by_k = defaultdict(int)

    for i, a in enumerate(agab_list, 1):
        La = len(a)
        anchors = tokens_pms4(a)

        for k in range(1, max_lev + 1):
            Lms = [L for L in range(La - k, La + k + 1) if L in ms_index.by_len]
            if not Lms:
                continue

            cand_map = ms_index.candidates(Lms, anchors, min_hits=min_hits)
            if not cand_map:
                continue

            local_rows: List[Tuple[str, str, int, bool]] = []
            for Lm, idxs in cand_map.items():
                B = ms_index.by_len[Lm]
                for j in idxs:
                    b = B[j]
                    d = edit_distance(a, b, k)
                    checked += 1
                    if d <= k:
                        # check vgene overlap
                        a_vs = agab_vmap.get(a, set())
                        b_vs = ms_vmap.get(b, set())
                        vgene_match = bool(a_vs and b_vs and (a_vs & b_vs))
                        local_rows.append((a, b, d,vgene_match,a_vs,b_vs))

            if not local_rows:
                continue

            # Partition by vgene match
            ok = [r for r in local_rows if r[3]]
            bad = [r for r in local_rows if not r[3]]

            # Count any hit for reporting
            hit_by_k[k] += 1

            if ok:
                # Prefer vgene-matching hits at this k and stop (like original behavior)
                rows.extend(ok)
                for (_, b, d, _, _, _) in ok:
                    if (b not in best_for_ms) or (d < best_for_ms[b][0]) or (d == best_for_ms[b][0] and a < best_for_ms[b][1]):
                        best_for_ms[b] = (d, a)
                break
            else:
                # No vgene-matching hit: keep non-matching hits (flagged) but continue to next k
                rows.extend(bad)
                for (_, b, d, _, _, _) in bad:
                    if (b not in best_for_ms) or (d < best_for_ms[b][0]) or (d == best_for_ms[b][0] and a < best_for_ms[b][1]):
                        best_for_ms[b] = (d, a)
                # continue to next k (do not break)

        if i % 1000 == 0:
            elapsed = time.time() - t0
            rate = i / max(elapsed, 1e-9)
            eta = (n_agab - i) / max(rate, 1e-9)
            k1, k2 = hit_by_k.get(1, 0), hit_by_k.get(2, 0)
            print(f"[INFO] processed {i}/{n_agab} AGAB | "
                  f"matches: {len(rows)} | checked: {checked} | "
                  f"elapsed: {elapsed:.1f}s | ETA: {eta/3600:.1f}h | "
                  f"k1={k1} k2={k2}",
                  flush=True)

    return rows, best_for_ms, checked
    # ms_index = MSIndex(ms_df)
    # rows: List[Tuple[str, str, int]] = []
    # best_for_ms: Dict[str, Tuple[int, str]] = {}
    # checked = 0
    # t0 = time.time()

    # agab_list = agab_df["cdr3_aa"].tolist()
    # n_agab = len(agab_list)
    # hit_by_k = defaultdict(int)

    # for i, a in enumerate(agab_list, 1):
    #     La = len(a)
    #     anchors = tokens_pms4(a)

    #     for k in range(1, max_lev + 1):
    #         Lms = [L for L in range(La - k, La + k + 1) if L in ms_index.by_len]
    #         if not Lms:
    #             continue

    #         cand_map = ms_index.candidates(Lms, anchors, min_hits=min_hits)
    #         if not cand_map:
    #             continue

    #         local_rows: List[Tuple[str, str, int]] = []
    #         for Lm, idxs in cand_map.items():
    #             B = ms_index.by_len[Lm]
    #             for j in idxs:
    #                 b = B[j]
    #                 d = edit_distance(a, b, k)
    #                 checked += 1
    #                 if d <= k:
    #                     local_rows.append((a, b, d))

    #         if local_rows:
    #             rows.extend(local_rows)
    #             hit_by_k[k] += 1
    #             for (_, b, d) in local_rows:
    #                 if (b not in best_for_ms) or (d < best_for_ms[b][0]) or (d == best_for_ms[b][0] and a < best_for_ms[b][1]):
    #                     best_for_ms[b] = (d, a)
    #             break

    #     if i % 1000 == 0:
    #         elapsed = time.time() - t0
    #         rate = i / max(elapsed, 1e-9)
    #         eta = (n_agab - i) / max(rate, 1e-9)
    #         k1, k2 = hit_by_k.get(1, 0), hit_by_k.get(2, 0)
    #         print(f"[INFO] processed {i}/{n_agab} AGAB | "
    #               f"matches: {len(rows)} | checked: {checked} | "
    #               f"elapsed: {elapsed:.1f}s | ETA: {eta/3600:.1f}h | "
    #               f"k1={k1} k2={k2}",
    #               flush=True)

    # return rows, best_for_ms, checked




def atomic_update_ms(ms_path: str,
                     best_for_ms: Dict[str, Tuple[int, str]],
                     ms_key_col: str = "cdr3_aa"):
    """
    Add 'agab_match' and 'lev' columns to MS-only CDR3 table and replace atomically.
    This is optional but mirrors the original behavior.
    """
    with FileLock(ms_path):
        # Read with explicit separator to avoid header being split
        if isinstance(ms_path, pd.DataFrame):
            ms_df_full = ms_path
        elif ms_path.endswith(".csv"):
            ms_df_full = pd.read_csv(ms_path, sep=",")
        elif ms_path.endswith(".tsv"):
            ms_df_full = pd.read_csv(ms_path, sep="\t")
        else:
            ms_df_full = pd.read_csv(ms_path, sep=None, engine="python")

        # Fix common header issue: 'cdr3_' instead of 'cdr3_aa'
        if ms_key_col not in ms_df_full.columns:
            if "cdr3_" in ms_df_full.columns:
                ms_df_full = ms_df_full.rename(columns={"cdr3_": ms_key_col})
            else:
                raise ValueError(
                    f"Column '{ms_key_col}' not found in {ms_path}. "
                    f"Columns: {list(ms_df_full.columns)}"
                )

        ms_norm = normalize_cdr3_series(ms_df_full[ms_key_col])

        agab_match_col, lev_col = [], []
        for s in ms_norm:
            if s in best_for_ms:
                lev, agab = best_for_ms[s]
                agab_match_col.append(agab)
                lev_col.append(lev)
            else:
                agab_match_col.append(pd.NA)
                lev_col.append(pd.NA)

        ms_df_full["agab_match"] = agab_match_col
        ms_df_full["lev"] = lev_col

        dir_ = os.path.dirname(ms_path) or "."
        base = os.path.basename(ms_path)
        ts = time.strftime("%Y%m%d_%H%M%S")
        backup = os.path.join(dir_, f".backup_{base}.{ts}")
        try:
            shutil.copy2(ms_path, backup)
            print(f"[INFO] backup created: {backup}")
        except Exception as e:
            print(f"[WARN] could not create backup: {e}")

        fd, tmp_path = tempfile.mkstemp(prefix=base + ".", suffix=".tmp", dir=dir_)
        os.close(fd)
        try:
            ms_df_full.to_csv(tmp_path, sep=",", index=False)
            os.replace(tmp_path, ms_path)
            print(f"[OK] updated MS-only CDR3 file in-place: {ms_path}")
        finally:
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except Exception:
                pass


# =========================
# Fisher helpers (antigen only)
# =========================

def _eta_msg(step_idx: int, total_steps: int, t0: float, step_name: str) -> None:
    """Print an ETA message after completing a step."""
    elapsed = time.perf_counter() - t0
    done = step_idx
    remain = max(total_steps - done, 0)
    avg = elapsed / done if done > 0 else float("nan")
    eta_sec = remain * avg if np.isfinite(avg) else float("nan")
    pct = 100.0 * done / total_steps

    def fmt(s):
        if not np.isfinite(s):
            return "?"
        m, s = divmod(int(s), 60)
        h, m = divmod(m, 60)
        return f"{h:02d}:{m:02d}:{s:02d}"

    print(f"[ETA] Step {done}/{total_steps} ({pct:4.1f}%) | {step_name:>28s} | "
          f"elapsed={fmt(elapsed)} | remaining≈{fmt(eta_sec)}")

def build_group_index(df: pd.DataFrame, group_col: str) -> Dict[str, Set[str]]:
    """
    Given df with columns [group_col, COLS_AGAB['cdr3']], return group -> set(cdr3_aa).
    """
    g = (
        df
        .dropna(subset=[group_col])
        .loc[:, [group_col, COLS_AGAB["cdr3"]]]
        .drop_duplicates()
    )
    out: Dict[str, Set[str]] = {}
    for name, sub in g.groupby(group_col):
        out[str(name)] = set(sub[COLS_AGAB["cdr3"]].tolist())
    return out

def fisher_and_residual(a: int, b: int, c: int, d: int) -> Tuple[float, float]:
    """
    Fisher exact test (one-sided 'greater') and standardized residual for cell a.
    """
    tbl = np.array([[a, b], [c, d]], dtype=int)
    _, p = fisher_exact(tbl, alternative="greater")
    n = tbl.sum()
    row_tot = tbl[0, :].sum()
    col_tot = tbl[:, 0].sum()
    if n == 0:
        return p, np.nan
    expected_a = (row_tot * col_tot) / n
    row_prop = row_tot / n
    col_prop = col_tot / n
    var_a = expected_a * (1 - row_prop) * (1 - col_prop)
    sr = (a - expected_a) / np.sqrt(var_a) if var_a > 0 else np.nan
    return p, float(sr)

def run_fisher_family(universe: Set[str],
                      in_data_set: Set[str],
                      group_index: Dict[str, Set[str]],
                      padj_method: str = "fdr_bh") -> pd.DataFrame:
    """
    Run Fisher enrichment for each group (e.g. antigen), given:
      - universe: all AGAB CDR3 eligible
      - in_data_set: AGAB CDR3 that matched MS at given lev
      - group_index: group -> set(CDR3)
    """
    rows: List[Dict] = []
    U, D = universe, in_data_set
    for group, G in group_index.items():
        a = len(D & G)
        b = len(D - G)
        c = len(G - D)
        d = len(U - (D | G))
        p, sr = fisher_and_residual(a, b, c, d)
        rows.append({
            "group": group,
            "a_in_data_in_group": a,
            "b_in_data_not_group": b,
            "c_not_data_in_group": c,
            "d_not_data_not_group": d,
            "p_raw": p,
            "std_residual": sr,
        })
    out = pd.DataFrame(rows).sort_values("p_raw").reset_index(drop=True)
    if len(out):
        out["p_adj"] = multipletests(out["p_raw"], method=padj_method)[1]
    else:
        out["p_adj"] = []
    return out.sort_values(["p_adj", "std_residual"], ascending=[True, False]).reset_index(drop=True)

def enrich_matches_with_metadata(matches: pd.DataFrame,
                                 ms_filtered: pd.DataFrame,
                                 agab: pd.DataFrame) -> pd.DataFrame:
    """
    Left-join MS and AgAb metadata into the matches table and keep all original match columns.
    Only antigen and origin_class are used from AgAb (no taxonomy).
    """
    m1 = matches.merge(
        ms_filtered[[COLS_MS["cdr3"], COLS_MS["subject"], COLS_MS["tissue"], COLS_MS["cluster"]]].drop_duplicates(),
        left_on=COLS_MATCH["cdr3_ms"],
        right_on=COLS_MS["cdr3"],
        how="left",
        suffixes=("", "_ms")
    )
    m2 = m1.merge(
        agab[[COLS_AGAB["cdr3"], COLS_AGAB["antigen"], COLS_AGAB["origin"]]].drop_duplicates(),
        left_on=COLS_MATCH["cdr3_agab"],
        right_on=COLS_AGAB["cdr3"],
        how="left",
        suffixes=("", "_agab")
    )
    drop_cols = [c for c in [COLS_MS["cdr3"], COLS_AGAB["cdr3"]] if c in m2.columns]
    m2 = m2.drop(columns=drop_cols, errors="ignore")
    front_cols = [
        COLS_MATCH["cdr3_agab"], COLS_MATCH["cdr3_ms"], COLS_MATCH["lev"],
        COLS_MS["subject"], COLS_MS["tissue"], COLS_MS["cluster"],
        COLS_AGAB["antigen"], COLS_AGAB["origin"],
    ]
    ordered_cols = front_cols + [c for c in m2.columns if c not in front_cols]
    return (m2.loc[:, ordered_cols])

def aggregate_tables(enriched: pd.DataFrame, out_prefix: str):
    """Create tabular summaries (antigen and tissue only)."""
    # Unique (per antigen × tissue) counts of unique AgAb CDR3
    antigen_by_tissue = (
        enriched
        .dropna(subset=[COLS_AGAB["antigen"], COLS_MS["tissue"], COLS_MATCH["cdr3_agab"]])
        .drop_duplicates(subset=[COLS_AGAB["antigen"], COLS_MS["tissue"], COLS_MATCH["cdr3_agab"]])
        .groupby([COLS_AGAB["antigen"], COLS_MS["tissue"]])
        .size()
        .reset_index(name="unique_agab_cdr3")
    )
    antigen_by_tissue.to_csv(f"{out_prefix}_antigen_by_tissue.tsv", sep="\t", index=False)

    # Tissue totals (unique AgAb CDR3 per tissue)
    tissue_summary = (
        enriched
        .dropna(subset=[COLS_MS["tissue"], COLS_MATCH["cdr3_agab"]])
        .drop_duplicates(subset=[COLS_MS["tissue"], COLS_MATCH["cdr3_agab"]])
        .groupby(COLS_MS["tissue"])
        .size()
        .reset_index(name="unique_agab_cdr3")
    )
    tissue_summary.to_csv(f"{out_prefix}_tissue_summary.csv", index=False)

    # Per-AgAb CDR3 roll-up with minimal lev and aggregated tissues/subjects
    rollup = (
        enriched
        .assign(lev_int=enriched[COLS_MATCH["lev"]].astype(int))
        .groupby(COLS_MATCH["cdr3_agab"])
        .agg(
            lev_min=("lev_int", "min"),
            n_ms_rows=("lev_int", "size"),
            n_ms_subjects=(COLS_MS["subject"], lambda x: x.dropna().nunique()),
            tissues=(COLS_MS["tissue"], lambda x: "|".join(sorted(pd.Series(x.dropna().unique(), dtype=str)))),
            subjects=(COLS_MS["subject"], lambda x: "|".join(sorted(pd.Series(x.dropna().unique(), dtype=str)))),
            clusters=(COLS_MS["cluster"], lambda x: "|".join(sorted(pd.Series(x.dropna().unique(), dtype=str)))),
            antigen=(COLS_AGAB["antigen"], lambda x: "|".join(sorted(pd.Series(x.dropna().unique(), dtype=str)))),
            origin_class=(COLS_AGAB["origin"], lambda x: "|".join(sorted(pd.Series(x.dropna().unique(), dtype=str)))),
        )
        .reset_index()
    )
    rollup.to_csv(f"{out_prefix}_agab_cdr3_rollup.tsv", sep="\t", index=False)

def compute_fisher_levels_antigen(agab_path: str, ms_df: str, matches_path: str, out_prefix: str):
    """
    Fisher test on antigen only, with lev 0/1/2 and enriched outputs.
    """
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    t0 = time.perf_counter()

    steps = [
        "load_agab", "prep_agab",
        "load_ms", "prep_ms", "filter_H_clusters",
        "load_matches", "prep_matches", "restrict_matches_to_ms",
        "index_antigen",
        "fisher_ant_lev0", "fisher_ant_lev1", "fisher_ant_lev2",
        "save_fisher",
        "enrich_matches", "save_enriched",
        "build_summaries", "save_summaries"
    ]
    step_i = 0
    total_steps = len(steps)
    COLS_AGAB["antigen"]='single_annotation'  # adjust for new column name
    # 1) Load AgAb annotated
    agab = pd.read_csv(agab_path, sep=None, engine="python")
    # agab= agab_df
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 1b) Prep AgAb
    for col in [COLS_AGAB["cdr3"], COLS_AGAB["antigen"], COLS_AGAB["origin"]]:
        if col not in agab.columns:
            raise ValueError(f"Column '{col}' not found in {agab_df}. Columns: {list(agab.columns)}")
    agab[COLS_AGAB["cdr3"]] = normalize_cdr3_series(agab[COLS_AGAB["cdr3"]])
    agab = agab[(agab[COLS_AGAB["cdr3"]].notna()) & (agab[COLS_AGAB["cdr3"]] != "")]
    agab["len"] = agab[COLS_AGAB["cdr3"]].str.len().astype("int32")
    agab = agab[agab["len"] >= MIN_LEN].copy()
    antigen_universe = set(
        agab.dropna(subset=[COLS_AGAB["antigen"]])[COLS_AGAB["cdr3"]].unique()
    )
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 2) Load MS/HC
    ms = pd.read_csv(ms_df, sep="\t")
    # ms=ms_df
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 2b) Prep MS
    for col in [COLS_MS["cdr3"], COLS_MS["subject"], COLS_MS["tissue"], COLS_MS["cluster"]]:
        if col not in ms.columns:
            raise ValueError(f"Column '{col}' not found in {ms_path}. Columns: {list(ms.columns)}")
    ms[COLS_MS["cdr3"]] = normalize_cdr3_series(ms[COLS_MS["cdr3"]])
    ms = ms[(ms[COLS_MS["cdr3"]].notna()) & (ms[COLS_MS["cdr3"]] != "")]
    ms["len"] = ms[COLS_MS["cdr3"]].str.len().astype("int32")
    ms = ms[ms["len"] >= MIN_LEN].copy()
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 2c) Filter H-clusters (same rule as MS_only1)
    ms_filt = filter_ms_exclude_H_clusters(ms)
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 3) Load matches
    matches = pd.read_csv(matches_path, sep="\t")
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 3b) Prep matches
    need = {COLS_MATCH["cdr3_agab"], COLS_MATCH["lev"]}
    if not need.issubset(matches.columns):
        raise ValueError(f"Matches file must have columns {need}. Columns: {list(matches.columns)}")
    if COLS_MATCH["cdr3_agab"] in matches.columns:
        matches[COLS_MATCH["cdr3_agab"]] = normalize_cdr3_series(matches[COLS_MATCH["cdr3_agab"]])
    if COLS_MATCH["cdr3_ms"] in matches.columns:
        matches[COLS_MATCH["cdr3_ms"]] = normalize_cdr3_series(matches[COLS_MATCH["cdr3_ms"]])
    matches = matches[matches[COLS_MATCH["lev"]].astype(int) <= MAX_LEV].copy()
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 3c) Restrict matches to filtered MS
    if COLS_MATCH["cdr3_ms"] in matches.columns:
        valid_ms = set(ms_filt[COLS_MS["cdr3"]].unique())
        matches = matches[matches[COLS_MATCH["cdr3_ms"]].isin(valid_ms)].copy()
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 4) Group indices (antigen) -- old, now Ags are from the matched cdr3s
    antigen_index = build_group_index(
        agab[[COLS_AGAB["cdr3"], COLS_AGAB["antigen"]]].copy(),
        COLS_AGAB["antigen"]
    )
# 4) Group indices (antigen) — restrict to matched AgAb CDR3s

    # matched_agab_cdr3s = set(matches[COLS_MATCH["cdr3_agab"]].unique())

    # agab_matched = agab[
    #     agab[COLS_AGAB["cdr3"]].isin(matched_agab_cdr3s)
    # ].copy()

    # antigen_universe = set(agab_matched[COLS_AGAB["cdr3"]].unique())

    # antigen_index = build_group_index(
    #     agab_matched[[COLS_AGAB["cdr3"], COLS_AGAB["antigen"]]],
    #     COLS_AGAB["antigen"]
    # )


    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # Build mapping antigen -> origin_class (pipe-separated if multiple)
    agab_antigen_origin = (
        agab
        .dropna(subset=[COLS_AGAB["antigen"], COLS_AGAB["origin"]])
        .drop_duplicates(subset=[COLS_AGAB["cdr3"], COLS_AGAB["antigen"], COLS_AGAB["origin"]])
        .groupby(COLS_AGAB["antigen"])[COLS_AGAB["origin"]]
        .apply(lambda x: "|".join(sorted(pd.Series(x.dropna().unique(), dtype=str))))
        .to_dict()
    )

    # 5) Fisher tests per lev
    def best_levels(sub: pd.DataFrame) -> pd.Series:
        return sub.groupby(COLS_MATCH["cdr3_agab"])[COLS_MATCH["lev"]].min().astype(int)

    best_antigen = best_levels(
        matches[matches[COLS_MATCH["cdr3_agab"]].isin(antigen_universe)]
    )
    in_ms_ant_lev0 = set(best_antigen.index[best_antigen == 0])
    in_ms_ant_lev1 = set(best_antigen.index[best_antigen == 1])
    in_ms_ant_lev2 = set(best_antigen.index[best_antigen == 2])

    r_ant_0 = run_fisher_family(antigen_universe, in_ms_ant_lev0, antigen_index, padj_method="fdr_bh")
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])
    r_ant_1 = run_fisher_family(antigen_universe, in_ms_ant_lev1, antigen_index, padj_method="fdr_bh")
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])
    r_ant_2 = run_fisher_family(antigen_universe, in_ms_ant_lev2, antigen_index, padj_method="fdr_bh")
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # Add origin_class annotation to Fisher tables (per antigen)
    for df in (r_ant_0, r_ant_1, r_ant_2):
        df["origin_class"] = df["group"].map(lambda g: agab_antigen_origin.get(g, pd.NA))

    # 6) Save Fisher tables (antigen only)
    r_ant_0.to_csv(f"{out_prefix}_antigen_lev0.csv", index=False)
    r_ant_1.to_csv(f"{out_prefix}_antigen_lev1.csv", index=False)
    r_ant_2.to_csv(f"{out_prefix}_antigen_lev2.csv", index=False)
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 7) Enrich matches and save
    enriched = enrich_matches_with_metadata(matches, ms_filt, agab)
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])
    enriched_out = f"{out_prefix}_matches_enriched.tsv"
    enriched.to_csv(enriched_out, sep="\t", index=False)
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    # 8) Summaries and save
    aggregate_tables(enriched, out_prefix)
    step_i += 1; _eta_msg(step_i, total_steps, t0, steps[step_i-1])

    print(f"[DONE] All Fisher/summary tables written with prefix: {out_prefix}")
    #print the significant antigens found at lev 0/1/2
    for df, lev in zip([r_ant_0, r_ant_1, r_ant_2], [0, 1, 2]):
        sig = df[df["p_adj"] < 0.05]
        n_sig = len(sig)
        print(f"[INFO] Lev {lev}: {n_sig} significant antigens (p_adj < 0.05)")
        if n_sig > 0:
            print(sig[["group", "p_adj", "std_residual"]].to_string(index=False))
    
    COLS_AGAB["antigen"]='antigen'  # reset to original


# =========================
# uniform antigen annotations
# =========================
import re
from typing import Iterable, List, Set, Dict, Any, Optional, Tuple
import pandas as pd


def normalize_antigen(text: str) -> str:
    """
    Lowercase, treat '_' and '-' as separators, remove punctuation except alphanumerics/spaces,
    and collapse whitespace.
    """
    if text is None:
        return ""
    text = text.lower()
    text = re.sub(r"[_\-]+", " ", text)          # separators -> space
    text = re.sub(r"[^a-z0-9 ]+", " ", text)     # keep only alnum + spaces
    text = re.sub(r"\s+", " ", text).strip()     # collapse whitespace
    return text


def tokenize_antigen(text: str, drop: Optional[Set[str]] = None) -> List[str]:
    """
    Normalize + split into tokens; drop generic tokens; return unique tokens (order not important).
    """
    if drop is None:
        drop = {"antigen", "protein", "proteins"}

    norm = normalize_antigen(text)
    toks = [t for t in norm.split(" ") if t]
    toks = [t for t in toks if t not in drop]
    # unique (stable-ish): use dict trick
    return list(dict.fromkeys(toks))


def jaccard_tokens(a_tokens: Iterable[str], b_tokens: Iterable[str]) -> float:
    """
    Jaccard similarity between two token sets.
    """
    a_set = set(a_tokens)
    b_set = set(b_tokens)
    if not a_set or not b_set:
        return 0.0
    inter = len(a_set.intersection(b_set))
    union = len(a_set.union(b_set))
    return inter / union


from itertools import combinations

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


def all_vs_all_jaccard(antigens: List[str], min_sim: float = 0.6, show_progress: bool = True) -> pd.DataFrame:
    """
    Compute all unordered pairs (i<j), compute Jaccard on token sets, and return pairs with sim>=min_sim.
    """
    # drop NA/empty and unique
    cleaned = []
    seen = set()
    for a in antigens:
        if a is None:
            continue
        a = str(a)
        if a == "":
            continue
        if a not in seen:
            cleaned.append(a)
            seen.add(a)

    toks = [tokenize_antigen(a) for a in cleaned]

    pairs = combinations(range(len(cleaned)), 2)
    if tqdm is not None and show_progress:
        pairs = tqdm(list(pairs), desc="All-vs-all Jaccard")  # materialize for tqdm
    else:
        pairs = list(pairs)

    rows = []
    for i, j in pairs:
        sim = jaccard_tokens(toks[i], toks[j])
        if sim >= min_sim:
            rows.append((cleaned[i], cleaned[j], sim))

    return pd.DataFrame(rows, columns=["a", "b", "sim"])


import networkx as nx


def canonicalize(strings: List[str]) -> str:
    """
    Pick a canonical representative for a cluster.
    Default rule: maximum token coverage (largest unique token set).
    Tie-breaker: longer normalized string.
    """
    if not strings:
        return ""
    token_lists = [tokenize_antigen(s) for s in strings]
    sizes = [len(set(t)) for t in token_lists]
    max_size = max(sizes)
    candidates = [s for s, sz in zip(strings, sizes) if sz == max_size]
    # tie-break by longer normalized length
    candidates.sort(key=lambda s: len(normalize_antigen(s)), reverse=True)
    return(candidates[0])


def build_clusters_with_canonical(pairs: pd.DataFrame, min_sim: float = 0.9) -> Dict[str, Any]:
    """
    - Build threshold graph from edges where sim>=min_sim
    - Connected components define clusters
    - Build canonical_map: canonical label and members per cluster
    - Compute avg_edge_sim (on threshold edges) and avg_allpairs_sim (within cluster, from full pairs)
    - Return graph + membership + canonical_map + stats
    """
    required_cols = {"a", "b", "sim"}
    if not required_cols.issubset(pairs.columns):
        raise ValueError(f"pairs must contain columns {required_cols}, got {pairs.columns.tolist()}")

    edges_thr = pairs.loc[pairs["sim"] >= min_sim, ["a", "b", "sim"]].copy()
    if edges_thr.empty:
        raise ValueError("No edges above min_sim. Lower min_sim or check your similarity computation.")

    # 1) Build graph from threshold edges
    g = nx.Graph()
    for a, b, sim in edges_thr.itertuples(index=False):
        g.add_edge(a, b, sim=float(sim))

    # 2) Connected components -> cluster ids
    components = list(nx.connected_components(g))
    membership_rows = []
    for cid, nodes in enumerate(components, start=1):
        for node in nodes:
            membership_rows.append((node, cid))
    membership = pd.DataFrame(membership_rows, columns=["antigen", "cluster"])

    # 3) Canonical map
    canonical_map_rows = []
    for cid, nodes in enumerate(components, start=1):
        members = sorted(nodes)
        canon = canonicalize(members)
        canonical_map_rows.append((cid, canon, members))
    canonical_map = pd.DataFrame(canonical_map_rows, columns=["cluster", "canonical", "members"])

    # Helper: map antigen -> cluster
    antigen_to_cluster = dict(membership[["antigen", "cluster"]].itertuples(index=False))

    # 4A) Avg EDGE sim within cluster (threshold edges only)
    edges_thr2 = edges_thr.copy()
    edges_thr2["cluster_a"] = edges_thr2["a"].map(antigen_to_cluster)
    edges_thr2["cluster_b"] = edges_thr2["b"].map(antigen_to_cluster)
    within_thr = edges_thr2.loc[edges_thr2["cluster_a"] == edges_thr2["cluster_b"], ["cluster_a", "sim"]]
    avg_edge_sim = (
        within_thr.groupby("cluster_a")
        .agg(n_edges=("sim", "size"), avg_edge_sim=("sim", "mean"))
        .reset_index()
        .rename(columns={"cluster_a": "cluster"})
    )

    # 4B) Avg ALL-PAIRS sim within cluster (from full pairs dataframe)
    pairs2 = pairs.loc[:, ["a", "b", "sim"]].copy()
    pairs2["cluster_a"] = pairs2["a"].map(antigen_to_cluster)
    pairs2["cluster_b"] = pairs2["b"].map(antigen_to_cluster)
    within_all = pairs2.loc[pairs2["cluster_a"] == pairs2["cluster_b"], ["cluster_a", "sim"]]
    avg_allpairs_sim = (
        within_all.groupby("cluster_a")
        .agg(n_pairs=("sim", "size"), avg_allpairs_sim=("sim", "mean"))
        .reset_index()
        .rename(columns={"cluster_a": "cluster"})
    )

    # Cluster sizes
    cluster_sizes = membership.groupby("cluster").size().reset_index(name="n_nodes")

    # Stats table
    stats = (
        cluster_sizes.merge(avg_edge_sim, on="cluster", how="left")
        .merge(avg_allpairs_sim, on="cluster", how="left")
        .merge(canonical_map[["cluster", "canonical"]], on="cluster", how="left")
        .sort_values(["n_nodes", "avg_allpairs_sim"], ascending=[False, False], na_position="last")
        .reset_index(drop=True)
    )

    return {
        "graph": g,
        "membership": membership,
        "canonical_map": canonical_map,
        "stats": stats,
    }

import pandas as pd
from typing import Dict

def build_antigen_to_canonical_map(res: dict) -> Dict[str, str]:
    """
    Build a dict mapping each antigen string to the canonical antigen of its cluster.
    """
    membership = res["membership"]          # columns: antigen, cluster
    canonical_map = res["canonical_map"]    # columns: cluster, canonical, members

    cluster_to_canonical = dict(
        canonical_map[["cluster", "canonical"]].itertuples(index=False, name=None)
    )

    antigen_to_cluster = dict(
        membership[["antigen", "cluster"]].itertuples(index=False, name=None)
    )

    antigen_to_canonical = {
        antigen: cluster_to_canonical[cluster]
        for antigen, cluster in antigen_to_cluster.items()
    }
    return antigen_to_canonical

def add_canonical_column(df0: pd.DataFrame, antigen_col: str, antigen_to_canonical: dict) -> pd.DataFrame:
    """
    Add a new column '<antigen_col>_canonical' with the canonical name from clustering.
    If an antigen was not in the graph, keep the original value.
    """
    df = df0.copy()
    new_col = f"{antigen_col}_canonical"
    df[new_col] = df[antigen_col].map(antigen_to_canonical).fillna(df[antigen_col])
    return df


# =========================
# Main
# =========================
MS_FULL_PATH = "/doctorai/chiarba/analysis/AT_THE_END/MS_HC_2026_h0h1.tsv"
# def main():
parser = argparse.ArgumentParser(description="Step 3: MS_only1 + AGAB–MS matching (lev<=2) + Fisher (antigen only).")
parser.add_argument("--agab", default=AGAB_PATH, help="Path to AGAB annotated table with 'cdr3_aa' and 'antigen'.")
parser.add_argument("--ms-full", default=MS_FULL_PATH, help="Full MS/HC table with 'cdr3_aa', 'subject_number', 'tissue', 'cluster_h1'.")
parser.add_argument("--ms-only", default=MS_ONLY_PATH, help="Output path for MS_only1 CDR3 table.")
parser.add_argument("--matches-out", default=MATCHES_PATH, help="TSV path for all retained matches (lev<=2).")
parser.add_argument("--out-prefix", default=OUT_PREFIX, help="Prefix for Fisher and summary outputs.")
parser.add_argument("--throttle-every", type=int, default=0, help="Sleep every N comparisons (0 = disable).")
parser.add_argument("--sleep-sec", type=float, default=0.005, help="Sleep seconds when throttling triggers.")
parser.add_argument("--min-hits", type=int, default=MIN_ANCHOR_HITS, help="Min shared 4-mer anchors to admit a candidate.")
args = parser.parse_args()

# 1) Build MS_only1
build_ms_only1(args.ms_full, args.ms_only)

# 2) Matching
setup_backend()
print(f"[INFO] Levenshtein backend = {BACKEND}")

ms_df_full = pd.read_csv(args.ms_full , sep="\t")



ms_only_df_full=filter_ms_exclude_H_clusters(ms_df_full)
agab_df_full = pd.read_csv(args.agab )

# agab_df = read_cdr3(args.agab, col=COLS_AGAB["cdr3"])
# ms_df   = read_cdr3(args.ms_only, col="cdr3_aa")

agab_df = read_cdr3(args.agab, col=COLS_AGAB["cdr3"],vgene_col='ighv_v_gene')
ms_df   = read_cdr3(ms_only_df_full, col="cdr3_aa",vgene_col='v_gene')
#generate a non-MS dataframe, the HC_df (healthy control), which is the difference between ms_df_full and ms_only_df_full
HC_df = filter_ms_exclude_H_clusters(ms_df_full, return_HC_df=True)
HC_df = read_cdr3(HC_df, col="cdr3_aa",vgene_col='v_gene')

print(f"[STEP 3.2] AGAB unique CDR3: {len(agab_df)} | MS_only1 unique CDR3: {len(ms_df)}")

t0 = time.time()
rows, best_for_ms, checked = match_all_ighv(
    agab_df, ms_df,
    max_lev=MAX_LEV,
    throttle_every=max(args.throttle_every, 0),
    sleep_sec=max(args.sleep_sec, 0.0),
    min_hits=max(args.min_hits, 1)
)
elapsed = time.time() - t0
print(f"[DONE MATCH] matches kept: {len(rows)} | pairs checked (after prefilter): {checked} | elapsed {elapsed/60:.1f} min")

res_df = (
    pd.DataFrame(rows, columns=["cdr3_agab", "cdr3_ms", "lev","vgene_match","vgene_agab","vgene_ms"])
    .sort_values(['vgene_match',"lev","cdr3_agab","cdr3_ms"])
    .reset_index(drop=True)
)
#cluster the antigen names
ag_jac=all_vs_all_jaccard(agab_df_full['single_annotation'], min_sim=0.3)
# ag_jac is your pairs dataframe with columns a,b,sim
res = build_clusters_with_canonical(ag_jac, min_sim=0.7)
antigen_to_canonical = build_antigen_to_canonical_map(res)
agab_df_full_clust = add_canonical_column(agab_df_full, antigen_col="single_annotation", antigen_to_canonical=antigen_to_canonical)
agab_df_full_clust['single_annotation_old']=agab_df_full_clust['single_annotation']
agab_df_full_clust['single_annotation']=agab_df_full_clust['single_annotation_canonical']

agab_df_full_clust.to_csv('/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/agab_cdr3_annotated3_disambiguated_CLUSTER.csv')



#match the cdr3 with the antigen name from agab full
agab_columns= ['cdr3_aa', 'antigen','single_annotation','patent_quality','dataset','ighv_v_gene','mut_stat']
res_df_metadata = res_df.merge(
    agab_df_full_clust[agab_columns].drop_duplicates(),
    left_on="cdr3_agab",
    right_on="cdr3_aa",
    how="left")


res_df.to_csv(args.matches_out, sep="\t", index=False)
print(f"[SAVED] all retained matches (lev<=2): {args.matches_out}")
res_df_metadata.to_csv(args.matches_out.replace(".tsv","_withAgabMetadata_IGHV_v2.tsv"), sep="\t", index=False)
print(f"[SAVED] file:  {args.matches_out.replace('.tsv','_withAgabMetadata_IGHV_v2.tsv')}")
      

# 3) Fisher (antigen only)
compute_fisher_levels_antigen(  
                                # args.agab,
                                '/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/agab_cdr3_annotated3_disambiguated_CLUSTER.csv',
                                args.ms_full,
                                "/doctorai/niccoloc/MS_db/MS_BCR/aggregating_test/agab_vs_ms_lev_matches1_withAgabMetadata_IGHV_v2.tsv",
                                args.out_prefix)

print("[ALL DONE]")
