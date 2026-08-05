"""
01_optimize_benchmark.py — Systematic benchmark diversity optimization.

Grid-searches over redundancy configs, sampling methods, feature dim combos,
and sample sizes N to identify the most diverse default benchmark configuration.
Run this first — 02_make_benchmark.py consumes its output (results.tsv) to
build the actual versioned release.

Sample sizes N are relative to each config's own post-redundancy pool size P:
N = round(f * P) for f in 0.25, 0.30, ..., 0.75 (11 points per config), instead
of a fixed absolute list shared across configs. This guarantees every config
gets the same number of N points regardless of pool size, and keeps N always
well inside the pool (max tested N = 0.75P).

This script only produces the raw grid (results.tsv/baseline_stages.tsv/plot);
it does not pick a winner. The actual redundancy-config -> method/dims -> N
selection is a 3-step hierarchical procedure implemented in
02_make_benchmark.py (see its module docstring and reports/methods.tex):
whole-pool diversity (both plain pairwise averages and local-redundancy
metrics like nearest-neighbor distance) was shown to just re-encode pool
size across the 27 configs, not an independent signal, so Step 1 fixes the
redundancy config by filter severity + a minimum pool size instead of by
D_composite. Steps 2-3 (method/dims, then N) do use D_composite, but
renormalized with bounds local to the winning config only, and Step 3 picks
the largest N within a small tolerance of the curve's peak rather than the
raw maximum (which collapses to the smallest N tested, an artifact of
M1_mnnd's mechanical N-dependence).

Outputs:
  pepxpro_benchmark/optimization/results.tsv
  pepxpro_benchmark/optimization/baseline_stages.tsv  — D_composite of the quality
    pool (pre-redundancy, pre-sampling) and of each config's post-redundancy pool
    (pre-sampling), normalized with the same bounds as results.tsv so the two
    files are directly comparable (e.g. to see how much redundancy reduction vs.
    sampling each contribute to the final diversity score).
  pepxpro_benchmark/optimization/plots/diversity_vs_N.png

Usage:
  python 01_optimize_benchmark.py            # 2 example configs (demo)
  python 01_optimize_benchmark.py --full     # all 27 configs

The grid itself (N fractions, sampling methods, feature-dimension combos,
redundancy presets, TM-score cutoffs) and the input summary TSV can all be
overridden from the command line without editing this file, e.g.:

  python 01_optimize_benchmark.py --full \\
      --summary-tsv my_summary.tsv \\
      --n-fractions 0.3,0.5,0.7 \\
      --methods stratified,kennard-stone \\
      --tm-thresholds 0.5,0.9 \\
      --rec-presets strict,medium \\
      --pep-presets strict \\
      --dims-combo dg,pep_len --dims-combo dg,pep_len,rec_len,n_contacts

Omitted flags fall back to the defaults used for the PepXPro Benchmark v1
optimization (see reports/Draft/figure1_distributions.ipynb, Figure 5).
"""

import sys
import math
import time
import csv
import argparse
from pathlib import Path

_here      = Path(__file__).resolve().parent
_root      = _here.parent
_gensample = _root.parent / "genSample"
_scrape_out = _root.parent / "scrape" / "outputs"
_fixed_cif_dir = _root.parent / "scrape" / "fixed_cif"
sys.path.insert(0, str(_gensample))

from modules.utils       import load_tsv, enrich_rows
from modules.filters     import filter_single_chain
from modules.sequence_redundancy   import (filter_receptor_redundancy,
                                           filter_peptide_redundancy, PRESETS)
from modules.structural_redundancy import (filter_tm_score_receptor_matrix,
                                           filter_tm_score_receptor)
from modules.sampling    import run_sampling, DIM_KEYS

# ---------------------------------------------------------------------------
# Parameter space
# ---------------------------------------------------------------------------

# Sample sizes N are relative to each config's own pool size: N = round(f * P).
N_FRACTIONS = [round(0.25 + 0.05 * i, 2) for i in range(11)]  # 0.25 .. 0.75


def n_values_for_pool(pool_size, fractions=N_FRACTIONS):
    """Sample sizes N for a pool of this size, as round(f * pool_size) for
    each fraction. Deduplicated and sorted (small pools can round-collide)."""
    return sorted({round(f * pool_size) for f in fractions})


METHODS = ["stratified", "kennard-stone", "kde", "cluster"]

DIM_COMBOS = {
    "D1": ["dg", "pep_len"],
    "D2": ["dg", "pep_len", "rec_len", "n_contacts"],
    "D3": ["dg", "pep_len", "pep_helix", "pep_coil", "rec_helix", "sasa_burial"],
}

# Methods that produce identical results regardless of dims beyond their effective
# dimensionality — only run them with D1 to avoid redundant computation.
SINGLE_DIM_METHODS = {"stratified", "kde"}

SEED = 42

# Step 1 (see 02_make_benchmark.py / reports/methods.tex): minimum pool size a
# redundancy config must retain to be eligible as the benchmark's default.
N_MIN_POOL = 150

# Fixed evaluation space for M1–M4 and M9 (always populated)
EVAL_DIMS = ["dg", "pep_len", "rec_len"]
EVAL_KEYS = [DIM_KEYS[d] for d in EVAL_DIMS]

OUT_DIR   = _root / "optimization"
PLOTS_DIR = OUT_DIR / "plots"

EXAMPLE_CONFIGS = [
    {"id": 1, "name": "strict_TM05", "rec": "strict", "pep": "strict", "tm": 0.5},
    {"id": 2, "name": "medium_TM07", "rec": "medium", "pep": "medium", "tm": 0.7},
]

TM_THRESHOLDS = [0.5, 0.7, 0.9]

METRIC_COLS = ["M1_mnnd", "M2_mpd", "M3_cov", "M4_cv_nnd",
               "M5_pep_div", "M6_rec_div", "M7_tm_div", "M8_struct_cov",
               "M9_wass", "M10_redund"]
LOWER_IS_BETTER = {"M4_cv_nnd", "M9_wass", "M10_redund"}

# Human-readable abbreviations for row labels, e.g. "RS-PS-TM0.5-StD1"
REC_ABBR    = {"strict": "RS", "medium": "RM", "lenient": "RL"}
PEP_ABBR    = {"strict": "PS", "medium": "PM", "lenient": "PL"}
METHOD_ABBR = {"stratified": "St", "kennard-stone": "KS", "kde": "Kd", "cluster": "Cl"}


def row_label(rec, pep, tm, method, dims_id):
    """e.g. row_label('strict', 'strict', 0.5, 'stratified', 'D1') -> 'RS-PS-TM0.5-StD1'"""
    return f"{REC_ABBR[rec]}-{PEP_ABBR[pep]}-TM{tm}-{METHOD_ABBR[method]}{dims_id}"


RESULT_COLS = (["config_id", "row_id", "label", "config_name", "rec_preset", "pep_preset",
                "tm_threshold", "pool_size", "method", "dims_id", "dims", "N"]
               + METRIC_COLS + ["D_composite", "Diversity_pool"])

# baseline_stages.tsv: pre-redundancy / pre-sampling reference points, scored with
# the exact same normalization bounds as results.tsv (see _composite_bounds) so
# D_composite is directly comparable across the two files.
BASELINE_COLS = (["config_id", "stage", "rec_preset", "pep_preset", "tm_threshold",
                  "ref_pool_size", "n"] + METRIC_COLS + ["D_composite"])


def all_configs(rec_presets=("strict", "medium", "lenient"),
                pep_presets=("strict", "medium", "lenient"),
                tm_thresholds=TM_THRESHOLDS):
    return [
        {"id": i, "name": f"rec{rec}_pep{pep}_TM{tm}", "rec": rec, "pep": pep, "tm": tm}
        for i, (rec, pep, tm) in enumerate(
            ((rec, pep, tm)
             for rec in rec_presets
             for pep in pep_presets
             for tm in tm_thresholds),
            start=1,
        )
    ]


# ---------------------------------------------------------------------------
# CLI grid-parameter parsing
# ---------------------------------------------------------------------------

def _parse_str_list(s):
    return [x.strip() for x in s.split(",") if x.strip()]


def _parse_float_list(s):
    return [float(x) for x in _parse_str_list(s)]


def _parse_dims_combos(combo_args):
    """Build the {"D1": [...], "D2": [...], ...} dict from repeated --dims-combo
    flags, or fall back to the built-in D1/D2/D3 combos if none were given."""
    if not combo_args:
        return dict(DIM_COMBOS)
    return {f"D{i+1}": _parse_str_list(combo) for i, combo in enumerate(combo_args)}


# ---------------------------------------------------------------------------
# Pool building
# ---------------------------------------------------------------------------

def build_pool_redundancy(pool_quality, config):
    rows = list(pool_quality)
    rows, _ = filter_receptor_redundancy(rows, PRESETS[config["rec"]]["receptor"])
    rows, _ = filter_peptide_redundancy(rows, PRESETS[config["pep"]]["peptide"])
    result = filter_tm_score_receptor_matrix(rows, config["tm"])
    if result is None:
        # Matrix incomplete for this pool (e.g. new entries since the matrix was
        # built) — fall back to on-the-fly TM-score computation, same as
        # 02_make_benchmark.py, instead of silently skipping the TM filter.
        result = filter_tm_score_receptor(rows, config["tm"], _fixed_cif_dir)
    return result


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def _eval_scaled(sample_rows, pool_rows):
    pool_v = [r for r in pool_rows   if all(r.get(k) is not None for k in EVAL_KEYS)]
    samp_v = [r for r in sample_rows if all(r.get(k) is not None for k in EVAL_KEYS)]
    if not pool_v:
        return [], []
    nd = len(EVAL_KEYS)
    means = [sum(r[k] for r in pool_v) / len(pool_v) for k in EVAL_KEYS]
    stds  = [max((sum((r[k] - means[i])**2 for r in pool_v) / len(pool_v))**0.5, 1e-8)
             for i, k in enumerate(EVAL_KEYS)]
    def normalize(rs):
        return [tuple((r[k] - means[i]) / stds[i] for i, k in enumerate(EVAL_KEYS))
                for r in rs if all(r.get(k) is not None for k in EVAL_KEYS)]
    return normalize(pool_v), normalize(samp_v)


def _d(a, b):
    return math.sqrt(sum((a[i] - b[i])**2 for i in range(len(a))))


def m1_mnnd(pts):
    n = len(pts)
    if n < 2: return 0.0
    return sum(min(_d(pts[i], pts[j]) for j in range(n) if j != i) for i in range(n)) / n


def m2_mpd(pts):
    n = len(pts)
    if n < 2: return 0.0
    dists = [_d(pts[i], pts[j]) for i in range(n) for j in range(i+1, n)]
    return sum(dists) / len(dists)


def m3_coverage(sample_rows, pool_rows, n_bins=5):
    pv = [r for r in pool_rows   if all(r.get(k) is not None for k in EVAL_KEYS)]
    sv = [r for r in sample_rows if all(r.get(k) is not None for k in EVAL_KEYS)]
    if not pv or not sv: return 0.0
    ranges = [(min(r[k] for r in pv), max(r[k] for r in pv)) for k in EVAL_KEYS]
    def cell(r):
        return tuple(min(int((r[k] - lo) / max(hi - lo, 1e-8) * n_bins), n_bins-1)
                     for k, (lo, hi) in zip(EVAL_KEYS, ranges))
    pool_cells = {cell(r) for r in pv}
    return len({cell(r) for r in sv} & pool_cells) / len(pool_cells)


def m4_cv_nnd(pts):
    n = len(pts)
    if n < 3: return 1.0
    nnds = [min(_d(pts[i], pts[j]) for j in range(n) if j != i) for i in range(n)]
    mu = sum(nnds) / n
    if mu == 0: return 1.0
    return (sum((d - mu)**2 for d in nnds) / n)**0.5 / mu


def _kmer(seq_a, seq_b, k):
    if not seq_a or not seq_b: return 0.0
    a = {seq_a[i:i+k] for i in range(max(1, len(seq_a)-k+1))}
    b = {seq_b[i:i+k] for i in range(max(1, len(seq_b)-k+1))}
    u = len(a | b)
    return len(a & b) / u if u > 0 else 0.0


def _biopython_identity(seq_a, seq_b, aligner):
    if not seq_a or not seq_b: return 0.0
    try:
        aln = next(iter(aligner.align(seq_a, seq_b)))
        aln_len = aln.shape[1]
        n_match = sum(seq_a[i] == seq_b[j] for i, j in zip(*aln.indices) if i >= 0 and j >= 0)
        return n_match / aln_len if aln_len > 0 else 0.0
    except Exception:
        return 0.0


def _make_aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices
    aln = PairwiseAligner()
    aln.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aln.open_gap_score   = -10
    aln.extend_gap_score = -0.5
    aln.mode = "global"
    return aln


def m5_pep_div(sample_rows, aligner):
    """1 - mean pairwise peptide identity (Biopython BLOSUM62, same as pipeline)."""
    seqs = [r.get("peptide_seq", "") for r in sample_rows if r.get("peptide_seq")]
    if len(seqs) < 2: return 0.0
    ids = [_biopython_identity(seqs[i], seqs[j], aligner)
           for i in range(len(seqs)) for j in range(i+1, len(seqs))]
    return 1.0 - sum(ids) / len(ids)


def _rec_seq(r):
    raw = r.get("receptor_seq", "")
    if not raw: return ""
    if "|" in raw: return "".join(p.split(":",1)[1] for p in raw.split("|") if ":" in p)
    return raw.split(":",1)[1] if ":" in raw else raw


def m6_rec_div(sample_rows):
    seqs = [_rec_seq(r) for r in sample_rows if _rec_seq(r)]
    if len(seqs) < 2: return 0.0
    sims = [_kmer(seqs[i], seqs[j], 5) for i in range(len(seqs)) for j in range(i+1, len(seqs))]
    return 1.0 - sum(sims) / len(sims)


def m7_tm_div(sample_rows, tm_mat, tm_idx):
    if tm_mat is None: return None
    import numpy as np
    idx = [tm_idx[r["pdb_id"]] for r in sample_rows if r["pdb_id"] in tm_idx]
    if len(idx) < 2: return None
    sub = tm_mat[np.ix_(idx, idx)]
    n   = len(idx)
    vals = [float(sub[i, j]) for i in range(n) for j in range(i+1, n)]
    return 1.0 - sum(vals) / len(vals)


def _pool_clusters(pool_rows, tm_mat, tm_idx, threshold=0.5):
    if tm_mat is None: return None
    import numpy as np
    pids = [r["pdb_id"] for r in pool_rows if r["pdb_id"] in tm_idx]
    if not pids: return None
    pidx = [tm_idx[p] for p in pids]
    mat  = tm_mat[np.ix_(pidx, pidx)]
    assigned, clusters = set(), []
    for i, pid in enumerate(pids):
        if pid in assigned: continue
        cl = [pid]; assigned.add(pid)
        for j, qid in enumerate(pids):
            if qid not in assigned and mat[i, j] >= threshold:
                cl.append(qid); assigned.add(qid)
        clusters.append(cl)
    return clusters


def m8_struct_cov(sample_rows, pool_clusters):
    if pool_clusters is None: return None
    sids = {r["pdb_id"] for r in sample_rows}
    covered = sum(1 for c in pool_clusters if any(p in sids for p in c))
    return covered / len(pool_clusters)


def m9_wasserstein(sample_rows, pool_rows):
    from scipy.stats import wasserstein_distance
    vals = []
    for k in EVAL_KEYS:
        pv = [r[k] for r in pool_rows   if r.get(k) is not None]
        sv = [r[k] for r in sample_rows if r.get(k) is not None]
        if not pv or not sv: continue
        rng = max(pv) - min(pv)
        vals.append(wasserstein_distance(pv, sv) / max(rng, 1e-8))
    return sum(vals) / len(vals) if vals else 0.0


def m10_redund(sample_rows, aligner, threshold=0.70):
    """Fraction of sample peptide pairs with identity > threshold (Biopython BLOSUM62)."""
    seqs = [r.get("peptide_seq", "") for r in sample_rows if r.get("peptide_seq")]
    if len(seqs) < 2: return 0.0
    pairs = sum(1 for i in range(len(seqs)) for j in range(i+1, len(seqs))
                if _biopython_identity(seqs[i], seqs[j], aligner) > threshold)
    total = len(seqs) * (len(seqs) - 1) // 2
    return pairs / total


def compute_metrics(sample, pool, pool_cl, tm_mat, tm_idx, aligner):
    _, spts = _eval_scaled(sample, pool)
    return {
        "M1_mnnd":       m1_mnnd(spts),
        "M2_mpd":        m2_mpd(spts),
        "M3_cov":        m3_coverage(sample, pool),
        "M4_cv_nnd":     m4_cv_nnd(spts),
        "M5_pep_div":    m5_pep_div(sample, aligner),
        "M6_rec_div":    m6_rec_div(sample),
        "M7_tm_div":     m7_tm_div(sample, tm_mat, tm_idx),
        "M8_struct_cov": m8_struct_cov(sample, pool_cl),
        "M9_wass":       m9_wasserstein(sample, pool),
        "M10_redund":    m10_redund(sample, aligner),
    }


# ---------------------------------------------------------------------------
# Composite score
# ---------------------------------------------------------------------------

def _composite_bounds(results):
    """Per-metric (lo, hi) used to min-max normalize D_composite. Computed once
    from the grid search rows so that any other row scored with these same
    bounds (e.g. the pre-redundancy / pre-sampling baseline stages) lands on
    an identical scale and is directly comparable to results.tsv."""
    bounds = {}
    for m in METRIC_COLS:
        vals = [r[m] for r in results if r.get(m) is not None]
        if not vals: continue
        lo, hi = min(vals), max(vals)
        bounds[m] = (lo, hi if hi != lo else lo + 1.0)
    return bounds


def apply_composite(results, bounds):
    for m in METRIC_COLS:
        if m not in bounds: continue
        lo, hi = bounds[m]
        span = hi - lo
        for r in results:
            v = r.get(m)
            if v is None:
                r[f"_n_{m}"] = None
            else:
                norm = (v - lo) / span
                r[f"_n_{m}"] = (1.0 - norm) if m in LOWER_IS_BETTER else norm
    for r in results:
        norms = [r[f"_n_{m}"] for m in METRIC_COLS if r.get(f"_n_{m}") is not None]
        r["D_composite"] = sum(norms) / len(norms) if norms else None


def add_composite(results):
    """Normalize D_composite using bounds derived from `results` itself."""
    apply_composite(results, _composite_bounds(results))


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(results, configs):
    print(f"\n{'='*72}")
    print("TOP 5 PER CONFIG  (by D_composite, best dims per N)")
    print(f"{'='*72}")
    for cfg in configs:
        sub = [r for r in results if r["config_name"] == cfg["name"]
               and r["D_composite"] is not None]
        if not sub: continue
        pool_n = sub[0]["pool_size"]
        print(f"\n[GC{cfg['id']}] {cfg['name']}  pool={pool_n}")
        print(f"  {'Method':<16} {'Dims':>4} {'N':>5}  {'D_comp':>6}  "
              f"{'M1_mnnd':>7}  {'M5_pep':>6}  {'M7_tm':>6}")
        for r in sorted(sub, key=lambda x: x["D_composite"], reverse=True)[:5]:
            tm = f"{r['M7_tm_div']:.3f}" if r["M7_tm_div"] is not None else "  N/A"
            print(f"  {r['method']:<16} {r['dims_id']:>4} {r['N']:>5}  "
                  f"{r['D_composite']:.3f}   {r['M1_mnnd']:.3f}    "
                  f"{r['M5_pep_div']:.3f}   {tm}")


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def _local_composite(rows):
    """Recompute D_composite using min-max bounds derived ONLY from `rows`
    (e.g. a single config's own grid rows), avoiding cross-config leakage
    from the global bounds baked into apply_composite()'s D_composite column
    (which spans all 27 configs). Adds a "D_local" key to each row in place
    and also returns the list for convenience."""
    bounds = {}
    for m in METRIC_COLS:
        vals = [r[m] for r in rows if r.get(m) is not None]
        if not vals:
            continue
        lo, hi = min(vals), max(vals)
        bounds[m] = (lo, hi if hi != lo else lo + 1.0)

    for r in rows:
        norms = []
        for m in METRIC_COLS:
            if m not in bounds or r.get(m) is None:
                continue
            lo, hi = bounds[m]
            norm = (r[m] - lo) / (hi - lo)
            norms.append((1 - norm) if m in LOWER_IS_BETTER else norm)
        r["D_local"] = sum(norms) / len(norms) if norms else None
    return rows


def _step1_winner(baseline_rows, n_min_pool=N_MIN_POOL):
    """Step 1: the most stringent (receptor, peptide, TM-score) combination
    whose non-redundant pool retains at least n_min_pool entries, ties broken
    by the larger margin above that floor. Operates on baseline_rows'
    "redundancy" stage (whole-pool diversity + pool size, computed directly
    on each config's pool with no sampling) -- NOT on the N-grid, and not on
    D_composite/whole-pool diversity, both shown to just re-encode pool size
    across the 27 configs (see reports/methods.tex). Lets main() resolve the
    winning config BEFORE running the expensive N-grid, so that grid only
    needs to run for one config instead of all of them.

    Returns the winning config_id (e.g. "GC2"), or None if no config
    qualifies."""
    strictness_rank  = {"strict": 0, "medium": 1, "lenient": 2}
    tm_severity_rank = {0.5: 0, 0.7: 1, 0.9: 2}

    qualified = [r for r in baseline_rows
                 if r["stage"] == "redundancy" and r["n"] >= n_min_pool]
    if not qualified:
        return None

    def severity(r):
        return (strictness_rank[r["rec_preset"]] + strictness_rank[r["pep_preset"]]
                + tm_severity_rank[r["tm_threshold"]])

    qualified.sort(key=lambda r: (severity(r), -r["n"]))
    return qualified[0]["config_id"]


def _select_winner(results, n_min_pool=N_MIN_POOL, tolerance=0.02):
    """Same 3-step hierarchical selection as 02_make_benchmark.py's
    _select_winning_config / _best_run_for_config, duplicated here only to
    annotate the winning point on this diagnostic plot. 02_make_benchmark.py
    remains the source of truth for the actual benchmark build; see its
    module docstring and reports/methods.tex for the full rationale."""
    strictness_rank  = {"strict": 0, "medium": 1, "lenient": 2}
    tm_severity_rank = {0.5: 0, 0.7: 1, 0.9: 2}

    by_cfg = {}
    for r in results:
        cid = r["config_id"]
        if cid not in by_cfg:
            by_cfg[cid] = {"config_id": cid, "config_name": r["config_name"],
                           "rec_preset": r["rec_preset"], "pep_preset": r["pep_preset"],
                           "tm_threshold": r["tm_threshold"], "pool_size": r["pool_size"]}

    qualified = [c for c in by_cfg.values() if c["pool_size"] >= n_min_pool]
    if not qualified:
        return None

    def severity(c):
        return (strictness_rank[c["rec_preset"]] + strictness_rank[c["pep_preset"]]
                + tm_severity_rank[c["tm_threshold"]])

    qualified.sort(key=lambda c: (severity(c), -c["pool_size"]))
    winner_cfg = qualified[0]

    sub = [r for r in results if r["config_id"] == winner_cfg["config_id"]
           and r["D_composite"] is not None]
    sub = _local_composite(sub)

    by_md = {}
    for r in sub:
        by_md.setdefault((r["method"], r["dims_id"]), []).append(r)

    means = {k: sum(r["D_local"] for r in v) / len(v) for k, v in by_md.items()}
    best_method, best_dims = max(means, key=means.get)

    curve = by_md[(best_method, best_dims)]
    d_max = max(r["D_local"] for r in curve)
    qualifying = [r for r in curve if r["D_local"] >= (1 - tolerance) * d_max]
    target = max(qualifying, key=lambda r: r["N"])

    return {"config_name": winner_cfg["config_name"], "method": best_method,
            "dims_id": best_dims, "N": target["N"], "D_composite": target["D_local"]}


def plot_results(results, configs):
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="ticks", context="paper", font_scale=0.85,
                  rc={"axes.spines.top": False, "axes.spines.right": False})

    colors = {
        "stratified":      "#4CAF50",
        "kennard-stone":   "#E53935",
        "kde":             "#FDD835",
        "cluster":         "#3F51B5",
    }
    # More dims = stronger/darker tone of the method's base color.
    dims_shade = {"D1": 0.4, "D2": 0.7, "D3": 1.0}
    # For methods with more than one dims combo (KS, cluster), also vary dash
    # length with dims: shorter dashes = fewer dims sampled, continuous = most.
    # stratified/kde only ever run D1 and keep a plain solid trend line
    # (no dims to contrast).
    dims_linestyle = {"D1": (0, (1, 1)), "D2": (0, (4, 2)), "D3": "-"}
    MULTI_DIM_METHODS = {"kennard-stone", "cluster"}

    def _shaded(hex_color, dims_id):
        factor = dims_shade.get(dims_id, 1.0)
        rgb = mcolors.to_rgb(hex_color)
        return tuple(1 - (1 - c) * factor for c in rgb)

    _abbr = {"recstrict": "RS", "recmedium": "RM", "reclenient": "RL",
             "pepstrict": "PS", "pepmedium": "PM", "peplenient": "PL"}

    def _short(name):
        parts = name.split("_")
        return " ".join(_abbr.get(p, p) for p in parts)

    winner = _select_winner(results)

    n_cfg = len(configs)
    n_cols = 5
    n_rows = math.ceil(n_cfg / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4.2 * n_cols, 3.8 * n_rows),
                             squeeze=False, constrained_layout=True)

    for idx, cfg in enumerate(configs):
        row, col = divmod(idx, n_cols)
        ax = axes[row][col]
        sub = [r for r in results if r["config_name"] == cfg["name"]
               and r["D_composite"] is not None]
        pool_n = int(float(str(sub[0]["pool_size"]))) if sub else "?"
        # Bounds LOCAL to this config only -- each subplot is independently
        # renormalized, so curve *shapes* within a config are comparable
        # (which method/dims/N performs best there) but the y-axis is no
        # longer on a common scale across subplots. That cross-config
        # comparison was already misleading (D_composite was shown to just
        # re-encode pool size across the 27 configs; see reports/methods.tex)
        # so nothing meaningful is lost, and this now matches Figure 5 Panel D.
        sub = _local_composite(sub)

        handles = {}
        for method in METHODS:
            dims_ids = sorted({r["dims_id"] for r in sub if r["method"] == method})
            for dims_id in dims_ids:
                pts = sorted((r["N"], r["D_local"]) for r in sub
                             if r["method"] == method and r["dims_id"] == dims_id)
                if not pts:
                    continue
                ns = [p[0] for p in pts]
                ds = [p[1] for p in pts]
                ls = dims_linestyle.get(dims_id, "-") if method in MULTI_DIM_METHODS else "-"
                line, = ax.plot(ns, ds, marker="o", linestyle=ls,
                                 color=_shaded(colors[method], dims_id),
                                 linewidth=1.3, markersize=3.5, alpha=0.9)
                handles[f"{method} ({dims_id})"] = line

        title = f"GC{cfg['id']}  {_short(cfg['name'])}\npool={pool_n}"
        if winner and cfg["name"] == winner["config_name"]:
            title += "  (winner)"
            star = ax.scatter([winner["N"]], [winner["D_composite"]], marker="*",
                               s=220, color="black", zorder=10, edgecolors="white",
                               linewidths=0.5)
            handles[f"winner: {winner['method']} ({winner['dims_id']}) N={winner['N']}"] = star

        ax.set_title(title, fontsize=7.5, fontweight="bold", pad=3)
        ax.set_xlabel("N", fontsize=7)
        ax.set_ylabel(r"$D_{\mathrm{composite}}$" if col == 0 else "", fontsize=7)
        ax.tick_params(labelsize=6.5)
        if idx == n_cfg - 1:
            ax.legend(handles.values(), handles.keys(), fontsize=5, loc="lower right",
                      framealpha=0.7, ncol=1)
        sns.despine(ax=ax)

    # hide unused axes in last row
    for empty in range(n_cfg, n_rows * n_cols):
        r, c = divmod(empty, n_cols)
        axes[r][c].set_visible(False)

    fig.suptitle(r"$D_{\mathrm{composite}}$ vs N — best dims per method per config"
                 " (bounds local to each config; black star = Step-1/2/3 winner)",
                 fontsize=9, fontweight="bold")
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out = PLOTS_DIR / "diversity_vs_N.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot → {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Grid-search benchmark redundancy/sampling configurations "
                     "for diversity (D_composite).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--full", action="store_true",
                        help="Consider the full 27-config redundancy grid instead of 2 "
                             "examples for Step 1. By default, only the Step-1 winner "
                             "(see reports/methods.tex) then gets its N x method x dims "
                             "grid searched -- pass --full-diagnostic to grid-search all "
                             "27 anyway.")
    parser.add_argument("--full-diagnostic", action="store_true",
                        help="Grid-search N x method x dims for ALL configs, not just "
                             "the Step-1 winner (~27x more runs). Only needed to "
                             "reproduce the full exploratory comparison across configs "
                             "(e.g. diversity_vs_N.png with every panel populated, or "
                             "the robustness check that most configs favor the same "
                             "method/dims as the winner). The actual benchmark build "
                             "(02_make_benchmark.py) only ever needs the winning "
                             "config's grid. Implies --full.")
    parser.add_argument("--plot-only", action="store_true",
                        help="Skip grid search; regenerate plot from existing results.tsv.")
    parser.add_argument("--summary-tsv", metavar="PATH", default="summary_merged.tsv",
                        help="Input summary TSV (default: summary_merged.tsv). Relative "
                             "paths are resolved against scrape/outputs/.")
    parser.add_argument("--n-fractions", metavar="F1,F2,...",
                        default=",".join(str(f) for f in N_FRACTIONS),
                        help="Comma-separated sample sizes N as a fraction of each "
                             "config's own pool size (default: 0.25-0.75 step 0.05).")
    parser.add_argument("--methods", metavar="M1,M2,...",
                        default="stratified,kennard-stone,kde,cluster",
                        help="Comma-separated sampling methods (subset of: "
                             "stratified, kennard-stone, kde, cluster).")
    parser.add_argument("--dims-combo", metavar="DIM1,DIM2,...", action="append",
                        help="Feature-dimension combo, e.g. 'dg,pep_len'. Repeat for "
                             "multiple combos (labeled D1, D2, ... in the order given). "
                             "Default: the built-in D1/D2/D3 combos.")
    parser.add_argument("--tm-thresholds", metavar="T1,T2,...", default="0.5,0.7,0.9",
                        help="Comma-separated receptor TM-score cutoffs.")
    parser.add_argument("--rec-presets", metavar="P1,P2,...", default="strict,medium,lenient",
                        help="Comma-separated receptor sequence-identity presets.")
    parser.add_argument("--pep-presets", metavar="P1,P2,...", default="strict,medium,lenient",
                        help="Comma-separated peptide sequence-identity presets.")
    args = parser.parse_args()

    n_fractions = _parse_float_list(args.n_fractions)
    methods     = _parse_str_list(args.methods)
    dim_combos  = _parse_dims_combos(args.dims_combo)
    tm_thresh   = _parse_float_list(args.tm_thresholds)
    rec_presets = _parse_str_list(args.rec_presets)
    pep_presets = _parse_str_list(args.pep_presets)

    summary_path = Path(args.summary_tsv)
    if not summary_path.is_absolute():
        summary_path = _scrape_out / summary_path

    if args.plot_only:
        results_path = OUT_DIR / "results.tsv"
        if not results_path.exists():
            print(f"ERROR: {results_path} not found — run without --plot-only first.")
            sys.exit(1)
        with open(results_path) as fh:
            results = list(csv.DictReader(fh, delimiter="\t"))
        for r in results:
            for k in ("N", "pool_size"):
                r[k] = int(r[k])
            for k in RESULT_COLS:
                if k not in ("row_id", "config_id", "label", "config_name", "rec_preset",
                             "pep_preset", "dims", "dims_id", "method"):
                    try: r[k] = float(r[k]) if r[k] not in ("", "None") else None
                    except: pass
        # Configs are derived from whatever grid actually produced results.tsv,
        # rather than re-run through all_configs(), so --plot-only stays correct
        # even after a custom (non-default) grid search.
        seen = {}
        for r in results:
            cid = r["config_id"]
            if cid not in seen:
                seen[cid] = {"id": int(cid[2:]), "name": r["config_name"],
                             "rec": r["rec_preset"], "pep": r["pep_preset"],
                             "tm": r["tm_threshold"]}
        configs = sorted(seen.values(), key=lambda c: c["id"])
        plot_results(results, configs)
        return

    args.full = args.full or args.full_diagnostic
    configs = all_configs(rec_presets, pep_presets, tm_thresh) if args.full else EXAMPLE_CONFIGS

    print("\nPepXPro benchmark optimization")
    print("=" * 40)
    rows = enrich_rows(load_tsv(summary_path))
    print(f"Loaded:          {len(rows)} entries  (from {summary_path.name})")
    pool_quality = filter_single_chain(rows)
    print(f"Single-chain:    {len(pool_quality)} entries\n")

    # Biopython aligner — same config as sequence_redundancy.py
    aligner = _make_aligner()
    print("Biopython aligner ready (BLOSUM62 global)\n")

    # Load TM matrix once
    tm_mat, tm_idx = None, {}
    if (_scrape_out / "tm_score_matrix.npz").exists():
        import numpy as np
        tm_idx = {pid: i for i, pid in
                  enumerate((_scrape_out / "tm_score_ids.txt").read_text().splitlines())}
        tm_mat = np.load(_scrape_out / "tm_score_matrix.npz")["matrix"]
        print(f"TM matrix:       {len(tm_idx)} entries\n")

    # Clusters of the quality pool itself — reference for the pre-redundancy /
    # pre-sampling baseline stages (see baseline_stages.tsv below).
    pc_quality = _pool_clusters(pool_quality, tm_mat, tm_idx)

    # Build redundancy pools (cheap -- filtering only, no sampling)
    print(f"Building {len(configs)} pool(s)...")
    pools, pool_cl = {}, {}
    for cfg in configs:
        print(f"\n  [{cfg['name']}]")
        pool = build_pool_redundancy(pool_quality, cfg)
        pools[cfg["name"]]   = pool
        pool_cl[cfg["name"]] = _pool_clusters(pool, tm_mat, tm_idx)
        print(f"  → pool: {len(pool)} entries")

    # Whole-pool baseline metrics (Step 1 candidate data): cheap, one row per
    # config, computed directly on the pool with no sampling. Computed BEFORE
    # the expensive N-grid so Step 1 can resolve the winning config first
    # (see reports/methods.tex) and the grid only needs to run for that one
    # config, instead of all of them.
    print("\nBaseline stages (pre-redundancy, pre-sampling)...")
    baseline_rows = []
    mets_quality = compute_metrics(pool_quality, pool_quality, pc_quality, tm_mat, tm_idx, aligner)
    baseline_rows.append({
        "config_id": "-", "stage": "quality", "rec_preset": "-", "pep_preset": "-",
        "tm_threshold": "-", "ref_pool_size": len(pool_quality), "n": len(pool_quality),
        **mets_quality,
    })
    for cfg in configs:
        pool = pools[cfg["name"]]
        mets_redund = compute_metrics(pool, pool_quality, pc_quality, tm_mat, tm_idx, aligner)
        baseline_rows.append({
            "config_id": f"GC{cfg['id']}", "stage": "redundancy",
            "rec_preset": cfg["rec"], "pep_preset": cfg["pep"], "tm_threshold": cfg["tm"],
            "ref_pool_size": len(pool_quality), "n": len(pool),
            **mets_redund,
        })

    # Step 1's whole-pool diversity, keyed by config_id -- raw (un-normalized)
    # mean of M5/M6/M7/(1-M10), each already a natural [0,1] fraction so none
    # needs min-max rescaling to its own tiny range across the 27 configs
    # (see reports/methods.tex). Reported in results.tsv too (one constant
    # value per config, repeated across its grid rows) so every number behind
    # the 3-step methodology lives in one file alongside D_composite.
    diversity_pool_by_cfg = {
        r["config_id"]: (r["M5_pep_div"] + r["M6_rec_div"] + r["M7_tm_div"]
                          + (1 - r["M10_redund"])) / 4
        for r in baseline_rows if r["stage"] == "redundancy"
    }

    # Step 1: resolve the winning config now, before the N-grid. Skipped
    # (grid all configs) for the 2-example demo or --full-diagnostic.
    if args.full and not args.full_diagnostic:
        winner_id = _step1_winner(baseline_rows)
        if winner_id is None:
            print(f"WARNING: no config retains >= {N_MIN_POOL} entries; "
                  f"falling back to grid-searching all {len(configs)} configs.")
            configs_to_grid = configs
        else:
            configs_to_grid = [cfg for cfg in configs if f"GC{cfg['id']}" == winner_id]
            print(f"\nStep 1 winner: {winner_id} — grid-searching only this config "
                  f"(pass --full-diagnostic to grid-search all {len(configs)} instead).")
    else:
        configs_to_grid = configs

    # Grid search — N is relative to each config's own pool size (0.25P-0.75P)
    def _active_combos(pool_size):
        for method in methods:
            dim_ids = ["D1"] if method in SINGLE_DIM_METHODS else list(dim_combos)
            for dims_id in dim_ids:
                for N in n_values_for_pool(pool_size, n_fractions):
                    yield method, dims_id, N

    total = sum(1 for cfg in configs_to_grid
                for _ in _active_combos(len(pools[cfg["name"]])))
    print(f"\nGrid search: {total} runs\n" + "-" * 60)

    results, n_done = [], 0
    for cfg in configs_to_grid:
        pool = pools[cfg["name"]]
        pc   = pool_cl[cfg["name"]]
        pool_n_values = n_values_for_pool(len(pool), n_fractions)
        for method in methods:
            dim_ids = ["D1"] if method in SINGLE_DIM_METHODS else list(dim_combos)
            for dims_id in dim_ids:
                dims = dim_combos[dims_id]
                eff_dims = [dims[0]] if method == "kde" else dims
                for N in pool_n_values:
                    t0 = time.time()
                    try:
                        sample = run_sampling(pool, N, method=method,
                                             dims=eff_dims, seed=SEED)
                    except Exception as e:
                        print(f"  SKIP {cfg['name']} {method} {dims_id} N={N}: {e}")
                        continue
                    mets = compute_metrics(sample, pool, pc, tm_mat, tm_idx, aligner)
                    n_done += 1
                    row = {"config_id": f"GC{cfg['id']}", "row_id": f"R{n_done}",
                           "label": row_label(cfg["rec"], cfg["pep"], cfg["tm"], method, dims_id),
                           "config_name": cfg["name"],
                           "rec_preset": cfg["rec"], "pep_preset": cfg["pep"],
                           "tm_threshold": cfg["tm"], "pool_size": len(pool),
                           "method": method, "dims_id": dims_id,
                           "dims": "+".join(eff_dims), "N": N,
                           "D_composite": None,
                           "Diversity_pool": diversity_pool_by_cfg[f"GC{cfg['id']}"]}
                    row.update(mets)
                    results.append(row)
                    elapsed = time.time() - t0
                    tm_str = f"{mets['M7_tm_div']:.3f}" if mets["M7_tm_div"] is not None else " N/A"
                    print(f"  [{n_done:3d}/{total}] {cfg['name']:<22} {method:<16} "
                          f"{dims_id} N={N:3d}  "
                          f"M1={mets['M1_mnnd']:.3f}  M5={mets['M5_pep_div']:.3f}  "
                          f"M7={tm_str}  ({elapsed:.1f}s)")

    bounds = _composite_bounds(results)
    apply_composite(results, bounds)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIR / "results.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=RESULT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(results)
    print(f"\nResults → {OUT_DIR / 'results.tsv'}  ({len(results)} rows)")

    # baseline_rows (whole-pool diversity, all configs) was already computed
    # above, before the N-grid -- just rescale it with the grid's bounds
    # (informational only: D_composite here isn't used by Step 1, which
    # works off the raw M5/M6/M7/M10 values, not this normalized score) and
    # write it out.
    apply_composite(baseline_rows, bounds)

    with open(OUT_DIR / "baseline_stages.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=BASELINE_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(baseline_rows)
    print(f"Baseline  → {OUT_DIR / 'baseline_stages.tsv'}  ({len(baseline_rows)} rows)")

    print_summary(results, configs_to_grid)
    plot_results(results, configs_to_grid)


if __name__ == "__main__":
    main()
