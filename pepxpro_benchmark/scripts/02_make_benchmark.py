"""
02_make_benchmark.py — Generate a versioned pepxpro_benchmark release.

Uses genSample modules directly.  No run_XX is created in genSample/outputs/.
Run 01_optimize_benchmark.py --full first — this script consumes its output
(pepxpro_benchmark/optimization/results.tsv).

Two ways to pick a configuration:

  --profile GC2          Frozen, citation-safe profile (see PROFILES below).
                          Reproduces an exact past release regardless of what
                          optimization/results.tsv currently contains. Once a
                          --from-config release is adopted as official, its
                          resolved parameters should be copied into PROFILES
                          so it stays reproducible even if the grid search is
                          re-run later with different data.

  --from-config auto     Derive the configuration live from results.tsv,
                          reproducing the Figure 5 hierarchical selection rule
                          (see reports/methods.tex):

                          Step 1 — redundancy config: whole-pool diversity
                          (both plain pairwise averages and local-redundancy
                          metrics like nearest-neighbor distance) was shown to
                          just re-encode pool size, not an independent
                          diversity signal, so the config is NOT picked by
                          D_composite at all. Instead: the most stringent
                          (receptor, peptide, TM-score) combination whose pool
                          retains at least N_MIN_POOL entries, ties broken by
                          the larger margin above that floor.

                          Step 2 — method + feature space: within that fixed
                          pool, the (method, dims) combination with the
                          highest mean D_composite across its own N curve,
                          renormalized with bounds local to this config only
                          (not the global grid bounds) to avoid cross-config
                          leakage.

                          Step 3 — sample size N: on the winning curve, the
                          largest N whose D_composite is within N_TOLERANCE of
                          that curve's own peak — not the raw maximum, which
                          collapses to the smallest N in the range due to the
                          mechanical M1_mnnd-vs-N artifact (see methods.tex).

  --from-config GC11     Target a specific grid config directly, regardless
                          of whether it "wins" Step 1 — e.g. to try an
                          alternative to the adopted default. Steps 2-3 still
                          apply to pick its method/dims/N.

Usage:
  python 02_make_benchmark.py --profile GC2                 # frozen v1 release
  python 02_make_benchmark.py --from-config auto --version v1
  python 02_make_benchmark.py --from-config GC11 --version v2
  python 02_make_benchmark.py --profile GC2 --dry-run        # show pool size only
"""

import argparse
import csv
import json
import shlex
import subprocess
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

_here       = Path(__file__).resolve().parent
_bench_root = _here.parent
_repo_root  = _bench_root.parent
_scrape_out = _repo_root / "scrape" / "outputs"
_fixed_cif  = _repo_root / "scrape" / "fixed_cif"
_filtered   = _scrape_out / "filtered_systems.tsv"

sys.path.insert(0, str(_repo_root / "genSample"))
from modules.utils import load_tsv, enrich_rows, save_tsv
from modules.filters import apply_quality_filters
from modules.sequence_redundancy import apply_redundancy_filters
from modules.structural_redundancy import apply_structural_redundancy_filters
from modules.sampling import run_sampling, method_label
from modules.export import write_sample, link_cifs, write_fasta


# ── Benchmark profiles ───────────────────────────────────────────────────────
# One entry per grid config selected in the optimization (Figures 5-6):
# rec/pep sequence-identity presets + TM-score cutoff, and the method/dims/N
# that scored best (D_composite) for that config.
#
# D2 = ["dg", "pep_len", "rec_len", "n_contacts"]  (pepxpro_benchmark/optimization)

SUMMARY_TSV = _scrape_out / "summary_merged.tsv"

PROFILES = {
    "GC2": {
        "version":         "v1",
        "seq_id_receptor": "strict",   # rec_preset — receptor 40% CD-HIT
        "seq_id_peptide":  "strict",   # pep_preset — peptide 70% Biopython pairwise
        "tm_score":        0.7,
        "sampling_method": "kennard-stone",
        "dims":            ["dg", "pep_len", "rec_len", "n_contacts"],  # D2
        "sample_n":        100,        # N=100, selected by the Step 3 tolerance rule
                                        # (largest N within 2% of the peak D_composite
                                        # on the winning kennard-stone/D2 curve for GC2;
                                        # see reports/methods.tex). Frozen here for
                                        # citation-safety independent of results.tsv.
                                        # Resolved under PDBbind v2020R1 (pool=167;
                                        # N_MIN_POOL raised to 150 for Step 1 after
                                        # v2020R1 made GC1 newly eligible at the
                                        # original 100-entry floor).
        "seed":            42,
    },
}

SINGLE_CHAIN  = True
FULL_METADATA = True
EXPORT_FASTA  = True

POOL_COLS = ["pdb_id", "receptor", "peptide", "dG_kcal_mol", "binding_type",
             "resolution", "peptide_ss_class", "receptor_func_class",
             "peptide_pct_helix", "peptide_pct_beta", "peptide_pct_coil",
             "receptor_pct_helix", "receptor_pct_beta", "receptor_pct_coil",
             "n_contacts", "n_hbonds", "pep_sasa_burial"]
# ─────────────────────────────────────────────────────────────────────────────

OPTIMIZATION_TSV = _bench_root / "optimization" / "results.tsv"
BASELINE_TSV     = _bench_root / "optimization" / "baseline_stages.tsv"
STRICTNESS_RANK  = {"strict": 0, "medium": 1, "lenient": 2}
TM_SEVERITY_RANK = {0.5: 0, 0.7: 1, 0.9: 2}

# Must match METRIC_COLS / LOWER_IS_BETTER in 01_optimize_benchmark.py exactly —
# results.tsv carries these raw metric columns alongside the pre-computed
# (globally-normalized) D_composite column.
METRIC_COLS = ["M1_mnnd", "M2_mpd", "M3_cov", "M4_cv_nnd",
               "M5_pep_div", "M6_rec_div", "M7_tm_div", "M8_struct_cov",
               "M9_wass", "M10_redund"]
LOWER_IS_BETTER = {"M4_cv_nnd", "M9_wass", "M10_redund"}

# Step 1: minimum pool size a redundancy config must retain to be eligible.
N_MIN_POOL = 150
# Step 3: N is the largest sample size within this fraction of the winning
# curve's peak D_composite (avoids the raw-maximum collapsing to the
# smallest N tested, an artifact of M1_mnnd's mechanical N-dependence).
N_TOLERANCE = 0.02


def _invoked_command():
    """The exact command line that produced this benchmark release, for reproducibility."""
    return "python " + Path(sys.argv[0]).name + " " + shlex.join(sys.argv[1:])


# ---------------------------------------------------------------------------
# Config resolution from optimization/results.tsv (--from-config auto|GCn)
# ---------------------------------------------------------------------------

def _load_optimization_results():
    if not OPTIMIZATION_TSV.exists():
        return None
    with open(OPTIMIZATION_TSV) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["N"] = int(r["N"])
        r["pool_size"] = int(r["pool_size"])
        r["tm_threshold"] = float(r["tm_threshold"])
        for m in METRIC_COLS + ["D_composite"]:
            v = r.get(m)
            r[m] = float(v) if v not in ("", "None", None) else None
    return rows


def _load_baseline_stages():
    """Loads baseline_stages.tsv's "redundancy" stage rows: whole-pool
    diversity + pool size for all 27 configs, computed directly on each
    config's pool with no sampling. This -- not results.tsv -- is what Step 1
    needs, and it stays complete for all 27 configs even when results.tsv
    itself only holds the winning config's grid (01_optimize_benchmark.py's
    default, non---full-diagnostic mode only grid-searches the Step-1
    winner; see its module docstring)."""
    if not BASELINE_TSV.exists():
        return None
    with open(BASELINE_TSV) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rows = [r for r in rows if r["stage"] == "redundancy"]
    for r in rows:
        r["n"] = int(r["n"])
        r["tm_threshold"] = float(r["tm_threshold"])
    return rows


def _select_winning_config(baseline_rows, n_min_pool=N_MIN_POOL):
    """Step 1: fix the redundancy config directly from filter severity + pool
    size, NOT from D_composite. Both the whole-pool pairwise-average
    diversity metrics (M5/M6/M7/M10) and local-redundancy metrics (nearest-
    neighbor distance, percentile distances) were shown to just re-encode
    pool size across the 27 configs (|Spearman| with pool_size up to 0.96),
    not an independent diversity signal -- see reports/methods.tex. Selects
    the most stringent (receptor, peptide, TM-score) combination whose pool
    retains at least n_min_pool entries, ties broken by the larger margin
    above that floor (i.e. the larger pool).

    Operates on baseline_stages.tsv's "redundancy" rows (see
    _load_baseline_stages), NOT on results.tsv -- the grid file may only
    contain the winning config's rows, so it can't be used to independently
    re-derive which config should win among all 27."""
    qualified = [r for r in baseline_rows if r["n"] >= n_min_pool]
    if not qualified:
        raise ValueError(f"No config retains >= {n_min_pool} entries in {BASELINE_TSV}")

    def severity(r):
        return (STRICTNESS_RANK[r["rec_preset"]] + STRICTNESS_RANK[r["pep_preset"]]
                + TM_SEVERITY_RANK[r["tm_threshold"]])

    qualified.sort(key=lambda r: (severity(r), -r["n"]))
    return qualified[0]["config_id"]


def _local_composite(rows):
    """Recompute D_composite using min-max bounds derived ONLY from `rows`
    (e.g. a single config's own grid rows), so Step 2/3 compare candidates
    without leakage from the global results.tsv bounds (which span all 27
    configs and would let other configs' extremes distort this comparison)."""
    bounds = {}
    for m in METRIC_COLS:
        vals = [r[m] for r in rows if r.get(m) is not None]
        if not vals:
            continue
        lo, hi = min(vals), max(vals)
        bounds[m] = (lo, hi if hi != lo else lo + 1.0)

    out = []
    for r in rows:
        norms = []
        for m in METRIC_COLS:
            if m not in bounds or r.get(m) is None:
                continue
            lo, hi = bounds[m]
            norm = (r[m] - lo) / (hi - lo)
            norms.append((1 - norm) if m in LOWER_IS_BETTER else norm)
        out.append({**r, "D_local": sum(norms) / len(norms) if norms else None})
    return out


def _best_run_for_config(rows, config_id, tolerance=N_TOLERANCE):
    """Steps 2-3, within one fixed config (from Step 1):

    Step 2 picks the (method, dims) combination with the highest mean
    D_local (bounds local to this config) across its own N curve.

    Step 3 picks the largest N on that winning curve whose D_local is within
    `tolerance` of that curve's own peak -- not the raw maximum, which
    collapses to the smallest N tested due to M1_mnnd's mechanical N-
    dependence (see reports/methods.tex)."""
    sub = [r for r in rows if r["config_id"] == config_id]
    if not sub:
        raise ValueError(f"No results found for {config_id} in {OPTIMIZATION_TSV}")

    sub = [r for r in _local_composite(sub) if r["D_local"] is not None]

    by_md = defaultdict(list)
    for r in sub:
        by_md[(r["method"], r["dims_id"])].append(r)

    means = {k: sum(r["D_local"] for r in v) / len(v) for k, v in by_md.items()}
    best_method, best_dims = max(means, key=means.get)

    curve = by_md[(best_method, best_dims)]
    d_max = max(r["D_local"] for r in curve)
    qualifying = [r for r in curve if r["D_local"] >= (1 - tolerance) * d_max]
    target = max(qualifying, key=lambda r: r["N"])

    return {
        "rec_preset": target["rec_preset"], "pep_preset": target["pep_preset"],
        "tm_threshold": target["tm_threshold"], "pool_size": target["pool_size"],
        "sampling_method": target["method"], "dims": target["dims"].split("+"),
        "sample_n": target["N"], "d_composite": target["D_local"],
    }


def _resolve_profile(args):
    """Returns (profile dict, profile_name, source) for either --profile or
    --from-config. Does not resolve --version — main() does that, since dry
    runs don't need a version directory."""
    if args.profile:
        profile = dict(PROFILES[args.profile])
        if args.seed is not None:
            profile["seed"] = args.seed
        return profile, args.profile, "profile"

    opt_rows = _load_optimization_results()
    if opt_rows is None:
        print(f"ERROR: {OPTIMIZATION_TSV} not found.")
        print("  Run 01_optimize_benchmark.py --full first, or use --profile "
              "for a frozen release that doesn't need it.")
        sys.exit(1)

    if args.from_config == "auto":
        baseline_rows = _load_baseline_stages()
        if baseline_rows is None:
            print(f"ERROR: {BASELINE_TSV} not found (needed for Step 1).")
            print("  Run 01_optimize_benchmark.py --full first.")
            sys.exit(1)
        config_id = _select_winning_config(baseline_rows)
        source = "auto"
    else:
        config_id = args.from_config
        source = "from-config"

    best = _best_run_for_config(opt_rows, config_id)
    print(f"  Resolved {config_id} from results.tsv: "
          f"rec={best['rec_preset']} pep={best['pep_preset']} tm={best['tm_threshold']} "
          f"| {best['sampling_method']} {'+'.join(best['dims'])} N={best['sample_n']} "
          f"(D_local={best['d_composite']:.3f}, pool={best['pool_size']})")

    profile = {
        "version":         None,
        "seq_id_receptor": best["rec_preset"],
        "seq_id_peptide":  best["pep_preset"],
        "tm_score":        best["tm_threshold"],
        "sampling_method": best["sampling_method"],
        "dims":            best["dims"],
        "sample_n":        best["sample_n"],
        "seed":            args.seed if args.seed is not None else 42,
    }
    return profile, config_id, source


def _filter_args(profile):
    """Build a Namespace with all fields expected by genSample filter modules."""
    return argparse.Namespace(
        single_chain=SINGLE_CHAIN,
        min_receptor=None,   max_receptor=None,
        min_peptide=None,    max_peptide=None,
        min_dg=None,         max_dg=None,
        min_year=None,       max_year=None,
        min_resolution=None, max_resolution=None,
        binding_type="both",
        max_added=None,      max_receptor_added=None,
        exclude_nmr=False,   exclude=None,
        seq_id_receptor=profile["seq_id_receptor"], seq_id_receptor_threshold=None,
        seq_id_peptide=profile["seq_id_peptide"],   seq_id_peptide_threshold=None,
        tm_score_receptor=profile["tm_score"],
        dims=profile["dims"],
    )


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def run_pipeline(profile_name, profile, source, version_dir, dry_run=False):
    a = _filter_args(profile)

    all_rows  = load_tsv(SUMMARY_TSV)
    full_cols = list(all_rows[0]) if all_rows else None
    all_rows  = enrich_rows(all_rows)
    print(f"  Loaded {len(all_rows)} entries from {SUMMARY_TSV.name}")

    pool_quality    = apply_quality_filters(all_rows, a, _filtered)
    pool_seq, _     = apply_redundancy_filters(pool_quality, a)
    pool_redundancy = apply_structural_redundancy_filters(pool_seq, a, _fixed_cif)

    if dry_run:
        dg  = [r["_dg"]      for r in pool_redundancy]
        pep = [r["_pep_len"] for r in pool_redundancy]
        rec = [r["_rec_len"] for r in pool_redundancy]
        print(f"\n  Pool after all filters: {len(pool_redundancy)} entries")
        if dg:
            print(f"  dG:         {min(dg):.2f} – {max(dg):.2f} kcal/mol")
            print(f"  Pep length: {min(pep)} – {max(pep)} residues")
            print(f"  Rec length: {min(rec)} – {max(rec)} residues")
        return None

    sampling_method = profile["sampling_method"]
    dims            = profile["dims"]
    sample_n        = profile["sample_n"]
    seed            = profile["seed"]

    selected = run_sampling(pool_redundancy, sample_n,
                            method=sampling_method, dims=dims, seed=seed)
    selected.sort(key=lambda r: r["_dg"])
    print(f"  Sampled {len(selected)} of {len(pool_redundancy)} entries "
          f"({method_label(sampling_method, dims)})")

    version_dir.mkdir(parents=True, exist_ok=True)

    write_sample(selected, version_dir, columns=full_cols if FULL_METADATA else None)

    save_tsv(pool_quality,    version_dir / "pool_quality.tsv",    POOL_COLS)
    save_tsv(pool_redundancy, version_dir / "pool_redundancy.tsv", POOL_COLS)
    print(f"  Written pool_quality.tsv ({len(pool_quality)}) and "
          f"pool_redundancy.tsv ({len(pool_redundancy)})")

    config = {
        "profile":         profile_name,
        "source":          source,   # "profile" (frozen) | "auto" | "from-config" (explicit GCn)
        "command":         _invoked_command(),
        "sample_size":     len(selected),
        "pool_size":       len(pool_redundancy),
        "sampling_method": method_label(sampling_method, dims),
        "dims":            dims,
        "seed":            seed,
        "filters": {
            "single_chain":      SINGLE_CHAIN,
            "seq_id_receptor":   profile["seq_id_receptor"],
            "seq_id_peptide":    profile["seq_id_peptide"],
            "tm_score_receptor": profile["tm_score"],
        },
        "sources": {
            "summary_tsv":   str(SUMMARY_TSV),
            "fixed_cif_dir": str(_fixed_cif),
        },
        "generated": date.today().isoformat(),
    }
    cfg_out = version_dir / "run_config.json"
    with open(cfg_out, "w", encoding="ascii") as fh:
        json.dump(config, fh, indent=2)
    print(f"  Written run_config.json")

    link_cifs(selected, version_dir, _fixed_cif, copy=True, fmt="cif")

    cif_dir = version_dir / "cif"
    if cif_dir.exists() and any(cif_dir.glob("*.cif")):
        _convert_cifs_to_pdb(cif_dir, version_dir / "pdb")

    if EXPORT_FASTA:
        seq_dir = version_dir / "sequences"
        seq_dir.mkdir(exist_ok=True)
        write_fasta(selected, seq_dir)

    return config, len(pool_quality), len(pool_redundancy), len(selected)


def _convert_cifs_to_pdb(cif_dir, pdb_dir):
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    pdb_dir.mkdir(parents=True, exist_ok=True)
    cif_files = list(cif_dir.glob("*.cif"))
    done, failed = 0, []
    for cif_path in cif_files:
        dst = pdb_dir / f"{cif_path.stem}.pdb"
        try:
            fixer = PDBFixer(filename=str(cif_path))
            with open(dst, "w", encoding="ascii") as fh:
                PDBFile.writeFile(fixer.topology, fixer.positions, fh)
            done += 1
        except Exception as exc:
            print(f"  Warning: CIF→PDB failed for {cif_path.name}: {exc}")
            failed.append(cif_path.stem)
    label = f" ({len(failed)} failed)" if failed else ""
    print(f"  Converted {done}/{len(cif_files)} CIFs → PDB{label}")
    return done


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def generate_plots(version_dir):
    plots_dst       = version_dir / "plots"
    plots_dst.mkdir(parents=True, exist_ok=True)
    _scrape_scripts = _repo_root / "scrape" / "scripts"
    _gen_scripts    = _repo_root / "genSample" / "scripts"

    print("  Generating summary_distributions.png ...")
    subprocess.run(
        [sys.executable, str(_scrape_scripts / "05_distributions.py"),
         "--summary-tsv", str(version_dir / "sample.tsv"),
         "--output-dir",  str(plots_dst)],
        check=True,
    )

    for dim_y, fname in [("ss_class",   "run_comparison_dg_ss_class.png"),
                         ("func_class", "run_comparison_dg_func_class.png")]:
        print(f"  Generating {fname} ...")
        subprocess.run(
            [sys.executable, str(_gen_scripts / "02_plot_sampleDist.py"),
             "--run-dir", str(version_dir),
             "--dims",    "dg", dim_y,
             "--output",  str(plots_dst / fname)],
            check=True,
        )

    print(f"  {len(list(plots_dst.glob('*.png')))} plots written → {plots_dst}")


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------

def write_readme(version_dir, n_quality, n_redundancy, n_sample, n_pdb):
    from modules.sequence_redundancy import PRESETS

    with open(version_dir / "run_config.json") as f:
        config = json.load(f)

    n_cif = len(list((version_dir / "cif").glob("*.cif"))) if (version_dir / "cif").exists() else 0

    filt = config["filters"]
    rec_pct = PRESETS[filt["seq_id_receptor"]]["receptor"] * 100
    pep_pct = PRESETS[filt["seq_id_peptide"]]["peptide"] * 100
    tm      = filt["tm_score_receptor"]

    lines = [
        f"PepXPro Benchmark — {version_dir.name}  (profile {config['profile']})",
        "=" * 42,
        f"Generated : {date.today().isoformat()}",
        "",
        "CONTENTS",
        "--------",
        f"sample.tsv                  {n_sample} entries, full metadata (incl. external dataset fields)",
        f"cif/                        {n_cif} fixed CIFs (AMBER14-compatible)",
        f"pdb/                        {n_pdb} PDB conversions",
        "sequences/",
        "  peptide_sequences.fasta",
        "  receptor_sequences.fasta",
        "plots/",
        "  summary_distributions.png",
        "  run_comparison_dg_ss_class.png",
        "  run_comparison_dg_func_class.png",
        "pool_quality.tsv            Entries after quality filters",
        "pool_redundancy.tsv         Entries after redundancy filters",
        "run_config.json             Full provenance",
        "",
        "PIPELINE",
        "--------",
        f"  Input:           {config['sources']['summary_tsv']}",
        f"  Quality pool:    {n_quality} entries  (--single-chain)",
        f"  Redundancy pool: {n_redundancy} entries  "
        f"(seq receptor {rec_pct:.0f}% [{filt['seq_id_receptor']}], "
        f"peptide {pep_pct:.0f}% [{filt['seq_id_peptide']}], TM-score {tm})",
        f"  Final sample:    {n_sample} entries  ({config['sampling_method']}, "
        f"N={config['sample_size']})",
        "",
        "SELECTION CRITERIA",
        "------------------",
        "Single receptor chain.",
        f"Sequence redundancy: receptor {rec_pct:.0f}% CD-HIT ({filt['seq_id_receptor']}), "
        f"peptide {pep_pct:.0f}% Biopython pairwise ({filt['seq_id_peptide']}).",
        f"Structural redundancy: receptor TM-score < {tm} (different fold).",
        f"Sampling: {config['sampling_method']} in {' x '.join(config['dims'])} space, "
        f"N chosen as the largest sample size within {N_TOLERANCE*100:.0f}% of the peak "
        "D_composite on the winning method/feature-space curve for this config "
        "(pepxpro_benchmark/optimization/results.tsv; see reports/methods.tex).",
        f"Source: PDBbind 2020 + external dataset (summary_merged.tsv), seed={config['seed']}.",
        "",
        "PARAMETERS",
        "----------",
        json.dumps(config, indent=2),
    ]

    out = version_dir / "README.txt"
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"  Written README.txt")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a versioned pepxpro_benchmark release.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--profile", choices=sorted(PROFILES),
                        help="Frozen, citation-safe benchmark profile (e.g. GC2 -> v1). "
                             "Reproduces an exact past release regardless of the current "
                             "optimization/results.tsv contents.")
    group.add_argument("--from-config", metavar="auto|GCn",
                        help="Derive redundancy/sampling parameters live from "
                             "optimization/results.tsv. 'auto' reproduces the Figure 5 "
                             "winner-selection rule; 'GCn' (e.g. GC11) targets that grid "
                             "config directly, regardless of whether it 'wins'.")
    parser.add_argument("--version", metavar="VX",
                        help="Version label, e.g. v3 (default for --profile: the profile's "
                             "fixed version, e.g. GC2 -> v1; required for --from-config).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Override the sampling seed (default: the profile's own seed, "
                             "or 42 for --from-config).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show filtered pool size without generating files.")
    args = parser.parse_args()

    profile, profile_name, source = _resolve_profile(args)

    if args.dry_run:
        print(f"PepXPro Benchmark — dry run ({profile_name}, source={source})")
        run_pipeline(profile_name, profile, source, None, dry_run=True)
        return

    version_label = args.version or profile.get("version")
    if not version_label:
        print("ERROR: --version is required when using --from-config.")
        sys.exit(1)
    version_dir = _bench_root / version_label
    if version_dir.exists():
        print(f"Error: {version_dir} already exists.")
        sys.exit(1)

    print(f"PepXPro Benchmark Builder")
    print(f"  Config  : {profile_name}  (source={source})")
    print(f"  Version : {version_dir.name}")
    print(f"  Output  : {version_dir}")
    print(f"  Command : {_invoked_command()}")
    print()

    result = run_pipeline(profile_name, profile, source, version_dir)
    if result is None:
        sys.exit(1)
    config, n_quality, n_redundancy, n_sample = result

    print(f"\n  Generating plots ...")
    generate_plots(version_dir)

    n_pdb = len(list((version_dir / "pdb").glob("*.pdb"))) if (version_dir / "pdb").exists() else 0
    write_readme(version_dir, n_quality, n_redundancy, n_sample, n_pdb)

    n_cif = len(list((version_dir / "cif").glob("*.cif"))) if (version_dir / "cif").exists() else 0
    print(f"\n  Done → {version_dir}")
    print(f"  {n_sample} entries | {n_cif} CIFs | {n_pdb} PDBs | sequences/ | 3 plots | README.txt")


if __name__ == "__main__":
    main()
