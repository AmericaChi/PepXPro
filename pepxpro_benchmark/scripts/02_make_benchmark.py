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
                          reproducing the Figure 5 winner-selection rule
                          (elbow_N/pool_size < 0.6, ties broken toward
                          stricter presets).
  --from-config GC11     Target a specific grid config directly, regardless
                          of whether it "wins" the comparison — e.g. to try
                          an alternative to the adopted default.

  In both --from-config modes, the method/dims/N for the chosen config are
  picked automatically: the (method, dims) combo with the highest
  D_composite at that config's own elbow N (falls back to the single
  best-scoring row if the config has fewer than 3 tested N values).

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
# that scored best (D_composite) for that config at its elbow N.
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
        "sample_n":        75,         # elbow N for GC2
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
STRICTNESS_RANK  = {"strict": 0, "medium": 1, "lenient": 2}


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
        d = r.get("D_composite")
        r["D_composite"] = float(d) if d not in ("", "None", None) else None
    return rows


def _elbow(ns, ds):
    """Point of maximum discrete curvature on (N, D_composite); needs >=3 points.

    Same rule used in reports/Draft/figure1_distributions.ipynb (Figure 5) —
    kept in sync manually since this script has no dependency on the notebook.
    """
    if len(ns) < 3:
        return None
    nr = max(ns) - min(ns); dr = max(ds) - min(ds)
    if nr == 0 or dr == 0:
        return None
    xn = [(n - min(ns)) / nr for n in ns]
    yn = [(d - min(ds)) / dr for d in ds]
    return max(
        ((abs(yn[i + 1] - 2 * yn[i] + yn[i - 1]), ns[i]) for i in range(1, len(ns) - 1)),
        key=lambda x: x[0]
    )[1]


def _best_run_for_config(rows, config_id):
    """Resolve one grid config to concrete build parameters: the (method, dims)
    combo with the highest D_composite at that config's own elbow N (envelope
    of best D_composite per N across all methods/dims). Falls back to the
    single globally best-scoring row if fewer than 3 N values were tested."""
    sub = [r for r in rows if r["config_id"] == config_id and r["D_composite"] is not None]
    if not sub:
        raise ValueError(f"No results found for {config_id} in {OPTIMIZATION_TSV}")

    by_n = {}
    for r in sub:
        if r["N"] not in by_n or r["D_composite"] > by_n[r["N"]]["D_composite"]:
            by_n[r["N"]] = r
    ns = sorted(by_n)
    elbow_n = _elbow(ns, [by_n[n]["D_composite"] for n in ns])
    target = by_n[elbow_n] if elbow_n is not None else max(sub, key=lambda r: r["D_composite"])

    return {
        "rec_preset": target["rec_preset"], "pep_preset": target["pep_preset"],
        "tm_threshold": target["tm_threshold"], "pool_size": target["pool_size"],
        "sampling_method": target["method"], "dims": target["dims"].split("+"),
        "sample_n": target["N"], "d_composite": target["D_composite"],
        "elbow_n": elbow_n,
    }


def _select_winning_config(rows):
    """Reproduce the Figure 5 winner-selection rule: among grid configs with a
    computable elbow that lands before 60% of the pool is exhausted, pick the
    highest mean D_composite, breaking ties toward stricter presets."""
    by_cfg = defaultdict(list)
    for r in rows:
        if r["D_composite"] is not None:
            by_cfg[r["config_id"]].append(r)

    summary = []
    for cfg_id, sub in by_cfg.items():
        by_n = {}
        for r in sub:
            if r["N"] not in by_n or r["D_composite"] > by_n[r["N"]]["D_composite"]:
                by_n[r["N"]] = r
        ns = sorted(by_n)
        elbow_n = _elbow(ns, [by_n[n]["D_composite"] for n in ns])
        pool = sub[0]["pool_size"]
        summary.append({
            "config_id": cfg_id,
            "rec_preset": sub[0]["rec_preset"], "pep_preset": sub[0]["pep_preset"],
            "mean_d": sum(r["D_composite"] for r in sub) / len(sub),
            "elbow_ratio": (elbow_n / pool) if elbow_n is not None else None,
        })

    qualified = [s for s in summary if s["elbow_ratio"] is not None and s["elbow_ratio"] < 0.6]
    if not qualified:
        raise ValueError("No grid config qualifies (elbow_ratio < 0.6) in results.tsv.")

    qualified.sort(key=lambda s: (-s["mean_d"], STRICTNESS_RANK[s["rec_preset"]],
                                   STRICTNESS_RANK[s["pep_preset"]]))
    return qualified[0]["config_id"]


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
        config_id = _select_winning_config(opt_rows)
        source = "auto"
    else:
        config_id = args.from_config
        source = "from-config"

    best = _best_run_for_config(opt_rows, config_id)
    elbow_note = f"elbow N={best['elbow_n']}" if best["elbow_n"] is not None else "no elbow — used best-scoring row"
    print(f"  Resolved {config_id} from results.tsv: "
          f"rec={best['rec_preset']} pep={best['pep_preset']} tm={best['tm_threshold']} "
          f"| {best['sampling_method']} {'+'.join(best['dims'])} N={best['sample_n']} "
          f"(D_composite={best['d_composite']:.3f}, {elbow_note}, pool={best['pool_size']})")

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
        f"elbow N={config['sample_size']})",
        "",
        "SELECTION CRITERIA",
        "------------------",
        "Single receptor chain.",
        f"Sequence redundancy: receptor {rec_pct:.0f}% CD-HIT ({filt['seq_id_receptor']}), "
        f"peptide {pep_pct:.0f}% Biopython pairwise ({filt['seq_id_peptide']}).",
        f"Structural redundancy: receptor TM-score < {tm} (different fold).",
        f"Sampling: {config['sampling_method']} in {' x '.join(config['dims'])} space, "
        f"N chosen at the D_composite-vs-N elbow for this grid config "
        "(pepxpro_benchmark/optimization/results.tsv).",
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
