"""
01_make_benchmark.py — Generate a versioned pepxpro_benchmark release.

Uses genSample modules directly.  No run_XX is created in genSample/outputs/.

Benchmark versions are generated from named PROFILES, each corresponding to a
grid configuration identified in pepxpro_benchmark/optimization/results.tsv
(see reports/Draft/figure1_distributions.ipynb, Figures 5-6). Every profile
uses its own elbow N (diminishing-returns benchmark size) and the method+dims
combo that scored best for that grid config.

Intended workflow:
  1. Explore and tune parameters with genSample/01_genSample.py
     (--dry-run, --sample, different --dims, etc.)
  2. Run pepxpro_benchmark/scripts/02_optimize_benchmark.py --full to grid-search
     redundancy configs and pick winners (see optimization/results.tsv).
  3. Run this script once per winning profile for a reproducible release.

Usage:
  python 01_make_benchmark.py --profile GC5            # -> v1
  python 01_make_benchmark.py --profile GC11           # -> v2
  python 01_make_benchmark.py --profile GC5 --dry-run  # show pool size only
"""

import argparse
import json
import shlex
import subprocess
import sys
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
    "GC5": {
        "version":         "v1",
        "seq_id_receptor": "strict",   # rec_preset — receptor 40% CD-HIT
        "seq_id_peptide":  "medium",   # pep_preset — peptide 80% Biopython pairwise
        "tm_score":        0.7,
        "sampling_method": "kennard-stone",
        "dims":            ["dg", "pep_len", "rec_len", "n_contacts"],  # D2
        "sample_n":        75,         # elbow N for GC5
        "seed":            42,
    },
    "GC11": {
        "version":         "v2",
        "seq_id_receptor": "medium",   # rec_preset — receptor 50% CD-HIT
        "seq_id_peptide":  "strict",   # pep_preset — peptide 70% Biopython pairwise
        "tm_score":        0.7,
        "sampling_method": "kennard-stone",
        "dims":            ["dg", "pep_len", "rec_len", "n_contacts"],  # D2
        "sample_n":        100,        # elbow N for GC11
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


def _invoked_command():
    """The exact command line that produced this benchmark release, for reproducibility."""
    return "python " + Path(sys.argv[0]).name + " " + shlex.join(sys.argv[1:])


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

def run_pipeline(profile_name, version_dir, dry_run=False):
    profile = PROFILES[profile_name]
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
    parser.add_argument("--profile", choices=sorted(PROFILES), required=True,
                        help="Named benchmark profile (grid config) to build.")
    parser.add_argument("--version", metavar="VX",
                        help="Version label override, e.g. v3 "
                             "(default: the profile's fixed version, e.g. GC5 -> v1).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show filtered pool size without generating files.")
    args = parser.parse_args()

    if args.dry_run:
        print(f"PepXPro Benchmark — dry run (profile {args.profile})")
        run_pipeline(args.profile, None, dry_run=True)
        return

    version_label = args.version or PROFILES[args.profile]["version"]
    version_dir = _bench_root / version_label
    if version_dir.exists():
        print(f"Error: {version_dir} already exists.")
        sys.exit(1)

    print(f"PepXPro Benchmark Builder")
    print(f"  Profile : {args.profile}")
    print(f"  Version : {version_dir.name}")
    print(f"  Output  : {version_dir}")
    print(f"  Command : {_invoked_command()}")
    print()

    result = run_pipeline(args.profile, version_dir)
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
