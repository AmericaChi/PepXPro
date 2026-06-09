"""
01_make_benchmark.py — Generate a versioned pepxpro_benchmark release.

Uses gen_sample modules directly.  No run_XX is created in gen_sample/outputs/.

Intended workflow:
  1. Explore and tune parameters with gen_sample/01_gen_sample.py
     (--dry-run, --sample, different --dims, etc.)
  2. Once parameters are final, run this script for a reproducible release.

Usage:
  python 01_make_benchmark.py              # auto-increment version (v1, v2, ...)
  python 01_make_benchmark.py --version v2 # explicit version label
  python 01_make_benchmark.py --dry-run    # show pool size without writing files
"""

import argparse
import json
import subprocess
import sys
from datetime import date
from pathlib import Path

_here       = Path(__file__).resolve().parent
_bench_root = _here.parent
_repo_root  = _bench_root.parent
_scrape_out = _repo_root / "scrape_pdb_bind" / "outputs"
_fixed_cif  = _repo_root / "scrape_pdb_bind" / "fixed_cif"
_filtered   = _scrape_out / "filtered_systems.tsv"

sys.path.insert(0, str(_repo_root / "gen_sample"))
from modules.utils import load_tsv, enrich_rows, save_tsv
from modules.filters import apply_quality_filters
from modules.sequence_redundancy import apply_redundancy_filters
from modules.structural_redundancy import apply_structural_redundancy_filters
from modules.sampling import run_sampling, method_label
from modules.export import write_sample, link_cifs, write_fasta


# ── Benchmark parameters ─────────────────────────────────────────────────────
# Edit these to define a new benchmark version.

SUMMARY_TSV     = _scrape_out / "summary_merged.tsv"

SAMPLE_N        = 89
SEED            = 42
SAMPLING_METHOD = "kennard-stone"
DIMS            = ["dg", "ss_class", "func_class"]

SINGLE_CHAIN    = True
SEQ_ID_RECEPTOR = "strict"   # 40% CD-HIT
SEQ_ID_PEPTIDE  = "strict"   # 70% Biopython pairwise
TM_SCORE        = 0.5

FULL_METADATA   = True
EXPORT_FASTA    = True

POOL_COLS = ["pdb_id", "receptor", "peptide", "dG_kcal_mol", "binding_type",
             "resolution", "peptide_ss_class", "receptor_func_class"]
# ─────────────────────────────────────────────────────────────────────────────


def _filter_args():
    """Build a Namespace with all fields expected by gen_sample filter modules."""
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
        seq_id_receptor=SEQ_ID_RECEPTOR, seq_id_receptor_threshold=None,
        seq_id_peptide=SEQ_ID_PEPTIDE,   seq_id_peptide_threshold=None,
        tm_score_receptor=TM_SCORE,
        dims=DIMS,
    )


def next_version_dir():
    for i in range(1, 1000):
        d = _bench_root / f"v{i}"
        if not d.exists():
            return d
    raise RuntimeError("No available version directory")


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def run_pipeline(version_dir, dry_run=False):
    a = _filter_args()

    all_rows  = load_tsv(SUMMARY_TSV)
    full_cols = list(all_rows[0]) if all_rows else None
    all_rows  = enrich_rows(all_rows)
    print(f"  Loaded {len(all_rows)} entries from {SUMMARY_TSV.name}")

    pool_quality    = apply_quality_filters(all_rows, a, _filtered)
    pool_seq        = apply_redundancy_filters(pool_quality, a)
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

    selected = run_sampling(pool_redundancy, SAMPLE_N,
                            method=SAMPLING_METHOD, dims=DIMS, seed=SEED)
    selected.sort(key=lambda r: r["_dg"])
    print(f"  Sampled {len(selected)} of {len(pool_redundancy)} entries "
          f"({method_label(SAMPLING_METHOD, DIMS)})")

    version_dir.mkdir(parents=True, exist_ok=True)

    write_sample(selected, version_dir, columns=full_cols if FULL_METADATA else None)

    save_tsv(pool_quality,    version_dir / "pool_quality.tsv",    POOL_COLS)
    save_tsv(pool_redundancy, version_dir / "pool_redundancy.tsv", POOL_COLS)
    print(f"  Written pool_quality.tsv ({len(pool_quality)}) and "
          f"pool_redundancy.tsv ({len(pool_redundancy)})")

    config = {
        "sample_size":     len(selected),
        "pool_size":       len(pool_redundancy),
        "sampling_method": method_label(SAMPLING_METHOD, DIMS),
        "dims":            DIMS,
        "seed":            SEED,
        "filters": {
            "single_chain":      SINGLE_CHAIN,
            "seq_id_receptor":   SEQ_ID_RECEPTOR,
            "seq_id_peptide":    SEQ_ID_PEPTIDE,
            "tm_score_receptor": TM_SCORE,
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
    _scrape_scripts = _repo_root / "scrape_pdb_bind" / "scripts"
    _gen_scripts    = _repo_root / "gen_sample" / "scripts"

    print("  Generating summary_distributions.png ...")
    subprocess.run(
        [sys.executable, str(_scrape_scripts / "04_distributions.py"),
         "--summary-tsv", str(version_dir / "sample.tsv"),
         "--output-dir",  str(plots_dst)],
        check=True,
    )

    for dim_y, fname in [("ss_class",   "run_comparison_dg_ss_class.png"),
                         ("func_class", "run_comparison_dg_func_class.png")]:
        print(f"  Generating {fname} ...")
        subprocess.run(
            [sys.executable, str(_gen_scripts / "03_compare_runs.py"),
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
    with open(version_dir / "run_config.json") as f:
        config = json.load(f)

    n_cif = len(list((version_dir / "cif").glob("*.cif"))) if (version_dir / "cif").exists() else 0

    lines = [
        f"PepXPro Benchmark — {version_dir.name}",
        "=" * 42,
        f"Generated : {date.today().isoformat()}",
        "",
        "CONTENTS",
        "--------",
        f"sample.tsv                  {n_sample} entries, full metadata",
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
        f"(seq receptor 40%, peptide 70%, TM-score 0.5)",
        f"  Final sample:    {n_sample} entries  ({config['sampling_method']})",
        "",
        "SELECTION CRITERIA",
        "------------------",
        "Single receptor chain.",
        "Sequence redundancy: receptor 40% CD-HIT, peptide 70% Biopython pairwise.",
        "Structural redundancy: receptor TM-score < 0.5 (different fold).",
        "Sampling: Kennard-Stone in dG x peptide-SS-class x receptor-function space.",
        f"Source: PDBbind 2020 + PEPBI (summary_merged.tsv), seed={SEED}.",
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
    parser.add_argument("--version", metavar="VX",
                        help="Version label, e.g. v2 (default: auto-increment).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show filtered pool size without generating files.")
    args = parser.parse_args()

    if args.dry_run:
        print("PepXPro Benchmark — dry run")
        run_pipeline(None, dry_run=True)
        return

    if args.version:
        version_dir = _bench_root / args.version
        if version_dir.exists():
            print(f"Error: {version_dir} already exists.")
            sys.exit(1)
    else:
        version_dir = next_version_dir()

    print(f"PepXPro Benchmark Builder")
    print(f"  Version : {version_dir.name}")
    print(f"  Output  : {version_dir}")
    print()

    result = run_pipeline(version_dir)
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
