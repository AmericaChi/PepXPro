"""
01_genSample.py — Generate a peptide-protein benchmark sample from PDBbind outputs.

Reads summary.tsv (and optionally filtered_systems.tsv) produced by scrape,
applies filters, samples entries, and writes outputs to genSample/outputs/run_XX/.
A log is written to genSample/logs/run_XX.log for every non-dry-run execution.

Usage examples:
  python 01_genSample.py --sample 20
  python 01_genSample.py --sample 20 --single-chain --max-peptide 15
  python 01_genSample.py --sample 20 --sampling kennard-stone --seed 42
  python 01_genSample.py --no-sample
  python 01_genSample.py --dry-run
  python 01_genSample.py --pick 1abc --protonate --ph 7.4
"""

import argparse
import sys
from pathlib import Path

# Allow running from scripts/ directly or from the genSample root.
_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, enrich_rows
from modules.filters import apply_quality_filters
from modules.sequence_redundancy import apply_redundancy_filters
from modules.structural_redundancy import apply_structural_redundancy_filters
from modules.sampling import run_sampling, SAMPLING_METHODS, method_label, AVAILABLE_DIMS
from modules.export import next_run_dir, write_sample, save_run_config, link_cifs, protonate_cifs, write_fasta
from modules.utils import save_tsv

POOL_COLS = ["pdb_id", "receptor", "peptide", "dG_kcal_mol", "binding_type",
             "resolution", "peptide_ss_class", "receptor_func_class"]

# ---------------------------------------------------------------------------
# Paths (relative to scrape, which lives next to genSample)
# ---------------------------------------------------------------------------

_SCRAPE_DIR    = _root.parent / "scrape"
SUMMARY_TSV    = _SCRAPE_DIR / "outputs" / "summary.tsv"
FILTERED_TSV   = _SCRAPE_DIR / "outputs" / "filtered_systems.tsv"
FIXED_CIF_DIR  = _SCRAPE_DIR / "fixed_cif"
OUTPUTS_DIR    = _root / "outputs"
LOGS_DIR       = _root / "logs"


class _Tee:
    """Write to both stdout and a log file simultaneously."""
    def __init__(self, log_path):
        self._file = open(log_path, "w", encoding="utf-8")
        self._stdout = sys.stdout

    def write(self, msg):
        self._stdout.write(msg)
        self._file.write(msg)

    def flush(self):
        self._stdout.flush()
        self._file.flush()

    def close(self):
        self._file.close()
        sys.stdout = self._stdout


# ---------------------------------------------------------------------------
# Pool summary (--dry-run)
# ---------------------------------------------------------------------------

def print_pool_summary(rows):
    dg_vals  = [r["_dg"]      for r in rows]
    pep_vals = [r["_pep_len"] for r in rows]
    rec_vals = [r["_rec_len"] for r in rows]

    bt_counts = {}
    for r in rows:
        bt = r.get("binding_type", "?")
        bt_counts[bt] = bt_counts.get(bt, 0) + 1

    print(f"\n  Pool: {len(rows)} entries")
    print(f"  {'':20}  {'min':>8} {'max':>8} {'avg':>8}")
    print(f"  {'dG (kcal/mol)':<20}  {min(dg_vals):>8.2f} {max(dg_vals):>8.2f} "
          f"{sum(dg_vals)/len(dg_vals):>8.2f}")
    print(f"  {'peptide length':<20}  {min(pep_vals):>8} {max(pep_vals):>8} "
          f"{sum(pep_vals)/len(pep_vals):>8.1f}")
    print(f"  {'receptor length':<20}  {min(rec_vals):>8} {max(rec_vals):>8} "
          f"{sum(rec_vals)/len(rec_vals):>8.1f}")
    bt_str = "  binding types: " + ", ".join(f"{k}={v}" for k, v in sorted(bt_counts.items()))
    print(bt_str)
    print()


# ---------------------------------------------------------------------------
# Sample summary (printed after selection)
# ---------------------------------------------------------------------------

def print_sample_summary(rows):
    dg_vals  = [r["_dg"]      for r in rows]
    pep_vals = [r["_pep_len"] for r in rows]
    rec_vals = [r["_rec_len"] for r in rows]

    print(f"\n  {'PDB':<8} {'dG (kcal/mol)':>15} {'type':<6}  receptor / peptide")
    print("  " + "-" * 62)
    for r in rows:
        print(f"  {r['pdb_id']:<8} {r['dG_kcal_mol']:>15} {r['binding_type']:<6} "
              f" {r['receptor']}  /  {r['peptide']}")
    print()
    print(f"  {'':20}  {'min':>8} {'max':>8} {'avg':>8}")
    print(f"  {'dG (kcal/mol)':<20}  {min(dg_vals):>8.2f} {max(dg_vals):>8.2f} "
          f"{sum(dg_vals)/len(dg_vals):>8.2f}")
    print(f"  {'peptide length':<20}  {min(pep_vals):>8} {max(pep_vals):>8} "
          f"{sum(pep_vals)/len(pep_vals):>8.1f}")
    print(f"  {'receptor length':<20}  {min(rec_vals):>8} {max(rec_vals):>8} "
          f"{sum(rec_vals)/len(rec_vals):>8.1f}")
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="Generate a peptide-protein benchmark sample from PDBbind outputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ---- Modes (mutually exclusive) ----
    mode_grp = p.add_argument_group("modes (one required)")
    mode = mode_grp.add_mutually_exclusive_group(required=True)
    mode.add_argument("--sample", metavar="N", type=int,
                      help="Sample N entries from the filtered pool.")
    mode.add_argument("--no-sample", action="store_true",
                      help="Take the full filtered pool without sampling.")
    mode.add_argument("--dry-run", action="store_true",
                      help="Print pool statistics after filtering and exit (no output written).")
    mode.add_argument("--pick", metavar="PDB_ID",
                      help="Extract a single system by PDB ID, bypassing all filters and sampling.")

    # ---- Sampling ----
    samp_grp = p.add_argument_group("sampling  (only with --sample)")
    samp_grp.add_argument("--sampling", choices=list(SAMPLING_METHODS), default="stratified",
                          metavar="METHOD",
                          help="Algorithm: stratified (default), kennard-stone, kde.")
    samp_grp.add_argument("--dims", nargs="+", default=["dg", "pep_len"], metavar="DIM",
                          help=f"Feature dimensions. Available: {', '.join(AVAILABLE_DIMS)}. "
                               "Default: dg pep_len. kde always uses only the first dim.")
    samp_grp.add_argument("--seed", type=int, default=None,
                          help="Random seed for reproducibility.")

    # ---- Filters ----
    filt_grp = p.add_argument_group("filters  (applied before sampling; all optional)")
    filt_grp.add_argument("--single-chain", action=argparse.BooleanOptionalAction, default=True,
                          help="Keep only single receptor-chain entries (default: on). Use --no-single-chain to include multi-chain receptors.")
    filt_grp.add_argument("--min-receptor", type=int, metavar="N",
                          help="Minimum receptor residue count.")
    filt_grp.add_argument("--max-receptor", type=int, metavar="N",
                          help="Maximum receptor residue count.")
    filt_grp.add_argument("--min-peptide", type=int, metavar="N",
                          help="Minimum peptide residue count.")
    filt_grp.add_argument("--max-peptide", type=int, metavar="N",
                          help="Maximum peptide residue count.")
    filt_grp.add_argument("--min-dg", type=float, metavar="F",
                          help="Minimum dG (kcal/mol), e.g. --min-dg -20.")
    filt_grp.add_argument("--max-dg", type=float, metavar="F",
                          help="Maximum dG (kcal/mol), e.g. --max-dg -5.")
    filt_grp.add_argument("--binding-type", choices=["kd", "ki", "both"], default="both",
                          help="Filter by measurement type: kd, ki, or both (default).")
    filt_grp.add_argument("--max-added", type=int, metavar="N",
                          help="Max residues PDBFixer added to peptide (requires filtered_systems.tsv).")
    filt_grp.add_argument("--max-receptor-added", type=int, metavar="N",
                          help="Max residues PDBFixer added to receptor (requires filtered_systems.tsv).")
    filt_grp.add_argument("--max-resolution", type=float, metavar="F",
                          help="Maximum resolution in Å (e.g. --max-resolution 2.5). "
                               "Excludes NMR/NA entries.")
    filt_grp.add_argument("--min-resolution", type=float, metavar="F",
                          help="Minimum resolution in Å.")
    filt_grp.add_argument("--exclude-nmr", action="store_true",
                          help="Exclude entries with no numeric resolution (NMR, NA).")
    filt_grp.add_argument("--min-year", type=int, metavar="N",
                          help="Minimum deposition year (e.g. --min-year 2000). "
                               "Excludes entries with NA year.")
    filt_grp.add_argument("--max-year", type=int, metavar="N",
                          help="Maximum deposition year.")
    filt_grp.add_argument("--exclude", metavar="FILE",
                          help="File with one pdb_id per line to exclude.")
    filt_grp.add_argument("--peptide-topology",
                          choices=["linear", "cyclic_disulfide", "both"], default="both",
                          help="Filter by peptide topology: linear, cyclic_disulfide, or both "
                               "(default: both).")

    # ---- Redundancy ----
    red_grp = p.add_argument_group(
        "redundancy  (applied after quality filters, before sampling)"
    )
    red_grp.add_argument(
        "--seq-id-receptor",
        choices=["strict", "medium", "lenient", "custom"],
        metavar="PRESET",
        help="Receptor sequence redundancy preset: strict (40%%), medium (50%%), "
             "lenient (70%%), custom. Uses CD-HIT (min threshold: 40%%). "
             "Requires receptor_seq in summary.tsv.",
    )
    red_grp.add_argument(
        "--seq-id-receptor-threshold", type=float, metavar="F",
        help="Receptor identity threshold for --seq-id-receptor custom (e.g. 0.40).",
    )
    red_grp.add_argument(
        "--seq-id-peptide",
        choices=["strict", "medium", "lenient", "custom"],
        metavar="PRESET",
        help="Peptide sequence redundancy preset: strict (70%%), medium (80%%), "
             "lenient (90%%), custom. Uses Biopython pairwise. Requires peptide_seq.",
    )
    red_grp.add_argument(
        "--seq-id-peptide-threshold", type=float, metavar="F",
        help="Peptide identity threshold for --seq-id-peptide custom (e.g. 0.80).",
    )

    # ---- Structural redundancy ----
    struct_grp = p.add_argument_group(
        "structural redundancy  (applied after seq redundancy; requires fixed_cif/)"
    )
    struct_grp.add_argument(
        "--tm-score-receptor", type=float, metavar="F",
        help="Receptor TM-score threshold (e.g. 0.5). Pairs >= threshold are redundant. "
             "Recommended: 0.5 = same fold. Requires tmtools (pip install tmtools).",
    )

    # ---- Input ----
    inp_grp = p.add_argument_group("input")
    inp_grp.add_argument("--summary-tsv", metavar="FILE", default=str(SUMMARY_TSV),
                         help="Path to summary TSV. Use summary_merged.tsv to include "
                              "external datasets (default: scrape/outputs/summary.tsv).")

    # ---- Output ----
    out_grp = p.add_argument_group("output")
    out_grp.add_argument("--no-save-pools", action="store_true",
                         help="Do not save pool_quality.tsv and pool_redundancy.tsv "
                              "in the run directory (saved by default).")
    out_grp.add_argument("--full-metadata", action="store_true",
                         help="Write all columns from summary.tsv to sample.tsv "
                              "(default: compact summary columns only).")

    # ---- Export ----
    exp_grp = p.add_argument_group("export  (all optional; can be combined)")
    exp_grp.add_argument("--link-cifs", action="store_true",
                         help="Symlink structure files into run_XX/cif/ (or pdb/).")
    exp_grp.add_argument("--copy-cifs", action="store_true",
                         help="Copy structure files into run_XX/cif/ (or pdb/).")
    exp_grp.add_argument("--protonate", action="store_true",
                         help="Add missing hydrogens at the given pH. Implies file copy.")
    exp_grp.add_argument("--ph", type=float, default=7.0,
                         help="pH for protonation (default: 7.0). Only used with --protonate.")
    exp_grp.add_argument("--format", choices=["cif", "pdb"], default="cif",
                         help="Output format: cif (default) or pdb.")
    exp_grp.add_argument("--export-fasta", action="store_true",
                         help="Write peptide_sequences.fasta and receptor_sequences.fasta "
                              "to the run directory. Requires peptide_seq / receptor_seq "
                              "columns in summary.tsv.")

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    for d in args.dims:
        if d not in AVAILABLE_DIMS:
            parser.error(f"Unknown dimension '{d}'. Available: {', '.join(AVAILABLE_DIMS)}")

    summary_tsv = Path(args.summary_tsv)
    if not summary_tsv.exists():
        print(f"  Error: {summary_tsv} not found — run scrape steps 1-3 first.")
        sys.exit(1)

    # Load and enrich
    all_rows = load_tsv(summary_tsv)
    # Capture public column order before enrich_rows adds _ private keys
    full_cols = [k for k in all_rows[0]] if all_rows else None
    all_rows = enrich_rows(all_rows)
    print(f"  Loaded {len(all_rows)} entries from {summary_tsv}")

    # --pick: bypass filters and sampling entirely
    if args.pick:
        pdb_id = args.pick.lower()
        match = [r for r in all_rows if r["pdb_id"] == pdb_id]
        if not match:
            print(f"  Error: '{pdb_id}' not found in {summary_tsv}.")
            sys.exit(1)
        selected = match
        print(f"  Picked {pdb_id}")
        print_sample_summary(selected)
        run_dir = next_run_dir(OUTPUTS_DIR)
        LOGS_DIR.mkdir(parents=True, exist_ok=True)
        tee = _Tee(LOGS_DIR / f"{run_dir.name}.log")
        sys.stdout = tee
        print(f"  Output → {run_dir}")
        write_sample(selected, run_dir,
                     columns=full_cols if args.full_metadata else None)
        config = {
            "mode": "pick",
            "pdb_id": pdb_id,
            "protonate": args.protonate,
            "ph": args.ph if args.protonate else None,
            "sources": {
                "summary_tsv":   str(summary_tsv),
                "fixed_cif_dir": str(FIXED_CIF_DIR),
            },
        }
        save_run_config(config, run_dir)
        _export_cifs(selected, run_dir, args)
        tee.close()
        return

    # Quality filters
    pool_quality = apply_quality_filters(all_rows, args, FILTERED_TSV)

    if not pool_quality:
        print("  No entries remain after quality filtering.")
        sys.exit(1)

    # Sequence redundancy
    pool_seq = apply_redundancy_filters(pool_quality, args)

    if not pool_seq:
        print("  No entries remain after sequence redundancy filtering.")
        sys.exit(1)

    # Structural redundancy
    pool_redundancy = apply_structural_redundancy_filters(pool_seq, args, FIXED_CIF_DIR)

    if not pool_redundancy:
        print("  No entries remain after structural redundancy filtering.")
        sys.exit(1)

    # Dry run — just show pool stats
    if args.dry_run:
        print_pool_summary(pool_redundancy)
        return

    # Set up output dir and logging
    run_dir = next_run_dir(OUTPUTS_DIR)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    tee = _Tee(LOGS_DIR / f"{run_dir.name}.log")
    sys.stdout = tee
    print(f"  Output → {run_dir}")

    # Save intermediate pools
    if not getattr(args, "no_save_pools", False):
        run_dir.mkdir(parents=True, exist_ok=True)
        save_tsv(pool_quality, run_dir / "pool_quality.tsv", POOL_COLS)
        print(f"  Written {run_dir / 'pool_quality.tsv'}  ({len(pool_quality)} entries)")
        save_tsv(pool_redundancy, run_dir / "pool_redundancy.tsv", POOL_COLS)
        print(f"  Written {run_dir / 'pool_redundancy.tsv'}  ({len(pool_redundancy)} entries)")

    # Sample or take all
    if args.no_sample:
        selected = pool_redundancy
        print(f"  --no-sample: keeping all {len(selected)} entries")
    else:
        selected = run_sampling(pool_redundancy, args.sample, method=args.sampling,
                                dims=args.dims, seed=args.seed)
        print(f"  Sampled {len(selected)} of {len(pool_redundancy)} entries "
              f"(method: {method_label(args.sampling, args.dims)})")

    selected.sort(key=lambda r: r["_dg"])
    print_sample_summary(selected)

    write_sample(selected, run_dir,
                 columns=full_cols if args.full_metadata else None)

    config = {
        "mode":            "pick" if args.pick else ("no-sample" if args.no_sample else "sample"),
        "sample_size":     len(selected),
        "pool_size":       len(pool_redundancy),
        "sampling_method": method_label(args.sampling, args.dims) if not args.no_sample else "none",
        "redundancy": {
            "seq_id_receptor":           args.seq_id_receptor,
            "seq_id_receptor_threshold": args.seq_id_receptor_threshold,
            "seq_id_peptide":            args.seq_id_peptide,
            "seq_id_peptide_threshold":  args.seq_id_peptide_threshold,
        },
        "dims":            args.dims if not args.no_sample else None,
        "seed":            args.seed,
        "protonate":       args.protonate,
        "ph":              args.ph if args.protonate else None,
        "filters": {
            "single_chain":        args.single_chain,
            "min_receptor":        args.min_receptor,
            "max_receptor":        args.max_receptor,
            "min_peptide":         args.min_peptide,
            "max_peptide":         args.max_peptide,
            "min_dg":              args.min_dg,
            "max_dg":              args.max_dg,
            "binding_type":        args.binding_type,
            "max_added":           args.max_added,
            "max_receptor_added":  args.max_receptor_added,
            "max_resolution":      args.max_resolution,
            "min_resolution":      args.min_resolution,
            "exclude_nmr":         args.exclude_nmr,
            "exclude":             args.exclude,
        },
        "sources": {
            "summary_tsv":   str(SUMMARY_TSV),
            "filtered_tsv":  str(FILTERED_TSV),
            "fixed_cif_dir": str(FIXED_CIF_DIR),
        },
    }
    save_run_config(config, run_dir)
    _export_cifs(selected, run_dir, args)
    if args.export_fasta:
        write_fasta(selected, run_dir)
    tee.close()


def _export_cifs(selected, run_dir, args):
    """Handle structure export: protonate (implies copy), copy, symlink, or skip."""
    fmt = args.format
    if args.protonate:
        protonate_cifs(selected, run_dir, FIXED_CIF_DIR, ph=args.ph, fmt=fmt)
    elif args.copy_cifs:
        link_cifs(selected, run_dir, FIXED_CIF_DIR, copy=True, fmt=fmt)
    elif args.link_cifs:
        link_cifs(selected, run_dir, FIXED_CIF_DIR, copy=False, fmt=fmt)


if __name__ == "__main__":
    main()
