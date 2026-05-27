"""
gen_sample.py — Generate a peptide-protein benchmark sample from PDBbind outputs.

Reads summary.tsv (and optionally filtered_systems.tsv) produced by scrape_pdb_bind,
applies filters, samples entries, and writes outputs to gen_sample/outputs/run_XX/.

Usage examples:
  python gen_sample.py --sample 20
  python gen_sample.py --sample 20 --single-chain --max-peptide 15
  python gen_sample.py --sample 20 --sampling kennard-stone --seed 42
  python gen_sample.py --sample 20 --sampling kde
  python gen_sample.py --no-sample
  python gen_sample.py --dry-run
  python gen_sample.py --sample 20 --link-cifs
  python gen_sample.py --sample 20 --copy-cifs
  python gen_sample.py --sample 20 --protonate
  python gen_sample.py --sample 20 --protonate --ph 6.5
  python gen_sample.py --pick 1abc
  python gen_sample.py --pick 1abc --protonate --ph 7.4
"""

import argparse
import sys
from pathlib import Path

# Allow running from scripts/ directly or from the gen_sample root.
_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, enrich_rows
from modules.filters import apply_all_filters
from modules.sampling import run_sampling, SAMPLING_METHODS
from modules.export import next_run_dir, write_sample, save_run_config, link_cifs, protonate_cifs

# ---------------------------------------------------------------------------
# Paths (relative to scrape_pdb_bind, which lives next to gen_sample)
# ---------------------------------------------------------------------------

_SCRAPE_DIR    = _root.parent / "scrape_pdb_bind"
SUMMARY_TSV    = _SCRAPE_DIR / "summary.tsv"
FILTERED_TSV   = _SCRAPE_DIR / "filtered_systems.tsv"
FIXED_CIF_DIR  = _SCRAPE_DIR / "fixed_cif"
OUTPUTS_DIR    = _root / "outputs"


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

    # Mode (mutually exclusive)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--sample", metavar="N", type=int,
                      help="Sample N entries from the filtered pool.")
    mode.add_argument("--no-sample", action="store_true",
                      help="Take the full filtered pool without sampling.")
    mode.add_argument("--dry-run", action="store_true",
                      help="Print pool statistics after filtering and exit (no output written).")
    mode.add_argument("--pick", metavar="PDB_ID",
                      help="Extract a single system by PDB ID, bypassing all filters and sampling.")

    # Sampling options
    p.add_argument("--sampling", choices=list(SAMPLING_METHODS), default="diverse",
                   metavar="METHOD",
                   help="Sampling algorithm: diverse (default), kennard-stone, kde.")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducibility.")

    # Filters
    p.add_argument("--single-chain", action=argparse.BooleanOptionalAction, default=True,
                   help="Keep only entries with a single receptor chain (default: on). Use --no-single-chain to include multi-chain receptors.")
    p.add_argument("--min-receptor", type=int, metavar="N",
                   help="Minimum receptor residue count.")
    p.add_argument("--max-receptor", type=int, metavar="N",
                   help="Maximum receptor residue count.")
    p.add_argument("--min-peptide", type=int, metavar="N",
                   help="Minimum peptide residue count.")
    p.add_argument("--max-peptide", type=int, metavar="N",
                   help="Maximum peptide residue count.")
    p.add_argument("--min-dg", type=float, metavar="F",
                   help="Minimum dG (kcal/mol), e.g. --min-dg -20.")
    p.add_argument("--max-dg", type=float, metavar="F",
                   help="Maximum dG (kcal/mol), e.g. --max-dg -5.")
    p.add_argument("--binding-type", choices=["kd", "ki", "both"], default="both",
                   help="Filter by binding measurement type (default: both).")
    p.add_argument("--max-added", type=int, metavar="N",
                   help="Max residues PDBFixer added to peptide (requires filtered_systems.tsv).")
    p.add_argument("--max-receptor-added", type=int, metavar="N",
                   help="Max residues PDBFixer added to receptor (requires filtered_systems.tsv).")
    p.add_argument("--exclude", metavar="FILE",
                   help="File with one pdb_id per line to exclude.")

    # Export options
    p.add_argument("--link-cifs", action="store_true",
                   help="Symlink selected CIF files into the run directory.")
    p.add_argument("--copy-cifs", action="store_true",
                   help="Copy (rather than symlink) selected CIF files into the run directory.")
    p.add_argument("--protonate", action="store_true",
                   help="Add missing hydrogens and protonate at the given pH. Implies --copy-cifs.")
    p.add_argument("--ph", type=float, default=7.0,
                   help="pH for protonation (default: 7.0). Only used with --protonate.")
    p.add_argument("--format", choices=["cif", "pdb"], default="cif",
                   help="Output structure format: cif (default) or pdb.")

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    if not SUMMARY_TSV.exists():
        print(f"  Error: {SUMMARY_TSV} not found — run scrape_pdb_bind steps 1-3 first.")
        sys.exit(1)

    # Load and enrich
    all_rows = load_tsv(SUMMARY_TSV)
    all_rows = enrich_rows(all_rows)
    print(f"  Loaded {len(all_rows)} entries from {SUMMARY_TSV}")

    # --pick: bypass filters and sampling entirely
    if args.pick:
        pdb_id = args.pick.lower()
        match = [r for r in all_rows if r["pdb_id"] == pdb_id]
        if not match:
            print(f"  Error: '{pdb_id}' not found in {SUMMARY_TSV}.")
            sys.exit(1)
        selected = match
        print(f"  Picked {pdb_id}")
        print_sample_summary(selected)
        run_dir = next_run_dir(OUTPUTS_DIR)
        print(f"  Output → {run_dir}")
        write_sample(selected, run_dir)
        config = {
            "mode": "pick",
            "pdb_id": pdb_id,
            "protonate": args.protonate,
            "ph": args.ph if args.protonate else None,
            "sources": {
                "summary_tsv":   str(SUMMARY_TSV),
                "fixed_cif_dir": str(FIXED_CIF_DIR),
            },
        }
        save_run_config(config, run_dir)
        _export_cifs(selected, run_dir, args)
        return

    # Filter
    rows = apply_all_filters(all_rows, args, FILTERED_TSV)

    if not rows:
        print("  No entries remain after filtering.")
        sys.exit(1)

    # Dry run — just show pool stats
    if args.dry_run:
        print_pool_summary(rows)
        return

    # Sample or take all
    if args.no_sample:
        selected = rows
        print(f"  --no-sample: keeping all {len(selected)} entries")
    else:
        selected = run_sampling(rows, args.sample, method=args.sampling, seed=args.seed)
        print(f"  Sampled {len(selected)} of {len(rows)} entries "
              f"(method: {args.sampling})")

    selected.sort(key=lambda r: r["_dg"])
    print_sample_summary(selected)

    # Write outputs
    run_dir = next_run_dir(OUTPUTS_DIR)
    print(f"  Output → {run_dir}")

    write_sample(selected, run_dir)

    config = {
        "mode":            "pick" if args.pick else ("no-sample" if args.no_sample else "sample"),
        "sample_size":     len(selected),
        "pool_size":       len(rows),
        "sampling_method": args.sampling if not args.no_sample else "none",
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
