"""
Step 3: Merge filtered_systems.tsv with peptide_entries.tsv (dG, binding type,
deposition year) and write a human-readable summary table.

Output columns:
    pdb_id          - PDB accession code
    receptor        - e.g. "AB:300" (chains A+B, 300 total residues)
    peptide         - e.g. "C:12" (post-fix residue count)
    dG_kcal_mol     - DeltaG = RT ln(Kd/Ki) (kcal/mol)
    temperature     - temperature used for ΔG calculation (K); from index notes or fallback
    binding_type    - Kd / Ki
    resolution      - crystallographic resolution in Å, or NA
    year            - deposition year

Input:  outputs/peptide_entries.tsv, outputs/filtered_systems.tsv
Output: outputs/summary.tsv, logs/03_summary.log
"""

import argparse
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, save_tsv, setup_logging
from modules.metadata import ss_percentages, classify_peptide_ss

OUTPUTS_DIR   = _root / "outputs"
LOGS_DIR      = _root / "logs"
ENTRIES_TSV   = OUTPUTS_DIR / "peptide_entries.tsv"
FILTERED_TSV  = OUTPUTS_DIR / "filtered_systems.tsv"
OUT_TSV       = OUTPUTS_DIR / "summary.tsv"
LOG_FILE      = LOGS_DIR / "03_summary.log"

OUT_COLS = [
    "pdb_id", "receptor", "peptide", "dG_kcal_mol", "binding_type",
    "resolution", "year", "temperature",
    "peptide_added", "peptide_seq",
    "peptide_ss", "peptide_ss_dssp",
    "peptide_pct_helix", "peptide_pct_beta", "peptide_pct_coil",
    "peptide_ss_class",
    "receptor_added", "receptor_seq",
    "receptor_ss", "receptor_ss_dssp",
    "receptor_pct_helix", "receptor_pct_beta", "receptor_pct_coil",
    "receptor_func_class",
    "n_contacts", "n_hbonds",
]


def main():
    parser = argparse.ArgumentParser(
        description="Step 3: Merge filtered_systems.tsv with peptide_entries.tsv "
                    "and write summary.tsv.")
    parser.add_argument("--clean", action="store_true",
                        help="Remove output summary TSV.")
    parser.add_argument("--entries", metavar="FILE", default=str(ENTRIES_TSV),
                        help=f"peptide_entries TSV (default: {ENTRIES_TSV.name}).")
    parser.add_argument("--filtered", metavar="FILE", default=str(FILTERED_TSV),
                        help=f"filtered_systems TSV (default: {FILTERED_TSV.name}).")
    parser.add_argument("--output", metavar="FILE", default=str(OUT_TSV),
                        help=f"Output summary TSV (default: {OUT_TSV.name}).")
    parser.add_argument("--log", metavar="FILE", default=str(LOG_FILE),
                        help=f"Log file path.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print a brief start/done line to stdout.")
    args = parser.parse_args()

    entries_tsv  = Path(args.entries)
    filtered_tsv = Path(args.filtered)
    out_tsv      = Path(args.output)
    log_file     = Path(args.log)

    if args.clean:
        if out_tsv.exists():
            out_tsv.unlink()
            print(f"Removed {out_tsv}")
        return

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    _logger, stdout = setup_logging(log_file, args.verbose)
    if args.verbose:
        stdout.write("[scrape/03-summary] starting ...\n")
        stdout.flush()

    entries_rows  = load_tsv(entries_tsv)
    entries_map   = {r["pdb_id"]: r for r in entries_rows}
    filtered_rows = load_tsv(filtered_tsv)
    accepted      = [r for r in filtered_rows if r.get("accepted") == "Y"]
    print(f"Accepted systems: {len(accepted)}")

    # Auto-detect extra columns from entries TSV beyond the standard set
    _standard_entry_cols = {"pdb_id", "n_mer", "year", "binding_type",
                             "dG_kcal", "dG_approx", "step1_reject", "resolution"}
    extra_entry_cols = [c for c in (entries_rows[0].keys() if entries_rows else [])
                        if c not in _standard_entry_cols and c not in OUT_COLS]

    out_rows = []
    for row in accepted:
        pdb_id = row["pdb_id"]
        entry  = entries_map.get(pdb_id, {})

        rec_chains = "".join(row.get("receptor_chains", "").split(","))
        receptor   = f"{rec_chains}:{row.get('receptor_length', '')}"
        peptide    = f"{row.get('peptide_chain', '')}:{row.get('peptide_length', '')}"

        dG_raw = entry.get("dG_kcal", "NA")
        try:
            dG_sort = float(dG_raw)
        except (ValueError, TypeError):
            dG_sort = float("inf")

        # Prefer CIF-based resolution; fall back to INDEX value
        resolution = row.get("resolution", "NA")
        if resolution == "NA":
            resolution = entry.get("resolution", "NA")

        out_rows.append((dG_sort, {
            "pdb_id":         pdb_id,
            "receptor":       receptor,
            "peptide":        peptide,
            "dG_kcal_mol":    dG_raw,
            "binding_type":   entry.get("binding_type", "NA"),
            "resolution":     resolution,
            "year":           entry.get("year", "NA"),
            "peptide_added":    row.get("peptide_added", ""),
            "peptide_seq":      row.get("peptide_seq", ""),
            "peptide_ss":       row.get("peptide_ss", "NA"),
            "peptide_ss_dssp":  row.get("peptide_ss_dssp", "NA"),
            "peptide_pct_helix": ss_percentages(row.get("peptide_ss_dssp", "NA"))[0],
            "peptide_pct_beta":  ss_percentages(row.get("peptide_ss_dssp", "NA"))[1],
            "peptide_pct_coil":  ss_percentages(row.get("peptide_ss_dssp", "NA"))[2],
            "peptide_ss_class":  classify_peptide_ss(
                ss_percentages(row.get("peptide_ss_dssp", "NA"))[0],
                ss_percentages(row.get("peptide_ss_dssp", "NA"))[1],
                ss_percentages(row.get("peptide_ss_dssp", "NA"))[2],
            ),
            "receptor_added":   row.get("receptor_added", ""),
            "receptor_seq":     row.get("receptor_seq", ""),
            "receptor_ss":      row.get("receptor_ss", "NA"),
            "receptor_ss_dssp": row.get("receptor_ss_dssp", "NA"),
            "receptor_pct_helix":  ss_percentages(row.get("receptor_ss_dssp", "NA"))[0],
            "receptor_pct_beta":   ss_percentages(row.get("receptor_ss_dssp", "NA"))[1],
            "receptor_pct_coil":   ss_percentages(row.get("receptor_ss_dssp", "NA"))[2],
            "receptor_func_class": row.get("receptor_func_class", "NA"),
            "n_contacts":          row.get("n_contacts", ""),
            "n_hbonds":            row.get("n_hbonds", ""),
            **{c: entry.get(c, "") for c in extra_entry_cols},
        }))

    out_rows.sort(key=lambda x: x[0])
    rows = [r for _, r in out_rows]

    all_cols = OUT_COLS + extra_entry_cols
    save_tsv(rows, out_tsv, all_cols)
    print(f"Written {out_tsv}")
    print(f"Accepted: {len(rows)}")
    print("\t".join(all_cols))
    for r in rows[:10]:
        print("\t".join(r.get(c, "") for c in all_cols))

    if args.verbose:
        stdout.write(f"[scrape/03-summary] done  —  {len(rows)} entries\n")
        stdout.flush()


if __name__ == "__main__":
    main()
