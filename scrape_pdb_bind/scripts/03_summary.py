"""
Step 3: Merge filtered_systems.tsv with peptide_entries.tsv (dG, binding type,
deposition year) and write a human-readable summary table.

Output columns:
    pdb_id          - PDB accession code
    receptor        - e.g. "AB:300" (chains A+B, 300 total residues)
    peptide         - e.g. "C:12" (post-fix residue count)
    dG_kcal_mol     - DeltaG at 300 K (kcal/mol)
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

from modules.utils import load_tsv, save_tsv
from modules.metadata import ss_percentages, classify_peptide_ss

OUTPUTS_DIR   = _root / "outputs"
LOGS_DIR      = _root / "logs"
ENTRIES_TSV   = OUTPUTS_DIR / "peptide_entries.tsv"
FILTERED_TSV  = OUTPUTS_DIR / "filtered_systems.tsv"
OUT_TSV       = OUTPUTS_DIR / "summary.tsv"
LOG_FILE      = LOGS_DIR / "03_summary.log"

OUT_COLS = [
    "pdb_id", "receptor", "peptide", "dG_kcal_mol", "binding_type",
    "resolution", "year",
    "peptide_added", "peptide_seq",
    "peptide_ss", "peptide_ss_dssp",
    "peptide_pct_helix", "peptide_pct_beta", "peptide_pct_coil",
    "peptide_ss_class",
    "receptor_added", "receptor_seq",
    "receptor_ss", "receptor_ss_dssp",
    "receptor_pct_helix", "receptor_pct_beta", "receptor_pct_coil",
    "receptor_func_class",
]


def main():
    parser = argparse.ArgumentParser(
        description="Step 3: Merge filtered_systems.tsv with peptide_entries.tsv "
                    "and write summary.tsv.")
    parser.add_argument("--clean", action="store_true",
                        help="Remove summary.tsv.")
    args = parser.parse_args()

    if args.clean:
        if OUT_TSV.exists():
            OUT_TSV.unlink()
            print(f"Removed {OUT_TSV}")
        return

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    entries_map   = {r["pdb_id"]: r for r in load_tsv(ENTRIES_TSV)}
    filtered_rows = load_tsv(FILTERED_TSV)
    accepted      = [r for r in filtered_rows if r.get("accepted") == "Y"]
    print(f"Accepted systems: {len(accepted)}")

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
            ),
            "receptor_added":   row.get("receptor_added", ""),
            "receptor_seq":     row.get("receptor_seq", ""),
            "receptor_ss":      row.get("receptor_ss", "NA"),
            "receptor_ss_dssp": row.get("receptor_ss_dssp", "NA"),
            "receptor_pct_helix":  ss_percentages(row.get("receptor_ss_dssp", "NA"))[0],
            "receptor_pct_beta":   ss_percentages(row.get("receptor_ss_dssp", "NA"))[1],
            "receptor_pct_coil":   ss_percentages(row.get("receptor_ss_dssp", "NA"))[2],
            "receptor_func_class": row.get("receptor_func_class", "NA"),
        }))

    out_rows.sort(key=lambda x: x[0])
    rows = [r for _, r in out_rows]

    save_tsv(rows, OUT_TSV, OUT_COLS)
    print(f"Written {OUT_TSV}")

    with open(LOG_FILE, "w", encoding="ascii") as fh:
        fh.write(f"Accepted: {len(rows)}\n")
        fh.write("\t".join(OUT_COLS) + "\n")
        for r in rows[:10]:
            fh.write("\t".join(r.get(c, "") for c in OUT_COLS) + "\n")


if __name__ == "__main__":
    main()
