"""
Step 2: For each downloaded CIF, select the best biological assembly,
identify the peptide and receptor chains, verify AMBER14 compatibility,
fill missing residues, extract sequences/resolution/secondary structure,
and write fixed CIFs to fixed_cif/.

Input:  outputs/peptide_entries.tsv, pdb_cache/<pdb_id>-<N>.cif
Output: outputs/filtered_systems.tsv, fixed_cif/<pdb_id>.cif
        logs/02_filter.log
"""

import argparse
import shutil
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv
from modules.filter import process_entry

OUTPUTS_DIR  = _root / "outputs"
LOGS_DIR     = _root / "logs"
PDB_CACHE    = _root / "pdb_cache"
FIXED_DIR    = _root / "fixed_cif"
IN_TSV       = OUTPUTS_DIR / "peptide_entries.tsv"
OUT_TSV      = OUTPUTS_DIR / "filtered_systems.tsv"
LOG_FILE     = LOGS_DIR / "02_filter.log"

OUT_COLS = [
    "pdb_id", "n_mer", "accepted", "assembly_num",
    "peptide_chain", "peptide_length", "peptide_added", "peptide_seq",
    "peptide_ss", "peptide_ss_dssp",
    "receptor_chains", "receptor_length", "receptor_added", "receptor_seq",
    "receptor_ss", "receptor_ss_dssp",
    "resolution",
    "reject_reason",
]


def load_entries(tsv_path):
    entries = {}
    for row in load_tsv(tsv_path):
        if row.get("step1_reject"):
            continue
        entries[row["pdb_id"]] = {
            "n_mer":      int(row["n_mer"]),
            "resolution": row.get("resolution", "NA"),
        }
    return entries


def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Filter entries, check AMBER14 compatibility, write fixed CIFs.")
    parser.add_argument("--clean", action="store_true",
                        help="Remove filtered_systems.tsv and fixed_cif/.")
    args = parser.parse_args()

    if args.clean:
        if OUT_TSV.exists():
            OUT_TSV.unlink()
            print(f"Removed {OUT_TSV}")
        if FIXED_DIR.exists():
            shutil.rmtree(FIXED_DIR)
            print(f"Removed {FIXED_DIR}/")
        return

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Reading {IN_TSV} ...")
    entries = load_entries(IN_TSV)
    print(f"  {len(entries)} entries to process")

    accepted_count = 0
    rejected_count = 0
    n_total = len(entries)

    with open(OUT_TSV, "w", encoding="ascii") as tsv_fh, \
         open(LOG_FILE, "w", encoding="ascii") as log_fh:

        tsv_fh.write("\t".join(OUT_COLS) + "\n")

        def log(msg):
            print(msg)
            log_fh.write(msg + "\n")

        for i, (pdb_id, entry_data) in enumerate(entries.items(), 1):
            n_mer           = entry_data["n_mer"]
            index_resolution = entry_data["resolution"]
            result = process_entry(pdb_id, n_mer, PDB_CACHE, FIXED_DIR)
            # Assembly CIFs lack experimental metadata — use INDEX resolution
            if result.get("resolution", "NA") == "NA":
                result["resolution"] = index_resolution

            if result["accepted"]:
                accepted_count += 1
                pep_id    = result["peptide_chain"]
                pep_added = result["residues_added"].get(pep_id, 0)
                rec_added = sum(result["residues_added"].get(cid, 0)
                                for cid in result["receptor_chains"])
                rec_chains_str = ",".join(result["receptor_chains"])

                row = {
                    "pdb_id":          pdb_id,
                    "n_mer":           str(n_mer),
                    "accepted":        "Y",
                    "assembly_num":    str(result["assembly_num"]),
                    "peptide_chain":   pep_id,
                    "peptide_length":  str(result["peptide_length"]),
                    "peptide_added":   str(pep_added),
                    "peptide_seq":     result.get("peptide_seq", ""),
                    "peptide_ss":      result.get("peptide_ss", "NA"),
                    "peptide_ss_dssp": result.get("peptide_ss_dssp", "NA"),
                    "receptor_chains": rec_chains_str,
                    "receptor_length": str(result["receptor_length"]),
                    "receptor_added":  str(rec_added),
                    "receptor_seq":    result.get("receptor_seq", ""),
                    "receptor_ss":     result.get("receptor_ss", "NA"),
                    "receptor_ss_dssp": result.get("receptor_ss_dssp", "NA"),
                    "resolution":      result.get("resolution", "NA"),
                    "reject_reason":   "",
                }

                pep_str = f"{pep_id}:{result['peptide_length'] - pep_added}+{pep_added}"
                rec_str = "  ".join(
                    f"{cid}:{result['receptor_chain_lengths'][cid] - result['residues_added'].get(cid,0)}"
                    f"+{result['residues_added'].get(cid,0)}"
                    for cid in result["receptor_chains"]
                )
                ions_all = sorted(set(
                    r for resnames in result["ion_chains"].values() for r in resnames
                ))
                ions_str = f"  ions=[{','.join(ions_all)}]" if ions_all else ""
                log(f"  [{i}/{n_total}] PASS   {pdb_id}  asm={result['assembly_num']}  "
                    f"pep={pep_str}  rec={rec_str}  res={result.get('resolution','NA')}"
                    f"{ions_str}  (accepted={accepted_count} rejected={rejected_count})")
                if result["peptide_length"] != n_mer:
                    log(f"    WARNING: peptide length {result['peptide_length']} != n_mer {n_mer}")
            else:
                rejected_count += 1
                row = {c: "" for c in OUT_COLS}
                row.update({
                    "pdb_id":       pdb_id,
                    "n_mer":        str(n_mer),
                    "accepted":     "N",
                    "assembly_num": str(result.get("assembly_num") or ""),
                    "reject_reason": result["reject_reason"],
                })
                log(f"  [{i}/{n_total}] REJECT {pdb_id}: {result['reject_reason']}  "
                    f"(accepted={accepted_count} rejected={rejected_count})")

            tsv_fh.write("\t".join(row.get(c, "") for c in OUT_COLS) + "\n")

    print(f"\nDone. accepted={accepted_count} rejected={rejected_count} "
          f"total={accepted_count + rejected_count}")
    print(f"Written {OUT_TSV}")
    print(f"Log: {LOG_FILE}")


if __name__ == "__main__":
    main()
