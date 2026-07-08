"""
Step 2: For each downloaded CIF, select the best biological assembly,
identify the peptide and receptor chains, verify AMBER14 compatibility,
fill missing residues, extract sequences/resolution/secondary structure,
and write fixed CIFs to fixed_cif/.

Input:  outputs/peptide_entries.tsv, pdb_cache/<pdb_id>-<N>.cif
Output: outputs/filtered_systems.tsv, fixed_cif/<pdb_id>.cif
        logs/02_filter.log

Cache behaviour (default):
  - Entries already in filtered_systems.tsv are reused without reprocessing.
  - If only lightweight columns (receptor_description, receptor_func_class) are
    missing from an accepted entry, they are patched cheaply from the raw CIF
    without re-running PDBFixer or DSSP.
  - Use --force to bypass the cache and reprocess all entries.
"""

import argparse
import shutil
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, setup_logging
from modules.filter import process_entry

OUTPUTS_DIR  = _root / "outputs"
LOGS_DIR     = _root / "logs"
PDB_CACHE    = _root / "pdb_cache"
FIXED_DIR    = _root / "fixed_cif"
IN_TSV       = OUTPUTS_DIR / "peptide_entries.tsv"
OUT_TSV      = OUTPUTS_DIR / "filtered_systems.tsv"
LOG_FILE     = LOGS_DIR / "02_filter.log"

OUT_COLS = [
    # Identity & raw metadata
    "pdb_id", "n_mer", "accepted", "assembly_num",
    # Peptide raw
    "peptide_chain", "peptide_length", "peptide_added", "peptide_seq",
    # Receptor raw
    "receptor_chains", "receptor_length", "receptor_added", "receptor_seq",
    # Crystallographic & functional metadata
    "resolution", "receptor_description", "receptor_func_class",
    # Processing flags
    "mse_converted", "caps_stripped", "cap_min_dist_A", "peptide_topology",
    # ── Computed / predicted ──────────────────────────────────────────────────
    "peptide_ss", "peptide_ss_dssp",
    "receptor_ss", "receptor_ss_dssp",
    "n_contacts", "n_hbonds", "pep_sasa_burial",
    # Status
    "reject_reason",
]

# Columns whose presence indicates expensive work (PDBFixer + DSSP) was done.
_HEAVY_COLS = ("peptide_seq", "receptor_seq")

# Columns that can be patched cheaply without re-running PDBFixer or DSSP.
# receptor_description/func_class: from raw pdb_cache CIF.
# n_contacts/n_hbonds/peptide_topology/pep_sasa_burial: from already-built fixed_cif/.
_PATCH_COLS = frozenset(("receptor_description", "receptor_func_class",
                          "n_contacts", "n_hbonds", "peptide_topology",
                          "pep_sasa_burial"))

# Columns that are legitimately empty for many accepted entries (not a cache miss).
# mse_converted = "" means no MSE was found; caps_stripped = "" means no caps stripped.
_NULLABLE_COLS = frozenset(("mse_converted", "caps_stripped", "cap_min_dist_A"))


def _load_cache(tsv_path):
    """Load existing filtered_systems.tsv as {pdb_id: row_dict}."""
    if not tsv_path.exists():
        return {}
    return {r["pdb_id"]: r for r in load_tsv(tsv_path)}


def _cache_status(cached):
    """
    Classify the state of a cached row.
    Returns:
      'full'  — row is up-to-date; no processing needed.
      'patch' — heavy work done; only _PATCH_COLS are missing; add cheaply.
      'miss'  — needs full reprocessing.
    """
    if cached is None:
        return "miss"
    # Rejected entries are always fully cached (no columns to patch)
    if cached.get("accepted") == "N":
        return "full"
    # Accepted entry: verify heavy sentinel columns are non-empty
    if not all(cached.get(c, "") not in ("", "NA") for c in _HEAVY_COLS):
        return "miss"
    # Find which required columns are absent or empty — skip reject_reason
    # (always empty for accepted entries) and _NULLABLE_COLS (legitimately empty).
    missing = [c for c in OUT_COLS
               if c != "reject_reason"
               and c not in _NULLABLE_COLS
               and not cached.get(c, "")]
    if not missing:
        return "full"
    if all(c in _PATCH_COLS for c in missing):
        return "patch"
    return "miss"


def _patch_row(pdb_id, cached):
    """Patch missing lightweight columns for a cached accepted row.

    receptor_description/func_class: reads raw pdb_cache CIF (no PDBFixer/DSSP).
    n_contacts/n_hbonds: reads already-built fixed_cif/<pdb_id>.cif.
    """
    row = dict(cached)

    if not row.get("receptor_description") or not row.get("receptor_func_class"):
        from modules.metadata import extract_receptor_description, classify_receptor_function
        asm_num  = cached.get("assembly_num", "1") or "1"
        cif_path = PDB_CACHE / f"{pdb_id}-{asm_num}.cif"
        if not cif_path.exists():
            cif_path = PDB_CACHE / f"{pdb_id}-1.cif"
        rec_chain_ids = [c.strip() for c in cached.get("receptor_chains", "").split(",")
                         if c.strip()]
        if cif_path.exists() and rec_chain_ids:
            desc = extract_receptor_description(cif_path, rec_chain_ids)
            row["receptor_description"] = desc
            row["receptor_func_class"]  = classify_receptor_function(desc)
        else:
            row.setdefault("receptor_description", "NA")
            row.setdefault("receptor_func_class",  "NA")

    if not row.get("n_contacts") or not row.get("n_hbonds"):
        from modules.filter import compute_interface_contacts
        fixed_cif = FIXED_DIR / f"{pdb_id}.cif"
        if fixed_cif.exists():
            pep_chain  = row.get("peptide_chain", "")
            rec_chains = [c.strip() for c in row.get("receptor_chains", "").split(",")
                          if c.strip()]
            n_c, n_h = compute_interface_contacts(fixed_cif, pep_chain, rec_chains)
            row["n_contacts"] = str(n_c)
            row["n_hbonds"]   = str(n_h)
        else:
            row.setdefault("n_contacts", "")
            row.setdefault("n_hbonds",   "")

    if not row.get("peptide_topology"):
        from modules.filter import _detect_peptide_topology
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        fixed_cif = FIXED_DIR / f"{pdb_id}.cif"
        if fixed_cif.exists():
            pep_chain = row.get("peptide_chain", "")
            try:
                d = MMCIF2Dict(str(fixed_cif))
                row["peptide_topology"] = _detect_peptide_topology(d, pep_chain)
            except Exception:
                row["peptide_topology"] = "linear"
        else:
            row["peptide_topology"] = "linear"

    if not row.get("pep_sasa_burial"):
        from modules.filter import compute_sasa_burial
        fixed_cif = FIXED_DIR / f"{pdb_id}.cif"
        if fixed_cif.exists():
            pep_chain = row.get("peptide_chain", "")
            burial = compute_sasa_burial(fixed_cif, pep_chain)
            row["pep_sasa_burial"] = str(burial) if burial is not None else ""
        else:
            row.setdefault("pep_sasa_burial", "")

    return row


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
                        help="Remove output TSV and fixed_cif/ (default output only).")
    parser.add_argument("--force", action="store_true",
                        help="Ignore the existing cache and reprocess all entries from scratch.")
    parser.add_argument("--input", metavar="FILE", default=str(IN_TSV),
                        help=f"Input peptide_entries TSV (default: {IN_TSV.name}).")
    parser.add_argument("--output", metavar="FILE", default=str(OUT_TSV),
                        help=f"Output filtered_systems TSV (default: {OUT_TSV.name}).")
    parser.add_argument("--log", metavar="FILE", default=str(LOG_FILE),
                        help=f"Log file path (default: {LOG_FILE.name}).")
    parser.add_argument("--no-convert-mse", dest="convert_mse",
                        action="store_false", default=True,
                        help="Disable MSE→MET conversion (on by default; MSE is a "
                             "crystallographic Met substitute with identical biology).")
    parser.add_argument("--strip-caps", action="store_true",
                        help="Remove ACE/NME/NH2 terminal caps from peptide if their "
                             "atoms are >4.5 Å from the receptor interface. "
                             "Saves cap_min_dist_A to TSV for transparency.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print a brief start/done line to stdout.")
    args = parser.parse_args()

    in_tsv   = Path(args.input)
    out_tsv  = Path(args.output)
    log_file = Path(args.log)

    if args.clean:
        if out_tsv.exists():
            out_tsv.unlink()
            print(f"Removed {out_tsv}")
        # Only remove fixed_cif/ when using the default output path
        if out_tsv == OUT_TSV and FIXED_DIR.exists():
            shutil.rmtree(FIXED_DIR)
            print(f"Removed {FIXED_DIR}/")
        return

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    _logger, stdout = setup_logging(log_file, args.verbose)
    if args.verbose:
        stdout.write("[scrape/02-filter] starting ...\n")
        stdout.flush()

    # Load existing results before opening out_tsv for writing
    if args.force:
        cache = {}
        print("  --force: ignoring existing cache, reprocessing all entries.")
    else:
        cache = _load_cache(out_tsv)
        if cache:
            n_acc = sum(1 for r in cache.values() if r.get("accepted") == "Y")
            n_rej = len(cache) - n_acc
            print(f"  Cache loaded: {n_acc} accepted, {n_rej} rejected ({len(cache)} total).")

    print(f"Reading {in_tsv} ...")
    entries = load_entries(in_tsv)
    print(f"  {len(entries)} entries to process")

    accepted_count = 0
    rejected_count = 0
    cached_count   = 0
    patched_count  = 0
    n_total = len(entries)

    with open(out_tsv, "w", encoding="ascii") as tsv_fh:

        tsv_fh.write("\t".join(OUT_COLS) + "\n")

        def _write_row(row):
            tsv_fh.write("\t".join(row.get(c, "") for c in OUT_COLS) + "\n")

        for i, (pdb_id, entry_data) in enumerate(entries.items(), 1):
            n_mer            = entry_data["n_mer"]
            index_resolution = entry_data["resolution"]

            # ── Cache check ───────────────────────────────────────────────
            cached = cache.get(pdb_id)
            status = _cache_status(cached)

            if status == "full":
                row = {c: cached.get(c, "") for c in OUT_COLS}
                _write_row(row)
                if cached.get("accepted") == "Y":
                    accepted_count += 1
                else:
                    rejected_count += 1
                cached_count += 1
                continue

            if status == "patch":
                patched = _patch_row(pdb_id, cached)
                row = {c: patched.get(c, "") for c in OUT_COLS}
                _write_row(row)
                accepted_count += 1
                patched_count  += 1
                missing = [c for c in OUT_COLS if c != "reject_reason" and not cached.get(c, "")]
                print(f"  [{i}/{n_total}] PATCH  {pdb_id}  ({', '.join(missing)} added)")
                continue

            # ── Full processing ───────────────────────────────────────────
            result = process_entry(pdb_id, n_mer, PDB_CACHE, FIXED_DIR,
                                   convert_mse=args.convert_mse,
                                   strip_caps=args.strip_caps)
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
                    "receptor_ss":          result.get("receptor_ss", "NA"),
                    "receptor_ss_dssp":     result.get("receptor_ss_dssp", "NA"),
                    "resolution":           result.get("resolution", "NA"),
                    "receptor_description": result.get("receptor_description", "NA"),
                    "receptor_func_class":  result.get("receptor_func_class", "NA"),
                    "n_contacts":           str(result.get("n_contacts", "")),
                    "n_hbonds":             str(result.get("n_hbonds", "")),
                    "pep_sasa_burial":      str(result["pep_sasa_burial"])
                                            if result.get("pep_sasa_burial") is not None else "",
                    "mse_converted":        result.get("mse_converted", ""),
                    "caps_stripped":        result.get("caps_stripped", ""),
                    "cap_min_dist_A":       result.get("cap_min_dist_A", ""),
                    "peptide_topology":     result.get("peptide_topology", "linear"),
                    "reject_reason":        "",
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
                print(f"  [{i}/{n_total}] PASS   {pdb_id}  asm={result['assembly_num']}  "
                    f"pep={pep_str}  rec={rec_str}  res={result.get('resolution','NA')}"
                    f"{ions_str}  (accepted={accepted_count} rejected={rejected_count})")
                if result["peptide_length"] != n_mer:
                    print(f"    WARNING: peptide length {result['peptide_length']} != n_mer {n_mer}")
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
                print(f"  [{i}/{n_total}] REJECT {pdb_id}: {result['reject_reason']}  "
                    f"(accepted={accepted_count} rejected={rejected_count})")

            _write_row(row)

    new_count = accepted_count + rejected_count - cached_count - patched_count
    print(f"\nDone. accepted={accepted_count} rejected={rejected_count} "
          f"total={accepted_count + rejected_count} "
          f"(cached={cached_count}  patched={patched_count}  new={new_count})")
    print(f"Written {out_tsv}")

    if args.verbose:
        stdout.write(
            f"[scrape/02-filter] done  —  {accepted_count} accepted"
            f"  {rejected_count} rejected\n"
        )
        stdout.flush()


if __name__ == "__main__":
    main()
