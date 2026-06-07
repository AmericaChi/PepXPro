from pathlib import Path
from .utils import load_tsv


def _report(label, before, after):
    print(f"  {label}: {after} of {before} entries kept")


def filter_single_chain(rows):
    before = len(rows)
    rows = [r for r in rows if r["_n_rec_chains"] == 1]
    _report("--single-chain", before, len(rows))
    return rows


def filter_min_receptor(rows, n):
    before = len(rows)
    rows = [r for r in rows if r["_rec_len"] >= n]
    _report(f"--min-receptor {n}", before, len(rows))
    return rows


def filter_max_receptor(rows, n):
    before = len(rows)
    rows = [r for r in rows if r["_rec_len"] <= n]
    _report(f"--max-receptor {n}", before, len(rows))
    return rows


def filter_min_peptide(rows, n):
    before = len(rows)
    rows = [r for r in rows if r["_pep_len"] >= n]
    _report(f"--min-peptide {n}", before, len(rows))
    return rows


def filter_max_peptide(rows, n):
    before = len(rows)
    rows = [r for r in rows if r["_pep_len"] <= n]
    _report(f"--max-peptide {n}", before, len(rows))
    return rows


def filter_min_dg(rows, val):
    before = len(rows)
    rows = [r for r in rows if r["_dg"] >= val]
    _report(f"--min-dg {val}", before, len(rows))
    return rows


def filter_max_dg(rows, val):
    before = len(rows)
    rows = [r for r in rows if r["_dg"] <= val]
    _report(f"--max-dg {val}", before, len(rows))
    return rows


def filter_binding_type(rows, btype):
    """btype: 'kd', 'ki', or 'both'"""
    if btype.lower() == "both":
        return rows
    before = len(rows)
    rows = [r for r in rows if r.get("binding_type", "").lower() == btype.lower()]
    _report(f"--binding-type {btype}", before, len(rows))
    return rows


def filter_max_added(rows, n, filtered_tsv_path):
    filtered_tsv_path = Path(filtered_tsv_path)
    if not filtered_tsv_path.exists():
        raise FileNotFoundError(f"--max-added requires {filtered_tsv_path}")
    pep_added = {}
    for row in load_tsv(filtered_tsv_path):
        if row.get("accepted") == "Y" and row.get("peptide_added"):
            try:
                pep_added[row["pdb_id"]] = int(row["peptide_added"])
            except ValueError:
                pass
    if not pep_added:
        raise ValueError("peptide_added column missing from filtered_systems.tsv")
    before = len(rows)
    rows = [r for r in rows if pep_added.get(r["pdb_id"], 0) <= n]
    _report(f"--max-added {n}", before, len(rows))
    return rows


def filter_max_receptor_added(rows, n, filtered_tsv_path):
    filtered_tsv_path = Path(filtered_tsv_path)
    if not filtered_tsv_path.exists():
        raise FileNotFoundError(f"--max-receptor-added requires {filtered_tsv_path}")
    rec_added = {}
    for row in load_tsv(filtered_tsv_path):
        if row.get("accepted") == "Y" and row.get("receptor_added"):
            try:
                rec_added[row["pdb_id"]] = int(row["receptor_added"])
            except ValueError:
                pass
    if not rec_added:
        raise ValueError("receptor_added column missing from filtered_systems.tsv")
    before = len(rows)
    rows = [r for r in rows if rec_added.get(r["pdb_id"], 0) <= n]
    _report(f"--max-receptor-added {n}", before, len(rows))
    return rows


def filter_max_resolution(rows, val):
    """Keep only entries with resolution <= val Å. NMR/NA entries are excluded."""
    before = len(rows)
    rows = [r for r in rows if r["_resolution"] is not None and r["_resolution"] <= val]
    _report(f"--max-resolution {val}", before, len(rows))
    return rows


def filter_min_resolution(rows, val):
    """Keep only entries with resolution >= val Å. NMR/NA entries are excluded."""
    before = len(rows)
    rows = [r for r in rows if r["_resolution"] is not None and r["_resolution"] >= val]
    _report(f"--min-resolution {val}", before, len(rows))
    return rows


def filter_exclude_nmr(rows):
    """Remove entries with no numeric resolution (NMR, NA)."""
    before = len(rows)
    rows = [r for r in rows if r["_resolution"] is not None]
    _report("--exclude-nmr", before, len(rows))
    return rows


def filter_exclude(rows, exclude_path):
    exclude_path = Path(exclude_path)
    if not exclude_path.exists():
        raise FileNotFoundError(f"--exclude file not found: {exclude_path}")
    with open(exclude_path, encoding="ascii") as fh:
        excluded = {line.strip() for line in fh if line.strip()}
    before = len(rows)
    rows = [r for r in rows if r["pdb_id"] not in excluded]
    _report(f"--exclude ({len(excluded)} ids)", before, len(rows))
    return rows


def apply_quality_filters(rows, args, filtered_tsv_path):
    """Apply all active filters from parsed argparse args."""
    if getattr(args, "single_chain", False):
        rows = filter_single_chain(rows)
    if getattr(args, "min_receptor", None) is not None:
        rows = filter_min_receptor(rows, args.min_receptor)
    if getattr(args, "max_receptor", None) is not None:
        rows = filter_max_receptor(rows, args.max_receptor)
    if getattr(args, "min_peptide", None) is not None:
        rows = filter_min_peptide(rows, args.min_peptide)
    if getattr(args, "max_peptide", None) is not None:
        rows = filter_max_peptide(rows, args.max_peptide)
    if getattr(args, "min_dg", None) is not None:
        rows = filter_min_dg(rows, args.min_dg)
    if getattr(args, "max_dg", None) is not None:
        rows = filter_max_dg(rows, args.max_dg)
    if getattr(args, "binding_type", "both") != "both":
        rows = filter_binding_type(rows, args.binding_type)
    if getattr(args, "max_added", None) is not None:
        rows = filter_max_added(rows, args.max_added, filtered_tsv_path)
    if getattr(args, "max_receptor_added", None) is not None:
        rows = filter_max_receptor_added(rows, args.max_receptor_added, filtered_tsv_path)
    if getattr(args, "max_resolution", None) is not None:
        rows = filter_max_resolution(rows, args.max_resolution)
    if getattr(args, "min_resolution", None) is not None:
        rows = filter_min_resolution(rows, args.min_resolution)
    if getattr(args, "exclude_nmr", False):
        rows = filter_exclude_nmr(rows)
    if getattr(args, "exclude", None) is not None:
        rows = filter_exclude(rows, args.exclude)
    return rows
