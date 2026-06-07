"""
structural_redundancy.py — Structural redundancy filtering.

Methods:
  filter_tm_score_receptor   Greedy clustering by receptor TM-score (tmtools).
  filter_irmsd_interface     Greedy clustering by binding interface RMSD.
  apply_structural_redundancy_filters  Entry point; enforces TM-score → iRMSD order.

Pipeline position: after sequence redundancy, before sampling.

Notes:
  - TM-score requires: pip install tmtools
  - iRMSD uses positional mapping for residue correspondence (v1 simplification).
  - iRMSD without TM-score pre-filtering produces a warning.
"""

import numpy as np
from pathlib import Path
from Bio.PDB import MMCIFParser

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H", "MSE": "M", "SEC": "C",
}


def _res_sort_key(row):
    try:
        return float(row.get("resolution", "99"))
    except (ValueError, TypeError):
        return 99.0


def _receptor_chains(row):
    """'AB:440' → ['A', 'B']"""
    return list(row["receptor"].split(":")[0])


def _peptide_chain(row):
    """'C:8' → 'C'"""
    return row["peptide"].split(":")[0]


def _load_struct(pdb_id, fixed_cif_dir):
    cif_path = Path(fixed_cif_dir) / f"{pdb_id}.cif"
    if not cif_path.exists():
        return None
    try:
        return MMCIFParser(QUIET=True).get_structure(pdb_id, str(cif_path))
    except Exception:
        return None


def _receptor_ca(struct, chains):
    """Return (ca_coords np.array, seq str) for receptor chains."""
    model = struct[0]
    coords, seq = [], []
    for cid in chains:
        if cid not in model:
            continue
        for res in model[cid].get_residues():
            if "CA" in res:
                coords.append(res["CA"].coord)
                seq.append(AA3_TO_1.get(res.resname, "X"))
    if not coords:
        return None, None
    return np.array(coords, dtype=float), "".join(seq)


def _peptide_all_coords(struct, chain):
    """Return (N, 3) array of all peptide heavy atom coordinates."""
    model = struct[0]
    coords = []
    if chain in model:
        for res in model[chain].get_residues():
            for atom in res.get_atoms():
                coords.append(atom.coord)
    return np.array(coords, dtype=float) if coords else None


def _kabsch(fixed, mobile):
    """Kabsch superimposition: returns (rot, tran) for row-vector form.

    Transforms mobile onto fixed:  mobile_new = mobile @ rot.T + tran
    Both arrays shape (N, 3).
    """
    cf = fixed.mean(axis=0)
    cm = mobile.mean(axis=0)
    F = fixed  - cf
    M = mobile - cm
    H = M.T @ F                        # (3, 3)
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, d])
    rot = Vt.T @ D @ U.T               # maps M columns → F columns
    tran = cf - cm @ rot.T
    return rot, tran


def _interface_mask(ca_coords, pep_coords, dist):
    """Boolean mask: True if CA is within dist Å of any peptide atom."""
    diff = ca_coords[:, np.newaxis, :] - pep_coords[np.newaxis, :, :]
    return np.sqrt((diff ** 2).sum(axis=2)).min(axis=1) < dist


# ---------------------------------------------------------------------------
# TM-score receptor
# ---------------------------------------------------------------------------

def filter_tm_score_receptor(rows, threshold, fixed_cif_dir):
    """Greedy TM-score clustering of receptor structures.

    threshold: TM-score >= threshold → redundant.  Standard: 0.5 = same fold.
    Keeps one representative per cluster (best crystallographic resolution first).
    Requires tmtools: pip install tmtools
    """
    try:
        from tmtools import tm_align
    except ImportError:
        raise ImportError("tmtools not installed. Run: pip install tmtools")

    before = len(rows)
    cache = {}

    def _get(row):
        pid = row["pdb_id"]
        if pid not in cache:
            struct = _load_struct(pid, fixed_cif_dir)
            if struct is None:
                cache[pid] = None
            else:
                coords, seq = _receptor_ca(struct, _receptor_chains(row))
                cache[pid] = (coords, seq) if coords is not None else None
        return cache[pid]

    kept = []
    for row in sorted(rows, key=_res_sort_key):
        data = _get(row)
        if data is None:
            kept.append(row)
            continue
        coords_q, seq_q = data

        redundant = False
        for rep in kept:
            rep_data = _get(rep)
            if rep_data is None:
                continue
            try:
                res = tm_align(coords_q, rep_data[0], seq_q, rep_data[1])
                if max(res.tm_norm_chain1, res.tm_norm_chain2) >= threshold:
                    redundant = True
                    break
            except Exception:
                continue

        if not redundant:
            kept.append(row)

    print(f"  TM-score receptor (threshold={threshold:.2f}): "
          f"{len(kept)} of {before} entries kept ({before - len(kept)} removed)")
    return kept


# ---------------------------------------------------------------------------
# iRMSD interface
# ---------------------------------------------------------------------------

def filter_irmsd_interface(rows, irmsd_threshold, fixed_cif_dir, interface_dist=4.0):
    """Greedy clustering by binding interface RMSD.

    Per-pair steps:
      1. Superimpose receptors globally (Kabsch/SVD).
      2. Find interface residues: CA within interface_dist Å of any peptide atom.
      3. Compute RMSD on common interface CAs (intersection of both interfaces).
      4. Redundant if n_common >= 5 and iRMSD < irmsd_threshold.

    Positional mapping is used (residues matched by array index after sorting).
    Designed to run after TM-score filtering.
    """
    before = len(rows)
    cache = {}

    def _get(row):
        pid = row["pdb_id"]
        if pid not in cache:
            struct = _load_struct(pid, fixed_cif_dir)
            if struct is None:
                cache[pid] = None
            else:
                ca_coords, _ = _receptor_ca(struct, _receptor_chains(row))
                pep_coords   = _peptide_all_coords(struct, _peptide_chain(row))
                cache[pid] = (ca_coords, pep_coords) \
                    if (ca_coords is not None and pep_coords is not None) else None
        return cache[pid]

    def _compute(data_a, data_b):
        ca_a, pep_a = data_a
        ca_b, pep_b = data_b
        n = min(len(ca_a), len(ca_b))
        if n < 10:
            return None, 0

        # Superimpose on first n CA atoms
        rot, tran = _kabsch(ca_a[:n], ca_b[:n])

        # Apply transform to ALL of ca_b and pep_b (no in-place mutation)
        ca_b_T   = ca_b   @ rot.T + tran
        pep_b_T  = pep_b  @ rot.T + tran

        # Interface masks
        mask_a = _interface_mask(ca_a[:n],   pep_a,   interface_dist)
        mask_b = _interface_mask(ca_b_T[:n], pep_b_T, interface_dist)
        common = mask_a & mask_b
        n_common = int(common.sum())

        if n_common < 5:
            return None, n_common

        irmsd = float(np.sqrt(
            ((ca_a[:n][common] - ca_b_T[:n][common]) ** 2).sum(axis=1).mean()
        ))
        return irmsd, n_common

    kept = []
    for row in sorted(rows, key=_res_sort_key):
        data = _get(row)
        if data is None:
            kept.append(row)
            continue

        redundant = False
        for rep in kept:
            rep_data = _get(rep)
            if rep_data is None:
                continue
            irmsd, _ = _compute(data, rep_data)
            if irmsd is not None and irmsd < irmsd_threshold:
                redundant = True
                break

        if not redundant:
            kept.append(row)

    print(f"  iRMSD interface ({interface_dist:.1f}Å def, threshold={irmsd_threshold:.1f}Å): "
          f"{len(kept)} of {before} entries kept ({before - len(kept)} removed)")
    return kept


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def apply_structural_redundancy_filters(rows, args, fixed_cif_dir):
    """Apply structural redundancy filters: TM-score receptor → iRMSD interface."""
    has_tm   = getattr(args, "tm_score_receptor", None) is not None
    has_irms = getattr(args, "irmsd", None) is not None

    if has_irms and not has_tm:
        print("  Warning: --irmsd without --tm-score-receptor may give poor results "
              "for structurally unrelated receptors.")

    if has_tm:
        rows = filter_tm_score_receptor(rows, args.tm_score_receptor, fixed_cif_dir)
    if has_irms:
        rows = filter_irmsd_interface(
            rows, args.irmsd, fixed_cif_dir,
            interface_dist=getattr(args, "irmsd_dist", 4.0),
        )
    return rows
