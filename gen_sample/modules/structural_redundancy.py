import concurrent.futures
import os

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


def filter_tm_score_receptor(rows, threshold, fixed_cif_dir):
    """Greedy TM-score clustering of receptor structures.

    threshold: TM-score >= threshold → redundant.  Standard: 0.5 = same fold.
    Keeps one representative per cluster (best crystallographic resolution first).
    Requires tmtools: pip install tmtools

    Inner loop is parallelized across CPU cores for speed.
    """
    try:
        from tmtools import tm_align
    except ImportError:
        raise ImportError("tmtools not installed. Run: pip install tmtools")

    before = len(rows)
    cache = {}

    def _load(row):
        pid = row["pdb_id"]
        if pid not in cache:
            struct = _load_struct(pid, fixed_cif_dir)
            if struct is None:
                cache[pid] = None
            else:
                coords, seq = _receptor_ca(struct, _receptor_chains(row))
                cache[pid] = (coords, seq) if coords is not None else None
        return cache[pid]

    sorted_rows = sorted(rows, key=_res_sort_key)

    # Pre-load all CIFs once (IO-bound, sequential with cache)
    print(f"  TM-score: loading {before} structures ...", end="", flush=True)
    for row in sorted_rows:
        _load(row)
    n_cif = sum(1 for r in sorted_rows if cache.get(r["pdb_id"]) is not None)
    print(f" {n_cif} CIFs found.")

    def _is_redundant(coords_q, seq_q, coords_r, seq_r):
        try:
            res = tm_align(coords_q, coords_r, seq_q, seq_r)
            return max(res.tm_norm_chain1, res.tm_norm_chain2) >= threshold
        except Exception:
            return False

    n_workers = min(8, os.cpu_count() or 4)
    kept = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        for i, row in enumerate(sorted_rows):
            data = cache.get(row["pdb_id"])
            if data is None:
                kept.append(row)
            else:
                coords_q, seq_q = data
                rep_datas = [cache[rep["pdb_id"]] for rep in kept
                             if cache.get(rep["pdb_id"]) is not None]

                redundant = False
                if rep_datas:
                    futs = [
                        executor.submit(_is_redundant, coords_q, seq_q, rd[0], rd[1])
                        for rd in rep_datas
                    ]
                    for fut in concurrent.futures.as_completed(futs):
                        if fut.result():
                            redundant = True
                            for f in futs:
                                f.cancel()
                            break
                    concurrent.futures.wait(futs)

                if not redundant:
                    kept.append(row)

            print(f"\r  TM-score: {i+1}/{before} checked, {len(kept)} kept ...",
                  end="", flush=True)

    print()
    print(f"  TM-score receptor (threshold={threshold:.2f}): "
          f"{len(kept)} of {before} entries kept ({before - len(kept)} removed)")
    return kept


def apply_structural_redundancy_filters(rows, args, fixed_cif_dir):
    if getattr(args, "tm_score_receptor", None) is not None:
        rows = filter_tm_score_receptor(rows, args.tm_score_receptor, fixed_cif_dir)
    return rows
