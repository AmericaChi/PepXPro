"""
06_tm_matrix.py — Pre-compute receptor TM-score matrix for all accepted entries.

Computes max(TM(i→j), TM(j→i)) for every receptor pair and saves a compressed
NumPy matrix that gen_sample loads instead of computing TM-scores on-the-fly.

Input:  outputs/filtered_systems.tsv  (accepted entries)
        fixed_cif/<pdb_id>.cif

Output: outputs/tm_score_matrix.npz   symmetric float32 N×N matrix
        outputs/tm_score_ids.txt       pdb_id per row/column (one per line)
        logs/06_tm_matrix.log

Usage:
  python 06_tm_matrix.py
  python 06_tm_matrix.py --force    # recompute even if matrix already exists
  python 06_tm_matrix.py -v         # verbose
"""

import argparse
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, setup_logging

OUTPUTS_DIR  = _root / "outputs"
FIXED_DIR    = _root / "fixed_cif"
FILTERED_TSV = OUTPUTS_DIR / "filtered_systems.tsv"
MATRIX_NPZ   = OUTPUTS_DIR / "tm_score_matrix.npz"
MATRIX_IDS   = OUTPUTS_DIR / "tm_score_ids.txt"
LOG_FILE     = _root / "logs" / "06_tm_matrix.log"

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H", "MSE": "M", "SEC": "C",
}

# Global storage for worker processes — avoids pickling structs on every call
_WORKER_STRUCTS: list = []


def _receptor_chain_ids(row):
    return [c.strip() for c in row["receptor_chains"].split(",")]


def _load_receptor(pdb_id, chain_ids):
    cif_path = FIXED_DIR / f"{pdb_id}.cif"
    if not cif_path.exists():
        return None
    try:
        struct = MMCIFParser(QUIET=True).get_structure(pdb_id, str(cif_path))
    except Exception:
        return None
    model = struct[0]
    coords, seq = [], []
    for cid in chain_ids:
        if cid not in model:
            continue
        for res in model[cid].get_residues():
            if "CA" in res:
                coords.append(res["CA"].coord)
                seq.append(AA3_TO_1.get(res.resname, "X"))
    if not coords:
        return None
    return np.array(coords, dtype=np.float64), "".join(seq)


def _worker_init(structs):
    """Load all structures into each worker process once at startup."""
    global _WORKER_STRUCTS
    _WORKER_STRUCTS = structs


def _compute_pair(ij):
    """Compute max(TM(i→j), TM(j→i)) for a single pair. Uses global structs."""
    i, j = ij
    di = _WORKER_STRUCTS[i]
    dj = _WORKER_STRUCTS[j]
    if di is None or dj is None:
        return i, j, 0.0
    try:
        from tmtools import tm_align
        res = tm_align(di[0], dj[0], di[1], dj[1])
        return i, j, float(max(res.tm_norm_chain1, res.tm_norm_chain2))
    except Exception:
        return i, j, 0.0


def build_matrix(structs, n_workers):
    """Compute symmetric N×N TM-score matrix using ProcessPoolExecutor."""
    n = len(structs)
    matrix = np.zeros((n, n), dtype=np.float32)
    np.fill_diagonal(matrix, 1.0)

    pairs = [(i, j) for i in range(n) for j in range(i)]
    total = len(pairs)

    print(f"  Submitting {total:,} pairs to {n_workers} processes ...")
    t0 = time.perf_counter()
    done = 0

    with mp.Pool(processes=n_workers, initializer=_worker_init,
                 initargs=(structs,)) as pool:
        for i, j, val in pool.imap_unordered(_compute_pair, pairs, chunksize=64):
            matrix[i, j] = val
            matrix[j, i] = val
            done += 1
            if done % 500 == 0 or done == total:
                elapsed = time.perf_counter() - t0
                pct     = done / total * 100
                eta_s   = (elapsed / done * (total - done)) if done else 0
                print(f"\r  {done:>6}/{total} pairs  {pct:5.1f}%  "
                      f"elapsed {int(elapsed//60)}m{int(elapsed%60):02d}s  "
                      f"ETA {int(eta_s//60)}m{int(eta_s%60):02d}s   ",
                      end="", flush=True)

    print()
    return matrix


def main():
    parser = argparse.ArgumentParser(
        description="Pre-compute receptor TM-score matrix for all accepted entries.")
    parser.add_argument("--force", action="store_true",
                        help="Recompute even if matrix already exists.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print brief start/done line to stdout.")
    args = parser.parse_args()

    _logger, stdout = setup_logging(LOG_FILE, args.verbose)
    if args.verbose:
        stdout.write("[scrape/06-tm-matrix] starting ...\n")
        stdout.flush()

    if MATRIX_NPZ.exists() and MATRIX_IDS.exists() and not args.force:
        ids = MATRIX_IDS.read_text().splitlines()
        print(f"Matrix already exists ({len(ids)} entries). Use --force to recompute.")
        return

    try:
        from tmtools import tm_align  # noqa: F401
    except ImportError:
        sys.stderr.write(
            "Error: tmtools not installed.\n"
            "Install with: python -m pip install tmtools\n"
        )
        sys.exit(1)

    rows = [r for r in load_tsv(FILTERED_TSV) if r.get("accepted") == "Y"]
    print(f"Accepted entries: {len(rows)}")

    print("Loading receptor structures ...")
    t_load = time.perf_counter()
    ids, structs = [], []
    for row in rows:
        data = _load_receptor(row["pdb_id"], _receptor_chain_ids(row))
        ids.append(row["pdb_id"])
        structs.append(data)
    n_loaded = sum(1 for s in structs if s is not None)
    print(f"  {n_loaded}/{len(ids)} CIFs loaded  "
          f"({len(ids) - n_loaded} missing)  "
          f"({time.perf_counter() - t_load:.1f}s)")

    n_workers = min(10, os.cpu_count() or 4)
    t_mat = time.perf_counter()
    matrix = build_matrix(structs, n_workers)
    elapsed = time.perf_counter() - t_mat
    mins, secs = divmod(int(elapsed), 60)
    print(f"Matrix built in {mins}m {secs:02d}s")

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(MATRIX_NPZ, matrix=matrix)
    MATRIX_IDS.write_text("\n".join(ids))
    size_kb = MATRIX_NPZ.stat().st_size / 1024
    print(f"Saved {MATRIX_NPZ.name} ({size_kb:.0f} KB)  +  {MATRIX_IDS.name}")
    print(f"  {len(ids)} entries  |  {matrix.shape}  |  float32")

    if args.verbose:
        stdout.write(
            f"[scrape/06-tm-matrix] done  —  {len(ids)} entries  {mins}m {secs:02d}s\n"
        )
        stdout.flush()


if __name__ == "__main__":
    # macOS Python 3.12+: default is spawn (safe with ObjC frameworks).
    # Linux: default is fork (faster). Both work with the initializer pattern.
    if sys.platform == "darwin":
        mp.set_start_method("spawn", force=True)
    main()
