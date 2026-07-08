"""
06_tm_matrix.py — Pre-compute receptor TM-score matrix for all accepted entries.

Computes max(TM(i->j), TM(j->i)) for every receptor pair and saves a compressed
NumPy matrix that genSample loads instead of computing TM-scores on-the-fly.

Reads all filtered_systems*.tsv files found in outputs/ so external datasets
(PPIKB, PEPBI, ...) are included automatically.

Incremental mode (default): if a matrix already exists, only new entries are
computed and the matrix is expanded. Use --force to recompute from scratch.

Input:  outputs/filtered_systems*.tsv  (accepted entries)
        fixed_cif/<pdb_id>.cif

Output: outputs/tm_score_matrix.npz   symmetric float32 N x N matrix
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
    global _WORKER_STRUCTS
    _WORKER_STRUCTS = structs


def _compute_pair(ij):
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


def _compute_pairs(pairs, structs, n_workers, label=""):
    total = len(pairs)
    results = []
    if total == 0:
        return results
    print(f"  Submitting {total:,} pairs to {n_workers} processes {label}...")
    t0 = time.perf_counter()
    done = 0
    with mp.Pool(processes=n_workers, initializer=_worker_init,
                 initargs=(structs,)) as pool:
        for i, j, val in pool.imap_unordered(_compute_pair, pairs, chunksize=64):
            results.append((i, j, val))
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
    return results


def _load_all_rows():
    """Load accepted rows from all filtered_systems*.tsv files, deduplicating by pdb_id."""
    skip = {"filtered_systems_merged.tsv"}
    tsvs = sorted(OUTPUTS_DIR.glob("filtered_systems*.tsv"))
    tsvs = [p for p in tsvs if p.name not in skip]
    seen, rows = set(), []
    for tsv in tsvs:
        for r in load_tsv(tsv):
            if r.get("accepted") == "Y" and r["pdb_id"] not in seen:
                seen.add(r["pdb_id"])
                rows.append(r)
    sources = ", ".join(p.name for p in tsvs)
    print(f"  Sources: {sources}")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Pre-compute receptor TM-score matrix for all accepted entries.")
    parser.add_argument("--force", action="store_true",
                        help="Recompute from scratch even if matrix already exists.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print brief start/done line to stdout.")
    args = parser.parse_args()

    _logger, stdout = setup_logging(LOG_FILE, args.verbose)
    if args.verbose:
        stdout.write("[scrape/06-tm-matrix] starting ...\n")
        stdout.flush()

    try:
        from tmtools import tm_align  # noqa: F401
    except ImportError:
        sys.stderr.write(
            "Error: tmtools not installed.\n"
            "Install with: python -m pip install tmtools\n"
        )
        sys.exit(1)

    rows = _load_all_rows()
    print(f"Accepted entries: {len(rows)}")

    # ── Load existing matrix for incremental update ───────────────────────────
    existing_ids = []
    existing_matrix = None
    if MATRIX_NPZ.exists() and MATRIX_IDS.exists() and not args.force:
        existing_ids = MATRIX_IDS.read_text().splitlines()
        existing_matrix = np.load(MATRIX_NPZ)["matrix"]
        print(f"Existing matrix: {len(existing_ids)} entries  {existing_matrix.shape}")

    all_ids = [r["pdb_id"] for r in rows]
    existing_set = set(existing_ids)
    new_rows = [r for r in rows if r["pdb_id"] not in existing_set]

    if existing_matrix is not None and not new_rows:
        print("Matrix is up to date — no new entries. Use --force to recompute.")
        return

    if new_rows:
        print(f"New entries to add: {len(new_rows)}")
    else:
        print("Building matrix from scratch ...")

    # ── Load all receptor structures ──────────────────────────────────────────
    print("Loading receptor structures ...")
    t_load = time.perf_counter()
    id_to_struct = {}
    load_rows = rows if existing_matrix is None else new_rows
    for row in load_rows:
        data = _load_receptor(row["pdb_id"], _receptor_chain_ids(row))
        id_to_struct[row["pdb_id"]] = data
    n_loaded = sum(1 for s in id_to_struct.values() if s is not None)
    print(f"  {n_loaded}/{len(id_to_struct)} CIFs loaded  "
          f"({time.perf_counter() - t_load:.1f}s)")

    n_workers = min(10, os.cpu_count() or 4)
    t_mat = time.perf_counter()

    if existing_matrix is None or args.force:
        # ── Full build ────────────────────────────────────────────────────────
        structs = [id_to_struct[pid] for pid in all_ids]
        n = len(structs)
        matrix = np.zeros((n, n), dtype=np.float32)
        np.fill_diagonal(matrix, 1.0)
        pairs = [(i, j) for i in range(n) for j in range(i)]
        for i, j, val in _compute_pairs(pairs, structs, n_workers):
            matrix[i, j] = val
            matrix[j, i] = val
    else:
        # ── Incremental expand ────────────────────────────────────────────────
        n_old = len(existing_ids)
        n_new = len(new_rows)
        n_total = n_old + n_new
        new_ids = [r["pdb_id"] for r in new_rows]

        # Also load existing structs needed for cross-pairs
        print(f"Loading existing receptor structures for cross-pairs ...")
        t2 = time.perf_counter()
        for row in rows:
            if row["pdb_id"] in existing_set and row["pdb_id"] not in id_to_struct:
                id_to_struct[row["pdb_id"]] = _load_receptor(
                    row["pdb_id"], _receptor_chain_ids(row))
        print(f"  Done ({time.perf_counter() - t2:.1f}s)")

        all_structs = [id_to_struct.get(pid) for pid in existing_ids + new_ids]

        # Expand matrix: copy old block, fill new rows/cols
        matrix = np.zeros((n_total, n_total), dtype=np.float32)
        matrix[:n_old, :n_old] = existing_matrix
        np.fill_diagonal(matrix, 1.0)

        # Pairs: new vs all (old + new, lower triangle only)
        pairs = []
        for ni, _ in enumerate(new_ids):
            i = n_old + ni
            for j in range(i):
                pairs.append((i, j))

        for i, j, val in _compute_pairs(pairs, all_structs, n_workers,
                                        label=f"(+{n_new} new vs {n_total-1} total)"):
            matrix[i, j] = val
            matrix[j, i] = val

        all_ids = existing_ids + new_ids

    elapsed = time.perf_counter() - t_mat
    mins, secs = divmod(int(elapsed), 60)
    print(f"Matrix built in {mins}m {secs:02d}s")

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(MATRIX_NPZ, matrix=matrix)
    MATRIX_IDS.write_text("\n".join(all_ids))
    size_kb = MATRIX_NPZ.stat().st_size / 1024
    print(f"Saved {MATRIX_NPZ.name} ({size_kb:.0f} KB)  +  {MATRIX_IDS.name}")
    print(f"  {len(all_ids)} entries  |  {matrix.shape}  |  float32")

    if args.verbose:
        stdout.write(
            f"[scrape/06-tm-matrix] done  --  {len(all_ids)} entries  {mins}m {secs:02d}s\n"
        )
        stdout.flush()


if __name__ == "__main__":
    if sys.platform == "darwin":
        mp.set_start_method("spawn", force=True)
    main()
