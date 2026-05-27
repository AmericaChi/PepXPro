"""
Step 4: Browse fixed_cif/*.cif structures in PyMOL one at a time.
Chains are coloured with util.cbc. Navigate with left/right arrow keys.

Flags:
  --sample N         Select N structures from summary.tsv stratified across ΔG
                     and peptide length (6 quantile bins each, 6×6 = 36 cells),
                     then iteratively replaces the most receptor-similar pair
                     (CA RMSD via PyMOL align) until all comparable pairs exceed
                     --receptor-rmsd. Prints the selection to screen and writes
                     sample.tsv (or the file given by --output). Does not open PyMOL.
                     Always prints a pairwise receptor RMSD cluster report at
                     0.2/0.5/1/2/3/5 Å thresholds.

  --expand N         Add N more structures to an existing sample.tsv (or the file
                     given by --output), chosen to be distributed across ΔG and
                     peptide length and structurally distinct from both each other
                     and the existing sample entries (CA RMSD via PyMOL align).
                     Overwrites the output file. Always prints a pairwise receptor
                     RMSD cluster report (new entries only) at 0.2/0.5/1/2/3/5 Å thresholds.

  --output FILE      TSV file to read/write for --sample, --expand, and
                     --view-sample. Defaults to sample.tsv. Useful for keeping
                     multiple independent samples (e.g. sample_sc.tsv).

  --receptor-rmsd F  RMSD threshold for the receptor diversity filter used with
                     --sample and --expand (default 0, disabled). Set > 0 to enable.

  --single-chain     With --sample/--expand: restrict pool to entries with a
                     single receptor chain before sampling.

  --max-receptor N   With --sample/--expand: restrict pool to entries whose
                     receptor residue count is ≤ N before sampling.

  --max-peptide N    With --sample/--expand: restrict pool to entries whose
                     peptide residue count is ≤ N before sampling.

  --max-added N      With --sample/--expand: restrict pool to entries where
                     PDBFixer added ≤ N residues to the peptide chain.
                     Reads peptide_added from filtered_systems.tsv directly
                     (step 2 must have been run with the current version).

  --view-sample      Open the PyMOL viewer restricted to the systems in sample.tsv
                     (or --output FILE), ordered by ΔG (tightest to weakest binder).

  --overwrite        Allow --sample to overwrite an existing output file. Without
                     this flag, --sample aborts if the output file already exists.
                     Not needed for --expand, which always rewrites its input file.

  --clean            Remove the output file (sample.tsv or --output FILE).

  (no flags)         Browse all structures in fixed_cif/.

Run: python 04_sample_and_view.py --sample 20
     python 04_sample_and_view.py --sample 20 --output sample_sc.tsv
     python 04_sample_and_view.py --sample 20 --receptor-rmsd 2.0
     python 04_sample_and_view.py --sample 20 --single-chain
     python 04_sample_and_view.py --sample 20 --max-receptor 200
     python 04_sample_and_view.py --sample 20 --max-peptide 15
     python 04_sample_and_view.py --sample 20 --max-added 3
     python 04_sample_and_view.py --expand 10
     python 04_sample_and_view.py --expand 10 --output sample_sc.tsv
     python 04_sample_and_view.py --expand 10 --single-chain
     python 04_sample_and_view.py --expand 10 --max-added 3
     python 04_sample_and_view.py --view-sample
     python 04_sample_and_view.py --view-sample --output sample_sc.tsv
     python 04_sample_and_view.py
"""

import contextlib
import os
import random
from pathlib import Path
import sys


@contextlib.contextmanager
def _suppress_fd_output():
    """Redirect C-level stdout/stderr to /dev/null (silences PyMOL banner/warnings)."""
    devnull = os.open(os.devnull, os.O_WRONLY)
    saved = [os.dup(1), os.dup(2)]
    os.dup2(devnull, 1)
    os.dup2(devnull, 2)
    try:
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        for fd in saved + [devnull]:
            os.close(fd)

FIXED_DIR    = Path("fixed_cif")
SUMMARY_TSV  = Path("summary.tsv")
SAMPLE_TSV   = Path("sample.tsv")
SUMMARY_COLS = ["pdb_id", "receptor", "peptide", "dG_kcal_mol", "binding_type"]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def load_tsv(path):
    rows = []
    with open(path, encoding="ascii") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.strip().split("\t"))))
    return rows


def diverse_sample(rows, n, dg_bins=6, pep_bins=6):
    """Return n rows with diverse coverage across ΔG and peptide length.

    ΔG and peptide length are each divided into dg_bins/pep_bins quantile-based
    bins (default 6×6 = 36 cells). Bins are equal-count so skewed distributions
    don't collapse into one bin. One random entry is drawn from each occupied cell
    in random order until n entries are selected; remaining slots are filled
    randomly from what's left.
    """
    from collections import defaultdict

    n = min(n, len(rows))
    if n == len(rows):
        return list(rows)

    def quantile_bins(vals, k):
        """Assign each value to a bin 0..k-1 based on its rank."""
        order = sorted(range(len(vals)), key=lambda i: vals[i])
        bins = [0] * len(vals)
        for rank, idx in enumerate(order):
            bins[idx] = min(int(rank * k / len(vals)), k - 1)
        return bins

    dg_b  = quantile_bins([r["_dg"]      for r in rows], dg_bins)
    pep_b = quantile_bins([r["_pep_len"] for r in rows], pep_bins)

    cells = defaultdict(list)
    for i in range(len(rows)):
        cells[(dg_b[i], pep_b[i])].append(i)

    # Draw one random entry per cell (cells visited in random order).
    cell_keys = list(cells.keys())
    random.shuffle(cell_keys)
    selected_idx = []
    for cell in cell_keys:
        if len(selected_idx) >= n:
            break
        selected_idx.append(random.choice(cells[cell]))

    # Fill any remaining slots from unselected entries.
    if len(selected_idx) < n:
        chosen = set(selected_idx)
        remaining = [i for i in range(len(rows)) if i not in chosen]
        random.shuffle(remaining)
        selected_idx.extend(remaining[:n - len(selected_idx)])

    return [rows[i] for i in selected_idx]


# ---------------------------------------------------------------------------
# Receptor RMSD helpers
# ---------------------------------------------------------------------------

def _redraw_similar_receptors(selected, pool, threshold, fixed_dir, max_iters=100,
                              anchors=None):
    """Iteratively replace the most receptor-similar pair (lowest PyMOL align
    RMSD) until all pairs exceed threshold or the pool is exhausted.

    anchors: optional list of fixed rows that participate in RMSD checks but
    are never evicted (used by --expand to protect the existing sample).

    Returns (updated_selected, n_replaced, min_final_rmsd_or_None).
    """
    import pymol2

    anchors = list(anchors) if anchors else []

    with _suppress_fd_output(), pymol2.PyMOL(["-c", "-Q"]) as p:
        cmd = p.cmd
        loaded = set()

        def load_cif(r):
            pid = r["pdb_id"]
            if pid not in loaded:
                cif = fixed_dir / f"{pid}.cif"
                if cif.exists():
                    cmd.load(str(cif), pid)
                    loaded.add(pid)

        def rec_sel(r):
            chains = "+".join(r["receptor"].split(":")[0])
            return f"{r['pdb_id']} and chain {chains} and name CA"

        def rmsd_pair(r1, r2):
            if r1["pdb_id"] not in loaded or r2["pdb_id"] not in loaded:
                return None
            n1, n2 = r1["_rec_len"], r2["_rec_len"]
            if min(n1, n2) / max(n1, n2) < 0.8:
                return None
            try:
                result = cmd.align(rec_sel(r1), rec_sel(r2), cycles=0, quiet=1)
                return result[0] if result[1] > 0 else None
            except Exception:
                return None

        selected = list(selected)
        for r in anchors + selected:
            load_cif(r)

        # Build initial RMSD cache: all selected-selected pairs plus
        # all selected-anchor pairs (anchor-anchor pairs are never needed).
        cache = {}
        for i in range(len(selected)):
            for j in range(i + 1, len(selected)):
                k = tuple(sorted([selected[i]["pdb_id"], selected[j]["pdb_id"]]))
                cache[k] = rmsd_pair(selected[i], selected[j])
            for a in anchors:
                k = tuple(sorted([selected[i]["pdb_id"], a["pdb_id"]]))
                cache[k] = rmsd_pair(selected[i], a)

        selected_ids = {r["pdb_id"] for r in selected}
        available = [r for r in pool if r["pdb_id"] not in selected_ids]

        n_replaced = 0
        for _ in range(max_iters):
            # Find the closest pair that involves at least one selected entry.
            # best_j == -1 means the partner is an anchor (selected[best_i] must be evicted).
            best_i = best_j = -1
            best_rmsd = threshold
            for i in range(len(selected)):
                for j in range(i + 1, len(selected)):
                    k = tuple(sorted([selected[i]["pdb_id"], selected[j]["pdb_id"]]))
                    r = cache.get(k)
                    if r is not None and r < best_rmsd:
                        best_rmsd, best_i, best_j = r, i, j
                for a in anchors:
                    k = tuple(sorted([selected[i]["pdb_id"], a["pdb_id"]]))
                    r = cache.get(k)
                    if r is not None and r < best_rmsd:
                        best_rmsd, best_i, best_j = r, i, -1

            if best_i == -1 or not available:
                break

            # If paired with an anchor we must evict best_i; otherwise pick randomly.
            evict_idx = best_i if best_j == -1 else random.choice([best_i, best_j])
            evicted = selected[evict_idx]
            random.shuffle(available)
            replacement = available.pop(0)
            load_cif(replacement)
            selected[evict_idx] = replacement
            available.append(evicted)
            selected_ids.add(replacement["pdb_id"])
            evicted_pid = evicted["pdb_id"]
            selected_ids.discard(evicted_pid)
            n_replaced += 1

            # Drop stale cache entries; compute only the new pairs.
            cache = {k: v for k, v in cache.items() if evicted_pid not in k}
            for m, other in enumerate(selected):
                if m == evict_idx:
                    continue
                k = tuple(sorted([replacement["pdb_id"], other["pdb_id"]]))
                cache[k] = rmsd_pair(replacement, other)
            for a in anchors:
                k = tuple(sorted([replacement["pdb_id"], a["pdb_id"]]))
                cache[k] = rmsd_pair(replacement, a)

        final_min = min((v for v in cache.values() if v is not None), default=None)

    return selected, n_replaced, final_min, cache


# ---------------------------------------------------------------------------
# RMSD matrix + cluster report
# ---------------------------------------------------------------------------

def compute_rmsd_matrix(selected, fixed_dir):
    """Compute all pairwise receptor CA-RMSD values for a list of rows.

    Returns a dict {(pdb_id_a, pdb_id_b): float|None} where None means PyMOL
    align found no matched CA atoms (incompatible receptor sizes).
    """
    import pymol2

    n = len(selected)
    n_pairs = n * (n - 1) // 2

    def rec_sel(r):
        chains = "+".join(r["receptor"].split(":")[0])
        return f"{r['pdb_id']} and chain {chains} and name CA"

    print(f"  Computing RMSD matrix ({n_pairs} pairs across {n} structures)...", flush=True)

    # Duplicate stdout before the suppressor redirects fd 1, so progress
    # writes survive inside the suppressed PyMOL session.
    _real_out = os.dup(1)

    def _progress(msg):
        os.write(_real_out, (msg + "\n").encode("ascii"))

    cache = {}
    done = 0
    with _suppress_fd_output(), pymol2.PyMOL(["-c", "-Q"]) as p:
        cmd = p.cmd
        for r in selected:
            cif = fixed_dir / f"{r['pdb_id']}.cif"
            if cif.exists():
                cmd.load(str(cif), r["pdb_id"])
        for i in range(n):
            for j in range(i + 1, n):
                key = tuple(sorted([selected[i]["pdb_id"], selected[j]["pdb_id"]]))
                n1 = selected[i]["_rec_len"]
                n2 = selected[j]["_rec_len"]
                if min(n1, n2) / max(n1, n2) < 0.8:
                    cache[key] = None
                else:
                    try:
                        result = cmd.align(rec_sel(selected[i]), rec_sel(selected[j]),
                                           cycles=0, quiet=1)
                        cache[key] = result[0] if result[1] > 0 else None
                    except Exception:
                        cache[key] = None
                done += 1
            _progress(f"  {done}/{n_pairs} pairs  (structure {i + 1}/{n} done)")

    os.close(_real_out)
    return cache


def print_rmsd_cluster_report(rmsd_cache, selected):
    """Print a pairwise RMSD cluster summary at several thresholds."""
    pids = [r["pdb_id"] for r in selected]
    n = len(pids)
    n_pairs = n * (n - 1) // 2

    valid = [v for v in rmsd_cache.values() if v is not None]
    n_skip = n_pairs - len(valid)

    print(f"  receptor RMSD matrix  ({n_pairs} pairs, {n} structures)")
    if not valid:
        print("  no comparable pairs (all receptor sizes incompatible)\n")
        return

    print(f"  all valid pairs: min {min(valid):.2f} Å  "
          f"mean {sum(valid)/len(valid):.2f} Å  "
          f"max {max(valid):.2f} Å"
          + (f"  ({len(valid)} comparable, {n_skip} skipped)" if n_skip else ""))
    print()

    def connected_components(pids, rmsd_cache, threshold):
        """Union-find connected components at given RMSD threshold."""
        parent = {pid: pid for pid in pids}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for (pa, pb), v in rmsd_cache.items():
            if v is not None and v <= threshold:
                ra, rb = find(pa), find(pb)
                if ra != rb:
                    parent[ra] = rb

        clusters = {}
        for pid in pids:
            root = find(pid)
            clusters.setdefault(root, []).append(pid)
        return list(clusters.values())

    print(f"  {'threshold':>10}   {'clusters':>8}   {'singletons':>10}   sizes")
    print(f"  {'─────────':>10}   {'────────':>8}   {'──────────':>10}   ─────")
    for thr in (0.2, 0.5, 1.0, 2.0, 3.0, 5.0):
        clusters = connected_components(pids, rmsd_cache, thr)
        sizes = sorted((len(c) for c in clusters), reverse=True)
        singletons = sum(1 for s in sizes if s == 1)
        multi = [str(s) for s in sizes if s > 1]
        sizes_str = ",".join(multi) if multi else "—"
        print(f"  {thr:>8.1f} Å   {len(clusters):>8}   {singletons:>10}   {sizes_str}")
    print()


# ---------------------------------------------------------------------------
# --sample N
# ---------------------------------------------------------------------------

def run_sample(n, single_chain=False, rmsd_threshold=3.0, max_receptor=None,
               max_peptide=None, max_added=None, overwrite=False, output_path=None):
    out = Path(output_path) if output_path else SAMPLE_TSV
    if out.exists() and not overwrite:
        print(f"  {out} already exists — use --overwrite to replace it.")
        return
    if not SUMMARY_TSV.exists():
        print(f"  {SUMMARY_TSV} not found — run step 3 first.")
        return

    rows = []
    for row in load_tsv(SUMMARY_TSV):
        try:
            row["_dg"]      = float(row.get("dG_kcal_mol", "NA"))
            row["_pep_len"] = int(row.get("peptide", ":0").split(":")[1])
            row["_rec_len"] = int(row.get("receptor", ":0").split(":")[1])
        except (ValueError, IndexError):
            continue
        rows.append(row)

    if not rows:
        print("  No entries with valid ΔG values in summary.tsv.")
        return

    if max_added is not None:
        filtered_tsv = Path("filtered_systems.tsv")
        if not filtered_tsv.exists():
            print(f"  --max-added requires filtered_systems.tsv — run step 2 first.")
            return
        pep_added_map = {}
        for frow in load_tsv(filtered_tsv):
            if frow.get("accepted") == "Y" and frow.get("peptide_added", ""):
                try:
                    pep_added_map[frow["pdb_id"]] = int(frow["peptide_added"])
                except ValueError:
                    pass
        if not pep_added_map:
            print("  --max-added: peptide_added column missing from filtered_systems.tsv — rerun step 2 first.")
            return
        before = len(rows)
        rows = [r for r in rows if pep_added_map.get(r["pdb_id"], 0) <= max_added]
        print(f"  --max-added {max_added}: {len(rows)} of {before} entries have ≤ {max_added} PDBFixer-added peptide residues")
        if not rows:
            print("  No entries within peptide-added limit.")
            return

    if single_chain:
        before = len(rows)
        # receptor field is "AB:300" — single chain has exactly one letter before ':'
        rows = [r for r in rows if len(r["receptor"].split(":")[0]) == 1]
        print(f"  --single-chain: {len(rows)} of {before} entries have a single receptor chain")
        if not rows:
            print("  No single-chain receptor entries found.")
            return

    if max_receptor is not None:
        before = len(rows)
        rows = [r for r in rows if r["_rec_len"] <= max_receptor]
        print(f"  --max-receptor {max_receptor}: {len(rows)} of {before} entries have receptor ≤ {max_receptor} residues")
        if not rows:
            print("  No entries within receptor size limit.")
            return

    if max_peptide is not None:
        before = len(rows)
        rows = [r for r in rows if r["_pep_len"] <= max_peptide]
        print(f"  --max-peptide {max_peptide}: {len(rows)} of {before} entries have peptide ≤ {max_peptide} residues")
        if not rows:
            print("  No entries within peptide size limit.")
            return

    selected = diverse_sample(rows, n)

    rmsd_cache = None
    if rmsd_threshold > 0:
        selected, n_replaced, min_rmsd, rmsd_cache = _redraw_similar_receptors(
            selected, rows, rmsd_threshold, FIXED_DIR
        )

    selected.sort(key=lambda r: r["_dg"])

    rows.sort(key=lambda r: r["_dg"])
    pool_dg_min  = rows[0]["_dg"]
    pool_dg_max  = rows[-1]["_dg"]
    pool_pep_min = min(r["_pep_len"] for r in rows)
    pool_pep_max = max(r["_pep_len"] for r in rows)

    with open(out, "w", encoding="ascii") as fh:
        fh.write("\t".join(SUMMARY_COLS) + "\n")
        for r in selected:
            fh.write("\t".join(r[c] for c in SUMMARY_COLS) + "\n")

    print(f"\n  Stratified sample: {len(selected)} of {len(rows)} structures  "
          f"(ΔG × peptide-length,  "
          f"ΔG {pool_dg_min:.2f}–{pool_dg_max:.2f} kcal/mol,  "
          f"peptide {pool_pep_min}–{pool_pep_max} aa)")
    print(f"  Written {out}\n")
    print(f"  {'PDB':<8} {'DeltaG_kcal_mol':>15} {'type':<6}  receptor / peptide")
    print("  " + "-" * 62)
    for r in selected:
        print(f"  {r['pdb_id']:<8} {r['dG_kcal_mol']:>15} {r['binding_type']:<6} "
              f" {r['receptor']}  /  {r['peptide']}")

    s_dg  = [r["_dg"]      for r in selected]
    s_pep = [r["_pep_len"] for r in selected]
    s_rec = [r["_rec_len"] for r in selected]
    print()
    print(f"  {'':8}  {'min':>8} {'max':>8} {'avg':>8}")
    print(f"  {'ΔG (kcal/mol)':<20} {min(s_dg):>8.2f} {max(s_dg):>8.2f} {sum(s_dg)/len(s_dg):>8.2f}")
    print(f"  {'peptide length':<20} {min(s_pep):>8} {max(s_pep):>8} {sum(s_pep)/len(s_pep):>8.1f}")
    print(f"  {'receptor length':<20} {min(s_rec):>8} {max(s_rec):>8} {sum(s_rec)/len(s_rec):>8.1f}")
    if rmsd_threshold > 0:
        if min_rmsd is None:
            rmsd_note = "no comparable pairs (lengths too different)"
        else:
            rmsd_note = f"min pairwise receptor RMSD {min_rmsd:.2f} Å"
        repl_note = f"{n_replaced} replacement(s)" if n_replaced else "no replacements needed"
        print(f"  {'receptor RMSD':<20}   {repl_note}, {rmsd_note}")
    print()

    if rmsd_cache is None:
        rmsd_cache = compute_rmsd_matrix(selected, FIXED_DIR)
    print_rmsd_cluster_report(rmsd_cache, selected)


# ---------------------------------------------------------------------------
# --expand N
# ---------------------------------------------------------------------------

def run_expand(n, single_chain=False, rmsd_threshold=3.0, max_receptor=None,
               max_peptide=None, max_added=None, output_path=None):
    out = Path(output_path) if output_path else SAMPLE_TSV
    if not out.exists():
        print(f"  {out} not found — run --sample N first.")
        return
    if not SUMMARY_TSV.exists():
        print(f"  {SUMMARY_TSV} not found — run step 3 first.")
        return

    existing = []
    for row in load_tsv(out):
        try:
            row["_dg"]      = float(row.get("dG_kcal_mol", "NA"))
            row["_pep_len"] = int(row.get("peptide",  ":0").split(":")[1])
            row["_rec_len"] = int(row.get("receptor", ":0").split(":")[1])
        except (ValueError, IndexError):
            continue
        existing.append(row)

    if not existing:
        print(f"  {out} is empty or unreadable.")
        return

    existing_ids = {r["pdb_id"] for r in existing}
    print(f"  Existing sample: {len(existing)} entries")

    pool = []
    for row in load_tsv(SUMMARY_TSV):
        if row["pdb_id"] in existing_ids:
            continue
        try:
            row["_dg"]      = float(row.get("dG_kcal_mol", "NA"))
            row["_pep_len"] = int(row.get("peptide",  ":0").split(":")[1])
            row["_rec_len"] = int(row.get("receptor", ":0").split(":")[1])
        except (ValueError, IndexError):
            continue
        pool.append(row)

    if max_added is not None:
        filtered_tsv = Path("filtered_systems.tsv")
        if not filtered_tsv.exists():
            print(f"  --max-added requires filtered_systems.tsv — run step 2 first.")
            return
        pep_added_map = {}
        for frow in load_tsv(filtered_tsv):
            if frow.get("accepted") == "Y" and frow.get("peptide_added", ""):
                try:
                    pep_added_map[frow["pdb_id"]] = int(frow["peptide_added"])
                except ValueError:
                    pass
        if not pep_added_map:
            print("  --max-added: peptide_added column missing from filtered_systems.tsv — rerun step 2 first.")
            return
        before = len(pool)
        pool = [r for r in pool if pep_added_map.get(r["pdb_id"], 0) <= max_added]
        print(f"  --max-added {max_added}: {len(pool)} of {before} entries have ≤ {max_added} PDBFixer-added peptide residues")

    if single_chain:
        before = len(pool)
        pool = [r for r in pool if len(r["receptor"].split(":")[0]) == 1]
        print(f"  --single-chain: {len(pool)} of {before} remaining entries have a single receptor chain")

    if max_receptor is not None:
        before = len(pool)
        pool = [r for r in pool if r["_rec_len"] <= max_receptor]
        print(f"  --max-receptor {max_receptor}: {len(pool)} of {before} entries within receptor size limit")

    if max_peptide is not None:
        before = len(pool)
        pool = [r for r in pool if r["_pep_len"] <= max_peptide]
        print(f"  --max-peptide {max_peptide}: {len(pool)} of {before} entries within peptide size limit")

    if not pool:
        print("  No additional entries available.")
        return

    new_entries = diverse_sample(pool, n)

    rmsd_cache = None
    if rmsd_threshold > 0:
        remaining = [r for r in pool if r["pdb_id"] not in {e["pdb_id"] for e in new_entries}]
        new_entries, n_replaced, min_rmsd, rmsd_cache = _redraw_similar_receptors(
            new_entries, remaining, rmsd_threshold, FIXED_DIR, anchors=existing
        )

    combined = existing + new_entries
    combined.sort(key=lambda r: r["_dg"])

    with open(out, "w", encoding="ascii") as fh:
        fh.write("\t".join(SUMMARY_COLS) + "\n")
        for r in combined:
            fh.write("\t".join(r[c] for c in SUMMARY_COLS) + "\n")

    print(f"\n  Added {len(new_entries)} entries  →  sample now {len(combined)} total")
    print(f"  Written {out}\n")
    print(f"  {'PDB':<8} {'DeltaG_kcal_mol':>15} {'type':<6}  receptor / peptide")
    print("  " + "-" * 62)
    for r in new_entries:
        print(f"  {r['pdb_id']:<8} {r['dG_kcal_mol']:>15} {r['binding_type']:<6} "
              f" {r['receptor']}  /  {r['peptide']}")

    s_dg  = [r["_dg"]      for r in new_entries]
    s_pep = [r["_pep_len"] for r in new_entries]
    s_rec = [r["_rec_len"] for r in new_entries]
    print()
    print(f"  {'':8}  {'min':>8} {'max':>8} {'avg':>8}  (new entries only)")
    print(f"  {'ΔG (kcal/mol)':<20} {min(s_dg):>8.2f} {max(s_dg):>8.2f} {sum(s_dg)/len(s_dg):>8.2f}")
    print(f"  {'peptide length':<20} {min(s_pep):>8} {max(s_pep):>8} {sum(s_pep)/len(s_pep):>8.1f}")
    print(f"  {'receptor length':<20} {min(s_rec):>8} {max(s_rec):>8} {sum(s_rec)/len(s_rec):>8.1f}")
    if rmsd_threshold > 0:
        if min_rmsd is None:
            rmsd_note = "no comparable pairs (lengths too different)"
        else:
            rmsd_note = f"min pairwise receptor RMSD {min_rmsd:.2f} Å (vs all {len(combined)} entries)"
        repl_note = f"{n_replaced} replacement(s)" if n_replaced else "no replacements needed"
        print(f"  {'receptor RMSD':<20}   {repl_note}, {rmsd_note}")
    print()

    if rmsd_cache is None:
        rmsd_cache = compute_rmsd_matrix(new_entries, FIXED_DIR)
    print_rmsd_cluster_report(rmsd_cache, new_entries)


# ---------------------------------------------------------------------------
# PyMOL viewer (shared by default mode and --view-sample)
# ---------------------------------------------------------------------------

def run_viewer(cif_files, label):
    import pymol
    from pymol import cmd
    launched_from_python = not __import__("sys").argv[0].endswith("pymol")
    if launched_from_python:
        pymol.finish_launching(["pymol", "-q"])
    if not cif_files:
        print(f"  No CIF files to display.")
        return

    idx = [0]

    def load_current():
        idx[0] = min(idx[0], len(cif_files) - 1)
        path = cif_files[idx[0]]
        pdb_id = path.stem
        cmd.delete("all")
        cmd.load(str(path), pdb_id)
        cmd.util.cbc()
        cmd.zoom()
        print(f"[{idx[0] + 1}/{len(cif_files)}]  {pdb_id}")

    def next_entry():
        idx[0] += 1
        load_current()

    def prev_entry():
        if idx[0] > 0:
            idx[0] -= 1
        load_current()

    cmd.set_key("right", next_entry)
    cmd.set_key("left", prev_entry)
    print(f"  {label}")
    print("  Left/right arrow keys to navigate.")
    load_current()
    if launched_from_python:
        import threading
        threading.Timer(0.5, lambda: cmd.viewport(1600, 1000)).start()
    else:
        cmd.viewport(1600, 1000)


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

_output_path = None
if "--output" in sys.argv:
    try:
        _output_path = sys.argv[sys.argv.index("--output") + 1]
    except IndexError:
        pass
_active_tsv = Path(_output_path) if _output_path else SAMPLE_TSV

if "-h" in sys.argv or "--help" in sys.argv:
    print(__doc__)

elif "--clean" in sys.argv:
    if _active_tsv.exists():
        _active_tsv.unlink()
        print(f"Removed {_active_tsv}")

elif "--expand" in sys.argv:
    idx_flag = sys.argv.index("--expand")
    try:
        _n = int(sys.argv[idx_flag + 1])
    except (IndexError, ValueError):
        _n = 10
    _rmsd = 0.0
    if "--receptor-rmsd" in sys.argv:
        try:
            _rmsd = float(sys.argv[sys.argv.index("--receptor-rmsd") + 1])
        except (IndexError, ValueError):
            pass
    _max_rec = None
    if "--max-receptor" in sys.argv:
        try:
            _max_rec = int(sys.argv[sys.argv.index("--max-receptor") + 1])
        except (IndexError, ValueError):
            pass
    _max_pep = None
    if "--max-peptide" in sys.argv:
        try:
            _max_pep = int(sys.argv[sys.argv.index("--max-peptide") + 1])
        except (IndexError, ValueError):
            pass
    _max_added = None
    if "--max-added" in sys.argv:
        try:
            _max_added = int(sys.argv[sys.argv.index("--max-added") + 1])
        except (IndexError, ValueError):
            pass
    run_expand(_n, single_chain="--single-chain" in sys.argv, rmsd_threshold=_rmsd,
               max_receptor=_max_rec, max_peptide=_max_pep, max_added=_max_added,
               output_path=_output_path)

elif "--sample" in sys.argv:
    idx_flag = sys.argv.index("--sample")
    try:
        _n = int(sys.argv[idx_flag + 1])
    except (IndexError, ValueError):
        _n = 20
    _rmsd = 0.0
    if "--receptor-rmsd" in sys.argv:
        try:
            _rmsd = float(sys.argv[sys.argv.index("--receptor-rmsd") + 1])
        except (IndexError, ValueError):
            pass
    _max_rec = None
    if "--max-receptor" in sys.argv:
        try:
            _max_rec = int(sys.argv[sys.argv.index("--max-receptor") + 1])
        except (IndexError, ValueError):
            pass
    _max_pep = None
    if "--max-peptide" in sys.argv:
        try:
            _max_pep = int(sys.argv[sys.argv.index("--max-peptide") + 1])
        except (IndexError, ValueError):
            pass
    _max_added = None
    if "--max-added" in sys.argv:
        try:
            _max_added = int(sys.argv[sys.argv.index("--max-added") + 1])
        except (IndexError, ValueError):
            pass
    run_sample(_n, single_chain="--single-chain" in sys.argv, rmsd_threshold=_rmsd,
               max_receptor=_max_rec, max_peptide=_max_pep, max_added=_max_added,
               overwrite="--overwrite" in sys.argv, output_path=_output_path)

elif "--view-sample" in sys.argv:
    if not _active_tsv.exists():
        print(f"  {_active_tsv} not found — run --sample N first.")
    else:
        rows = load_tsv(_active_tsv)
        for r in rows:
            try:
                r["_dg"] = float(r.get("dG_kcal_mol", "0"))
            except ValueError:
                r["_dg"] = 0.0
        rows.sort(key=lambda r: r["_dg"])
        cif_files = [FIXED_DIR / f"{r['pdb_id']}.cif" for r in rows
                     if (FIXED_DIR / f"{r['pdb_id']}.cif").exists()]
        run_viewer(cif_files, f"Viewing sample ({len(cif_files)} structures from {_active_tsv})")

else:
    import pymol
    from pymol import cmd
    launched_from_python = not sys.argv[0].endswith("pymol")
    if launched_from_python:
        pymol.finish_launching(["pymol", "-q"])

    def get_all_cifs():
        return sorted(FIXED_DIR.glob("*.cif"))

    idx = [0]

    def load_current():
        cif_files = get_all_cifs()
        if not cif_files:
            print(f"  No CIF files found in {FIXED_DIR}/ yet.")
            return
        idx[0] = min(idx[0], len(cif_files) - 1)
        path = cif_files[idx[0]]
        pdb_id = path.stem
        cmd.delete("all")
        cmd.load(str(path), pdb_id)
        cmd.util.cbc()
        cmd.zoom()
        print(f"[{idx[0] + 1}/{len(cif_files)}]  {pdb_id}")

    def next_entry():
        idx[0] += 1
        load_current()

    def prev_entry():
        if idx[0] > 0:
            idx[0] -= 1
        load_current()

    cmd.set_key("right", next_entry)
    cmd.set_key("left", prev_entry)
    print(f"  Browsing {FIXED_DIR}/ (rescanned on each keypress).")
    print("  Left/right arrow keys to navigate.")
    load_current()
    if launched_from_python:
        import threading
        threading.Timer(0.5, lambda: cmd.viewport(1600, 1000)).start()
    else:
        cmd.viewport(1600, 1000)
