"""
sequence_redundancy.py — Sequence-based redundancy filtering.

Receptor: CD-HIT clustering on concatenated receptor chain sequences.
Peptide:  Biopython global pairwise alignment (Needleman-Wunsch, BLOSUM62),
          greedy clustering (best resolution kept as representative).

Tools: cd-hit (external binary), biopython (already in env).

Presets:
  strict  → receptor 40%, peptide 70%   (CD-HIT minimum: 40%)
  medium  → receptor 50%, peptide 80%
  lenient → receptor 70%, peptide 90%
  custom  → user-defined thresholds (receptor threshold must be >= 0.40)
"""

import re
import subprocess
import tempfile
import os
from pathlib import Path

PRESETS = {
    "strict":  {"receptor": 0.40, "peptide": 0.70},
    "medium":  {"receptor": 0.50, "peptide": 0.80},
    "lenient": {"receptor": 0.70, "peptide": 0.90},
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _receptor_seq(row):
    """Concatenate all receptor chain sequences from 'A:SEQ|B:SEQ' format."""
    raw = row.get("receptor_seq", "")
    if not raw:
        return ""
    if "|" in raw:
        return "".join(p.split(":", 1)[1] for p in raw.split("|") if ":" in p)
    if ":" in raw:
        return raw.split(":", 1)[1]
    return raw


def _peptide_seq(row):
    """Return the peptide sequence string."""
    return row.get("peptide_seq", "")


def _check_seq_column(rows, col, method_name):
    """Raise if the required sequence column is missing from summary.tsv."""
    if not rows:
        return
    if col not in rows[0] or not rows[0][col]:
        raise ValueError(
            f"{method_name} requires '{col}' column in summary.tsv. "
            "Re-run scrape_pdb_bind steps 2-3 to generate it."
        )


def _res_sort_key(r):
    """Sort key: best resolution first, NMR/NA last."""
    res = r.get("_resolution")
    return (0, res) if res is not None else (1, 999)


# ---------------------------------------------------------------------------
# CD-HIT word size (required by cd-hit based on threshold)
# ---------------------------------------------------------------------------

def _cdhit_word_size(threshold):
    if threshold >= 0.7:
        return 5
    if threshold >= 0.6:
        return 4
    if threshold >= 0.5:
        return 3
    return 2  # 0.4–0.5; CD-HIT minimum supported threshold is 0.4


# ---------------------------------------------------------------------------
# Receptor filtering — CD-HIT
# ---------------------------------------------------------------------------

def filter_receptor_redundancy(rows, threshold):
    """Remove receptor-redundant entries using CD-HIT.

    Keeps one representative per cluster (best resolution first).
    Returns filtered rows and prints stats.
    threshold: identity fraction [0, 1], e.g. 0.50 for 50%.
    """
    _check_seq_column(rows, "receptor_seq", "receptor sequence redundancy")

    before = len(rows)
    rows_sorted = sorted(rows, key=_res_sort_key)

    fasta_lines = []
    for r in rows_sorted:
        seq = _receptor_seq(r)
        if seq:
            fasta_lines.append(f">{r['pdb_id']}\n{seq}")

    if not fasta_lines:
        print(f"  --seq-id-receptor: no receptor sequences found, skipping")
        return rows

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_in  = os.path.join(tmpdir, "receptors.fasta")
        fasta_out = os.path.join(tmpdir, "receptors_nr")

        with open(fasta_in, "w") as fh:
            fh.write("\n".join(fasta_lines))

        n = _cdhit_word_size(threshold)
        cmd = [
            "cd-hit",
            "-i", fasta_in,
            "-o", fasta_out,
            "-c", str(threshold),
            "-n", str(n),
            "-d", "0",        # full sequence name in .clstr
            "-T", "0",        # use all threads
            "-M", "0",        # unlimited memory
        ]
        if threshold < 0.40:
            raise ValueError(
                f"cd-hit does not support thresholds below 0.40 (got {threshold:.2f}). "
                "Use --seq-id-receptor-threshold 0.40 or higher, or use a preset "
                "(strict=0.40, medium=0.50, lenient=0.70)."
            )
        try:
            subprocess.run(cmd, capture_output=True, check=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"cd-hit failed: {exc.stderr.decode()}")
        except FileNotFoundError:
            raise RuntimeError(
                "cd-hit binary not found. Install with: "
                "conda install -c conda-forge cd-hit"
            )

        rep_ids = _parse_cdhit_clstr(fasta_out + ".clstr")

    filtered = [r for r in rows if r["pdb_id"] in rep_ids]
    print(f"  --seq-id-receptor {threshold:.0%}: {len(filtered)} of {before} entries kept "
          f"({before - len(filtered)} removed)")
    return filtered


def _parse_cdhit_clstr(clstr_path):
    """Return set of pdb_ids that are cluster representatives (lines ending with *)."""
    reps = set()
    with open(clstr_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(">Cluster"):
                continue
            if line.rstrip().endswith("*"):
                m = re.search(r">([^.]+)\.\.\.", line)
                if m:
                    reps.add(m.group(1).lower())
    return reps


# ---------------------------------------------------------------------------
# Peptide filtering — Biopython pairwise
# ---------------------------------------------------------------------------

def filter_peptide_redundancy(rows, threshold):
    """Remove peptide-redundant entries using Biopython global pairwise alignment.

    Greedy clustering: rows sorted by best resolution first. For each
    representative, all sequences with identity > threshold are removed.
    threshold: identity fraction [0, 1], e.g. 0.80 for 80%.
    """
    _check_seq_column(rows, "peptide_seq", "peptide sequence redundancy")

    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score    = -10
    aligner.extend_gap_score  = -0.5
    aligner.mode = "global"

    def identity(seq_a, seq_b):
        if not seq_a or not seq_b:
            return 0.0
        try:
            aln = next(iter(aligner.align(seq_a, seq_b)))
            aln_len = aln.shape[1]
            n_match = sum(
                seq_a[i] == seq_b[j]
                for i, j in zip(*aln.indices)
                if i >= 0 and j >= 0
            )
            return n_match / aln_len if aln_len > 0 else 0.0
        except Exception:
            return 0.0

    before = len(rows)
    rows_sorted = sorted(rows, key=_res_sort_key)

    selected = []
    removed  = set()

    for i, r in enumerate(rows_sorted):
        if r["pdb_id"] in removed:
            continue
        selected.append(r)
        seq_i = _peptide_seq(r)
        for j in range(i + 1, len(rows_sorted)):
            rj = rows_sorted[j]
            if rj["pdb_id"] in removed:
                continue
            if identity(seq_i, _peptide_seq(rj)) > threshold:
                removed.add(rj["pdb_id"])

    print(f"  --seq-id-peptide {threshold:.0%}: {len(selected)} of {before} entries kept "
          f"({before - len(selected)} removed)")
    return selected


# ---------------------------------------------------------------------------
# Combined entry point
# ---------------------------------------------------------------------------

def apply_redundancy_filters(rows, args):
    """Apply active sequence redundancy filters in order: receptor then peptide."""
    rec_preset  = getattr(args, "seq_id_receptor", None)
    pep_preset  = getattr(args, "seq_id_peptide",  None)

    if rec_preset:
        if rec_preset == "custom":
            threshold = args.seq_id_receptor_threshold
        else:
            threshold = PRESETS[rec_preset]["receptor"]
        rows = filter_receptor_redundancy(rows, threshold)

    if pep_preset:
        if pep_preset == "custom":
            threshold = args.seq_id_peptide_threshold
        else:
            threshold = PRESETS[pep_preset]["peptide"]
        rows = filter_peptide_redundancy(rows, threshold)

    return rows
