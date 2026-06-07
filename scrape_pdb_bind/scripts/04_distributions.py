"""
Step 4: Generate box plot distributions for all numeric columns in summary.tsv
and write a statistics TSV with min, max, mean, std, quartiles, etc.

Outputs: outputs/plots/summary_distributions.png
         outputs/plots/summary_statistics.tsv
"""

import argparse
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv

OUTPUTS_DIR = _root / "outputs"
SUMMARY_TSV = OUTPUTS_DIR / "summary.tsv"
PLOTS_DIR   = OUTPUTS_DIR / "plots"


def parse_numeric_columns(rows):
    """Auto-detect numeric columns and return {col: [values]}."""
    if not rows:
        return {}
    numeric = {}
    for col in rows[0].keys():
        vals = []
        for r in rows:
            try:
                vals.append(float(r[col]))
            except (ValueError, TypeError):
                pass
        if len(vals) >= len(rows) * 0.5:  # at least 50% parseable
            numeric[col] = vals
    return numeric


# Human-readable axis labels
COL_LABELS = {
    "dG_kcal_mol":      "ΔG (kcal/mol)",
    "resolution":       "Resolution (Å)",
    "year":             "Deposition year",
}


def compute_stats(vals):
    """Compute descriptive statistics for a list of floats."""
    n = len(vals)
    if n == 0:
        return {k: "NA" for k in ("N", "mean", "std", "min", "Q1", "median", "Q3", "max")}
    s = sorted(vals)
    mean = sum(s) / n
    std  = (sum((v - mean) ** 2 for v in s) / n) ** 0.5

    def percentile(p):
        idx = (n - 1) * p / 100
        lo, hi = int(idx), min(int(idx) + 1, n - 1)
        return s[lo] + (s[hi] - s[lo]) * (idx - lo)

    return {
        "N":      str(n),
        "mean":   f"{mean:.4f}",
        "std":    f"{std:.4f}",
        "min":    f"{s[0]:.4f}",
        "Q1":     f"{percentile(25):.4f}",
        "median": f"{percentile(50):.4f}",
        "Q3":     f"{percentile(75):.4f}",
        "max":    f"{s[-1]:.4f}",
    }


def write_stats_tsv(numeric, out_path):
    stat_keys = ["N", "mean", "std", "min", "Q1", "median", "Q3", "max"]
    with open(out_path, "w", encoding="ascii") as fh:
        fh.write("column\t" + "\t".join(stat_keys) + "\n")
        for col, vals in numeric.items():
            stats = compute_stats(vals)
            fh.write(col + "\t" + "\t".join(stats[k] for k in stat_keys) + "\n")
    print(f"  Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Step 4: Box plot distributions and statistics for numeric columns in summary.tsv.")
    parser.add_argument("--clean", action="store_true",
                        help="Remove summary_distributions.png and summary_statistics.tsv.")
    args = parser.parse_args()

    out_path  = PLOTS_DIR / "summary_distributions.png"
    stats_path = PLOTS_DIR / "summary_statistics.tsv"

    if args.clean:
        for p in (out_path, stats_path):
            if p.exists():
                p.unlink()
                print(f"Removed {p}")
        return

    if not SUMMARY_TSV.exists():
        print(f"Error: {SUMMARY_TSV} not found — run steps 1-3 first.")
        sys.exit(1)

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    rows = load_tsv(SUMMARY_TSV)
    print(f"  Loaded {len(rows)} entries from {SUMMARY_TSV}")

    numeric = parse_numeric_columns(rows)
    if not numeric:
        print("  No numeric columns found.")
        sys.exit(1)

    cols = list(numeric.keys())
    n    = len(cols)
    print(f"  Numeric columns: {', '.join(cols)}")

    import math
    n_cols  = 6
    n_rows  = math.ceil(n / n_cols)
    fig = plt.figure(figsize=(4 * n_cols, 5 * n_rows))
    gs  = gridspec.GridSpec(n_rows, n_cols, figure=fig, wspace=0.5, hspace=0.7)

    for i, col in enumerate(cols):
        vals = numeric[col]
        ax   = fig.add_subplot(gs[i // n_cols, i % n_cols])
        bp   = ax.boxplot(vals, patch_artist=True, widths=0.5,
                          medianprops=dict(color="black", linewidth=2))
        bp["boxes"][0].set_facecolor("#90CAF9")

        label = COL_LABELS.get(col, col)
        ax.set_title(label, fontsize=10, fontweight="bold", pad=6)
        ax.set_ylabel(label, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.set_xticks([])


    # Hide unused cells in the last row
    for j in range(n, n_rows * n_cols):
        fig.add_subplot(gs[j // n_cols, j % n_cols]).set_visible(False)

    fig.suptitle(
        f"Summary distributions  |  n={len(rows)} entries",
        fontsize=13, y=1.01,
    )

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  Saved {out_path}")
    plt.close(fig)

    write_stats_tsv(numeric, stats_path)


if __name__ == "__main__":
    main()
