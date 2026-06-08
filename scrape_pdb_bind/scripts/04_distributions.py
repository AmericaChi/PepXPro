"""
Step 4: Generate box plot distributions for all numeric columns in summary.tsv,
plus bar charts for peptide secondary structure class and receptor functional class.

Outputs: outputs/plots/summary_distributions.png
         outputs/plots/summary_statistics.tsv
"""

import argparse
import math
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv

OUTPUTS_DIR = _root / "outputs"
SUMMARY_TSV = OUTPUTS_DIR / "summary.tsv"
PLOTS_DIR   = OUTPUTS_DIR / "plots"

COL_LABELS = {
    "dG_kcal_mol": "ΔG (kcal/mol)",
    "resolution":  "Resolution (Å)",
    "year":        "Deposition year",
}

SS_ORDER   = ["helix", "partial-helix", "beta", "mixed", "coil", "other", "NA"]
FUNC_ORDER = ["enzyme", "chaperone", "structural", "signal",
              "regulatory", "peptide-binding", "other", "NA"]

SS_COLORS   = ["#EF5350", "#FF8A65", "#42A5F5", "#AB47BC", "#26A69A", "#78909C", "#BDBDBD"]
FUNC_COLORS = ["#EF5350", "#FF7043", "#FFA726", "#FFCA28",
               "#66BB6A", "#26C6DA", "#78909C", "#BDBDBD"]


def parse_numeric_columns(rows):
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
        if len(vals) >= len(rows) * 0.5:
            numeric[col] = vals
    return numeric


def classify_peptide_ss(pct_helix_str, pct_beta_str, pct_coil_str=""):
    try:
        h = float(pct_helix_str)
        e = float(pct_beta_str)
    except (ValueError, TypeError):
        return "NA"
    if h >= 0.70:               return "helix"
    if e >= 0.50:               return "beta"
    if h >= 0.20 and e >= 0.20: return "mixed"
    if h >= 0.30:               return "partial-helix"
    try:
        c = float(pct_coil_str)
    except (ValueError, TypeError):
        c = 1.0 - h - e
    if c > 0.50:                return "coil"
    return "other"


def count_categories(rows, col, order):
    counts = {k: 0 for k in order}
    for r in rows:
        val = r.get(col, "") or "NA"
        if val not in counts:
            val = "NA"
        counts[val] += 1
    return [counts[k] for k in order]


def plot_bar(ax, labels, counts, colors, title):
    bars = ax.bar(range(len(labels)), counts, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
    ax.set_ylabel("count", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    total = sum(counts)
    for bar, cnt in zip(bars, counts):
        if cnt > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + total * 0.01,
                    str(cnt), ha="center", va="bottom", fontsize=7)


def compute_stats(vals):
    n = len(vals)
    if n == 0:
        return {k: "NA" for k in ("N", "mean", "std", "min", "Q1", "median", "Q3", "max")}
    s = sorted(vals)
    mean = sum(s) / n
    std  = (sum((v - mean) ** 2 for v in s) / n) ** 0.5

    def pct(p):
        idx = (n - 1) * p / 100
        lo, hi = int(idx), min(int(idx) + 1, n - 1)
        return s[lo] + (s[hi] - s[lo]) * (idx - lo)

    return {"N": str(n), "mean": f"{mean:.4f}", "std": f"{std:.4f}",
            "min": f"{s[0]:.4f}", "Q1": f"{pct(25):.4f}",
            "median": f"{pct(50):.4f}", "Q3": f"{pct(75):.4f}", "max": f"{s[-1]:.4f}"}


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
        description="Step 4: Distributions and statistics for summary.tsv.")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--summary-tsv", metavar="FILE", default=None,
                        help="Path to summary TSV (default: outputs/summary.tsv).")
    parser.add_argument("--output-dir", metavar="DIR", default=None,
                        help="Directory for output plots/stats (default: inferred from summary-tsv).")
    args = parser.parse_args()

    summary_tsv = Path(args.summary_tsv) if args.summary_tsv else SUMMARY_TSV
    if args.output_dir:
        plots_dir = Path(args.output_dir)
    else:
        plots_dir = summary_tsv.parent / "plots" if args.summary_tsv else PLOTS_DIR

    out_path   = plots_dir / "summary_distributions.png"
    stats_path = plots_dir / "summary_statistics.tsv"

    if args.clean:
        for p in (out_path, stats_path):
            if p.exists():
                p.unlink()
                print(f"Removed {p}")
        return

    if not summary_tsv.exists():
        print(f"Error: {summary_tsv} not found.")
        sys.exit(1)

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    rows = load_tsv(summary_tsv)
    print(f"  Loaded {len(rows)} entries from {summary_tsv}")

    # Numeric box plots
    numeric = parse_numeric_columns(rows)
    if not numeric:
        print("  No numeric columns found.")
        sys.exit(1)
    num_cols = list(numeric.keys())
    print(f"  Numeric columns: {', '.join(num_cols)}")

    # Categorical bar charts
    # Derive peptide SS class from pct columns
    for r in rows:
        r["_peptide_ss_class"] = classify_peptide_ss(
            r.get("peptide_pct_helix", ""), r.get("peptide_pct_beta", ""),
            r.get("peptide_pct_coil", ""))

    ss_counts   = count_categories(rows, "_peptide_ss_class", SS_ORDER)
    func_counts = count_categories(rows, "receptor_func_class", FUNC_ORDER)

    # Grid: numeric box plots + 2 bar charts
    n_total = len(num_cols) + 2
    N_COLS  = 6
    n_rows  = math.ceil(n_total / N_COLS)

    fig = plt.figure(figsize=(4 * N_COLS, 5 * n_rows))
    gs  = gridspec.GridSpec(n_rows, N_COLS, figure=fig, wspace=0.5, hspace=0.8)

    # Box plots
    for i, col in enumerate(num_cols):
        ax = fig.add_subplot(gs[i // N_COLS, i % N_COLS])
        bp = ax.boxplot(numeric[col], patch_artist=True, widths=0.5,
                        medianprops=dict(color="black", linewidth=2))
        bp["boxes"][0].set_facecolor("#90CAF9")
        label = COL_LABELS.get(col, col)
        ax.set_title(label, fontsize=10, fontweight="bold", pad=6)
        ax.set_ylabel(label, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.set_xticks([])

    # Bar chart — peptide SS class
    i_ss = len(num_cols)
    ax_ss = fig.add_subplot(gs[i_ss // N_COLS, i_ss % N_COLS])
    plot_bar(ax_ss, SS_ORDER, ss_counts, SS_COLORS, "Peptide SS class")

    # Bar chart — receptor functional class
    i_fn = len(num_cols) + 1
    ax_fn = fig.add_subplot(gs[i_fn // N_COLS, i_fn % N_COLS])
    plot_bar(ax_fn, FUNC_ORDER, func_counts, FUNC_COLORS, "Receptor function")

    # Hide unused grid cells
    for j in range(n_total, n_rows * N_COLS):
        fig.add_subplot(gs[j // N_COLS, j % N_COLS]).set_visible(False)

    fig.suptitle(f"Summary distributions  |  n={len(rows)} entries",
                 fontsize=13, y=1.01)

    plots_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  Saved {out_path}")
    plt.close(fig)

    write_stats_tsv(numeric, stats_path)


if __name__ == "__main__":
    main()
