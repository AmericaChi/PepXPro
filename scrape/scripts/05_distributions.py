"""
Step 5: Generate box plot distributions for all numeric columns in summary.tsv,
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

from modules.utils import load_tsv, setup_logging

OUTPUTS_DIR = _root / "outputs"
LOGS_DIR    = _root / "logs"
SUMMARY_TSV = OUTPUTS_DIR / "summary.tsv"
PLOTS_DIR   = OUTPUTS_DIR / "plots"
LOG_FILE    = LOGS_DIR / "05_distributions.log"

COL_LABELS = {
    "dG_kcal_mol":     "ΔG (kcal/mol)",
    "resolution":      "Resolution (Å)",
    "year":            "Deposition year",
    "_pep_len":        "Peptide length",
    "pep_sasa_burial": "Buried SASA",
    "n_contacts":      "Interface contacts",
    "n_hbonds":        "H-bond proxies",
}

YLABEL_OVERRIDES = {
    "pep_sasa_burial": r"$f_\mathrm{buried}$",
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
    import seaborn as sns
    nonzero = [(l, c, col) for l, c, col in zip(labels, counts, colors) if c > 0]
    if not nonzero:
        ax.set_visible(False)
        return
    lbls, cnts, cols = zip(*nonzero)
    bars = ax.bar(range(len(lbls)), cnts, color=cols, width=0.65,
                  edgecolor="white", linewidth=0.6)
    ax.set_xticks(range(len(lbls)))
    ax.set_xticklabels(lbls, rotation=40, ha="right", fontsize=6.5)
    ax.set_title(title, fontsize=8, fontweight="bold", pad=4)
    ax.set_ylabel("count", fontsize=7)
    ax.tick_params(axis="y", labelsize=7)
    total = sum(cnts)
    for bar, cnt in zip(bars, cnts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + total * 0.015,
                str(cnt), ha="center", va="bottom", fontsize=6)
    sns.despine(ax=ax)


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
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print a brief start/done line to stdout.")
    parser.add_argument("--log", metavar="FILE", default=str(LOG_FILE),
                        help=f"Log file path (default: {LOG_FILE.name}).")
    args = parser.parse_args()

    summary_tsv = Path(args.summary_tsv) if args.summary_tsv else SUMMARY_TSV
    if args.output_dir:
        plots_dir = Path(args.output_dir)
    else:
        stem = summary_tsv.stem  # e.g. "summary", "summary_merged", "summary_PPIKB"
        dir_name = "plots" + stem[len("summary"):]  # -> "plots", "plots_merged", "plots_PPIKB"
        plots_dir = summary_tsv.parent / dir_name

    out_path   = plots_dir / "summary_distributions.png"
    stats_path = plots_dir / "summary_statistics.tsv"

    if args.clean:
        for p in (out_path, stats_path):
            if p.exists():
                p.unlink()
                print(f"Removed {p}")
        return

    if not summary_tsv.exists():
        sys.stderr.write(f"Error: {summary_tsv} not found.\n")
        sys.exit(1)

    log_file = Path(args.log)
    _logger, stdout = setup_logging(log_file, args.verbose)
    if args.verbose:
        stdout.write("[scrape/05-distributions] starting ...\n")
        stdout.flush()

    import matplotlib.pyplot as plt
    import seaborn as sns

    import pandas as pd

    sns.set_theme(style="ticks", context="paper", font_scale=0.9,
                  rc={"axes.spines.top": False, "axes.spines.right": False})

    rows = load_tsv(summary_tsv)
    print(f"  Loaded {len(rows)} entries from {summary_tsv}")

    for r in rows:
        val = r.get("peptide", "")
        if ":" in val:
            try:
                r["_pep_len"] = str(int(val.split(":")[1]))
            except (ValueError, IndexError):
                pass

    numeric = parse_numeric_columns(rows)
    if not numeric:
        print("  No numeric columns found.")
        sys.exit(1)

    # Separate pct columns from scalar columns
    PCT_PEP  = ["peptide_pct_helix",   "peptide_pct_beta",   "peptide_pct_coil"]
    PCT_REC  = ["receptor_pct_helix",  "receptor_pct_beta",  "receptor_pct_coil"]
    PCT_ALL  = set(PCT_PEP + PCT_REC)
    scalar_cols = [c for c in numeric if c not in PCT_ALL]
    print(f"  Numeric columns: {', '.join(scalar_cols)}")

    for r in rows:
        r["_peptide_ss_class"] = classify_peptide_ss(
            r.get("peptide_pct_helix", ""), r.get("peptide_pct_beta", ""),
            r.get("peptide_pct_coil", ""))

    ss_counts   = count_categories(rows, "_peptide_ss_class", SS_ORDER)
    func_counts = count_categories(rows, "receptor_func_class", FUNC_ORDER)

    N_COLS  = 4
    n_rows  = 3   # 12 panels fixed → 3 rows × 4 cols
    palette = sns.color_palette("muted")
    box_color = palette[0]

    fig, axes = plt.subplots(n_rows, N_COLS,
                             figsize=(N_COLS * 2.6, n_rows * 2.8),
                             squeeze=False)

    def _boxstrip(ax, data, label, integer_y=False, ylabel=None):
        """Narrow boxplot + jittered strip for a single numeric series."""
        sns.boxplot(y=data, ax=ax, color=box_color, linewidth=0.8, width=0.18,
                    fliersize=0)
        sns.stripplot(y=data, ax=ax, color=box_color, size=1.4, alpha=0.25,
                      jitter=True)
        ax.set_title(label, fontsize=8, fontweight="bold", pad=3)
        ax.set_ylabel(ylabel if ylabel is not None else label, fontsize=7)
        ax.set_xlabel("")
        ax.tick_params(labelsize=7)
        ax.set_xticks([])
        if integer_y:
            import matplotlib.ticker as ticker
            ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        sns.despine(ax=ax)

    def _grouped_num(ax, col_map, title, ylabel="residues"):
        """Grouped boxplot + strip for two numeric columns on the same axis."""
        frames = []
        for lbl, col in col_map.items():
            vals = numeric.get(col, [])
            frames.append(pd.DataFrame({"val": vals, "type": lbl}))
        df = pd.concat(frames, ignore_index=True)
        sns.boxplot(data=df, x="type", y="val", ax=ax,
                    palette="muted", linewidth=0.8, width=0.45, fliersize=0)
        sns.stripplot(data=df, x="type", y="val", ax=ax,
                      palette="muted", size=1.2, alpha=0.2, jitter=True)
        ax.set_title(title, fontsize=8, fontweight="bold", pad=3)
        ax.set_ylabel(ylabel, fontsize=7)
        ax.set_xlabel("")
        ax.tick_params(labelsize=7)
        import matplotlib.ticker as ticker
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        sns.despine(ax=ax)

    def _grouped_pct(ax, cols, title):
        """Grouped boxplot + strip for a set of pct columns (converted to %)."""
        ss_labels = [c.split("_pct_")[1] for c in cols]
        frames = []
        for col, lbl in zip(cols, ss_labels):
            vals = [v * 100 for v in numeric.get(col, [])]
            frames.append(pd.DataFrame({"pct": vals, "type": lbl}))
        df = pd.concat(frames, ignore_index=True)
        sns.boxplot(data=df, x="type", y="pct", ax=ax,
                    palette="muted", linewidth=0.8, width=0.45, fliersize=0)
        sns.stripplot(data=df, x="type", y="pct", ax=ax,
                      palette="muted", size=1.2, alpha=0.2, jitter=True)
        ax.set_title(title, fontsize=8, fontweight="bold", pad=3)
        ax.set_ylabel("%", fontsize=7)
        ax.set_xlabel("")
        ax.tick_params(labelsize=7)
        ax.set_ylim(-2, 102)
        sns.despine(ax=ax)

    INTEGER_COLS = {"_pep_len"}

    # Explicit panel order
    PANELS = [
        ("scalar",      "dG_kcal_mol"),
        ("scalar",      "resolution"),
        ("scalar",      "year"),
        ("scalar",      "_pep_len"),
        ("pct",         PCT_PEP,   "Peptide SS (%)"),
        ("bar",         "_peptide_ss_class", SS_ORDER, SS_COLORS, "Peptide SS class"),
        ("grouped_num", {"Peptide": "peptide_added", "Receptor": "receptor_added"},
                        "Added residues"),
        ("pct",         PCT_REC,   "Receptor SS (%)"),
        ("bar",         "receptor_func_class", FUNC_ORDER, FUNC_COLORS, "Receptor function"),
        ("scalar",      "pep_sasa_burial"),
        ("scalar",      "n_contacts"),
        ("scalar",      "n_hbonds"),
    ]

    for i, spec in enumerate(PANELS):
        ax = axes[i // N_COLS, i % N_COLS]
        kind = spec[0]
        if kind == "scalar":
            col = spec[1]
            label = COL_LABELS.get(col, col.replace("_", " "))
            if numeric.get(col):
                _boxstrip(ax, numeric[col], label,
                          integer_y=(col in INTEGER_COLS),
                          ylabel=YLABEL_OVERRIDES.get(col))
            else:
                ax.set_visible(False)
        elif kind == "pct":
            _grouped_pct(ax, spec[1], spec[2])
        elif kind == "grouped_num":
            _grouped_num(ax, spec[1], spec[2])
        elif kind == "bar":
            _, col, order, colors, title = spec
            plot_bar(ax, order, count_categories(rows, col, order), colors, title)

    for j in range(len(PANELS), n_rows * N_COLS):
        axes[j // N_COLS, j % N_COLS].set_visible(False)

    fig.suptitle(f"PepXPro dataset  ·  n = {len(rows)}",
                 fontsize=9, fontweight="bold", y=1.01)
    fig.tight_layout(w_pad=1.5, h_pad=2.0)

    plots_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"  Saved {out_path}")
    plt.close(fig)

    write_stats_tsv(numeric, stats_path)

    if args.verbose:
        stdout.write(f"[scrape/05-distributions] done  —  {out_path.name}\n")
        stdout.flush()


if __name__ == "__main__":
    main()
