"""
03_compare_runs.py — Plot filtering and sampling distribution for a single run.

Reads pool_quality.tsv, pool_redundancy.tsv, and sample.tsv from a run directory
and produces a 1x3 panel plot:
  Panel 1: ΔG histogram (quality → redundancy → sample)
  Panel 2: peptide length histogram
  Panel 3: scatter ΔG × peptide length

Output: run_XX/plots/run_comparison.png

Usage:
  python 03_compare_runs.py --run run_01
"""

import argparse
import json
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, enrich_rows, SS_CLASS_ORDER, FUNC_CLASS_ORDER

OUTPUTS_DIR = _root / "outputs"

C_QUALITY  = "#AAAAAA"   # gray
C_REDUN    = "#5C85D6"   # blue
C_SAMPLE   = "#E8602C"   # orange

DIM_KEY = {
    "dg":           "_dg",
    "pep_len":      "_pep_len",
    "rec_len":      "_rec_len",
    "n_rec_chains": "_n_rec_chains",
    "ss_class":     "_ss_class",
}

DIM_LABEL = {
    "dg":           "ΔG (kcal/mol)",
    "pep_len":      "Peptide length (residues)",
    "rec_len":      "Receptor length (residues)",
    "n_rec_chains": "N receptor chains",
    "ss_class":     "SS class",
}


def _dim_vals(rows, dim):
    key = DIM_KEY.get(dim, f"_{dim}")
    return [r.get(key) for r in rows]


def _apply_cat_ticks(ax, order, axis="y"):
    """Set integer ticks with categorical labels on an axis."""
    ticks = range(len(order))
    if axis == "y":
        ax.set_yticks(ticks)
        ax.set_yticklabels(order, fontsize=7)
    else:
        ax.set_xticks(ticks)
        ax.set_xticklabels(order, rotation=30, ha="right", fontsize=7)


CAT_DIMS = {
    "ss_class":   SS_CLASS_ORDER,
    "func_class": FUNC_CLASS_ORDER,
}


def _load(path, label):
    if not path.exists():
        print(f"  Error: {path} not found.")
        print(f"  Re-run 01_gen_sample.py without --no-save-pools to generate {label}.")
        sys.exit(1)
    rows = enrich_rows(load_tsv(path))
    print(f"  Loaded {len(rows):>4} entries from {path.name}")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Plot filtering and sampling distribution for a single run.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--run", metavar="RUN_DIR",
                        help="Run directory name inside gen_sample/outputs/, e.g. run_01")
    parser.add_argument("--run-dir", metavar="PATH",
                        help="Full path to any directory containing pool_quality.tsv, "
                             "pool_redundancy.tsv, and sample.tsv.")
    parser.add_argument("--dims", nargs=2, metavar=("DIM_X", "DIM_Y"),
                        help="Override sampling dims from run_config (e.g. --dims dg ss_class).")
    parser.add_argument("--output", metavar="FILE",
                        help="Override output path (default: run_XX/plots/run_comparison[_dimx_dimy].png).")
    args = parser.parse_args()

    if args.run_dir:
        run_dir = Path(args.run_dir)
    elif args.run:
        run_dir = OUTPUTS_DIR / args.run
    else:
        parser.error("one of --run or --run-dir is required")

    if not run_dir.exists():
        print(f"  Error: {run_dir} not found.")
        sys.exit(1)

    quality_rows    = _load(run_dir / "pool_quality.tsv",    "pool_quality.tsv")
    redundancy_rows = _load(run_dir / "pool_redundancy.tsv", "pool_redundancy.tsv")
    sample_rows     = _load(run_dir / "sample.tsv",          "sample.tsv")

    # Read sampling method and dims from run_config.json (--dims overrides)
    config_path = run_dir / "run_config.json"
    sampling_method = "unknown"
    dims_config = ["dg", "pep_len"]
    if config_path.exists():
        with open(config_path, encoding="ascii") as fh:
            cfg = json.load(fh)
        sampling_method = cfg.get("sampling_method", "unknown")
        dims_config = cfg.get("dims") or ["dg", "pep_len"]

    dims  = list(args.dims) if args.dims else dims_config
    dim_x = dims[0] if len(dims) >= 1 else "dg"
    dim_y = dims[1] if len(dims) >= 2 else "pep_len"

    def _extract(rows, dx, dy):
        """Return (dg, x, y) lists with paired filtering — only rows where both x and y are non-None."""
        dg, x, y = [], [], []
        for r in rows:
            xv = _dim_vals([r], dx)[0]
            yv = _dim_vals([r], dy)[0]
            dg.append(r["_dg"])
            if xv is not None and yv is not None:
                x.append(xv)
                y.append(yv)
        return dg, x, y

    dg_q, x_q, y_q = _extract(quality_rows,    dim_x, dim_y)
    dg_r, x_r, y_r = _extract(redundancy_rows, dim_x, dim_y)
    dg_s, x_s, y_s = _extract(sample_rows,     dim_x, dim_y)

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    BINS = 20
    label_x = DIM_LABEL.get(dim_x, dim_x)
    label_y = DIM_LABEL.get(dim_y, dim_y)

    # Panel 1: dG histogram (always)
    ax = axes[0]
    ax.hist(dg_q, bins=BINS, color=C_QUALITY, alpha=0.6, label=f"quality ({len(dg_q)})")
    ax.hist(dg_r, bins=BINS, color=C_REDUN,   alpha=0.6, label=f"redundancy ({len(dg_r)})")
    ax.hist(dg_s, bins=BINS, color=C_SAMPLE,  alpha=0.9, label=f"sample ({len(dg_s)})")
    ax.set_xlabel("ΔG (kcal/mol)")
    ax.set_ylabel("Count")
    ax.set_title("Binding affinity distribution")
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.legend(fontsize=8)

    # Panel 2: dim_y histogram
    ax = axes[1]
    if dim_y in CAT_DIMS and y_q:
        cat_order = CAT_DIMS[dim_y]
        bins_y = range(len(cat_order) + 1)
        ax.hist(y_q, bins=bins_y, color=C_QUALITY, alpha=0.6, label=f"quality ({len(y_q)})")
        ax.hist(y_r, bins=bins_y, color=C_REDUN,   alpha=0.6, label=f"redundancy ({len(y_r)})")
        ax.hist(y_s, bins=bins_y, color=C_SAMPLE,  alpha=0.9, label=f"sample ({len(y_s)})")
        _apply_cat_ticks(ax, cat_order, axis="x")
    elif y_q:
        bins_y = range(int(min(y_q)), int(max(y_q)) + 2)
        ax.hist(y_q, bins=bins_y, color=C_QUALITY, alpha=0.6, label=f"quality ({len(y_q)})")
        ax.hist(y_r, bins=bins_y, color=C_REDUN,   alpha=0.6, label=f"redundancy ({len(y_r)})")
        ax.hist(y_s, bins=bins_y, color=C_SAMPLE,  alpha=0.9, label=f"sample ({len(y_s)})")
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.set_xlabel(label_y)
    ax.set_ylabel("Count")
    ax.set_title(f"{label_y} distribution")
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.legend(fontsize=8)

    # Panel 3: scatter dim_x × dim_y
    ax = axes[2]
    ax.scatter(x_q, y_q, s=10, color=C_QUALITY, alpha=0.3, label=f"quality ({len(x_q)})")
    ax.scatter(x_r, y_r, s=10, color=C_REDUN,   alpha=0.5, label=f"redundancy ({len(x_r)})")
    ax.scatter(x_s, y_s, s=30, color=C_SAMPLE,  alpha=1.0, label=f"sample ({len(x_s)})",
               zorder=5)
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    ax.set_title(f"Coverage: {label_x} × {label_y}  [{sampling_method}]")
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    if dim_y in CAT_DIMS:
        _apply_cat_ticks(ax, CAT_DIMS[dim_y], axis="y")
    else:
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.legend(fontsize=8)

    fig.suptitle(f"{args.run} — pool filtering and sampling", fontsize=12)
    fig.tight_layout()

    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
    else:
        plots_dir = run_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        suffix = f"_{dim_x}_{dim_y}" if args.dims else ""
        out = plots_dir / f"run_comparison{suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
