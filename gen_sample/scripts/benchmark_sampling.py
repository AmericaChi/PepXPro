"""
benchmark_sampling.py — Compare available sampling methods visually.

Runs all three sampling methods on the same filtered pool and produces a single
figure showing selected entries in a 2D scatter (dims[0] x dims[1]).
KDE always samples on dims[0] only; the label reflects this (kde-1D).

Output: gen_sample/outputs/plots/sampling_benchmark.png

Usage:
  python benchmark_sampling.py
  python benchmark_sampling.py --n 30
  python benchmark_sampling.py --n 20 --seed 42
  python benchmark_sampling.py --dims rec_len pep_len
"""

import argparse
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, enrich_rows
from modules.filters import filter_single_chain
from modules.sampling import SAMPLING_METHODS, method_label, DIM_KEYS, AVAILABLE_DIMS

_SCRAPE_DIR = _root.parent / "scrape_pdb_bind"
SUMMARY_TSV = _SCRAPE_DIR / "summary.tsv"
PLOTS_DIR   = _root / "outputs" / "plots"

DIM_LABELS = {
    "dg":           "ΔG (kcal/mol)",
    "pep_len":      "Peptide length (aa)",
    "rec_len":      "Receptor length (aa)",
    "n_rec_chains": "N receptor chains",
}

COLORS = {
    "stratified":       "#2196F3",
    "kennard-stone":    "#E91E63",
    "kde":              "#4CAF50",
    "isolation-forest": "#FF9800",
    "cluster":          "#9C27B0",
}


def main():
    parser = argparse.ArgumentParser(description="Benchmark and compare sampling methods.")
    parser.add_argument("--n",    type=int, default=20,              help="Number of entries to sample (default: 20).")
    parser.add_argument("--seed", type=int, default=None,            help="Random seed for stochastic methods.")
    parser.add_argument("--dims", nargs="+", default=["dg", "pep_len"], metavar="DIM",
                        help=f"Two dimensions for the plot axes. Available: {', '.join(AVAILABLE_DIMS)}. "
                             "Default: dg pep_len. KDE always uses only the first dim.")
    args = parser.parse_args()

    for d in args.dims:
        if d not in AVAILABLE_DIMS:
            parser.error(f"Unknown dimension '{d}'. Available: {', '.join(AVAILABLE_DIMS)}")
    if len(args.dims) < 2:
        parser.error("--dims requires at least two dimensions for the 2D benchmark plot.")

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    if not SUMMARY_TSV.exists():
        print(f"  Error: {SUMMARY_TSV} not found.")
        sys.exit(1)

    rows = load_tsv(SUMMARY_TSV)
    rows = enrich_rows(rows)
    pool = filter_single_chain(rows)
    print(f"  Pool: {len(pool)} entries after single-chain filter")

    x_key = DIM_KEYS[args.dims[0]]
    y_key = DIM_KEYS[args.dims[1]]
    pool_x = [r[x_key] for r in pool]
    pool_y = [r[y_key] for r in pool]

    methods = list(SAMPLING_METHODS.keys())
    results = {}
    for method in methods:
        fn   = SAMPLING_METHODS[method]
        dims = [args.dims[0]] if method == "kde" else args.dims
        selected = fn(pool, args.n, dims=dims, seed=args.seed)
        results[method] = (selected, dims)
        print(f"  {method_label(method, dims)}: {len(selected)} selected")

    # --- Plot ---
    n_methods = len(methods)
    fig = plt.figure(figsize=(5 * n_methods, 5))
    gs  = gridspec.GridSpec(1, n_methods, figure=fig, wspace=0.35)

    xlabel = DIM_LABELS.get(args.dims[0], args.dims[0])
    ylabel = DIM_LABELS.get(args.dims[1], args.dims[1])

    for i, method in enumerate(methods):
        selected, dims = results[method]
        label    = method_label(method, dims)
        dims_str = " × ".join(dims)

        sel_x = [r[x_key] for r in selected]
        sel_y = [r[y_key] for r in selected]

        ax = fig.add_subplot(gs[i])
        ax.scatter(pool_x, pool_y, c="#cccccc", s=12, alpha=0.4, linewidths=0,
                   zorder=1, label=f"Pool (n={len(pool)})")
        ax.scatter(sel_x, sel_y, c=COLORS[method], s=45, alpha=0.85,
                   linewidths=0.4, edgecolors="white", zorder=2,
                   label=f"Selected (n={len(selected)})")

        ax.set_title(f"{label}\n({dims_str})", fontsize=11, fontweight="bold", pad=8)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel if i == 0 else "", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7, loc="upper right", framealpha=0.7)
        ax.set_xlim(min(pool_x) - 0.5, max(pool_x) + 0.5)
        ax.set_ylim(min(pool_y) - 0.5, max(pool_y) + 0.5)

    seed_str = f"seed={args.seed}" if args.seed is not None else "no seed"
    fig.suptitle(
        f"Sampling method comparison  |  N={args.n}, {seed_str}, pool={len(pool)}",
        fontsize=12, y=1.02,
    )

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out = PLOTS_DIR / "sampling_benchmark.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
