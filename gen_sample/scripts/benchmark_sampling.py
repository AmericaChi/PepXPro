"""
benchmark_sampling.py — Compare available sampling methods visually.

Runs all three sampling methods on the same filtered pool and produces a single
figure showing selected entries in (dG x peptide_length) space.

Output: gen_sample/outputs/plots/sampling_benchmark.png

Usage:
  python benchmark_sampling.py
  python benchmark_sampling.py --n 30
  python benchmark_sampling.py --n 20 --seed 42
"""

import argparse
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_root = _here.parent
sys.path.insert(0, str(_root))

from modules.utils import load_tsv, enrich_rows
from modules.filters import filter_single_chain
from modules.sampling import SAMPLING_METHODS

_SCRAPE_DIR   = _root.parent / "scrape_pdb_bind"
SUMMARY_TSV   = _SCRAPE_DIR / "summary.tsv"
PLOTS_DIR     = _root / "outputs" / "plots"


def main():
    parser = argparse.ArgumentParser(description="Benchmark and compare sampling methods.")
    parser.add_argument("--n",    type=int, default=20,   help="Number of entries to sample (default: 20).")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for stochastic methods.")
    args = parser.parse_args()

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    if not SUMMARY_TSV.exists():
        print(f"  Error: {SUMMARY_TSV} not found.")
        sys.exit(1)

    rows = load_tsv(SUMMARY_TSV)
    rows = enrich_rows(rows)
    pool = filter_single_chain(rows)
    print(f"  Pool: {len(pool)} entries after single-chain filter")

    pool_dg  = [r["_dg"]      for r in pool]
    pool_pep = [r["_pep_len"] for r in pool]

    methods = list(SAMPLING_METHODS.keys())
    results = {}
    for method in methods:
        fn = SAMPLING_METHODS[method]
        selected = fn(pool, args.n, seed=args.seed)
        results[method] = selected
        print(f"  {method}: {len(selected)} selected")

    # --- Plot ---
    fig = plt.figure(figsize=(15, 5))
    gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35)

    titles = {
        "diverse":       "Diverse\n(quantile grid)",
        "kennard-stone": "Kennard-Stone\n(max-min distance)",
        "kde":           "KDE Inverse Density\n(tail-favoring)",
    }
    colors = {
        "diverse":       "#2196F3",
        "kennard-stone": "#E91E63",
        "kde":           "#4CAF50",
    }

    for i, method in enumerate(methods):
        ax = fig.add_subplot(gs[i])
        selected = results[method]
        sel_dg  = [r["_dg"]      for r in selected]
        sel_pep = [r["_pep_len"] for r in selected]

        # Pool background
        ax.scatter(pool_dg, pool_pep, c="#cccccc", s=12, alpha=0.4, linewidths=0, zorder=1,
                   label=f"Pool (n={len(pool)})")
        # Selected
        ax.scatter(sel_dg, sel_pep, c=colors[method], s=45, alpha=0.85, linewidths=0.4,
                   edgecolors="white", zorder=2, label=f"Selected (n={len(selected)})")

        ax.set_title(titles[method], fontsize=11, fontweight="bold", pad=8)
        ax.set_xlabel("ΔG (kcal/mol)", fontsize=9)
        ax.set_ylabel("Peptide length (aa)" if i == 0 else "", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7, loc="upper right", framealpha=0.7)

        # Axis limits consistent across panels
        ax.set_xlim(min(pool_dg) - 0.5, max(pool_dg) + 0.5)
        ax.set_ylim(min(pool_pep) - 0.5, max(pool_pep) + 0.5)

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
