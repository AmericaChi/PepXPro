"""
02_plot_sampleDist.py — Plot filtering and sampling distribution for a single run.

Reads pool_quality.tsv, pool_redundancy.tsv, and sample.tsv from a run directory
and produces a 2x2 panel plot:
  Panel 1: dG histogram (quality -> redundancy -> sample)
  Panel 2: dim_y histogram
  Panel 3: scatter dim_x x dim_y
  Panel 4: PCA coverage in feature space (diversity score)

Output: run_XX/plots/run_comparison.png

Usage:
  python 02_plot_sampleDist.py --run run_01
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

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
    "dg":           "dG (kcal/mol)",
    "pep_len":      "Peptide length (residues)",
    "rec_len":      "Receptor length (residues)",
    "n_rec_chains": "N receptor chains",
    "ss_class":     "SS class",
}

# All numeric fields available from enrich_rows for PCA
PCA_FIELDS = [
    "_dg", "_pep_len", "_rec_len", "_n_rec_chains",
    "_pep_helix", "_pep_beta", "_pep_coil",
    "_rec_helix", "_rec_beta", "_rec_coil",
    "_n_contacts", "_n_hbonds",
    "_resolution", "_binding_kd",
]


def _dim_vals(rows, dim):
    key = DIM_KEY.get(dim, f"_{dim}")
    return [r.get(key) for r in rows]


def _apply_cat_ticks(ax, order, axis="y"):
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


def _pca_matrix(rows, fields):
    """Build float matrix from rows, keeping only cols with >50% non-None coverage."""
    valid_fields = [f for f in fields
                    if sum(1 for r in rows if r.get(f) is not None) > len(rows) * 0.5]
    mat = []
    for r in rows:
        row_vals = [r.get(f) if r.get(f) is not None else 0.0 for f in valid_fields]
        mat.append(row_vals)
    return np.array(mat, dtype=float), valid_fields


def _internal_diversity(X_scaled):
    """Mean pairwise Euclidean distance in standardized feature space."""
    if len(X_scaled) < 2:
        return 0.0
    diffs = X_scaled[:, None, :] - X_scaled[None, :, :]
    dists = np.sqrt((diffs ** 2).sum(axis=-1))
    n = len(X_scaled)
    return dists.sum() / (n * (n - 1))


def _plot_pca(ax, quality_rows, redundancy_rows, sample_rows):
    """Returns list of PCA feature names used, or [] on failure."""
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        ax.text(0.5, 0.5, "scikit-learn required for PCA panel",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_title("PCA coverage (unavailable)")
        return []

    mat_q, fields = _pca_matrix(quality_rows, PCA_FIELDS)
    if mat_q.shape[0] < 3 or mat_q.shape[1] < 2:
        ax.text(0.5, 0.5, "Not enough data for PCA",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        return []

    sc  = StandardScaler()
    Xq  = sc.fit_transform(mat_q)
    pca = PCA(n_components=2, random_state=42)
    Xq2 = pca.fit_transform(Xq)
    v1, v2 = pca.explained_variance_ratio_[:2] * 100

    def _project(rows):
        # Use the same fields as quality pool — don't re-filter per subset
        mat = np.array([[r.get(f) if r.get(f) is not None else 0.0 for f in fields]
                        for r in rows], dtype=float)
        if len(mat) == 0:
            return np.empty((0, 2))
        return pca.transform(sc.transform(mat))

    Xr2 = _project(redundancy_rows)
    Xs2 = _project(sample_rows)

    mat_s = np.array([[r.get(f) if r.get(f) is not None else 0.0 for f in fields]
                      for r in sample_rows], dtype=float)
    Xs_scaled = sc.transform(mat_s) if len(mat_s) > 0 else np.empty((0, len(fields)))
    diversity  = _internal_diversity(Xs_scaled)

    ax.scatter(Xq2[:, 0], Xq2[:, 1], s=10, color=C_QUALITY, alpha=0.35,
               zorder=1, label=f"quality ({len(Xq2)})")
    if len(Xr2):
        ax.scatter(Xr2[:, 0], Xr2[:, 1], s=12, color=C_REDUN, alpha=0.5,
                   zorder=2, label=f"redundancy ({len(Xr2)})")
    if len(Xs2):
        ax.scatter(Xs2[:, 0], Xs2[:, 1], s=35, color=C_SAMPLE, alpha=0.9,
                   zorder=3, label=f"sample ({len(Xs2)})")
    ax.set_xlabel(f"PC1 ({v1:.0f}%)", fontsize=8)
    ax.set_ylabel(f"PC2 ({v2:.0f}%)", fontsize=8)
    ax.set_title(f"PCA coverage  |  diversity={diversity:.2f}  [{len(fields)}D]", fontsize=9)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7)
    return fields


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

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    BINS = 20
    label_x = DIM_LABEL.get(dim_x, dim_x)
    label_y = DIM_LABEL.get(dim_y, dim_y)

    # Panel 1: dG histogram
    ax1.hist(dg_q, bins=BINS, color=C_QUALITY, alpha=0.6, label=f"quality ({len(dg_q)})")
    ax1.hist(dg_r, bins=BINS, color=C_REDUN,   alpha=0.6, label=f"redundancy ({len(dg_r)})")
    ax1.hist(dg_s, bins=BINS, color=C_SAMPLE,  alpha=0.9, label=f"sample ({len(dg_s)})")
    ax1.set_xlabel("dG (kcal/mol)")
    ax1.set_ylabel("Count")
    ax1.set_title("Binding affinity distribution")
    ax1.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax1.legend(fontsize=8)

    # Panel 2: dim_y histogram
    if dim_y in CAT_DIMS and y_q:
        cat_order = CAT_DIMS[dim_y]
        bins_y = range(len(cat_order) + 1)
        ax2.hist(y_q, bins=bins_y, color=C_QUALITY, alpha=0.6, label=f"quality ({len(y_q)})")
        ax2.hist(y_r, bins=bins_y, color=C_REDUN,   alpha=0.6, label=f"redundancy ({len(y_r)})")
        ax2.hist(y_s, bins=bins_y, color=C_SAMPLE,  alpha=0.9, label=f"sample ({len(y_s)})")
        _apply_cat_ticks(ax2, cat_order, axis="x")
    elif y_q:
        bins_y = range(int(min(y_q)), int(max(y_q)) + 2)
        ax2.hist(y_q, bins=bins_y, color=C_QUALITY, alpha=0.6, label=f"quality ({len(y_q)})")
        ax2.hist(y_r, bins=bins_y, color=C_REDUN,   alpha=0.6, label=f"redundancy ({len(y_r)})")
        ax2.hist(y_s, bins=bins_y, color=C_SAMPLE,  alpha=0.9, label=f"sample ({len(y_s)})")
        ax2.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax2.set_xlabel(label_y)
    ax2.set_ylabel("Count")
    ax2.set_title(f"{label_y} distribution")
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax2.legend(fontsize=8)

    # Panel 3: scatter dim_x x dim_y
    ax3.scatter(x_q, y_q, s=10, color=C_QUALITY, alpha=0.3, label=f"quality ({len(x_q)})")
    ax3.scatter(x_r, y_r, s=10, color=C_REDUN,   alpha=0.5, label=f"redundancy ({len(x_r)})")
    ax3.scatter(x_s, y_s, s=30, color=C_SAMPLE,  alpha=1.0, label=f"sample ({len(x_s)})",
                zorder=5)
    ax3.set_xlabel(label_x)
    ax3.set_ylabel(label_y)
    ax3.set_title(f"Coverage: {label_x} x {label_y}  [{sampling_method}]")
    ax3.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    if dim_y in CAT_DIMS:
        _apply_cat_ticks(ax3, CAT_DIMS[dim_y], axis="y")
    else:
        ax3.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax3.legend(fontsize=8)

    # Panel 4: PCA coverage
    pca_fields = _plot_pca(ax4, quality_rows, redundancy_rows, sample_rows)

    run_label = args.run or Path(args.run_dir).name if args.run_dir else "run"
    fig.suptitle(f"{run_label} — pool filtering and sampling", fontsize=12)
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

    if pca_fields and config_path.exists():
        with open(config_path, encoding="ascii") as fh:
            cfg = json.load(fh)
        cfg["pca_features"] = [f.lstrip("_") for f in pca_fields]
        with open(config_path, "w", encoding="ascii") as fh:
            json.dump(cfg, fh, indent=2)
        print(f"  Updated {config_path.name} with pca_features ({len(pca_fields)})")


if __name__ == "__main__":
    main()
