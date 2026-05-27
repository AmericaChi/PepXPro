import random
from collections import defaultdict


def diverse_sample(rows, n, dg_bins=6, pep_bins=6, seed=None):
    """6x6 quantile-bin grid on dG x peptide_length. One random entry per cell,
    remaining slots filled randomly. Adapted from scrape_pdb_bind/04_sample_and_view.py."""
    if seed is not None:
        random.seed(seed)

    n = min(n, len(rows))
    if n == len(rows):
        return list(rows)

    def quantile_bins(vals, k):
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

    cell_keys = list(cells.keys())
    random.shuffle(cell_keys)
    selected_idx = []
    for cell in cell_keys:
        if len(selected_idx) >= n:
            break
        selected_idx.append(random.choice(cells[cell]))

    if len(selected_idx) < n:
        chosen = set(selected_idx)
        remaining = [i for i in range(len(rows)) if i not in chosen]
        random.shuffle(remaining)
        selected_idx.extend(remaining[:n - len(selected_idx)])

    return [rows[i] for i in selected_idx]


def kennard_stone_sample(rows, n, seed=None):
    """Select n entries maximizing minimum pairwise distance in (dG, pep_len) space.

    Features are normalized to [0, 1] before computing distances. The first
    entry is the one closest to the centroid; subsequent entries maximize the
    minimum distance to already-selected entries (Kennard-Stone algorithm).
    """
    import math

    if seed is not None:
        random.seed(seed)

    n = min(n, len(rows))
    if n == len(rows):
        return list(rows)

    dg_vals  = [r["_dg"]      for r in rows]
    pep_vals = [r["_pep_len"] for r in rows]

    dg_min,  dg_max  = min(dg_vals),  max(dg_vals)
    pep_min, pep_max = min(pep_vals), max(pep_vals)

    def norm_dg(v):
        return (v - dg_min)  / (dg_max  - dg_min)  if dg_max  != dg_min  else 0.0

    def norm_pep(v):
        return (v - pep_min) / (pep_max - pep_min) if pep_max != pep_min else 0.0

    pts = [(norm_dg(r["_dg"]), norm_pep(r["_pep_len"])) for r in rows]

    def dist(i, j):
        return math.sqrt((pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2)

    # Start with the point closest to the centroid.
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    first = min(range(len(rows)), key=lambda i: (pts[i][0] - cx) ** 2 + (pts[i][1] - cy) ** 2)

    selected = [first]
    min_dists = [dist(first, i) for i in range(len(rows))]

    while len(selected) < n:
        next_idx = max(
            (i for i in range(len(rows)) if i not in selected),
            key=lambda i: min_dists[i],
        )
        selected.append(next_idx)
        for i in range(len(rows)):
            if i not in selected:
                d = dist(next_idx, i)
                if d < min_dists[i]:
                    min_dists[i] = d

    return [rows[i] for i in selected]


def kde_inverse_sample(rows, n, seed=None):
    """Sample n entries with probability inversely proportional to local dG density.

    Uses a Gaussian KDE on dG values. Entries in sparse regions (tails) are
    favored over entries in the crowded intermediate-affinity region.
    """
    import math

    if seed is not None:
        random.seed(seed)

    n = min(n, len(rows))
    if n == len(rows):
        return list(rows)

    dg_vals = [r["_dg"] for r in rows]
    m = len(dg_vals)

    # Silverman bandwidth
    std = (sum((v - sum(dg_vals) / m) ** 2 for v in dg_vals) / m) ** 0.5
    h = 1.06 * std * m ** (-0.2) if std > 0 else 1.0

    def kde_density(x):
        return sum(
            math.exp(-0.5 * ((x - xi) / h) ** 2) for xi in dg_vals
        ) / (m * h * math.sqrt(2 * math.pi))

    densities = [kde_density(v) for v in dg_vals]
    weights = [1.0 / d if d > 0 else 0.0 for d in densities]
    total = sum(weights)
    weights = [w / total for w in weights]

    indices = random.choices(range(m), weights=weights, k=n)
    # Deduplicate while preserving order; fill from unused if needed.
    seen = set()
    selected = []
    for i in indices:
        if i not in seen:
            seen.add(i)
            selected.append(i)

    if len(selected) < n:
        remaining = [i for i in range(m) if i not in seen]
        random.shuffle(remaining)
        selected.extend(remaining[:n - len(selected)])

    return [rows[i] for i in selected[:n]]


SAMPLING_METHODS = {
    "diverse":       diverse_sample,
    "kennard-stone": kennard_stone_sample,
    "kde":           kde_inverse_sample,
}


def run_sampling(rows, n, method="diverse", seed=None):
    if method not in SAMPLING_METHODS:
        raise ValueError(f"Unknown sampling method '{method}'. Choose from: {', '.join(SAMPLING_METHODS)}")
    return SAMPLING_METHODS[method](rows, n, seed=seed)
