import random
from collections import defaultdict

DIM_KEYS = {
    "dg":           "_dg",
    "pep_len":      "_pep_len",
    "rec_len":      "_rec_len",
    "n_rec_chains": "_n_rec_chains",
}

AVAILABLE_DIMS = list(DIM_KEYS.keys())


def _quantile_bins(vals, k):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    bins = [0] * len(vals)
    for rank, idx in enumerate(order):
        bins[idx] = min(int(rank * k / len(vals)), k - 1)
    return bins


def _resolve_vals(rows, dims):
    result = []
    for d in dims:
        key = DIM_KEYS.get(d)
        if key is None:
            raise ValueError(f"Unknown dimension '{d}'. Available: {', '.join(AVAILABLE_DIMS)}")
        result.append([r[key] for r in rows])
    return result


def stratified_sample(rows, n, dims=None, seed=None):
    """Quantile-bin stratified sampling.

    1D: n equal-quantile bins on the single dimension, one entry per bin.
    2D: 6x6 quantile grid; one random entry per cell until n is reached.
    Only the first two dims are used if more are provided.
    Default dims: ["dg", "pep_len"].
    """
    if dims is None:
        dims = ["dg", "pep_len"]
    if len(dims) > 2:
        print(f"  Note: stratified uses at most 2 dims — ignoring {dims[2:]}.")
        dims = dims[:2]
    if seed is not None:
        random.seed(seed)

    n = min(n, len(rows))
    if n == len(rows):
        return list(rows)

    val_lists = _resolve_vals(rows, dims)

    if len(dims) == 1:
        bin_arrs = [_quantile_bins(val_lists[0], n)]
    else:
        bin_arrs = [_quantile_bins(vl, 6) for vl in val_lists]

    cells = defaultdict(list)
    for i in range(len(rows)):
        key = tuple(b[i] for b in bin_arrs)
        cells[key].append(i)

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


def kennard_stone_sample(rows, n, dims=None, seed=None):
    """Maximizes minimum pairwise distance in the given feature space (N-dimensional).

    Features are normalized to [0, 1]. The first entry is the one closest to
    the centroid; subsequent entries maximize the minimum distance to already-
    selected entries (Kennard-Stone algorithm).
    Default dims: ["dg", "pep_len"].
    """
    import math

    if dims is None:
        dims = ["dg", "pep_len"]
    if seed is not None:
        random.seed(seed)

    n = min(n, len(rows))
    if n == len(rows):
        return list(rows)

    val_lists = _resolve_vals(rows, dims)
    nd = len(dims)

    normed = []
    for vals in val_lists:
        vmin, vmax = min(vals), max(vals)
        span = vmax - vmin if vmax != vmin else 1.0
        normed.append([(v - vmin) / span for v in vals])

    pts = [tuple(normed[d][i] for d in range(nd)) for i in range(len(rows))]

    def dist(i, j):
        return math.sqrt(sum((pts[i][d] - pts[j][d]) ** 2 for d in range(nd)))

    centroid = [sum(pts[i][d] for i in range(len(pts))) / len(pts) for d in range(nd)]
    first = min(range(len(rows)),
                key=lambda i: sum((pts[i][d] - centroid[d]) ** 2 for d in range(nd)))

    selected = [first]
    selected_set = {first}
    min_dists = [dist(first, i) for i in range(len(rows))]

    while len(selected) < n:
        next_idx = max(
            (i for i in range(len(rows)) if i not in selected_set),
            key=lambda i: min_dists[i],
        )
        selected.append(next_idx)
        selected_set.add(next_idx)
        for i in range(len(rows)):
            if i not in selected_set:
                d = dist(next_idx, i)
                if d < min_dists[i]:
                    min_dists[i] = d

    return [rows[i] for i in selected]


def kde_inverse_sample(rows, n, dims=None, seed=None):
    """Inverse KDE density sampling (always 1D).

    Uses a Gaussian KDE with Silverman bandwidth on a single dimension.
    Entries in sparse regions (tails) are favored over dense intermediate
    regions. Only dims[0] is used; additional dims are ignored with a warning.
    Default dim: ["dg"].
    """
    import math

    if dims is None:
        dims = ["dg"]
    if len(dims) > 1:
        print(f"  Note: kde is 1D — using only '{dims[0]}', ignoring {dims[1:]}.")
        dims = [dims[0]]
    if seed is not None:
        random.seed(seed)

    n = min(n, len(rows))
    if n == len(rows):
        return list(rows)

    vals = _resolve_vals(rows, dims)[0]
    m = len(vals)

    std = (sum((v - sum(vals) / m) ** 2 for v in vals) / m) ** 0.5
    h = 1.06 * std * m ** (-0.2) if std > 0 else 1.0

    def kde_density(x):
        return sum(
            math.exp(-0.5 * ((x - xi) / h) ** 2) for xi in vals
        ) / (m * h * math.sqrt(2 * math.pi))

    densities = [kde_density(v) for v in vals]
    weights = [1.0 / d if d > 0 else 0.0 for d in densities]
    total = sum(weights)
    weights = [w / total for w in weights]

    indices = random.choices(range(m), weights=weights, k=n)
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


def isolation_forest_sample(rows, n, dims=None, seed=None):
    """Sample n entries weighted by isolation-forest anomaly score (N-dimensional).

    Uses sklearn IsolationForest. Entries in sparse/unusual regions of feature
    space receive higher anomaly scores and are favored. N-dimensional analogue
    of kde_inverse_sample.
    Default dims: ["dg", "pep_len"].
    """
    try:
        from sklearn.ensemble import IsolationForest
    except ImportError:
        raise ImportError("isolation-forest requires scikit-learn: conda install scikit-learn")

    if dims is None:
        dims = ["dg", "pep_len"]
    if seed is not None:
        random.seed(seed)

    n = min(n, len(rows))
    if n == len(rows):
        return list(rows)

    val_lists = _resolve_vals(rows, dims)
    X = list(zip(*val_lists))  # shape (n_rows, n_dims)

    clf = IsolationForest(random_state=seed if seed is not None else 0, contamination="auto")
    clf.fit(X)
    # score_samples: more negative = more anomalous
    scores = clf.score_samples(X)
    # invert so anomalous entries get higher weight
    weights = [-s for s in scores]
    total = sum(weights)
    weights = [w / total for w in weights]

    indices = random.choices(range(len(rows)), weights=weights, k=n)
    seen = set()
    selected = []
    for i in indices:
        if i not in seen:
            seen.add(i)
            selected.append(i)

    if len(selected) < n:
        remaining = [i for i in range(len(rows)) if i not in seen]
        random.shuffle(remaining)
        selected.extend(remaining[:n - len(selected)])

    return [rows[i] for i in selected[:n]]


def cluster_sample(rows, n, dims=None, seed=None):
    """K-means clustering with k=n; selects the entry nearest each cluster centroid.

    Guarantees exactly n entries with uniform coverage of the feature space.
    Deterministic for a fixed seed. Works in N dimensions.
    Default dims: ["dg", "pep_len"].
    """
    try:
        from sklearn.cluster import KMeans
    except ImportError:
        raise ImportError("cluster requires scikit-learn: conda install scikit-learn")

    if dims is None:
        dims = ["dg", "pep_len"]

    n = min(n, len(rows))
    if n == len(rows):
        return list(rows)

    val_lists = _resolve_vals(rows, dims)
    nd = len(dims)

    # Normalize to [0, 1] before clustering
    normed = []
    for vals in val_lists:
        vmin, vmax = min(vals), max(vals)
        span = vmax - vmin if vmax != vmin else 1.0
        normed.append([(v - vmin) / span for v in vals])

    X = list(zip(*normed))

    km = KMeans(n_clusters=n, random_state=seed if seed is not None else 0, n_init="auto")
    km.fit(X)

    import math
    # For each cluster, find the entry closest to its centroid
    selected = []
    seen = set()
    for label, centroid in enumerate(km.cluster_centers_):
        cluster_idx = [i for i, lbl in enumerate(km.labels_) if lbl == label]
        if not cluster_idx:
            continue
        closest = min(cluster_idx,
                      key=lambda i: math.sqrt(sum((X[i][d] - centroid[d]) ** 2
                                                  for d in range(nd))))
        if closest not in seen:
            seen.add(closest)
            selected.append(closest)

    # Fill if any cluster was empty or duplicate
    if len(selected) < n:
        remaining = [i for i in range(len(rows)) if i not in seen]
        random.shuffle(remaining)
        selected.extend(remaining[:n - len(selected)])

    return [rows[i] for i in selected[:n]]


SAMPLING_METHODS = {
    "stratified":        stratified_sample,
    "kennard-stone":     kennard_stone_sample,
    "kde":               kde_inverse_sample,
    "isolation-forest":  isolation_forest_sample,
    "cluster":           cluster_sample,
}


def effective_ndims(method, dims):
    """Number of dimensions actually used by the method."""
    if method == "kde":
        return 1
    if method == "stratified":
        return min(len(dims), 2)
    return len(dims)  # kennard-stone, isolation-forest, cluster: full N-D


def method_label(method, dims):
    """Display label including dimensionality, e.g. 'stratified-2D' or 'kde-1D'."""
    return f"{method}-{effective_ndims(method, dims)}D"


def run_sampling(rows, n, method="stratified", dims=None, seed=None):
    if method not in SAMPLING_METHODS:
        raise ValueError(f"Unknown method '{method}'. Choose from: {', '.join(SAMPLING_METHODS)}")
    return SAMPLING_METHODS[method](rows, n, dims=dims, seed=seed)
