"""Which GEOMETRY of the grid predicts best?

Every cell gets a region id; v5 uses the regions as context. Only the map
changes between rows, so the compressed size measures how much each geometry
knows about the data (map cost included: log2(R) bits per cell).

  kmeans    : spherical k-means on PCA coordinates (v5's map)
  pca_grid  : a literal coordinate grid in PCA space: top axes cut into bins
  laplacian : Laplacian eigenmaps of the cell-transition graph (manifold
              learning), then k-means on those coordinates
  random    : control
"""
import json, math, os, sys
import numpy as np
from gridzip import learn_grid, map_bytes
import v5


def features(tokens, m):
    t = np.asarray(tokens, dtype=np.int64)
    A = np.zeros((m, m), dtype=np.float64)
    np.add.at(A, (t[:-1], t[1:]), 1.0)
    return A


def pca_coords(A, dims=16):
    F = np.log1p(np.hstack([A, A.T]))
    F /= np.linalg.norm(F, axis=1, keepdims=True) + 1e-9
    U, S, _ = np.linalg.svd(F - F.mean(0), full_matrices=False)
    return U[:, :dims] * S[:dims]


def kmeans(X, R, seed=0, iters=30):
    X = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-9)
    rng = np.random.default_rng(seed)
    C = X[rng.choice(len(X), R, replace=False)]
    for _ in range(iters):
        lab = (X @ C.T).argmax(1)
        for r in range(R):
            p = X[lab == r]
            if len(p): c = p.sum(0); C[r] = c / (np.linalg.norm(c) + 1e-9)
    return lab


def pca_grid(X, axes, bins):
    """quantile bins on each of the top `axes` principal axes -> one id"""
    lab = np.zeros(len(X), dtype=np.int64)
    for a in range(axes):
        q = np.quantile(X[:, a], np.linspace(0, 1, bins + 1)[1:-1])
        lab = lab * bins + np.searchsorted(q, X[:, a])
    return lab


def laplacian_coords(A, dims=16):
    W = np.log1p(A + A.T)
    d = W.sum(1) + 1e-9
    Dm = 1 / np.sqrt(d)
    L = np.eye(len(W)) - (Dm[:, None] * W * Dm[None, :])
    vals, vecs = np.linalg.eigh(L)
    return vecs[:, 1:dims + 1] * Dm[:, None]   # skip the trivial eigenvector


if __name__ == "__main__":
    path, n, cells, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    data = open(path, "rb").read()[:n]
    M, toks = learn_grid(data[:min(n, 2_000_000)], cells)
    table = [bytes([i]) for i in range(256)]
    for a, b in M: table.append(table[a] + table[b])
    m = 256 + len(M)
    A = features(toks, m)
    P = pca_coords(A)
    Lc = laplacian_coords(A)
    maps = {
        "kmeans16 (v5)": kmeans(P, 16),
        "pca_grid 2 axes x 4 bins": pca_grid(P, 2, 4),
        "pca_grid 4 axes x 2 bins": pca_grid(P, 4, 2),
        "pca_grid 3 axes x 3 bins (27)": pca_grid(P, 3, 3),
        "laplacian kmeans16": kmeans(Lc, 16),
        "laplacian grid 4 axes x 2": pca_grid(Lc, 4, 2),
        "random16": np.random.default_rng(1).integers(0, 16, m),
    }
    res = []
    for name, reg in maps.items():
        R = int(reg.max()) + 1
        nb, _ = v5.compress(data, table, True, True, use_space=True, reg=reg)
        size = nb + map_bytes(M) + m * math.log2(max(R, 2)) / 8
        res.append(dict(n=n, map=name, regions=R, size=size))
        print(json.dumps(res[-1]), flush=True)
        json.dump(res, open(out, "w"), indent=1)
