"""EM (Gaussian mixture) maps vs k-means maps.

EM gives every cell a SOFT membership over regions, with each region its own
shape (diagonal covariance). Two ways to use it:
  em_hard  : most likely region only (drop-in replacement for k-means)
  em_soft  : most likely region AND second most likely region, fed to v6's
             two map slots, so cells "between" regions carry both
Compared with k-means16 and k-means16 + k-means128 (v6) on the same grid.
"""
import json, math, sys
import numpy as np
from gridzip import learn_grid, map_bytes
from geometry import features, pca_coords, kmeans
import v6


def gmm_em(X, R, iters=60, seed=0):
    n, d = X.shape
    rng = np.random.default_rng(seed)
    mu = X[rng.choice(n, R, replace=False)].copy()
    var = np.tile(X.var(0) + 1e-6, (R, 1))
    pi = np.full(R, 1.0 / R)
    for _ in range(iters):
        # E step (log domain)
        ll = (-0.5 * (((X[:, None, :] - mu[None]) ** 2) / var[None]).sum(2)
              - 0.5 * np.log(var).sum(1)[None] + np.log(pi)[None])
        ll -= ll.max(1, keepdims=True)
        g = np.exp(ll); g /= g.sum(1, keepdims=True)
        # M step
        nk = g.sum(0) + 1e-9
        mu = (g.T @ X) / nk[:, None]
        var = (g.T @ (X ** 2)) / nk[:, None] - mu ** 2 + 1e-4
        pi = nk / n
    return g  # responsibilities, n x R


if __name__ == "__main__":
    path, n, cells, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    data = open(path, "rb").read()[:n]
    M, toks = learn_grid(data[:min(n, 2_000_000)], cells)
    table = [bytes([i]) for i in range(256)]
    for a, b in M: table.append(table[a] + table[b])
    m = 256 + len(M)
    P = pca_coords(features(toks, m))
    Xn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-9)
    g = gmm_em(Xn, 16)
    order = np.argsort(-g, axis=1)
    em1, em2 = order[:, 0], order[:, 1]
    soft = em1 * 16 + em2                      # 256 ids: (best, second) region pair
    km16, km128 = kmeans(P, 16), kmeans(P, 128, seed=1)
    print("cells with a real second region (2nd prob > 0.2):", int((np.sort(g, 1)[:, -2] > 0.2).sum()), "of", m, flush=True)
    runs = {
        "kmeans16 + kmeans128 (v6)": (km16, km128, 4 + 7),
        "em_hard16 + kmeans128": (em1, km128, 4 + 7),
        "em_soft: best16 + (best,second) pair": (em1, soft, 4 + 4),
        "kmeans16 only": (km16, km16, 4),
        "em_hard16 only": (em1, em1, 4),
    }
    res = []
    for name, (ra, rb, bits) in runs.items():
        nb, _ = v6.compress(data, table, True, True, use_space=True, reg=ra, regB=rb)
        res.append(dict(n=n, map=name, size=nb + map_bytes(M) + m * bits / 8))
        print(json.dumps(res[-1]), flush=True)
        json.dump(res, open(out, "w"), indent=1)
