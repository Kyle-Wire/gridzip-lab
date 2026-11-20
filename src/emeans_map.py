"""e-means (arXiv 2611.09501) vs Lloyd k-means as the grid's spatial map.

Implemented from the paper's description (Sec. 3.1-3.2): per-cluster EMA
count and sum, centroid = sum / count; mini-batch updates; step size
annealed on a cosine schedule from a peak to a floor; dead clusters revived
every 5 batches by splitting a donor; last-epoch-mean finalisation.
The paper fits its peak step and batch size with scaling laws; those
equations did not render in the HTML we read, so PEAK/FLOOR/BATCH below are
OUR choices, not the paper's.

Question: does a lower clustering loss give a better compression map?
"""
import json, math, sys
import numpy as np
from gridzip import learn_grid, map_bytes
from geometry import features, pca_coords, kmeans
import v6

PEAK, FLOOR, BATCH, EPOCHS, REVIVE_FRAC = 0.5, 0.01, 64, 3, 0.1


def assign(X, C):
    return ((X[:, None, :] - C[None]) ** 2).sum(2).argmin(1)


def emeans(X, R, seed=0):
    rng = np.random.default_rng(seed)
    n = len(X)
    C0 = X[rng.choice(n, R, replace=False)]
    cnt = np.full(R, 1e-3); sm = C0 * 1e-3            # tiny weight on the init
    steps = EPOCHS * math.ceil(n / BATCH); t = 0
    for ep in range(EPOCHS):
        perm = rng.permutation(n)
        last_lab = np.empty(n, dtype=np.int64)
        for b0 in range(0, n, BATCH):
            idx = perm[b0:b0 + BATCH]
            C = sm / cnt[:, None]
            lab = assign(X[idx], C); last_lab[idx] = lab
            alpha = FLOOR + 0.5 * (PEAK - FLOOR) * (1 + math.cos(math.pi * t / max(1, steps - 1)))
            bc = np.bincount(lab, minlength=R).astype(float)
            bs = np.zeros_like(sm); np.add.at(bs, lab, X[idx])
            cnt = (1 - alpha) * cnt + alpha * bc
            sm = (1 - alpha) * sm + alpha * bs
            t += 1
            if t % 5 == 0:                               # dead-cluster revival
                dead = np.nonzero(cnt < REVIVE_FRAC * cnt.mean())[0]
                for d in dead:
                    donor = int(cnt.argmax())
                    c = sm[donor] / cnt[donor]
                    eps = 1e-3 * rng.standard_normal(X.shape[1])
                    cnt[donor] /= 2; cnt[d] = cnt[donor]
                    sm[donor] = (c + eps) * cnt[donor]; sm[d] = (c - eps) * cnt[d]
    C = sm / cnt[:, None]
    for r in range(R):                                   # last-epoch-mean finalisation
        p = X[last_lab == r]
        if len(p): C[r] = p.mean(0)
    return assign(X, C)


def loss(X, lab):
    return float(sum(((X[lab == r] - X[lab == r].mean(0)) ** 2).sum() for r in np.unique(lab)) / len(X))


if __name__ == "__main__":
    path, n, cells, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    data = open(path, "rb").read()[:n]
    M, toks = learn_grid(data[:min(n, 2_000_000)], cells)
    table = [bytes([i]) for i in range(256)]
    for a, b in M: table.append(table[a] + table[b])
    m = 256 + len(M)
    P = pca_coords(features(toks, m))
    X = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-9)   # same space k-means uses
    maps = {"kmeans16 + kmeans128": (kmeans(P, 16), kmeans(P, 128, seed=1)),
            "emeans16 + emeans128": (emeans(X, 16), emeans(X, 128, seed=1))}
    res = []
    for name, (a, b) in maps.items():
        r = dict(n=n, map=name, loss16=loss(X, a), loss128=loss(X, b),
                 used128=int(len(np.unique(b))))
        nb, _ = v6.compress(data, table, True, True, use_space=True, reg=a, regB=b)
        r["size"] = nb + map_bytes(M) + m * (4 + 7) / 8
        res.append(r); print(json.dumps(r), flush=True)
        json.dump(res, open(out, "w"), indent=1)
