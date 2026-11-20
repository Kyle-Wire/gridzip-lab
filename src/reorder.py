"""Curriculum for a compressor: reorder enwik8's articles so related ones sit
together, then compress. (cmix ships an article order for enwik9; this tests
the idea with our own topic geometry.)

Cost of the reordering = the permutation, which the decoder needs:
log2(N!) bits for N articles (an upper bound: any ordering coded as a Lehmer
code). Prefix (before the first <page>) and suffix stay in place.

Orders tried:
  original   : as in the file
  topic      : articles sorted by their dominant topic, then by 2nd PCA axis
  chain      : greedy nearest-neighbour chain through article vectors
  random     : control
Article vector = histogram of grid cells in the article -> log/TF-IDF -> SVD.
"""
import math, os, re, sys
import numpy as np
from gridzip import learn_grid, apply_grid


def split_pages(data):
    starts = [m.start() for m in re.finditer(rb"  <page>", data)]
    prefix = data[: starts[0]]
    pages = [data[s:e] for s, e in zip(starts, starts[1:] + [len(data)])]
    # the last "page" may be cut by the 100 MB boundary: keep it last, unmoved
    return prefix, pages[:-1], pages[-1]


from numba import njit


@njit(cache=True)
def _page_cell_counts(buf, starts, G, out, m):
    V = np.zeros((len(starts), m), dtype=np.float32)
    st = 0; pg = -1
    for i in range(len(buf)):
        while pg + 1 < len(starts) and i >= starts[pg + 1]:
            pg += 1
        st = G[st, buf[i]]
        if pg >= 0:
            V[pg, out[st]] += 1.0
    return V


def article_vectors(pages, merges, m, dims=32):
    """bag of grid cells per article (longest cell ending at each byte)"""
    from v4 import build_automaton
    table = [bytes([i]) for i in range(256)]
    for a, b in merges: table.append(table[a] + table[b])
    G, out, _ = build_automaton(table)
    buf = np.frombuffer(b"".join(pages), dtype=np.uint8)
    starts = np.cumsum([0] + [len(p) for p in pages[:-1]]).astype(np.int64)
    V = _page_cell_counts(buf, starts, G, out, m)
    df = (V > 0).sum(0) + 1
    np.log1p(V, out=V); V *= np.log(len(pages) / df).astype(np.float32)[None]
    V /= np.linalg.norm(V, axis=1, keepdims=True) + 1e-9
    mu = V.mean(0)
    # top principal axes from the small m x m covariance (same subspace as an SVD, far less memory)
    Cv = (V.T @ V).astype(np.float64) / len(V) - np.outer(mu, mu)
    w, E = np.linalg.eigh(Cv)
    E = E[:, ::-1][:, :dims].astype(np.float32)
    X = V @ E - (mu @ E)[None]
    del V
    return X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-9)


def chain_order(X, start=0):
    n = len(X); left = np.ones(n, dtype=bool); order = [start]; left[start] = False
    cur = start
    for _ in range(n - 1):
        sims = X @ X[cur]; sims[~left] = -9
        cur = int(sims.argmax()); order.append(cur); left[cur] = False
    return order


if __name__ == "__main__":
    src, outdir = sys.argv[1], sys.argv[2]
    data = open(src, "rb").read()
    prefix, pages, tail = split_pages(data)
    n = len(pages)
    perm_bytes = math.lgamma(n + 1) / math.log(2) / 8
    print(f"articles={n:,}  permutation cost <= {perm_bytes:,.0f} B", flush=True)
    merges, _ = learn_grid(data[:2_000_000], 1024)
    X = article_vectors(pages, merges, 256 + len(merges))
    from geometry import kmeans
    topic = kmeans(X, 64)
    want = sys.argv[3].split(",") if len(sys.argv) > 3 else ["original", "topic", "chain", "random"]
    orders = {
        "original": list(range(n)),
        "topic": sorted(range(n), key=lambda i: (topic[i], X[i, 1])),
        "chain": chain_order(X) if "chain" in want else None,
        "random": list(np.random.default_rng(0).permutation(n)),
    }
    os.makedirs(outdir, exist_ok=True)
    for name, order in orders.items():
        if name not in want: continue
        with open(f"{outdir}/enwik8.{name}", "wb") as f:   # stream, no full copy in memory
            f.write(prefix)
            for i in order: f.write(pages[i])
            f.write(tail)
        assert os.path.getsize(f"{outdir}/enwik8.{name}") == len(data)
        np.save(f"{outdir}/order.{name}.npy", np.array(order))
    print("written", list(orders), flush=True)
