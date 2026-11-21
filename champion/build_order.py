"""Step 1 of the champion: the article order and the header transform.

Articles are chained by word similarity: each page becomes a hashed tf-idf vector over its words of 3+ letters
(2^14 buckets, weight 1 + log count, unit length), and a greedy nearest-neighbour chain runs from the file's first
page. The prefix (before the first <page>) and the cut last page stay in place. The decoder needs the permutation,
which costs log2(n!) bits = 18,750 B for 12,346 pages (the same for any order).
Then headers.py cuts page, revision and contributor ids that the decoder can rebuild (side stream 16,736 B), and the
rebuild is checked against enwik8.

Writes OUT/enwik8.wordchain, OUT/order.wordchain.npy (order[k] = original index of the k-th page),
OUT/enwik8.word.hdr (what lpaq1g compresses) and OUT/enwik8.word.hdr.side.
PYTHONHASHSEED must be 0: the word hashing uses Python's hash(). About 15 s and 2 GB.
Usage (repo root):  PYTHONHASHSEED=0 python champion/build_order.py data/enwik8 champion/out"""
import math, os, re, sys, time
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import headers


def word_vectors(pages, D=1 << 14):
    W = np.zeros((len(pages), D), np.float32)
    for i, p in enumerate(pages):
        for w, c in Counter(re.findall(rb"[a-z]{3,}", p.lower())).items(): W[i, hash(w) % D] += 1 + math.log(c)
    W *= np.log(len(pages) / ((W > 0).sum(0) + 1))[None].astype(np.float32)
    return W / (np.linalg.norm(W, axis=1, keepdims=True) + 1e-9)


def greedy_chain(S):
    n = len(S); left = np.ones(n, bool); cur = 0; order = [0]; left[0] = False
    for _ in range(n - 1):
        s = S[cur].copy(); s[~left] = -9; cur = int(s.argmax()); order.append(cur); left[cur] = False
    return np.array(order)


def main(src, out):
    assert os.environ.get("PYTHONHASHSEED") == "0", "run with PYTHONHASHSEED=0 (the word hashing must be reproducible)"
    t0 = time.time(); os.makedirs(out, exist_ok=True)
    data = open(src, "rb").read()
    prefix, pages, tail = headers.split_pages(data); n = len(pages)
    W = word_vectors(pages); order = greedy_chain(W @ W.T); del W
    assert sorted(order.tolist()) == list(range(n))
    chain = prefix + b"".join(pages[i] for i in order) + tail; assert len(chain) == len(data)
    open(os.path.join(out, "enwik8.wordchain"), "wb").write(chain); np.save(os.path.join(out, "order.wordchain.npy"), order)
    s, side = headers.strip(chain, order)
    assert headers.rebuild(s, order, side) == data, "header rebuild failed"
    open(os.path.join(out, "enwik8.word.hdr"), "wb").write(s); open(os.path.join(out, "enwik8.word.hdr.side"), "wb").write(side)
    print(f"{n:,} pages, permutation {math.lgamma(n + 1) / math.log(2) / 8:,.0f} B; stripped {len(s):,} B + side "
          f"{len(side):,} B; rebuild ok ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
