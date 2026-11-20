"""Can WHERE things sit hold data? A direct test.

Take N distinct items (words). Any ordering of them is a valid file, and
there are N! orderings, so the order alone can carry log2(N!) bits.

1. Hide a payload purely in the ORDER of the words (no bit of any word
   changes), then recover it exactly.
2. Measure what that costs: compress the sorted list (order carries nothing)
   and the payload-ordered list. If position were "free" storage, the
   second file would compress to the same size as the first.

Blocks of B words each carry floor(log2(B!)) bits via the factorial
number system (Lehmer code).
"""
import lzma, math, os, re, sys
import numpy as np

def perm_from_int(x, items):
    pool = list(items); out = []
    for k in range(len(pool), 0, -1):
        x, r = divmod(x, k)
        out.append(pool.pop(r))
    return out

def int_from_perm(perm, items):
    pool = list(items); digits = []
    for w in perm:
        r = pool.index(w); digits.append(r); pool.pop(r)
    x = 0
    for k, r in zip(range(1, len(digits) + 1), reversed(digits)):
        x = x * k + r
    return x

if __name__ == "__main__":
    src = sys.argv[1]
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 256
    text = open(src, "rb").read()[:20_000_000].decode("latin-1")
    words = sorted(set(re.findall(r"[a-z]{3,}", text.lower())))
    nblocks = len(words) // B
    words = words[: nblocks * B]
    cap = math.floor(math.log2(math.factorial(B)))     # bits per block
    payload = os.urandom(nblocks * cap // 8)
    P = int.from_bytes(payload, "big")
    mask = (1 << cap) - 1
    ordered = []
    for j in range(nblocks):
        chunk = (P >> (cap * (nblocks - 1 - j))) & mask
        ordered += perm_from_int(chunk, words[j * B:(j + 1) * B])
    # recover
    Q = 0
    for j in range(nblocks):
        Q = (Q << cap) | int_from_perm(ordered[j * B:(j + 1) * B], words[j * B:(j + 1) * B])
    ok = Q == P
    enc = lambda ws: "\n".join(ws).encode()
    xs = lambda b: len(lzma.compress(b, preset=9 | lzma.PRESET_EXTREME))
    s_sorted, s_payload = xs(enc(words)), xs(enc(ordered))
    rng = np.random.default_rng(0); shuffled = list(words); rng.shuffle(shuffled)
    s_random = xs(enc(shuffled))
    def front(ws):  # set-aware coding: shared-prefix length + new suffix
        out, prev = [], ""
        for w in ws:
            k = 0
            while k < min(len(w), len(prev)) and w[k] == prev[k]: k += 1
            out.append(chr(48 + k) + w[k:]); prev = w
        return enc(out)
    s_front = xs(front(words))
    total_bits = math.log2(math.factorial(len(words)))
    print(f"words={len(words):,} blocks={nblocks} of {B}, raw list={len(enc(words)):,} B")
    print(f"order can hold: {len(payload):,} B in blocks (whole-list maximum log2(N!)/8 = {total_bits/8:,.0f} B)")
    print(f"payload recovered exactly from word ORDER alone: {ok}")
    print(f"xz size, sorted (order holds nothing)   : {s_sorted:,} B")
    print(f"xz size, ordered to carry the payload   : {s_payload:,} B  (+{s_payload - s_sorted:,})")
    print(f"xz size, randomly shuffled              : {s_random:,} B  (+{s_random - s_sorted:,})")
    print(f"payload size                            : {len(payload):,} B")
    print(f"set-aware coding of the sorted set      : {s_front:,} B")
    print(f"  -> best way to store set + payload    : {s_front + len(payload):,} B (payload kept as plain bytes)")
    print(f"  -> set + payload hidden in the order  : {s_payload:,} B")
