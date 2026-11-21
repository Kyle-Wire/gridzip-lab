"""Header fields the decoder can work out for itself: cut them out of the chain, rebuild them after decoding.
(Extends pageid.py, which does page ids only, to three fields.)

Every enwik8 page header reads
    <title>T</title> / <id>PAGE</id> / [<restrictions>] / <revision> / <id>REV</id> / <timestamp>TS</timestamp> /
    <contributor> / <username>U</username> / <id>USER</id> (or <ip>...</ip>) ...
  page id   rises with the page's original position: rebuilt from the article permutation (which the decoder needs
            anyway) plus the gaps between consecutive ids.
  rev id    MediaWiki numbers edits in time order: predicted from the timestamp by interpolating between the pages
            already rebuilt (original order) with the nearest earlier and later timestamps; the residual is sent.
  user id   a username always has the same id: cut when the username already appeared earlier (original order) with
            that id; the decoder keeps a username -> id table. Nothing sent.
Side stream: page-id gaps and rev-id residuals, coded with a small adaptive binary arithmetic coder (Elias-gamma
shape, adaptive bucket and top mantissa bits).

  strip(chain, order)            -> (chain without those lines, side stream)
  rebuild(stripped, order, side) -> the original file
order[i] = original index of the i-th page in the chain (work/order.chain.npy). The prefix and the cut last page stay
in place (reorder.py). CLI (WSL, repo root):  headers.py test CHAIN ORDER.npy ORIGINAL"""
import calendar, re, sys, time
import numpy as np

PAGE = re.compile(rb"(  <page>\n    <title>[^\n]*</title>\n)    <id>(\d+)</id>\n")
TITLE = re.compile(rb"  <page>\n    <title>[^\n]*</title>\n")
REV = re.compile(rb"(<revision>\n)      <id>(\d+)</id>\n(      <timestamp>([^<\n]+)</timestamp>)")
REV_CUT = re.compile(rb"<revision>\n      <timestamp>([^<\n]+)</timestamp>")
USER = re.compile(rb"(        <username>([^<\n]*)</username>\n)        <id>(\d+)</id>\n")
USER_ANY = re.compile(rb"        <username>([^<\n]*)</username>\n(        <id>(\d+)</id>\n)?")


def split_pages(data):                                      # as reorder.py: prefix, whole pages, cut last page
    starts = [m.start() for m in re.finditer(rb"  <page>", data)]
    pages = [data[s:e] for s, e in zip(starts, starts[1:] + [len(data)])]
    return data[:starts[0]], pages[:-1], pages[-1]


def secs(ts):
    t = time.strptime(ts.decode(), "%Y-%m-%dT%H:%M:%SZ"); return calendar.timegm(t)


class Model:                                                # adaptive bits for Elias-gamma-shaped integers
    def __init__(self): self.p = {}
    def bits(self, slot, n):                                # n >= 0 -> list of (context, bit) and raw low bits
        v = n + 1; b = v.bit_length(); out = [((slot, "u", i), 1) for i in range(b - 1)] + [((slot, "u", b - 1), 0)]
        for j in range(b - 2, -1, -1):
            k = b - 2 - j
            out.append(((slot, "m", b, k, v >> (j + 1)) if k < 2 else None, v >> j & 1))
        return out


class Encoder:
    def __init__(self): self.x1, self.x2, self.out, self.m = 0, 0xFFFFFFFF, bytearray(), Model()
    def bit(self, ctx, y):
        p = self.m.p.get(ctx, 2048) if ctx is not None else 2048
        xm = self.x1 + ((self.x2 - self.x1) >> 12) * p
        if y: self.x2 = xm
        else: self.x1 = xm + 1
        if ctx is not None: self.m.p[ctx] = min(4095, max(1, p + ((((y << 12) - p) * 3) >> 5)))
        while (self.x1 ^ self.x2) & 0xFF000000 == 0:
            self.out.append(self.x2 >> 24); self.x1 = self.x1 << 8 & 0xFFFFFFFF; self.x2 = (self.x2 << 8 & 0xFFFFFFFF) | 255
    def int(self, slot, n):
        for ctx, y in self.m.bits(slot, n): self.bit(ctx, y)
    def finish(self): self.out += self.x1.to_bytes(4, "big"); return bytes(self.out)


class Decoder:
    def __init__(self, b): self.b, self.i, self.x1, self.x2, self.m = b + b"\0" * 4, 4, 0, 0xFFFFFFFF, Model(); self.x = int.from_bytes(self.b[:4], "big")
    def bit(self, ctx):
        p = self.m.p.get(ctx, 2048) if ctx is not None else 2048
        xm = self.x1 + ((self.x2 - self.x1) >> 12) * p
        y = int(self.x <= xm)
        if y: self.x2 = xm
        else: self.x1 = xm + 1
        if ctx is not None: self.m.p[ctx] = min(4095, max(1, p + ((((y << 12) - p) * 3) >> 5)))
        while (self.x1 ^ self.x2) & 0xFF000000 == 0:
            self.x1 = self.x1 << 8 & 0xFFFFFFFF; self.x2 = (self.x2 << 8 & 0xFFFFFFFF) | 255
            self.x = (self.x << 8 & 0xFFFFFFFF) | (self.b[self.i] if self.i < len(self.b) else 0); self.i += 1
        return y
    def int(self, slot):
        b = 1
        while self.bit((slot, "u", b - 1)): b += 1
        v = 1
        for k in range(b - 1): v = v << 1 | self.bit((slot, "m", b, k, v) if k < 2 else None)
        return v - 1


class RevPredictor:                                         # rev id from timestamp, from the pages seen so far
    def __init__(self): self.t, self.r = [], []
    def predict(self, t):
        import bisect
        if not self.t: return 0
        i = bisect.bisect_left(self.t, t); lo, hi = max(i - 1, 0), min(i, len(self.t) - 1)
        if lo == hi or self.t[hi] == self.t[lo]: return self.r[lo] if abs(self.t[lo] - t) <= abs(self.t[hi] - t) else self.r[hi]
        return self.r[lo] + (self.r[hi] - self.r[lo]) * (t - self.t[lo]) // (self.t[hi] - self.t[lo])
    def add(self, t, r):
        import bisect
        j = bisect.bisect_left(self.t, t); self.t.insert(j, t); self.r.insert(j, r)


def zz(e): return 2 * e if e >= 0 else -2 * e - 1
def unzz(n): return n // 2 if n % 2 == 0 else -(n + 1) // 2


def strip(chain, order):
    prefix, pages, tail = split_pages(chain)
    n = len(pages); order = [int(x) for x in order]
    assert len(order) == n and sorted(order) == list(range(n)), "order does not match the chain's pages"
    orig = [None] * n
    for i, p in enumerate(pages): orig[order[i]] = p
    orig.append(tail)
    enc, rp, users, prev_id, cut = Encoder(), RevPredictor(), {}, 0, [None] * (n + 1)
    for k, p in enumerate(orig):                            # original order, as the decoder will rebuild
        m = PAGE.match(p); assert m, f"page {k}: no id line after the title"
        d = m.group(2); assert str(int(d)).encode() == d and int(d) > prev_id, f"page {k}: id {d!r}"
        enc.int(0, int(d) - prev_id - 1); prev_id = int(d); q = m.group(1) + p[m.end():]
        assert not REV_CUT.search(q), f"page {k}: a revision without an id"
        r = REV.search(q)
        if r:
            rv = int(r.group(2)); assert str(rv).encode() == r.group(2)
            t = secs(r.group(4)); enc.int(1, zz(rv - rp.predict(t))); rp.add(t, rv)
            q = q[:r.start()] + r.group(1) + r.group(3) + q[r.end():]
        u = USER_ANY.search(q)
        if u:
            assert u.group(2), f"page {k}: username without an id"
            name, uid = u.group(1), u.group(3)
            if users.get(name) == uid: q = q[:u.start(2)] + q[u.end(2):]
            users[name] = uid
        cut[k] = q
    out = [prefix] + [cut[order[i]] for i in range(n)] + [cut[n]]
    return b"".join(out), enc.finish()


def rebuild(stripped, order, side):
    prefix, pages, tail = split_pages(stripped)
    n = len(pages); order = [int(x) for x in order]
    assert len(order) == n
    orig = [None] * n
    for i, p in enumerate(pages): orig[order[i]] = p
    orig.append(tail)
    dec, rp, users, prev_id, out = Decoder(side), RevPredictor(), {}, 0, [prefix]
    for k, q in enumerate(orig):
        prev_id += dec.int(0) + 1
        m = TITLE.match(q); q = q[:m.end()] + b"    <id>%d</id>\n" % prev_id + q[m.end():]
        r = REV_CUT.search(q)
        if r:
            t = secs(r.group(1)); rv = rp.predict(t) + unzz(dec.int(1)); rp.add(t, rv)
            a = r.start() + len(b"<revision>\n"); q = q[:a] + b"      <id>%d</id>\n" % rv + q[a:]
        u = USER_ANY.search(q)
        if u:
            name = u.group(1)
            if u.group(2): users[name] = u.group(3)
            else: q = q[:u.end()] + b"        <id>" + users[name] + b"</id>\n" + q[u.end():]
        out.append(q)
    return b"".join(out)


if __name__ == "__main__":
    cmd, a = sys.argv[1], sys.argv[2:]
    if cmd == "test":
        t0 = time.time()
        chain, order, orig = open(a[0], "rb").read(), np.load(a[1]), open(a[2], "rb").read()
        s, side = strip(chain, order); t1 = time.time()
        ok = rebuild(s, order, side) == orig
        print(f"chain {len(chain):,} B -> stripped {len(s):,} B (-{len(chain) - len(s):,}), side {len(side):,} B, "
              f"rebuild == original: {ok} (strip {t1 - t0:.0f} s, rebuild {time.time() - t1:.0f} s)")
        if len(a) > 3: open(a[3], "wb").write(s); open(a[3] + ".side", "wb").write(side)
        sys.exit(0 if ok else 1)
