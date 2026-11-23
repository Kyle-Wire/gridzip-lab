"""Self-tests for codelen.lawsearch: python -m codelen.test_lawsearch (from the repo root)."""
import numpy as np
import codelen
from codelen import lawsearch as L


def test_bits_match_codelen():
    rng = np.random.default_rng(1); y = rng.normal(0.2, 0.3, 300); mu = rng.normal(0.2, 0.1, 300); sd = 0.25
    g = L.Grid(0.01, -2.5, 2.5)
    ours = L.bits_gaussian(y, mu, sd, g)
    ref = codelen.score_gaussian(y, mu, np.full(300, sd), 0.01, -2.5, 2.5, per_reading=True)["bits"]
    assert np.max(np.abs(ours - ref)) < 1e-4, np.max(np.abs(ours - ref))


def test_context_matches_codelen():
    rng = np.random.default_rng(2); n = 200
    pos = rng.uniform(0, 1, (n, 2)); y = np.sin(4 * pos[:, 0]) + rng.normal(0, 0.05, n)
    order = rng.permutation(n); g = L.Grid(0.01, -2.5, 2.5)
    ctx = np.full((n, 4), -1)
    for i, j in enumerate(order):
        prev = order[:i]
        if len(prev):
            d = np.hypot(*(pos[prev] - pos[j]).T); k = np.argsort(d)[:4]; k = k[d[k] < 0.2]
            ctx[i, :len(k)] = prev[k]
    ours = L.context_bits(y, order, ctx, g).sum()
    rank = np.empty(n, int); rank[order] = np.arange(n)       # ctx in coding-order rows for codelen
    ctx_rows = np.where(ctx >= 0, rank[np.maximum(ctx, 0)], -1)
    gg = codelen.Grid(0.01, -2.5, 2.5)
    ref = codelen.score(y[order], codelen.ContextMean(gg, ctx_rows, s0=0.3, w0=2.0), gg)["total_bits"]
    assert abs(ours - ref) < 0.05, (ours, ref)


def test_recovers_planted_family():
    rng = np.random.default_rng(3); n = 240; x = rng.uniform(0, 4, n); grp = (rng.uniform(size=n) < 0.5).astype(float)
    y = -0.1 - 0.35 * grp * np.exp(-x / 1.5) + rng.normal(0, 0.05, n)
    lib = {"G": grp, "x": x, "Gexp": L.Family(lambda a: grp * np.exp(-x / max(a, 1e-3)), 0.2, 10, log=True),
           "noise": rng.normal(size=n)}
    g = L.Grid(0.01, -2.5, 2.5)
    top, _ = L.search(lib, y, g, max_terms=2)
    assert top[0]["terms"] == ["Gexp"], top[0]["terms"]
    assert abs(top[0]["fit_alpha"]["Gexp"] - 1.5) < 0.45 and abs(top[0]["fit_coef"]["Gexp"] + 0.35) < 0.07
    v = L.verify(top[0], y, g); assert v["roundtrip_ok"]


def test_null_picks_bias_only():
    rng = np.random.default_rng(4); n = 240
    lib = {f"f{i}": rng.normal(size=n) for i in range(6)}
    y = 0.12 + rng.normal(0, 0.1, n); g = L.Grid(0.01, -2.5, 2.5)
    top, _ = L.search(lib, y, g, max_terms=2)
    assert top[0]["terms"] == [], top[0]["terms"]


def test_prescreen_matches_exhaustive():
    """The prescreened search returns the exhaustive search's top formulas on a 12-feature library with a planted law."""
    rng = np.random.default_rng(5); n = 600
    lib = {f"x{i}": rng.normal(0, 1, n) for i in range(12)}
    y = 0.1 + 0.4 * lib["x3"] - 0.25 * lib["x7"] + rng.normal(0, 0.2, n)
    g = L.Grid(0.01, -3.0, 3.0)
    full, nf = L.search(lib, y, g, max_terms=3, keep=5)
    fast, nq = L.search(lib, y, g, max_terms=3, keep=5, prescreen=40)
    assert nf == nq
    assert [r["terms"] for r in full] == [r["terms"] for r in fast], ([r["terms"] for r in full], [r["terms"] for r in fast])
    assert sorted(full[0]["terms"]) == ["x3", "x7"]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"): fn(); print("ok", name)
