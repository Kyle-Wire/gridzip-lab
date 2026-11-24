"""Self-tests for physlib: python -m physlib.test_physlib (from the repo root, ~30 s)."""
import warnings
import numpy as np
from physlib import engine as E, synth as S, envtests as V
warnings.filterwarnings("ignore")


def _run(dom, law, seed, noise="iid"):
    terms, ctx, y, g, planted = S.case(dom, law, seed, noise=noise)
    return E.run_case(terms, ctx, y, E.Grid(*g), prev=S.prev_of(dom, len(y))), planted


def test_planted_found():
    for dom, law in (("ground", "sinkhole onset 2022"), ("flood", "weir drowns to orifice"), ("wind", "contrast + entrance decay")):
        r, planted = _run(dom, law, 0)
        assert sorted(r["winner"]["terms"]) == sorted(planted), (dom, r["winner"]["terms"])
    r, _ = _run("ground", "sinkhole onset 2022", 0)
    assert abs(r["winner"]["alpha"]["onset"] - 2022.0) < 0.15


def test_noise_silent():
    for dom in ("heat", "ground", "flood"):
        for noise in ("iid", "ar1"):
            r, _ = _run(dom, None, 2000, noise)
            assert r["winner"]["terms"] == [], (dom, noise, r["winner"]["terms"])


def test_coder_roundtrip():
    r, _ = _run("flood", "weir drowns to orifice", 1)
    assert r["coder_check"]["roundtrip_ok"]


def test_env_farmland():
    out = V.one(("Irrigated farmland (heat)", "flag", 0)); assert out["named"]
    out = V.one(("Irrigated farmland (heat)", "control", 0)); assert not out["flagged"]


def test_v7_grammar_finds_product():
    from physlib import grammar as G, terms as T
    rng = np.random.default_rng(90123); ctx = S.wind_ctx(rng); n = len(ctx["lp"])
    y = 0.1 + 0.3 * (ctx["hc"] * ctx["lp"] - (ctx["hc"] * ctx["lp"]).mean()) + rng.normal(0, 0.15, n)
    r = E.run_case(list(T.WIND) + G.candidates(ctx, n), ctx, y, E.Grid(0.01, -4, 4), robust=True)
    assert any("hc" in t and "lp" in t for t in r["winner"]["terms"]), r["winner"]["terms"]


def test_v7_grammar_noise_silent():
    from physlib import grammar as G
    terms, ctx, y, g, _ = S.case("wind", None, 91000)
    r = E.run_case(list(terms) + G.candidates(ctx, len(y)), ctx, y, E.Grid(*g), robust=True)
    assert r["winner"]["terms"] == [], r["winner"]["terms"]


def test_v7_structure_bits_unchanged_without_grammar():
    terms, ctx, y, g, _ = S.case("wind", None, 1)
    b = E.assemble(terms, ctx, len(y)); import math
    assert abs(E.structure_bits(b[:2], len(b)) - (E.M.delta_bits(6) + 5 * math.log2(len(b) + 2))) < 1e-9


if __name__ == "__main__":
    for k, f in list(globals().items()):
        if k.startswith("test_"): f(); print("ok", k)
