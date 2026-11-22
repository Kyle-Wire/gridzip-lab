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


if __name__ == "__main__":
    for k, f in list(globals().items()):
        if k.startswith("test_"): f(); print("ok", k)
