"""Surprise score self-tests: python -m surprise.test_surprise (from the repo root, CPU, about half a minute)."""
import numpy as np
from surprise import core


def test_null_is_calibrated_and_roundtrips():
    rng = np.random.default_rng(1); U, S = 200, 4
    xy = rng.uniform(0, 300, (U, 2)); y = 0.1 * rng.standard_normal((U, S))
    r = core.score(y, xy, 60, 0.02, -5, 5, verify=True)
    assert r["roundtrip_ok"]
    p95 = np.percentile(r["surprise"], 95)
    assert -2 < np.median(r["surprise"]) < 1 and p95 < 7, (np.median(r["surprise"]), p95)


def test_planted_house_stands_out():
    rng = np.random.default_rng(2); U, S = 200, 4
    xy = rng.uniform(0, 300, (U, 2)); y = (rng.random((U, 1)) < 0.3) * 1.0 + 0.1 * rng.standard_normal((U, S)); y[5] += 2.0   # wetter than either mode
    r = core.score(y, xy, 60, 0.02, -5, 5, verify=False)
    assert r["surprise"][5] > np.percentile(r["surprise"], 95) and r["sign"][5] > 0   # a bimodal block has other odd houses too


def test_height_explains_a_low_house():
    rng = np.random.default_rng(3); U, S = 200, 3
    xy = rng.uniform(0, 300, (U, 2)); g = rng.normal(0, 0.3, U)
    y = -0.8 * g[:, None] + 0.05 * rng.standard_normal((U, S))   # lower ground, wetter
    r = core.score(y, xy, 60, 0.02, -5, 5, descriptor=g, verify=False)
    low = np.argmin(g)
    assert r["surprise_given"][low] < r["surprise"][low]


if __name__ == "__main__":
    for f in (test_null_is_calibrated_and_roundtrips, test_planted_house_stands_out, test_height_explains_a_low_house):
        f(); print("PASS", f.__name__)
