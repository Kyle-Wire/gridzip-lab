"""Rescore the ground chain after CHI A11 (commit c78b4cc): same towns, cells, vocabulary and grid as ground.py.

Per cell: fix 2 = reference_check on the series, then each cell's straight-line rate re-fitted; fix 1 (sinkhole) = the
rate after the detected bend (hinge_at at the chain's own break date) for cells in the house patch (within 150 m),
standing in for the patch reading the chain reports. Post hoc: run after the fixes, against the PROTOCOL-2 baseline.
python physics/ground_rescore.py <CHI checkout at c78b4cc>
"""
import json, math, sys
from pathlib import Path
import numpy as np
CHI = Path(sys.argv[1]) / "tools" / "gate4"; sys.path.insert(0, str(CHI)); sys.path.insert(0, str(CHI.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import synthetic as S, ground_core as gc
import runner
from codelen import lawsearch as L
G = L.Grid(0.01, -40.0, 40.0); OUT = Path(__file__).resolve().parents[1] / "results" / "physics"
BREAK = 2021.9866529774038     # the chain's own break date for the sinkhole house (change.json, A11 run)


def rate(t, Y):
    A = np.column_stack([np.ones(t.size), t - t[0], np.sin(2 * math.pi * t), np.cos(2 * math.pi * t)])
    ok = np.all(np.isfinite(Y), axis=1) if Y.ndim == 2 else np.isfinite(Y)
    return np.linalg.lstsq(A[ok], Y[ok], rcond=None)[0][1]


out = {}
for name in S.SCENARIOS:
    tw = S.town(name); t = tw["t"]; valid = tw["valid"]
    ser, info = gc.reference_check(tw["series"], t, valid)
    ok = np.where(valid)[0]; idx = np.sort(np.random.default_rng(0).choice(ok, 3000, replace=False))
    Y = np.asarray(ser[idx], float).T
    vel = rate(t, Y)
    x, y = tw["x"][idx], tw["y"][idx]; dist = np.hypot(x - S.HOUSE[0], y - S.HOUSE[1]) / 1000
    if name == "sinkhole":
        for k in np.where(dist <= 0.15)[0]: vel[k] = gc.hinge_at(t, Y[:, k], BREAK)["rate_after_mm_yr"]
    k3 = t >= t[0] + 2 * (t[-1] - t[0]) / 3; h = t >= t[0] + (t[-1] - t[0]) / 2
    late = rate(t[k3], Y[k3]) - vel; accel = rate(t[h], Y[h]) - rate(t[~h], Y[~h])
    r = tw["truth_rate"][idx] - vel; c = lambda a: a - a.mean()
    lib = {"vel": c(vel), "std": c(tw["std"][idx]), "amp": c(tw["amp"][idx]), "late": late, "accel": accel, "dist": c(dist), "near": (dist <= 0.25).astype(float)}
    west, east = np.where(x < 0)[0], np.where(x >= 0)[0]
    res = runner.case(lib, r, G, [("west->east", west, east), ("east->west", east, west)])
    near = dist <= 0.1
    out[name] = {k: v for k, v in res.items() if not k.startswith("_")}
    out[name]["reference_check"] = {k: info.get(k) for k in ("result", "steps", "bad_dates")}
    out[name]["house_100m"] = {"truth": float(tw["truth_rate"][idx][near].mean()), "chain": float(vel[near].mean()), "cells": int(near.sum())}
    print(name, "one", round(res["one_correction"]["bits_per_reading"], 3), {k: round(v, 2) for k, v in res["one_correction"]["fit_coef"].items()},
          "| best", res["top"][0]["terms"], round(res["top"][0]["bits_per_reading"], 3), "| house", out[name]["house_100m"], info.get("result"))
(OUT / "ground_after_a11.json").write_text(json.dumps(out, indent=1, default=float))
