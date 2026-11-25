"""Protocol 12: the compressor on FloodSim's planted depth cases (benchmarks/physics/depth/cases.mjs output).

Usage:
  python physics/flood_depth.py <cases.json> [--out results/physics/flood_depth.json]

Families, each searched alone (and the instrument checked against the protocol's expected terms):
  C  walled channel: y = ln(h_model / h_true_walls); context W, Wm, dx, h (model depth), S, n, Fr
  B  backwater, free edge and held edge apart: y = ln(h_model / h_true); context x_e, Lb, Fr, S, n
  R  rain sheet flow: y = ln(h_model / h_true); context x, S, n, dx
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from physlib import engine as E, guards as Q, terms as T  # noqa: E402

GRID = E.Grid(0.005, -3.0, 3.0)
G = 9.81


def fam_C(cases):
    rows = [c for c in cases if "h_model" in c]
    y = np.array([np.log(c["h_model"] / c["h_true_walls"]) for c in rows])
    ctx = {k: np.array([float(c[s]) for c in rows]) for k, s in (("W", "W"), ("Wm", "Wm"), ("dx", "dx"), ("h", "h_model"), ("S", "S"), ("n", "n"))}
    ctx["Fr"] = np.array([c["Q"] / (c["Wm"] * c["h_model"]) / np.sqrt(G * c["h_model"]) for c in rows])
    setting = np.array([f'S{c["S"]}' for c in rows])
    return ctx, y, setting, rows


def fam_B(cases, mode):
    ctx = {k: [] for k in ("x_e", "Lb", "Fr", "S", "n")}; y = []; setting = []
    for c in cases:
        if c["mode"] != mode:
            continue
        for r in c["rows"]:
            if r["h_model"] <= 1e-3:
                continue
            y.append(np.log(r["h_model"] / r["h_true"]))
            ctx["x_e"].append(r["x_from_edge"]); ctx["Lb"].append(c["Lb"]); ctx["S"].append(c["S"]); ctx["n"].append(c["n"])
            ctx["Fr"].append(c["q"] / r["h_model"] / np.sqrt(G * r["h_model"])); setting.append(f'S{c["S"]}r{c["rise"]}')
    return {k: np.array(v, float) for k, v in ctx.items()}, np.array(y), np.array(setting), None


def fam_R(cases):
    ctx = {k: [] for k in ("x", "S", "n", "dx")}; y = []; setting = []
    for c in cases:
        for r in c["rows"]:
            if r["x"] < 20 or r["h_model"] <= 0:
                continue
            y.append(np.log(r["h_model"] / r["h_true"]))
            ctx["x"].append(r["x"]); ctx["S"].append(c["S"]); ctx["n"].append(c["n"]); ctx["dx"].append(c["dx"])
            setting.append(f'S{c["S"]}n{c["n"]}i{c["mmh"]}')
    return {k: np.array(v, float) for k, v in ctx.items()}, np.array(y), np.array(setting), None


def search(name, ctx, y, setting):
    terms = [t for t in T.FLOOD_DEPTH if all(k in ctx for k in t.needs) and all(np.ptp(np.atleast_1d(ctx[k])) > 0 for k in t.needs)]
    folds = [(f"leave out {s}", np.where(setting != s)[0], np.where(setting == s)[0]) for s in sorted(set(setting.tolist()))]
    res = Q.analyse(terms, ctx, y, GRID, folds=folds, groups=setting, bootstrap_B=30, null_N=19, name=name, neighbour=False)
    gd = res["guards"]
    return {"name": name, "n": res["n"], "terms_offered": [t.name for t in terms], "winner": res["winner"]["terms"],
            "coef": res["winner"]["coef"], "bits_saved_per_reading": res["bits_saved_per_reading"],
            "sd_before": res["sd_before"], "sd_after": res["sd_after"], "mean_before": float(np.mean(y)),
            "null_p": gd.get("null_bits_saved", {}).get("p"), "verdict": Q.verdict(res), "transfers": res["transfers"],
            "_y": res["_y"].tolist(), "_mu": res["_mu"].tolist()}


EXPECT = {"C": {"raster_width", "wall_friction"}, "B_free": {"edge_drawdown"}, "B_held": set()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cases", type=Path)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "physics" / "flood_depth.json")
    a = ap.parse_args()
    D = json.loads(a.cases.read_text())
    out = {}
    if "C" in D:
        ctx, y, s, _ = fam_C(D["C"]); out["C"] = search("C walled channel", ctx, y, s)
    if "B" in D:
        for mode in ("free", "held"):
            ctx, y, s, _ = fam_B(D["B"], mode)
            if len(y) > 8:
                out[f"B_{mode}"] = search(f"B backwater, {mode} edge", ctx, y, s)
    if "R" in D:
        ctx, y, s, _ = fam_R(D["R"]); out["R"] = search("R rain sheet flow", ctx, y, s)
    for k, v in out.items():
        if k in EXPECT:
            v["instrument_check"] = {"expected": sorted(EXPECT[k]), "named": v["winner"],
                                     "pass": EXPECT[k] <= set(v["winner"]) if EXPECT[k] else (not v["winner"] or (v["null_p"] or 1) > 0.05)}
        print(k, v["winner"], v["coef"],
              "bits", round(v["bits_saved_per_reading"], 3), "null p", v["null_p"], "sd", round(v["sd_before"], 4), "->", round(v["sd_after"], 4),
              v.get("instrument_check", {}).get("pass"), flush=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
