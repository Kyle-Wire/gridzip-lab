"""Protocol 16: name the LBM's Case D error from the district's geometry (WindSim F3 steps 8 and 9).

python physics/district_fields.py [--seeds 10] [--runs f3s8_D,f3s9_DW,f3s9_D16]
Reads results/physics/caseD_heights.npz and results/physics/<run>_points.csv (WindSim Case D points files, raw).
Writes results/physics/district_fields.json.
Target y = ln(model / measured), mirror pairs folded (mean over +y and -y), grid 0.01.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "physics"))
from physlib import engine as E, terms as T  # noqa: E402
from physlib.engine import Term  # noqa: E402
import lateral_fields as LF  # noqa: E402

PHYS = ROOT / "results" / "physics"
B, H_LOW, U_H = 0.0625, 0.025, 6.65
RUN_SETTINGS = {"f3s8_D": {"cells_low": 4.0, "walls": 0.0}, "f3s9_DW": {"cells_low": 4.0, "walls": 1.0},
                "f3s9_D16": {"cells_low": 6.4, "walls": 0.0}}

DIST = [
    Term("d_tower", "wind", "b · d / B", "the miss grows with distance from the tower", ("d_tower",), lambda c: c["d_tower"]),
    Term("wake_side", "wind", "b · cos φ", "the miss differs downstream (+1) and upstream (−1) of the tower", ("wake_side",), lambda c: c["wake_side"]),
    Term("shelter_up", "wind", "b · H_up / H_low", "a low block just upwind shelters the point (under-resolved block drag)", ("shelter_up",), lambda c: c["shelter_up"]),
    Term("gap", "wind", "b · w_gap / B", "the open width across the flow at the point (channelling or blockage)", ("gap",), lambda c: c["gap"]),
    Term("in_canyon", "wind", "δ · [blocks both sides within 1 B]", "the point sits in a street canyon between low blocks", ("in_canyon",), lambda c: c["in_canyon"]),
    Term("tower_view", "wind", "δ · [clear line to the tower]", "the tower's downwash reaches the point unobstructed", ("tower_view",), lambda c: c["tower_view"]),
]
RUN_TERMS = [
    Term("cells_low", "wind", "b · cells per low-block height", "the low blocks are under-resolved", ("cells_low",), lambda c: c["cells_low"]),
    Term("walls", "wind", "δ · [wall model on]", "the missing wall shear on the blocks", ("walls",), lambda c: c["walls"]),
]
GEO_NAMES = [t.name for t in DIST]


def heights():
    h = np.load(PHYS / "caseD_heights.npz"); return h["heights"], h["origin_xy"], float(h["dx"])


def geometry(xy):
    hm, o, dx = heights()
    low = np.where(hm > H_LOW * 1.5, 0.0, hm)            # low-rise only (the tower is excluded)
    ix = lambda x: int(np.clip(np.floor((x - o[0]) / dx), 0, hm.shape[0] - 1))
    iy = lambda y: int(np.clip(np.floor((y - o[1]) / dx), 0, hm.shape[1] - 1))
    out = {k: [] for k in GEO_NAMES}
    for x, y in xy:
        d = np.hypot(x, y); out["d_tower"].append(d / B); out["wake_side"].append(x / d)
        i, j = ix(x), iy(y)
        i0 = ix(x - 2 * H_LOW); out["shelter_up"].append(float(low[i0:i, j].max() / H_LOW) if i > i0 else 0.0)
        col = hm[i]; up = next((k for k in range(j, hm.shape[1]) if col[k] > 0), hm.shape[1] - 1)
        dn = next((k for k in range(j, -1, -1) if col[k] > 0), 0)
        w = (up - dn) * dx; out["gap"].append(min(w, 0.5) / B)
        out["in_canyon"].append(float((up - j) * dx <= B and (j - dn) * dx <= B))
        s = np.linspace(0, 1, 200)[:-1]; px, py = x * (1 - s), y * (1 - s)
        far = np.maximum(np.abs(px), np.abs(py)) > B / 2 + dx           # outside the tower footprint
        hit = any(low[ix(a), iy(b)] > 0 for a, b in zip(px[far], py[far]))
        out["tower_view"].append(0.0 if hit else 1.0)
    return {k: np.array(v, float) for k, v in out.items()}


def load_run(name):
    p = np.genfromtxt(PHYS / f"{name}_points.csv", delimiter=",", names=True)
    key = np.round(p["x"], 4) * 1e4 + np.round(np.abs(p["y"]), 4)
    keys = np.unique(key); idx = [np.where(key == k)[0] for k in keys]
    x = np.array([p["x"][i[0]] for i in idx]); y = np.array([abs(p["y"][i[0]]) for i in idx])
    r = np.array([np.log(p["vs_pred"][i] / p["vs_obs"][i]).mean() for i in idx])
    lp = np.array([np.log(p["vs_pred"][i].mean() / U_H) for i in idx])
    return x, y, r, lp


def folds(x, d):
    near = d < 1.6
    return [("upstream -> downstream", np.where(x < 0)[0], np.where(x >= 0)[0]),
            ("downstream -> upstream", np.where(x >= 0)[0], np.where(x < 0)[0]),
            ("near -> far", np.where(near)[0], np.where(~near)[0]),
            ("far -> near", np.where(~near)[0], np.where(near)[0])]


def pair_transfers(res):
    h = {f["fold"]: f["transfers"] for f in res["holdout"]}
    return {"x_pair": bool(h.get("upstream -> downstream") and h.get("downstream -> upstream")),
            "distance_pair": bool(h.get("near -> far") and h.get("far -> near"))}


def pays(res, boot, keys):
    inc = res["guards"].get("inclusion_by_input", {}); pt = pair_transfers(res); p = res["guards"].get("null_bits_saved", {}).get("p", 1.0)
    return {k: {"inclusion": inc.get(k, 0), "bootstrap": boot.get(k, 0), "folds": pt, "null_p": p,
                "pays": bool(inc.get(k, 0) >= 0.8 and boot.get(k, 0) >= 0.6 and all(pt.values()) and p <= 0.05)} for k in keys}


def one_run(name, seeds):
    x, y, r, lp = load_run(name); g = geometry(np.column_stack([x, y])); ctx = dict(g, lp=lp)
    lib = T.WIND + DIST; fo = folds(x, g["d_tower"]); n = len(r)
    syn = []
    for s in range(seeds):
        w = LF.search(lib, ctx, np.random.default_rng(1600 + s).normal(0, r.std(), n))["winner"]["terms"]
        syn.append([t for t in w if any(k in t for k in GEO_NAMES)])
    silent = sum(not f for f in syn)
    res = LF.search(lib, ctx, r, fo, full=True, name=f"{name}: geometry terms")
    boot = LF.boot_inputs(res, r); rule = pays(res, boot, GEO_NAMES)
    named = [k for k, v in rule.items() if v["pays"]]
    reading = "untrustworthy" if silent < 9 else ("named: " + ", ".join(named) if named else "not named")
    print(name, "n", n, "syn silent", silent, "winner", res["winner"]["terms"], reading, flush=True)
    return {"n": n, "level": float(np.exp(r.mean())), "synthetic": {"silent": silent, "runs": syn, "passed": silent >= 9},
            "search": LF._sum(res), "pair_transfers": pair_transfers(res), "rule": rule, "reading": reading,
            "geometry": {k: v.tolist() for k, v in g.items()}, "x": x.tolist(), "y": y.tolist(), "r": r.tolist(),
            "mu": res["_mu"].tolist()}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--runs", default="f3s8_D,f3s9_DW,f3s9_D16"); a = ap.parse_args()
    runs = [r for r in a.runs.split(",") if (PHYS / f"{r}_points.csv").exists()]
    out = {"protocol": "physics/PROTOCOL-16.md", "runs": {}}
    for r in runs: out["runs"][r] = one_run(r, a.seeds)
    if len(runs) > 1:                                  # pooled: geometry plus run terms, stacked
        parts = [load_run(r) for r in runs]
        x = np.concatenate([p[0] for p in parts]); yy = np.concatenate([p[1] for p in parts])
        r = np.concatenate([p[2] for p in parts]); lp = np.concatenate([p[3] for p in parts])
        g = geometry(np.column_stack([x, yy])); ctx = dict(g, lp=lp)
        for k in ("cells_low", "walls"):
            ctx[k] = np.concatenate([np.full(len(p[0]), RUN_SETTINGS[n][k]) for p, n in zip(parts, runs)])
        res = LF.search(T.WIND + DIST + RUN_TERMS, ctx, r, folds(x, g["d_tower"]), full=True, name="pooled")
        boot = LF.boot_inputs(res, r)
        out["pooled"] = {"runs": runs, "search": LF._sum(res), "rule": pays(res, boot, GEO_NAMES + ["cells_low", "walls"])}
        print("pooled", res["winner"]["terms"], flush=True)
    (PHYS / "district_fields.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
