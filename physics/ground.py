"""PROTOCOL-2 case G: the ground chain's straight-line velocity against the planted truth in CHI's four synthetic towns."""
import json, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import runner
from codelen import lawsearch as L
CHI = Path(__file__).resolve().parents[2] / "climate-home-inspection" / "tools" / "gate4"
sys.path.insert(0, str(CHI)); sys.path.insert(0, str(CHI.parent))
import synthetic as S

OUT = Path(__file__).resolve().parents[1] / "results" / "physics"
G = L.Grid(0.01, -40.0, 40.0)


def rate(t, Y):
    A = np.column_stack([np.ones(t.size), t - t[0], np.sin(2 * math.pi * t), np.cos(2 * math.pi * t)])
    return np.linalg.lstsq(A, Y, rcond=None)[0][1]


def one(name, n=3000):
    tw = S.town(name); t = tw["t"]; ok = np.where(tw["valid"])[0]
    idx = np.sort(np.random.default_rng(0).choice(ok, n, replace=False))
    ser = tw["series"][idx].T.astype(np.float64); vel = tw["vel"][idx]
    k3 = t >= t[0] + 2 * (t[-1] - t[0]) / 3; h = t >= t[0] + (t[-1] - t[0]) / 2
    late = rate(t[k3], ser[k3]) - vel; accel = rate(t[h], ser[h]) - rate(t[~h], ser[~h])
    x, y = tw["x"][idx], tw["y"][idx]; dist = np.hypot(x - S.HOUSE[0], y - S.HOUSE[1]) / 1000
    r = tw["truth_rate"][idx] - vel
    c = lambda a: a - a.mean()
    lib = {"vel": c(vel), "std": c(tw["std"][idx]), "amp": c(tw["amp"][idx]), "late": late, "accel": accel,
           "dist": c(dist), "near": (dist <= 0.25).astype(float)}
    west, east = np.where(x < 0)[0], np.where(x >= 0)[0]
    res = runner.case(lib, r, G, [("west->east", west, east), ("east->west", east, west)])
    order = np.random.default_rng(0).permutation(n); ctx = np.full((n, 4), -1); pos = np.c_[x, y]
    for i, j in enumerate(order):
        prev = order[:i]
        if len(prev):
            d = np.hypot(*(pos[prev] - pos[j]).T); kk = np.argsort(d)[:4]; kk = kk[d[kk] <= 100]; ctx[i, :len(kk)] = prev[kk]
    res["neighbour"], _ = runner.neighbour_test(r, order, ctx, G, [np.arange(n)], n_perm=99)
    res["data"] = {"x": x.tolist(), "y": y.tolist(), "r": r.tolist(), "mu_best": res["_mu"].tolist(), "late": late.tolist(),
                   "truth": tw["truth_rate"][idx].tolist(), "vel": vel.tolist()}
    return res


if __name__ == "__main__":
    out = {}
    for name in S.SCENARIOS:
        res = one(name); out[name] = {k: v for k, v in res.items() if not k.startswith("_")}
        print(name, "one", round(res["one_correction"]["bits_per_reading"], 3), {k: round(v, 3) for k, v in res["one_correction"]["fit_coef"].items()}, "sd", res["one_correction"]["sd"])
        for t in res["top"][:4]: print("   ", t["terms"], round(t["bits_per_reading"], 3), {k: round(v, 3) for k, v in t["fit_coef"].items()}, "sd", t["sd"])
        for h in res["holdout"][:3]: print("    HO", h["terms"], h["transfers"], [(round(f["formula_bpr"], 2), round(f["one_correction_bpr"], 2)) for f in h["folds"]])
        print("    neighbour", res["neighbour"], res["coder_check"]["roundtrip_ok"])
    (OUT / "ground.json").write_text(json.dumps(out, indent=1, default=float))
