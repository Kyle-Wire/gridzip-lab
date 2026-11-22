"""PROTOCOL-2 case F: FloodSim inlet exchange. F1 Sheffield rig (PV-SS10), F2 grate on a grade (PV-SS13)."""
import json, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import runner
from codelen import lawsearch as L

FS = Path(__file__).resolve().parents[2] / "FloodSim" / "benchmarks" / "physics"
OUT = Path(__file__).resolve().parents[1] / "results" / "physics"; OUT.mkdir(parents=True, exist_ok=True)
G = L.Grid(0.005, -3.0, 3.0)


def f1():
    s = json.loads((FS / "spec/PV-SS10.json").read_text())["parameters"]
    T = s["tests"]["value"]; D = s["inletDiameterM"]["value"]; Cw = s["appCoefficients"]["value"]["Cw"]
    h = np.array([t["hsTableM"] for t in T]); Q = np.array([t["QeRawLs"] for t in T]) / 1000; Q1 = np.array([t["Q1RawLs"] for t in T]) / 1000
    model = Cw * math.pi * D * math.sqrt(2 * 9.81) * h ** 1.5
    r = np.log(Q / model)
    Fr = (Q1 / (1.0 * h)) / np.sqrt(9.81 * h)          # rig 1 m wide (protocol); approach Froude number
    lib = {"lnh": np.log(h / h.mean()), "lnQ1": np.log(Q1 / Q1.mean()), "Fr": Fr - Fr.mean(),
           "datum": L.Family(lambda d: 1.5 * np.log((h + d) / h), -0.005, 0.010, fixed=True)}
    odd = np.arange(len(h))[::2]; even = np.arange(len(h))[1::2]
    res = runner.case(lib, r, G, [("odd->even", odd, even), ("even->odd", even, odd)], max_terms=2)
    res["data"] = {"h_mm": (1000 * h).tolist(), "measured_ls": (1000 * Q).tolist(), "floodsim_ls": (1000 * model).tolist(),
                   "Q1_ls": (1000 * Q1).tolist(), "r": r.tolist(), "fr": Fr.tolist(), "corr_lnh_lnQ1": float(np.corrcoef(np.log(h), np.log(Q1))[0, 1])}
    return res


def f2():
    d = json.loads((FS / "results/2026-11-17T22-43-34-587Z/PV-SS13.json").read_text())
    spec = json.loads((FS / "spec/PV-SS13.json").read_text())["parameters"]["cases"]["value"]
    ids = list(spec); hand = d["hand"]
    rep = {k.split(" ")[0]: v for k, v in d["cases"].items() if "[cpu]" in k}
    V = np.array([hand[i]["velocityMs"] for i in ids]); Tt = np.array([hand[i]["spreadM"] for i in ids])
    Qg = np.array([hand[i]["flowReachingTheGrateLs"] for i in ids]); Eo = np.array([hand[i]["frontalShare"] for i in ids])
    E = np.array([hand[i]["efficiency"] for i in ids]); share = np.array([rep[i]["share"] for i in ids])
    dg = np.array([rep[i]["depthInGrateCellM"] for i in ids]); SL = np.array([spec[i]["sceneArgs"]["SL"] for i in ids])
    rain = np.array([spec[i]["rain"]["mmhr"] for i in ids])
    r = np.log(share / E)
    c = lambda a: a - a.mean()
    lib = {"V": c(V), "T": c(Tt), "lnQ": c(np.log(Qg)), "Eo": c(Eo), "SL": c(SL), "d": c(1000 * dg)}
    folds = [(f"leave out {R} mm/h", np.where(rain != R)[0], np.where(rain == R)[0]) for R in (50, 100, 200)]
    res = runner.case(lib, r, G, folds, max_terms=2)
    res["data"] = {"id": ids, "V": V.tolist(), "spread": Tt.tolist(), "Q_ls": Qg.tolist(), "Eo": Eo.tolist(), "hec22": E.tolist(),
                   "floodsim": share.tolist(), "SL": SL.tolist(), "rain": rain.tolist(), "r": r.tolist(),
                   "corr": {k: float(np.corrcoef(v, r)[0, 1]) for k, v in lib.items()}}
    return res


if __name__ == "__main__":
    out = {}
    for name, fn in (("F1_sheffield", f1), ("F2_grate_on_grade", f2)):
        res = fn(); out[name] = {k: v for k, v in res.items() if not k.startswith("_")}
        print(name, "one", round(res["one_correction"]["bits_per_reading"], 3), res["one_correction"]["fit_coef"])
        for t in res["top"][:5]: print("   ", t["terms"], round(t["bits_per_reading"], 3), {k: round(v, 3) for k, v in t["fit_coef"].items()}, t["fit_alpha"])
        for h in res["holdout"]: print("    HO", h["terms"], h["transfers"], [(round(f["formula_bpr"], 2), round(f["one_correction_bpr"], 2)) for f in h["folds"]])
        print("   data", {k: v for k, v in res["data"].items() if k in ("corr", "corr_lnh_lnQ1")})
    (OUT / "flood.json").write_text(json.dumps(out, indent=1, default=float))
