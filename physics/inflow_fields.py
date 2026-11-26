"""Protocol 13: the LBM's Case C errors with an inflow that carries shear stress (WindSim F3 step 7), roof terms guarded.

python physics/inflow_fields.py [--seeds 10]
Reads results/physics/f3s7_{C1,W1}_fields.npz + points, f3s7_N1_points.csv and f3s6_R0_points.csv.
Writes results/physics/inflow_fields.json and inflow_fields_arrays.npz.
"""
import argparse, json, math, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "physics"))
from physlib import engine as E, terms as T  # noqa: E402
import lateral_fields as LF  # noqa: E402
import roof_fields as RF  # noqa: E402

PHYS = ROOT / "results" / "physics"


def pts(name):
    p = np.genfromtxt(PHYS / name, delimiter=",", names=True)
    return p[np.argsort(p["id"])]


def rule(res, boot, keys):
    inc = res["guards"].get("inclusion_by_input", {}); q = LF._quad_transfers(res)
    return {k: {"inclusion": inc.get(k, 0), "bootstrap": boot.get(k, 0), "quadrant_folds": q,
                "pays": bool(inc.get(k, 0) >= 0.8 and boot.get(k, 0) >= 0.6 and q >= 3)} for k in keys}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seeds", type=int, default=10); a = ap.parse_args()
    cC, yC, x, y, idC, ledC = RF.ctx_build("C1", "f3s7"); cW, yW, _, _, idW, ledW = RF.ctx_build("W1", "f3s7")
    R0, N1 = pts("f3s6_R0_points.csv"), pts("f3s7_N1_points.csv")
    assert (R0["id"] == idC).all() and (N1["id"] == idC).all() and (idW == idC).all()
    lr0, ln1 = np.log(R0["vs_pred"] / R0["vs_obs"]), np.log(N1["vs_pred"] / N1["vs_obs"])
    dN, dC = lr0 - ln1, lr0 - yC; s_pt = float(dN.std()); n = len(yC); fo = LF.folds(x, y)
    out = {"protocol": "physics/PROTOCOL-13.md", "n": n,
           "level": {"R0": math.exp(lr0.mean()), "N1": math.exp(ln1.mean()), "C1": math.exp(yC.mean()), "W1": math.exp(yW.mean())},
           "Q0": {"s_pt": s_pt, "noise_level_ln": float(abs(dN.mean()))}, "ledger": {"C1": ledC, "W1": ledW}}
    print("levels", out["level"], "Q0", out["Q0"], flush=True)
    FIELD, ROOF = LF.FIELD, RF.ROOF
    ctxA = {k: cC[k] for k in LF.BASE_KEYS + tuple(LF.FIELD_NAMES)}
    t = time.time(); A = LF.search(T.WIND + FIELD, ctxA, yC, fo, full=True, name="C1: FINDINGS-7 inputs + protocol 9 terms")
    print("A", A["winner"]["terms"], round(A["sd_after"], 4), f"{time.time() - t:.0f}s", flush=True)
    res = yC - A["_mu"]; sdA = float(res.std())
    roof_ctx = {k: cC[k] for k in RF.ROOF_NAMES}
    def named(yy):
        w = LF.search(ROOF, roof_ctx, yy)["winner"]["terms"]
        return w, [t for t in w if any(f in t for f in RF.ROOF_NAMES)]
    syn = {"noise_sd": sdA, "runs": []}
    for s in range(a.seeds):
        w, f = named(np.random.default_rng(1300 + s).normal(0, sdA, n)); syn["runs"].append({"winner": w, "roof": f})
    syn["silent"] = sum(not r["roof"] for r in syn["runs"]); syn["passed"] = syn["silent"] >= 9
    print("syn silent", syn["silent"], flush=True)
    syn["power"] = {}
    for nm in ("roof_stress", "roof_speed"):
        z = (cC[nm] - cC[nm].mean()) / cC[nm].std(); row = {}
        for sd in (0.02, 0.04, 0.06, 0.08, 0.10):
            row[str(sd)] = sum(any(nm in t for t in named(np.random.default_rng(1400 + s).normal(0, sdA, n) + sd * z)[1]) for s in range(a.seeds))
        syn["power"][nm] = {"found_by_sd": row, "mde_sd": next((float(k) for k, h in row.items() if h >= 8), None)}
        print("power", nm, row, flush=True)
    out["synthetic"] = syn; out["A"] = LF._sum(A)
    incA = A["guards"].get("inclusion_by_input", {})
    out["Q1"] = {"lateral_recurs": bool(incA.get("lateral", 0) >= 0.8 and LF._quad_transfers(A) >= 3),
                 "lateral_inclusion": incA.get("lateral", 0), "quadrant_folds_A": LF._quad_transfers(A)}
    if syn["passed"]:
        t = time.time(); B = LF.search(ROOF, roof_ctx, res, fo, full=True, name="C1 residual after A: roof terms")
        print("B", B["winner"]["terms"], round(B["sd_after"], 4), f"{time.time() - t:.0f}s", flush=True)
        bB = LF.boot_inputs(B, res); out["B"] = LF._sum(B); out["Q1"]["roof"] = rule(B, bB, RF.ROOF_NAMES)
        out["Q1"]["roof_pays"] = [k for k, v in out["Q1"]["roof"].items() if v["pays"]]
    else:
        out["Q1"]["roof_pays"] = None; out["Q1"]["note"] = "guarded synthetic check failed: roof terms dropped (protocol 13)"
    # Q2: the inflow change, point by point
    R0c = RF.ctx_build("R0", "f3s6")[0]
    cd = dict(cC); cd.update({f"d_{k}": cC[k] - R0c[k] for k in RF.DIFF_KEYS})
    lib = T.WIND + FIELD + (ROOF if syn["passed"] else []) + [t for t in RF.DIFF if syn["passed"] or not any(r in t.name for r in RF.ROOF_NAMES)]
    t = time.time(); C = LF.search(lib, cd, dC, fo, full=True, name="dC = ln R0 - ln C1")
    print("C", C["winner"]["terms"], round(C["sd_after"], 4), f"{time.time() - t:.0f}s", flush=True)
    bC = LF.boot_inputs(C, dC); keys = sorted({k for t in lib for k in t.needs}); nr = rule(C, bC, keys)
    p = C["guards"].get("null_bits_saved", {}).get("p", 1.0); named_C = [k for k, v in nr.items() if v["pays"]]
    structured = bool(dC.std() > 1.5 * s_pt and p <= 0.05 and named_C)
    out["Q2"] = {"dC_mean": float(dC.mean()), "dC_sd": float(dC.std()), "s_pt": s_pt, "null_p": p, "named": named_C,
                 "search": LF._sum(C), "reading": "structured" if structured else f"level shift {dC.mean():+.3f} in ln"}
    (PHYS / "inflow_fields.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez(PHYS / "inflow_fields_arrays.npz", yC=yC, yW=yW, lr0=lr0, ln1=ln1, x=x, y=y, muA=A["_mu"], muC=C["_mu"])
    print("done", json.dumps({"Q1": out["Q1"].get("roof_pays"), "Q2": out["Q2"]["reading"]}), flush=True)


if __name__ == "__main__":
    main()
