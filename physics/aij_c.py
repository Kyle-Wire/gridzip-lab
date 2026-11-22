"""Missing physics in WindSim's AIJ Case C street speeds, by code length (protocol: physics/PROTOCOL.md).

python physics/aij_c.py [--windsim ../WindSim] [--out results/physics]
Reads WindSim's committed per-point files; writes aij_c.json and the figures. CPU only, a few minutes.
"""
import argparse, json, math, subprocess, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from codelen import lawsearch as L

D = 0.2
U_APP = 2.43           # measured approach speed at z = 0.02 m, AF_caseC.csv (WindSim docs/validation.md)
G = L.Grid(0.01, -4.0, 4.0)
HC = {"0D": 0.0, "1D": 1.0, "2D": 2.0}


def load(path, col):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["x"], d["y"], d["vs_obs"], d[col], (d["pred_with_k"] if "pred_with_k" in d.dtype.names else None)


def street(v): a = round(abs(v), 3); return 0.1 < a < 0.3
def edge(v): return round(abs(v), 3) in (0.1, 0.3)


def features(x, y, hc, pred, predk=None):
    """The library fixed in the protocol, from geometry and the model's output only."""
    sx = np.array([street(v) for v in x]); sy = np.array([street(v) for v in y]); ex = np.array([edge(v) for v in x])
    A = (sy & ~sx).astype(float); X = (sx & ~sy).astype(float); I = (sx & sy).astype(float)
    mouth = (sy & ex).astype(float)
    xi = (x + 0.5) / D
    lee = np.where(sx & ~sy, (x - np.sign(x) * 0.2) / 0.1, 0.0)
    off = np.where(sy & ~sx, np.abs(np.abs(y) - 0.2) / 0.1, np.where(sx & ~sy, np.abs(np.abs(x) - 0.2) / 0.1, 0.0))
    dx = np.maximum(np.abs(x) - 0.1, 0); dy = np.maximum(np.abs(y) - 0.1, 0)
    near_c = ((np.hypot(dx, dy) <= 0.1 + 1e-9) & (hc > 0)).astype(float)
    lp = np.log(pred / U_APP)
    ex_ = lambda m: (lambda a: m * np.exp(-xi / max(a, 1e-3)))
    lib = {"A": A, "X": X, "I": I, "mouth": mouth, "xi": xi, "exp_xi": L.Family(ex_(1.0), 0.25, 16, log=True),
           "lee": lee, "off": off, "hc": hc, "near_c": near_c, "lp": lp}
    for nm, m in (("A", A), ("X", X), ("I", I), ("near_c", near_c)):
        lib[f"{nm}*xi"] = m * xi; lib[f"{nm}*exp_xi"] = L.Family(ex_(m), 0.25, 16, log=True)
        lib[f"{nm}*hc"] = m * hc; lib[f"{nm}*lp"] = m * lp
    if predk is not None:
        q = (predk / pred) ** 2 - 1
        lib["kfrac"] = L.Family(lambda a: 0.5 * np.log1p(a * q), 0.0, 3.0, fixed=True)
    exclusive = [{"xi", "exp_xi"}] + [{f"{m}*xi", f"{m}*exp_xi"} for m in ("A", "X", "I", "near_c")]
    cls = np.where(A > 0, "along", np.where(X > 0, "cross", "intersection"))
    return lib, exclusive, cls


def dataset(root, files, col="vs_pred"):
    xs, ys, ms, ps, pks, hcs, lay = [], [], [], [], [], [], []
    for layout, f in files.items():
        x, y, m, p, pk = load(root / f, col)
        xs += [x]; ys += [y]; ms += [m]; ps += [p]; hcs += [np.full(len(x), HC[layout])]; lay += [np.full(len(x), layout)]
        if pk is not None: pks += [pk]
    cat = np.concatenate
    return dict(x=cat(xs), y=cat(ys), m=cat(ms), p=cat(ps), pk=cat(pks) if pks else None, hc=cat(hcs), layout=cat(lay))


def mirror_floor(ds):
    """sd of ln m(x, y) - ln m(x, -y) / sqrt 2 over every layout: measurement repeatability + tunnel asymmetry."""
    d = []
    for lay in np.unique(ds["layout"]):
        s = ds["layout"] == lay; key = {(round(a, 3), round(b, 3)): v for a, b, v in zip(ds["x"][s], ds["y"][s], ds["m"][s])}
        d += [math.log(key[(a, b)] / key[(a, -b)]) for (a, b) in key if b > 0 and (a, -b) in key]
    sd = float(np.std(d) / math.sqrt(2))
    ent = math.log2(sd * math.sqrt(2 * math.pi * math.e) / G.res)      # bits per reading of N(0, sd) on the grid
    return sd, ent, len(d)


def spatial_test(ds, r, n_perm=199, seed=0, radius=0.075):
    """Model-free: are residuals predictable from already-coded neighbours? Against shuffles over the positions."""
    rng = np.random.default_rng(seed); out = {}
    tot_obs, tot_null = 0.0, np.zeros(n_perm)
    for lay in np.unique(ds["layout"]):
        idx = np.where(ds["layout"] == lay)[0]; pos = np.c_[ds["x"][idx], ds["y"][idx]]
        order = rng.permutation(len(idx)); ctx = np.full((len(idx), 4), -1)
        for i, j in enumerate(order):
            prev = order[:i]
            if len(prev):
                dd = np.hypot(*(pos[prev] - pos[j]).T); k = np.argsort(dd)[:4]; k = k[dd[k] <= radius + 1e-9]
                ctx[i, :len(k)] = prev[k]
        rr = r[idx]
        obs = L.context_bits(rr, order, ctx, G).sum(); tot_obs += obs
        for t in range(n_perm): tot_null[t] += L.context_bits(rng.permutation(rr), order, ctx, G).sum()
    n = len(r)
    return {"bits_per_reading": tot_obs / n, "null_mean": float(tot_null.mean() / n), "null_sd": float(tot_null.std() / n),
            "p": float((1 + np.sum(tot_null <= tot_obs)) / (n_perm + 1)), "null": (tot_null / n).tolist()}


def bias_only_bits(r): return L.score({"_": np.zeros(len(r))}, [], r, G)


def run_engine(name, ds, max_terms=3, keep=25):
    r = np.log(ds["m"] / ds["p"]); lib, excl, cls = features(ds["x"], ds["y"], ds["hc"], ds["p"], ds["pk"])
    t0 = time.time(); top, n_cand = L.search(lib, r, G, max_terms=max_terms, exclusive=excl, keep=keep)
    none = L.score(lib, [], r, G)
    zero_bits = float(L.bits_gaussian(r, 0.0, max(float(np.sqrt(np.mean(r ** 2))), G.res / 2), G).sum())
    best = top[0]; ver = L.verify(best, r, G)
    # hold-out by layout, top 5
    hold = []
    for res in top[:5]:
        row = {"terms": res["terms"], "layouts": {}}
        for lay in ("0D", "1D", "2D"):
            te = np.where(ds["layout"] == lay)[0]; tr = np.where(ds["layout"] != lay)[0]
            if len(te) == 0 or len(tr) == 0: continue
            b_f = L.holdout(lib, res["terms"], r, G, tr, te); b_0 = L.holdout(lib, [], r, G, tr, te)
            row["layouts"][lay] = {"formula_bits_per_reading": b_f / len(te), "bias_only_bits_per_reading": b_0 / len(te)}
        row["transfers"] = bool(row["layouts"]) and all(v["formula_bits_per_reading"] < v["bias_only_bits_per_reading"]
                                                        for v in row["layouts"].values())
        hold.append(row)
    strip = lambda d: {k: v for k, v in d.items() if not k.startswith("_") and k not in ("rounded", "alphas")}
    return {"engine": name, "n": len(r), "candidates": n_cand, "search_seconds": round(time.time() - t0, 1),
            "no_correction_bits_per_reading": zero_bits / len(r), "bias_only": strip(none), "top": [strip(t) for t in top],
            "winner_coder_check": ver, "holdout": hold, "_r": r, "_mu": best["_mu"], "_mu0": none["_mu"], "_cls": cls,
            "_lib": lib}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--windsim", default=str(Path(__file__).resolve().parents[2] / "WindSim"))
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "results" / "physics"))
    ap.add_argument("--perms", type=int, default=199)
    a = ap.parse_args(); ws = Path(a.windsim); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    v = ws / "docs" / "validation"
    commit = subprocess.run(["git", "-C", str(ws), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    lbm = dataset(v / "lbm_case_c_d10_blind", {k: f"points_{k}.csv" for k in HC})
    rans = dataset(v / "rans_case_c_d10", {k: f"points_{k}.csv" for k in HC}, col="pred_vector")
    res = {"windsim_commit": commit, "grid": {"res": G.res, "lo": G.lo, "hi": G.hi, "eps": G.eps}}

    # step 0: sanity on planted laws, same geometry
    rng = np.random.default_rng(0)
    lib, excl, _ = features(lbm["x"], lbm["y"], lbm["hc"], lbm["p"])
    planted = -0.1 - 0.35 * lib["A"] * np.exp(-lib["xi"] / 1.5) + rng.normal(0, 0.08, len(lbm["x"]))
    t, _ = L.search(lib, planted, G, 3, exclusive=excl, keep=3)
    null = 0.12 + rng.normal(0, 0.15, len(lbm["x"])); tn, _ = L.search(lib, null, G, 3, exclusive=excl, keep=3)
    U, W = rans["p"], rans["pk"]; synth_m = np.sqrt(U ** 2 + (2 / 3) * (W ** 2 - U ** 2)) * np.exp(rng.normal(0, 0.05, len(U)))
    libk, exk, _ = features(rans["x"], rans["y"], rans["hc"], U, W)
    tk, _ = L.search(libk, np.log(synth_m / U), G, 2, exclusive=exk, keep=3)
    s0 = {"planted": {"top": t[0]["terms"], "fit_coef": t[0]["fit_coef"], "fit_alpha": t[0]["fit_alpha"]},
          "null": {"top": tn[0]["terms"]},
          "kfrac": {"top": tk[0]["terms"], "fit_alpha": tk[0]["fit_alpha"]}}
    ok_a = t[0]["terms"] == ["A*exp_xi"] and abs(t[0]["fit_alpha"]["A*exp_xi"] / 1.5 - 1) < 0.3 \
        and abs(t[0]["fit_coef"]["A*exp_xi"] / -0.35 - 1) < 0.2
    ok_b = tn[0]["terms"] == []
    ok_c = "kfrac" in tk[0]["terms"] and abs(tk[0]["fit_alpha"]["kfrac"] - 2 / 3) < 0.15
    s0["pass"] = {"planted": bool(ok_a), "null": bool(ok_b), "kfrac": bool(ok_c)}
    res["step0_sanity"] = s0; print("step 0", json.dumps(s0["pass"]), s0["planted"], s0["kfrac"], flush=True)

    # step 1: structure and the noise floor
    floor = mirror_floor(lbm)
    res["noise_floor"] = {"sd_log": floor[0], "bits_per_reading": floor[1], "pairs": floor[2]}
    for nm, ds in (("lbm", lbm), ("rans", rans)):
        r = np.log(ds["m"] / ds["p"]); st = spatial_test(ds, r, a.perms)
        res[f"step1_{nm}"] = {k: v for k, v in st.items() if k != "null"}; res[f"_null_{nm}"] = st["null"]
        print("step 1", nm, res[f"step1_{nm}"], flush=True)

    # steps 2-3: search + hold-out, both engines
    eng = {}
    for nm, ds in (("lbm", lbm), ("rans", rans)):
        e = run_engine(nm, ds); eng[nm] = e
        print("step 2", nm, e["top"][0]["terms"], round(e["top"][0]["bits_per_reading"], 3),
              "bias-only", round(e["bias_only"]["bits_per_reading"], 3), flush=True)
    # RANS read as sqrt(U^2 + 2k) (the bracket's top), for comparison
    rk = dict(rans); rk["p"] = rans["pk"]; rk["pk"] = None
    eng["rans_with_k"] = run_engine("rans_with_k", rk)

    # step 4: resolution, LBM layout 1D: refit the LBM winner's terms at each N
    rdir = v / "lbm_case_c_1D_resolution"
    series = [("D/10", v / "lbm_case_c_d10_blind" / "points_1D.csv")] + [
        (f"D/{n}{' z0' if 'z0' in d else ''}", rdir / d / "points_1D.csv") for d, n in
        (("d16", 16), ("d16_z0tunnel", 16), ("d20", 20), ("d24", 24), ("d32_z0tunnel", 32))]
    win = eng["lbm"]["top"][0]["terms"]; resol = []
    for lab, f in series:
        x, y, m, p, _ = load(f, "vs_pred"); hc = np.ones(len(x)); lb, _, cls = features(x, y, hc, p)
        r = np.log(m / p); s = L.score(lb, win, r, G); s_b = L.score(lb, [], r, G)
        resol.append({"run": lab, "fit_coef": s["fit_coef"], "fit_alpha": s["fit_alpha"],
                      "bits_per_reading": s["bits_per_reading"], "bias_only_bits_per_reading": s_b["bits_per_reading"],
                      "mean_r_by_street": {c: float(np.mean(r[cls == c])) for c in ("along", "cross", "intersection")}})
    res["step4_resolution"] = resol

    for nm, e in eng.items(): res[f"step2_{nm}"] = {k: v for k, v in e.items() if not k.startswith("_")}
    (out / "aij_c.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez(out / "aij_c_arrays.npz", **{f"{nm}_{k[1:]}": e[k] for nm, e in eng.items() for k in ("_r", "_mu", "_mu0")},
             **{f"{nm}_{c}": d[c] for nm, d in (("lbm", lbm), ("rans", rans)) for c in ("x", "y", "m", "p", "hc")},
             lbm_layout=lbm["layout"], rans_layout=rans["layout"], cls=eng["lbm"]["_cls"],
             null_lbm=np.array(res.pop("_null_lbm")), null_rans=np.array(res.pop("_null_rans")))
    (out / "aij_c.json").write_text(json.dumps(res, indent=1, default=float))
    print("wrote", out)


if __name__ == "__main__":
    main()
