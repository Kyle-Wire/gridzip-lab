"""PROTOCOL-4 part B: the physics compressor on every model's errors, before and after the 22 Nov fixes.

python -m physlib.cases <case ...>      cases: w1 w2 heat flood ground (default: all)
Inputs (paths default to sibling checkouts; override with env vars):
  results/physics/*.json, aij_c_arrays.npz  the PROTOCOL 1-3 data (WindSim AIJ C points, Layer 1 pairs, flood before)
  HEAT_NC   folder with fixes2_phl_noanchor_{stamp,hour_ending}_*.nc (HeatSim M1, Kyle's PC)
  HEAT_ISD  folder with the three stations' 2025 ISD csvs (AWS noaa-global-hourly-pds)
  FLOOD_AFTER  FloodSim results dir of the run at 8f525e9 (PV-SS10, PV-SS13)
  CHI       climate-home-inspection checkout with A11 (main)
Writes results/physlib/<case>.json; arrays for figures in results/physlib/<case>_arrays.npz.
"""
import json, math, os, sys
from pathlib import Path

import numpy as np

from . import engine as E, terms as T

ROOT = Path(__file__).resolve().parents[1]
PHYS = ROOT / "results" / "physics"
OUT = ROOT / "results" / "physlib"
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "physics"))
from codelen import lawsearch as LS  # noqa: E402


def neighbour(y, order, ctx, g, groups, n_perm=99, seed=0):
    rng = np.random.default_rng(seed)
    obs = LS.context_bits(y, order, ctx, g).sum() / len(y); null = np.zeros(n_perm)
    for k in range(n_perm):
        yp = y.copy()
        for grp in groups: yp[grp] = rng.permutation(y[grp])
        null[k] = LS.context_bits(yp, order, ctx, g).sum() / len(y)
    return {"bits_per_reading": float(obs), "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "structure_left_bits": float(null.mean() - obs), "p": float((1 + np.sum(null <= obs)) / (n_perm + 1))}


def spatial_ctx(pos, radius, k=4, seed=0, blocks=None):
    n = len(pos); order = np.random.default_rng(seed).permutation(n); ctx = np.full((n, k), -1)
    for i, j in enumerate(order):
        prev = order[:i]
        if blocks is not None: prev = prev[blocks[prev] == blocks[j]]
        if len(prev):
            d = np.hypot(*(pos[prev] - pos[j]).T); kk = np.argsort(d)[:k]; kk = kk[d[kk] <= radius + 1e-9]; ctx[i, :len(kk)] = prev[kk]
    return order, ctx


def finish(res, extra=None, remainder=None):
    out = E.public(res)
    if remainder is not None: out["remainder_structure"] = remainder
    if extra: out.update(extra)
    return out


def arrays(res, **kw):
    return {"y": res["_y"], "mu": res["_mu"], "mu0": res["_mu0"], **kw}


# ================================================================================================ W1 AIJ Case C
def w1():
    import aij_c as A
    d = np.load(PHYS / "aij_c_arrays.npz", allow_pickle=True); out, arr = {}, {}
    g = E.Grid(0.01, -4.0, 4.0)
    for eng in ("lbm", "rans"):
        x, y, m, p, hc, lay = (d[f"{eng}_{k}"] for k in ("x", "y", "m", "p", "hc", "layout"))
        r = np.log(m / p)
        sx = np.array([A.street(v) for v in x]); sy = np.array([A.street(v) for v in y])
        dx = np.maximum(np.abs(x) - 0.1, 0); dy = np.maximum(np.abs(y) - 0.1, 0)
        ctx = {"lp": np.log(p / A.U_APP), "fetch": (x + 0.5) / A.D, "along": (sy & ~sx).astype(float), "cross": (sx & ~sy).astype(float),
               "hc": hc.astype(float), "corner": ((np.hypot(dx, dy) <= 0.1 + 1e-9) & (hc > 0)).astype(float)}
        folds = [(f"leave out {L}", np.where(lay != L)[0], np.where(lay == L)[0]) for L in ("0D", "1D", "2D")]
        res = E.run_case(T.WIND, ctx, r, g, folds, name=f"AIJ Case C {eng.upper()} D/10")
        blocks = np.array([{"0D": 0, "1D": 1, "2D": 2}[v] for v in lay])
        order, cx = spatial_ctx(np.c_[x, y], 0.075, blocks=blocks)
        rem = res["_y"] - res["_mu"]
        out[eng] = finish(res, remainder={"before": neighbour(r - res["_mu0"], order, cx, g, [np.where(blocks == b)[0] for b in range(3)]),
                                          "after": neighbour(rem, order, cx, g, [np.where(blocks == b)[0] for b in range(3)])})
        out[eng]["noise_floor_sd"] = float(json.loads((PHYS / "aij_c.json").read_text())["noise_floor"]["sd_log"])
        arr[eng] = arrays(res, x=x, yy=y, layout=lay, fetch=ctx["fetch"], along=ctx["along"])
        print("W1", eng, out[eng]["winner"]["terms"], round(out[eng]["bits_saved_per_reading"], 3), out[eng]["transfers"])
    return out, arr


# ================================================================================================ W2 Layer 1
def w2():
    L1 = json.loads((PHYS / "layer1.json").read_text()); out, arr = {}, {}
    g = E.Grid(0.005, -3.0, 3.0)
    for q, lab in (("Q1_thunderstorm", "thunderstorm"), ("Q2_ordinary", "ordinary")):
        D = L1[q]["data"]; y = np.array(D["y"], float); n = len(y)
        st = np.array(D["station"]); pl = np.array(D["place"])
        key = np.char.add(np.char.add(pl, "|"), st)
        loo = np.zeros(n)
        for k in set(key.tolist()):
            ii = np.where(key == k)[0]
            if len(ii) > 1: loo[ii] = (y[ii].sum() - y[ii]) / (len(ii) - 1)
        ctx = {"stn_loo": loo, "dz": np.array(D["dz_elev"], float), "dirn": np.array(D["dirn"], float),
               "hour": np.array(D["hl"], float), "month": np.array(D["month"], float)}
        if "z0" in D: ctx.update(z=np.full(n, 10.0), z0=np.array(D["z0"], float))
        folds = [(f"leave out {p}", np.where(pl != p)[0], np.where(pl == p)[0]) for p in sorted(set(pl.tolist()))]
        res = E.run_case(T.WIND, ctx, y, g, folds, name=f"Layer 1 {lab}")
        out[lab] = finish(res); arr[lab] = arrays(res, place=pl, station=st)
        print("W2", lab, out[lab]["winner"]["terms"], round(out[lab]["bits_saved_per_reading"], 3), out[lab]["transfers"])
    out["after"] = "pending: WindSim F5 (station z0 from its own gusts) has a protocol and code (WindSim 01d1864) but no result files yet"
    return out, arr


# ================================================================================================ H HeatSim
STATIONS = {"72030464752": "Wings Field", "72408013739": "PHL airport", "72408594732": "NE PHL airport"}
HDAYS = ("2025-07-08", "2025-07-25", "2025-07-28", "2025-07-29", "2025-07-30")


def heat_rows(ncdir, isddir, variant):
    import pandas as pd, xarray as xr
    from pyproj import Transformer
    import heat as H
    to_utm = Transformer.from_crs("EPSG:4326", "EPSG:32618", always_xy=True).transform
    obs = {s: H.isd(Path(isddir) / f"{s}.csv") for s in STATIONS}
    rows = []
    for day in HDAYS:
        ds = xr.open_dataset(Path(ncdir) / f"fixes2_phl_noanchor_{variant}_{day}.nc"); assert int(ds.attrs["anchored"]) == 0
        times = pd.DatetimeIndex(ds.time.values)
        for sid in STATIONS:
            f, lat, lon = obs[sid]; sx, sy = to_utm(lon, lat)
            ix = int(np.argmin(np.abs(ds.x.values - sx))); iy = int(np.argmin(np.abs(ds.y.values - sy)))
            tm = ds.open_site_temperature.values[:, iy, ix].astype(float); o = f.reindex(times)
            fr = {k: float(ds[k].values[iy, ix]) for k in ("frac_low", "frac_roof", "frac_canyon", "frac_tree")}
            dTm = np.gradient(tm)
            for k, t in enumerate(times):
                rows.append(dict(day=day, sid=sid, t=t, Tm=tm[k], dTm=dTm[k], To=o["T"].iloc[k], Td=o["Td"].iloc[k], W=o["W"].iloc[k],
                                 cloud=o["cloud"].iloc[k], lat=lat, lon=lon, **fr))
        ver = ds.attrs["heatsim_version"]; ds.close()
    df = pd.DataFrame(rows); df["cloud"] = df["cloud"].fillna(0.5)
    df = df.dropna(subset=["To", "Tm"]).reset_index(drop=True)
    df["W"] = df["W"].fillna(df["W"].median()); df["Td"] = df["Td"].fillna(df["Td"].median())
    return df, ver


def heat():
    import pandas as pd
    import heat as H
    ncdir = os.environ.get("HEAT_NC", "/mnt/user-data/uploads/Climate Inspection/runs/heatsim")
    isddir = os.environ.get("HEAT_ISD", "/tmp/claude-0/isd")
    g = E.Grid(0.1, -20.0, 20.0); out, arr = {}, {}
    for variant, lab in (("stamp", "before"), ("hour_ending", "after")):
        df, ver = heat_rows(ncdir, isddir, variant)
        r = (df["To"] - df["Tm"]).values
        lat, lon = df["lat"].values, df["lon"].values
        sun_at = lambda tau: np.maximum(H.sun_sin(df["t"] - pd.to_timedelta(tau, unit="h"), lat, lon), 0)
        hl = ((df["t"].dt.hour + df["t"].dt.minute / 60 - 4) % 24).values
        ctx = {"sun": sun_at(0.0), "sun_at": sun_at, "Tm": df["Tm"].values, "dTm": df["dTm"].values,
               "dayid": (df["day"] + df["sid"]).values, "W": np.maximum(df["W"].values, 0.5), "cloud": df["cloud"].values,
               "Td": df["Td"].values, "fveg": (df["frac_tree"] + df["frac_low"]).values, "fbuilt": (df["frac_roof"] + df["frac_canyon"]).values,
               "station": df["sid"].values, "hour": hl}
        folds = [(f"leave out {d}", np.where(df["day"] != d)[0], np.where(df["day"] == d)[0]) for d in HDAYS] + \
                [(f"leave out {STATIONS[s]}", np.where(df["sid"] != s)[0], np.where(df["sid"] == s)[0]) for s in STATIONS]
        res = E.run_case(T.HEAT, ctx, r, g, folds, name=f"HeatSim {ver} {variant}")
        keyt = df["t"].astype(str).values; order = np.random.default_rng(0).permutation(len(r)); cx = np.full((len(r), 2), -1)
        for i, j in enumerate(order):
            same = [q for q in np.where(keyt == keyt[j])[0] if q != j][:2]; cx[i, :len(same)] = same
        groups = [np.where(df["sid"] == s)[0] for s in STATIONS]
        rem = {"before": neighbour(r - res["_mu0"], order, cx, g, groups), "after": neighbour(r - res["_mu"], order, cx, g, groups)}
        dr = df.groupby(["day", "sid"]).agg(mo=("Tm", lambda v: v.max() - v.min()), ob=("To", lambda v: v.max() - v.min()))
        out[lab] = finish(res, extra={"heatsim_version": ver, "sun_timing": variant, "rmse_K": float(np.sqrt(np.mean(r ** 2))),
                                      "bias_K": float(np.mean(r)), "range_ratio": float((dr["mo"] / dr["ob"]).mean())}, remainder=rem)
        arr[lab] = arrays(res, hour=hl, sid=df["sid"].values, day=df["day"].values, Tm=df["Tm"].values, To=df["To"].values)
        print("H", lab, out[lab]["winner"]["terms"], {k: v for k, v in out[lab]["winner"]["coef"].items()}, out[lab]["winner"]["alpha"],
              round(out[lab]["bits_saved_per_reading"], 3), out[lab]["transfers"], "rmse", round(out[lab]["rmse_K"], 2))
    return out, arr


# ================================================================================================ F FloodSim
FS = Path(os.environ.get("FLOODSIM", ROOT.parent / "FloodSim")) / "benchmarks" / "physics"


def _ss13(path):
    d = json.loads(Path(path).read_text()); spec = json.loads((FS / "spec/PV-SS13.json").read_text())["parameters"]["cases"]["value"]
    ids = list(spec); hand = d["hand"]; rep = {k.split(" ")[0]: v for k, v in d["cases"].items() if "[cpu]" in k}
    col = lambda f: np.array([f(i) for i in ids], float)
    return dict(ids=ids, V=col(lambda i: hand[i]["velocityMs"]), T=col(lambda i: hand[i]["spreadM"]), Q=col(lambda i: hand[i]["flowReachingTheGrateLs"]),
                Eo=col(lambda i: hand[i]["frontalShare"]), E=col(lambda i: hand[i]["efficiency"]), share=col(lambda i: rep[i]["share"]),
                h=col(lambda i: rep[i]["depthInGrateCellM"]), SL=col(lambda i: spec[i]["sceneArgs"]["SL"]), rain=col(lambda i: spec[i]["rain"]["mmhr"]))


def flood():
    g = E.Grid(0.005, -3.0, 3.0); out, arr = {}, {}
    F1 = json.loads((PHYS / "flood.json").read_text())["F1_sheffield"]["data"]
    h = np.array(F1["h_mm"]) / 1000; meas = np.array(F1["measured_ls"]); Q1 = np.array(F1["Q1_ls"]) / 1000
    after_dir = Path(os.environ.get("FLOOD_AFTER", sorted((FS / "results").iterdir())[-1]))
    rows = [r for r in json.loads((after_dir / "PV-SS10.json").read_text())["rows"] if r["quantity"].startswith("exchange")]
    sim_after = np.array([r["sim"] for r in rows]); ref_after = np.array([r["ref"] for r in rows])
    assert np.allclose(ref_after, meas, atol=1e-3), "PV-SS10 tests do not line up"
    ctx = {"h": h, "Fr": (Q1 / h) / np.sqrt(9.81 * h), "Q": Q1}
    odd, even = np.arange(len(h))[::2], np.arange(len(h))[1::2]
    folds = [("odd -> even", odd, even), ("even -> odd", even, odd)]
    for lab, sim in (("before", np.array(F1["floodsim_ls"])), ("after", sim_after)):
        r = np.log(meas / sim)
        res = E.run_case(T.flood_terms(ctx), ctx, r, g, folds, name=f"Sheffield rig {lab}")
        out[f"F1_{lab}"] = finish(res, extra={"cw": 0.5 if lab == "before" else 0.32, "median_ratio": float(np.exp(np.median(r)))})
        arr[f"F1_{lab}"] = arrays(res, h=h, meas=meas, sim=sim)
        print("F1", lab, out[f"F1_{lab}"]["winner"]["terms"], out[f"F1_{lab}"]["winner"]["coef"], round(out[f"F1_{lab}"]["bits_saved_per_reading"], 3))
    before_dir = FS / "results" / "2026-11-17T22-43-34-587Z"
    for lab, path in (("before", before_dir / "PV-SS13.json"), ("after", after_dir / "PV-SS13.json")):
        D = _ss13(path); r = np.log(D["share"] / D["E"])
        ctx = {k: D[k] for k in ("V", "T", "Q", "Eo", "SL", "rain", "h")}
        folds = [(f"leave out {R:.0f} mm/h", np.where(D["rain"] != R)[0], np.where(D["rain"] == R)[0]) for R in (50, 100, 200)]
        res = E.run_case(T.flood_terms(ctx), ctx, r, g, folds, max_terms=3, name=f"Grate on a grade {lab}")
        out[f"F2_{lab}"] = finish(res, extra={"run": path.parent.name, "median_share_ratio": float(np.exp(np.median(r)))})
        arr[f"F2_{lab}"] = arrays(res, SL=D["SL"], V=D["V"], rain=D["rain"])
        print("F2", lab, out[f"F2_{lab}"]["winner"]["terms"], out[f"F2_{lab}"]["winner"]["coef"], round(out[f"F2_{lab}"]["bits_saved_per_reading"], 3))
    out["F2_note"] = ("the PV-SS13 rerun at FloodSim 8f525e9 (capture 'hec22' default) gave the same shares as the 17 Nov run to 3 "
                      "decimals, so 'after' does not yet carry the grate fix; reported to the fixes thread")
    return out, arr


# ================================================================================================ G ground
def ground():
    chi = Path(os.environ.get("CHI", ROOT.parent / "climate-home-inspection")) / "tools" / "gate4"
    sys.path.insert(0, str(chi)); sys.path.insert(0, str(chi.parent))
    import synthetic as S, ground_core as gc
    BREAK = 2021.9866529774038       # the chain's own break date for the sinkhole house (CHI change.json, A11)
    gr = E.Grid(0.01, -40.0, 40.0); gs = E.Grid(0.01, -200.0, 200.0); out, arr = {}, {}

    def fitline(t, Y):
        A = np.column_stack([np.ones(t.size), t - t[0], np.sin(2 * math.pi * t), np.cos(2 * math.pi * t)])
        ok = np.all(np.isfinite(Y), axis=1) if Y.ndim == 2 else np.isfinite(Y)
        c = np.linalg.lstsq(A[ok], np.nan_to_num(Y[ok]), rcond=None)[0]; return c[1], A @ c

    for name in S.SCENARIOS:
        tw = S.town(name); t = tw["t"]; valid = tw["valid"]
        ok = np.where(valid)[0]; idx = np.sort(np.random.default_rng(0).choice(ok, 3000, replace=False))
        x, y = tw["x"][idx], tw["y"][idx]; dist = np.hypot(x - S.HOUSE[0], y - S.HOUSE[1]) / 1000
        truth = tw["truth_rate"][idx]
        for lab in ("before", "after"):
            ser = tw["series"]
            if lab == "after": ser, info = gc.reference_check(tw["series"], t, valid)
            Y = np.asarray(ser[idx], float).T
            vel, fitY = fitline(t, Y)
            if lab == "after" and name == "sinkhole":
                for k in np.where(dist <= 0.15)[0]:
                    hb = gc.hinge_at(t, Y[:, k], BREAK); vel[k] = hb["rate_after_mm_yr"]
                    A = np.column_stack([np.ones(t.size), np.minimum(t - BREAK, 0), np.maximum(t - BREAK, 0), np.sin(2 * math.pi * t), np.cos(2 * math.pi * t)])
                    okk = np.isfinite(Y[:, k]); c = np.linalg.lstsq(A[okk], Y[okk, k], rcond=None)[0]; fitY[:, k] = A @ c
            # (i) rate level: truth - chain rate
            k3 = t >= t[0] + 2 * (t[-1] - t[0]) / 3; h2 = t >= t[0] + (t[-1] - t[0]) / 2
            late = fitline(t[k3], Y[k3])[0] - vel; accel = fitline(t[h2], Y[h2])[0] - fitline(t[~h2], Y[~h2])[0]
            ctx = {"late": late, "accel": accel, "vel": vel, "std": tw["std"][idx], "amp": tw["amp"][idx],
                   "near": (dist <= 0.25).astype(float), "dist": dist}
            r = truth - vel
            folds = [("west -> east", np.where(x < 0)[0], np.where(x >= 0)[0]), ("east -> west", np.where(x >= 0)[0], np.where(x < 0)[0])]
            res = E.run_case(T.GROUND_RATE, ctx, r, gr, folds, name=f"ground {name} rate {lab}")
            order, cx = spatial_ctx(np.c_[x, y], 100.0)
            rem = {"before": neighbour(r - res["_mu0"], order, cx, gr, [np.arange(len(r))], n_perm=49),
                   "after": neighbour(r - res["_mu"], order, cx, gr, [np.arange(len(r))], n_perm=49)}
            near = dist <= 0.1
            key = f"{name}_{lab}"
            out[key] = {"rate": finish(res, remainder=rem, extra={"house_100m": {"truth": float(truth[near].mean()), "chain": float(vel[near].mean())}})}
            # (ii) series level: what the chain's per-cell model leaves in the data
            R = Y - fitY
            net = np.nanmedian(R, axis=1)
            sctx = {"t": t}
            half = np.arange(len(t)) < len(t) // 2
            sf = [("first half -> second", np.where(half)[0], np.where(~half)[0]), ("second half -> first", np.where(~half)[0], np.where(half)[0])]
            rn = E.run_case(T.ground_terms(sctx), sctx, net, gs, sf, name=f"ground {name} network {lab}")
            out[key]["network"] = finish(rn)
            rng = np.random.default_rng(0)
            close = np.where(dist <= 0.3)[0]; close = rng.choice(close, min(60, len(close)), replace=False)
            rest = rng.choice(np.setdiff1d(np.arange(len(idx)), close), 60, replace=False)
            cells = np.r_[close, rest]; wins = []
            for k in cells:
                yk = R[:, k]; okk = np.isfinite(yk)
                if okk.sum() < 50: continue
                c2 = {"t": t[okk]}
                rk = E.run_case(T.ground_terms(c2), c2, yk[okk], gs, (), max_terms=3, beam=4)
                wins.append({"cell": int(k), "dist_km": float(dist[k]), "terms": rk["winner"]["terms"], "alpha": rk["winner"]["alpha"],
                             "bits_saved": rk["bits_saved_per_reading"]})
            out[key]["cells"] = wins
            arr[key] = {"t": t, "net": net, "net_mu": rn["_mu"], "x": x, "y": y, "r": r, "mu": res["_mu"], "dist": dist}
            print("G", key, "rate", out[key]["rate"]["winner"]["terms"], "| network", rn["winner"]["terms"], rn["winner"]["alpha"],
                  "| cells near with onset", sum("onset" in w["terms"] for w in wins if w["dist_km"] <= 0.15), "/", sum(w["dist_km"] <= 0.15 for w in wins))
    return out, arr


CASES = {"w1": w1, "w2": w2, "heat": heat, "flood": flood, "ground": ground}

if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    for c in (sys.argv[1:] or list(CASES)):
        res, arr = CASES[c]()
        (OUT / f"{c}.json").write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
        flat = {f"{k}__{kk}": np.asarray(vv) for k, v in arr.items() for kk, vv in v.items()}
        np.savez_compressed(OUT / f"{c}_arrays.npz", **flat)
