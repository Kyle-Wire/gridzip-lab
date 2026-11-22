"""PROTOCOL-3: WindSim Layer 1 (regional storm-gust chain) errors, after and before the F2 hurricane split.

python physics/layer1.py <WindSim PR 10 checkout> <dir with a3_blind/*.json> <dir with ISD header csvs>
"""
import json, math, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import runner
from codelen import lawsearch as L
from codelen import mdl as M

OUT = Path(__file__).resolve().parents[1] / "results" / "physics"
G = L.Grid(0.005, -3.0, 3.0)
EXCLUDE = ()
CA = {"los-angeles", "riverside", "san-francisco"}


def stations_meta(hdr):
    m = {}
    for f in Path(hdr).glob("*.csv"):
        d = pd.read_csv(f, nrows=1, usecols=["STATION", "LATITUDE", "LONGITUDE", "ELEVATION"], dtype={"STATION": str})
        m[f.stem] = (float(d.LATITUDE[0]), float(d.LONGITUDE[0]), float(d.ELEVATION[0]))
    return m


def km(a, b):
    p = np.radians; la1, lo1, la2, lo2 = a[0], a[1], b[..., 0], b[..., 1]
    h = np.sin(p(la2 - la1) / 2) ** 2 + np.cos(p(la1)) * np.cos(p(la2)) * np.sin(p(lo2 - lo1) / 2) ** 2
    return 6371 * 2 * np.arcsin(np.sqrt(h))


def load(files, group, meta):
    rows = []
    for f in files:
        d = json.loads(Path(f).read_text())
        if not d.get("pairs"): continue
        ids = {s["icao"]: s["id"] for s in d["stations"]}
        reg = {s["time"]: s["regional_gust_ms"] for s in d["storms"]}
        for p in d["pairs"]:
            la, lo, el = meta[ids[p["station"]]]
            rows.append(dict(place=d["region"], group=group, station=p["station"], storm=p["storm_time"], type=p["type"],
                             obs=p["observed_ms"], pred=p["predicted_ms"], dirn=p["direction_deg"], z0=p["z0_m"],
                             kzt=p["kzt"], near_km=p["nearest_km"], n_nb=p["n_neighbours"], reg=reg.get(p["storm_time"], np.nan),
                             lat=la, lon=lo, elev=el / 1000))
    return pd.DataFrame(rows)


def features(df):
    df = df.reset_index(drop=True).copy()
    df["y"] = np.log(df.obs / df.pred)
    t = pd.to_datetime(df.storm)
    df["month"] = t.dt.month; df["year"] = t.dt.year
    hl = (t.dt.hour + t.dt.minute / 60 + df.lon / 15) % 24
    dz = np.zeros(len(df)); sl = np.zeros(len(df)); tl = np.zeros(len(df))
    for _, g in df.groupby(["place", "storm", "type"]):
        ix = g.index.values; P = g[["lat", "lon"]].values
        for i in ix:
            o = ix[ix != i]
            if len(o) == 0: continue
            dd = np.maximum(km(df.loc[i, ["lat", "lon"]].values.astype(float), df.loc[o, ["lat", "lon"]].values.astype(float)), 1.0)
            w = 1 / dd ** 2; dz[i] = df.elev[i] - np.sum(w * df.elev[o]) / w.sum()
            sl[i] = df.y[o].mean()
    for _, g in df.groupby(["place", "station", "type"]):
        ix = g.index.values
        for i in ix:
            o = ix[ix != i]; tl[i] = df.y[o].mean() if len(o) else 0.0
    df["dz_elev"] = dz; df["storm_loo"] = sl; df["stn_loo"] = tl
    df["hl"] = hl
    medreg = df.groupby("place").reg.transform("median"); df["lnreg"] = np.log(df.reg / medreg).fillna(0)
    return df


def library(df):
    c = lambda a: np.asarray(a, float) - np.mean(a)
    d = np.radians(df.dirn.values)
    ca = df.place.isin(CA).values.astype(float)
    off = ca * ((df.dirn.values >= 0) & (df.dirn.values <= 135)).astype(float)
    return {"lnz0": c(np.log(df.z0)), "kzt": c(df.kzt), "near_km": c(df.near_km), "n_nb": c(df.n_nb),
            "dz_elev": c(df.dz_elev), "elev": c(df.elev),
            "dir_c": L.Family(lambda a: np.cos(d - np.radians(a)), 0.0, 345.0, n=24),
            "offshore_ca": c(off), "lnpred": c(np.log(df.pred)),
            "msin": np.sin(2 * np.pi * df.month / 12), "mcos": np.cos(2 * np.pi * df.month / 12),
            "hsin": np.sin(2 * np.pi * df.hl / 24), "hcos": np.cos(2 * np.pi * df.hl / 24),
            "lnreg": c(df.lnreg), "yr": c(df.year - 2012), "storm_loo": c(df.storm_loo), "stn_loo": c(df.stn_loo),
            "ca": c(ca), "desert": c((df.place == "phoenix").values.astype(float))}


def folds(df, cal_blind=True):
    f = [(f"leave out {p}", np.where(df.place != p)[0], np.where(df.place == p)[0]) for p in sorted(df.place.unique())]
    if cal_blind and (df.group == "blind").any():
        f.append(("calibration -> blind", np.where(df.group == "cal")[0], np.where(df.group == "blind")[0]))
    return f


def spread_test(y, mu, groups, g):
    """Bits of the residuals under one pooled sd against one sd per group (each sd paid as codelen.mdl does)."""
    r = y - mu
    def sd_bits(sd0):
        p = max(0, -math.floor(math.log2(max(sd0, g.res)))) + 3; b, s = M.param_bits(sd0, p); return b, max(s, g.res / 2)
    b, s = sd_bits(float(np.sqrt(np.mean(r ** 2)))); pooled = b + float(L.bits_gaussian(y, mu, s, g).sum())
    tot = 0.0; per = {}
    for k in np.unique(groups):
        m = groups == k; b, s = sd_bits(float(np.sqrt(np.mean(r[m] ** 2))))
        tot += b + float(L.bits_gaussian(y[m], mu[m], s, g).sum()); per[str(k)] = {"n": int(m.sum()), "sd": s}
    return {"pooled_bpr": pooled / len(y), "per_group_bpr": tot / len(y), "saving_bpr": (pooled - tot) / len(y), "groups": per}


def neighbours(df, key, n_perm=99):
    """Context = the earlier-coded pairs sharing `key` (same storm, or same station); shuffles within each place."""
    y = df.y.values; n = len(y); order = np.random.default_rng(0).permutation(n); ctx = np.full((n, 8), -1)
    grp = df[key].astype(str).agg("|".join, axis=1).values
    for i, j in enumerate(order):
        same = [q for q in np.where(grp == grp[j])[0] if q != j][:8]; ctx[i, :len(same)] = same
    res, _ = runner.neighbour_test(y, order, ctx, G, [np.where(df.place == p)[0] for p in df.place.unique()], n_perm=n_perm)
    return res


def run_case(df, name, extra_terms=None, exclude=None):
    exclude = () if "--first" in sys.argv else ("storm_loo",) if exclude is None else exclude
    lib = library(df)
    for t in EXCLUDE + tuple(exclude): lib.pop(t, None)
    if extra_terms: lib.update(extra_terms)
    fo = folds(df)
    res = runner.case(lib, df.y.values, G, fo, exclusive=[{"ca", "offshore_ca"}, {"desert", "ca"}], n_hold=8)
    for h in res["holdout"]:
        pl = [r for r in h["folds"] if r["fold"].startswith("leave out")]
        h["share_place_folds"] = float(np.mean([r["formula_bpr"] < r["one_correction_bpr"] for r in pl]))
        cb = [r for r in h["folds"] if r["fold"] == "calibration -> blind"]
        h["cal_to_blind"] = (cb[0]["formula_bpr"] < cb[0]["one_correction_bpr"]) if cb else None
        h["counts"] = h["share_place_folds"] >= 0.8 and h["cal_to_blind"] is not False
    res["neighbour_same_storm"] = neighbours(df, ["place", "storm", "type"])
    res["neighbour_same_station"] = neighbours(df, ["place", "station", "type"])
    res["by_place"] = {p: {"n": int((df.place == p).sum()), "mean_y": float(df.y[df.place == p].mean()),
                           "sd_y": float(df.y[df.place == p].std()),
                           "share_within_25": float(np.mean(np.abs(df.pred[df.place == p] / df.obs[df.place == p] - 1) <= 0.25)),
                           "share_obs_high": float(np.mean(df.y[df.place == p] > 0))} for p in sorted(df.place.unique())}
    print(f"\n== {name}: n {len(df)}  one {res['one_correction']['bits_per_reading']:.3f}  sd {res['one_correction']['sd']:.3f}")
    for t in res["top"][:6]:
        print("  ", t["terms"], round(t["bits_per_reading"], 3), {k: round(v, 3) for k, v in t["fit_coef"].items()},
              {k: round(v, 1) for k, v in t["fit_alpha"].items()}, "sd", round(t["sd"], 3))
    for h in res["holdout"]:
        print("   HO", h["terms"], "places", round(h["share_place_folds"], 2), "cal->blind", h["cal_to_blind"], "counts", h["counts"])
    print("   nb storm", res["neighbour_same_storm"], "\n   nb station", res["neighbour_same_station"])
    return res


def strip(res): return {k: v for k, v in res.items() if not k.startswith("_")}


def main(ws, a3dir, hdr):
    meta = stations_meta(hdr); ws = Path(ws); l1 = ws / "docs" / "validation" / "layer1"
    after = pd.concat([load(sorted((l1 / "f2" / "cal").glob("*.json")), "cal", meta),
                       load(sorted((l1 / "f2" / "blind").glob("*.json")), "blind", meta)])
    before = pd.concat([load(sorted((l1 / "run4" / "checks").glob("*.json")), "cal", meta),
                        load(sorted(Path(a3dir).glob("*.json")), "cal", meta)])
    out = {"sources": {"windsim_pr10": "d36598b", "f2_commit": "63dd5c5", "run4": "1297da7", "a3": "622f9e6"}}

    # Q1 thunderstorms
    th = features(after[after.type == "thunderstorm"])
    q1 = run_case(th, "Q1 thunderstorm, after F2")
    places = th.place.values
    q1["spread_per_place"] = spread_test(th.y.values, q1["_mu0"], places, G)
    q1["spread_desert"] = spread_test(th.y.values, q1["_mu0"], np.where(places == "phoenix", "phoenix", "other"), G)
    q1["data"] = th[["place", "group", "station", "storm", "y", "dirn", "hl", "month", "lnreg", "dz_elev", "storm_loo", "stn_loo"]].to_dict("list")
    print("  spread per place", round(q1["spread_per_place"]["saving_bpr"], 4), " desert", round(q1["spread_desert"]["saving_bpr"], 4))
    out["Q1_thunderstorm"] = strip(q1)

    # Q2 ordinary storms
    nt = features(after[after.type == "non-thunderstorm"])
    q2 = run_case(nt, "Q2 ordinary storms, after F2")
    ca_folds = []
    for t in q2["top"][:8]:
        rows = []
        for p in sorted(CA & set(nt.place)):
            tr, te = np.where(nt.place != p)[0], np.where(nt.place == p)[0]
            lib = library(nt); lib.pop("storm_loo")
            rows.append({"fold": f"leave out {p}", "formula_bpr": L.holdout(lib, t["terms"], nt.y.values, G, tr, te) / len(te),
                         "one_correction_bpr": L.holdout(lib, [], nt.y.values, G, tr, te) / len(te)})
        ca_folds.append({"terms": t["terms"], "folds": rows, "all_ca_hold": all(r["formula_bpr"] < r["one_correction_bpr"] for r in rows)})
    q2["ca_holdout"] = ca_folds
    for t in q2["holdout"]: pass
    q2["named_ca_terms"] = {t: L.score(library(nt), [t], nt.y.values, G)["bits_per_reading"] for t in
                            ("ca", "offshore_ca", "dz_elev", "elev", "msin", "mcos", "lnz0", "near_km")}
    q2["spread_per_place"] = spread_test(nt.y.values, q2["_mu0"], nt.place.values, G)
    q2["spread_ca"] = spread_test(nt.y.values, q2["_mu0"], np.where(nt.place.isin(CA), "california", "other"), G)
    for c in ca_folds: print("   CA-HO", c["terms"], c["all_ca_hold"], [(round(r["formula_bpr"], 3), round(r["one_correction_bpr"], 3)) for r in c["folds"]])
    print("   one-term", {k: round(v, 3) for k, v in q2["named_ca_terms"].items()})
    print("  spread per place", round(q2["spread_per_place"]["saving_bpr"], 4), " CA", round(q2["spread_ca"]["saving_bpr"], 4))
    q2["data"] = nt[["place", "group", "station", "storm", "y", "dirn", "hl", "month", "dz_elev", "elev", "z0"]].to_dict("list")
    out["Q2_ordinary"] = strip(q2)

    # Q3 before vs after the hurricane split, 10 calibration places
    trop = after[after.type == "tropical"]; trop_keys = set(zip(trop.place, trop.storm))
    b = features(before[before.type == "non-thunderstorm"])
    b_tr = np.array([(p, s) in trop_keys for p, s in zip(b.place, b.storm)], float)
    a = features(after[(after.group == "cal") & (after.type == "non-thunderstorm")])
    q3b = run_case(b, "Q3 before F2 (tropical inside ordinary)", extra_terms={"tropical": b_tr - b_tr.mean()})
    q3a = run_case(a, "Q3 after F2 (ordinary only)")
    q3 = {"before": strip(q3b), "after": strip(q3a), "tropical_pairs_in_before": int(b_tr.sum()),
          "tropical_term_alone": L.score({"tropical": b_tr - b_tr.mean()}, ["tropical"], b.y.values, G)["bits_per_reading"]}
    keys_a = set(zip(a.place, a.storm, a.station)); shared = np.array([(p, s, t) in keys_a for p, s, t in zip(b.place, b.storm, b.station)])
    keys_b = set(zip(b.place, b.storm, b.station)); shared_a = np.array([(p, s, t) in keys_b for p, s, t in zip(a.place, a.storm, a.station)])
    q3["shared_pairs"] = int(shared.sum())
    q3["shared_one_correction_before_after"] = (L.score(library(b), [], b.y.values, G, rows=np.where(shared)[0])["bits_per_reading"],
                                                L.score(library(a), [], a.y.values, G, rows=np.where(shared_a)[0])["bits_per_reading"])
    tr_all = features(after[(after.group == "cal") & (after.type == "tropical")])
    q3["tropical_after"] = {"n": len(tr_all), "one_correction_bpr": L.score(library(tr_all), [], tr_all.y.values, G)["bits_per_reading"],
                            "sd": float(tr_all.y.std())}
    print("\nQ3 tropical pairs in before", q3["tropical_pairs_in_before"], "tropical alone", round(q3["tropical_term_alone"], 3),
          "shared", q3["shared_pairs"], q3["shared_one_correction_before_after"], "tropical after", q3["tropical_after"])
    out["Q3_hurricane_split"] = q3
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / ("layer1_firstrun.json" if "--first" in sys.argv else "layer1.json")).write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main(*sys.argv[1:4])
