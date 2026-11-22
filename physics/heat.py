"""PROTOCOL-2 case H: HeatSim 0.6.1 open-site air temperature (station correction off) against 3 Philadelphia ISD stations.

python physics/heat.py <dir with phl_noanchor_*.nc> <dir with ISD 2025 csvs>
"""
import json, math, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from pyproj import Transformer
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import runner
from codelen import lawsearch as L

OUT = Path(__file__).resolve().parents[1] / "results" / "physics"
G = L.Grid(0.1, -20.0, 20.0)
STATIONS = {"72030464752": "Wings Field", "72408013739": "PHL airport", "72408594732": "NE PHL airport"}
DAYS = ("2025-07-08", "2025-07-25", "2025-07-28")


def isd(path):
    d = pd.read_csv(path, low_memory=False, usecols=["DATE", "LATITUDE", "LONGITUDE", "WND", "TMP", "DEW", "GF1", "GA1"])
    t = pd.to_datetime(d["DATE"]).dt.floor("h")
    def num(col, i, scale, bad):
        v = d[col].astype(str).str.split(",").str[i]; v = pd.to_numeric(v, errors="coerce"); return v.where(v != bad) / scale
    T = num("TMP", 0, 10, 9999); Td = num("DEW", 0, 10, 9999); W = num("WND", 3, 10, 9999)
    gf = num("GF1", 0, 1, 99) if "GF1" in d else None; ga = num("GA1", 0, 1, 99)
    cl = (gf.where(gf <= 8) if gf is not None else np.nan); cl = cl.fillna(ga.where(ga <= 8)) / 8
    f = pd.DataFrame({"t": t, "T": T, "Td": Td, "W": W, "cloud": cl}).groupby("t").mean()
    return f, float(d["LATITUDE"].iloc[0]), float(d["LONGITUDE"].iloc[0])


def sun_sin(times, lat, lon):
    """sin(solar elevation) from the NOAA approximation; times UTC."""
    ts = pd.DatetimeIndex(times); doy = ts.dayofyear.values; hr = ts.hour.values + ts.minute.values / 60
    g = 2 * math.pi / 365 * (doy - 1 + (hr - 12) / 24)
    decl = 0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g) + 0.000907 * np.sin(2 * g)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    tst = hr * 60 + eqt + 4 * lon; ha = np.radians(tst / 4 - 180); la = np.radians(lat)
    return np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)


def main(ncdir, isddir):
    to_utm = Transformer.from_crs("EPSG:4326", "EPSG:32618", always_xy=True).transform
    rows = []
    for day in DAYS:
        ds = xr.open_dataset(Path(ncdir) / f"phl_noanchor_{day}.nc"); assert int(ds.attrs["anchored"]) == 0
        ver = ds.attrs["heatsim_version"]
        for sid, nm in STATIONS.items():
            f, lat, lon = isd(Path(isddir) / f"{sid}.csv")
            sx, sy = to_utm(lon, lat); ix = int(np.argmin(np.abs(ds.x.values - sx))); iy = int(np.argmin(np.abs(ds.y.values - sy)))
            tm = ds.open_site_temperature.values[:, iy, ix]; times = pd.DatetimeIndex(ds.time.values)
            o = f.reindex(times)
            fr = {k: float(ds[k].values[iy, ix]) for k in ("frac_low", "frac_water", "frac_roof", "frac_canyon", "frac_tree")}
            for k, t in enumerate(times):
                rows.append(dict(day=day, sid=sid, t=t, Tm=float(tm[k]), To=o["T"].iloc[k], Td=o["Td"].iloc[k], W=o["W"].iloc[k],
                                 cloud=o["cloud"].iloc[k], lat=lat, lon=lon, **fr))
    df = pd.DataFrame(rows)
    df["cloud"] = df["cloud"].fillna(0.5)
    df = df.dropna(subset=["To", "Tm"]).reset_index(drop=True)
    df["W"] = df["W"].fillna(df["W"].median()); df["Td"] = df["Td"].fillna(df["Td"].median())
    r = (df["To"] - df["Tm"]).values
    sun = lambda lag: np.maximum(sun_sin(df["t"] - pd.to_timedelta(lag, unit="h"), df["lat"].values, df["lon"].values), 0)
    hl = ((df["t"].dt.hour + df["t"].dt.minute / 60 - 4) % 24).values          # local daylight time, for the harmonics
    W = np.maximum(df["W"].values, 0.5); c = lambda a: a - np.mean(a)
    lib = {"sun": sun(0), "sun_lag": L.Family(lambda a: sun(a), 0.0, 6.0, n=25), "night": (sun(0) <= 0).astype(float),
           "hsin": np.sin(2 * np.pi * hl / 24), "hcos": np.cos(2 * np.pi * hl / 24), "wind": c(W), "1/wind": c(1 / W),
           "sun/wind": sun(0) / W, "dew": c(df["Td"].values), "cloud": c(df["cloud"].values), "Tm": c(df["Tm"].values),
           "stn_PHL": (df["sid"] == "72408013739").astype(float).values, "stn_NE": (df["sid"] == "72408594732").astype(float).values}
    ex = [{"sun", "sun_lag"}]
    folds = [(f"leave out {d}", np.where(df["day"] != d)[0], np.where(df["day"] == d)[0]) for d in DAYS] + \
            [(f"leave out {STATIONS[s]}", np.where(df["sid"] != s)[0], np.where(df["sid"] == s)[0]) for s in STATIONS]
    res = runner.case(lib, r, G, folds, exclusive=ex)
    # neighbour test: same hour at the other stations; shuffles over time within each station
    key = df["day"] + df["t"].astype(str); order = np.random.default_rng(0).permutation(len(r)); ctx = np.full((len(r), 2), -1)
    for i, j in enumerate(order):
        same = [q for q in np.where(key.values == key.values[j])[0] if q != j][:2]; ctx[i, :len(same)] = same
    nb, null = runner.neighbour_test(r, order, ctx, G, [np.where(df["sid"] == s)[0] for s in STATIONS])
    res["neighbour"] = nb
    res["version"] = ver; res["stations"] = STATIONS
    res["data"] = {"r": r.tolist(), "hour_local": hl.tolist(), "sid": df["sid"].tolist(), "day": df["day"].tolist(),
                   "Tm": df["Tm"].tolist(), "To": df["To"].tolist(), "W": W.tolist(), "cloud": df["cloud"].tolist(),
                   "mu_best": res["_mu"].tolist(), "cell_fractions": df.groupby("sid")[["frac_low", "frac_water", "frac_roof", "frac_canyon", "frac_tree"]].first().to_dict()}
    out = {k: v for k, v in res.items() if not k.startswith("_")}
    (OUT / "heat.json").write_text(json.dumps(out, indent=1, default=float))
    print("n", len(r), "one", round(res["one_correction"]["bits_per_reading"], 3), res["one_correction"]["fit_coef"], "sd", res["one_correction"]["sd"])
    for t in res["top"][:8]: print("  ", t["terms"], round(t["bits_per_reading"], 3), {k: round(v, 3) for k, v in t["fit_coef"].items()}, {k: round(v, 2) for k, v in t["fit_alpha"].items()}, "sd", t["sd"])
    for h in res["holdout"]: print("  HO", h["terms"], h["transfers"], [(round(f["formula_bpr"], 2), round(f["one_correction_bpr"], 2)) for f in h["folds"]])
    print("  neighbour", nb, "coder", res["coder_check"])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
