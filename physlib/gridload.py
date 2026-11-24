"""PROTOCOL-8: the physics compressor on electric grid load (ISO New England, GEFCom2017 data). Written fresh for
SPARC / gridzip-lab (no lpaq1, cmix or other GPL code).

Data: the GEFCom2017 "defined data" track, ISO New England hourly demand with dry-bulb and dew-point temperature for
eight load zones and two aggregates, March 2003 to April 2017, as packaged by camroach87/gefcom2017data
(data/gefcom.rda). The underlying numbers are ISO-NE's public SMD hourly files.

python -m physlib.gridload fetch            download data/gefcom.rda into GRIDLOAD_DATA (default ~/.cache/gridload)
python -m physlib.gridload timing           engine time on pure noise (no demand read), to fix the subsample step
python -m physlib.gridload planted          part A: planted laws on real weather inputs, synthetic outcome
python -m physlib.gridload explore          part B: what the vanilla benchmark misses, 2011-2013 (exploration years)
python -m physlib.gridload sealed           part C: the head-to-head on 2014-2016, read once
Results go to results/gridload/.
"""
from __future__ import annotations

import json, math, os, sys, time
from pathlib import Path

import numpy as np

from . import engine as E, guards as Q
from .engine import Term

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "gridload"
DATA = Path(os.environ.get("GRIDLOAD_DATA", Path.home() / ".cache" / "gridload"))
URL = "https://raw.githubusercontent.com/camroach87/gefcom2017data/master/data/gefcom.rda"
ZONES = ("TOTAL", "ME", "NH", "VT", "CT", "RI", "SEMASS", "WCMASS", "NEMASSBOST")
TRAIN = (2004, 2010)        # vanilla fitted here for the exploration (part B)
EXPLORE = (2011, 2012)      # compressor searches here
CHECK = (2013, 2013)        # its formula must transfer here
FIT_ALL = (2004, 2013)      # every head-to-head model is fitted here
SEALED = (2014, 2016)       # read once (part C)
STEP = int(os.environ.get("GRIDLOAD_STEP", "1"))   # compressor reads every STEP-th hour (fixed in PROTOCOL-8)
LAT, LON = 42.3, -71.8      # centre of New England for the sun terms
G = E.Grid(0.002, -1.0, 1.0)


# ------------------------------------------------------------------------------------------------------- data
def fetch():
    import urllib.request
    DATA.mkdir(parents=True, exist_ok=True); p = DATA / "gefcom.rda"
    if not p.exists(): urllib.request.urlretrieve(URL, p)
    return p


_FRAME = None


def frame():
    global _FRAME
    if _FRAME is None:
        import pyreadr
        _FRAME = pyreadr.read_r(str(fetch()))["gefcom"]
    return _FRAME


def zone(name):
    """One zone as a dict of hourly arrays on a gap-free clock (missing hours NaN), with the derived inputs."""
    import pandas as pd
    df = frame(); z = df[df["zone"].str.upper() == name].copy()
    z = z.groupby("ts", as_index=False).agg(demand=("demand", "mean"), drybulb=("drybulb", "mean"), dewpnt=("dewpnt", "mean"),
                                            holiday=("holiday", "max"))
    z = z.set_index("ts").reindex(pd.date_range(z["ts"].min(), z["ts"].max(), freq="h"))
    return derive(z)


def _rh(T_F, Td_F):
    T = (T_F - 32) / 1.8; Td = (Td_F - 32) / 1.8
    es = lambda t: np.exp(17.625 * t / (243.04 + t))
    return np.clip(100 * es(Td) / es(T), 1, 100)


def _sun_sin(ts):
    """sin(solar elevation) at LAT, LON; timestamps are local clock time (ISO-NE prevailing time, hour beginning)."""
    import pandas as pd
    loc = pd.DatetimeIndex(ts).tz_localize("America/New_York", ambiguous="NaT", nonexistent="shift_forward")
    utc = loc.tz_convert("UTC"); ok = ~utc.isna()
    doy = np.where(ok, utc.dayofyear, 1).astype(float); hr = np.where(ok, utc.hour + utc.minute / 60 + 0.5, 0.0)
    dec = np.radians(23.44) * np.sin(2 * np.pi * (284 + doy) / 365)
    ha = np.radians(15 * (hr + LON / 15 - 12))
    s = np.sin(np.radians(LAT)) * np.sin(dec) + np.cos(np.radians(LAT)) * np.cos(dec) * np.cos(ha)
    return np.where(ok, s, np.nan)


def derive(z):
    import pandas as pd
    ts = z.index; T = z["drybulb"].to_numpy(float); Td = z["dewpnt"].to_numpy(float)
    Ts = pd.Series(T, index=ts)
    d = {"ts": ts.to_numpy(), "y": z["demand"].to_numpy(float), "T": T, "Td": Td, "RH": _rh(T, Td),
         "hour": ts.hour.to_numpy(float), "dow": ts.dayofweek.to_numpy(float), "month": ts.month.to_numpy(float),
         "doy": ts.dayofyear.to_numpy(float), "year": ts.year.to_numpy(float),
         "trend": (ts - pd.Timestamp("2004-01-01")).total_seconds().to_numpy() / 3600 / 8766.0}
    hol = z["holiday"].astype(float).fillna(0).to_numpy()
    day = ts.normalize(); dh = pd.Series(hol, index=ts).groupby(day).max()
    d["holiday"] = hol
    d["hol_adj"] = np.clip(day.map((dh.shift(1, freq="D").reindex(dh.index).fillna(0) + dh.shift(-1, freq="D").reindex(dh.index).fillna(0))).to_numpy(float), 0, 1) * (1 - hol)
    d["T_lag3"] = Ts.shift(3).to_numpy(); d["T_lag1"] = Ts.shift(1).to_numpy(); d["T_lag2"] = Ts.shift(2).to_numpy()
    d["Tavg24"] = Ts.shift(1).rolling(24, min_periods=20).mean().to_numpy()
    d["Tavg48"] = Ts.shift(1).rolling(48, min_periods=40).mean().to_numpy()
    dmax = Ts.groupby(day).max(); d["Tmax_prev"] = day.map(dmax.shift(1, freq="D")).to_numpy(float)
    hot = (dmax >= 85).astype(int); run = hot.groupby((hot != hot.shift()).cumsum()).cumsum() * hot
    d["heatrun"] = day.map(run.shift(1, freq="D")).to_numpy(float)
    s = _sun_sin(ts); d["sun"] = np.maximum(0, s)
    d["dark_awake"] = ((s < 0) & (d["hour"] >= 6) & (d["hour"] <= 22)).astype(float)
    d["weekend"] = (d["dow"] >= 5).astype(float)
    return d


def rows(d, years, step=1, need=("y",)):
    """Indices in the year range with every needed input finite, every step-th hour of the clock."""
    yr = d["year"]; m = (yr >= years[0]) & (yr <= years[1])
    for k in need + ("T", "Td", "RH", "T_lag1", "T_lag2", "T_lag3", "Tavg24", "Tavg48", "Tmax_prev", "heatrun", "sun"):
        m &= np.isfinite(d[k])
    idx = np.flatnonzero(m)
    return idx[(idx % step) == 0] if step > 1 else idx


# ---------------------------------------------------------------------------------------------------- models
def _oh(v, n):
    v = np.asarray(v, int); out = np.zeros((len(v), n)); out[np.arange(len(v)), v] = 1; return out


def vanilla_X(d, idx):
    """Hong's GEFCom2012 'vanilla' benchmark: Trend + M + W×H + T + T² + T³ + T·M + T²·M + T³·M + T·H + T²·H + T³·H
    (temperature scaled to T/100 for conditioning; M and H one-hot, so T, T², T³ main effects sit in the T·M blocks)."""
    t = d["T"][idx] / 100; M_ = _oh(d["month"][idx] - 1, 12); H = _oh(d["hour"][idx], 24); W = _oh(d["dow"][idx], 7)
    WH = (W[:, :, None] * H[:, None, :]).reshape(len(idx), -1)
    cols = [d["trend"][idx][:, None], M_, WH]
    for p in (1, 2, 3): cols += [M_ * t[:, None] ** p, H * t[:, None] ** p]
    return np.column_stack(cols)


def expert_X(d, idx):
    """Pre-registered expert model (PROTOCOL-8 M2), written from memory of the load-forecasting literature (recency
    effect, relative humidity, holidays; papers not opened): vanilla + the same polynomial-by-hour blocks for Tavg24 and
    for the 1-3 h lags + RH, RH·T, RH·T² by month + holiday and day-next-to-holiday by hour."""
    X = [vanilla_X(d, idx)]; H = _oh(d["hour"][idx], 24); M_ = _oh(d["month"][idx] - 1, 12)
    for k in ("Tavg24", "T_lag1", "T_lag2", "T_lag3"):
        t = d[k][idx] / 100
        for p in (1, 2, 3): X.append(H * t[:, None] ** p)
    rh = d["RH"][idx] / 100; t = d["T"][idx] / 100
    for c in (rh, rh * t, rh * t * t): X.append(M_ * c[:, None])
    X += [H * d["holiday"][idx][:, None], H * d["hol_adj"][idx][:, None]]
    return np.column_stack(X)


def ols(X, y):
    sc = np.sqrt((X ** 2).mean(0)) + 1e-12
    th, *_ = np.linalg.lstsq(X / sc, y, rcond=None); return th / sc


GBM_FEATS = ("hour", "dow", "month", "doy", "trend", "holiday", "hol_adj", "T", "Td", "RH", "T_lag1", "T_lag2", "T_lag3",
             "Tavg24", "Tavg48", "Tmax_prev", "heatrun", "sun")


def gbm(d, tr, te, seed=0, extra=None):
    from sklearn.ensemble import HistGradientBoostingRegressor
    F = lambda i: np.column_stack([d[k][i] for k in GBM_FEATS] + ([extra(i)] if extra else []))
    m = HistGradientBoostingRegressor(max_iter=800, learning_rate=0.05, max_leaf_nodes=63, l2_regularization=1.0,
                                      early_stopping=False, random_state=seed)
    m.fit(F(tr), np.log(d["y"][tr])); return np.exp(m.predict(F(te)))


def mape(y, p): return float(np.mean(np.abs(p - y) / y) * 100)


# ----------------------------------------------------------------------------------------------- the library
def _hinge(x, q): return np.maximum(0.0, x - q)


LOAD = [
    Term("recency_cool", "load", "max(0, T̄24 − 65 °F)", "buildings stay hot after a hot day: demand follows yesterday's heat, not only this hour's",
         ("Tavg24",), lambda c: _hinge(c["Tavg24"], 65.0), group="recency"),
    Term("recency_heat", "load", "max(0, 50 °F − T̄24)", "cold soaks in: heating demand follows the last day's cold",
         ("Tavg24",), lambda c: _hinge(50.0 - c["Tavg24"], 0.0), group="recency_h"),
    Term("lag_gap", "load", "T(t−3 h) − T(t)", "thermal lag: demand answers temperature a few hours late",
         ("T_lag3", "T"), lambda c: c["T_lag3"] - c["T"], group="lag"),
    Term("prev_day_heat", "load", "max(0, Tmax(yesterday) − 80 °F)", "a hot yesterday raises today's base (pre-cooled mass, behaviour)",
         ("Tmax_prev",), lambda c: _hinge(c["Tmax_prev"], 80.0), group="recency"),
    Term("heat_wave", "load", "consecutive days with Tmax ≥ 85 °F", "demand climbs with each day of a heat wave",
         ("heatrun",), lambda c: c["heatrun"], group="heatwave"),
    Term("humid_cool", "load", "max(0, T − 70 °F) · RH/100", "humid heat: air conditioners also remove moisture",
         ("T", "RH"), lambda c: _hinge(c["T"], 70.0) * c["RH"] / 100, group="humidity"),
    Term("dewpoint", "load", "max(0, T_d − 60 °F)", "muggy air raises cooling load beyond dry-bulb temperature",
         ("Td",), lambda c: _hinge(c["Td"], 60.0), group="humidity"),
    Term("holiday", "load", "[holiday]", "the model treats holidays as ordinary days", ("holiday",), lambda c: c["holiday"], group="holiday"),
    Term("holiday_adj", "load", "[day next to a holiday]", "bridge days: people take the day off around holidays",
         ("hol_adj",), lambda c: c["hol_adj"], group="holiday_adj"),
    Term("lighting", "load", "[sun down, 06-22 h]", "lighting load follows sunset, which clock-hour terms miss across the year and DST",
         ("dark_awake",), lambda c: c["dark_awake"], group="light"),
    Term("btm_solar", "load", "sin⁺(elev) · max(0, year − 2010)", "rooftop solar growing behind the meter hides midday demand",
         ("sun", "year"), lambda c: c["sun"] * _hinge(c["year"], 2010.0), group="solar"),
    Term("weekend_temp", "load", "[weekend] · (T − 60 °F)", "weather sensitivity differs on weekends (homes vs offices)",
         ("weekend", "T"), lambda c: c["weekend"] * (c["T"] - 60.0), group="behaviour"),
    Term("trend_curve", "load", "(years since 2004)²", "growth is not linear (efficiency, economy)", ("trend",), lambda c: c["trend"] ** 2, group="trend"),
    Term("daily_shape", "load", "sin, cos(2π·hour/24)", "shape only: the daily curve is off", ("hour",),
         lambda c: np.column_stack([np.sin(2 * np.pi * c["hour"] / 24), np.cos(2 * np.pi * c["hour"] / 24)]), physics=False, group="shape"),
]
GRAMMAR_INPUTS = ("T", "Td", "Tavg24", "Tmax_prev", "hour", "doy", "trend", "sun")


def ctx_of(d, idx):
    c = {k: d[k][idx].astype(float) for k in ("T", "Td", "RH", "Tavg24", "T_lag3", "Tmax_prev", "heatrun", "holiday", "hol_adj",
                                                "dark_awake", "sun", "year", "weekend", "trend", "hour", "doy")}
    return c


def prev_of(idx, step):
    p = np.full(len(idx), -1); pos = {v: i for i, v in enumerate(idx)}
    for i, v in enumerate(idx):
        j = pos.get(v - step)
        if j is not None: p[i] = j
    return p


def analyse(c, r, idx, step, null_N=0, bootstrap_B=0, folds=(), seed=0, grammar=True, name=""):
    prev = prev_of(idx, step)
    excl = tuple(k for k in c if k not in GRAMMAR_INPUTS)
    weeks = (idx // (24 * 7)).astype(int); weeks = [np.flatnonzero(weeks == w) for w in np.unique(weeks)]
    return Q.analyse(LOAD, c, r, G, grammar=grammar, robust=False, prev=prev, folds=folds, bootstrap_B=bootstrap_B, null_N=null_N,
                     neighbour=False, name=name, exclude=excl, blocks=weeks, seed=seed, boot_by_group=False)


def summary(res):
    gd = res["guards"]
    out = {"winner": res["winner"]["terms"], "coef": res["winner"]["coef"], "alpha": res["winner"]["alpha"],
           "bits_saved_per_reading": res["bits_saved_per_reading"], "sd_before": res["sd_before"], "sd_after": res["sd_after"],
           "transfers": res["transfers"], "holdout": res["holdout"], "verdict": Q.verdict(res), "n": res["n"],
           "inclusion": gd.get("inclusion"), "lookalikes": gd.get("lookalikes")}
    if "null_bits_saved" in gd: out["null_bits_saved"] = gd["null_bits_saved"]
    if "bootstrap" in gd: out["bootstrap"] = gd["bootstrap"]
    return out


def _dump(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(obj, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))


# ------------------------------------------------------------------------------------------------- part 0: timing
def timing():
    d = zone("TOTAL"); idx = rows(d, EXPLORE, STEP, need=()); c = ctx_of(d, idx)
    rng = np.random.default_rng(0); r = rng.normal(0, 0.03, len(idx))
    t = time.time(); res = analyse(c, r, idx, STEP, name="timing"); dt = time.time() - t
    print(json.dumps({"n": len(idx), "step": STEP, "seconds": round(dt, 1), "winner": res["winner"]["terms"],
                      "grammar_size": res["grammar_size"]}))


# ------------------------------------------------------------------------------------------------- part A: planted
PLANTS = {
    "P1_recency": ["recency_cool"],
    "P2_humidity": ["humid_cool"],
    "P3_holiday": ["holiday"],
    "P4_solar": ["btm_solar"],
    "P5_combo": ["recency_cool", "humid_cool", "holiday"],
    "P6_grammar": ["g:max(0, Tmax_prev−q50) · hour"],
}


def _noise(rng, idx, step, sd, rho):
    prev = prev_of(idx, step); e = rng.normal(0, sd * math.sqrt(1 - rho * rho), len(idx)); out = np.empty(len(idx))
    for i in range(len(idx)):
        out[i] = rho * out[prev[i]] + e[i] if prev[i] >= 0 else rng.normal(0, sd)
    return out


def planted(seeds=10, null_seeds=20, years=EXPLORE, zone_name="TOTAL"):
    d = zone(zone_name); idx = rows(d, years, STEP, need=()); c = ctx_of(d, idx)
    by = {t.name: t for t in LOAD}; SD, RHO = 0.03, 0.9
    res = {}
    for pname, terms in list(PLANTS.items()) + [("N_null", [])]:
        rows_ = []
        for s in range(null_seeds if not terms else seeds):
            rng = np.random.default_rng(1000 * len(res) + s)
            sig = np.zeros(len(idx)); planted_groups = []
            for tn in terms:
                if tn.startswith("g:"):
                    x = c["Tmax_prev"]; v = _hinge(x, float(np.quantile(x, 0.5))) * c["hour"]; grp = "grammar"
                else:
                    v = E._as2d(by[tn].build(c))[:, 0]; grp = by[tn].group
                v = v - v.mean(); v = v / (v.std() + 1e-300)
                sig += v * (0.5 * SD / math.sqrt(len(terms))) * (rng.choice([-1, 1]) if tn != "holiday" else -1)
                planted_groups.append(grp)
            r = sig + _noise(rng, idx, STEP, SD, RHO)
            t0 = time.time(); a = analyse(c, r, idx, STEP, seed=s, name=f"{pname} seed {s}")
            win = a["winner"]["terms"]; gw = [next((t.group for t in LOAD if t.name == w), "grammar") for w in win]
            named = all(g_ in gw for g_ in planted_groups) if terms else (len(win) == 0)
            extra = [w for w, g_ in zip(win, gw) if g_ not in planted_groups]
            mu = a["_mu"] - a["_mu0"]; r2 = float(np.corrcoef(mu, sig)[0, 1] ** 2) if terms and np.std(mu) > 0 else None
            rows_.append({"seed": s, "winner": win, "named": named, "extra": extra, "r2": r2, "sec": round(time.time() - t0, 1)})
            print(pname, s, win, named, f"{rows_[-1]['sec']}s", flush=True)
        ok = sum(x["named"] for x in rows_); ex = sum(1 for x in rows_ if x["extra"])
        res[pname] = {"planted": terms, "named": ok, "of": len(rows_), "with_extra_terms": ex, "rows": rows_}
    _dump("A_planted", res); return res


# ------------------------------------------------------------------------------------------------- part B: explore
def explore(zones=ZONES, null_N=9, bootstrap_B=10):
    out = {}
    for zn in zones:
        d = zone(zn)
        tr = rows(d, TRAIN); ex = rows(d, (EXPLORE[0], CHECK[1]), STEP)
        th = ols(vanilla_X(d, tr), d["y"][tr]); pred = vanilla_X(d, ex) @ th
        r = np.log(d["y"][ex] / pred); c = ctx_of(d, ex)
        fo = [("check year 2013", np.flatnonzero(d["year"][ex] <= EXPLORE[1]), np.flatnonzero(d["year"][ex] == CHECK[0]))]
        sel = np.flatnonzero(d["year"][ex] <= EXPLORE[1])
        cc = {k: v[sel] for k, v in c.items()}
        t0 = time.time()
        a = analyse(cc, r[sel], ex[sel], STEP, null_N=null_N if zn == "TOTAL" else 0, bootstrap_B=bootstrap_B if zn == "TOTAL" else 0,
                    name=f"ISO-NE {zn} vanilla residuals 2011-12")
        s = summary(a)
        # transfer to 2013: fit the formula on 2011-12, code 2013 with it (AR(1)), against one correction
        built = E.assemble(list(LOAD) + (__import__("physlib.grammar", fromlist=["x"]).candidates(c, len(ex), exclude=tuple(k for k in c if k not in GRAMMAR_INPUTS))), c, len(ex))
        names = [b.name for b in built]
        terms_ok = [t for t in s["winner"] if t in names]
        s["transfer_2013"] = E.holdout(built, terms_ok, r, G, fo[0][1], fo[0][2], E.library_info(built), prev=prev_of(ex, STEP)) if terms_ok else None
        s["vanilla_mape_2011_13"] = mape(d["y"][ex], pred); s["sec"] = round(time.time() - t0, 1)
        out[zn] = s; print(zn, s["winner"], round(s["bits_saved_per_reading"], 3), s["transfer_2013"], f"{s['sec']}s", flush=True)
        _dump("B_explore", out)
    return out


# ------------------------------------------------------------------------------------------------- part C: sealed
def found_X(d, idx, winner, alphas=None):
    """The compressor's found terms as extra regression columns (library terms rebuilt; grammar terms rebuilt from their
    label with the quantile constant fixed at its exploration value, carried in the label)."""
    by = {t.name: t for t in LOAD}; c = ctx_of(d, idx); cols = []
    for w in winner:
        if w in by: cols.append(E._as2d(by[w].build(c)))
        else: cols.append(E._as2d(grammar_col(w, c)))
    return np.column_stack(cols) if cols else np.zeros((len(idx), 0))


def grammar_col(label, c):
    """Rebuild a grammar term from its label (labels carry the constant to 3 significant figures)."""
    import re
    lab = label[2:]
    ops = [(r"^max\(0, (\w+)−([-\d.e]+)\)$", lambda x, q: np.maximum(0, x - q)), (r"^\[(\w+)>([-\d.e]+)\]$", lambda x, q: (x > q) * 1.0),
           (r"^e\^\(−(\w+)/([-\d.e]+)\)$", lambda x, q: np.exp(-x / q)), (r"^(\w+)²$", lambda x, q: x ** 2), (r"^(\w+)$", lambda x, q: x),
           (r"^ln (\w+)$", lambda x, q: np.log(x)), (r"^√(\w+)$", lambda x, q: np.sqrt(x)), (r"^1/(\w+)$", lambda x, q: 1 / x)]
    def one(s):
        for pat, f in ops:
            m = re.match(pat, s.strip())
            if m: return f(c[m.group(1)], float(m.group(2)) if m.lastindex and m.lastindex > 1 else None)
        raise ValueError(f"can't rebuild {s!r}")
    return np.prod([one(p) for p in lab.split(" · ")], axis=0)


def sealed(zones=ZONES):
    B = json.loads((OUT / "B_explore.json").read_text()); out = {}
    for zn in zones:
        d = zone(zn); tr = rows(d, FIT_ALL); te = rows(d, SEALED); y = d["y"]
        win = B[zn]["winner"]
        p0 = vanilla_X(d, te) @ ols(vanilla_X(d, tr), y[tr])
        X1 = lambda i: np.column_stack([vanilla_X(d, i), found_X(d, i, win)])
        p1 = X1(te) @ ols(X1(tr), y[tr])
        p2 = expert_X(d, te) @ ols(expert_X(d, tr), y[tr])
        p3 = gbm(d, tr, te)
        o = {"found_terms": win, "n_test_hours": len(te)}
        for k, p in (("M0_vanilla", p0), ("M1_vanilla_plus_compressor", p1), ("M2_expert", p2), ("M3_gbm", p3)):
            o[k] = {"mape": mape(y[te], p), "by_year": {int(Y): mape(y[te][d["year"][te] == Y], p[d["year"][te] == Y]) for Y in range(SEALED[0], SEALED[1] + 1)}}
            top = np.argsort(-y[te])[: max(1, len(te) // 100)]   # the top 1 % of hours (peaks)
            o[k]["peak_mape"] = mape(y[te][top], p[top])
        g0 = o["M0_vanilla"]["mape"]; o["gain_M1"] = g0 - o["M1_vanilla_plus_compressor"]["mape"]; o["gain_M2"] = g0 - o["M2_expert"]["mape"]
        o["share_of_expert_gain"] = o["gain_M1"] / o["gain_M2"] if o["gain_M2"] > 0 else None
        out[zn] = o; print(zn, {k: round(v["mape"], 3) for k, v in o.items() if k.startswith("M")}, flush=True)
        _dump("C_sealed", out)
    return out


# --------------------------------------------------------------------------------- A1 parts D and E: audit strong models
FRESH = (2017, 2017)


def audit(zones=ZONES, null_N=9):
    out = {}
    for zn in zones:
        d = zone(zn); tr = rows(d, TRAIN); ex = rows(d, (EXPLORE[0], CHECK[1])); y = d["y"]
        preds = {"M2_expert": expert_X(d, ex) @ ols(expert_X(d, tr), y[tr]), "M3_gbm": gbm(d, tr, ex)}
        o = {}
        for m, pred in preds.items():
            r = np.log(y[ex] / pred); c = ctx_of(d, ex)
            sel = np.flatnonzero(d["year"][ex] <= EXPLORE[1]); cc = {k: v[sel] for k, v in c.items()}
            t0 = time.time()
            a = analyse(cc, r[sel], ex[sel], STEP, null_N=null_N if zn == "TOTAL" else 0, name=f"ISO-NE {zn} {m} residuals 2011-12")
            s = summary(a)
            from . import grammar as GR
            built = E.assemble(list(LOAD) + GR.candidates(c, len(ex), exclude=tuple(k for k in c if k not in GRAMMAR_INPUTS)), c, len(ex))
            names = [b.name for b in built]; ok = [t for t in s["winner"] if t in names]
            fo = (np.flatnonzero(d["year"][ex] <= EXPLORE[1]), np.flatnonzero(d["year"][ex] == CHECK[0]))
            s["transfer_2013"] = E.holdout(built, ok, r, G, fo[0], fo[1], E.library_info(built), prev=prev_of(ex, STEP)) if ok else None
            s["mape_2011_13"] = mape(y[ex], pred); s["sec"] = round(time.time() - t0, 1); o[m] = s
            print(zn, m, s["winner"], round(s["bits_saved_per_reading"], 3), s["transfer_2013"], s.get("null_bits_saved", {}).get("max"), f"{s['sec']}s", flush=True)
        out[zn] = o; _dump("D_audit", out)
    return out


def fresh(zones=ZONES):
    D = json.loads((OUT / "D_audit.json").read_text()); out = {}
    for zn in zones:
        d = zone(zn); tr = rows(d, FIT_ALL); y = d["y"]; o = {}
        ok = lambda w: [t for t in w if np.all(np.isfinite(found_X(d, tr, [t])))]   # A1 note: a term that is not finite over
        w2 = ok(D[zn]["M2_expert"]["winner"]); w3 = ok(D[zn]["M3_gbm"]["winner"])  # the fit years (1/trend at 2004) is dropped
        for lab, yrs in (("2014_2016_post_hoc", SEALED), ("jan_apr_2017_fresh", FRESH)):
            te = rows(d, yrs)
            p2 = expert_X(d, te) @ ols(expert_X(d, tr), y[tr])
            X2 = lambda i: np.column_stack([expert_X(d, i), found_X(d, i, w2)])
            p2p = X2(te) @ ols(X2(tr), y[tr]) if w2 else p2
            p3 = gbm(d, tr, te)
            p3p = gbm(d, tr, te, extra=(lambda i: found_X(d, i, w3))) if w3 else p3
            o[lab] = {"n_hours": len(te), "M2": mape(y[te], p2), "M2_plus": mape(y[te], p2p), "M3": mape(y[te], p3), "M3_plus": mape(y[te], p3p)}
        o["terms_M2"] = w2; o["terms_M3"] = w3; out[zn] = o
        print(zn, {k: {kk: round(vv, 3) for kk, vv in v.items()} for k, v in o.items() if isinstance(v, dict)}, flush=True)
        _dump("E_fresh", out)
    return out


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    cmd = sys.argv[1] if len(sys.argv) > 1 else "timing"
    {"fetch": lambda: print(fetch()), "timing": timing, "planted": planted, "explore": explore, "sealed": sealed, "audit": audit, "fresh": fresh}[cmd]()
