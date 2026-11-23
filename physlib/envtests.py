"""PROTOCOL-4 part C (A3): hard environments as tell-tale tests. python -m physlib.envtests

For each environment a model is likely to get wrong, plant the error it would make, then ask the compressor whether
it flags it, names the right term, stays quiet on a control, and what a missing term looks like (unnamed runs).
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from . import engine as E, terms as T, synth as S
from .engine import Term

OUT = Path(__file__).resolve().parents[1] / "results" / "physlib"
NSTN = 6


# ------------------------------------------------------------------------------------------------ inputs
def heat_ctx(rng, days=5):
    """6 stations, each with its own surroundings (crop and irrigation share, distance to the coast, height
    relative to the valley floor, canyon H/W, runway beside it), hourly over 5 days."""
    n1 = days * 24; hrs = np.arange(n1, dtype=float)
    hour = np.tile(hrs % 24, NSTN); day = np.tile(hrs // 24, NSTN); stn = np.repeat(np.arange(NSTN), n1)
    site = lambda v: np.asarray(v, float)[stn]
    base = S._sun(hour)
    Tm = 25 + rng.normal(0, 1.5, days)[day.astype(int)] + 8 * S._sun(hour, 3.0) ** 0.7 + 1.2 * np.cos(2 * np.pi * (hour - 16) / 24) \
        + rng.normal(0, 0.3, len(hour))
    dTm = np.gradient(Tm.reshape(NSTN, -1), axis=1).ravel()
    W = np.exp(rng.normal(1.0, 0.45, len(hour)))
    wdir = (rng.normal(200, 60, days)[day.astype(int)] + rng.normal(0, 25, len(hour))) % 360
    return {"sun": base, "sun_at": (lambda tau: S._sun(hour, tau)), "Tm": Tm, "dTm": dTm, "dayid": day * 10 + stn, "W": W,
            "cloud": rng.uniform(0, 1, len(hour)), "Td": rng.normal(20, 2, len(hour)), "station": stn, "hour": hour,
            "fveg": site(rng.uniform(0.1, 0.6, NSTN)), "fbuilt": site(rng.uniform(0.1, 0.6, NSTN)),
            "firrig": site([0.0, 0.7, 0.0, 0.3, 0.0, 0.0]), "coast_km": site([1.0, 30, 4.0, 60, 15, 45]),
            "onshore": ((np.abs(((wdir - 120) + 180) % 360 - 180) < 60)).astype(float),
            "rel_elev": site([0.0, -60, 40, -25, 120, 0.0]), "hw": site([0.3, 0.0, 2.5, 0.2, 0.0, 1.2]), "runway": site([0, 0, 0, 0, 0, 1])}


def flood_ctx(rng, n=60):
    c = S.flood_ctx(rng, n)
    SL = np.exp(rng.uniform(np.log(0.0003), np.log(0.05), n)); c["SL"] = SL
    c.update(tide=rng.normal(0.4, 0.5, n), coastal=(rng.random(n) < 0.4).astype(float), soilT=rng.normal(3, 5, n),
             burn=np.clip(rng.normal(0.0, 0.4, n), 0, 1) * (rng.random(n) < 0.35), karst=(rng.random(n) < 0.25).astype(float))
    return c


def wind_ctx(rng, n=500):
    c = S.wind_ctx(rng, n)
    c.update(valley=(rng.random(n) < 0.4).astype(float), water_frac=np.clip(rng.beta(0.6, 2.0, n), 0, 1))
    return c


# ------------------------------------------------------------------------------------------------ environment terms
HEAT_ENV = [
    Term("irrigation_cooling", "heat", "a · sin(e)⁺ · f_irrigated", "irrigated fields evaporate and cool the day (farmland bias)",
         ("sun", "firrig"), lambda c: c["sun"] * c["firrig"]),
    Term("sea_breeze", "heat", "a · [12–19 h] · [onshore] · e^(−d/L)", "an afternoon sea breeze cools the coast",
         ("hour", "onshore", "coast_km"), lambda c: (lambda L: ((c["hour"] >= 12) & (c["hour"] <= 19)) * c["onshore"] * np.exp(-c["coast_km"] / L)),
         family=(1.0, 40.0, True), unit="km"),
    Term("cold_pool", "heat", "a · [night] / U · (−Δz)⁺", "calm nights pool cold air in low ground", ("sun", "W", "rel_elev"),
         lambda c: (c["sun"] <= 0) / c["W"] * np.maximum(-c["rel_elev"], 0) / 50),
    Term("canyon_trapping", "heat", "a · [night] · H/W", "deep street canyons hold heat at night", ("sun", "hw"), lambda c: (c["sun"] <= 0) * c["hw"]),
    Term("runway_siting", "heat", "a · sin(e)⁺ · [runway]", "a station beside hot pavement reads warm by day", ("sun", "runway"), lambda c: c["sun"] * c["runway"]),
]
FLOOD_ENV = [
    Term("flat_ground", "flood", "b · max(0, ln(S₀ / S))", "on very flat ground water spreads and ponds more than modelled", ("SL",),
         lambda c: (lambda s0: np.maximum(0.0, np.log(s0 / c["SL"]))), family=(0.0005, 0.01, True), unit="slope"),
    Term("tide_backwater", "flood", "b · tide · [coastal outfall]", "high tide backs water up the pipes", ("tide", "coastal"), lambda c: c["tide"] * c["coastal"]),
    Term("frozen_soil", "flood", "b · [soil T < 0]", "frozen ground lets no water in", ("soilT",), lambda c: (c["soilT"] < 0).astype(float)),
    Term("burned_soil", "flood", "b · burn severity", "burned soil repels water", ("burn",), lambda c: c["burn"]),
    Term("karst_drain", "flood", "b · [karst]", "sinkholes and fractures drain water away", ("karst",), lambda c: c["karst"]),
]
WIND_ENV = [
    Term("valley_channel", "wind", "a · [valley] · cos(θ − φ)", "a valley steers and speeds wind along its axis", ("valley", "dirn"),
         lambda c: (lambda phi: c["valley"] * np.cos(np.radians(c["dirn"] - phi))), family=(0.0, 175.0, False), unit="°"),
    Term("water_fetch", "wind", "b · share of fetch over water", "wind off open water is less slowed than modelled", ("water_frac",), lambda c: c["water_frac"]),
]

ENVS = {  # name: (domain, ctx, base library, env term, planted law [(term, coef, alpha)], noise sd, grid)
    "Irrigated farmland (heat)": ("heat", heat_ctx, "irrigation_cooling", [("irrigation_cooling", -2.5, None)], 0.8),
    "Coastal sea breeze (heat)": ("heat", heat_ctx, "sea_breeze", [("sea_breeze", -3.0, 6.0)], 0.8),
    "Valley cold-air pools (heat)": ("heat", heat_ctx, "cold_pool", [("cold_pool", -1.0, None)], 0.8),
    "High-rise canyons at night (heat)": ("heat", heat_ctx, "canyon_trapping", [("canyon_trapping", 0.8, None)], 0.8),
    "Station beside a runway (heat)": ("heat", heat_ctx, "runway_siting", [("runway_siting", 2.0, None)], 0.8),
    "Very flat floodplain (flood)": ("flood", flood_ctx, "flat_ground", [("flat_ground", 0.15, 0.003)], 0.08),
    "Tidal backwater (flood)": ("flood", flood_ctx, "tide_backwater", [("tide_backwater", 0.25, None)], 0.08),
    "Frozen soil (flood)": ("flood", flood_ctx, "frozen_soil", [("frozen_soil", 0.2, None)], 0.08),
    "Burned soil after wildfire (flood)": ("flood", flood_ctx, "burned_soil", [("burned_soil", 0.4, None)], 0.08),
    "Karst sinks (flood)": ("flood", flood_ctx, "karst_drain", [("karst_drain", -0.15, None)], 0.08),
    "Valley channelling (wind)": ("wind", wind_ctx, "valley_channel", [("valley_channel", 0.25, 40.0)], 0.15),
    "Onshore wind off water (wind)": ("wind", wind_ctx, "water_fetch", [("water_fetch", 0.3, None)], 0.15),
}
GRID = {"heat": (0.1, -20, 20), "flood": (0.005, -3, 3), "wind": (0.01, -4, 4)}


def library(domain, ctx, drop=None):
    base = {"heat": T.HEAT, "flood": T.flood_terms(ctx), "wind": T.WIND}[domain]
    env = {"heat": HEAT_ENV, "flood": FLOOD_ENV, "wind": WIND_ENV}[domain]
    return [t for t in list(base) + env if t.name != drop]


def plant(terms, ctx, law):
    sig = 0.0
    for name, coef, alpha in law:
        t = next(x for x in terms if x.name == name); f = t.build(ctx); v = np.asarray(f(alpha) if alpha is not None else f, float)
        sig = sig + (v if t.fixed else coef * (v - v.mean()))
    return sig


def one(args):
    env, mode, seed = args
    dom, mk, term, law, sd = ENVS[env]; rng = np.random.default_rng(seed)
    ctx = mk(rng); n = len(ctx["sun"] if dom == "heat" else ctx["h"] if dom == "flood" else ctx["lp"])
    full = library(dom, ctx)
    y = 0.1 + rng.normal(0, sd, n) + (plant(full, ctx, law) if mode != "control" else 0)
    terms = library(dom, ctx, drop=term if mode == "unnamed" else None)
    prev = None
    if dom == "heat":
        prev = np.arange(n) - 1; prev[np.arange(0, n, 24)] = -1
    res = E.run_case(terms, ctx, y, E.Grid(*GRID[dom]), prev=prev)
    rem = y - res["_mu"]
    return {"env": env, "mode": mode, "seed": seed, "found": res["winner"]["terms"], "alpha": res["winner"]["alpha"],
            "flagged": bool(res["winner"]["terms"]), "named": term in res["winner"]["terms"],
            "bits_saved": res["bits_saved_per_reading"], "sd_before": res["sd_before"], "sd_after": res["sd_after"]}


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    J = [(e, m, s) for e in ENVS for m in ("flag", "control", "unnamed") for s in range(10)]
    t0 = time.time()
    with Pool() as p: rows = p.map(one, J, chunksize=2)
    summ = {}
    for e in ENVS:
        r = {m: [x for x in rows if x["env"] == e and x["mode"] == m] for m in ("flag", "control", "unnamed")}
        stand = {}
        for x in r["unnamed"]:
            for t in x["found"]: stand[t] = stand.get(t, 0) + 1
        summ[e] = {"term": ENVS[e][2], "flagged": np.mean([x["flagged"] for x in r["flag"]]), "named": np.mean([x["named"] for x in r["flag"]]),
                   "false_flag": np.mean([x["flagged"] for x in r["control"]]),
                   "unnamed_flagged": np.mean([x["flagged"] for x in r["unnamed"]]), "unnamed_stand_ins": dict(sorted(stand.items(), key=lambda kv: -kv[1])),
                   "bits_saved_flag": float(np.mean([x["bits_saved"] for x in r["flag"]])),
                   "sd_explained_flag": float(np.mean([1 - x["sd_after"] / x["sd_before"] for x in r["flag"]])),
                   "sd_explained_unnamed": float(np.mean([1 - x["sd_after"] / x["sd_before"] for x in r["unnamed"]]))}
        print(f"{e:38s} flagged {summ[e]['flagged']:.1f} named {summ[e]['named']:.1f} false {summ[e]['false_flag']:.1f} | unnamed flagged {summ[e]['unnamed_flagged']:.1f} via {list(stand)[:3]}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "envtests.json").write_text(json.dumps({"summary": summ, "runs": rows}, indent=1, default=float))
    print(f"{len(J)} runs {time.time() - t0:.0f} s")
