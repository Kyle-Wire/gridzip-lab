"""Planted-law synthetic cases for each domain: realistic inputs, a known law added to noise, so the search can be
scored on whether it names what was planted (and only that), and on pure noise, where it must name nothing."""
import numpy as np

from . import terms as T


# ---------------------------------------------------------------------------------------------- inputs per domain
def wind_ctx(rng, n=360):
    cls = rng.choice(3, n, p=[0.45, 0.35, 0.2])
    return {"lp": rng.normal(0, 0.35, n), "z": rng.uniform(4, 30, n), "z0": np.exp(rng.uniform(np.log(0.02), np.log(1.0), n)),
            "zref": 10.0, "fetch": rng.uniform(0, 12, n), "along": (cls == 0).astype(float), "cross": (cls == 1).astype(float),
            "hc": rng.choice([0.0, 1.0, 2.0], n), "corner": (rng.random(n) < 0.12).astype(float), "kq": rng.uniform(0, 0.6, n),
            "stn_loo": rng.normal(0, 0.12, n), "dz": rng.normal(0, 40, n), "dirn": rng.uniform(0, 360, n),
            "hour": rng.uniform(0, 24, n), "month": rng.integers(1, 13, n).astype(float)}


def _sun(hour, tau=0.0):
    """A clear July day at 40° N, local solar time: sin(elevation)+, with the curve shifted tau hours later."""
    h = np.asarray(hour) - tau
    return np.maximum(0.0, 0.2 + 0.75 * np.cos(2 * np.pi * (h - 12) / 24))


def heat_ctx(rng, days=5, stations=3):
    hrs = np.arange(days * 24, dtype=float); hour = np.tile(hrs % 24, stations); day = np.tile(hrs // 24, stations)
    stn = np.repeat(np.arange(stations), days * 24)
    sun = _sun(hour)
    lvl = rng.normal(0, 1.5, days)[day.astype(int)]
    Tm = 25 + lvl + 8 * _sun(hour, 3.0) ** 0.7 + 1.2 * np.cos(2 * np.pi * (hour - 16) / 24) + rng.normal(0, 0.3, len(hour))
    dTm = np.gradient(Tm.reshape(stations, -1), axis=1).ravel()
    fveg = np.array([0.6, 0.2, 0.35])[:stations][stn]; fbuilt = np.array([0.1, 0.5, 0.3])[:stations][stn]
    return {"sun": sun, "sun_at": (lambda tau: _sun(hour, tau)), "Tm": Tm, "dTm": dTm, "dayid": day * 10 + stn,
            "W": np.exp(rng.normal(1.0, 0.45, len(hour))), "cloud": rng.uniform(0, 1, len(hour)), "Td": rng.normal(20, 2, len(hour)),
            "fveg": fveg, "fbuilt": fbuilt, "station": stn, "hour": hour}


def flood_ctx(rng, n=40):
    h = rng.uniform(0.005, 0.06, n); V = rng.uniform(0.2, 2.2, n); SL = np.exp(rng.uniform(np.log(0.004), np.log(0.08), n))
    return {"h": h, "V": V, "SL": SL, "Fr": V / np.sqrt(9.81 * h), "Eo": rng.uniform(0.3, 0.9, n), "T": rng.uniform(0.5, 3, n),
            "Q": np.exp(rng.normal(np.log(0.01), 0.7, n)), "rain": rng.choice([50.0, 100.0, 200.0], n)}


def ground_ctx(rng=None):
    t = np.arange(2016.5, 2025.9, 12 / 365.25)
    return {"t": t}


def ar1(rng, n, sd, rho, blocks=1):
    """AR(1) noise with the given marginal sd, restarted in each block (station, cell)."""
    e = rng.normal(0, sd, n); out = np.empty(n); m = n // blocks
    for b in range(blocks):
        s = slice(b * m, (b + 1) * m if b < blocks - 1 else n); x = e[s]; z = np.empty_like(x); z[0] = x[0]
        for i in range(1, len(x)): z[i] = rho * z[i - 1] + np.sqrt(1 - rho ** 2) * x[i]
        out[s] = z
    return out


# ---------------------------------------------------------------------------------------------- planted laws
def _lib(domain, ctx):
    if domain == "flood": return T.flood_terms(ctx)
    if domain == "ground": return T.ground_terms(ctx)
    return T.DOMAINS[domain]


def _term(domain, ctx, name):
    return next(t for t in _lib(domain, ctx) if t.name == name)


def plant(domain, ctx, law):
    """law: list of (term name, coefficient or list, family parameter or None). Returns the planted signal."""
    sig = 0.0
    for name, coef, alpha in law:
        t = _term(domain, ctx, name); f = t.build(ctx)
        v = np.asarray(f(alpha) if alpha is not None else f, float)
        if v.ndim == 1: v = v - (0 if t.fixed or t.family else v.mean())
        if v.ndim == 2: v = (v - v.mean(0)) @ np.atleast_1d(coef)
        elif not t.fixed: v = coef * v
        sig = sig + v
    return sig


SCENARIOS = {
    "wind": dict(ctx=wind_ctx, sd=0.15, grid=(0.01, -4, 4), laws={
        "contrast + entrance decay": [("contrast", -0.3, None), ("fetch_decay", 0.35, 2.0)],
        "roughness length x3": [("log_law_z0", None, 3.0)],
        "station siting + direction": [("station_bias", 1.0, None), ("direction", [0.12, 0.0], None)],
    }),
    "heat": dict(ctx=heat_ctx, sd=0.8, grid=(0.1, -20, 20), laws={
        "swing too large + runs late": [("amplitude", -0.22, None), ("timing_shift", 0.8, None)],
        "sunlight 1 h late": [("radiation_lag", 2.5, 1.0)],
        "calm nights + station siting": [("night_decoupling", -2.0, None), ("station_offsets", [0.8, -0.6], None)],
    }),
    "flood": dict(ctx=flood_ctx, sd=0.05, grid=(0.005, -3, 3), laws={
        "weir drowns to orifice": [("orifice_transition", None, 0.025)],
        "splash-over on steep grades": [("splash_over", -0.45, 1.1), ("grade", -3.0, None)],
        "datum offset 3 mm": [("datum", None, 0.003)],
    }),
    "ground": dict(ctx=ground_ctx, sd=3.0, grid=(0.01, -200, 200), laws={
        "sinkhole onset 2022": [("onset", -8.0, 2022.0)],
        "two reference jumps": [("jump_a", 6.0, 2019.3), ("jump_b", -5.0, 2023.1)],
        "seasonal + consolidation": [("seasonal", [3.0, 1.0], None), ("consolidation", 12.0, 2.0)],
    }),
}


def case(domain, law_name, seed, scale=1.0, noise="iid", rho=0.8):
    """A synthetic case: (terms, ctx, y, grid, planted term names)."""
    sc = SCENARIOS[domain]; rng = np.random.default_rng(seed)
    ctx = sc["ctx"](rng); n = len(next(v for v in ctx.values() if isinstance(v, np.ndarray)))
    law = sc["laws"][law_name] if law_name else []
    sig = plant(domain, ctx, law) * scale if law else np.zeros(n)
    blocks = 3 if domain == "heat" else 1
    e = rng.normal(0, sc["sd"], n) if noise == "iid" else ar1(rng, n, sc["sd"], rho, blocks)
    y = 0.1 + sig + e
    return _lib(domain, ctx), ctx, y, sc["grid"], [nm for nm, _, _ in law]


def prev_of(domain, n):
    """A1: ordered domains (heat: hours within each station; ground: dates) get the AR(1) residual code."""
    if domain not in ("heat", "ground"): return None
    blocks = 3 if domain == "heat" else 1; m = n // blocks; p = np.arange(n) - 1
    p[np.arange(0, n, m)] = -1
    return p


# ---------------------------------------------------------------------------------------------- v2 (PROTOCOL-7)
def conv_ctx(rng, n=400, stations=8):
    """Synthetic thunderstorm gust rows with physically linked inputs: the pressure jump and rain rise with the
    temperature drop, DCAPE with dry air aloft, the stability with the hour."""
    stn = rng.integers(0, stations, n); site = rng.normal(0, 0.12, stations)
    dT = rng.gamma(2.0, 2.0, n); Tk = rng.normal(300, 4, n)
    dp = 100 * np.maximum(0.0, 0.25 * dT + rng.normal(0, 0.8, n))
    dd700 = rng.gamma(2.0, 4.0, n); dcape = np.exp(rng.normal(6.2, 0.4, n) + 0.03 * dd700)
    hour = rng.uniform(0, 24, n); zL = np.clip(-0.4 * np.cos(2 * np.pi * (hour - 15) / 24) + rng.normal(0, 0.25, n), -2, 1)
    return {"dT": dT, "Tk": Tk, "dp": dp, "dcape": dcape, "cape": np.exp(rng.normal(7.0, 0.7, n)), "lapse": rng.normal(7.0, 1.0, n),
            "dd700": dd700, "rain_mmh": np.exp(rng.normal(1.5, 1.0, n) + 0.1 * dT), "u_model": rng.uniform(12, 30, n),
            "zL": zL, "z": np.full(n, 10.0), "z0": np.exp(rng.uniform(np.log(0.01), np.log(0.5), n)), "ri_b": rng.normal(0, 0.3, n),
            "u_storm": rng.uniform(3, 20, n), "w_star": rng.uniform(0.5, 3.0, n), "hour": hour, "month": rng.integers(4, 10, n).astype(float),
            "dirn": rng.uniform(0, 360, n), "station": stn, "stn_loo": site[stn] + rng.normal(0, 0.03, n)}


def trop_ctx(rng, n=400):
    """Synthetic landfall rows: hours since landfall (negative offshore), radius over R_max, azimuth from the motion."""
    return {"t_land": rng.uniform(-12, 48, n), "r_rmax": np.exp(rng.uniform(np.log(0.3), np.log(8.0), n)),
            "theta_rel": rng.uniform(0, 360, n), "u_storm": rng.uniform(2, 12, n), "u_model": rng.uniform(15, 60, n),
            "land_frac": np.clip(rng.beta(1.2, 1.2, n), 0, 1), "z": np.full(n, 10.0), "z0": np.exp(rng.uniform(np.log(0.0002), np.log(0.5), n)),
            "hour": rng.uniform(0, 24, n), "dirn": rng.uniform(0, 360, n)}
