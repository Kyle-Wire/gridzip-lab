"""The curated term library: known physical terms per domain, each a formula a model's errors might follow.

A term reads named inputs from a case's context (ctx) and is skipped when they are missing, so one library serves
every case of its domain. Forms are the standard textbook ones (neutral log law, power law, weir and orifice laws,
OHM storage hysteresis, hinge and step models for InSAR series); constants that a model might have wrong are left
free as a coefficient or a one-parameter family. Residual conventions (fixed per domain):
  wind    r = ln(measured / model speed)
  heat    r = measured - model air temperature (K)
  flood   r = ln(measured or reference / model) for an exchange flow or capture share
  ground  r = measured series - the chain's fitted model (mm), one cell's time series, or truth - chain rate (mm/yr)
"""
import numpy as np

from .engine import Term

G0 = 9.81


def _c(a): a = np.asarray(a, np.float64); return a - a.mean()


def _harm(x, period):
    w = 2 * np.pi * np.asarray(x, np.float64) / period; return np.column_stack([np.sin(w), np.cos(w)])


def _onehot(labels):
    labels = np.asarray(labels); u = sorted(set(labels.tolist()))
    return np.column_stack([(labels == v).astype(float) for v in u[1:]]) if len(u) > 1 else np.zeros((len(labels), 0))


# ============================================================================================================ wind
WIND = [
    Term("contrast", "wind", "β · ln(U_model / U_ref)", "the model exaggerates (β<0) or flattens (β>0) its own speed-ups",
         ("lp",), lambda c: c["lp"]),
    Term("log_law_z0", "wind", "ln[ ln(z / (a·z0)) / ln(z / z0) ]", "the roughness length is wrong by a factor a",
         ("z", "z0"), lambda c: (lambda a: np.log(np.log(c["z"] / (a * c["z0"])) / np.log(c["z"] / c["z0"])))
         , family=(0.1, 10.0, True), fixed=True, group="roughness", unit="× z0"),
    Term("power_law", "wind", "Δα · ln(z / z_ref)", "the profile exponent is wrong", ("z", "zref"),
         lambda c: np.log(c["z"] / c["zref"]), group="profile"),
    Term("displacement", "wind", "ln[ ln((z − d) / z0) / ln(z / z0) ]", "a missing zero-plane displacement d",
         ("z", "z0"), lambda c: (lambda d: np.log(np.log(np.maximum(c["z"] - d, 1.01 * c["z0"]) / c["z0"]) / np.log(c["z"] / c["z0"]))),
         family=(0.0, 6.0, False), fixed=True, unit="m"),
    Term("gust_factor", "wind", "g · 1 / ln(z / z0)", "gusts scale with turbulence intensity I ≈ 1/ln(z/z0) (peak factor g wrong)",
         ("z", "z0"), lambda c: 1.0 / np.log(c["z"] / c["z0"]), group="roughness"),
    Term("fetch_decay", "wind", "A · exp(−x / L)", "an entrance or internal boundary layer: the miss decays with distance x",
         ("fetch",), lambda c: (lambda L: np.exp(-c["fetch"] / L)), family=(0.25, 16.0, True), group="fetch", unit="L"),
    Term("fetch_linear", "wind", "b · x", "the miss grows steadily with distance into the array", ("fetch",), lambda c: c["fetch"], group="fetch"),
    Term("canyon_along", "wind", "δ · [street along the wind]", "channelling along streets is off", ("along",), lambda c: c["along"]),
    Term("canyon_cross", "wind", "δ · [street across the wind]", "sheltering in cross streets is off", ("cross",), lambda c: c["cross"]),
    Term("height_ratio", "wind", "b · H_c / H", "a taller neighbour changes the flow more than the model says", ("hc",), lambda c: c["hc"]),
    Term("corner", "wind", "δ · [near a tall corner]", "corner jets / downwash are off", ("corner",), lambda c: c["corner"]),
    Term("contrast_along", "wind", "β · [along] · ln(U/U_ref)", "speed-up contrast differs in along-wind streets",
         ("lp", "along"), lambda c: c["lp"] * c["along"]),
    Term("fetch_decay_along", "wind", "A · [along] · exp(−x/L)", "entrance decay only in along-wind streets",
         ("fetch", "along"), lambda c: (lambda L: c["along"] * np.exp(-c["fetch"] / L)), family=(0.25, 16.0, True), group="fetch_along", unit="L"),
    Term("turbulence_k", "wind", "½ ln(1 + a·(k-ratio))", "the measured speed includes turbulence the model's mean leaves out",
         ("kq",), lambda c: (lambda a: 0.5 * np.log1p(a * c["kq"])), family=(0.0, 3.0, False), fixed=True),
    Term("station_bias", "wind", "δ · (station's mean miss in its other events)", "a per-station siting or instrument offset",
         ("stn_loo",), lambda c: c["stn_loo"]),
    Term("elevation", "wind", "b · Δz_station", "station height above its surroundings (hill speed-up or mountain sheltering)",
         ("dz",), lambda c: c["dz"]),
    Term("direction", "wind", "a·sin θ + b·cos θ", "the miss depends on wind direction (upwind terrain or exposure)",
         ("dirn",), lambda c: _harm(c["dirn"], 360.0)),
    Term("stability", "wind", "a·sin(2πh/24) + b·cos(2πh/24)", "day/night stability changes gusts the neutral model misses",
         ("hour",), lambda c: _harm(c["hour"], 24.0)),
    Term("season", "wind", "a·sin(2πm/12) + b·cos(2πm/12)", "a seasonal miss (storm type, foliage roughness)",
         ("month",), lambda c: _harm(c["month"], 12.0)),
]

# ============================================================================================================ heat
HEAT = [
    Term("radiation", "heat", "a · sin(e)⁺", "daytime heating by sunlight is too strong or weak", ("sun",), lambda c: c["sun"], group="sun"),
    Term("radiation_lag", "heat", "a · sin(e)⁺ at t − τ", "sunlight's heating arrives late or early by τ hours", ("sun_at",),
         lambda c: (lambda tau: c["sun_at"](tau)), family=(-3.0, 6.0, False), group="sun", unit="h"),
    Term("storage_hysteresis", "heat", "a · d sin(e)⁺ / dt", "heat storage releases with a phase lag (OHM a2 term)", ("sun_at",),
         lambda c: (c["sun_at"](-0.5) - c["sun_at"](0.5))),
    Term("amplitude", "heat", "a · (T_model − daily mean)", "the model's daily swing is too large (a<0) or small (a>0)", ("Tm", "dayid"),
         lambda c: c["Tm"] - _groupmean(c["Tm"], c["dayid"])),
    Term("timing_shift", "heat", "δ · dT_model / dt", "the model's daily cycle runs late (δ>0, h) or early", ("dTm",), lambda c: c["dTm"]),
    Term("wind_mixing", "heat", "a · sin(e)⁺ / U", "daytime mixing: light-wind days heat more than modelled", ("sun", "W"), lambda c: c["sun"] / c["W"]),
    Term("night_decoupling", "heat", "a · [night] / U", "calm nights decouple and cool more", ("sun", "W"), lambda c: (c["sun"] <= 0) / c["W"]),
    Term("sky_cooling", "heat", "a · [night] · (1 − cloud)", "clear nights radiate more heat away", ("sun", "cloud"),
         lambda c: (c["sun"] <= 0) * (1 - c["cloud"])),
    Term("cloud_shading", "heat", "a · sin(e)⁺ · cloud", "clouds cut sunlight more or less than modelled", ("sun", "cloud"), lambda c: c["sun"] * c["cloud"]),
    Term("humidity", "heat", "a · (T_dew − mean)", "moisture: evaporation or latent heat off", ("Td",), lambda c: _c(c["Td"])),
    Term("vegetation_et", "heat", "a · sin(e)⁺ · f_veg", "evapotranspiration by trees and grass off by day", ("sun", "fveg"), lambda c: c["sun"] * c["fveg"]),
    Term("built_storage", "heat", "a · [night] · f_built", "roofs and paving release stored heat at night", ("sun", "fbuilt"),
         lambda c: (c["sun"] <= 0) * c["fbuilt"]),
    Term("station_offsets", "heat", "δ_s per station", "a siting offset at each station", ("station",), lambda c: _onehot(c["station"])),
    Term("diurnal_shape", "heat", "a·sin(2πh/24) + b·cos(2πh/24)", "a daily shape with no physics name (a gap the library can't name)",
         ("hour",), lambda c: _harm(c["hour"], 24.0), physics=False),
]


def _groupmean(v, gid):
    v = np.asarray(v, np.float64); gid = np.asarray(gid); out = np.empty_like(v)
    for g in set(gid.tolist()): m = gid == g; out[m] = v[m].mean()
    return out


# ============================================================================================================ flood
FLOOD = [
    Term("weir_exponent", "flood", "β · ln(h / h̄)", "flow follows a different power of depth than h^1.5", ("h",), lambda c: np.log(c["h"] / np.mean(c["h"]))),
    Term("orifice_transition", "flood", "−ln(h / h_t) for h > h_t", "the inlet drowns: weir becomes orifice (Q ∝ h^0.5) above h_t",
         ("h",), lambda c: (lambda ht: -np.log(np.maximum(c["h"], ht) / ht)), family=None, fixed=True, unit="m"),
    Term("datum", "flood", "1.5 · ln((h + δ) / h)", "the depth is measured from the wrong level by δ", ("h",),
         lambda c: (lambda d: 1.5 * np.log(np.maximum(c["h"] + d, 1e-6) / c["h"])), family=None, fixed=True, unit="m"),
    Term("approach_froude", "flood", "b · Fr", "approach velocity carries water past (or into) the inlet", ("Fr",), lambda c: c["Fr"]),
    Term("velocity_head", "flood", "b · ln(1 + V² / 2gh)", "the inlet sees the energy head, not only the depth", ("V", "h"),
         lambda c: np.log1p(c["V"] ** 2 / (2 * G0 * c["h"]))),
    Term("grade", "flood", "b · S_L", "the miss grows with street grade (shallower, faster flow)", ("SL",), lambda c: c["SL"], group="grade"),
    Term("manning_slope", "flood", "b · ln S_L", "the depth-slope exponent of Manning's law is off", ("SL",), lambda c: np.log(c["SL"]), group="grade"),
    Term("splash_over", "flood", "b · max(0, V − V0)", "above a splash-over speed V0 water jumps the grate", ("V",),
         lambda c: (lambda v0: np.maximum(0.0, c["V"] - v0)), family=None, unit="m/s"),
    Term("frontal_share", "flood", "b · E_o", "the frontal share of gutter flow (HEC-22) drives capture", ("Eo",), lambda c: c["Eo"]),
    Term("spread", "flood", "b · T", "flow spread across the street", ("T",), lambda c: c["T"]),
    Term("flow", "flood", "b · ln Q", "the miss scales with flow (a rating-curve shape error)", ("Q",), lambda c: np.log(c["Q"])),
    Term("rain_intensity", "flood", "b · ln i", "the miss depends on rain rate", ("rain",), lambda c: np.log(c["rain"])),
    Term("manning_n", "flood", "0.6 · ln(n_true / n_model)", "roughness n wrong: depth ∝ n^0.6 (wide channel)", ("n_ratio_log",),
         lambda c: 0.6 * c["n_ratio_log"]),
]


def flood_terms(ctx):
    """Families whose range depends on the data (transition depth, datum, splash-over speed)."""
    out = []
    for t in FLOOD:
        t = Term(**{**t.__dict__})
        if t.name == "orifice_transition" and "h" in ctx:
            h = np.asarray(ctx["h"]); t.family = (float(np.quantile(h, 0.1)), float(np.quantile(h, 0.9)), True)
        if t.name == "datum" and "h" in ctx:
            h = np.asarray(ctx["h"]); t.family = (-0.5 * float(h.min()), 1.0 * float(h.min()), False)
        if t.name == "splash_over" and "V" in ctx:
            v = np.asarray(ctx["V"]); t.family = (float(np.quantile(v, 0.1)), float(np.quantile(v, 0.9)), False)
        out.append(t)
    return out


# ------------------------------------------------------------------------------------------- flood at real marks
# Protocol 11 (flood wall check). Residual r = model minus mark (m) at USGS high-water marks, rule D's error. Context
# comes from the run and its inputs only (CHI tools/gate2/walls_terms.py), never from other marks.
FLOOD_WALL = [
    Term("volume_scale", "flood", "b · h", "the miss grows with the model's depth: too much or too little water (rain "
         "total, losses) or a roughness that scales every depth", ("h_model",), lambda c: c["h_model"]),
    Term("depth_power", "flood", "b · ln(h + 0.1)", "shallow and deep water are off by different factors (a depth-"
         "dependent roughness or infiltration error)", ("h_model",), lambda c: np.log(c["h_model"] + 0.1)),
    Term("floodplain_height", "flood", "b · HAND", "the miss grows with height above the local low: the channel carries "
         "too little (missing bathymetry under lidar water) or too much", ("hand",), lambda c: c["hand"], group="valley"),
    Term("channel_distance", "flood", "b · ln(1 + d_ch / 10 m)", "the miss changes away from the main flow path "
         "(overbank exchange, floodplain storage)", ("d_channel",), lambda c: np.log1p(c["d_channel"] / 10.0), group="valley"),
    Term("edge_loss", "flood", "b · exp(−d_edge / L)", "water leaves or piles up at the model's edge within L of it "
         "(boundary condition, box too small)", ("d_edge",), lambda c: (lambda L: np.exp(-c["d_edge"] / L)),
         family=(30.0, 600.0, True), unit="m"),
    Term("blockage", "flood", "b · building share within 30 m", "dense buildings raise or lower water more than the "
         "model's walls do (conveyance, storage lost)", ("bld_frac",), lambda c: c["bld_frac"]),
    Term("ground_slope", "flood", "b · S", "the miss tracks ground slope (steep, fast flow: friction or timing)",
         ("slope",), lambda c: c["slope"]),
    Term("velocity_head_wall", "flood", "b · V² / 2g", "marks on upstream faces record run-up the model's level omits",
         ("V",), lambda c: c["V"] ** 2 / (2 * G0)),
    Term("datum_kind", "flood", "b · [height-above-ground mark]", "marks read as height above ground and marks with a "
         "NAVD88 level disagree: a datum or ground-reading error", ("kind_hag",), lambda c: c["kind_hag"]),
    Term("on_building", "flood", "b · [mark on a building]", "marks on walls differ from open-ground marks (wall "
         "treatment, local run-up or shelter)", ("at_building",), lambda c: c["at_building"]),
    # protocol 15 (27 Nov 2026): after the thin-film fix, the Lumber marks still read ~0.6 m low
    Term("storage_fill", "flood", "b · d_ch / (d_ch + L)", "the floodplain far from the river has not filled to the "
         "river's level when the run stops (a steady run cut short, or too little volume reaching the swamp)",
         ("d_channel",), lambda c: (lambda L: c["d_channel"] / (c["d_channel"] + L)), family=(30.0, 1000.0, True), unit="m",
         group="valley"),
    Term("spill_rise", "flood", "b · (spill − ground)", "water must climb over a ridge (road, railway, levee) to reach "
         "the mark; the real flood got through culverts the model does not have", ("spill_rise",),
         lambda c: c["spill_rise"]),
    Term("town_offset", "flood", "per-town offset", "a whole town is off by one amount (its rain product or inflow)",
         ("site_onehot",), lambda c: c["site_onehot"], group="offset", physics=False),
]

# Protocol 11, buildings: target = the model's depth at an observed-dry building's wall (should be 0). CHI
# tools/gate2/walls_bld.py writes the context.
FLOOD_WALL_BLD = [
    Term("contributing_area", "flood", "b · ln(1 + A / 100 m²)", "uphill runoff piles against the walls: the model's "
         "walls block flow the real ground lets past (gaps, yards, drains)", ("contrib",), lambda c: np.log1p(c["contrib"] / 100.0)),
    Term("roof_drop", "flood", "b · footprint / ring area", "too much roof rain lands at the walls (FloodSim spreads "
         "each roof's rain over the cells touching it)", ("roof_drop",), lambda c: c["roof_drop"]),
    Term("bld_slope", "flood", "b · S", "water banks against walls on steep ground", ("slope",), lambda c: c["slope"]),
    Term("bld_hand", "flood", "b · HAND", "houses high above the stream get wet (rain on the ground, not the river)",
         ("hand",), lambda c: c["hand"]),
    Term("bld_edge", "flood", "b · exp(−d_edge / L)", "water piles up or leaves near the model's edge", ("d_edge",),
         lambda c: (lambda L: np.exp(-c["d_edge"] / L)), family=(30.0, 600.0, True), unit="m"),
    Term("bld_footprint", "flood", "b · ln(footprint / 100 m²)", "big buildings collect more water (their own size)",
         ("footprint",), lambda c: np.log(np.maximum(c["footprint"], 1.0) / 100.0)),
]


# ------------------------------------------------------------------------------------------- depth (protocol 12)
# Planted FloodSim depth cases (FloodSim benchmarks/physics/depth). Target y = ln(h_model / h_true).
FLOOD_DEPTH = [
    Term("raster_width", "flood", "−0.6 · ln(Wm / W)", "the channel the model sees is wider or narrower than the real one "
         "(cells snap to the footprint); depth ∝ width^−0.6 at a given flow", ("Wm", "W"),
         lambda c: -0.6 * np.log(c["Wm"] / c["W"]), fixed=True),
    Term("wall_friction", "flood", "b · ln(1 + 2h/W)", "walls drag on the flow and the 2D model has bed friction only "
         "(expected b ≈ −0.4)", ("h", "W"), lambda c: np.log1p(2 * c["h"] / c["W"])),
    Term("cells_across", "flood", "b · dx / Wm", "too few cells across the channel to carry the flow", ("Wm", "dx"),
         lambda c: c["dx"] / c["Wm"]),
    Term("edge_drawdown", "flood", "b · exp(−x_e / (λ·Lb))", "a free edge drains the river within λ backwater lengths",
         ("x_e", "Lb"), lambda c: (lambda lam: np.exp(-c["x_e"] / (lam * c["Lb"]))), family=(0.05, 2.0, True)),
    Term("froude", "flood", "b · Fr", "the miss grows with flow speed", ("Fr",), lambda c: c["Fr"]),
    Term("sheet_scale", "flood", "b · ln x", "sheet-flow depth grows along the slope with the wrong power", ("x",),
         lambda c: np.log(c["x"])),
    Term("slope_term", "flood", "b · ln S", "the miss tracks slope", ("S",), lambda c: np.log(c["S"])),
    Term("cell_size", "flood", "b · ln dx", "the miss tracks cell size", ("dx",), lambda c: np.log(c["dx"])),
    Term("roughness", "flood", "b · ln n", "the miss tracks roughness", ("n",), lambda c: np.log(c["n"])),
]


# ============================================================================================================ ground
def _span(c): t = np.asarray(c["t"]); return float(t.min()), float(t.max())


GROUND = [
    Term("linear", "ground", "v · (t − t̄)", "the rate is wrong (reference offset or real motion)", ("t",), lambda c: c["t"]),
    Term("seasonal", "ground", "a·sin 2πt + b·cos 2πt", "an annual cycle (groundwater, thermal) is in the series", ("t",), lambda c: _harm(c["t"], 1.0)),
    Term("semiannual", "ground", "a·sin 4πt + b·cos 4πt", "a half-year cycle", ("t",), lambda c: _harm(c["t"], 0.5)),
    Term("onset", "ground", "v₂ · max(0, t − t₀)", "motion starts mid-record at t₀ (a sinkhole or new load)", ("t",),
         lambda c: (lambda t0: np.maximum(0.0, c["t"] - t0)), family=None, unit="year"),
    Term("jump_a", "ground", "s · H(t − t_j)", "a step on one date (reference pixel jump, unwrapping error)", ("t",),
         lambda c: (lambda tj: (c["t"] >= tj).astype(float)), family=None, unit="year"),
    Term("jump_b", "ground", "s · H(t − t_k)", "a second step on another date", ("t",),
         lambda c: (lambda tj: (c["t"] >= tj).astype(float)), family=None, unit="year"),
    Term("consolidation", "ground", "a · exp(−(t − t_start) / τ)", "settling that slows with time (consolidation after loading)", ("t",),
         lambda c: (lambda tau: np.exp(-(c["t"] - c["t"].min()) / tau)), family=(0.3, 10.0, True), unit="year"),
    Term("acceleration", "ground", "a · (t − t̄)²", "the rate changes steadily", ("t",), lambda c: (c["t"] - c["t"].mean()) ** 2),
    Term("thermal", "ground", "a · T_air(t)", "thermal expansion of structures follows air temperature", ("temp",), lambda c: c["temp"]),
]


def ground_terms(ctx):
    out = []
    lo, hi = _span(ctx)
    for t in GROUND:
        t = Term(**{**t.__dict__})
        if t.name in ("onset", "jump_a", "jump_b"): t.family = (lo + 0.3, hi - 0.3, False)
        out.append(t)
    return out


# Rate-level ground terms (one value per cell: truth - chain rate), from the chain's own outputs only
GROUND_RATE = [
    Term("late_rate", "ground", "b · (rate over last third − rate)", "the straight line misses a late change (onset)", ("late",), lambda c: c["late"]),
    Term("rate_change", "ground", "b · (rate 2nd half − rate 1st half)", "the rate changed mid-record", ("accel",), lambda c: c["accel"]),
    Term("shrink", "ground", "b · (rate − area mean)", "cell-to-cell detail is partly noise: shrink it", ("vel",), lambda c: c["vel"]),
    Term("noise_level", "ground", "b · series sd", "noisier cells are biased", ("std",), lambda c: c["std"]),
    Term("seasonal_amp", "ground", "b · seasonal amplitude", "seasonal motion leaks into the rate", ("amp",), lambda c: c["amp"]),
    Term("near_house", "ground", "δ · [within 250 m of the house]", "a local feature at the property", ("near",), lambda c: c["near"]),
    Term("distance", "ground", "b · distance to the house", "a smooth gradient across the area", ("dist",), lambda c: c["dist"]),
]

DOMAINS = {"wind": WIND, "heat": HEAT, "flood": FLOOD, "flood_wall": FLOOD_WALL, "flood_wall_bld": FLOOD_WALL_BLD, "flood_depth": FLOOD_DEPTH, "ground": GROUND, "ground_rate": GROUND_RATE}


# ============================================================================================================ v2 (PROTOCOL-7)
# Convective and stability terms for storm gusts, r = ln(measured gust / model gust). Standard forms written from
# memory (density-current speed, Bernoulli pressure-jump speed, downdraft speed from DCAPE, Businger-Dyer psi_m);
# the source papers were not opened, so every constant is a free coefficient.
G_H = G0 * 1000.0        # g times a nominal 1 km cold-pool depth; the depth's error goes into the coefficient
RHO = 1.2


def psi_m(zl):
    """Businger-Dyer stability correction for momentum: unstable (z/L < 0) and stable (z/L >= 0, -5 z/L)."""
    zl = np.asarray(zl, np.float64); x = np.power(np.maximum(1 - 16 * np.minimum(zl, 0), 1e-12), 0.25)
    un = 2 * np.log((1 + x) / 2) + np.log((1 + x ** 2) / 2) - 2 * np.arctan(x) + np.pi / 2
    return np.where(zl < 0, un, -5 * np.minimum(zl, 2.0))


CONVECTIVE = [
    Term("cold_pool", "wind", "b · √(g·h·ΔT⁺/T) / U_model", "the storm's cold pool pushes the gust beyond the model (density current)",
         ("dT", "Tk", "u_model"), lambda c: np.sqrt(G_H * np.maximum(c["dT"], 0) / c["Tk"]) / c["u_model"], group="outflow"),
    Term("pressure_jump", "wind", "b · √(2Δp⁺/ρ) / U_model", "a mesohigh pressure jump drives the outflow", ("dp", "u_model"),
         lambda c: np.sqrt(2 * np.maximum(c["dp"], 0) / RHO) / c["u_model"], group="mesohigh"),
    Term("downdraft", "wind", "b · √(2·DCAPE) / U_model", "downdraft strength (downdraft CAPE) adds to the gust", ("dcape", "u_model"),
         lambda c: np.sqrt(2 * np.maximum(c["dcape"], 0)) / c["u_model"], group="downdraft"),
    Term("instability", "wind", "b · √(2·CAPE) / U_model", "updraft instability (CAPE) organises stronger storms", ("cape", "u_model"),
         lambda c: np.sqrt(2 * np.maximum(c["cape"], 0)) / c["u_model"], group="instability"),
    Term("lapse_rate", "wind", "b · Γ₀₋₃ₖₘ", "steep low-level lapse rates mix momentum down", ("lapse",), lambda c: c["lapse"]),
    Term("dry_aloft", "wind", "b · (T − T_d)₇₀₀", "dry air aloft evaporates rain and drives the downdraft", ("dd700",), lambda c: c["dd700"]),
    Term("rain_loading", "wind", "b · ln(1 + rain rate)", "precipitation drag loads the downdraft", ("rain_mmh",), lambda c: np.log1p(c["rain_mmh"])),
    Term("mo_stability", "wind", "b · ψ_m(z/L) / ln(z/z0)", "surface-layer stability (Monin-Obukhov) changes the profile", ("zL", "z", "z0"),
         lambda c: psi_m(c["zL"]) / np.log(c["z"] / c["z0"]), group="stability"),
    Term("richardson", "wind", "b · Ri_b", "bulk stability (Richardson number) damps or frees the gust", ("ri_b",), lambda c: c["ri_b"], group="stability"),
    Term("storm_motion", "wind", "b · U_storm / U_model", "the storm's own motion adds to the outflow", ("u_storm", "u_model"),
         lambda c: c["u_storm"] / c["u_model"]),
    Term("convective_velocity", "wind", "b · w* / U_model", "boundary-layer convection (Deardorff w*) drives the gusts", ("w_star", "u_model"),
         lambda c: c["w_star"] / c["u_model"]),
]

# Tropical cyclone terms, r = ln(measured / model wind). Kaplan-DeMaria-style decay after landfall, written from memory.
TROPICAL = [
    Term("land_decay", "wind", "b · e^(−α·t_land)", "the storm weakens after landfall at a rate α the model has wrong", ("t_land",),
         lambda c: (lambda a: np.exp(-a * np.maximum(c["t_land"], 0.0))), family=(0.01, 0.5, True), group="decay", unit="1/h"),
    Term("radial_shape", "wind", "b · ln(r / R_max) for r > R_max", "the wind falls off with radius faster or slower than modelled",
         ("r_rmax",), lambda c: np.log(np.maximum(c["r_rmax"], 1.0)), group="radial"),
    Term("inner_core", "wind", "δ · [r < R_max]", "inside the eyewall radius the model is off", ("r_rmax",), lambda c: (c["r_rmax"] < 1).astype(float)),
    Term("motion_asymmetry", "wind", "b · cos(θ − θ_motion) · U_storm / U_model", "the right-front quadrant is stronger (storm motion)",
         ("theta_rel", "u_storm", "u_model"), lambda c: np.cos(np.radians(c["theta_rel"])) * c["u_storm"] / c["u_model"]),
    Term("sea_land", "wind", "b · f_land upwind", "friction over land upwind cuts the wind more than modelled", ("land_frac",), lambda c: c["land_frac"]),
]

HEAT_V2 = [
    Term("mixed_layer", "heat", "a · sin(e)⁺ / z_i", "a shallow daytime mixed layer heats more", ("sun", "zi"), lambda c: c["sun"] / c["zi"]),
    Term("advection", "heat", "a · U·∂T/∂x", "warm or cold air carried in by the wind", ("adv",), lambda c: c["adv"]),
    Term("soil_moisture", "heat", "a · sin(e)⁺ · (θ − θ̄)", "wet soil evaporates and cools the day", ("sun", "swc"), lambda c: c["sun"] * _c(c["swc"])),
]

WIND_FULL = WIND + CONVECTIVE + TROPICAL
DOMAINS.update({"wind_full": WIND_FULL, "heat_full": HEAT + HEAT_V2})

# ============================================================================================================ nwp
# Errors of a learned (or any gridded) weather model, r = (forecast − analysis) / spread for one variable at one point
# and lead. Inputs are known at the start or fixed: lat, land (0-1), height_km, sun (cos zenith at valid time),
# msl_anom (hPa, from the zonal mean), msl_curv (hPa per cell², discrete Laplacian), spin850 (1e-5 s⁻¹),
# wind850 (m/s), q850 (g/kg), tend (the variable's last 6 h change). Added 24 Nov 2026 for WindSim Part F.
OMEGA = 7.2921e-5
NWP = [
    Term("diffusion", "nwp", "κ · ∇²x₀", "the model smooths highs and lows like extra diffusion (κ>0 fills lows, flattens highs)",
         ("msl_curv",), lambda c: c["msl_curv"], group="smoothing"),
    Term("vortex_spin_down", "nwp", "b · ζ₈₅₀ · max(0, U₈₅₀ − U₀)", "strong, windy vortices (storms) are spun down",
         ("spin850", "wind850"), lambda c: (lambda u0: c["spin850"] * np.maximum(c["wind850"] - u0, 0.0)),
         family=(0.0, 15.0, False), group="vortex", unit="m/s"),
    Term("vortex_abs", "nwp", "b · |ζ₈₅₀|", "errors grow with spin of either sign (cyclones and anticyclones)", ("spin850",),
         lambda c: np.abs(c["spin850"]), group="vortex"),
    Term("coriolis", "nwp", "b · f,  f = 2Ω sin φ", "wind-pressure balance is off by latitude", ("lat",),
         lambda c: 2 * OMEGA * np.sin(np.radians(c["lat"])) * 1e4, group="latitude"),
    Term("polar", "nwp", "b · max(0, |φ| − 60°)", "the lat-lon grid converges near the poles", ("lat",),
         lambda c: np.maximum(np.abs(c["lat"]) - 60.0, 0.0), group="latitude"),
    Term("terrain_reduction", "nwp", "b · [land] · z_s · p′", "sea-level pressure reduced over high ground is off", ("land", "height_km", "msl_anom"),
         lambda c: c["land"] * c["height_km"] * c["msl_anom"], group="land"),
    Term("land_anomaly", "nwp", "b · [land] · p′", "highs and lows behave differently over land", ("land", "msl_anom"),
         lambda c: c["land"] * c["msl_anom"], group="land"),
    Term("extrapolation", "nwp", "b · Δx₋₆ₕ", "the model carries the last change on too far (b>0) or too little (b<0)", ("tend",),
         lambda c: c["tend"]),
    Term("diurnal", "nwp", "b · cos Z⁺", "daytime heating is out of phase or too strong", ("sun",), lambda c: c["sun"]),
    Term("latent_heating", "nwp", "b · q₈₅₀ · max(ζ₈₅₀ · sign φ, 0)", "moist cyclones deepen from latent heat the model lacks",
         ("q850", "spin850", "lat"), lambda c: c["q850"] * np.maximum(c["spin850"] * np.sign(c["lat"] + 1e-9), 0.0)),
]
DOMAINS.update({"nwp": NWP})


# ============================================================================================================ InSAR (PROTOCOL-13)
from .insar import INSAR_RATE  # noqa: E402
DOMAINS["insar_rate"] = INSAR_RATE
