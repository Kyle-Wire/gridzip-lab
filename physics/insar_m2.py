"""Protocol 13: what method 2 (insar_pipeline's own MiaplPy -> snaphu -> MintPy route) is missing on stable ground, V2 Pawnee.

    python -m physics.insar_m2 [synth|A|B|C|all]

Reads data/insar_m2/v2_pawnee_points.npz (README there), writes results/physics/insar_m2.json.
Written fresh for SPARC / gridzip-lab.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from physlib import engine as E, guards as Q
from physlib.engine import Term

ROOT = Path(__file__).resolve().parents[1]
NPZ = ROOT / "data/insar_m2/v2_pawnee_points.npz"
OUT = ROOT / "results/physics/insar_m2.json"
GRID = E.Grid(0.05, -60.0, 60.0)
TILE_KM = 3.0


def _c(a): a = np.asarray(a, np.float64); return a - a.mean()


# ------------------------------------------------------------------------------------------------------- data
def load():
    d = np.load(NPZ)
    D = {k: d[k] for k in d.files}
    lat0, lon0 = float(np.nanmean(D["lat"])), float(np.nanmean(D["lon"]))
    D["east_km"] = (D["lon"] - lon0) * 111.32 * np.cos(np.radians(lat0))
    D["north_km"] = (D["lat"] - lat0) * 110.57
    D["t_yr"] = np.array([(np.datetime64(f"{s[:4]}-{s[4:6]}-{s[6:]}") - np.datetime64("2021-10-01")).astype(int) / 365.25
                          for s in D["dates"]])
    D["series_mm"] = D["series_0p1mm"].astype(np.float64) / 10.0
    D["series_mm"][D["series_0p1mm"] == np.iinfo(np.int16).min] = np.nan
    D["valid"] = np.isfinite(D["rate_mm_yr"]) & D["in_conncomp_mask"] & np.all(np.isfinite(D["series_mm"]), axis=1)
    tx = np.floor((D["east_km"] - D["east_km"].min()) / TILE_KM).astype(int)
    ty = np.floor((D["north_km"] - D["north_km"].min()) / TILE_KM).astype(int)
    D["tile"] = tx * 1000 + ty
    return D


def sample(D, seed, n=20_000, exclude=None):
    idx = np.where(D["valid"])[0]
    if exclude is not None: idx = np.setdiff1d(idx, exclude)
    return np.sort(np.random.default_rng(seed).choice(idx, size=min(n, len(idx)), replace=False))


# ----------------------------------------------------------------------------------------------- part A terms
def ctx_for(D, rows):
    g = lambda k: np.asarray(D[k][rows], np.float64)
    c = {"coherence": g("temporal_coherence"), "ps_share": g("ps_fraction"), "amp_disp": g("amp_dispersion_median"),
         "outside_comp": g("n_ifgrams_outside_component") / 375.0, "closure_amb": np.log1p(g("n_triplets_nonzero_ambiguity")),
         "jumps": g("n_jump_dates"), "height": _c(g("dem_height_m")), "incidence": g("incidence_deg"),
         "ref_dist": g("dist_to_reference_m") / 1000.0, "east": g("east_km"), "north": g("north_km"),
         "annual_amp": g("annual_amplitude_mm")}
    for k, v in list(c.items()):                     # a NaN covariate (amp_disp where no PS pixel) -> its median
        if not np.all(np.isfinite(v)): v = v.copy(); v[~np.isfinite(v)] = np.nanmedian(v); c[k] = v
    return c


from physlib.insar import INSAR_RATE as TERMS  # the library lives in physlib/insar.py


def folds_for(c):
    e, n = c["east"], c["north"]; me, mn = np.median(e), np.median(n)
    W, Ea, S, N = np.where(e < me)[0], np.where(e >= me)[0], np.where(n < mn)[0], np.where(n >= mn)[0]
    return [("west -> east", W, Ea), ("east -> west", Ea, W), ("south -> north", S, N), ("north -> south", N, S)]


def blocks_for(tiles):
    return [np.where(tiles == t)[0] for t in np.unique(tiles)]


def summary(res):
    gd = res["guards"]
    return {"winner": res["winner"]["terms"], "sd_before": res["sd_before"], "sd_after": res["sd_after"],
            "bits_saved_per_reading": res["bits_saved_per_reading"], "transfers": res["transfers"],
            "holdout": res["holdout"], "inclusion": gd.get("inclusion"), "bootstrap": gd.get("bootstrap"),
            "null": {k: v for k, v in gd.get("null_bits_saved", {}).items() if k != "all"},
            "neighbour": gd.get("neighbour"), "verdict": Q.verdict(res), "ladder": res["ladder"],
            "coefs": res["winner"].get("theta")}


def run_A(D, seed, rows=None, y=None, B=10, N=9, grammar=True):
    rows = sample(D, seed) if rows is None else rows
    c = ctx_for(D, rows); y = np.asarray(D["rate_mm_yr"][rows], np.float64) if y is None else y
    tiles = D["tile"][rows]
    t = time.time()
    res = Q.analyse(TERMS, c, y, GRID, grammar=grammar, folds=folds_for(c), groups=tiles, bootstrap_B=B, null_N=N,
                    blocks=blocks_for(tiles), name=f"V2 rate field seed {seed}", seed=seed, boot_by_group=True)
    out = summary(res); out["seconds"] = time.time() - t; out["n"] = len(rows)
    return out, res


def bootshare(s, term):
    bs = s.get("bootstrap") or {}
    sel = bs.get("selected") or bs.get("share") or {}
    return float(sel.get(term, 0.0)) if isinstance(sel, dict) else 0.0


def named_A(s1, s2):
    out = []
    for t in s1["winner"]:
        if t not in s2["winner"]: continue
        ok_folds = sum(1 for h in s1["holdout"] if h.get("transfers"))
        if ok_folds >= 3 and bootshare(s1, t) >= 0.6: out.append(t)
    return out


# --------------------------------------------------------------------------------------------- part B series
def design(t, extra=None):
    X = [np.ones_like(t), t, np.sin(2 * np.pi * t), np.cos(2 * np.pi * t)]
    if extra is not None: X += list(extra)
    return np.column_stack(X)


def step_cols(D):
    ms = D["ministack"].astype(int)
    return [(ms >= k).astype(float) for k in range(1, ms.max() + 1)]


def fit_rates(Y, X):
    """Y: points x dates. Least squares per point (shared X). Returns rates (col 1), residual sd per point, rss."""
    beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None)
    R = Y.T - X @ beta
    return beta[1], R, beta


def rsd(v): v = v[np.isfinite(v)]; return float(1.4826 * np.median(np.abs(v - np.median(v))))


def gauss_bits(R, k_extra, n_dates):
    """Total bits per point of the residuals under N(0, sd per point), plus 0.5 log2(n) per extra coefficient."""
    sd = np.sqrt(np.mean(R ** 2, axis=0)) + 1e-9
    return 0.5 * n_dates * np.log2(2 * np.pi * np.e * sd ** 2) + k_extra * 0.5 * np.log2(n_dates)


def run_B(D, rows, Y=None):
    t = D["t_yr"]; Y = D["series_mm"][rows] if Y is None else Y; n = len(t)
    X0 = design(t); r0, R0, _ = fit_rates(Y, X0); b0 = gauss_bits(R0, 0, n)
    out = {"plain": {"rate_rsd": rsd(r0), "rate_sd": float(np.std(r0)), "rate_median": float(np.median(r0))}}
    # B1 ministack steps
    S = step_cols(D); X1 = design(t, S); r1, R1, be1 = fit_rates(Y, X1); b1 = gauss_bits(R1, len(S), n)
    steps = be1[4:]
    out["B1_ministack_steps"] = {"rate_rsd": rsd(r1), "rate_sd": float(np.std(r1)), "bits_saved_per_point": float(np.mean(b0 - b1)),
                                 "step_rsd_mm": [rsd(s) for s in steps], "step_median_mm": [float(np.median(s)) for s in steps],
                                 "corr_rate_plain_vs_sum_steps": float(np.corrcoef(r0, steps.sum(0))[0, 1])}
    # B2 common mode per tile per date
    tiles = D["tile"][rows]; Yc = Y.copy()
    for tl in np.unique(tiles):
        m = tiles == tl
        Yc[m] -= np.median(D["series_mm"][D["valid"] & (D["tile"] == tl)], axis=0)
    r2, R2, _ = fit_rates(Yc, X0); b2 = gauss_bits(R2, 0, n)   # per-tile medians come from many points: not charged per point
    out["B2_common_mode"] = {"rate_rsd": rsd(r2), "rate_sd": float(np.std(r2)), "bits_saved_per_point": float(np.mean(b0 - b2)),
                             "tile_rate_rsd": rsd(r0 - r2)}
    # B3 fading: short rate per ministack (line on its own 10 dates)
    ms = D["ministack"].astype(int); short = []
    for k in range(ms.max() + 1):
        m = ms == k
        if m.sum() < 5: continue
        Xk = np.column_stack([np.ones(m.sum()), t[m]]); bk, *_ = np.linalg.lstsq(Xk, Y[:, m].T, rcond=None); short.append(bk[1])
    short = np.array(short)                                     # ministacks x points
    within = short.mean(0)                                      # mean of within-ministack rates
    out["B3_fading"] = {"median_short_rate_minus_whole": float(np.median(within - r0)),
                        "short_rate_median_per_ministack": [float(np.median(s)) for s in short],
                        "rsd_mean_short_rate": rsd(within), "corr_whole_vs_mean_short": float(np.corrcoef(r0, within)[0, 1])}
    # B4 late start
    m = t >= 1.0; r4, R4, _ = fit_rates(Y[:, m], design(t[m]))
    out["B4_late_start"] = {"rate_rsd": rsd(r4), "rate_sd": float(np.std(r4)), "corr_with_plain": float(np.corrcoef(r0, r4)[0, 1])}
    # combined B1 + B2
    r12, R12, _ = fit_rates(Yc, X1)
    out["B1_plus_B2"] = {"rate_rsd": rsd(r12), "rate_sd": float(np.std(r12)), "bits_saved_per_point": float(np.mean(b0 - gauss_bits(R12, len(S), n)))}
    for k in ("B1_ministack_steps", "B2_common_mode"):
        o = out[k]; o["rsd_cut"] = 1 - o["rate_rsd"] / out["plain"]["rate_rsd"]
        o["named"] = bool(o["bits_saved_per_point"] > 0 and o["rsd_cut"] >= 0.20)
    out["_rates"] = {"plain": r0, "B1": r1, "B2": r2, "B12": r12, "B4": r4, "within": within}
    return out


# ------------------------------------------------------------------------------------------------ part C GNSS
def run_C(D):
    t = D["t_yr"]; X0 = design(t); X1 = design(t, step_cols(D)); out = {}
    for i, s in enumerate(D["stations"]):
        rad, gn = D["station_radar_mean_200m_mm"][i], D["station_gnss_los_mm"][i]
        m = np.isfinite(rad) & np.isfinite(gn)
        if m.sum() < 20: out[str(s)] = {"n": int(m.sum())}; continue
        dlt = rad[m] - gn[m]
        rate = lambda X: float(np.linalg.lstsq(X[m], dlt, rcond=None)[0][1])
        r_rad = float(np.linalg.lstsq(X0[m], rad[m], rcond=None)[0][1]); r_gn = float(np.linalg.lstsq(X0[m], gn[m], rcond=None)[0][1])
        out[str(s)] = {"n": int(m.sum()), "radar_rate": r_rad, "gnss_rate": r_gn, "diff_rate_plain": rate(X0), "diff_rate_steps": rate(X1),
                       "diff_rate_late": float(np.linalg.lstsq(X0[m & (t >= 1)], (rad - gn)[m & (t >= 1)], rcond=None)[0][1]) if (m & (t >= 1)).sum() > 20 else None,
                       "diff_rms_mm": float(np.sqrt(np.mean((dlt - dlt.mean()) ** 2)))}
    return out


# ------------------------------------------------------------------------------------------------ synthetic
def run_synth(D, seeds=(0, 1, 2, 3, 4)):
    rows = sample(D, 10, n=3000); c = ctx_for(D, rows)
    base = np.asarray(D["rate_mm_yr"][rows], float)
    # noise: a smooth field + white, of the plain residual's size, so spatial correlation is kept
    out = {"a_null": [], "b_coherence": [], "c_steps": []}
    for s in seeds:
        rng = np.random.default_rng(100 + s)
        y0 = rng.permutation(base)            # covariates' link broken, marginal kept
        sa, _ = run_A(D, s, rows=rows, y=y0, B=0, N=0)
        out["a_null"].append(sa["winner"])
        coh = _c(c["coherence"]); k = np.sqrt(0.10 / 0.90) * np.std(y0) / np.std(coh)
        sb, _ = run_A(D, s, rows=rows, y=y0 + k * coh, B=0, N=0)
        out["b_coherence"].append(sb["winner"])
        Y = D["series_mm"][rows[:3000]].copy(); S = np.array(step_cols(D))
        Y += (rng.normal(0, 2.0, (3000, S.shape[0])) @ S)
        b = run_B(D, rows[:3000], Y=Y); b0 = run_B(D, rows[:3000])
        out["c_steps"].append({"cut": b["B1_ministack_steps"]["rsd_cut"], "cut_without_plant": b0["B1_ministack_steps"]["rsd_cut"],
                               "plain_rsd_with_plant": b["plain"]["rate_rsd"], "plain_rsd_without": b0["plain"]["rate_rsd"]})
    out["pass_a"] = sum(1 for w in out["a_null"] if not [t for t in w if not t.startswith("g:")] and not w) >= 4
    out["pass_b"] = sum(1 for w in out["b_coherence"] if "coherence" in w or any("coherence" in t for t in w)) >= 4
    out["pass_c"] = sum(1 for z in out["c_steps"] if z["cut"] >= 0.20) >= 4
    return out


def _clean(o):
    if isinstance(o, dict): return {str(k): _clean(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, (list, tuple)): return [_clean(v) for v in o]
    if isinstance(o, np.ndarray): return _clean(o.tolist())
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.bool_,)): return bool(o)
    return o


def main(part="all"):
    D = load(); res = json.loads(OUT.read_text()) if OUT.exists() else {}
    if part in ("synth", "all"):
        t = time.time(); res["synthetic"] = run_synth(D); print("synth", {k: v for k, v in res["synthetic"].items() if k.startswith("pass")}, f"{time.time()-t:.0f}s", flush=True)
    if part in ("B", "all"):
        rows = sample(D, 10); b = run_B(D, rows); res["B"] = b; print("B", json.dumps(_clean(b))[:1500], flush=True)
        np.savez_compressed(ROOT / "results/physics/insar_m2_rates.npz", rows=rows, **b["_rates"])
    if part in ("C", "all"):
        res["C"] = run_C(D); print("C", json.dumps(_clean(res["C"]), indent=0)[:2000], flush=True)
    if part in ("A", "all"):
        r1 = sample(D, 10, n=3000); r2 = sample(D, 11, n=3000, exclude=r1)
        s1, _ = run_A(D, 10, rows=r1); print("A1", s1["winner"], s1["sd_before"], s1["sd_after"], f"{s1['seconds']:.0f}s", flush=True)
        s2, _ = run_A(D, 11, rows=r2); print("A2", s2["winner"], s2["sd_before"], s2["sd_after"], f"{s2['seconds']:.0f}s", flush=True)
        res["A"] = {"seed10": s1, "seed11": s2, "named": named_A(s1, s2)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(_clean(res), indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
