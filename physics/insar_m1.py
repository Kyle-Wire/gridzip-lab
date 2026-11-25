"""Protocol 10 A3, part D: method 1 (OPERA DISP-S1, GNSS tie off) against method 2 (V2) at Pawnee.

    python -m physics.insar_m1            # writes results/physics/insar_m1_vs_m2.json
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
from physlib import guards as Q, engine as E
from physlib.insar import INSAR_RATE
from physlib.engine import Term
from physics.insar_m2 import load as load_m2, design, rsd, folds_for, blocks_for, summary, _clean, GRID

ROOT = Path(__file__).resolve().parents[1]
M1 = ROOT / "data/insar_m1/m1_pawnee_points.npz"
OUT = ROOT / "results/physics/insar_m1_vs_m2.json"


def load_m1(lat0, lon0):
    d = np.load(M1); D = {k: d[k] for k in d.files}
    D["east_km"] = (D["lon"] - lon0) * 111.32 * np.cos(np.radians(lat0)); D["north_km"] = (D["lat"] - lat0) * 110.57
    s = D["series_0p1mm"].astype(float) / 10; s[D["series_0p1mm"] == -32768] = np.nan; D["series_mm"] = s
    D["t_yr"] = np.array([(np.datetime64(f"{x[:4]}-{x[4:6]}-{x[6:]}") - np.datetime64("2021-10-01")).astype(int) / 365.25 for x in D["dates"]])
    D["valid"] = np.isfinite(D["rate_mm_yr"])
    return D


def tilemap(e, n, r, step=1.0):
    ex = np.arange(-16, 16.01, step); ny = np.arange(-16.5, 16.01, step)
    M = np.full((len(ny) - 1, len(ex) - 1), np.nan)
    ix = np.digitize(e, ex) - 1; iy = np.digitize(n, ny) - 1; key = iy * 1000 + ix
    for k in np.unique(key):
        m = key == k
        if m.sum() >= 30 and 0 <= k // 1000 < M.shape[0] and 0 <= k % 1000 < M.shape[1]: M[k // 1000, k % 1000] = np.median(r[m])
    return ex, ny, M


def shared_rate(D, mask_dates):
    """Per-point rate (line + annual) on the given dates, each point using its own finite dates (>= 20)."""
    t = D["t_yr"][mask_dates]; Y = D["series_mm"][:, mask_dates]; W = np.isfinite(Y).astype(float); Y0 = np.where(W > 0, Y, 0.0)
    X = design(t)
    A = np.einsum("pd,di,dj->pij", W, X, X); b = np.einsum("pd,di->pi", W * Y0, X)
    r = np.full(len(Y), np.nan); good = W.sum(1) >= 20
    r[good] = np.linalg.solve(A[good], b[good][..., None])[:, 1, 0]
    return r


def spread_short(e, n, r, a=0.2, b=0.4, N=400000, seed=0):
    rng = np.random.default_rng(seed); i = rng.integers(0, len(r), N); j = rng.integers(0, len(r), N)
    d = np.hypot(e[i] - e[j], n[i] - n[j]); m = (d >= a) & (d < b)
    return rsd(r[i][m] - r[j][m]) / np.sqrt(2)


def lagrms(Y, lags=(1, 10, 40, 120)):
    Y = Y - np.nanmedian(Y, 0); return {int(k): float(np.sqrt(np.nanmean((Y[:, k:] - Y[:, :-k]) ** 2))) for k in lags}


def station_rates(D, mask_dates, gnss_demean=True):
    t = D["t_yr"][mask_dates]; X = design(t); st = [str(s) for s in D["stations"]]; out = {}
    for i, s in enumerate(st):
        rad = D["station_radar_mean_200m_mm"][i][mask_dates]; g = D["station_gnss_los_mm"][i][mask_dates]
        m = np.isfinite(rad) & np.isfinite(g)
        out[s] = {"radar": float(np.linalg.lstsq(X[m], rad[m], rcond=None)[0][1]), "gnss": float(np.linalg.lstsq(X[m], g[m], rcond=None)[0][1]), "n": int(m.sum())}
    return out


def pair_share(sr, tol=2.0):
    st = list(sr); ok = tot = 0; rows = []
    for a in range(len(st)):
        for b in range(a + 1, len(st)):
            err = (sr[st[a]]["radar"] - sr[st[b]]["radar"]) - (sr[st[a]]["gnss"] - sr[st[b]]["gnss"])
            rows.append((st[a], st[b], err)); tot += 1; ok += abs(err) <= tol
    return ok, tot, rows


def main():
    D2 = load_m2(); lat0, lon0 = float(np.nanmean(D2["lat"])), float(np.nanmean(D2["lon"]))
    D1 = load_m1(lat0, lon0)
    d2 = list(D2["dates"]); sh1 = np.isin(D1["dates"], d2); sh2 = np.isin(D2["dates"], D1["dates"][sh1])
    out = {"shared_dates": int(sh1.sum())}
    v1, v2 = D1["valid"], D2["valid"]
    r1 = shared_rate({**D1, "series_mm": D1["series_mm"][v1]}, sh1); r2 = shared_rate({**D2, "series_mm": D2["series_mm"][v2]}, sh2)
    e1, n1, e2, n2 = D1["east_km"][v1], D1["north_km"][v1], D2["east_km"][v2], D2["north_km"][v2]
    f1 = np.isfinite(r1); r1, e1, n1 = r1[f1], e1[f1], n1[f1]; v1 = np.where(v1)[0][f1]
    # D1 corner
    def corner(e, n, r):
        nw = np.median(r[(e < -4) & (n > 4)]); rest = np.median(r[(e > -2) & (n < 4)]); return float(nw), float(rest), float(nw - rest)
    c1, c2 = corner(e1, n1, r1), corner(e2, n2, r2)
    verdict = "ground or common physics" if c1[2] <= -2 else ("our processing" if abs(c1[2]) < 1 else "unclear")
    _, _, M1t = tilemap(e1, n1, r1); _, _, M2t = tilemap(e2, n2, r2); m = np.isfinite(M1t) & np.isfinite(M2t)
    out["D1"] = {"m1_nw_rest": c1, "m2_nw_rest": c2, "verdict": verdict, "tile_corr": float(np.corrcoef(M1t[m], M2t[m])[0, 1]), "tiles": int(m.sum()),
                 "m1_tile_rsd": rsd(M1t[np.isfinite(M1t)] - np.nanmedian(M1t)), "m2_tile_rsd": rsd(M2t[np.isfinite(M2t)] - np.nanmedian(M2t))}
    # D2 stations
    s1, s2 = station_rates(D1, sh1), station_rates(D2, sh2)
    p1, p2 = pair_share(s1), pair_share(s2)
    out["D2"] = {"m1": s1, "m2": s2, "m1_pairs_ok": [p1[0], p1[1]], "m2_pairs_ok": [p2[0], p2[1]],
                 "m1_pairs": [(a, b, round(x, 2)) for a, b, x in p1[2]], "m2_pairs": [(a, b, round(x, 2)) for a, b, x in p2[2]]}
    # D3, D4
    out["D3"] = {"m1_sigma_median": float(np.nanmedian(D1["sigma_mm_yr"][v1])), "m1_spread_0p2_0p4": spread_short(e1, n1, r1),
                 "m2_sigma_median": float(np.nanmedian(D2["sigma_mm_yr"][v2])), "m2_spread_0p2_0p4": spread_short(e2, n2, r2),
                 "m1_rate_rsd": rsd(r1), "m2_rate_rsd": rsd(r2)}
    out["D4"] = {"m1": lagrms(D1["series_mm"][v1][::4][:, sh1]), "m2": lagrms(D2["series_mm"][v2][::4][:, sh2])}
    print(json.dumps(_clean(out), indent=0)[:3000], flush=True)
    # D5 compressor
    rows = np.sort(np.random.default_rng(10).choice(len(r1), 3000, replace=False))
    ex, ny, M2full = tilemap(e2, n2, r2)
    ix = np.clip(np.digitize(e1[rows], ex) - 1, 0, M2full.shape[1] - 1); iy = np.clip(np.digitize(n1[rows], ny) - 1, 0, M2full.shape[0] - 1)
    m2tile = M2full[iy, ix]; m2tile = np.where(np.isfinite(m2tile), m2tile, np.nanmedian(M2full))
    g = lambda k: np.asarray(D1[k][v1][rows], float)
    ctx = {"coherence": g("temporal_coherence"), "incidence": g("incidence_deg"), "east": e1[rows].astype(float), "north": n1[rows].astype(float),
           "annual_amp": np.nan_to_num(g("annual_amplitude_mm"), nan=float(np.nanmedian(g("annual_amplitude_mm")))), "jumps": g("n_jump_dates"), "m2_tile": m2tile}
    keep = [t for t in INSAR_RATE if all(k in ctx for k in t.needs)]
    keep.append(Term("m2_tile", "insar", "b · method 2's 1 km tile rate here", "the two methods share this structure", ("m2_tile",), lambda c: c["m2_tile"] - c["m2_tile"].mean()))
    y = r1[rows].astype(float); tiles = (np.floor((ctx["east"] + 16) / 3) * 1000 + np.floor((ctx["north"] + 16.5) / 3)).astype(int)
    t0 = time.time()
    res = Q.analyse(keep, ctx, y, E.Grid(0.05, -150, 150), grammar=True, folds=folds_for(ctx), groups=tiles, bootstrap_B=10, null_N=9,
                    blocks=blocks_for(tiles), name="Pawnee method 1 rate field", seed=10)
    s = summary(res); s["seconds"] = time.time() - t0; s["library"] = [t.name for t in keep]
    out["D5"] = s; print("D5", s["winner"], s["sd_before"], s["sd_after"], s["bootstrap"]["selected"] if s.get("bootstrap") else None, s["null"], flush=True)
    OUT.write_text(json.dumps(_clean(out), indent=1))
    np.savez_compressed(ROOT / "results/physics/insar_m1_tiles.npz", M1=M1t, M2=M2t)


if __name__ == "__main__":
    main()
