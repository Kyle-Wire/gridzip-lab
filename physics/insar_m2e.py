"""Protocol 10 A4, part E: which processing step makes the OK02 corner (V2 Pawnee)."""
import json, time, numpy as np
from pathlib import Path
from physlib import guards as Q
from physlib.insar import INSAR_RATE
from physlib.engine import Term
from physics.insar_m2 import load, design, rsd, ctx_for, folds_for, blocks_for, summary, _clean, GRID, sample

ROOT = Path(__file__).resolve().parents[1]; DD = ROOT / "data/insar_m2"


def rates(Y, t):
    W = np.isfinite(Y).astype(float); Y0 = np.where(W > 0, Y, 0.0); X = design(t)
    A = np.einsum("pd,di,dj->pij", W, X, X); b = np.einsum("pd,di->pi", W * Y0, X)
    return np.linalg.solve(A, b[..., None])[:, 1, 0]


def main():
    D = load(); v = D["valid"]; t = D["t_yr"]; e, n = D["east_km"], D["north_km"]
    mm = lambda a: np.where(a == -32768, np.nan, a.astype(float) / 10)
    pre = mm(np.load(DD / "v2_series_before_hrrr.npz")["series_0p1mm"]); hr = mm(np.load(DD / "v2_tropo.npz")["hrrr_correction_0p1mm"])
    C = np.load(DD / "v2_components.npz"); fin = D["series_mm"]
    dem = pre - fin - hr
    r = {k: rates(Y[v], t) for k, Y in (("before", pre), ("hrrr", hr), ("dem", dem), ("final", fin))}
    ev, nv = e[v], n[v]; NW = (ev < -4) & (nv > 4); RE = (ev > -2) & (nv < 4)
    cor = {k: float(np.median(x[NW]) - np.median(x[RE])) for k, x in r.items()}
    out = {"E1": {"nw_minus_rest": cor, "makes_corner": [k for k in ("hrrr", "dem") if np.sign(cor[k]) == np.sign(cor["final"]) and abs(cor[k]) >= 0.5 * abs(cor["final"])],
                  "before_has_corner": bool(np.sign(cor["before"]) == np.sign(cor["final"]) and abs(cor["before"]) >= 0.5 * abs(cor["final"]))}}
    apart = C["n_ifgrams_apart"][v].astype(float) / 375
    bits = np.unpackbits(C["apart_from_reference_bits"][v], axis=1)[:, :375]
    late_if = np.array([d[0] >= "20231101" for d in C["ifgram_dates"]])
    na = bits.sum(1); apart_late = np.where(na > 0, bits[:, late_if].sum(1) / np.maximum(na, 1), 0.5)
    never = apart == 0
    out["E2"] = {"share_never_apart": float(never.mean()), "nw_apart_median": float(np.median(apart[NW])), "rest_apart_median": float(np.median(apart[RE])),
                 "nw_minus_rest_never_apart": float(np.median(r["final"][NW & never]) - np.median(r["final"][RE & never])),
                 "n_nw_never": int((NW & never).sum()), "n_rest_never": int((RE & never).sum())}
    out["E2"]["corner_survives_never_apart"] = bool(np.sign(out["E2"]["nw_minus_rest_never_apart"]) == np.sign(cor["final"]) and abs(out["E2"]["nw_minus_rest_never_apart"]) >= 0.5 * abs(cor["final"]))
    print(json.dumps(out, indent=1), flush=True)
    # E3
    idx_valid = np.where(v)[0]; rows = sample(D, 10, n=3000); pos = np.searchsorted(idx_valid, rows)
    c = ctx_for(D, rows); c.update({"hrrr_rate": r["hrrr"][pos], "dem_rate": r["dem"][pos], "apart": apart[pos], "apart_late": apart_late[pos]})
    extra = [Term(k, "insar", f"b · {k}", m, (k,), (lambda kk: lambda cc: cc[kk] - cc[kk].mean())(k)) for k, m in
             (("hrrr_rate", "the weather correction's own rate leaks"), ("dem_rate", "the DEM-error step's own rate leaks"),
              ("apart", "unwrapping islands apart from the reference"), ("apart_late", "islands concentrated late in the record"))]
    y = np.asarray(D["rate_mm_yr"][rows], float); tiles = D["tile"][rows]; t0 = time.time()
    res = Q.analyse(INSAR_RATE + extra, c, y, GRID, grammar=True, folds=folds_for(c), groups=tiles, bootstrap_B=10, null_N=9, blocks=blocks_for(tiles), name="V2 corner part E", seed=10)
    s = summary(res); s["seconds"] = time.time() - t0; out["E3"] = s
    print("E3", s["winner"], s["sd_before"], s["sd_after"], s["bootstrap"]["selected"], s["null"], flush=True)
    (ROOT / "results/physics/insar_m2_partE.json").write_text(json.dumps(_clean(out), indent=1))


if __name__ == "__main__":
    main()
