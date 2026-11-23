"""PROTOCOL-4: PHL-POC-001 house patch − wide ring (DISP-S1 LOS, after A11). python physics/phl_patch.py <series.json>"""
import json, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import runner
from codelen import lawsearch as L

OUT = Path(__file__).resolve().parents[1] / "results" / "physics"
G = L.Grid(0.05, -60.0, 60.0)


def runmed(x, k=2):
    return np.array([np.median(x[max(0, i - k):i + k + 1]) for i in range(len(x))])


def main(path):
    s = json.loads(Path(path).read_text())["routes"]["disp"]
    t = np.array(s["t"]); P, N, W = (np.array(s[k], float) for k in ("patch", "near ring", "wide ring"))
    ok = ~(np.isnan(P) | np.isnan(W)); dropped = t[~ok].tolist()          # A11 bad dates are blanked in the file
    t, P, N, W = t[ok], P[ok], np.where(np.isnan(N), W, N)[ok], W[ok]
    y = P - W; n = len(y); c = lambda a: a - a.mean(); idx = np.arange(n)
    wdev = W - runmed(W); wlag = np.r_[0.0, wdev[:-1]]; ar1 = np.r_[y[0], y[:-1]]
    lib = {"t": c(t), "step": L.Family(lambda a: c((t >= a).astype(float)), 2017.0, 2025.5, n=69),
           "annual": L.Family(lambda a: np.cos(2 * np.pi * (t - a)), 0.0, 0.975, n=40),
           "semi": L.Family(lambda a: np.cos(4 * np.pi * (t - a)), 0.0, 0.4875, n=40),
           "wdev": c(wdev), "wdev_lag": c(wlag), "ar1": c(ar1), "nr": c(N - W),
           "spike": L.Family(lambda a: c((idx == int(round(a))).astype(float)), 0, n - 1, n=n)}
    h = n // 2
    res = runner.case(lib, y, G, [("first -> second", idx[:h], idx[h:]), ("second -> first", idx[h:], idx[:h])],
                      exclusive=[{"annual", "semi"}], n_hold=8)
    order = idx.copy(); ctx = np.full((n, 2), -1); ctx[1:, 0] = idx[:-1]; ctx[2:, 1] = idx[:-2]
    nb, _ = runner.neighbour_test(y, order, ctx, G, [idx], n_perm=99)
    d = y[1:-1] - (y[:-2] + y[2:]) / 2; rsd = 1.4826 * np.median(np.abs(d - np.median(d))); z = d / rsd
    j = int(np.argmin(np.abs(t - 2019.51))); k = j - 1
    spike = {"robust_sd_mm": float(rsd), "z_2019_51": float(z[k]), "date": float(t[j]), "rel_mm": float(y[j]),
             "n_abs_z_over_5": int(np.sum(np.abs(z) > 5)), "dates_over_5": [(float(t[i + 1]), float(z[i])) for i in np.where(np.abs(z) > 5)[0]],
             "z_sorted_top": sorted([(float(abs(v)), float(t[i + 1])) for i, v in enumerate(z)], reverse=True)[:6],
             "gaussian_expected_max_abs_z": float(math.sqrt(2 * math.log(2 * len(z))))}
    for name, terms in (("spike only", ["spike"]), ("annual only", ["annual"]), ("wdev only", ["wdev"]), ("t only", ["t"]), ("step only", ["step"]), ("ar1 only", ["ar1"])):
        r = L.score(lib, terms, y, G); spike.setdefault("single_terms", {})[name] = {"bpr": r["bits_per_reading"], "coef": r["fit_coef"], "alpha": r["fit_alpha"]}
    res["neighbour_prev2"] = nb; res["dropped_dates"] = dropped; res["spike_test"] = spike
    res["data"] = {"t": t.tolist(), "y": y.tolist(), "wdev": wdev.tolist(), "mu_best": res["_mu"].tolist(), "W": W.tolist(), "P": P.tolist(), "N": N.tolist()}
    print("n", n, "one", round(res["one_correction"]["bits_per_reading"], 3), "sd", res["one_correction"]["sd"])
    for r in res["top"][:8]: print("  ", r["terms"], round(r["bits_per_reading"], 3), {k: round(v, 3) for k, v in r["fit_coef"].items()}, {k: round(v, 3) for k, v in r["fit_alpha"].items()}, "sd", round(r["sd"], 3))
    for hh in res["holdout"]: print("  HO", hh["terms"], hh["transfers"], [(round(f["formula_bpr"], 3), round(f["one_correction_bpr"], 3)) for f in hh["folds"]])
    print("  neighbour", nb); print("  spike", {k: v for k, v in spike.items() if k != "single_terms"})
    for k, v in spike["single_terms"].items(): print("   ", k, round(v["bpr"], 3), {a: round(b, 3) for a, b in v["coef"].items()}, {a: round(b, 3) for a, b in v["alpha"].items()})
    (OUT / "phl_patch.json").write_text(json.dumps({k: v for k, v in res.items() if not k.startswith("_")}, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1])
