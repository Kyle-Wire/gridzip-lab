"""PROTOCOL-FIRE amendment A14: the compressor on the errors of FireSim runs of real fires.

python physics/fire_real.py <name>=<table.npz> <name>=<table.npz> [--out results/physics/fire_real.json]

Tables come from FireSim scripts/error_table.py (one per fire). Three error tables are searched with the fire
formulas in physlib.terms.FIRE: T1 arrival, T2 reach, T3 house blocks. Each fire's own mean error is taken out first
and reported as that fire's level. Planted checks run on the real term columns before the real errors are read.
Written fresh.
"""
import json, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from physlib import engine as E, guards as Q, terms as T            # noqa: E402

CONST = ("ember_L", "gap_d", "slow_v", "held_v")
SLOW_V = (0.5, 1.0, 2.0, 5.0, 10.0)
A18 = "--a18" in sys.argv                                           # amendment A18: slow_edge on T2o only
THIN = [1]                                                          # amendment A21: lattice thinning step (set in main)
A16 = "--a16" in sys.argv                                           # amendment A16: tables T2 and T2o with the overrun inputs
PLANT = {"wind": ("lnu", None), "lee": ("lee5", None), "calm": ("calm", None), "ember": ("ember", 1000.0), "houses": ("gap", 10.0), "flank": ("flank", None),
         "late": ("late", None), "slow": ("slow", 2.0), "held": ("held", 2.44)}
GRIDS = {"T1": E.Grid(0.05, -8.0, 8.0), "T2": E.Grid(0.25, -3.0, 3.0), "T3": E.Grid(0.05, -2.0, 2.0), "T2o": E.Grid(0.25, -3.0, 3.0)}


def load(path):
    d = np.load(path, allow_pickle=False); meta = json.loads(str(d["meta"]))
    return {k: d[k] for k in d.files if k != "meta"}, meta


def ros(d, te, step_m=300.0):
    """A18: the model's spread rate (m/min) at each lattice point: the step over the largest arrival-time difference to
    a 4-neighbour the fire reached earlier. Infinite where there is no such neighbour (a start, or a spot fire)."""
    R = d["t12_row"] // 10; C = d["t12_col"] // 10; t = np.where(d["t12_t_mod"] <= te, d["t12_t_mod"], np.inf)
    G = np.full((R.max() + 3, C.max() + 3), np.inf); G[R + 1, C + 1] = t
    dt = np.zeros(len(t))
    for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        nb = G[R + 1 + dr, C + 1 + dc]; x = np.where(np.isfinite(nb) & (nb < t), t - nb, 0.0); dt = np.maximum(dt, np.where(np.isfinite(x), x, 0.0))
    return np.where(dt > 0, step_m / np.maximum(dt, 1e-9), np.inf)


def table(d, meta, which):
    """ctx and y of one fire for one table. Hourly inputs: T1 over the hours before the observed arrival, T2 and T3
    over the whole window (a point only the model burned has no observed time)."""
    te = meta["t_end_min"]
    if which == "T3":
        p, suf = "t3_", "_all"; rows = np.ones(len(d["t3_n"]), bool); y = d["t3_obs"] - d["t3_mod"]
        blk = d["t3_block3km"]
    else:
        p = "t12_"; ob = d["t12_t_obs"] <= te; mo = d["t12_t_mod"] <= te
        if which == "T1":
            suf = "_to"; rows = ob & mo
            with np.errstate(all="ignore"): y = np.log(d["t12_t_obs"] / np.maximum(d["t12_t_mod"], 10.0))
        elif which == "T2o":                                         # A16: points the model burned; 1 = the record does not show them burned
            suf = "_all"; rows = mo; y = 1.0 - ob.astype(float)
        else:
            suf = "_all"; rows = ob | mo; y = mo.astype(float) - ob.astype(float)
        blk = (d["t12_row"] // 100) * 1000 + d["t12_col"] // 100
    if which != "T3" and THIN[0] > 1:                            # A21: every THIN-th lattice point in both directions
        rows = rows & ((d["t12_row"] // 10) % THIN[0] == 0) & ((d["t12_col"] // 10) % THIN[0] == 0)
    g = lambda k: d[p + k][rows] if p + k in d else None
    c = {"lnu": np.log(np.maximum(g("u10" + suf), 0.3)), "run": g("run" + suf), "calm": g("calm" + suf), "night": g("night" + suf),
         "crew": g("crew" + suf), "moist": g("moist" + suf), "rh": g("rh" + suf), "lndist": np.log(g("dist")),
         "ember": g("ember"), "ember_L": np.array(meta["ember_L"]), "bld100": g("bld100")}
    if g("u850" + suf) is not None:
        c.update(aloft=np.log(np.maximum(g("u850" + suf), 0.3)) - c["lnu"], lee1=g("lee1" + suf) / 100.0, lee5=g("lee5" + suf) / 100.0)
    if g("dz_model") is not None: c["dz_model"] = g("dz_model") / 1000.0
    rd = g("road_dist")
    if rd is not None and np.all(np.isfinite(rd)): c["lnroad"] = np.log(np.maximum(rd, 15.0))
    for k in ("grass", "shrub", "timber", "nonburn", "developed", "canopy", "slope_run", "slope_any"):
        c["line_" + k] = g("line_" + k)
    if which == "T1": c["lnhours"] = np.log(d["t12_t_obs"][rows] / 60.0)
    if which in ("T2", "T2o") and A16:
        c.update(flank=np.sqrt(np.maximum(0.0, 1.0 - c["run"] ** 2)), back=(c["run"] < 0).astype(float))
    if which == "T2o": c.update(late=d["t12_t_mod"][rows] / te, edge_age=(te - d["t12_t_mod"][rows]) / 60.0)
    if which == "T2o" and A18:
        r_ = ros(d, te)[rows]; c.update(slow=np.stack([(r_ < v).astype(float) for v in SLOW_V], 1), slow_v=np.array(SLOW_V))
    if which == "T2o" and "t12_flame_arr" in d:                  # A22
        f_ = d["t12_flame_arr"][rows]; c.update(held=np.stack([(f_ < v).astype(float) for v in (1.22, 2.44)], 1), held_v=np.array((1.22, 2.44)))
    if which == "T3": c.update(gap=g("gap"), gap_d=np.array(meta["gap_d"]), mod=g("mod"))
    return c, y[rows], blk[rows]


def pool(fires, which):
    """Stack the fires: only inputs every fire has; each fire's mean error removed and returned as its level."""
    parts = [table(d, m, which) for d, m in fires.values()]
    keys = set.intersection(*[{k for k, v in c.items() if v is not None} for c, _, _ in parts])
    ctx = {k: (parts[0][0][k] if k in CONST else np.concatenate([c[k] for c, _, _ in parts])) for k in sorted(keys)}
    fid = np.concatenate([np.full(len(y), i) for i, (_, y, _) in enumerate(parts)])
    y = np.concatenate([y for _, y, _ in parts]); blk = np.concatenate([b + 10 ** 7 * i for i, (_, _, b) in enumerate(parts)])
    level = {n: float(y[fid == i].mean()) for i, n in enumerate(fires)}; spread = {n: float(y[fid == i].std()) for i, n in enumerate(fires)}
    for i in range(len(parts)): y[fid == i] -= y[fid == i].mean()
    return ctx, y, fid, blk, level, spread


def block_shuffle(y, blk, rng):
    """Noise with the errors' own 3 km structure but no tie to place: block means dealt to other blocks, and the
    within-block remainders shuffled inside each block."""
    ids, inv = np.unique(blk, return_inverse=True)
    mean = np.bincount(inv, y) / np.bincount(inv); out = mean[rng.permutation(len(ids))][inv]
    rem = y - mean[inv]
    for j in range(len(ids)):
        ii = np.where(inv == j)[0]; out[ii] += rem[rng.permutation(ii)]
    return out


def groups_of(terms): return sorted({next(t for t in T.FIRE if t.name == n).group or n for n in terms})


def planted(ctx, y, blk, g, seed=0):
    """Each planted term's group must be named; a noise-only twin must name nothing."""
    rng = np.random.default_rng(seed); n = len(y); built = E.assemble(T.FIRE, ctx, n); L = E.library_info(built)
    by = {b.name: b for b in built}; out = {}
    names = {"wind": "wind_speed", "lee": "lee_descent_5km", "calm": "calm_share", "ember": "ember", "houses": "gap", "flank": "flank_spread",
             "late": "late_arrival", "slow": "slow_edge", "held": "held_edge"}
    noise = block_shuffle(y, blk, rng)
    w, _, _ = E.search(built, noise, g, library_size=L, robust=True)
    out["noise_only"] = {"named": w["terms"], "pass": not w["terms"]}
    for grp, name in names.items():
        if name not in by: continue
        b = by[name]; col = (b.fn(PLANT[grp][1]) if b.fn else b.cols)[:, 0]; col = col - col.mean()
        if col.std() < 1e-12: continue
        yp = 0.5 * y.std() * col / col.std() + block_shuffle(y, blk, rng)
        w, _, _ = E.search(built, yp, g, library_size=L, robust=True)
        got = groups_of(w["terms"])
        out[grp] = {"planted": name, "named": w["terms"], "pass": grp in got, "extra_groups": [k for k in got if k != grp]}
    return out


def block_null(res, y, blk, g, N=19, seed=1):
    rng = np.random.default_rng(seed); built = res["_built"]; L = E.library_info(built); out = []
    for _ in range(N):
        yp = block_shuffle(y, blk, rng)
        w, ranked, _ = E.search(built, yp, g, library_size=L, robust=True)
        none = next((r for r in ranked if not r["terms"]), None) or E.exact_bits([], yp, g, L, robust=True)
        out.append((none["total_bits"] - w["total_bits"]) / len(y))
    return out


def one(fires, which):
    g = GRIDS[which]; ctx, y, fid, blk, level, spread = pool(fires, which); names = list(fires); t = time.time()
    out = {"n": int(len(y)), "n_by_fire": {n: int((fid == i).sum()) for i, n in enumerate(names)}, "level": level, "spread": spread,
           "inputs": sorted(k for k in ctx if k not in CONST)}
    if len(y) < 40:
        out["skipped"] = "fewer than 40 readings"; return out
    out["planted"] = planted(ctx, y, blk, g)
    folds = [(f"fit {names[1 - i]}, test {names[i]}" if len(names) == 2 else f"fit the others, test {names[i]}", np.where(fid != i)[0], np.where(fid == i)[0])
             for i in range(len(names))] if len(names) >= 2 else []
    res = Q.analyse(T.FIRE, ctx, y, g, grammar=False, folds=folds, groups=blk, bootstrap_B=30, null_N=0, name=which, blocks=None,
                    neighbour=False)
    nb = block_null(res, y, blk, g)
    gd = res["guards"]
    # neighbour test over all inputs at once (the library routine builds n x n tables for every pair of inputs,
    # which does not fit in memory at several thousand points)
    from scipy.spatial import cKDTree
    from physlib import grammar as G
    prims = G.primitives(ctx, len(y), ())
    X = np.column_stack([np.argsort(np.argsort(v)) / (len(v) - 1) for v in prims.values()])
    idx = cKDTree(X).query(X, k=9)[1][:, 1:]
    gd["neighbour"] = Q.neighbour_test(prims, y - res["_mu"], nb=[("all inputs", idx)])
    out.update(winner=res["winner"]["terms"], coef=res["winner"]["coef"], alpha=res["winner"]["alpha"],
               bits_saved_per_reading=res["bits_saved_per_reading"], sd_before=res["sd_before"], sd_after=res["sd_after"],
               holdout=res["holdout"], transfers=res["transfers"], ladder=res["ladder"], ranked=res["ranked"][:5],
               verdict=Q.verdict(res), bootstrap=gd.get("bootstrap"), lookalikes=gd.get("lookalikes"), neighbour=gd.get("neighbour"),
               null_bits_saved={"max": max(nb), "mean": float(np.mean(nb)), "all": nb}, terms_meta=res["terms_meta"])
    per = {}
    for i, n in enumerate(names):                                # each fire alone: a lead for that fire, drives nothing
        ii = np.where(fid == i)[0]
        if len(ii) >= 40:
            c1 = {k: (v if k in CONST else v[ii]) for k, v in ctx.items()}
            r1 = E.run_case(T.FIRE, c1, y[ii], g, robust=True, name=f"{which} {n}")
            per[n] = {"winner": r1["winner"]["terms"], "coef": r1["winner"]["coef"], "alpha": r1["winner"]["alpha"],
                      "bits_saved_per_reading": r1["bits_saved_per_reading"], "sd_before": r1["sd_before"], "sd_after": r1["sd_after"]}
    out["each_fire_alone"] = per
    # the four rules of A14, per winning term
    sel = (gd.get("bootstrap") or {}).get("selected", {}); find = {}
    for tn in out["winner"]:
        grp = groups_of([tn])[0]; meta = res["terms_meta"].get(tn, {})
        pl = out["planted"].get(grp)
        find[tn] = {"group": grp, "physics": meta.get("physics", True), "transfers_both_ways": bool(res["transfers"]),
                    "bootstrap_share": sel.get(grp, sel.get(tn)), "beats_every_null": out["bits_saved_per_reading"] > max(nb),
                    "planted_check": None if pl is None else pl["pass"]}
        f = find[tn]
        f["finding"] = bool(f["physics"] and f["transfers_both_ways"] and (f["bootstrap_share"] or 0) >= 0.5 and f["beats_every_null"]
                            and f["planted_check"] is not False and out["planted"]["noise_only"]["pass"])
    out["rules"] = find; out["seconds"] = round(time.time() - t, 1)
    print(which, "n", out["n"], "level", {k: round(v, 3) for k, v in level.items()}, "winner", out["winner"],
          "saved", round(out["bits_saved_per_reading"], 3), "null max", round(max(nb), 3), "transfers", res["transfers"], flush=True)
    return out


def twins(fires, which, seeds=5):
    """A21: how many of `seeds` noise-only twins name a term at the current thinning."""
    g = GRIDS[which]; ctx, y, fid, blk, _, _ = pool(fires, which); built = E.assemble(T.FIRE, ctx, len(y)); L = E.library_info(built); bad = []
    for sd in range(seeds):
        w, _, _ = E.search(built, block_shuffle(y, blk, np.random.default_rng(sd)), g, library_size=L, robust=True)
        bad.append(w["terms"])
    return len(y), bad


def main():
    args = [a for a in sys.argv[1:] if "=" in a and not a.startswith("--")]
    outp = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else ROOT / "results" / "physics" / "fire_real.json"
    fires = {a.split("=", 1)[0]: load(a.split("=", 1)[1]) for a in args}
    tabs = next((a.split("=")[1].split(",") for a in sys.argv if a.startswith("--tables=")), None)
    if "--thin=auto" in sys.argv:                                # A21: chosen on noise alone, before the real errors are read
        choice = {}
        for s_ in (1, 2, 3, 5, 10):
            THIN[0] = s_; out_ = {w_: twins(fires, w_) for w_ in tabs}
            choice[s_] = {w_: {"n": n_, "twins_naming": sum(bool(b) for b in bad), "named": bad} for w_, (n_, bad) in out_.items()}
            print("thin", s_, {w_: (v["n"], v["twins_naming"]) for w_, v in choice[s_].items()}, flush=True)
            if all(v["twins_naming"] == 0 for v in choice[s_].values()): break
        else:
            outp.parent.mkdir(parents=True, exist_ok=True); outp.write_text(json.dumps({"protocol": "PROTOCOL-FIRE A21", "thinning": choice, "result": "no thinning step silences the noise twins; no real search run"}, indent=1))
            print("no thinning step silences the twins; wrote", outp); return
    res = {"thin": THIN[0], "thinning": locals().get("choice"), "protocol": "PROTOCOL-FIRE A18" if A18 else "PROTOCOL-FIRE A16 part 2" if A16 else "PROTOCOL-FIRE A14", "tables": {n: m for n, (_, m) in fires.items()}}
    for which in (tuple(a.split("=")[1].split(",") for a in sys.argv if a.startswith("--tables="))[0] if any(a.startswith("--tables=") for a in sys.argv) else ("T2o",) if A18 else ("T2", "T2o") if A16 else ("T1", "T2", "T3")):
        res[which] = one(fires, which)
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    print("wrote", outp)


if __name__ == "__main__":
    main()
