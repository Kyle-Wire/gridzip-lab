"""PROTOCOL-FIRE amendment A4: planted-mechanism checks for the fire model.

python physics/fire_planted.py        (needs a FireSim checkout: FIRESIM=/path, default ../FireSim)
Writes results/physics/fire_planted.json. Written fresh; reads FireSim as a library, changes nothing in it.
"""
import json, os, subprocess, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FS = Path(os.environ.get("FIRESIM", ROOT.parent / "FireSim"))
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(FS))
from physlib import engine as E                      # noqa: E402
from physlib.engine import Term                      # noqa: E402
from firesim import synth                            # noqa: E402
from firesim.landscape import fields                 # noqa: E402
from firesim.levelset import arrival                 # noqa: E402

D2L2 = (.06, .07, .08, .60, .90)
CELL, WIND, T, NPTS, CR = 10.0, 6.0, 90.0, 150, 0.4
GROUPS = {"grass": (102, 104), "shrub": (122, 145), "timber": (165, 183, 188)}


def one(seed, off=None):
    ls = synth.landscape(seed); g = synth.on_grid(ls, CELL)
    wto = np.random.default_rng(1000 + seed).uniform(0, 2 * np.pi)
    tan = np.zeros_like(g["tan"]) if off == "slope" else g["tan"]
    cc = np.zeros_like(g["cc"]) if off == "shelter" else g["cc"]
    F = fields(g["fuel"], D2L2, WIND, wto, tan, g["upslope"], cc, g["ch"], g["cbh"], g["cbd"], crown_ratio=CR, crown=off != "crown")
    b = ls["block"]; nb = ls["fuel"].shape[0]; by, bx = np.mgrid[:nb, :nb]
    tx = ls["size"] / 2 - 1200 * np.cos(wto); ty = ls["size"] / 2 - 1200 * np.sin(wto)
    d = np.hypot((bx + .5) * b - tx, (by + .5) * b - ty); d[ls["fuel"] == 91] = 1e9
    iy, ix = np.unravel_index(d.argmin(), d.shape)
    r, c = int((iy + .5) * b / CELL), int((ix + .5) * b / CELL)
    a, _ = arrival(F["rate"], F["heading"], F["lb"], CELL, [(r, c, 0, 15.0)], T)
    return a, g, (r, c), wto


def line_terms(g, ig, pts, wto):
    """Each vocabulary term averaged along the straight line from the ignition to each point."""
    out = {k: [] for k in ("slope_run", "slope_any", "crown_fuel", "canopy", "shelter", "grass", "shrub", "timber", "wind_run", "ln_dist")}
    for py, px in pts:
        n = int(max(abs(py - ig[0]), abs(px - ig[1]))) + 1
        yy = np.rint(np.linspace(ig[0], py, n)).astype(int); xx = np.rint(np.linspace(ig[1], px, n)).astype(int)
        ang = np.arctan2(py - ig[0], px - ig[1])
        t2 = g["tan"][yy, xx] ** 2
        out["slope_run"].append(np.mean(t2 * np.cos(g["upslope"][yy, xx] - ang)))
        out["slope_any"].append(np.mean(t2))
        can = g["cbd"][yy, xx] > 0
        out["crown_fuel"].append(np.mean(can & (g["cbh"][yy, xx] <= 3.0)))
        out["canopy"].append(np.mean(can))
        out["shelter"].append(np.mean(g["cc"][yy, xx] * CR / 3.0))
        f = g["fuel"][yy, xx]
        for k, nums in GROUPS.items():
            out[k].append(np.mean(np.isin(f, nums)))
        out["wind_run"].append(np.cos(wto - ang))
        out["ln_dist"].append(np.log(np.hypot(py - ig[0], px - ig[1]) * CELL))
    return {k: np.array(v) for k, v in out.items()}


TERMS = [
    Term("slope_run", "fire", "b · mean[tan²(slope) · cos(upslope − path)]", "fire runs faster upslope than the model says: slope is missing or weak", ("slope_run",), lambda c: c["slope_run"], group="slope"),
    Term("slope_any", "fire", "b · mean[tan²(slope)]", "steep ground is off whatever its direction", ("slope_any",), lambda c: c["slope_any"], group="slope"),
    Term("crown_fuel", "fire", "b · share of path under low canopy", "crown fire is missing or weak", ("crown_fuel",), lambda c: c["crown_fuel"]),
    Term("canopy", "fire", "b · share of path under any canopy", "something about forest canopy is off", ("canopy",), lambda c: c["canopy"]),
    Term("shelter", "fire", "b · mean[cover · crown ratio / 3]", "the wind under canopy is off (sheltering)", ("shelter",), lambda c: c["shelter"]),
    Term("grass", "fire", "b · share of path in grass", "grass spread is off", ("grass",), lambda c: c["grass"], group="fuel"),
    Term("shrub", "fire", "b · share of path in shrub", "shrub spread is off", ("shrub",), lambda c: c["shrub"], group="fuel"),
    Term("timber", "fire", "b · share of path in timber", "timber spread is off", ("timber",), lambda c: c["timber"], group="fuel"),
    Term("wind_run", "fire", "b · cos(wind − path)", "head against flank and back is off (fire shape or wind response)", ("wind_run",), lambda c: c["wind_run"]),
    Term("ln_dist", "fire", "b · ln(distance)", "the miss grows or shrinks with distance (start-up, acceleration)", ("ln_dist",), lambda c: c["ln_dist"]),
]
PLANT = {"slope": ("slope_run", "slope_any"), "crown": ("crown_fuel",), "shelter": ("shelter",), "twin": ()}
OTHER = {"slope": ("crown_fuel", "canopy", "shelter"), "crown": ("slope_run", "slope_any", "shelter"), "shelter": ("slope_run", "slope_any", "crown_fuel"), "twin": ()}


def main():
    seeds = list(range(20, 32)); truth = {s: one(s) for s in seeds}
    res = {"firesim_commit": subprocess.run(["git", "-C", str(FS), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
           "settings": dict(cell=CELL, wind=WIND, minutes=T, points=NPTS, noise_sd=0.05, seeds=seeds)}
    g01 = E.Grid(0.01, -4.0, 4.0)
    for off in ("slope", "crown", "shelter", "twin"):
        ys, ctxs, sd_id = [], [], []
        for s in seeds:
            a, g, ig, wto = truth[s]
            c = a if off == "twin" else one(s, off)[0]
            ok = np.isfinite(a) & np.isfinite(c) & (a > 5) & (c > 5)
            py, px = np.nonzero(ok)
            rng = np.random.default_rng(5000 + s)
            pick = rng.choice(len(py), min(NPTS, len(py)), replace=False)
            pts = list(zip(py[pick], px[pick]))
            y = np.log(a[py[pick], px[pick]] / c[py[pick], px[pick]]) + rng.normal(0, 0.05, len(pick))
            ys.append(y); ctxs.append(line_terms(g, ig, pts, wto)); sd_id.append(np.full(len(pick), s))
        y = np.concatenate(ys); sid = np.concatenate(sd_id)
        ctx = {k: np.concatenate([c[k] for c in ctxs]) for k in ctxs[0]}
        A = np.where(sid <= 25)[0]; B = np.where(sid > 25)[0]
        r = E.run_case(TERMS, ctx, y, g01, [("fit 20-25, test 26-31", A, B), ("fit 26-31, test 20-25", B, A)], max_terms=4, name=f"fire planted: {off} off")
        win = r["winner"]["terms"]
        if off == "twin":
            ok = len(win) == 0
        else:
            ok = any(t in win for t in PLANT[off]) and bool(r["transfers"]) and not any(t in win for t in OTHER[off])
        out = E.public(r); out.update(passed=bool(ok), n=int(len(y)), sd_y=float(y.std()), planted=list(PLANT[off]))
        res[off] = out
        print(off, "PASS" if ok else "FAIL", "winner", win, "bits saved/reading", round(r["bits_saved_per_reading"], 3), "transfers", r["transfers"], "sd", round(float(y.std()), 3), "n", len(y), flush=True)
    out_dir = ROOT / "results" / "physics"; out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "fire_planted.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
