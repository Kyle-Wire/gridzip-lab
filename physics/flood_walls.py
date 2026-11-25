"""Protocol 11: what FloodSim's misses at real high-water marks follow (flood wall check, tuning towns only).

Usage:
  python physics/flood_walls.py synth                      # planted-law check and null on made-up marks (CPU, ~1 min)
  python physics/flood_walls.py run terms1.csv [terms2.csv ...] --out results/physics/flood_walls.json

Each CSV is one run's per-mark context from Climate Home Inspection `tools/gate2/walls_terms.py` (one tuning town or a
spent case). The target is rule D's error, model minus mark (m). Model-dry marks are tallied apart, by term value.
Folds: leave one town out when there are two towns or more, else leave one waterbody out.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from physlib import engine as E, guards as Q, terms as T  # noqa: E402

GRID = E.Grid(0.01, -5.0, 5.0)


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def load(paths):
    rows = []
    for p in paths:
        with open(p, newline="", encoding="utf-8") as f:
            rows += list(csv.DictReader(f))
    return rows


def build(rows):
    """ctx, y, groups for the wet marks; the dry tally apart. Inputs with any blank are left out of ctx (terms skip)."""
    wet = [r for r in rows if r["model_dry"] in ("0", 0) and r["err_m"] not in ("", None)]
    dry = [r for r in rows if r["model_dry"] in ("1", 1)]
    y = np.array([_f(r["err_m"]) for r in wet])
    col = lambda k: np.array([_f(r.get(k)) for r in wet])  # noqa: E731
    ctx = {}
    for key, src in (("h_model", "h_model_m"), ("hand", "hand_m"), ("d_channel", "d_channel_m"), ("d_edge", "d_edge_m"),
                     ("bld_frac", "bld_frac_30"), ("slope", "slope"), ("V", "v_max_ms")):
        v = col(src)
        if len(v) and np.isfinite(v).all():
            ctx[key] = v
    ctx["kind_hag"] = np.array([float(r["kind"] == "depth") for r in wet])
    ctx["at_building"] = np.array([float(r["at_building"] in ("1", 1)) for r in wet])
    sites = np.array([r["site"] for r in wet])
    if len(set(sites.tolist())) > 1:
        ctx["site_onehot"] = T._onehot(sites)
    for k in ("kind_hag", "at_building"):
        if np.ptp(ctx[k]) == 0:
            del ctx[k]
    wb = np.array([f'{r["site"]}:{r["waterbody"]}' for r in wet])
    return ctx, y, sites, wb, dry


def folds_for(sites, wb):
    key = sites if len(set(sites.tolist())) > 1 else wb
    return [(f"leave out {k}", np.where(key != k)[0], np.where(key == k)[0]) for k in sorted(set(key.tolist()))
            if (key == k).sum() >= 2 and (key != k).sum() >= 4], key


def analyse(ctx, y, sites, wb, name, B=30, N=19):
    folds, groups = folds_for(sites, wb)
    terms = [t for t in T.FLOOD_WALL if all(k in ctx for k in t.needs)]
    res = Q.analyse(terms, ctx, y, GRID, folds=folds, groups=groups, bootstrap_B=B, null_N=N, name=name)
    gd = res["guards"]
    return {"name": name, "n": res["n"], "terms_offered": [t.name for t in terms], "winner": res["winner"]["terms"],
            "coef": res["winner"]["coef"], "bits_saved_per_reading": res["bits_saved_per_reading"],
            "sd_before": res["sd_before"], "sd_after": res["sd_after"], "holdout": res["holdout"],
            "transfers": res["transfers"], "verdict": Q.verdict(res),
            "null_p": gd.get("null_bits_saved", {}).get("p"), "inclusion": gd.get("inclusion"),
            "neighbour": gd.get("neighbour"), "_y": res["_y"], "_mu": res["_mu"]}


def dry_tally(dry):
    if not dry:
        return {"n": 0}
    out = {"n": len(dry)}
    for k in ("hand_m", "d_channel_m", "d_edge_m", "bld_frac_30", "slope"):
        v = np.array([_f(r.get(k)) for r in dry]); v = v[np.isfinite(v)]
        if len(v):
            out[k] = {"median": float(np.median(v)), "min": float(v.min()), "max": float(v.max())}
    return out


def synth(seed=0):
    """Three made-up towns, 25 marks each, noise sd 0.15 m. Planted: (a) edge loss −0.4·exp(−d_edge/120 m);
    (b) volume scale +0.25·h; (c) nothing. Pass: (a) and (b) named with no false term, (c) names nothing (null p > 0.05)."""
    rng = np.random.default_rng(seed)
    n = 75
    sites = np.repeat(["A", "B", "C"], 25)
    rows = []
    for i in range(n):
        rows.append({"site": sites[i], "hwm_id": str(i), "kind": "level" if rng.random() < 0.6 else "depth",
                     "quality": "Good", "type": "Seed line", "waterbody": f"w{i % 2}", "at_building": str(int(rng.random() < 0.6)),
                     "model_dry": "0", "h_model_m": f"{rng.uniform(0.1, 2.5):.3f}", "hand_m": f"{rng.uniform(0, 6):.3f}",
                     "d_channel_m": f"{rng.uniform(1, 300):.1f}", "d_edge_m": f"{rng.uniform(5, 800):.1f}",
                     "bld_frac_30": f"{rng.uniform(0, 0.4):.4f}", "slope": f"{rng.uniform(0, 0.08):.4f}", "v_max_ms": ""})
    out = {}
    for lab, law in (("edge_loss", lambda r: -0.4 * np.exp(-_f(r["d_edge_m"]) / 120)),
                     ("volume_scale", lambda r: 0.25 * _f(r["h_model_m"])), ("null", lambda r: 0.0)):
        rr = [dict(r, err_m=f"{law(r) + rng.normal(0, 0.15):.3f}") for r in rows]
        ctx, y, s, wb, _ = build(rr)
        a = analyse(ctx, y, s, wb, f"synth {lab}", B=10, N=19)
        a.pop("_y"); a.pop("_mu")
        ok = (a["winner"] == [] or a["null_p"] > 0.05) if lab == "null" else (lab in a["winner"] and all(
            t == lab or t == "town_offset" for t in a["winner"]))
        out[lab] = {"winner": a["winner"], "bits": a["bits_saved_per_reading"], "null_p": a["null_p"], "pass": bool(ok)}
        print(lab, out[lab], flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("synth")
    r = sub.add_parser("run")
    r.add_argument("csv", nargs="+", type=Path)
    r.add_argument("--out", type=Path, default=ROOT / "results" / "physics" / "flood_walls.json")
    a = ap.parse_args()
    if a.cmd == "synth":
        res = synth()
        p = ROOT / "results" / "physics" / "flood_walls_synth.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(res, indent=1, default=float))
        sys.exit(0 if all(v["pass"] for v in res.values()) else 1)
    rows = load(a.csv)
    ctx, y, s, wb, dry = build(rows)
    res = analyse(ctx, y, s, wb, "flood walls: " + ", ".join(sorted(set(s.tolist()))))
    arr = {"y": res.pop("_y"), "mu": res.pop("_mu")}
    res["dry"] = dry_tally(dry)
    res["sources"] = [p.name for p in a.csv]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1, default=float))
    np.savez(a.out.with_suffix(".npz"), **arr)
    print(json.dumps({k: res[k] for k in ("n", "winner", "bits_saved_per_reading", "sd_before", "sd_after", "null_p", "verdict")}, indent=1, default=str))


if __name__ == "__main__":
    main()
