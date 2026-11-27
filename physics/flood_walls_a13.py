"""Protocol 15: exploratory compressor read of flood depth misses after the thin-film fix (A13 runs, tuning towns).

Usage: python physics/flood_walls_a13.py <CHI repo> [--spill lumber_barrier.json] --out results/physics/flood_walls_a13.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "physics"))
from physlib import guards as Q, terms as T  # noqa: E402
import flood_walls as FW  # noqa: E402

RUNS = (("lumber-sc-2016", "lumber-terms.csv", "lumber-sc-2016-edge-thinfilm-A13.json"),
        ("highland-falls-ny-2023", "highland-falls-terms-A9.csv", "highland-falls-ny-2023-thinfilm-A13.json"))


def rows_a13(chi: Path, spill: dict):
    out = []
    for site, ctx_csv, score in RUNS:
        ctx = {r["hwm_id"]: r for r in csv.DictReader(open(chi / "runs-summary/walls" / ctx_csv, encoding="utf-8"))}
        for m in json.load(open(chi / "docs/validation/walls-tuning" / score))["depth"]:
            r = dict(ctx[str(m["hwm_id"])])
            if m.get("err_m") is None:
                r.update(model_dry="1", err_m="")
            else:
                r.update(model_dry="0", err_m=f'{m["err_m"]:.3f}')
                if m["kind"] == "depth":
                    r["h_model_m"] = f'{m["model_m"]:.3f}'
            sp = spill.get(str(m["hwm_id"]))
            r["spill_rise_m"] = f'{sp["spill_from_channel"] - sp["ground"]:.3f}' if sp else ""
            out.append(r)
    return out


def explore(rows, name, B=30, N=19):
    ctx, y, sites, wb, dry = FW.build(rows)
    lum = np.where(sites == "lumber-sc-2016")[0]
    folds = [(f"leave out {i}", np.setdiff1d(np.arange(len(y)), [i]), np.array([i])) for i in lum]
    terms = [t for t in T.FLOOD_WALL if all(k in ctx for k in t.needs)]
    res = Q.analyse(terms, ctx, y, FW.GRID, folds=folds, groups=sites, bootstrap_B=B, null_N=N, name=name)
    gd = res["guards"]
    return {"name": name, "n": res["n"], "exploratory": True, "terms_offered": [t.name for t in terms],
            "winner": res["winner"]["terms"], "coef": res["winner"]["coef"], "bits_saved_per_reading": res["bits_saved_per_reading"],
            "sd_before": res["sd_before"], "sd_after": res["sd_after"], "null_p": gd.get("null_bits_saved", {}).get("p"),
            "inclusion": gd.get("inclusion"), "bootstrap": gd.get("bootstrap"), "verdict": Q.verdict(res),
            "dry": FW.dry_tally(dry), "_y": res["_y"].tolist(), "_mu": res["_mu"].tolist(),
            "_ids": [r["hwm_id"] for r in rows if r["model_dry"] == "0"]}


def synth(seed=1):
    rng = np.random.default_rng(seed); out = {}
    base = [{"site": "S", "hwm_id": str(i), "kind": "level", "quality": "Good", "type": "Seed line", "waterbody": "w",
             "at_building": str(int(rng.random() < 0.6)), "model_dry": "0", "h_model_m": f"{rng.uniform(0.2, 2.5):.3f}",
             "hand_m": f"{rng.uniform(0, 4):.3f}", "d_channel_m": f"{rng.uniform(5, 600):.1f}", "d_edge_m": f"{rng.uniform(500, 5000):.1f}",
             "bld_frac_30": f"{rng.uniform(0, 0.3):.4f}", "slope": f"{rng.uniform(0, 0.05):.4f}", "v_max_ms": "",
             "spill_rise_m": f"{rng.uniform(0, 0.4):.3f}"} for i in range(25)]
    for lab, law in (("storage_fill", lambda r: -0.8 * FW._f(r["d_channel_m"]) / (FW._f(r["d_channel_m"]) + 150)), ("null", lambda r: 0.0)):
        rr = [dict(r, err_m=f"{law(r) + rng.normal(0, 0.15):.3f}") for r in base]
        ctx, y, s, wb, _ = FW.build(rr)
        folds = [(f"leave out {i}", np.setdiff1d(np.arange(len(y)), [i]), np.array([i])) for i in range(len(y))]
        terms = [t for t in T.FLOOD_WALL if all(k in ctx for k in t.needs)]
        res = Q.analyse(terms, ctx, y, FW.GRID, folds=folds, bootstrap_B=10, null_N=19, name=f"synth {lab}")
        p = res["guards"]["null_bits_saved"]["p"]; w = res["winner"]["terms"]
        ok = (w == [] or p > 0.05) if lab == "null" else (lab in w)
        out[lab] = {"winner": w, "bits": res["bits_saved_per_reading"], "null_p": p, "pass": bool(ok)}
        print(lab, out[lab], flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("chi", type=Path)
    ap.add_argument("--spill", type=Path)
    ap.add_argument("--out", type=Path, default=ROOT / "results/physics/flood_walls_a13.json")
    a = ap.parse_args()
    spill = {o["hwm_id"]: o for o in json.load(open(a.spill))} if a.spill else {}
    res = {"synth": synth()}
    rows = rows_a13(a.chi, spill)
    res["all"] = explore(rows, "A13 both towns (exploratory, n<20)")
    keep = [r for r in rows if r["hwm_id"] not in ("44719", "44724")]
    res["no_step_pair"] = explore(keep, "A13 without HF 44719/44724 (exploratory)")
    res["lumber_only"] = explore([r for r in rows if r["site"] == "lumber-sc-2016"], "A13 Lumber only (exploratory)")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1, default=float))
    for k in ("all", "no_step_pair", "lumber_only"):
        r = res[k]
        print(k, json.dumps({x: r[x] for x in ("n", "winner", "coef", "bits_saved_per_reading", "sd_before", "sd_after", "null_p", "verdict")}, default=str))


if __name__ == "__main__":
    main()
