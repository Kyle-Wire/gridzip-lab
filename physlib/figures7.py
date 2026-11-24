"""PROTOCOL-7 figures. python -m physlib.figures7  -> physics/figures/library/19-22_*.png"""
import json, warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import engine as E, synth as S, terms as T, grammar as G, guards as Q
from .validate7 import plant_snr

ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results" / "physlib"; FIG = ROOT / "physics" / "figures" / "library"
SURF, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
GOOD, CRIT = "#0ca30c", "#d03b3b"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "axes.titlecolor": INK, "axes.titleweight": "bold", "axes.titlesize": 11, "axes.titlelocation": "left"})


def load(s):
    p = RES / f"v7_{s}.json"; return json.loads(p.read_text()) if p.exists() else None


def rules():
    """(group, label, observed, threshold, direction, pass)."""
    R = []
    A = load("A")
    if A:
        for k, v in A["summary"].items(): R.append(("A  bigger library", k[3:], min(v["named"], v["explained"]), 0.8, ">=", v["pass"]))
    Ar = load("Areg")
    if Ar:
        worst = min(Ar["summary"].items(), key=lambda kv: kv[1]["named"] - (kv[1]["protocol4"] or 0))
        allp = all(v["pass"] for v in Ar["summary"].values())
        d = worst[1]["named"] - (worst[1]["protocol4"] or 0)
        R.append(("A  bigger library", f"old 12 laws, worst change ({worst[0]})", d, -0.10, ">=", allp))
    B = load("B")
    if B:
        for k, v in B["summary"].items():
            obs = v.get("shape_found", v.get("explained_0.8", v.get("named"))); th = 0.7 if "shape_found" in v or "explained_0.8" in v else 0.8
            R.append(("B  new term shapes", k[3:], obs, th, ">=", v["pass"]))
    C1 = load("C1")
    if C1:
        for k, v in C1["summary"].items(): R.append(("C  guards: noise", f"noise {k}", v["any_term"], 0.05, "<=", v["pass"]))
    C2 = load("C2")
    if C2:
        s = C2["summary"]; R.append(("C  guards", "twins r 0.998: says it can't tell", s["twins r=0.998"]["says_cant_tell"], 0.8, ">=", s["twins r=0.998"]["pass"]))
        R.append(("C  guards", "pair r 0.6: names A, sure", s["pair r=0.6"]["a_named_sure"], 0.8, ">=", s["pair r=0.6"]["pass"]))
    C3 = load("C3")
    if C3:
        v = C3["summary"]["group confounding"]
        R.append(("C  guards", f"station-level term flagged when named (named in {v['attr_named_runs']} of 20)", v["flagged_when_named"], 0.9, ">=", v["flagged_when_named"] >= 0.9))
        R.append(("C  guards", "unflagged false physics", v["unflagged_false_share"], 0.05, "<=", v["unflagged_false_share"] <= 0.05))
    C4 = load("C4")
    if C4:
        s = C4["summary"]
        for n in (30, 60, 120): R.append(("C  guards", f"false terms per run, n = {n}", s[f"n={n}"]["false_per_run"], 0.2, "<=", s[f"n={n}"]["pass_false"]))
        if s["n=120"].get("coverage") == s["n=120"].get("coverage"):
            R.append(("C  guards", "bootstrap 90 % interval covers truth (n 120)", s["n=120"]["coverage"], 0.8, ">=", bool(s["n=120"].get("pass_coverage"))))
    C5 = load("C5")
    if C5:
        s = C5["summary"]
        R.append(("C  guards", "no false 'structure left' after the law", s["a law named, no structure left"]["flagged"], 0.10, "<=", s["a law named, no structure left"]["pass"]))
        R.append(("C  guards", "flags a shape outside the grammar", s["b form outside grammar"]["flagged"], 0.7, ">=", s["b form outside grammar"]["pass"]))
        R.append(("C  guards", "noise level read within 10 %", s["c noise: sd within 10 %"]["share"], 0.9, ">=", s["c noise: sd within 10 %"]["pass"]))
    return R


def fig_scoreboard():
    R = rules(); n = len(R)
    pos, cur, prev = [], 0.0, None
    for r in R:
        if r[0] != prev: cur += 1.0; prev = r[0]
        pos.append(cur); cur += 1.0
    tot = cur; y = tot - np.array(pos); n = tot
    fig, ax = plt.subplots(figsize=(10.5, 0.3 * n + 1.6))
    for yi, (grp, lab, obs, th, dirn, ok) in zip(y, R):
        if lab.startswith("old 12"):
            ax.text(0.5, yi, f"{'✓ pass' if ok else '✗ fail'}  every law within 10 points of protocol 4 (worst {obs:+.2f})", va="center", ha="center",
                    color=INK2, fontsize=9); ax.text(-0.02, yi, lab, ha="right", va="center", fontsize=9, color=INK); continue
        ax.plot([0, 1], [yi, yi], color=GRID, lw=6, solid_capstyle="round", zorder=1)
        ax.plot([th, th], [yi - 0.3, yi + 0.3], color=INK2, lw=1.5, zorder=2)
        ax.scatter([obs], [yi], s=70, color=GOOD if ok else CRIT, zorder=3, edgecolor=SURF, linewidth=2)
        ax.text(1.03, yi, ("✓ pass" if ok else "✗ fail") + f"  {obs:.0%}" + (f"  (rule {dirn} {th:.0%})"), va="center", fontsize=9, color=INK)
        ax.text(-0.02, yi, lab, ha="right", va="center", fontsize=9, color=INK)
    prev = None
    for yi, r in zip(y, R):
        if r[0] != prev: ax.text(-0.02, yi + 0.75, r[0], ha="right", va="bottom", fontsize=9.5, fontweight="bold", color=INK); prev = r[0]
    ax.set_xlim(0, 1); ax.set_ylim(-0.3, n + 0.3); ax.set_yticks([]); ax.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax.set_xticklabels(["0", "25 %", "50 %", "75 %", "100 %"])
    ax.spines["left"].set_visible(False)
    npass = sum(r[5] for r in R)
    ax.set_title(f"The stronger compressor against its own rules: {npass} of {len(R)} pass\n(dot = what happened, tick = the rule, written before any run)", pad=14, loc="center")
    fig.text(0.01, 0.005, "Synthetic cases with a known answer, 20 seeds each (noise sets 30-50). Rules were written before any run (PROTOCOL-7.md, A1).",
             fontsize=8, color=MUTED)
    fig.subplots_adjust(left=0.37, right=0.79, top=0.93, bottom=0.05); fig.savefig(FIG / "19_v7_scoreboard.png", dpi=150); plt.close(fig)


def fig_discovery():
    warnings.filterwarnings("ignore")
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    cases = [("Contrast that changes with neighbour height", lambda c: c["hc"] * c["lp"], "hc · lp"),
             ("Corner jets fading with distance", lambda c: c["corner"] * np.exp(-c["fetch"] / np.quantile(c["fetch"][c["fetch"] > 0], 0.5)), "corner · e^(−fetch/L)"),
             ("A threshold in height", lambda c: np.maximum(0, c["z"] - np.quantile(c["z"], 0.5)), "max(0, z − z̄)")]
    for ax, (title, f, lab) in zip(axs, cases):
        rng = np.random.default_rng(0); ctx = S.wind_ctx(rng); n = len(ctx["lp"]); col = f(ctx)
        sig, _ = plant_snr([col - col.mean()], [1.0], 0.15); y = 0.1 + sig + rng.normal(0, 0.15, n)
        r = E.run_case(list(T.WIND) + G.candidates(ctx, n), ctx, y, E.Grid(0.01, -4, 4), robust=True)
        fit = r["_mu"] - r["_mu"].mean()
        ax.scatter(sig, y - y.mean(), s=9, color=MUTED, alpha=0.5, label="readings (law + noise)", linewidths=0)
        o = np.argsort(sig); ax.plot(sig[o], fit[o], color=BLUE, lw=2, label="what the compressor found")
        ax.plot([sig.min(), sig.max()], [sig.min(), sig.max()], color=ORANGE, lw=1.5, ls="--", label="the planted law")
        ax.set_title(title, fontsize=10); ax.set_xlabel(f"planted: {lab}  (not in the library)"); ax.set_ylabel("model error (log)")
        ax.text(0.02, 0.97, "found: " + " + ".join(t.replace("g:", "") for t in r["winner"]["terms"]), transform=ax.transAxes, va="top", fontsize=8.5, color=INK,
                wrap=True)
    axs[0].legend(loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("Finding term shapes the library doesn't have: the grammar builds ~2,000 candidates and pays for each one's description",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout(); fig.savefig(FIG / "20_v7_discovery.png", dpi=150); plt.close(fig)


def fig_guards():
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
    C2 = load("C2"); ax = axs[0]
    if C2:
        for i, rho in enumerate((0.6, 0.97, 0.998)):
            rr = [r for r in C2["runs"] if r["rho"] == rho]; v = [r["inclusion_a"] for r in rr]
            ax.scatter(np.full(len(v), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(v)), v, s=22, color=BLUE, edgecolor=SURF, linewidth=1)
            ct = np.mean([r["cant_tell"] for r in rr])
            ax.text(i, 1.08, f"says can't tell: {ct:.0%}", ha="center", fontsize=8.5, color=INK2)
        ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["r 0.6\n(different)", "r 0.97\n(close)", "r 0.998\n(twins)"])
        ax.axhline(0.9, color=MUTED, lw=1, ls=":"); ax.set_ylim(0, 1.15); ax.set_ylabel("odds it's the planted input (A)")
        ax.set_title("Look-alike terms: sure only\nwhen the data can tell", fontsize=10)
    C4 = load("C4"); ax = axs[1]
    if C4:
        ns = [30, 60, 120]; s = C4["summary"]
        ax.plot(ns, [s[f"n={n}"]["named"] for n in ns], "-o", color=BLUE, lw=2, ms=8, label="planted cold-pool term found")
        ax.plot(ns, [s[f"n={n}"]["false_per_run"] for n in ns], "-o", color=ORANGE, lw=2, ms=8, label="false terms per run")
        ax.axhline(0.2, color=MUTED, lw=1, ls=":"); ax.text(32, 0.23, "false-term limit 0.2", fontsize=8, color=MUTED)
        cv = s["n=120"].get("coverage")
        if cv == cv and cv is not None: ax.text(120, s["n=120"]["named"] - 0.12, f"90 % interval\ncovers truth {cv:.0%}", ha="right", fontsize=8.5, color=INK2)
        ax.set_xticks(ns); ax.set_xlabel("readings"); ax.set_ylim(0, 1.1); ax.legend(fontsize=8, frameon=False, loc="center right")
        ax.set_title("Small samples: quiet\nrather than guessing", fontsize=10)
    C1 = load("C1"); ax = axs[2]
    if C1:
        items = sorted(C1["summary"].items()); y = np.arange(len(items))
        ax.barh(y, [v["any_term"] for _, v in items], color=[GOOD if v["pass"] else CRIT for _, v in items], height=0.6)
        ax.set_yticks(y); ax.set_yticklabels([k.replace("t3", "heavy tails").replace("ar1", "correlated") for k, _ in items], fontsize=8.5)
        for yi, (_, v) in zip(y, items): ax.text(0.002, yi, f"{round(v['any_term'] * v['runs'])} of {v['runs']} runs", va="center", fontsize=8, color=INK2)
        ax.axvline(0.05, color=INK2, lw=1.2); ax.text(0.052, len(items) - 0.4, "5 % limit", fontsize=8, color=INK2)
        ax.set_xlim(0, max(0.1, max(v["any_term"] for _, v in items) * 1.3)); ax.set_xlabel("share of pure-noise runs where it named anything")
        ax.set_title("Pure noise, grammar on:\nfalse finds", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "21_v7_guards.png", dpi=150); plt.close(fig)


def fig_real():
    cases = []
    for f, keys in (("D_w1", ("lbm", "rans")), ("D_w2", ("thunderstorm", "ordinary")), ("D_storm", ("storm",))):
        d = load(f)
        if d:
            for k in keys:
                if k in d: cases.append((f, k, d[k]))
    if not cases: return
    lab = {"lbm": "AIJ Case C, LBM", "rans": "AIJ Case C, RANS", "thunderstorm": "Layer 1, thunderstorm gusts", "ordinary": "Layer 1, ordinary gusts",
           "storm": "Thunderstorm table with cold pool\n(post hoc, spent data)"}
    fig, ax = plt.subplots(figsize=(12, 1.25 * len(cases) + 1.2)); y = np.arange(len(cases))[::-1]
    for yi, (_, k, d) in zip(y, cases):
        nb = d.get("null_bits_saved", {}).get("all", [])
        ax.scatter(nb, np.full(len(nb), yi), s=16, color=MUTED, alpha=0.7, linewidths=0)
        ax.scatter([d["bits_saved_per_reading"]], [yi], s=90, color=BLUE, edgecolor=SURF, linewidth=2, zorder=3)
        ax.text(-0.02, yi, lab[k], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=9.5, color=INK)
        sel = d["guards"].get("bootstrap", {}).get("selected", {})
        grp = {t.name: (t.group or t.name) for t in T.WIND_FULL}
        pick = lambda t: sel.get(grp.get(t, t), 0.0)
        words = "  ·  ".join(f"{t.replace('g:', '').replace('−-', '+')} ({pick(t):.0%})" for t, v in d["verdict"])
        nbd = d["guards"].get("neighbour", {})
        left = " · structure still unnamed" if nbd.get("structure_left") else ""
        ax.text(0.99, yi + 0.28, (words or "nothing pays") + left, transform=ax.get_yaxis_transform(), ha="right", fontsize=8, color=INK2)
    ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xlabel("bits saved per reading (blue: the real errors; grey: the same search on shuffled errors, 19 times)\nText: what it named, with the share of 30 resamples that pick the same term")
    ax.set_ylim(-0.7, len(cases) - 0.3)
    ax.set_title("The real cases again, grammar and guards on: is what it finds bigger than what it finds in noise?", pad=10)
    fig.subplots_adjust(left=0.24, right=0.98, top=0.9, bottom=0.14); fig.savefig(FIG / "22_v7_real_cases.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    for f in (fig_scoreboard, fig_guards, fig_real, fig_discovery): f(); print("ok", f.__name__)
