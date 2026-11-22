"""Storytelling figures for the physics compressor (PROTOCOL-4): python -m physlib.figures [out dir]."""
import json, math, sys
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results" / "physlib"
FIG = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "physics" / "figures" / "library"; FIG.mkdir(parents=True, exist_ok=True)
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"
BEFORE, AFTER = ORANGE, BLUE
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "savefig.facecolor": SURF, "axes.titlesize": 10.5, "axes.titlelocation": "left", "legend.frameon": False})
NICE = {"contrast": "speed-up contrast", "fetch_decay": "entrance decay", "fetch_decay_along": "entrance decay (along streets)",
        "height_ratio": "taller neighbour", "canyon_cross": "cross streets", "canyon_along": "along streets", "contrast_along": "contrast (along streets)",
        "station_bias": "station's own offset", "amplitude": "daily swing too big", "timing_shift": "cycle runs late",
        "sky_cooling": "clear-night cooling", "station_offsets": "station siting", "radiation_lag": "sunlight timing",
        "level": "one correction", "one correction": "one correction", "onset": "onset (motion starts)", "jump_a": "reference jump",
        "jump_b": "second jump", "seasonal": "seasonal cycle", "late_rate": "late rate change", "rate_change": "rate change",
        "shrink": "shrink to area mean", "near_house": "near the house", "distance": "distance to house", "grade": "street grade",
        "splash_over": "splash-over", "weir_exponent": "depth exponent", "orifice_transition": "weir→orifice", "datum": "depth datum",
        "consolidation": "consolidation", "linear": "rate offset", "acceleration": "acceleration", "gust_factor": "gust factor / roughness",
        "log_law_z0": "roughness length", "direction": "wind direction", "diurnal_shape": "unnamed daily shape", "noise_level": "noise level",
        "seasonal_amp": "seasonal amplitude", "flow": "flow size", "velocity_head": "velocity head", "frontal_share": "frontal share"}
nice = lambda t: NICE.get(t, t.replace("_", " "))


def load(name):
    p = RES / f"{name}.json"; a = RES / f"{name}_arrays.npz"
    if not p.exists(): return None, None
    j = json.loads(p.read_text()); arr = {}
    if a.exists():
        z = np.load(a, allow_pickle=True)
        for k in z.files:
            c, kk = k.split("__", 1); arr.setdefault(c, {})[kk] = z[k]
    return j, arr


def ladder_ax(ax, res, title=None, color=AFTER):
    lad = res["ladder"]; n = res["n"]
    b = [s["bits_per_reading"] for s in lad]; labels = [nice(s["added"]) for s in lad]
    x = np.arange(len(b))
    ax.bar(x[:1], b[:1], color=MUTED, width=0.62)
    for i in range(1, len(b)):
        ax.bar(i, b[i - 1] - b[i], bottom=b[i], color=color, width=0.62)
        ax.text(i, b[i - 1] + 0.01 * max(b), f"−{(b[i-1]-b[i]):.2f}", ha="center", va="bottom", fontsize=8, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(["one\ncorrection"] + ["+ " + l for l in labels[1:]], fontsize=7.5, rotation=0)
    lo = min(b) - 0.25 * (max(b) - min(b) + 0.05); ax.set_ylim(max(0, lo), max(b) + 0.12 * (max(b) - min(b) + 0.05))
    ax.set_ylabel("bits per reading"); ax.grid(axis="y", color=GRID, lw=0.5)
    ax.axhline(b[-1], color=INK2, lw=0.8, ls=":")
    if title: ax.set_title(title)


def hour_mean(h, v):
    hh = np.round(h).astype(int) % 24; return np.array([np.mean(v[hh == k]) if np.any(hh == k) else np.nan for k in range(24)])


# ------------------------------------------------------------------------------------------------ 11 how it works
def fig_how():
    from physlib import synth as S, engine as E
    terms, ctx, y, g, planted = S.case("ground", "sinkhole onset 2022", 4)
    res = E.run_case(terms, ctx, y, E.Grid(*g), prev=S.prev_of("ground", len(y)))
    t = ctx["t"]
    fig = plt.figure(figsize=(13, 4.3)); gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.25], wspace=0.32)
    ax = fig.add_subplot(gs[0]); ax.plot(t, y - y.mean(), ".", ms=3, color=MUTED)
    ax.set_title("1. A model's errors (here: a made-up InSAR series)"); ax.set_xlabel("year"); ax.set_ylabel("measured − model (mm)")
    ax.grid(color=GRID, lw=0.5)
    ax = fig.add_subplot(gs[1]); ax.axis("off"); ax.set_title("2. Try terms from the physics library")
    cards = [("onset  v·max(0, t − t₀)", AFTER, "kept: saves bits"), ("reference jump  s·H(t − tⱼ)", MUTED, "dropped: costs more than it saves"),
             ("seasonal  a·sin 2πt + b·cos 2πt", MUTED, "dropped"), ("consolidation  a·e^(−t/τ)", MUTED, "dropped"), ("… 5 more ground terms", MUTED, "")]
    for i, (txt, col, note) in enumerate(cards):
        yy = 0.86 - i * 0.18
        ax.add_patch(FancyBboxPatch((0.02, yy - 0.06), 0.96, 0.12, boxstyle="round,pad=0.01,rounding_size=0.02", fc=SURF, ec=col, lw=1.6 if col == AFTER else 0.8))
        ax.text(0.05, yy + 0.005, txt, fontsize=9, color=INK, va="center"); ax.text(0.05, yy - 0.04, note, fontsize=7.5, color=col if col == AFTER else INK2, va="center")
    ax.text(0.0, -0.06, "Each term must pay for its own description\n(name + numbers) in bits, or it is dropped.", fontsize=8, color=INK2, transform=ax.transAxes)
    ax = fig.add_subplot(gs[2])
    ax.plot(t, y - y.mean(), ".", ms=3, color=MUTED, label="errors")
    ax.plot(t, res["_mu"] - res["_mu"].mean(), color=AFTER, lw=2.2, label=f"found: onset at {res['winner']['alpha']['onset']:.2f} (planted 2022.0)")
    rem = y - res["_mu"]; ax.plot(t, rem - rem.mean() - 40, ".", ms=3, color=YELLOW, label="what's left (shifted down)")
    ax.set_title("3. What the term explains, and the remainder"); ax.set_xlabel("year"); ax.legend(fontsize=8, loc="lower left"); ax.grid(color=GRID, lw=0.5)
    b0, b1 = res["ladder"][0]["bits_per_reading"], res["winner"]["bits_per_reading"]
    fig.suptitle(f"How the physics compressor works: storing the errors took {b0:.2f} bits per reading; with the onset term, {b1:.2f}", x=0.01, ha="left", fontsize=12, color=INK)
    fig.savefig(FIG / "11_how_it_works.png", dpi=170, bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------------------------------------ 12 validation
def fig_validation():
    A = json.loads((RES / "synthetic.json").read_text())["summary"]
    A0 = json.loads((RES / "synthetic_independent.json").read_text())["summary"] if (RES / "synthetic_independent.json").exists() else None
    fig = plt.figure(figsize=(13.5, 5.6)); gs = fig.add_gridspec(2, 3, width_ratios=[1.5, 1, 1.15], height_ratios=[1, 1], hspace=0.55, wspace=0.42)
    ax = fig.add_subplot(gs[:, 0])
    keys = list(A["planted"]); yy = np.arange(len(keys))[::-1]
    nm = [A["planted"][k]["named_exactly"] for k in keys]; ex = [A["planted"][k]["law_explained"] for k in keys]
    ax.barh(yy + 0.18, ex, height=0.34, color=AFTER, label="law explained (R² ≥ 0.9)")
    ax.barh(yy - 0.18, nm, height=0.34, color=AQUA, label="named exactly")
    ax.axvline(0.8, color=INK2, ls=":", lw=1); ax.axvline(0.7, color=INK2, ls="--", lw=0.8)
    ax.set_yticks(yy); ax.set_yticklabels([k.replace(": ", ":  ") for k in keys], fontsize=8)
    ax.set_xlim(0, 1.02); ax.set_xticks([0, 0.5, 0.7, 0.8, 1]); ax.set_xticklabels(["0", "50 %", "70", "80", "100 %"])
    ax.set_title("Planted laws: found in 20 seeds each?"); ax.legend(fontsize=8, loc="lower left", bbox_to_anchor=(0, -0.16), ncol=2)
    ax = fig.add_subplot(gs[0, 1])
    doms = ["wind", "heat", "flood", "ground"]; x = np.arange(4)
    iid = [A["noise"][f"{d}: iid"]["any_term"] for d in doms]; ar = [A["noise"][f"{d}: ar1"]["any_term"] for d in doms]
    ar0 = [A0["noise"][f"{d}: ar1"]["any_term"] for d in doms] if A0 else None
    if ar0: ax.bar(x - 0.27, ar0, 0.26, color=BEFORE, label="correlated noise, first run")
    ax.bar(x, ar, 0.26, color=AFTER, label="correlated noise, after A1")
    ax.bar(x + 0.27, iid, 0.26, color=AQUA, label="independent noise")
    for i in range(4):
        ax.text(x[i] + 0.27, iid[i] + 0.01, f"{iid[i]:.0%}", ha="center", fontsize=7.5, color=INK2)
        ax.text(x[i], ar[i] + 0.01, f"{ar[i]:.0%}", ha="center", fontsize=7.5, color=INK2)
    ax.axhline(0.05, color=INK2, ls=":", lw=1); ax.set_xticks(x); ax.set_xticklabels(doms)
    ax.set_ylim(0, 0.42); ax.set_yticks([0, 0.05, 0.2, 0.4]); ax.set_yticklabels(["0", "5 %", "20 %", "40 %"])
    ax.set_title("Pure noise: how often it 'finds' a law"); ax.legend(fontsize=7, loc="upper left")
    ax = fig.add_subplot(gs[1, 1])
    for d, c in zip(doms, (BLUE, ORANGE, AQUA, VIOLET)):
        pts = A["sweep"][d]["points"]; ax.plot([p["scale"] for p in pts], [p["explained"] for p in pts], "-o", color=c, ms=4, lw=2, label=d)
    ax.set_xlabel("planted law's size (1 = as stated)"); ax.set_ylabel("explained"); ax.set_ylim(-0.03, 1.05)
    ax.set_title("How small a law can it see?"); ax.legend(fontsize=7.5, ncol=2); ax.grid(color=GRID, lw=0.5)
    ax = fig.add_subplot(gs[:, 2]); ax.axis("off")
    pr = sum(1 for k in keys if A["planted"][k]["law_explained"] >= 0.8 and A["planted"][k]["named_exactly"] >= 0.7)
    lines = [("Pass rules, fixed before the run", INK, 11),
             (f"1. Planted laws: {pr} of {len(keys)} pass", INK, 10),
             ("   (explained ≥ 80 % and named ≥ 70 %)", INK2, 8.5),
             (f"2. Independent noise: {max(iid):.0%} false finds", INK, 10),
             ("   in 400 runs (rule ≤ 5 %): pass", INK2, 8.5),
             (f"3. Correlated noise: {max(ar):.0%} worst domain", INK, 10),
             (f"   (first run: {max(ar0):.0%}; reported, no gate)" if ar0 else "", INK2, 8.5), ("", INK, 6)]
    fails = [k for k in keys if not (A["planted"][k]["law_explained"] >= 0.8 and A["planted"][k]["named_exactly"] >= 0.7)]
    if fails:
        lines.append(("Misses (reported, not tuned away):", INK, 9.5))
        for k in fails:
            v = A["planted"][k]; top = list(v["found_sets"].items())[0]
            lines.append((f"• {k.split(': ')[1]}: named {v['named_exactly']:.0%}", INK2, 8.5))
            lines.append((f"   most often found: {top[0] or '(nothing)'}", MUTED, 8))
    for i, (t, c, s) in enumerate(lines): ax.text(0, 0.98 - i * 0.062, t, color=c, fontsize=s, va="top", transform=ax.transAxes)
    fig.suptitle("Does it find real laws, and stay quiet on noise? 1,080 synthetic cases, truth known", x=0.01, ha="left", fontsize=12.5)
    fig.savefig(FIG / "12_validation.png", dpi=170, bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------------------------------------ 13 wind
def fig_wind():
    W, Wa = load("w1"); L1, _ = load("w2")
    if W is None: return
    fig, axs = plt.subplots(2, 3, figsize=(14, 7.2), gridspec_kw={"width_ratios": [1.2, 1.2, 1.0], "hspace": 0.5, "wspace": 0.3})
    for row, eng in enumerate(("lbm", "rans")):
        r, a = W[eng], Wa[eng]; f = a["fetch"]; y, mu, mu0 = a["y"], a["mu"], a["mu0"]
        ax = axs[row, 0]; ax.scatter(f, y - mu0, s=9, color=BEFORE, alpha=0.6, ec="none")
        ax.axhline(0, color=INK2, lw=0.8, ls=":"); ax.set_ylabel("ln(measured / model)"); ax.set_xlabel("distance into the array (building widths)")
        ax.set_title(f"{eng.upper()}: errors before (sd {r['sd_before']:.2f})"); ax.grid(color=GRID, lw=0.5); ax.set_ylim(-1.4, 1.4)
        ax = axs[row, 1]; ax.scatter(f, y - mu, s=9, color=AFTER, alpha=0.6, ec="none")
        ax.axhline(0, color=INK2, lw=0.8, ls=":"); ax.axhspan(-r["noise_floor_sd"], r["noise_floor_sd"], color=GRID, alpha=0.6, lw=0)
        ax.set_title(f"after the winning terms (sd {r['sd_after']:.2f}; grey = tunnel repeatability)"); ax.set_xlabel("distance into the array (building widths)")
        ax.grid(color=GRID, lw=0.5); ax.set_ylim(-1.4, 1.4)
        ladder_ax(axs[row, 2], r, f"bits saved ({'transfers' if r['transfers'] else 'does not transfer'} across layouts)")
    fig.suptitle("WindSim, AIJ Case C wind tunnel: both engines exaggerate speed-ups and miss the entrance decay", x=0.01, ha="left", fontsize=12.5)
    fig.savefig(FIG / "13_wind_aij.png", dpi=170, bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------------------------------------ 14 heat
def fig_heat():
    H, Ha = load("heat")
    if H is None: return
    fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.3), gridspec_kw={"width_ratios": [1.3, 1.3, 1.0], "wspace": 0.3})
    hrs = np.arange(24)
    for ax, lab, title in ((axs[0], "before_round", "Before the fix (sunlight stamped as instant)"),
                           (axs[1], "after_round", "After the fix (sunlight as hour-ending mean)")):
        a, r = Ha[lab], H[lab]; h = a["hour"]; y, mu = a["y"], a["mu"]
        ax.scatter(h + np.random.default_rng(0).uniform(-0.2, 0.2, len(h)), y, s=8, color=MUTED, alpha=0.45, ec="none")
        ax.plot(hrs, hour_mean(h, y), color=INK, lw=2, label="mean error by hour")
        ax.plot(hrs, hour_mean(h, mu), color=BEFORE if "before" in lab else AFTER, lw=2.2, ls="--",
                label="found: " + " + ".join(nice(t) for t in r["winner"]["terms"]))
        ax.axhline(0, color=INK2, lw=0.8, ls=":"); ax.set_xticks(range(0, 24, 3)); ax.set_xlabel("local hour (EDT)")
        ax.set_ylim(-9, 6); ax.grid(color=GRID, lw=0.5); ax.legend(fontsize=8, loc="lower left")
        c = r["winner"]["coef"]; extra = f"swing ×{1/(1 + c['amplitude']):.2f}" if "amplitude" in c else ""
        if "timing_shift" in c: extra += f", {c['timing_shift']:+.1f} h late"
        ax.set_title(f"{title}\nRMSE {r['rmse_K']:.2f} K; {extra}")
    axs[0].set_ylabel("station − model air temperature (K)")
    ax = axs[2]; labs = [("before", "before,\nfloor clock"), ("after", "after,\nfloor clock"), ("before_round", "before"), ("after_round", "after")]
    x = np.arange(len(labs))
    for i, (k, _) in enumerate(labs):
        c = H[k]["winner"]["coef"]
        ax.bar(i - 0.2, -c.get("amplitude", 0), 0.38, color=VIOLET); ax.bar(i + 0.2, c.get("timing_shift", 0), 0.38, color=YELLOW)
    ax.set_xticks(x); ax.set_xticklabels([l for _, l in labs], fontsize=8)
    ax.bar([], [], color=VIOLET, label="swing too big (share)"); ax.bar([], [], color=YELLOW, label="runs late (hours)")
    ax.axvline(1.5, color=GRID, lw=1); ax.legend(fontsize=8); ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_title("What each run's winning formula says\n(left pair: stations on the wrong clock)")
    fig.suptitle("HeatSim, Philadelphia, 5 hot days, 3 stations: the timing gap is gone after the fix; a 20 % too-large daily swing is left",
                 x=0.01, ha="left", fontsize=12.5)
    fig.savefig(FIG / "14_heat.png", dpi=170, bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------------------------------------ 15 flood
def fig_flood():
    F, Fa = load("flood")
    if F is None: return
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw={"wspace": 0.32})
    ax = axs[0]
    for lab, col in (("before", BEFORE), ("after", AFTER)):
        a = Fa[f"F1_{lab}"]; ax.plot(a["h"] * 1000, a["y"], "o", color=col, ms=6, mec=SURF,
                                      label=f"Cw {F[f'F1_{lab}']['cw']}: median ratio ×{F[f'F1_{lab}']['median_ratio']:.2f}")
    ax.axhline(0, color=INK2, lw=0.8, ls=":"); ax.set_xlabel("depth over the inlet (mm)"); ax.set_ylabel("ln(measured / FloodSim)")
    ax.legend(fontsize=8); ax.grid(color=GRID, lw=0.5); ax.set_title("Sheffield rig: one factor, then nothing")
    ax = axs[1]
    for lab, col in (("before", BEFORE), ("after", AFTER)):
        k = f"F2_{lab}"
        if k in Fa: a = Fa[k]; ax.plot(a["SL"] * 100, a["y"], "o", color=col, ms=6, mec=SURF, label=f"{lab}: " + (" + ".join(nice(t) for t in F[k]["winner"]["terms"]) or "nothing"))
    ax.axhline(0, color=INK2, lw=0.8, ls=":"); ax.set_xlabel("street grade (%)"); ax.set_ylabel("ln(FloodSim share / HEC-22 share)")
    ax.legend(fontsize=8); ax.grid(color=GRID, lw=0.5); ax.set_title("Grate on a grade: the miss tracks the grade")
    ax = axs[2]; keys = [k for k in ("F1_before", "F1_after", "F2_before", "F2_after") if k in F]
    ax.bar(range(len(keys)), [F[k]["sd_before"] for k in keys], 0.38, color=MUTED, label="error sd, raw")
    ax.bar(np.arange(len(keys)) + 0.4, [F[k]["sd_after"] for k in keys], 0.38, color=AFTER, label="after winning terms")
    ax.set_xticks(np.arange(len(keys)) + 0.2); ax.set_xticklabels([k.replace("_", "\n") for k in keys], fontsize=8)
    ax.legend(fontsize=8); ax.set_ylabel("sd of ln error"); ax.grid(axis="y", color=GRID, lw=0.5); ax.set_title("What is left")
    fig.suptitle("FloodSim inlets: the weir coefficient fix removes the only term on the rig", x=0.01, ha="left", fontsize=12.5)
    fig.savefig(FIG / "15_flood.png", dpi=170, bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------------------------------------ 16 ground
def fig_ground():
    Gj, Ga = load("ground")
    if Gj is None: return
    fig, axs = plt.subplots(2, 3, figsize=(14.5, 8), gridspec_kw={"hspace": 0.45, "wspace": 0.3})
    for col, lab in enumerate(("before", "after")):
        a = Ga[f"bad-reference_{lab}"]; n = Gj[f"bad-reference_{lab}"]["network"]
        ax = axs[0, col]; ax.plot(a["t"], a["net"], ".", ms=3, color=MUTED)
        ax.plot(a["t"], a["net_mu"], color=BEFORE if lab == "before" else AFTER, lw=2,
                label="found: " + (" + ".join(nice(t) for t in n["winner"]["terms"]) or "nothing"))
        ax.set_title(f"Faulty reference town, {lab} A11: what every cell shares"); ax.set_xlabel("year"); ax.set_ylabel("median of (series − chain fit), mm")
        ax.legend(fontsize=8); ax.grid(color=GRID, lw=0.5)
        cells = Gj[f"sinkhole_{lab}"]["cells"]; ax = axs[1, col]
        d = np.array([c["dist_km"] for c in cells]) * 1000; on = np.array(["onset" in c["terms"] for c in cells])
        other = np.array([bool(c["terms"]) and "onset" not in c["terms"] for c in cells])
        ax.scatter(d[~on & ~other], np.zeros((~on & ~other).sum()), s=14, color=MUTED, label="nothing left to name")
        ax.scatter(d[on], np.ones(on.sum()), s=14, color=BEFORE if lab == "before" else AFTER, label="onset term wins")
        ax.scatter(d[other], 2 * np.ones(other.sum()), s=14, color=YELLOW, label="another term wins")
        ax.axvline(150, color=INK2, ls=":", lw=1); ax.text(155, 2.3, "A11 bend fit applied\nwithin 150 m", fontsize=7.5, color=INK2)
        ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["none", "onset", "other"]); ax.set_ylim(-0.5, 2.9); ax.set_xlabel("distance from the house (m)")
        ax.set_title(f"Sinkhole town, {lab} A11: what each cell's series still holds"); ax.legend(fontsize=7.5, loc="center right")
    ax = axs[0, 2]; towns = ["quiet", "sinkhole", "regional", "bad-reference"]; x = np.arange(4)
    ax.bar(x - 0.2, [Gj[f"{t}_before"]["rate"]["sd_before"] for t in towns], 0.38, color=BEFORE, label="before A11")
    ax.bar(x + 0.2, [Gj[f"{t}_after"]["rate"]["sd_before"] for t in towns], 0.38, color=AFTER, label="after A11")
    ax.set_xticks(x); ax.set_xticklabels(towns, fontsize=8); ax.set_ylabel("sd of (truth − chain rate), mm/yr"); ax.legend(fontsize=8)
    ax.set_title("Rate error per cell, every town"); ax.grid(axis="y", color=GRID, lw=0.5)
    ax = axs[1, 2]
    hb = [Gj[f"sinkhole_{l}"]["rate"]["house_100m"] for l in ("before", "after")]
    ax.bar([0, 1], [h["chain"] for h in hb], 0.5, color=[BEFORE, AFTER]); ax.axhline(hb[0]["truth"], color=INK, lw=1.5, ls="--")
    ax.text(1.3, hb[0]["truth"], f"truth {hb[0]['truth']:.1f}", va="center", fontsize=8.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["before A11", "after A11"]); ax.set_ylabel("rate at the house (mm/yr)"); ax.set_xlim(-0.5, 1.9)
    ax.set_title("Sinkhole town: the house's rate"); ax.grid(axis="y", color=GRID, lw=0.5)
    fig.suptitle("Ground (InSAR) chain on four made-up towns: A11 removes the jumps and the missed onset", x=0.01, ha="left", fontsize=12.5)
    fig.savefig(FIG / "16_ground.png", dpi=170, bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------------------------------------ 17 scoreboard
def fig_scoreboard():
    rows = []
    W, _ = load("w1"); L1, _ = load("w2"); H, _ = load("heat"); F, _ = load("flood"); Gj, _ = load("ground")
    if W: rows += [("Wind: AIJ tunnel, LBM", W["lbm"], None), ("Wind: AIJ tunnel, RANS", W["rans"], None)]
    if L1: rows += [("Wind: storm gusts, thunder", L1["thunderstorm"], None), ("Wind: storm gusts, ordinary", L1["ordinary"], None)]
    if H: rows += [("Heat: Philadelphia", H["before_round"], H["after_round"])]
    if F: rows += [("Flood: Sheffield rig", F["F1_before"], F["F1_after"]), ("Flood: grate on a grade", F["F2_before"], F.get("F2_after"))]
    if Gj:
        for t in ("sinkhole", "bad-reference"):
            rows += [(f"Ground: {t} town (rate)", Gj[f"{t}_before"]["rate"], Gj[f"{t}_after"]["rate"])]
    fig, ax = plt.subplots(figsize=(13, 0.62 * len(rows) + 1.6))
    yy = np.arange(len(rows))[::-1]
    for y, (lab, b, a) in zip(yy, rows):
        sb = 1 - b["sd_after"] / b["sd_before"]
        ax.barh(y + 0.18, sb, 0.34, color=BEFORE if a else MUTED)
        ax.text(sb + 0.01, y + 0.18, " + ".join(nice(t) for t in b["winner"]["terms"]) or "nothing", va="center", fontsize=8, color=INK2)
        if a:
            sa = 1 - a["sd_after"] / a["sd_before"]
            ax.barh(y - 0.18, sa, 0.34, color=AFTER)
            ax.text(sa + 0.01, y - 0.18, " + ".join(nice(t) for t in a["winner"]["terms"]) or "nothing left to name", va="center", fontsize=8, color=INK2)
    ax.set_yticks(yy); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlim(0, 1.15); ax.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax.set_xticklabels(["0", "25 %", "50 %", "75 %", "100 %"])
    ax.set_xlabel("share of the error spread the winning physics terms explain")
    ax.barh([], [], color=BEFORE, label="before today's fix"); ax.barh([], [], color=AFTER, label="after the fix"); ax.barh([], [], color=MUTED, label="no fix today")
    ax.legend(fontsize=8, loc="lower right"); ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_title("Every model: how much of its error the library can name, before and after today's fixes (the rest is remainder)", fontsize=11.5)
    fig.savefig(FIG / "17_scoreboard.png", dpi=170, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    sys.path.insert(0, str(ROOT))
    for f in (fig_how, fig_validation, fig_wind, fig_heat, fig_flood, fig_ground, fig_scoreboard):
        try: f(); print("ok", f.__name__)
        except FileNotFoundError as e: print("skip", f.__name__, e)
