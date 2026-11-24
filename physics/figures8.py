"""Figures for protocol 8 (physics/FINDINGS-8.md): python physics/figures8.py"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm
ROOT = Path(__file__).resolve().parents[1]; F = ROOT / "physics/figures"
d = json.loads((ROOT / "results/physics/lateral_fields.json").read_text()); a = np.load(ROOT / "results/physics/lateral_fields_arrays.npz")
INK, MUTE, ACC, WARN = "#1f2933", "#7b8794", "#2f6fdf", "#d9822b"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTE})


def blocks(ax):
    for cx in (-0.4, 0.0, 0.4):
        for cy in (-0.4, 0.0, 0.4):
            tall = cx == 0 and cy == 0
            ax.add_patch(Rectangle((cx - 0.1, cy - 0.1), 0.2, 0.2, color="#52606d" if tall else "#cbd2d9", zorder=0))
    ax.set_xlim(-0.55, 0.55); ax.set_ylim(-0.55, 0.55); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.annotate("wind", xy=(-0.33, 0.52), xytext=(-0.53, 0.52), arrowprops=dict(arrowstyle="->", color=INK), va="center", color=INK)


# 1: the error at each point, before and after the terms
x, y, r = a["x"], a["y"], a["r"]
fig, axs = plt.subplots(1, 3, figsize=(14, 4.8)); fig.subplots_adjust(right=0.86, wspace=0.08, top=0.82)
nrm = TwoSlopeNorm(vcenter=0, vmin=-0.3, vmax=0.3)
panels = [(r - r.mean(), f"Model vs tunnel, level removed (sd {r.std():.3f})"),
          (r - a["muA"], f"After the step on speed-up (sd {d['A']['sd_after']:.3f})"),
          (r - a["muB"], f"After adding downwash × lateral (sd {d['B']['sd_after']:.3f})")]
for ax, (v, t) in zip(axs, panels):
    blocks(ax); sc = ax.scatter(x, y, c=v, cmap="RdBu_r", norm=nrm, s=38, edgecolor="white", lw=0.4, zorder=2)
    ax.set_title(t, loc="left", fontsize=10.5, color=INK)
cb = fig.colorbar(sc, cax=fig.add_axes([0.875, 0.15, 0.012, 0.6])); cb.set_ticks(np.log([0.75, 0.87, 1, 1.15, 1.33]) - 0)
cb.set_ticklabels(["model 25 % low", "13 % low", "on the level", "15 % high", "33 % high"])
fig.suptitle(f"AIJ Case C 1D, LBM D/32: every point reads about ×{d['level_factor']:.2f} high; what's left after that, point by point "
             "\n(dark block = the tall centre block; red = model faster than the tunnel; tunnel repeat noise ≈ 0.05)", x=0.01, ha="left", fontsize=11, color=INK)
fig.savefig(F / "lateral_01_points.png", dpi=160, bbox_inches="tight")

# 2: which flow words pay, and the level ledger
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]; bi = d["bootstrap_by_input"]["B"]
lab = {"lp": "model's own speed-up", "lateral": "distance from centreline", "downwash": "vertical flow at 5 cm",
       "ti_local": "local turbulence", "fetch": "distance into array", "side_flux": "outward flow", "street_flux": "air through the street",
       "turn": "flow turning", "reverse": "reversed flow", "x_pos": "x position", "corner": "near tall corner", "cross": "cross street"}
field = {"downwash", "ti_local", "side_flux", "street_flux", "turn", "reverse"}
keys = [k for k in lab if k in bi or k in field]; keys = sorted(keys, key=lambda k: -bi.get(k, 0))
ax.barh(range(len(keys)), [bi.get(k, 0) for k in keys], color=[WARN if k in field else MUTE for k in keys], height=0.6)
ax.axvline(0.6, color=INK, lw=0.8, ls="--"); ax.text(0.61, len(keys) - 0.6, "protocol mark 60 %", fontsize=8.5, color=INK)
ax.set_yticks(range(len(keys)), [lab[k] for k in keys]); ax.invert_yaxis(); ax.set_xlim(0, 1)
ax.set_xlabel("share of 30 resamples whose best formula uses it")
ax.set_title("Vertical flow joins distance from centreline; it doesn't replace it\n(orange = from the model's own flow field)", loc="left", fontsize=10.5, color=INK)
ax = axs[1]; led = d["ledger"]["Q_over_Q_upstream_by_x"]; xs = np.array([float(k) for k in led]); qs = np.array(list(led.values()))
ax.axvspan(-0.5, 0.5, color="#e4e7eb", zorder=0); ax.text(0, 1.03, "block array", ha="center", color=MUTE, fontsize=9)
ax.plot(xs, qs, "o-", color=ACC); ax.axhline(1, color=INK, lw=0.8)
ax.set_ylim(0.4, 1.1); ax.set_xlabel("x (m), wind from the left"); ax.set_ylabel("air flow below 10 cm / upstream")
ax.set_title(f"Inside the array the model's lowest 10 cm carries ~⅔ of the air;\nnet flow at 10 cm is upward (+{d['ledger']['plan_mean_w_top_over_Uapp']*100:.1f} % of approach speed)",
             loc="left", fontsize=10.5, color=INK)
fig.tight_layout(); fig.savefig(F / "lateral_02_terms_ledger.png", dpi=160); print("ok")
