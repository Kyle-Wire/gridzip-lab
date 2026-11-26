"""Figure for protocol 10 (physics/FINDINGS-10.md): python physics/figures10.py"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
ROOT = Path(__file__).resolve().parents[1]; F = ROOT / "physics/figures/roof"; F.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "physics"))
import lateral_fields as LF  # noqa: E402
d = json.loads((ROOT / "results/physics/roof_fields.json").read_text()); a = np.load(ROOT / "results/physics/roof_fields_arrays.npz")
INK, MUTE, ACC, WARN = "#1f2933", "#7b8794", "#2f6fdf", "#d9822b"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTE})


def blocks(ax):
    from matplotlib.patches import Rectangle
    for cx in (-0.4, 0.0, 0.4):
        for cy in (-0.4, 0.0, 0.4):
            ax.add_patch(Rectangle((cx - 0.1, cy - 0.1), 0.2, 0.2, color="#cbd2d9", zorder=0))
    ax.set_xlim(-0.55, 0.55); ax.set_ylim(-0.55, 0.55); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.annotate("wind", xy=(-0.33, 0.52), xytext=(-0.53, 0.52), arrowprops=dict(arrowstyle="->", color=INK), va="center", color=INK)
C0, C1 = ACC, WARN

fig, axs = plt.subplots(1, 3, figsize=(15, 4.9)); fig.subplots_adjust(wspace=0.42, top=0.8)

# 1: approach flow at x = -0.8, y = 0, in both sections
ax = axs[0]; ax2 = ax.twiny()
for run, c in (("R0", C0), ("R1", C1)):
    f = np.load(ROOT / f"results/physics/f3s6_{run}_fields.npz"); st = f["mean_sums"].astype(float); o = f["origin"].astype(float); dx = float(f["dx_m"])
    z = np.linspace(0.01, 0.44, 60); P = np.c_[np.full(60, -0.8), np.zeros(60), z]
    v = LF.interp(st[..., [0, 2, 8]], o, dx, P); u, w, uw = v[:, 0], v[:, 1], v[:, 2]
    uref = float(LF.interp(st[..., [0]], o, dx, np.array([[-0.8, 0, 0.2]]))[0, 0])
    lab = "standard section" if run == "R0" else "wider, taller section"
    ax.plot(u, z / 0.2, color=c, lw=2, label=f"speed, {lab}")
    ax2.plot(-(uw - u * w) / uref ** 2 * 1000, z / 0.2, color=c, lw=1.4, ls="--", label=f"stress ×1000, {lab}")
ax.axhline(1, color=MUTE, lw=0.6); ax.text(2.3, 1.03, "roof height", color=MUTE, fontsize=8.5)
ax.set_xlabel("mean speed, m/s (solid)"); ax2.set_xlabel("−u'w' / U_ref² × 1000 (dashed)", color=INK); ax.set_ylabel("height / block height D")
ax.legend(loc="upper left", fontsize=8, frameon=False); ax2.legend(loc="lower right", fontsize=8, frameon=False)
ax.set_title("1. The wind coming in is not the same in the two runs\n(5.1 % apart, limit 3 %), so blockage can't be judged",
             loc="left", fontsize=10.5, color=INK, pad=34)

# 2: per-point change between the two runs
ax = axs[1]; blocks(ax); dy = a["y0"] - a["y1"]
sc = ax.scatter(a["x"], a["y"], c=-dy, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-0.2, vmax=0.2), s=40, edgecolor="white", lw=0.4, zorder=2)
cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02); cb.set_ticks(np.log([0.85, 1, 1.18])); cb.set_ticklabels(["15 % lower", "same", "18 % higher"])
ax.set_title(f"2. Overall level barely moves (×{d['level']['R0']:.3f} → ×{d['level']['R1']:.3f}),\n"
             f"but points move ±{dy.std() * 100:.0f} %: wider vs standard, per point", loc="left", fontsize=10.5, color=INK)

# 3: the compressor's false-alarm test and its power
ax = axs[2]; syn = d["synthetic"]; sds = [0.02, 0.04, 0.06, 0.08, 0.10]
for nm, c in (("roof_speed", C0), ("roof_stress", C1)):
    ax.plot(sds, [syn["power"][nm]["found_by_sd"][str(s)] for s in sds], "o-", color=c, label=f"planted {nm.replace('_', ' ')}")
ax.axhline(8, color=MUTE, lw=0.8, ls=":"); ax.text(0.021, 8.2, "found in 8 of 10 = detectable", fontsize=8.5, color=MUTE)
fa = syn["a_none"]["of"] - syn["a_none"]["no_roof_term"]
ax.bar([0.0], [fa], width=0.012, color="#c0392b"); ax.text(0.0, fa + 0.3, f"{fa} false\nalarms\n(limit 1)", ha="center", fontsize=8.5, color="#c0392b")
ax.set_xticks([0.0] + sds); ax.set_xticklabels(["nothing\nplanted"] + [f"{s:.2f}" for s in sds]); ax.set_ylim(0, 10.5)
ax.set_xlabel("planted effect size (sd of ln speed)"); ax.set_ylabel("runs out of 10"); ax.legend(fontsize=8, frameon=False, loc="center right")
ax.set_title("3. Self-test failed: with nothing planted it named\n'roof shear' twice, so its roof answer isn't trusted", loc="left", fontsize=10.5, color=INK)

fig.suptitle("AIJ Case C 1D, LBM D/16 (WindSim F3 step 6, 24 Nov 2026): wider tunnel and flow over the roofs, read by the compressor (protocol 10)",
             x=0.01, y=1.06, ha="left", fontsize=11.5, color=INK)
fig.savefig(F / "roof_01_step6.png", dpi=160, bbox_inches="tight")
print(F / "roof_01_step6.png")
