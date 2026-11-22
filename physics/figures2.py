"""Figures for physics/FINDINGS-2.md (flood, heat, ground) and the cross-model scoreboard."""
import json, math
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results" / "physics"; FIG = ROOT / "physics" / "figures"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
S1, S2, S3, S4, S5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "savefig.facecolor": SURF})
fl = json.loads((RES / "flood.json").read_text()); ht = json.loads((RES / "heat.json").read_text())
gd = json.loads((RES / "ground.json").read_text()); w = json.loads((RES / "aij_c.json").read_text())

# 06 heat: error by local hour, per station, and the best formula
d = ht["data"]; hl = np.array(d["hour_local"]); r = np.array(d["r"]); mu = np.array(d["mu_best"]); sid = np.array(d["sid"])
fig, ax = plt.subplots(figsize=(8.5, 3.8))
ax.axhline(0, color=INK2, lw=1, ls=":")
for (s, nm), col in zip(ht["stations"].items(), (S1, S2, S3)):
    m = sid == s
    ax.scatter(hl[m] + (np.random.default_rng(1).uniform(-0.15, 0.15, m.sum())), r[m], s=14, color=col, alpha=0.55, ec="none", label=nm)
hrs = np.arange(24); mm = [np.mean(mu[np.round(hl) % 24 == h]) for h in hrs]; rr = [np.mean(r[np.round(hl) % 24 == h]) for h in hrs]
ax.plot(hrs, rr, color=INK, lw=2, label="mean error, all stations")
ax.plot(hrs, mm, color=S4, lw=2, ls="--", label="best formula (sun 2 h earlier + hour + light wind)")
ax.set_xlabel("local hour (EDT)"); ax.set_ylabel("station − model (K)"); ax.set_xticks(range(0, 24, 3))
ax.legend(frameon=False, fontsize=8, ncol=2, loc="lower left"); ax.grid(axis="y", color=GRID, lw=0.5)
ax.set_title("HeatSim open-site air temperature, station correction off: the daily swing is too large and runs late",
             loc="left", fontsize=10.5, color=INK)
fig.tight_layout(); fig.savefig(FIG / "06_heat_diurnal.png", dpi=170); plt.close(fig)

# 07 flood: Sheffield rig and grate on a grade
f1, f2 = fl["F1_sheffield"]["data"], fl["F2_grate_on_grade"]["data"]
fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
ax = axs[0]; h = np.array(f1["h_mm"])
ax.plot(h, f1["floodsim_ls"], "o", color=S2, ms=6, mec=SURF, label="FloodSim (weir, coefficient 0.5)")
ax.plot(h, f1["measured_ls"], "o", color=S1, ms=6, mec=SURF, label="measured, Sheffield rig")
c = math.exp(fl["F1_sheffield"]["one_correction"]["fit_coef"]["c"])
ax.plot(h, np.array(f1["floodsim_ls"]) * c, color=INK, lw=1.5, label=f"FloodSim × {c:.2f} (the only term that pays)")
ax.set_xlabel("depth over the inlet (mm)"); ax.set_ylabel("flow into the sewer (l/s)"); ax.legend(frameon=False, fontsize=8)
ax.set_title("Sheffield rig: right shape, coefficient 55 % too high", loc="left", fontsize=10, color=INK); ax.grid(color=GRID, lw=0.5)
ax = axs[1]; SL = np.array(f2["SL"]); rain = np.array(f2["rain"])
for R, col in zip((50, 100, 200), (S1, S3, S2)):
    m = rain == R; ax.plot(100 * SL[m], np.array(f2["floodsim"])[m] / np.array(f2["hec22"])[m], "o-", color=col, ms=7, mec=SURF, lw=1.5, label=f"{R} mm/h")
ax.set_xlabel("street grade (%)"); ax.set_ylabel("FloodSim share / HEC-22 share"); ax.legend(frameon=False, fontsize=8, title="rain", title_fontsize=8)
ax.set_title("Grate on a grade: the miss tracks the grade, not the flow", loc="left", fontsize=10, color=INK); ax.grid(color=GRID, lw=0.5)
ax.set_xticks([1, 2, 4])
fig.tight_layout(); fig.savefig(FIG / "07_flood_exchange.png", dpi=170); plt.close(fig)

# 08 ground: sinkhole town, error before and after the formula
g = gd["sinkhole"]["data"]; x, y = np.array(g["x"]), np.array(g["y"]); r = np.array(g["r"]); mu = np.array(g["mu_best"])
fig, axs = plt.subplots(1, 2, figsize=(10, 4.4))
for ax, e, ttl in zip(axs, (r, r - mu), ("truth − straight-line velocity", "left after: velocity + acceleration + near the house")):
    sc = ax.scatter(x / 1000, y / 1000, c=np.clip(e, -4, 4), cmap="RdBu", vmin=-4, vmax=4, s=40, ec=SURF, lw=0.5)
    ax.plot(0.32, -0.21, marker="*", color=INK, ms=12); ax.set_aspect("equal"); ax.set_xlim(-0.9, 1.5); ax.set_ylim(-1.4, 1.0); ax.set_title(ttl, loc="left", fontsize=10, color=INK)
    ax.set_xlabel("km"); ax.text(0.40, -0.35, "house", fontsize=8, color=INK)
cb = fig.colorbar(sc, ax=axs, shrink=0.8); cb.set_label("mm/yr (red: chain reads too little sinking)")
fig.suptitle("Ground chain, synthetic sinkhole town: a straight line over 2016-2025 misses a bowl that started in 2022", x=0.02, ha="left", fontsize=11, color=INK)
fig.savefig(FIG / "08_ground_sinkhole.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# 09 scoreboard: bits per point, one correction vs best formula vs neighbours, every case
cases = [("Wind: LBM D/10", w["step2_lbm"]["bias_only"]["bits_per_reading"], w["step2_lbm"]["top"][0]["bits_per_reading"], w["step1_lbm"]["bits_per_reading"], True),
         ("Wind: RANS D/10", w["step2_rans"]["bias_only"]["bits_per_reading"], w["step2_rans"]["top"][0]["bits_per_reading"], w["step1_rans"]["bits_per_reading"], True),
         ("Heat: Philadelphia, 5 days", ht["one_correction"]["bits_per_reading"], ht["top"][0]["bits_per_reading"], ht["neighbour"]["bits_per_reading"], ht["holdout"][0]["transfers"]),
         ("Flood: Sheffield rig", fl["F1_sheffield"]["one_correction"]["bits_per_reading"], fl["F1_sheffield"]["top"][0]["bits_per_reading"], None, False),
         ("Flood: grate on a grade", fl["F2_grate_on_grade"]["one_correction"]["bits_per_reading"], fl["F2_grate_on_grade"]["top"][0]["bits_per_reading"], None, fl["F2_grate_on_grade"]["holdout"][0]["transfers"]),
         ("Ground: sinkhole town", gd["sinkhole"]["one_correction"]["bits_per_reading"], gd["sinkhole"]["top"][0]["bits_per_reading"], gd["sinkhole"]["neighbour"]["bits_per_reading"], True),
         ("Ground: regional town", gd["regional"]["one_correction"]["bits_per_reading"], gd["regional"]["top"][0]["bits_per_reading"], gd["regional"]["neighbour"]["bits_per_reading"], True)]
fig, ax = plt.subplots(figsize=(9, 4.2)); yy = np.arange(len(cases))[::-1]
for (nm, b0, bf, bn, tr), y in zip(cases, yy):
    save = b0 - bf
    ax.barh(y, save, color=S1 if tr else GRID, height=0.55, ec=SURF)
    if bn is not None: ax.plot([b0 - bn], [y], marker="D", color=S4, ms=7, mec=SURF)
    ax.text(max(save, 0) + 0.03, y, f"{save:+.2f}" + ("" if tr else "  (does not hold on every hold-out)" if save > 0 else "  (one factor is best)"),
            va="center", fontsize=8.5, color=INK2)
ax.set_yticks(yy, [c[0] for c in cases]); ax.axvline(0, color=INK2, lw=1)
ax.set_xlabel("bits per point saved over 'model × one correction' (more = more missing structure named)")
ax.plot([], [], "D", color=S4, label="saved by guessing from neighbours (no formula)"); ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.set_title("What a short formula names in each model's errors", loc="left", fontsize=11, color=INK); ax.set_xlim(-0.1, 2.6)
fig.tight_layout(); fig.savefig(FIG / "09_scoreboard.png", dpi=170); plt.close(fig)
print("ok")
