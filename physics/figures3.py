"""Figures for physics/FINDINGS-3.md (WindSim Layer 1)."""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results" / "physics"; FIG = ROOT / "physics" / "figures"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"; S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})
ph = json.loads((RES / "layer1_posthoc.json").read_text()); r = json.loads((RES / "layer1.json").read_text())
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.7, 1]})
ax = axs[0]
for k, typ, col, off in ((0, "non-thunderstorm", S1, -0.15), (1, "thunderstorm", S2, 0.15)):
    pl = sorted(ph[typ]["places"]); allp = sorted(set(ph["non-thunderstorm"]["places"]) | set(ph["thunderstorm"]["places"]))
    for p in pl:
        x = allp.index(p) + off
        for s, v in ph[typ]["places"][p]["worst_stations"].items():
            if v["count"] < 4: continue
            pct = 100 * (np.exp(v["mean"]) - 1); ax.plot(x, pct, "o", color=col, ms=5, alpha=0.8, mec=SURF)
            if abs(pct) > 24: ax.text(x + 0.12, pct, s, fontsize=7, color=INK2, va="center")
ax.axhspan(-25, 25, color=GRID, alpha=0.5, lw=0); ax.axhline(0, color=INK2, lw=0.8)
ax.set_xticks(range(len(allp)), [p.replace("south-street", "phila").replace("tract-353-02", "ne-phila").replace("san-francisco", "sf").replace("los-angeles", "la") for p in allp], rotation=40, ha="right")
ax.set_ylabel("station's average miss (%), observed vs chain"); ax.plot([], [], "o", color=S1, label="ordinary storms"); ax.plot([], [], "o", color=S2, label="thunderstorms")
ax.legend(frameon=False, fontsize=8, loc="lower left"); ax.set_title("Each airport keeps its own miss, storm after storm (grey: ±25 %)", loc="left", fontsize=10.5, color=INK)
ax = axs[1]
rows = [("Thunderstorms", r["Q1_thunderstorm"]), ("Ordinary storms", r["Q2_ordinary"])]
lab, sv, nb = [], [], []
for nm, q in rows:
    b0 = q["one_correction"]["bits_per_reading"]; lab.append(nm); sv.append(b0 - q["top"][0]["bits_per_reading"])
    nb.append(q["neighbour_same_station"]["null_mean"] - q["neighbour_same_station"]["bits_per_reading"])
y = np.arange(len(lab))[::-1]
ax.barh(y + 0.18, sv, 0.34, color=S1, label="formula: station's own offset + gust size")
ax.barh(y - 0.18, nb, 0.34, color=S4, label="model-free: same station's other storms")
ca = r["Q2_ordinary"]["named_ca_terms"]; b0 = r["Q2_ordinary"]["one_correction"]["bits_per_reading"]
best_phys = max(b0 - ca[k] for k in ("offshore_ca", "dz_elev", "msin", "lnz0"))
ax.barh([-1], [max(best_phys, 0.002)], 0.34, color=GRID); ax.text(0.01, -1, f"best California physics term: {best_phys:+.3f}", va="center", fontsize=8, color=INK2)
ax.set_yticks(list(y) + [-1], lab + ["California (ordinary)"]); ax.set_xlabel("bits per gust saved over 'chain × one factor'")
ax.set_ylim(-1.6, 2.0); ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=1); ax.set_title("What names the misses", loc="left", fontsize=10.5, color=INK)
fig.tight_layout(); fig.savefig(FIG / "10_layer1_stations.png", dpi=170); plt.close(fig); print("ok")
