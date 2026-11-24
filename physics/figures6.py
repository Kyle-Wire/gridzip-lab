"""Figures for protocol 6 (physics/FINDINGS-6.md): python physics/figures6.py"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; F = ROOT / "physics/figures"
d = json.loads((ROOT / "results/physics/storm_search.json").read_text()); a = np.load(ROOT / "results/physics/storm_arrays.npz", allow_pickle=True)
INK, MUTE, ACC, WARN = "#1f2933", "#7b8794", "#2f6fdf", "#d9822b"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTE})
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.1, 1]})
ax = axs[0]
lab = {"cold_pool": "Cold-pool push", "gust_front": "Gust front passed", "pressure_jump": "Pressure jump", "dew_drop": "Dewpoint drop",
       "dryness": "Dry air before", "heat": "Afternoon heat", "elevation": "Station height", "direction": "Wind direction",
       "stability": "Time of day", "season": "Season"}
s = d["savings_beyond_station_bits"]; names = sorted(s, key=lambda k: -s[k])
storm = {"cold_pool", "gust_front", "pressure_jump", "dew_drop", "dryness", "heat"}
ax.barh(range(len(names)), [s[n] for n in names], color=[WARN if n in storm else MUTE for n in names], height=0.6)
ax.axvline(0, color=INK, lw=0.8); ax.set_yticks(range(len(names)), [lab[n] for n in names]); ax.invert_yaxis()
ax.set_xlabel("bits saved beyond each station's own offset (below 0 = doesn't pay)")
ax.set_title("No storm term pays its way (orange = new)", loc="left", fontsize=11, color=INK)
ax = axs[1]; y = a["y"]; gf = a["gf"] > 0.5
b = np.linspace(-1.0, 0.8, 37)
ax.hist(np.exp(-0) * y[gf], bins=b, alpha=0.75, color=ACC, label=f"gust front passed ({gf.sum()})", density=True)
ax.hist(y[~gf], bins=b, alpha=0.6, color=MUTE, label=f"no front seen ({(~gf).sum()})", density=True)
for v, c in ((np.median(y[gf]), ACC), (np.median(y[~gf]), INK)): ax.axvline(v, color=c, lw=1.2, ls="--")
t = np.log([0.5, 0.7, 1, 1.4, 2]); ax.set_xticks(t, ["0.5×", "0.7×", "1×", "1.4×", "2×"])
ax.set_xlabel("measured gust / chain's gust"); ax.set_yticks([]); ax.legend(frameon=False, fontsize=9)
ax.set_title("Behind a front: ~10 % higher, lost in the spread", loc="left", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(F / "storm_01_terms.png", dpi=160); print("ok")
