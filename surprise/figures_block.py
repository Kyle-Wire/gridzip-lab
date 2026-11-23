"""Figures for the South Street surprise run: python surprise/figures_block.py <CHI folder> <work dir> <fig dir> [--mark]
--mark outlines the sample home (for the private copy only; the committed figures never mark it)."""
import json, sys, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from pyproj import Transformer
chi, work, figd = sys.argv[1:4]; mark = "--mark" in sys.argv; os.makedirs(figd, exist_ok=True)
SEQ = ["#f0efec", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
cmap = LinearSegmentedColormap.from_list("seq", SEQ); INK, MUTED, SURF = "#1f1f1e", "#6b6b68", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": SURF, "axes.facecolor": SURF})
d = np.load(f"{work}/flood_block_south_street.npz"); res = json.load(open(f"{work}/flood_block_south_street.json"))
fp = json.load(open(f"{chi}/data/raw/south-street/li_footprints.p0.geojson"))["features"]
t = Transformer.from_crs("EPSG:4326", "EPSG:32618", always_xy=True)
polys = []
for i in d["idx"]:
    g = fp[i]["geometry"]; ring = np.array(g["coordinates"][0]); x, y = t.transform(ring[:, 0], ring[:, 1]); polys.append(np.column_stack([x, y]))
s = d["surprise"]; x0, y0 = min(p[:, 0].min() for p in polys), min(p[:, 1].min() for p in polys)
polys = [p - [x0, y0] for p in polys]
fig, ax = plt.subplots(figsize=(9, 6.4))
norm = Normalize(0, 16)
pc = PolyCollection(polys, array=np.clip(s, 0, 16), cmap=cmap, norm=norm, edgecolors=SURF, linewidths=0.4); ax.add_collection(pc)
if mark and res.get("homes"):
    ky = json.load(open(f"{chi}/private/property-key.json")); hx, hy = t.transform(ky["lon"], ky["lat"])
    ax.plot(hx - x0, hy - y0, "o", ms=14, mfc="none", mec="#e34948", mew=2); ax.annotate("PHL-POC-001", (hx - x0, hy - y0), xytext=(12, 12), textcoords="offset points", color=INK, fontsize=10, weight="bold")
ax.set_aspect("equal"); ax.autoscale_view(); ax.set_xlabel("metres east"); ax.set_ylabel("metres north")
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
cb = fig.colorbar(pc, ax=ax, shrink=0.6, pad=0.02); cb.set_label("surprise (extra bits vs neighbours within 30 m)"); cb.outline.set_visible(False)
ax.set_title("South Street: how unlike its neighbours each house's modelled flooding is\nFloodSim runs of 16 Nov 2026, 6 design storms; 0 = as its neighbours suggest", loc="left", fontsize=11, color=INK)
fig.text(0.01, 0.01, "Screening model output, not a measurement. Grey/light = typical. Each extra bit = half as likely for a typical house.", color=MUTED, fontsize=8)
fig.savefig(f"{figd}/01_block_map{'_marked' if mark else ''}.png", dpi=160, bbox_inches="tight"); plt.close(fig)
if mark: sys.exit()
# 2: real vs shuffled
fig, ax = plt.subplots(figsize=(8, 4.2)); bins = np.arange(-8, 32, 1)
ax.hist(np.clip(d["shuffled"], -8, 31), bins=bins, color="#eb6834", alpha=0.85, label="readings shuffled between houses (location ignored)", rwidth=0.9)
ax.hist(np.clip(s, -8, 31), bins=bins, color="#2a78d6", alpha=0.75, label="real block", rwidth=0.9)
for p, c in ((95, "#2a78d6"),):
    v = res["percentiles"]["95"]; ax.axvline(v, color=c, lw=1.5, ls="--"); ax.text(v + 0.4, ax.get_ylim()[1] * 0.85, f"real 95th pct {v:.1f} bits", color=INK, fontsize=9)
v = res["shuffle"]["percentiles"]["95"]; ax.axvline(v, color="#eb6834", lw=1.5, ls="--"); ax.text(v + 0.4, ax.get_ylim()[1] * 0.7, f"shuffled 95th {v:.1f}", color=INK, fontsize=9)
ax.set_xlabel("surprise (bits)"); ax.set_ylabel("houses"); ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1, 0.6))
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.set_title("Where a house sits matters: the real block has far more odd houses than chance", loc="left", fontsize=11, color=INK)
fig.savefig(f"{figd}/02_real_vs_shuffled.png", dpi=160, bbox_inches="tight"); plt.close(fig)
# 3: planted sensitivity
pl = res["planted"]; ks = list(pl)
fig, ax = plt.subplots(figsize=(7, 3.8))
ax.bar(ks, [pl[k]["median_extra_bits"] for k in ks], color="#2a78d6", width=0.55)
for i, k in enumerate(ks): ax.text(i, pl[k]["median_extra_bits"] + 0.15, f"{pl[k]['median_extra_bits']:.1f} bits\n{100*pl[k]['detected_share_mean']:.0f}% in top 5%", ha="center", fontsize=9, color=INK)
ax.set_ylabel("median extra bits"); ax.set_xlabel("water added around 30 random houses, every storm")
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.set_ylim(0, max(pl[k]["median_extra_bits"] for k in ks) * 1.35)
ax.set_title("Planted wet houses on the real block: the score rises, but slowly", loc="left", fontsize=11, color=INK)
fig.savefig(f"{figd}/03_planted.png", dpi=160, bbox_inches="tight"); plt.close(fig)
print("ok")
