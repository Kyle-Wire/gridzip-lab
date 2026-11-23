"""Figure for the synthetic street test: python surprise/figures_street.py <street dir> <fig dir>"""
import json, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sd, fd = sys.argv[1:3]; d = json.load(open(f"{sd}/street.json")); r = json.load(open(f"{sd}/street_result.json"))
nx, ny, W = d["nx"], d["ny"], d["W"]; a = np.fromfile(f"{sd}/street_100yr.f32", np.float32).reshape(ny, nx); house = np.array(d["house"]).reshape(ny, nx)
INK, MUTED, SURF = "#1f1f1e", "#6b6b68", "#fcfcfb"; plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
cmap = LinearSegmentedColormap.from_list("s", ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6.6), gridspec_kw={"height_ratios": [1.1, 1]})
img = np.ma.masked_where(house > 0, np.clip(a, 0, 0.3)); ax1.imshow(img, origin="lower", cmap=cmap, vmin=0, vmax=0.3, extent=[0, nx, 0, ny])
ax1.imshow(np.ma.masked_where(house == 0, house), origin="lower", cmap=LinearSegmentedColormap.from_list("h", ["#3b423e", "#3b423e"]), extent=[0, nx, 0, ny])
P = d["planted"]; names = {P["low_spot"]: "low spot in back yard", P["neighbour_fill"]: "neighbour's fill next door", P["low_lot"]: "whole lot 25 cm lower"}
for q in range(d["lots"]):
    ax1.axvline(q * W, color="#dde2de", lw=0.6); ax1.text(q * W + W / 2, 60, str(q), ha="center", color=MUTED, fontsize=8)
for q, n in names.items(): ax1.text(q * W + W / 2, 46, n, ha="center", color="#e34948", fontsize=8, weight="bold", wrap=True)
ax1.set_xticks([]); ax1.set_yticks([]); ax1.set_title("Synthetic street, 100-yr 1-hr storm: deepest water (blue, up to 30 cm); houses dark; street at the bottom", loc="left", color=INK, fontsize=10)
s = np.array(r["surprise"]); cols = ["#e34948" if q in names else "#2a78d6" for q in range(len(s))]
ax2.bar(range(len(s)), s, color=cols, width=0.6); ax2.axhline(0, color=MUTED, lw=0.8)
for q in range(len(s)): ax2.text(q, s[q] + (0.3 if s[q] >= 0 else -0.9), f"{s[q]:.1f}", ha="center", fontsize=8, color=INK)
ax2.set_xticks(range(len(s))); ax2.set_xlabel("lot (red = planted)"); ax2.set_ylabel("surprise, bits")
for sp in ("top", "right"): ax2.spines[sp].set_visible(False)
ax2.set_title("The score flags the low-spot yard (11.8 bits, first by far).\nThe other two plants barely changed the modelled water at the house, so there was nothing to flag.", loc="left", color=INK, fontsize=10)
fig.tight_layout(); fig.savefig(f"{fd}/04_synthetic_street.png", dpi=160); print("ok")
