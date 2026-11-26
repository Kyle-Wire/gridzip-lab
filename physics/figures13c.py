import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from physics.insar_m2 import load
from physics.figures13 import INK, INK2, BLUE, ORANGE

def main(out):
    D = load(); C = np.load("data/insar_m2/v2_components.npz"); a = C["n_ifgrams_apart"]; v = D["valid"]; e, n = D["east_km"], D["north_km"]; r = D["rate_mm_yr"]
    main_c = v & (a == 0); br = v & (a == 375)
    fig = plt.figure(figsize=(14, 7), facecolor="#fcfcfb")
    fig.suptitle("Found it: the false-sinking corner is one piece of the map that the processing stitched on separately", fontsize=14, color=INK, x=0.02, ha="left")
    ax = fig.add_axes([0.04, 0.08, 0.5, 0.8]); rng = np.random.default_rng(0)
    for m, col, lab in [(main_c, "#bdbcb7", "connected directly to the reference point"), (br, ORANGE, "a separate piece, stitched on afterwards")]:
        idx = np.where(m)[0]; idx = rng.choice(idx, min(len(idx), 15000), replace=False); ax.scatter(e[idx], n[idx], s=1.2, c=col, label=lab, rasterized=True)
    for i, s in enumerate(D["stations"]):
        sx = (D["station_lon"][i] - D["lon"].mean()) * 111.32 * np.cos(np.radians(D["lat"].mean())); sy = (D["station_lat"][i] - D["lat"].mean()) * 110.57
        ax.plot(sx, sy, "^", ms=11, mfc="white", mec=INK, mew=1.5); ax.text(sx + 0.6, sy + 0.5, s, fontsize=10, color=INK, weight="bold")
    ax.set_aspect("equal"); ax.set_xlabel("km east"); ax.set_ylabel("km north"); ax.legend(frameon=False, loc="lower left", markerscale=8, fontsize=9.5)
    ax.set_title("Where each point's reading comes from (method 2, Pawnee)", loc="left", fontsize=11, color=INK)
    ax2 = fig.add_axes([0.63, 0.25, 0.33, 0.5]); vals = [float(np.median(r[main_c])), float(np.median(r[br]))]
    ax2.bar([0, 1], vals, color=["#8f8e89", ORANGE], width=0.55); ax2.axhline(0, color=INK2, lw=0.8)
    for x, val in zip([0, 1], vals): ax2.text(x, val - 0.35, f"{val:+.1f} mm/yr", ha="center", fontsize=10.5, color=INK)
    ax2.set_xticks([0, 1], ["connected\ndirectly", "stitched-on\npiece"]); ax2.set_ylim(-4.8, 0.8)
    for s in ("top", "right"): ax2.spines[s].set_visible(False)
    ax2.set_ylabel("typical false movement, mm/yr (should be 0)")
    ax2.set_title("The stitched piece reads 3 mm/yr of\nsinking that the ground isn't doing", loc="left", fontsize=11, color=INK)
    fig.text(0.02, 0.01, "Grey = in the reference point's unwrapping piece in all 375 radar pairs; orange = apart from it in all 375, joined by MintPy's bridging step. "
             "Weather and height corrections don't cause it (the corner is there before them). gridzip-lab protocol 13, part E.", fontsize=8, color=INK2)
    fig.savefig(out, dpi=150); print(out)

if __name__ == "__main__":
    import sys; main(sys.argv[1])
