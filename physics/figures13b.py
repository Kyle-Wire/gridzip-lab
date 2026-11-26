"""Protocol 13 part D figure: method 1 vs method 2 at Pawnee, for a reader with no InSAR background."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from physics.figures13 import cmap, INK, INK2, GRID_C, BLUE, ORANGE

def main(out):
    r = json.load(open("results/physics/insar_m1_vs_m2.json")); T = np.load("results/physics/insar_m1_tiles.npz")
    ex = np.arange(-16, 16.01, 1.0); ny = np.arange(-16.5, 16.01, 1.0)
    fig = plt.figure(figsize=(15, 9.6), facecolor="#fcfcfb")
    fig.suptitle("Pawnee, Oklahoma (ground that isn't moving): the ready-made product (method 1) against our own processing (method 2)",
                 fontsize=14, color=INK, x=0.02, ha="left", y=0.98)
    for k, (key, lab, col) in enumerate([("M1", "Method 1 (ready-made OPERA product, what reports use)", ORANGE), ("M2", "Method 2 (our own processing)", BLUE)]):
        ax = fig.add_axes([0.03 + k * 0.31, 0.36, 0.29, 0.54]); M = T[key]
        im = ax.pcolormesh(ex, ny, np.ma.masked_invalid(M - np.nanmedian(M)), cmap=cmap, vmin=-10, vmax=10)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(lab, loc="left", fontsize=11, color=INK)
        sp = r["D3"][f"m{k+1}_rate_rsd"]; ax.text(0.02, 0.02, f"typical false movement: {sp:.1f} mm/yr", transform=ax.transAxes, fontsize=10.5,
                                                   color=INK, weight="bold", bbox=dict(fc="white", ec="none", alpha=0.85))
    cb = fig.colorbar(im, cax=fig.add_axes([0.05, 0.32, 0.55, 0.015]), orientation="horizontal")
    cb.set_label("mm per year it says the ground moves (should be 0 everywhere). Blue = sinking, red = rising, grey = no trusted reading. 32 km box.", fontsize=8.5)
    # right: three comparisons
    def bars(ax, vals, title, unit, better="lower"):
        ax.barh([1, 0], vals, color=[ORANGE, BLUE], height=0.55)
        for i, v in zip([1, 0], vals): ax.text(v * 1.02 + 0.05, i, f"{v:g} {unit}", va="center", fontsize=9.5, color=INK)
        ax.set_yticks([1, 0], ["method 1", "method 2"]); ax.set_xlim(0, max(vals) * 1.45)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        ax.set_title(title, loc="left", fontsize=10.5, color=INK); ax.set_xticks([])
    bars(fig.add_axes([0.71, 0.74, 0.26, 0.14]), [round(r["D3"]["m1_rate_rsd"], 1), round(r["D3"]["m2_rate_rsd"], 1)], "False movement across the area (lower is better)", "mm/yr")
    bars(fig.add_axes([0.71, 0.52, 0.26, 0.14]), [round(r["D3"]["m1_spread_0p2_0p4"], 1), round(r["D3"]["m2_spread_0p2_0p4"], 1)], "Disagreement between points 300 m apart", "mm/yr")
    p1, p2 = r["D2"]["m1_pairs_ok"][0], r["D2"]["m2_pairs_ok"][0]
    bars(fig.add_axes([0.71, 0.30, 0.26, 0.14]), [p1, p2], "GPS station pairs matched within 2 mm/yr (of 15)", "of 15")
    ax = fig.add_axes([0.06, 0.07, 0.9, 0.17]); lags = [1, 10, 40, 120]; days = np.array(lags) * 12 / 365.25
    ax.plot(days, [r["D4"]["m1"][str(l)] for l in lags], color=ORANGE, lw=2, marker="o", ms=6, label="method 1")
    ax.plot(days, [r["D4"]["m2"][str(l)] for l in lags], color=BLUE, lw=2, marker="o", ms=6, label="method 2")
    ax.set_xlabel("years between two dates"); ax.set_ylabel("typical false change, mm"); ax.grid(color=GRID_C); ax.legend(frameon=False, loc="upper left")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title("How fast false movement builds up over time (lower is better)", loc="left", fontsize=10.5, color=INK)
    fig.text(0.02, 0.01, "Same satellite track, box and 126 dates. Method 1: OPERA DISP-S1, no weather correction, GNSS tie off. Method 2: insar_pipeline V2 (MiaplPy) with HRRR weather correction. "
             "Maps shifted so their middle reads 0. gridzip-lab protocol 13, part D.", fontsize=8, color=INK2)
    fig.savefig(out, dpi=150); print(out)

if __name__ == "__main__":
    import sys; main(sys.argv[1])
