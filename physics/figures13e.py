"""Pawnee V2 rerun with the bridged-cells fix (insar_pipeline PR #22): where the stitched ground is, before/after grades."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from physics.figures13 import cmap, INK, INK2, GRID_C, BLUE, ORANGE

def main(out):
    d = np.load("data/insar_m2/fix_rerun/v2_fix_points.npz"); lat0, lon0 = float(np.nanmean(d["lat"])), float(np.nanmean(d["lon"]))
    e = (d["lon"] - lon0) * 111.32 * np.cos(np.radians(lat0)); n = (d["lat"] - lat0) * 110.57; r = d["rate_mm_yr"]; br = d["bridged"]
    fig = plt.figure(figsize=(15, 8.6), facecolor="#fcfcfb")
    fig.suptitle("Pawnee rerun with the stitched-piece fix: the false patch at OK02 is caught, but two other checks still fail",
                 fontsize=14, color=INK, x=0.02, ha="left", y=0.975)
    for k, (title, show) in enumerate((("Ground movement it measured, mm/yr", "rate"), ("Where the ground is a separate stitched piece", "br"))):
        ax = fig.add_axes([0.03 + k * 0.245, 0.14, 0.23, 0.74]); ok = np.isfinite(r)
        if show == "rate":
            sc = ax.scatter(e[ok], n[ok], c=r[ok] - np.nanmedian(r), s=0.5, cmap=cmap, vmin=-8, vmax=8, rasterized=True)
            cb = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.04, pad=0.08)
            cb.set_label("blue = away from satellite, red = toward", fontsize=8.5)
        else:
            ax.scatter(e[ok & ~br], n[ok & ~br], s=0.5, color="#d9d8d3", rasterized=True)
            ax.scatter(e[ok & br], n[ok & br], s=0.5, color=ORANGE, rasterized=True)
            ax.text(0.02, -0.1, "orange = stitched on by the software (12% of points):\nnow labelled, wider error bars, not graded",
                    transform=ax.transAxes, fontsize=8.5, color=INK2, va="top")
        for i, s in enumerate(d["stations"]):
            sx = (d["station_lon"][i] - lon0) * 111.32 * np.cos(np.radians(lat0)); sy = (d["station_lat"][i] - lat0) * 110.57
            ax.plot(sx, sy, "^", ms=9, mfc="white", mec=INK, mew=1.4); ax.text(sx + 0.6, sy + 0.5, s, fontsize=8.5, color=INK, weight="bold")
        ax.set_aspect("equal"); ax.set_facecolor("#bdbcb7"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title, loc="left", fontsize=10.5, color=INK)
    ax = fig.add_axes([0.52, 0.30, 0.46, 0.58]); ax.axis("off")
    rows = [("Check", "Before the fix", "With the fix"),
            ("Pair rate error", "5.7 mm/yr  caveat", "4.8 mm/yr  PASS"),
            ("Station pairs agreeing", "8 of 15  FAIL", "6 of 10  FAIL"),
            ("Noise on still ground", "7.75 mm/yr  FAIL", "7.75 mm/yr  FAIL"),
            ("Bias", "+0.13 mm/yr  PASS", "+0.33 mm/yr  PASS"),
            ("Error bars hold the truth", "6 of 6  PASS", "5 of 5  PASS"),
            ("OK02 (the false patch)", "off by 3 to 10 mm/yr", "caught: labelled, not graded")]
    for i, row in enumerate(rows):
        y = 1 - i * 0.125
        for j, (txt, x) in enumerate(zip(row, (0.0, 0.42, 0.70))):
            col = INK if i == 0 or j == 0 else (ORANGE if "FAIL" in txt or "off by" in txt else (INK2 if "caveat" in txt else BLUE))
            ax.text(x, y, txt, fontsize=10.5, color=col, weight="bold" if i == 0 or "PASS" in txt or "FAIL" in txt else "normal", va="top", transform=ax.transAxes)
        ax.plot([0, 1], [y - 0.09, y - 0.09], color=GRID_C, lw=0.8, transform=ax.transAxes)
    fig.text(0.52, 0.27, "What still fails is not stitched ground. The pairs still out all involve OK04, OK05 or OK07:\n"
             "OK04-OK05 -7.9, OK05-OK06 +6.7, OK04-OK07 -6.6, OK03-OK05 -5.5 mm/yr.\n"
             "OK05 has GPS on only 57 of 127 dates. The station series also jump together\n"
             "at a few of the 13 processing batches, which points at the batch joins next.",
             fontsize=9.5, color=INK, va="top")
    fig.text(0.02, 0.01, "Method 2 = insar_pipeline's own processing (MiaplPy, snaphu, MintPy), Sentinel-1 ascending track 34, 127 dates, 6 public NGL GPS stations. "
             "Same rates as before (MintPy kept its inversion); what changed is the stitched-piece label, its error bars and the grading.",
             fontsize=7.8, color=INK2)
    fig.savefig(out, dpi=150); print(out)

if __name__ == "__main__":
    import sys; main(sys.argv[1])
