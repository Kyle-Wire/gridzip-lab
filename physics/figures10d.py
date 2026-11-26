"""Las Vegas BR1 run A (method 2 with the bridged-cells fix): map, station check, scorecard against V2 Pawnee."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from physics.figures10 import cmap, INK, INK2, GRID_C, BLUE, ORANGE
from physics.insar_m2 import design

def rate(t, y):
    m = np.isfinite(y); return np.linalg.lstsq(design(t[m]), y[m], rcond=None)[0][1]

def main(out):
    d = np.load("data/insar_lv/br1_a_points.npz"); lat0, lon0 = float(np.nanmean(d["lat"])), float(np.nanmean(d["lon"]))
    e = (d["lon"] - lon0) * 111.32 * np.cos(np.radians(lat0)); n = (d["lat"] - lat0) * 110.57; r = d["rate_mm_yr"]
    t = np.array([(np.datetime64(f"{x[:4]}-{x[4:6]}-{x[6:]}") - np.datetime64(f"{d['dates'][0][:4]}-{d['dates'][0][4:6]}-{d['dates'][0][6:]}")).astype(int) / 365.25 for x in d["dates"]])
    fig = plt.figure(figsize=(15, 8.6), facecolor="#fcfcfb")
    fig.suptitle("Our own InSAR processing (method 2) passes its first fair test: Las Vegas, with the stitched-piece fix", fontsize=14, color=INK, x=0.02, ha="left", y=0.975)
    ax = fig.add_axes([0.05, 0.12, 0.38, 0.76]); ok = np.isfinite(r)
    sc = ax.scatter(e[ok], n[ok], c=r[ok] - np.nanmedian(r), s=0.6, cmap=cmap, vmin=-8, vmax=8, rasterized=True)
    for i, s in enumerate(d["stations"]):
        sx = (d["station_lon"][i] - lon0) * 111.32 * np.cos(np.radians(lat0)); sy = (d["station_lat"][i] - lat0) * 110.57
        ax.plot(sx, sy, "^", ms=10, mfc="white", mec=INK, mew=1.4); ax.text(sx + 0.5, sy + 0.4, s, fontsize=9, color=INK, weight="bold")
    ax.set_aspect("equal"); ax.set_xlabel("km east"); ax.set_ylabel("km north"); ax.set_facecolor("#bdbcb7")
    ax.set_title("Ground movement it measured, mm/yr (real movement, not noise)", loc="left", fontsize=11, color=INK)
    cb = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.04, pad=0.09); cb.set_label("blue = away from the satellite (mostly sinking), red = toward it (rising); grey = no trusted reading", fontsize=8.5)
    # station pairs
    ax = fig.add_axes([0.5, 0.52, 0.22, 0.36]); st = list(d["stations"])
    R = [rate(t, d["station_radar_mean_200m_mm"][i]) for i in range(len(st))]; G = [rate(t, d["station_gnss_los_mm"][i]) for i in range(len(st))]
    xs, ys = [], []
    for a in range(len(st)):
        for b in range(a + 1, len(st)): xs.append(G[a] - G[b]); ys.append(R[a] - R[b])
    lim = max(map(abs, xs + ys)) * 1.15
    ax.plot([-lim, lim], [-lim, lim], color=INK2, lw=1, ls="--"); ax.fill_between([-lim, lim], [-lim - 2, lim - 2], [-lim + 2, lim + 2], color=GRID_C, alpha=0.6, lw=0)
    ax.scatter(xs, ys, s=40, color=BLUE, edgecolor="white", linewidth=1.5, zorder=3)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.set_xlabel("GPS says, mm/yr"); ax.set_ylabel("satellite says, mm/yr")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title("All 21 GPS station pairs: the satellite\nmatches the GPS (grey band = within 2 mm/yr)", loc="left", fontsize=10.5, color=INK)
    # scorecard
    ax = fig.add_axes([0.5, 0.06, 0.48, 0.36]); ax.axis("off")
    rows = [("Check", "Pawnee (V2, earlier site)", "Las Vegas (with the fix)"),
            ("Station pairs agreeing (pipeline test)", "8 of 15  FAIL", "21 of 21  PASS"),
            ("Pair rate error", "5.7 mm/yr  caveat", "1.1 mm/yr  PASS"),
            ("Noise on still ground", "7.75 mm/yr  FAIL", "2.05 mm/yr  PASS"),
            ("Bias", "+0.13 mm/yr  PASS", "+0.04 mm/yr  PASS"),
            ("Error bars hold the truth", "6 of 6  PASS", "7 of 7  PASS"),
            ("Our two-station GPS check", "not run", "1.54 vs 1.64 allowed  PASS")]
    for i, row in enumerate(rows):
        y = 1 - i * 0.135
        for j, (txt, x) in enumerate(zip(row, (0.0, 0.36, 0.69))):
            col = INK if i == 0 or j == 0 else (ORANGE if "FAIL" in txt else (INK2 if "caveat" in txt or "not run" in txt else BLUE))
            ax.text(x, y, txt, fontsize=10.5 if i else 10, color=col, weight="bold" if i == 0 or "PASS" in txt or "FAIL" in txt else "normal", va="top", transform=ax.transAxes)
        ax.plot([0, 1], [y - 0.11, y - 0.11], color=GRID_C, lw=0.8, transform=ax.transAxes)
    fig.text(0.02, 0.005, "Method 2 = insar_pipeline's own processing (MiaplPy, snaphu, MintPy), Sentinel-1 descending track 173, 148 dates; 7 public NGL GPS stations. "
             "\nBox from public facts only. Las Vegas has almost no stitched ground (38 of 537,526 points): it shows the fix does no harm, not yet that it helps. Different sites, not a before/after of one place.", fontsize=7.8, color=INK2)
    fig.savefig(out, dpi=150); print(out)

if __name__ == "__main__":
    import sys; main(sys.argv[1])
