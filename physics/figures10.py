"""Protocol 10 figure: method 2 (V2 Pawnee) on stable ground, for a reader with no InSAR background."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from physics.insar_m2 import load, design, rsd

INK, INK2, GRID_C = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
cmap = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#ee9a99", "#a12a2a"])
cmap.set_bad("#bdbcb7")

def main(out):
    D = load(); v = D["valid"]; e, n = D["east_km"], D["north_km"]; t = D["t_yr"]; st = list(D["stations"])
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
    fig = plt.figure(figsize=(15, 10.5), facecolor="#fcfcfb")
    fig.suptitle("Our own InSAR processing (method 2) on ground that isn't moving: Pawnee, Oklahoma test site, Oct 2021 to Nov 2025",
                 fontsize=14, color=INK, x=0.02, ha="left", y=0.985)
    # A map
    ax = fig.add_axes([0.05, 0.1, 0.4, 0.8]); ax.set_facecolor("#fcfcfb")
    step = 1.0; ex = np.arange(-16, 16.01, step); ny = np.arange(-16.5, 16.01, step)
    M = np.full((len(ny) - 1, len(ex) - 1), np.nan); r = D["rate_mm_yr"]
    ix = np.digitize(e[v], ex) - 1; iy = np.digitize(n[v], ny) - 1
    for a in range(M.shape[1]):
        for b in range(M.shape[0]):
            m = (ix == a) & (iy == b)
            if m.sum() >= 30: M[b, a] = np.median(r[v][m])
    off = np.nanmedian(M)
    im = ax.pcolormesh(ex, ny, np.ma.masked_invalid(M - off), cmap=cmap, vmin=-6, vmax=6)
    for i, s in enumerate(st):
        sx = (D["station_lon"][i] - D["lon"].mean()) * 111.32 * np.cos(np.radians(D["lat"].mean())); sy = (D["station_lat"][i] - D["lat"].mean()) * 110.57
        ax.plot(sx, sy, "^", ms=11, mfc="white", mec=INK, mew=1.5); ax.text(sx + 0.6, sy + 0.5, s, fontsize=10, color=INK, weight="bold")
    ax.set_aspect("equal"); ax.set_xlabel("km east"); ax.set_ylabel("km north")
    ax.set_title("A. What method 2 says the ground is doing (should be about 0 everywhere)", loc="left", fontsize=11, color=INK)
    ax.annotate("this corner reads 4 to 7 mm/yr sinking,\nbut the GPS station there (OK02) says it isn't",
                xy=(-9, 8.5), xytext=(-15.5, 13.6), fontsize=9.5, color=INK, arrowprops=dict(arrowstyle="-", color=INK2))
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.04, pad=0.08)
    cb.set_label("mm per year: blue = moving away from the satellite (mostly sinking), red = toward it (rising)\ngrey = no trusted reading (the processing could not follow the signal there)", fontsize=8.5)
    # B stations
    ax = fig.add_axes([0.53, 0.58, 0.44, 0.31]); ax.set_facecolor("#fcfcfb")
    G = D["station_gnss_los_mm"]; R = D["station_radar_mean_200m_mm"]; X = design(t)
    def rate(y):
        m = np.isfinite(y); return np.linalg.lstsq(X[m], y[m], rcond=None)[0][1]
    names = [s for s in st if s != "OK03"]; i0 = st.index("OK03")
    rr = [rate(R[st.index(s)] - R[i0]) for s in names]; gg = [rate(G[st.index(s)] - G[i0]) for s in names]
    x = np.arange(len(names)); w = 0.36
    ax.bar(x - w / 2 - 0.01, rr, w, color=BLUE, label="satellite (method 2)"); ax.bar(x + w / 2 + 0.01, gg, w, color=ORANGE, label="GPS station")
    ax.axhline(0, color=INK2, lw=0.8); ax.set_xticks(x, [s + ("\n(GPS has gaps)" if s == "OK05" else "") for s in names])
    ax.set_ylabel("mm/yr, compared with OK03"); ax.grid(axis="y", color=GRID_C); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    ax.set_title("B. Station check: the satellite and GPS agree at OK04 and OK06, not at OK02 and OK07", loc="left", fontsize=11, color=INK)
    for xi, a, b in zip(x, rr, gg):
        ax.text(xi - w / 2, a + (0.25 if a >= 0 else -0.55), f"{a:+.1f}", ha="center", fontsize=8.5, color=INK)
        ax.text(xi + w / 2, b + (0.25 if b >= 0 else -0.55), f"{b:+.1f}", ha="center", fontsize=8.5, color=INK)
    # C error bars
    ax = fig.add_axes([0.68, 0.33, 0.29, 0.15]); ax.set_facecolor("#fcfcfb")
    labs = ["error it states", "real spread, points 0.3 km apart", "real spread, places 10+ km apart"]; vals = [0.72, 2.4, 3.0]
    ax.barh(range(3)[::-1], vals, color=[ORANGE, BLUE, BLUE], height=0.6)
    for k, (l, vv) in enumerate(zip(labs, vals)): ax.text(vv + 0.05, 2 - k, f"{vv:.1f} mm/yr", va="center", fontsize=9, color=INK)
    ax.set_yticks(range(3)[::-1], labs); ax.set_xlim(0, 3.8); ax.set_xlabel("mm/yr")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title("C. Its error bars are about 3 times too narrow", loc="left", fontsize=11, color=INK)
    # D lag
    ax = fig.add_axes([0.53, 0.07, 0.44, 0.17]); ax.set_facecolor("#fcfcfb")
    Y = D["series_mm"][v][::4]; Y = Y - np.median(Y, 0); lags = [1, 2, 3, 5, 8, 12, 20, 30, 40, 60, 80, 100, 120]
    rms = [np.sqrt(np.mean((Y[:, k:] - Y[:, :-k]) ** 2)) for k in lags]; days = np.array(lags) * 12 / 365.25
    ax.plot(days, rms, color=BLUE, lw=2, marker="o", ms=5, label="what the data does")
    ax.plot(days, [rms[0]] * len(days), color=ORANGE, lw=2, ls="--", label="what its error bars assume")
    ax.set_xlabel("years between the two dates"); ax.set_ylabel("typical change, mm"); ax.grid(color=GRID_C); ax.set_ylim(0, 16)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    ax.set_title("D. Why: the noise slowly drifts instead of averaging out, so it builds up into fake trends", loc="left", fontsize=11, color=INK)
    fig.text(0.02, 0.012, "Source: insar_pipeline V2 run (MiaplPy, snaphu, MintPy), 24 Nov 2026; GPS: NGL IGS20 daily, along the satellite's view. "
             "The map is shifted so its middle reads 0. Test site, not a customer property. gridzip-lab protocol 10.", fontsize=8, color=INK2)
    fig.savefig(out, dpi=150); print(out)

if __name__ == "__main__":
    import sys; main(sys.argv[1])
