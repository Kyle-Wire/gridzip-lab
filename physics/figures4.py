"""Figure for physics/FINDINGS-4.md (PHL house patch vs rings)."""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; d = json.loads((ROOT / "results/physics/phl_patch.json").read_text())["data"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"; S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})
t, P, N, W = (np.array(d[k]) for k in ("t", "P", "N", "W")); rel, nr, pn = P - W, N - W, P - N
fig = plt.figure(figsize=(12, 6.4)); gs = fig.add_gridspec(2, 3)
ax = fig.add_subplot(gs[0, :2])
ax.plot(t, rel, color=S2, lw=1.1, label="house patch − wide ring (what A11 tests)"); ax.plot(t, nr, color=S1, lw=1.1, alpha=0.8, label="near ring − wide ring")
ax.axvline(2019.5096, color=INK2, ls=":", lw=1); ax.text(2019.56, -11, "2019.51: patch (−10) and near ring (−8) dip together\nagainst the wide ring for one date", fontsize=8, color=INK2)
ax.set_ylabel("mm (line of sight)"); ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=2); ax.grid(color=GRID, lw=0.5)
ax.set_title("The house patch moves with its near ring, including the 2019.51 dip", loc="left", fontsize=10.5, color=INK)
ax = fig.add_subplot(gs[0, 2])
ax.plot(nr, rel, "o", color=S1, ms=3.5, alpha=0.6, mec="none"); x = np.linspace(nr.min(), nr.max(), 2); b = np.polyfit(nr, rel, 1)
ax.plot(x, np.polyval(b, x), color=INK, lw=1.2); i = int(np.argmin(np.abs(t - 2019.5096))); ax.plot(nr[i], rel[i], "o", color=S2, ms=7)
ax.set_xlabel("near ring − wide ring (mm)"); ax.set_ylabel("patch − wide ring (mm)"); ax.grid(color=GRID, lw=0.5)
ax.set_title(f"r = {np.corrcoef(nr, rel)[0, 1]:.2f}, slope {b[0]:.2f}", loc="left", fontsize=10, color=INK)
ax = fig.add_subplot(gs[1, :2])
ax.plot(t, pn, color=S3, lw=1.1, label="house patch − near ring")
dd = pn[1:-1] - (pn[:-2] + pn[2:]) / 2; rs = 1.4826 * np.median(np.abs(dd - np.median(dd))); z = dd / rs
for k in np.where(np.abs(z) > 5)[0]: ax.plot(t[k + 1], pn[k + 1], "o", color=S2, ms=6); ax.text(t[k + 1] + 0.05, pn[k + 1], f"{t[k + 1]:.2f}", fontsize=7, color=INK2)
ax.set_ylabel("mm"); ax.set_ylim(-11, 11); ax.grid(color=GRID, lw=0.5); ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.set_title(f"Against its near ring the house is {rel.std() / pn.std():.1f}× quieter (sd {pn.std():.2f} mm); orange: single-date spikes, |z| > 5", loc="left", fontsize=10.5, color=INK)
ax = fig.add_subplot(gs[1, 2]); r = json.loads((ROOT / "results/physics/phl_patch.json").read_text())
st = r["spike_test"]["single_terms"]; b0 = r["one_correction"]["bits_per_reading"]
names = ["near ring", "spike (one date)", "annual cycle", "rate", "AR(1)", "wide-ring swing", "step"]
vals = [b0 - [t_ for t_ in r["top"] if t_["terms"] == ["nr"]][0]["bits_per_reading"], b0 - st["spike only"]["bpr"], b0 - st["annual only"]["bpr"], b0 - st["t only"]["bpr"], b0 - st["ar1 only"]["bpr"], b0 - st["wdev only"]["bpr"], b0 - st["step only"]["bpr"]]
ax.barh(range(len(names))[::-1], vals, color=[S1 if v > 0 else GRID for v in vals]); ax.axvline(0, color=INK2, lw=0.8)
ax.set_yticks(range(len(names))[::-1], names); ax.set_xlabel("bits per date saved (one term)"); ax.set_title("Only the near ring pays", loc="left", fontsize=10, color=INK)
fig.tight_layout(); fig.savefig(ROOT / "physics/figures/11_phl_patch.png", dpi=170); print("ok")
