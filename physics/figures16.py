"""Protocol 16 figure: Case D's per-point miss on the district map, and the compressor's distance term."""
import json, sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "physics")); sys.path.insert(0, str(ROOT))
import district_fields as D

hm, o, dx = D.heights(); res = json.loads((D.PHYS / "district_fields.json").read_text())
out = ROOT / "physics" / "figures" / "district"; out.mkdir(parents=True, exist_ok=True)
for run, d in res["runs"].items():
    x, y, r, mu = map(np.array, (d["x"], d["y"], d["r"], d["mu"])); dt = np.array(d["geometry"]["d_tower"])
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.2), dpi=150)
    ext = [o[0], o[0] + hm.shape[0] * dx, o[1], o[1] + hm.shape[1] * dx]
    ax[0].imshow(np.where(hm > 0, np.where(hm > 0.1, 2, 1), 0).T, origin="lower", extent=ext, cmap="Greys", vmax=2.5)
    sc = ax[0].scatter(np.r_[x, x], np.r_[y, -y], c=np.exp(np.r_[r, r]), cmap="RdBu", vmin=0.3, vmax=1.7, s=34, edgecolors="k", lw=0.5)
    ax[0].set_xlim(-0.25, 0.25); ax[0].set_ylim(-0.2, 0.2); ax[0].set_aspect("equal")
    ax[0].annotate("wind", xy=(-0.17, 0.17), xytext=(-0.24, 0.17), arrowprops=dict(arrowstyle="->"), va="center")
    ax[0].set_title("Simulated ÷ measured at 2 m (red = too slow)\nblack = tower, grey = low blocks"); fig.colorbar(sc, ax=ax[0], shrink=0.8)
    o_ = np.argsort(dt)
    ax[1].axhspan(0.75, 1.25, color="#d9ead3", label="hit band (±25 %)")
    ax[1].scatter(dt, np.exp(r), c="k", s=18, label="probe (mirror pairs averaged)")
    ax[1].plot(dt[o_], np.exp(mu[o_]), "r-", lw=2, label="compressor's best formula: distance from tower")
    ax[1].set_xlabel("distance from the tower, in tower widths"); ax[1].set_ylabel("simulated ÷ measured"); ax[1].set_ylim(0, 1.5)
    ax[1].set_title(f"The only word the compressor finds is distance from the tower\n(sheltering, street-canyon and line-of-sight words add nothing) — {run}")
    ax[1].legend(fontsize=8, loc="lower left"); fig.tight_layout(); fig.savefig(out / f"district_{run}.png")
