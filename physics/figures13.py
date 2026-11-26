"""Figure for protocol 13 (physics/FINDINGS-13.md): python physics/figures13.py"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; F = ROOT / "physics/figures/inflow"; F.mkdir(parents=True, exist_ok=True)
INK, MUTE, ACC, WARN = "#1f2933", "#7b8794", "#2f6fdf", "#d9822b"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTE})
d = json.loads((ROOT / "results/physics/inflow_fields.json").read_text())
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8)); fig.subplots_adjust(wspace=0.3, top=0.8)
L = d["level"]; names = [("R0", "standard"), ("N1", "standard,\nrepeat"), ("C1", "standard +\ninflow stress"), ("W1", "wide +\ninflow stress")]
vals = [L[k] for k, _ in names]
a1.bar(range(4), vals, color=[MUTE, MUTE, ACC, WARN], width=0.6)
for i, v in enumerate(vals): a1.text(i, v + 0.01, f"×{v:.3f}", ha="center", color=INK)
a1.axhline(1, color=INK, lw=0.8); a1.text(3.35, 1.01, "tunnel", color=INK, fontsize=9)
a1.set_xticks(range(4)); a1.set_xticklabels([n for _, n in names]); a1.set_ylim(0.9, 1.35); a1.set_ylabel("street speed, model ÷ tunnel (all 120 points)")
a1.set_title("1. Adding shear stress to the incoming wind barely moves\nthe 1.25× overshoot (a repeat run moves it by 0.02 %)", loc="left", color=INK)
for f, lab, c in (("f3s6_R0", "standard inflow (R0)", MUTE), ("f3s7_C1", "inflow with stress (C1)", ACC)):
    z = np.load(ROOT / f"results/physics/{f}_fields.npz"); st = z["mean_sums"].astype(float); o = z["origin"]; dx = float(z["dx_m"])
    x = o[0] + (np.arange(st.shape[0]) + .5) * dx; j = st.shape[1] // 2; k = int((0.2 - o[2]) / dx); m = x < -0.5
    a2.plot(x[m], -(st[m, j, k, 8] - st[m, j, k, 0] * st[m, j, k, 2]), color=c, lw=2, label=lab)
a2.axhline(0.0265, color=INK, ls=":", lw=1); a2.text(-0.79, 0.0465, "target for a tunnel boundary layer (dotted)", fontsize=8.5, color=INK)
a2.axvline(-0.8, color=MUTE, lw=0.7); a2.text(-0.795, -0.03, "where step 6\nread it", fontsize=8.5, color=MUTE)
a2.axvline(-0.5, color=MUTE, lw=0.7); a2.text(-0.545, -0.03, "first\nblocks", fontsize=8.5, color=MUTE, ha="right")
a2.set_xlabel("distance along the tunnel, m (inlet at −1.0)"); a2.set_ylabel("−u'w' at roof height, m²/s²"); a2.legend(frameon=False, fontsize=8.5, loc="lower left")
a2.set_title("2. Both runs build the right stress by the time the wind reaches\nthe blocks; step 6's 'no stress' was read too close to the inlet", loc="left", color=INK)
fig.suptitle("AIJ Case C 1D, LBM D/16, WindSim F3 step 7 (25 Nov 2026): the incoming wind's turbulence is not the cause", x=0.01, ha="left", fontsize=11.5, color=INK, y=1.0)
fig.savefig(F / "inflow_01_step7.png", dpi=160, bbox_inches="tight"); print(F / "inflow_01_step7.png")
