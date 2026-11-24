"""Figures for protocol 5 (physics/FINDINGS-5.md): python physics/figures5.py"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; F = ROOT / "physics" / "figures"
d = json.loads((ROOT / "results/physics/aij_c_bulk.json").read_text()); a = np.load(ROOT / "results/physics/aij_c_bulk_arrays.npz")
INK, MUTE, ACC, WARN, OK = "#1f2933", "#7b8794", "#2f6fdf", "#d9822b", "#3f9142"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTE})

# 1. ledger: share of ln 1.30
L = d["ledgers"]; tot = -L["level_ln"]
rows = [("Approach flow near the ground", L["S1_inflow"]["share"], "measured"),
        ("Evaluation height (bound)", L["S2_height"]["share_bound"], "bound"),
        ("Domain blockage (bound)", L["S3_domain"]["share_blockage_bound"], "bound"),
        ("Made inside the array (what is left)", L["made_inside_array_share"], "rest")]
fig, ax = plt.subplots(figsize=(7.2, 2.9))
for i, (lab, s, kind) in enumerate(rows):
    ax.barh(i, s * 100, color={"measured": ACC, "bound": "none", "rest": WARN}[kind], edgecolor=ACC if kind == "bound" else "none",
            hatch="//" if kind == "bound" else None, height=0.6)
    ax.text(s * 100 + 1.5, i, f"{s*100:.0f} %", va="center", color=INK)
ax.set_yticks(range(len(rows)), [r[0] for r in rows]); ax.invert_yaxis(); ax.set_xlim(0, 100)
ax.set_xlabel("share of the 1.30× excess (log terms)")
ax.set_title("Most of the 1.30× is made inside the block array", loc="left", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(F / "bulk_01_ledger.png", dpi=160); plt.close(fig)

# 2. run levels vs settings
lv = d["run_levels"]; names = list(lv)
fig, ax = plt.subplots(figsize=(7.2, 3.2))
col = [OK if v["z0"] else MUTE for v in lv.values()]
ax.bar(range(len(names)), [v["median_ratio_model_over_measured"] for v in lv.values()], color=col, width=0.65)
for i, v in enumerate(lv.values()): ax.text(i, v["median_ratio_model_over_measured"] + 0.01, f"{v['median_ratio_model_over_measured']:.2f}", ha="center", fontsize=9)
ax.axhline(1.0, color=INK, lw=0.8); ax.set_ylim(0.9, 1.55)
lab = {"d16": "D/16", "d16_z0tunnel": "D/16\ntunnel z0", "d20": "D/20", "d24": "D/24", "d32_z0tunnel": "D/32", "d32_walls": "D/32\nwalls",
       "d32_L0.5D": "D/32\nL 0.5D", "d32_cs0.17": "D/32\nCs 0.17"}
ax.set_xticks(range(len(names)), [lab[n] for n in names], fontsize=8.5)
ax.set_ylabel("model / measured (median)")
ax.set_title("Eight runs, one excess: only the ground roughness setting (green) moves it", loc="left", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(F / "bulk_02_runs.png", dpi=160); plt.close(fig)

# 3. the D/32 pattern on the array, and against RANS
fig, axs = plt.subplots(1, 2, figsize=(9.4, 4.1))
ax = axs[0]
for bx in (-0.4, 0, 0.4):
    for by in (-0.4, 0, 0.4): ax.add_patch(plt.Rectangle((bx - 0.1, by - 0.1), 0.2, 0.2, color="#cbd2d9"))
pat = -(a["r32"] - a["r32"].mean()) * 100
sc = ax.scatter(a["x"], a["y"], c=pat, cmap="RdBu_r", vmin=-25, vmax=25, s=46, edgecolor="white", lw=0.4)
ax.set_aspect("equal"); ax.set_xlim(-0.55, 0.55); ax.set_ylim(-0.55, 0.55); ax.set_xlabel("x (m), wind →"); ax.set_ylabel("y (m)")
plt.colorbar(sc, ax=ax, shrink=0.8, label="model faster than its own 1.30× (%)")
ax.set_title("What is left after the 1.30×: no nameable shape", loc="left", fontsize=10.5, color=INK)
ax = axs[1]; rp = -(a["rans"] - a["rans"].mean())
ax.scatter(-(a["r32"] - a["r32"].mean()), rp, s=16, color=ACC, alpha=0.75)
ax.axhline(0, color=MUTE, lw=0.6); ax.axvline(0, color=MUTE, lw=0.6)
ax.set_xlabel("LBM D/32 pattern, ln"); ax.set_ylabel("RANS D/20 pattern, ln")
s6 = d["S6_rans"]
ax.set_title(f"LBM and RANS miss in different places (r = {s6['corr_lbm_rans']:.2f})", loc="left", fontsize=10.5, color=INK)
fig.tight_layout(); fig.savefig(F / "bulk_03_pattern.png", dpi=160); plt.close(fig)
print("ok")
