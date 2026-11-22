"""Figures for physics/FINDINGS.md from results/physics/*.json and WindSim's per-point files."""
import json, math, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import aij_c as A
from codelen import lawsearch as L

ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results" / "physics"; FIG = ROOT / "physics" / "figures"
FIG.mkdir(exist_ok=True)
V = ROOT.parent / "WindSim" / "docs" / "validation"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
DIV = LinearSegmentedColormap.from_list("div", ["#d03b3b", "#ec9a8f", "#f0efec", "#86b6ef", "#1c5cab"])
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})
main = json.loads((RES / "aij_c.json").read_text()); post = json.loads((RES / "aij_c_posthoc.json").read_text())
floor = main["noise_floor"]


def one(path, col):
    ds = A.dataset(path.parent, {"1D": path.name}, col); r = np.log(ds["m"] / ds["p"])
    lib, ex, cls = A.features(ds["x"], ds["y"], ds["hc"], ds["p"], ds["pk"]); return ds, r, lib, cls


RUNS = {"LBM D/10": (V / "lbm_case_c_d10_blind/points_1D.csv", "vs_pred", "lbm_d10"),
        "LBM D/32": (V / "lbm_case_c_1D_resolution/d32_z0tunnel/points_1D.csv", "vs_pred", "lbm_d32_z0"),
        "RANS D/10": (V / "rans_case_c_d10/points_1D.csv", "pred_vector", "rans_d10_1D"),
        "RANS D/20": (V / "rans_case_c_d20/results_points.csv", "pred_vector", "rans_d20")}


# 1. Where the bits go (layout 1D)
fig, ax = plt.subplots(figsize=(8.2, 3.6))
names = list(RUNS); yy = np.arange(len(names))[::-1]
for k, nm in enumerate(names):
    p = post[RUNS[nm][2]]; b0, bf, bc = p["bias_only_bpr"], p["top"][0]["bpr"], p["spatial"]["bits_per_reading"]
    y = yy[k]
    ax.plot([floor["bits_per_reading"], b0], [y, y], color=GRID, lw=6, solid_capstyle="round", zorder=1)
    ax.scatter([b0], [y], s=90, color=S2, zorder=3, ec=SURF, lw=2, label="model + one correction factor" if k == 0 else None)
    ax.scatter([bc], [y], s=90, color=S4, zorder=3, ec=SURF, lw=2, label="+ each point guessed from its neighbours" if k == 0 else None)
    ax.scatter([bf], [y], s=90, color=S1, zorder=4, ec=SURF, lw=2, label="+ best short formula" if k == 0 else None)
    ax.text(b0 + 0.08, y, f"{b0:.2f}", va="center", fontsize=8.5, color=INK2)
    ax.text(bf - 0.08, y + 0.22, f"{bf:.2f}", va="bottom", ha="center", fontsize=8.5, color=INK2)
ax.axvline(floor["bits_per_reading"], color=INK2, lw=1, ls="--")
ax.text(floor["bits_per_reading"] + 0.05, yy[0] + 0.55, f"measurement noise floor {floor['bits_per_reading']:.2f}\n(mirror pairs, sd {100*floor['sd_log']:.1f} %)",
        fontsize=8, color=INK2, va="top")
ax.set_yticks(yy, names); ax.set_xlabel("bits per measured point to store the model's error (fewer = more explained)")
ax.set_xlim(4.4, 8.8); ax.set_ylim(-0.6, len(names) - 0.1); ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.legend(loc="lower left", frameon=False, fontsize=8.5, bbox_to_anchor=(0, 1.0), ncol=3)
ax.set_title("AIJ Case C, layout 1D: how much of each engine's error is structure", loc="left", fontsize=11, color=INK, pad=26)
fig.tight_layout(); fig.savefig(FIG / "01_bits_ladder.png", dpi=170); plt.close(fig)


# 2. Residual maps
def blocks(ax, hc=1):
    for bx in (-0.4, 0, 0.4):
        for by in (-0.4, 0, 0.4):
            if (bx, by) == (0, 0) and hc == 0: continue
            ax.add_patch(Rectangle((bx - 0.1, by - 0.1), 0.2, 0.2, fc="#d9d8d3", ec=INK2, lw=0.6))
    ax.set_xlim(-0.52, 0.52); ax.set_ylim(-0.36, 0.36); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    ax.annotate("", xy=(-0.42, 0.33), xytext=(-0.52, 0.33), arrowprops=dict(arrowstyle="->", color=INK2))
    ax.text(-0.52, 0.345, "wind", fontsize=7, color=INK2)


panels = []
for nm in ("LBM D/10", "LBM D/32", "RANS D/10"):
    ds, r, lib, cls = one(RUNS[nm][0], RUNS[nm][1])
    b0 = L.score(lib, [], r, A.G); terms = post[RUNS[nm][2]]["top"][0]["terms"]
    terms = [t for t in terms if "hc" not in t]; bf = L.score(lib, terms, r, A.G)
    panels.append((nm, ds, r - b0["_mu"], r - bf["_mu"], terms, b0["fit_coef"]["c"]))
fig, axs = plt.subplots(2, 3, figsize=(12, 6.4), gridspec_kw={"hspace": 0.35})
norm = TwoSlopeNorm(0, -0.6, 0.6)
for j, (nm, ds, e0, ef, terms, c) in enumerate(panels):
    for i, (e, lab) in enumerate(((e0, f"after one factor (model x {math.exp(-c):.2f} too {'fast' if c < 0 else 'slow'})"),
                                  (ef, "after best formula: " + (" + ".join(terms) if terms else "none pays for itself")))):
        ax = axs[i, j]; blocks(ax)
        sc = ax.scatter(ds["x"], ds["y"], c=np.clip(e, -0.6, 0.6), cmap=DIV, norm=norm, s=58, ec=SURF, lw=0.8, zorder=3)
        ax.set_title(f"{nm}\n{lab}" if i == 0 else lab, fontsize=8.5, color=INK, loc="left")
        ax.text(0.52, -0.355, f"sd {100*np.std(e):.0f} %", fontsize=8, color=INK2, ha="right")
cb = fig.colorbar(sc, ax=axs, shrink=0.7, pad=0.02); cb.set_label("ln(measured / model) left over\nred: model too fast, blue: too slow", fontsize=8.5)
fig.suptitle("Where each engine's error sits in the streets (z = 0.02 m, layout 1D)", x=0.02, ha="left", fontsize=11, color=INK)
fig.savefig(FIG / "02_residual_maps.png", dpi=170, bbox_inches="tight"); plt.close(fig)


# 3. The recovered law: measured ~ model^b
fig, axs = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
for ax, nm, col in zip(axs, ("LBM D/10", "LBM D/32", "RANS D/20"), (S1, S3, S2)):
    ds, r, lib, cls = one(RUNS[nm][0], RUNS[nm][1]); lp = lib["lp"]
    keep = ds["p"] > 0.25                      # drawn: the 4 RANS points near |U| = 0 sit off the chart
    for c_, mk in (("along", "o"), ("cross", "s"), ("intersection", "^")):
        s = (cls == c_) & keep
        ax.scatter(np.exp(lp[s]), ds["m"][s] / A.U_APP, s=26, marker=mk, color=col, alpha=0.8, ec=SURF, lw=0.6, label=c_)
    k = np.polyfit(lp[keep], np.log(ds["m"][keep] / A.U_APP), 1); xx = np.linspace(0.1, 1.6, 50)
    ax.plot(xx, np.exp(np.polyval(k, np.log(xx))), color=INK, lw=2)
    ax.plot(xx, xx, color=INK2, lw=1, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.12, 1.7); ax.set_ylim(0.12, 1.7)
    tk = [0.15, 0.25, 0.5, 1.0, 1.5]
    ax.set_xticks(tk, [str(t) for t in tk]); ax.set_yticks(tk, [str(t) for t in tk]); ax.minorticks_off()
    ax.set_title(f"{nm}: measured ∝ model^{k[0]:.2f}", loc="left", fontsize=10, color=INK)
    ax.set_xlabel("model speed / approach speed"); ax.grid(color=GRID, lw=0.5, which="both")
    ax.text(0.13, 1.45, "dotted: model right", fontsize=7.5, color=INK2)
axs[0].set_ylabel("measured speed / approach speed"); axs[0].legend(frameon=False, fontsize=8, loc="lower right")
fig.suptitle("The term both engines need: their street-to-street contrasts are too strong (exponent 1 = right; RANS fit leaves out 4 points where |U| ~ 0)",
             x=0.02, ha="left", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(FIG / "03_contrast_law.png", dpi=170); plt.close(fig)


# 4. Resolution: the shape fades, a uniform factor stays
rs = main["step4_resolution"]; N = [10, 16, 16, 20, 24, 32]; z0 = ["z0" in r["run"] for r in rs]
b = [1 + r["fit_coef"]["lp"] for r in rs]
ent = [100 * (math.exp(r["fit_coef"]["exp_xi"] * math.exp(-1 / r["fit_alpha"]["exp_xi"])) - 1) for r in rs]
save = [r["bias_only_bits_per_reading"] - r["bits_per_reading"] for r in rs]
fac = []
for r_, (lab, f) in zip(rs, [("D/10", V / "lbm_case_c_d10_blind/points_1D.csv")] + [(None, V / "lbm_case_c_1D_resolution" / d / "points_1D.csv") for d in ("d16", "d16_z0tunnel", "d20", "d24", "d32_z0tunnel")]):
    d = np.genfromtxt(f, delimiter=",", names=True); fac.append(float(np.exp(np.mean(np.log(d["vs_pred"] / d["vs_obs"])))))
fig, axs = plt.subplots(1, 4, figsize=(12.5, 3.3))
for ax, val, ttl, ref in zip(axs, (b, ent, fac, save),
                             ("contrast exponent b\n(1 = right)", "entrance excess, first street\n(% too slow in the model)",
                              "uniform factor left\n(model / measured, geometric mean)", "bits saved by the formula\nover one factor (per point)"),
                             (1, 0, 1, 0)):
    ax.axhline(ref, color=INK2, lw=1, ls=":")
    ax.plot(N, val, color=S1, lw=2, zorder=2)
    for n_, v_, z in zip(N, val, z0):
        ax.scatter([n_], [v_], s=60, color=S2 if z else S1, ec=SURF, lw=2, zorder=3)
    ax.set_title(ttl, fontsize=9.5, loc="left", color=INK); ax.set_xticks([10, 16, 20, 24, 32]); ax.set_xlabel("cells per block width")
    ax.grid(axis="y", color=GRID, lw=0.5)
axs[0].scatter([], [], color=S1, s=40, label="log-fit ground z0"); axs[0].scatter([], [], color=S2, s=40, label="tunnel-rule ground z0")
axs[0].legend(frameon=False, fontsize=8, loc="lower right")
fig.suptitle("LBM layout 1D: the shape of the error fades with resolution; one uniform factor (about 1.3x) does not",
             x=0.02, ha="left", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(FIG / "04_resolution.png", dpi=170); plt.close(fig)


# 5. Is there structure at all? neighbours vs shuffles
arr = np.load(RES / "aij_c_arrays.npz")
fig, axs = plt.subplots(1, 2, figsize=(9, 3))
for ax, nm, key in zip(axs, ("LBM-LES D/10, 360 points", "RANS k-eps D/10, 360 points"), ("lbm", "rans")):
    nul = arr[f"null_{key}"]; obs = main[f"step1_{key}"]["bits_per_reading"]
    ax.hist(nul, bins=25, color=GRID, ec=SURF)
    ax.axvline(obs, color=S1, lw=2); ax.text(obs, ax.get_ylim()[1] * 0.9, f" real layout\n {obs:.2f} b, p = {main[f'step1_{key}']['p']:.3f}", color=INK, fontsize=8.5)
    ax.set_title(nm, loc="left", fontsize=10, color=INK); ax.set_xlabel("bits per point, each guessed from its neighbours")
axs[0].set_ylabel("199 shuffles")
fig.suptitle("Do errors at neighbouring points predict each other? LBM: clearly. RANS pooled: borderline (protocol mark p < 0.01)", x=0.02, ha="left", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(FIG / "05_structure_test.png", dpi=170); plt.close(fig)
print("figures in", FIG)
