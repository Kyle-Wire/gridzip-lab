"""PROTOCOL-8 figures (grid load). python -m physlib.figures8  -> physics/figures/gridload/*.png"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results" / "gridload"; F = ROOT / "physics" / "figures" / "gridload"
SURF, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
GOOD, CRIT = "#0ca30c", "#d03b3b"
NAMES = {"TOTAL": "All New England", "ME": "Maine", "NH": "New Hampshire", "VT": "Vermont", "CT": "Connecticut", "RI": "Rhode Island",
         "SEMASS": "SE Mass.", "WCMASS": "W/Central Mass.", "NEMASSBOST": "Boston / NE Mass."}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF})


def _title(fig, t, sub):
    fig.text(0.02, 0.965, t, fontsize=16, weight="bold", color=INK, va="top")
    fig.text(0.02, 0.915, sub, fontsize=11, color=INK2, va="top")


def head_to_head():
    C = json.loads((R / "C_sealed.json").read_text()); zs = list(C)
    labs = [("M0_vanilla", "Standard utility model"), ("M1_vanilla_plus_compressor", "Standard + 5 terms the compressor found"),
            ("M2_expert", "Expert hand-built model"), ("M3_gbm", "Machine learning (boosted trees)")]
    fig, axs = plt.subplots(1, 2, figsize=(14, 7.2), gridspec_kw={"width_ratios": [1.7, 1]})
    _title(fig, "Adding what the compressor found makes the standard model better everywhere",
           "Average hourly forecast error, 2014-2016 (years no model or search had seen). ISO New England, actual weather. Lower is better.")
    ax = axs[0]; w = 0.2; x = np.arange(len(zs))
    for k, (m, lab) in enumerate(labs):
        v = [C[z][m]["mape"] for z in zs]
        ax.bar(x + (k - 1.5) * w, v, w * 0.92, color=S[k], label=lab, zorder=3)
    ax.set_xticks(x); ax.set_xticklabels([NAMES[z] for z in zs], rotation=30, ha="right"); ax.set_ylabel("error, % of load (MAPE)")
    ax.yaxis.grid(True, color=GRID, zorder=0); ax.set_ylim(0, 6)
    ax.legend(loc="upper left", frameon=False, fontsize=10, ncol=2)
    ax = axs[1]; t = C["TOTAL"]
    for k, (m, lab) in enumerate(labs):
        for j, key in enumerate(("mape", "peak_mape")):
            v = t[m][key]; ax.barh(j * 5 + (3 - k), v, 0.85, color=S[k], zorder=3)
            ax.text(v + 0.05, j * 5 + (3 - k), f"{v:.2f} %", va="center", fontsize=10, color=INK)
    ax.set_yticks([1.5, 6.5]); ax.set_yticklabels(["every hour", "peak hours\n(top 1 %)"]); ax.set_xlim(0, 5)
    ax.set_title("All New England", loc="left", color=INK, fontsize=12); ax.xaxis.grid(True, color=GRID, zorder=0); ax.set_xlabel("error, %")
    fig.text(0.02, 0.015, "Each bar group: standard (blue), + compressor terms (orange), expert (aqua), machine learning (yellow). "
             "The compressor's 5 terms close 74 % of the gap the expert's ~370 extra inputs open.", fontsize=9.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 0.88)); fig.savefig(F / "01_head_to_head.png", dpi=150); plt.close(fig)


GROUPS = [("recency", "Yesterday's heat\n(lagged temperature)"), ("humidity", "Humidity /\ndew point"),
          ("sun", "Sun & daylight\n(solar, lighting, hour)"), ("season", "Season shape\n(day of year)"), ("temp", "This hour's\ntemperature (shape)"),
          ("trend", "Growth trend"), ("holiday", "Holidays")]
LIB = {"recency_cool": "recency", "recency_heat": "recency", "prev_day_heat": "recency", "lag_gap": "recency", "heat_wave": "recency",
       "humid_cool": "humidity", "dewpoint": "humidity", "btm_solar": "sun", "lighting": "sun", "trend_curve": "trend",
       "holiday": "holiday", "holiday_adj": "holiday", "daily_shape": "sun"}
TOK = {"Tavg24": "recency", "Tmax_prev": "recency", "Td": "humidity", "sun": "sun", "hour": "sun", "doy": "season", "trend": "trend", "T": "temp"}


def _groups_of(term):
    import re
    if term in LIB: return {LIB[term]}
    return {TOK[t] for t in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", term[2:]) if t in TOK}


def what_it_named(fname="B_explore.json", out="02_what_it_named.png", title=None, model=None):
    B = json.loads((R / fname).read_text()); zs = list(B)
    M = np.zeros((len(GROUPS), len(zs)))
    for j, z in enumerate(zs):
        win = B[z][model]["winner"] if model else B[z]["winner"]
        for t in win:
            for g in _groups_of(t): M[[k for k, _ in GROUPS].index(g), j] += 1
    fig, ax = plt.subplots(figsize=(13, 6.4))
    _title(fig, title or "What the standard model misses, in plain physics",
           "Terms the compressor named in each region's errors, 2011-12 (it searched 14 named terms and ~900 shapes it built itself)")
    for i in range(len(GROUPS)):
        for j in range(len(zs)):
            v = M[i, j]
            ax.add_patch(plt.Rectangle((j + 0.05, i + 0.05), 0.9, 0.9, color=S[0] if v else "#efeeea", alpha=min(1, 0.35 + 0.3 * v) if v else 1))
            if v: ax.text(j + 0.5, i + 0.5, f"{int(v)}", ha="center", va="center", color="white", fontsize=12, weight="bold")
    ax.set_xlim(0, len(zs)); ax.set_ylim(len(GROUPS), 0)
    ax.set_xticks(np.arange(len(zs)) + 0.5); ax.set_xticklabels([NAMES[z] for z in zs], rotation=30, ha="right")
    ax.set_yticks(np.arange(len(GROUPS)) + 0.5); ax.set_yticklabels([g[1] for g in GROUPS]); ax.tick_params(length=0)
    for s_ in ax.spines.values(): s_.set_visible(False)
    fig.text(0.02, 0.015, "Number = how many of the region's winning terms involve that input. Holidays were never named on the real data, "
             "though it found a planted holiday dip in 10 of 10 tests.", fontsize=9.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 0.88)); fig.savefig(F / out, dpi=150); plt.close(fig)


def trust():
    A = json.loads((R / "A_planted.json").read_text())
    lab = {"P1_recency": "Yesterday's heat", "P2_humidity": "Humid heat", "P3_holiday": "Holiday dip", "P4_solar": "Rooftop solar",
           "P5_combo": "Three at once", "P6_grammar": "A shape not in its library", "N_null": "Pure noise (must find nothing)"}
    fig, ax = plt.subplots(figsize=(11, 5.4))
    _title(fig, "Can it be trusted? Hidden effects found, and silence on noise",
           "Real New England weather 2011-12, made-up demand errors with a known answer. Rules written before any run.")
    ks = list(A); y = np.arange(len(ks))[::-1]
    for yi, k in zip(y, ks):
        v = A[k]["named"] / A[k]["of"] * 100
        ax.barh(yi, v, 0.6, color=GOOD if v >= 80 else CRIT, zorder=3)
        ax.text(v + 1, yi, f"{A[k]['named']} of {A[k]['of']}", va="center", color=INK, fontsize=11)
    ax.axvline(80, color=MUTED, ls="--", lw=1); ax.text(80, len(ks) - 0.4, "pass line 80 % (noise: 95 %)", color=MUTED, fontsize=9, ha="center")
    ax.set_yticks(y); ax.set_yticklabels([lab[k] for k in ks]); ax.set_xlim(0, 115); ax.set_xlabel("% of runs answered correctly")
    ax.xaxis.grid(True, color=GRID, zorder=0)
    fig.tight_layout(rect=(0, 0.02, 1, 0.86)); fig.savefig(F / "00_trust.png", dpi=150); plt.close(fig)


def audit_fig():
    E_ = json.loads((R / "E_fresh.json").read_text()); zs = list(E_)
    fig, axs = plt.subplots(1, 2, figsize=(14, 6.6), sharey=True)
    _title(fig, "Auditing strong models: does what it finds still help?",
           "Change in hourly error after adding the compressor's terms. Fresh test: Jan-Apr 2017, never read before. Below zero = better.")
    for ax, (m, lab) in zip(axs, (("M2", "Expert hand-built model"), ("M3", "Machine learning (boosted trees)"))):
        d = [E_[z]["jan_apr_2017_fresh"][m + "_plus"] - E_[z]["jan_apr_2017_fresh"][m] for z in zs]
        x = np.arange(len(zs)); ax.bar(x, d, 0.6, color=[GOOD if v < 0 else CRIT for v in d], zorder=3)
        for xi, v in zip(x, d): ax.text(xi, v + (0.01 if v >= 0 else -0.01), f"{v:+.2f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=9, color=INK)
        ax.axhline(0, color=INK2, lw=1); ax.set_xticks(x); ax.set_xticklabels([NAMES[z] for z in zs], rotation=30, ha="right")
        ax.set_title(lab, loc="left", color=INK, fontsize=12); ax.yaxis.grid(True, color=GRID, zorder=0)
    axs[0].set_ylabel("change in error, percentage points")
    fig.text(0.02, 0.015, "The two worst results (SE Mass. expert, Boston machine learning) are the two audits that failed the compressor's own "
             "check on 2013: its guard flagged them before the test. 2,876 hours per region; winter and spring only.", fontsize=9.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 0.86)); fig.savefig(F / "03_audit_fresh.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    F.mkdir(parents=True, exist_ok=True); trust(); head_to_head(); what_it_named()
    if (R / "E_fresh.json").exists(): audit_fig()
