"""InSAR rate-field terms for the physics compressor (PROTOCOL-13), reusable on any InSAR stack.

Target: a per-point LOS rate on ground whose true rate is known (stable ground, or GNSS-tied), or a rate error.
ctx keys (all (n,) arrays): coherence (temporal coherence 0-1), ps_share (share of PS pixels in the cell),
amp_disp (amplitude dispersion), outside_comp (share of interferograms outside the unwrapping component),
closure_amb (ln(1 + closure triplets with a non-zero ambiguity)), jumps (date-to-date steps > lambda/4),
height (DEM height, m), incidence (deg), ref_dist (km to the reference), east, north (km), annual_amp (mm).
Written fresh for SPARC / gridzip-lab.
"""
from __future__ import annotations

import numpy as np

from .engine import Term


def _c(a): a = np.asarray(a, np.float64); return a - a.mean()


def _t(name, formula, meaning, key, group=None):
    return Term(name, "ground_m2", formula, meaning, (key,), (lambda k: lambda c: _c(c[k]))(key), group=group)


INSAR_RATE = [
    _t("coherence", "b · temporal coherence", "noisier (low-coherence) cells carry a rate bias", "coherence", "quality"),
    _t("ps_share", "b · PS share of the cell", "distributed scatterers (phase-linked) and PS disagree in rate", "ps_share"),
    _t("amp_disp", "b · amplitude dispersion", "unstable scatterers carry a rate bias", "amp_disp"),
    _t("outside_comp", "b · share of ifgrams outside snaphu's component", "unwrapping islands", "outside_comp", "unwrap"),
    _t("closure_amb", "b · ln(1 + closure triplets with a non-zero ambiguity)", "unwrapping errors (closure)", "closure_amb", "unwrap"),
    _t("jumps", "b · number of >λ/4 date-to-date jumps", "cycle slips left in the series", "jumps", "unwrap"),
    _t("height", "b · DEM height", "stratified troposphere or DEM error left in the rate", "height"),
    _t("incidence", "b · incidence angle", "a look-geometry term (range ramp)", "incidence"),
    _t("ref_dist", "b · distance to the reference cell", "error grows away from the reference", "ref_dist"),
    Term("ramp", "ground_m2", "b1 · east + b2 · north", "a long-wavelength tilt (orbit, ionosphere, troposphere)", ("east", "north"),
         lambda c: np.column_stack([_c(c["east"]), _c(c["north"])]), group="ramp"),
    Term("ramp2", "ground_m2", "quadratic in east, north", "a long-wavelength bowl (troposphere, ionosphere)", ("east", "north"),
         lambda c: np.column_stack([_c(c["east"]), _c(c["north"]), _c(c["east"] ** 2), _c(c["north"] ** 2), _c(c["east"] * c["north"])]), group="ramp"),
    _t("annual_amp", "b · annual amplitude", "seasonal ground (soil moisture) leaks into the rate", "annual_amp"),
]

# Fading-signal term (PROTOCOL-13 A5): ctx key "sawtooth" = the point's amplitude on the scene's per-ministack ramp
# pattern, fitted after the point's own [1, t, sin, cos] fit (see sawtooth_amplitude). Not in INSAR_RATE by default,
# because it needs the series.
SAWTOOTH = _t("sawtooth", "b · a_saw,  a_saw = point's amplitude on the scene ministack-ramp pattern",
              "fading-signal bias of phase-linked short-baseline stacks: a rate bias that scales with the sawtooth", "sawtooth", "fading")


def sawtooth_amplitude(Y, t, ministack, ref_rows=None):
    """Y: points x dates (mm), t: years, ministack: (dates,) ints. Returns (amplitude per point, pattern P)."""
    X = np.column_stack([np.ones_like(t), t, np.sin(2 * np.pi * t), np.cos(2 * np.pi * t)])
    beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None); R = Y.T - X @ beta            # dates x points
    med = np.median(R if ref_rows is None else R[:, ref_rows], axis=1)
    P = np.zeros_like(t)
    for k in np.unique(ministack):
        m = ministack == k
        if m.sum() < 2: continue
        A = np.column_stack([np.ones(m.sum()), t[m]]); b, *_ = np.linalg.lstsq(A, med[m], rcond=None)
        P[m] = t[m] * b[1] - (t[m] * b[1]).mean()
    P = P - X @ np.linalg.lstsq(X, P, rcond=None)[0]
    amp = (P @ R) / (P @ P)
    return amp, P


