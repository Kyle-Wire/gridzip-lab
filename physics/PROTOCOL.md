# Finding missing physics by code length: protocol (written 22 Nov 2026, before any search on real residuals)

## The question
Kyle asked whether the compressor and codelen can pull formulas and physical laws out of data, so we can see
what our models are missing. The honest version of the claim is narrower than "discover physics":

> If a model's errors against measurement can be stored in fewer bits once a short formula is added, and the
> saving pays for writing the formula down, the errors contain structure, and the formula says what shape it has.

This is the minimum description length (MDL) principle. It is the scoring rule behind symbolic regression and
SINDy-style sparse term selection; codelen makes every bit checkable with a real arithmetic coder.

## Case
WindSim on AIJ Case C (3 x 3 cubes, D = 0.2 m, 0 deg wind, scalar speed at z = 0.02 m, 120 points per layout,
centre block 0D / 1D / 2D). Per-point files committed in WindSim at `3c780ff`, `docs/validation/`:
- LBM-LES blind run at D/10, three layouts (`lbm_case_c_d10_blind/`)
- LBM-LES layout 1D at D/16, D/16 tunnel z0, D/20, D/24, D/32 tunnel z0 (`lbm_case_c_1D_resolution/`)
- RANS k-eps at D/10, three layouts, |U| and sqrt(|U|^2 + 2k) (`rans_case_c_d10/`), and D/20 layout 1D

**What had been seen before writing this:** WindSim's `docs/validation.md` summaries (hit rates, biases, mean
speed by street type for LBM 1D). No per-point residual had been fitted or plotted by this work.

## Target and grid
r = ln(measured / predicted) for each point: the factor the model is missing. Coded on a grid of 0.01 (1 %),
range -2.5 to 2.5, escape eps 1e-3 (codelen's pre-registered value). Bits compare only on this grid.

## Feature library (fixed now; geometry and model output only, never the measurement)
Wind is +x; the front row spans x = -0.5 to -0.3 m.
- `A`, `X`, `I`: point in an along-wind street, a cross street, an intersection (indicator).
- `mouth`: along-wind street point on the line of a block face (x = +-0.1, +-0.3 m).
- `xi`: streamwise distance from the array's front face, in D.
- `exp_xi(L)`: exp(-xi / L), L fitted (the decay of an entrance excess along a street).
- `lee`: in a cross street, signed position from the upstream block's leeward face (-0.5) to the next block's
  windward face (+0.5).
- `off`: lateral offset from the street centreline, in street half-widths (0 or 0.5 or 1).
- `hc`: centre block height in D (0, 1, 2); `near_c`: point touching the centre block's footprint (within 0.1 m).
- `lp`: ln(predicted / 2.43 m/s), the model's own speed-up over the approach flow at z = 0.02 m. A term in `lp`
  is the law measured = a * predicted^b, i.e. the model exaggerating (b < 1) or damping (b > 1) contrasts.
- RANS only, `kfrac(alpha)`: 0.5 ln(1 + alpha (W^2/U^2 - 1)), W = sqrt(U^2 + 2k). This is measured =
  sqrt(U^2 + alpha * 2k): how much of the turbulence a thermistor sees. Theory for isotropic turbulence and a
  strong mean flow gives alpha near 2/3 (the two cross-flow variances add to |u|).
- Products of an indicator with `xi`, `exp_xi`, `hc`, `lp`.

## Steps and readings (fixed now)
0. **Sanity, synthetic, before real data.** On Case C's own geometry and the LBM D/10 predictions:
   (a) a planted law r = -0.1 - 0.35 A exp(-xi/1.5) + noise must be ranked first, with L within 30 % and the
   coefficient within 20 %; (b) pure noise around a constant must rank "bias only" first; (c) RANS data with
   measured = sqrt(U^2 + (2/3) 2k) x noise must give alpha within 0.15 of 2/3. If any fails, the tool is fixed
   and the sanity re-run before step 1.
1. **Is there structure at all?** (model-free) Code r with each point predicted from its already-coded spatial
   neighbours (codelen `ContextMean`, nearest earlier points within 0.075 m, random order, seed 0) and compare
   with 199 shuffles of r over the same positions. Structure is called present when p < 0.01.
   **Noise floor:** the measurements are mirror-symmetric across y = 0 in a symmetric flow, so
   sd_noise = sd(ln m(x, y) - ln m(x, -y)) / sqrt(2). The entropy of N(0, sd_noise) on the grid is the floor; bits
   above it are structure still to find (or asymmetry of the tunnel, which the floor also contains).
2. **Search.** Every formula with up to 3 terms from the library plus an intercept, least-squares fit, each
   parameter's precision chosen by MDL (codelen.mdl conventions). Ranked by total bits (formula + parameters +
   residuals) on the pooled data of one engine. The winner is verified with the arithmetic coder.
3. **Hold-out.** For the top 5 formulas: fit on two layouts, code the third (all three ways). A formula
   **transfers** only if it beats the bias-only fit on every held-out layout.
4. **Resolution (LBM, layout 1D).** Fit the step-2 winner's terms at each resolution. A term whose coefficient
   shrinks by more than half from D/10 to D/32 is labelled *resolution*; otherwise *persistent*.
5. **Engines.** Same search on RANS |U|. A term found in both with the same sign points at something the
   engines share (inflow, measurement, set-up); opposite signs point at each engine's own physics.

## What this can and cannot claim
- It finds the *shape* of what is missing, in the vocabulary we give it. It cannot name a mechanism the
  vocabulary has no word for; a large gap to the noise floor after the search says that is the case.
- Case C has been seen by WindSim's own diagnosis, and here it is fitted. Nothing found here is validated: a
  formula becomes a model change only through WindSim's own protocol and a held-out case (B or D).
- An empirical correction fitted to Case C must not be applied to reports. Its use is to point at the physics.

## Amendment A1 (22 Nov 2026, before any real-data search; nothing but the planted-law unit tests had run)
`kfrac(alpha)` enters with its coefficient fixed at 1, as the law measured = sqrt(U^2 + alpha 2k) says; only alpha
is fitted. A free slope on it would trade against alpha (slope x alpha is what a small k sees), so alpha could not
be read. Reported parameters are the least-squares values; the MDL-rounded values the code length uses are given
beside them.

## Amendment A2 (22 Nov 2026, after the first run stopped on an error, before any search result was read)
The grid range is widened from -2.5..2.5 to -4..4. RANS |U| falls to 0.04-0.06 m/s at four points (recirculation
centres) where the thermistors read 1.2-1.6 m/s, so r reaches 3.75 and the coder refused it. Nothing else changes;
every number below is on the -4..4 grid.
