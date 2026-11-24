# Protocol 7: a stronger physics compressor (written 24 Nov 2026, before any run of the new code)

Ask (Kyle, 24 Nov): keep building the compressor. Make it more general and harder to fool: a bigger term library
(atmospheric stability, convective and tropical terms), a way to find term forms that are not in the library, and
guards against its known failure modes (look-alike terms, small samples, confounding by group, noise that no term
can explain), with honest uncertainty on what it reports. Internal research engine; nothing here goes on the website.

**Seen before this was written:** `physlib/` as merged (main 931685c) and its PROTOCOL-4 results
(`results/physlib/synthetic.json`: 10 of 12 planted laws pass, 0 false finds on noise). The WindSim thread's
protocol 6 (gridzip-lab branch `claude/project-thread-fnbdyl`, 75980a8) has already fitted cold-pool, pressure-jump
and dryness terms to 903 real thunderstorm gusts and found nothing beyond station offsets; I have read its summary
note, not its fit. Those data are therefore spent: part D runs them only as a post-hoc check of the engine.

## What is built (physlib v2)
1. **Library growth** (`physlib/terms.py`). Wind gains convective and stability terms (cold-pool density-current
   speed √(g·h·ΔT/T)/U, pressure-jump speed √(2Δp/ρ)/U, downdraft speed √(2·DCAPE)/U, CAPE, low-level lapse rate,
   dry air aloft, rain loading, Monin-Obukhov ψ_m(z/L)/ln(z/z0), bulk Richardson number, storm motion, convective
   velocity w*) and tropical-cyclone terms (Kaplan-DeMaria style decay after landfall e^(−α·t) with α a family,
   radius-to-maximum-wind shape, motion asymmetry cos(θ−θ_motion)·U_storm/U, sea-to-land roughness). Heat gains
   mixed-layer depth, advection and soil moisture. Forms are standard textbook ones written from memory; none of the
   papers they come from was opened for this work, so no constant is quoted from them: every constant is free.
2. **Discovery beyond the library** (`physlib/grammar.py`). A small grammar builds candidate terms from a case's raw
   inputs: one input under one operation (identity, log, square root, square, inverse, decay e^(−x/L), hinge
   max(0, x − q), step [x > q], with L and q at the 25/50/75 % quantiles of x), or a product of two such pieces.
   Each candidate costs its own description in bits (operation, input, constant choice, product flag), so a
   grammar term pays more to be named than a library term and the library's physics names win ties. One flag bit
   per term says library or grammar. A prescreen (correlation with the current remainder) keeps the beam affordable.
3. **Guards** (`physlib/guards.py`):
   - *Look-alikes:* swap each winning term for every candidate that resembles it (|r| > 0.8) and re-score exactly; any
     swap within 3.3 bits (odds under 10:1) is listed as "can't rule out". The search's scored formulas give
     approximate posterior odds (weight 2^−bits) and an inclusion probability per term.
   - *Group confounding:* a winning term that is nearly constant within groups (stations, places) is compared with
     plain group offsets; if the offsets code as well (within 3.3 bits), the term is flagged "group-level, can't be
     told from siting".
   - *Small samples:* a bootstrap (by group when groups exist) reports how often each term is chosen and a 90 %
     interval for its coefficient.
   - *Is the gain real:* a permutation null (y shuffled within blocks, the whole search re-run) gives the bits a
     search of this size finds in noise.
   - *Irreducible noise and unnamed structure:* a neighbour test in input space (each reading's k nearest neighbours
     over the case's standardised inputs predict its remainder; 99 shuffles) says whether structure is left that no
     term named; the remainder sd is reported next to the replicate sd where a case has one.
   - *Heavy tails:* residuals may be coded under a two-scale Gaussian mixture (5, 10 or 20 % of readings at 4× the
     sd) when that is shorter, at 2 bits for the choice, so a few wild readings can't buy a term.

## Tests with known answers (synthetic; 20 seeds each unless stated; noise iid Gaussian unless stated)
Planted signal size: each planted law's signal sd equals the noise sd (signal-to-noise 1), unless stated.

**A. The bigger library still finds what is planted.** Full wind library (old + convective + tropical).
On synthetic thunderstorm rows (n = 400, inputs drawn with physical correlations, e.g. Δp rises with ΔT):
A1 cold pool; A2 cold pool + pressure jump; A3 MO stability; A4 downdraft + station offset. On synthetic
landfall rows (n = 400): A5 decay after landfall (α = 0.095 per hour); A6 motion asymmetry.
- **Rule A:** each is named exactly (a term of the same declared group counts) in ≥ 80 % of seeds and explained
  (R² of fitted against planted signal ≥ 0.9) in ≥ 80 %.
- **Rule A-reg:** the 12 PROTOCOL-4 laws, rerun with the v2 engine (wind on the full library), each named at a rate no
  more than 10 points below its PROTOCOL-4 rate.

**B. Discovery: forms the library does not have** (grammar on, the planted form removed from the library).
B1 an interaction hc · lp (contrast that changes with neighbour height); B2 a gated decay corner · e^(−fetch/L) with
L at the median of fetch; B3 a heat moisture-under-cloud product T_dew · cloud; B4 a hinge max(0, z − median z).
B5 off-grid: e^(−fetch/3.1), a length the grammar does not offer. B6 control: the PROTOCOL-4 law "contrast +
entrance decay", which the library has, with the grammar on.
- **Rule B:** B1 to B4: the winner holds a term whose column has |r| ≥ 0.95 with the planted one in ≥ 70 % of seeds,
  and explained in ≥ 80 %. B5: R² ≥ 0.8 in ≥ 70 %. B6: the library terms named exactly in ≥ 80 % (the grammar must
  not displace the physics names).

**C. Guards.**
- C1 noise, grammar on: wind, heat, flood, ground on iid noise (50 seeds each); heat and ground on AR(1) ρ 0.8 (30);
  wind and flood on Student-t ν = 3 (50); wind at n = 40 (50). **Rule:** any term named in ≤ 5 % of runs, each set.
- C2 look-alikes: plant term A; the library holds A and a twin B with r(A, B) = 0.97. **Rule:** the report says it
  can't tell (B listed in "can't rule out", or A's inclusion probability < 0.9) in ≥ 80 % of seeds. Separable pair
  r = 0.6: A named with inclusion ≥ 0.9 in ≥ 80 %.
- C3 group confounding: 6 stations with siting offsets (sd = noise sd), plus a station-level attribute term with no
  effect. **Rule:** whenever the attribute term is named, the guard flags it, in ≥ 90 % of those seeds; unflagged
  false physics in ≤ 5 % of all seeds.
- C4 small samples: A1 at n = 30, 60, 120, grammar on. **Rule:** false terms per run ≤ 0.2 at every n. At n = 120,
  the bootstrap 90 % interval (30 resamples) covers the planted coefficient in ≥ 80 % of seeds where it is named.
  Detection rate at each n is reported (no rule).
- C5 unnamed structure: (a) law planted and named: the neighbour test reports structure left (p ≤ 0.05) in ≤ 10 %;
  (b) a planted form outside the grammar (sin(3·lp) · cos(fetch)): structure left flagged in ≥ 70 %; (c) noise
  only: the reported remainder sd within 10 % of the true noise sd in ≥ 90 %.

If running time forces fewer seeds, that is written as an amendment before the run, never after.

## D. The real cases again (exploratory; after A to C are scored; nothing here is a test)
W1 AIJ Case C (LBM D/10 and RANS, `results/physics/aij_c_arrays.npz`) with the grammar and guards on: does
discovery name any of the ~0.5 bit per point of neighbour structure FINDINGS-4 left unnamed? W2 Layer 1 storm gusts
(`results/physics/layer1.json`) with the group guard. The WindSim protocol 6 thunderstorm table (903 gusts,
spent) through the new search, post hoc. Each is reported next to its FINDINGS-4 (or FINDINGS-6) result; agreement
and disagreement both go in FINDINGS-7. No gate, rule or seed count above changes after any result is seen.

## Amendment A1 (24 Nov 2026, after development smoke runs on seeds 90000 and up, before any protocol seed is run)
Seen: smoke runs of C2 and C3 (4 seeds each, seeds ≥ 90000), and while debugging the search, one run each of the
B6 law on seeds 0 and 9000-9004 (development; seed 0 is also a B6 protocol seed, so B6 seed 0 is reported as seen).
Changes, each from a failure of the code or of a test's premise, none from a protocol result:
1. *Search:* approx bits undercount parameters, so with thousands of grammar terms the exact finalists were all
   5-term formulas and a shorter true formula was never scored exactly. With the grammar on, the best 3 formulas of
   each size also go to exact scoring. A non-finite score is now +∞ (it scrambled the beam's sort).
2. *Grammar:* a grammar shape that the library already has (|r| > 0.995 with a plain library term or with a family
   member anywhere on its grid) is dropped, so a family's free parameter isn't undercut by the grammar's fixed
   constant. Grammar columns are scaled to unit sd.
3. *Inclusion is also reported per input* (a grammar shape of x counts for x, as a library term reading x does):
   shapes of one input were splitting its mass. C2's "A's inclusion" means A's input's inclusion.
4. *C2's twin correlation was wrong for its purpose.* At n = 360 and signal-to-noise 1, r = 0.97 twins cost about
   n/2 · log2(1 + 1 − r²) ≈ 15 bits to swap, so the data can tell them apart and the guard rightly says so. The rule
   now uses r = 0.998 (expected swap cost under 1 bit), with the same pass rule; r = 0.97 is reported, no rule.
5. *C3's group guard as written misses a one-parameter station-level term that beats five offsets on parameter
   count.* A group-level term is now also flagged when it doesn't pay on the group means (k means of the remainder
   without it, coded with and without it: k/2 · log2 of the SSE ratio against naming plus one parameter). "The
   attribute is named" includes any grammar shape of it.
6. The neighbour test (C5) takes the best of the full input space and of every pair of inputs, the same maximum in
   each shuffle, so a pattern in two inputs among twenty isn't washed out.
