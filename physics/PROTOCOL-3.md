# Protocol 3: what WindSim's storm-gust chain (Layer 1) is missing (written 22 Nov 2026, before any fit)

Same method as [PROTOCOL.md](PROTOCOL.md) and [PROTOCOL-2.md](PROTOCOL-2.md): an MDL formula search over a named vocabulary
(`codelen/lawsearch.py`, `physics/runner.py`), hold-out folds, and a model-free neighbour test.

**What has been seen before this was written.** The per-place shares within 25 %, medians and verdicts in WindSim's
`docs/validation.md` (run 4, A3, T1, T2, F2) and the PR 10 description. No per-pair error has been fitted, plotted or
grouped by any feature below. Station coordinates and elevations were fetched (first lines of each station's ISD file,
AWS `noaa-global-hourly-pds`) and not yet joined to errors.

## Data
WindSim PR 10 (branch `claude/project-thread-fnbdyl`, d36598b) check files, every station-storm pair:
- **after F2** (63dd5c5, tropical storms separated): `docs/validation/layer1/f2/{cal,blind}/*.json`, 10 calibration and
  4 blind places (Boston cannot run).
- **before F2** (tropical storms inside "non-thunderstorm"): run 4 (`run4/checks/*.json`, 1297da7) plus A3
  (`a3_blind/*.json`, 622f9e6, from branch `claude/project-thread-er7q2w`): the same 10 calibration places.

The quantity is **y = ln(observed / predicted)** for each pair (the chain's own error, log form so that ±25 % is
symmetric). Grid: 0.005 in y, −3..3.

## Vocabulary (every term centred; fixed now)
Pair or station level:
- `lnz0`: ln of the target's upwind z0 for the storm direction; `kzt`.
- `near_km`: distance to the nearest station used; `n_nb`: number of stations used.
- `dz_elev`: the target's elevation minus the inverse-distance-squared mean elevation of the other stations in that
  storm (km). Terrain between stations.
- `elev`: the target's elevation (km).
- `dir_c(a)`: cos(direction − a), a fitted over 0..345° (one parameter). `offshore_ca`: in California only,
  1 if the direction is from 0° to 135° (Santa Ana and other offshore winds), else 0.
- `lnpred`: ln of the predicted gust (does the miss depend on gust size).
Storm level:
- `msin`, `mcos`: season (month harmonics). `hsin`, `hcos`: local solar hour harmonics (from the station's longitude).
- `lnreg`: ln of the storm's regional gust over its place's median regional gust (storm size).
- `yr`: storm year − 2012 (ASOS sensor changes).
Context, leave-one-out (computed without the pair itself, so they cannot copy it):
- `storm_loo`: mean y of the other stations in the same storm.
- `stn_loo`: mean y of the same station in the place's other storms of the same type.
Place:
- `ca`: 1 at a California place (Los Angeles, Riverside, San Francisco).
- `desert`: 1 at Phoenix. Used only in the spread test (Q1b) and as a labelled term, never as a "law".

## Questions and readings, fixed now
**Q1, thunderstorm gusts (every thunderstorm pair, after F2).**
- (a) Formula search, max 3 terms. Folds: leave out each place; and calibration → blind. A term counts only if it
  holds on 80 % of place folds and on calibration → blind.
- (b) Spread: code the pairs with one pooled sd against one sd per place (each extra sd paid at MDL cost), and
  against pooled + a Phoenix sd. If per-place spread pays, the ±30 % range should vary by place, and which places
  are wide is reported.
- Neighbour tests (99 shuffles): context = the same storm's other stations; and separately the same station's other
  storms. If neither beats its null, what is left is point noise, and a range (T2) is the right form. If the storm
  context pays, the miss is storm-wide (the regional transfer is off), which a better regional scale could fix. If the
  station context pays, it is the station (exposure, sensor), which the Verkaik gust-z0 check could fix.
- Phoenix and Atlanta are reported apart: which side the misses are on, and whether any term above explains them.

**Q2, ordinary storms (non-thunderstorm pairs, after F2).**
- Formula search as Q1a, with `ca` and `offshore_ca` in the vocabulary. **Reading:** a California pattern is named
  by a term only if (i) it beats one correction on all pairs, (ii) it holds when each California place is left out in
  turn (fit on the other 13 places), and (iii) a physical term (`dz_elev`, `offshore_ca`, `dir_c`, season) does at least
  as well as the bare `ca` label. If only `ca` pays, the pattern is real but unnamed. If nothing pays, the California
  misses are wider scatter, and Q1b's spread test is run on ordinary storms too (per-place sd, and `ca` sd).
- Neighbour tests as Q1.

**Q3, the hurricane split (F2).** At the 10 calibration places:
- bits per pair of "one correction" and of the best formula for ordinary storms before F2 against after F2 (on the
  pairs both sets share and on all pairs), and the neighbour tests before and after;
- in the before set, a `tropical` term (1 if the storm time matches a storm F2 typed tropical at that place): if the
  split named real structure, this term pays in the before set.
- Reading: F2 helped by this measure if the after set's ordinary storms code in fewer bits per pair than the before
  set's, and the `tropical` term pays before.

Nothing here changes the chain or a gate. Findings go to the WindSim thread as leads for its own protocols.

## Amendment A1 (22 Nov 2026, after the first run, before any write-up)
The first run's winner in every set used `storm_loo` with a coefficient of −2 to −2.7, and the same-storm neighbour
test coded *worse* than its shuffles (p = 1.0). That is the check's own construction, not physics: each station is
predicted from the other stations of the same storm, so a station that reads high pushes its neighbours' predictions
up and their errors down. Errors within a storm are anti-correlated by design. `storm_loo` is therefore dropped from
the vocabulary and the search is rerun unchanged otherwise. The same-storm neighbour test is kept in the record and is
read only as a sign of this coupling; storm-wide misses cannot be tested from these files. Both runs are reported.
