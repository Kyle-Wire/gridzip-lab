# codelen: score forecasts by code length

The score for a forecast is the number of bits needed to store each real reading, quantised to a stated resolution, under the forecast's predicted distribution. Fewer bits means a better forecast. Unlike RMSE, this score also rewards honest uncertainty: an overconfident sd is punished hard, and a too-wide one is punished gently.

Every number can be checked. With `verify=True`, the readings are arithmetic-coded with integer tables, decoded again, and the round trip is asserted. The coded length comes within a few bits of the ideal, Σ −log2 p.

The code is written fresh, with no lpaq1, cmix or other GPL code. It needs numpy; scipy is optional and makes it faster. It runs on CPU only.

## Conventions
These follow the pre-registered SPARC field-model scores (sparc-transformers `docs/design/field-model.md`, "Scores for real tests", 21 Nov 2026). They agree with the analytic reference `sparc_tf/bits.py` to within 1e-5 bits per reading, and the tests check this.

- **Bins:** a reading y is snapped to r = round(y / res) · res and coded as the bin [r − res/2, r + res/2], for r from `lo` to `hi`. Every output states res, lo, hi and eps. Bits only compare at the same resolution: making `res` 10 times finer adds 3.32 bits to every reading of a smooth forecast.
- **Escape:** each bin gets (1 − eps) · model mass + eps · res / (hi − lo). The model mass is the exact bin mass and is not renormalised: mass the forecast puts outside [lo, hi] is lost, just as in the analytic score. The worst case per reading is −log2(eps · res / (hi − lo)), which is 20.2 bits at 0.1 K, eps 1e-3 over −60..60 °C.
- **Gaussian:** a difference of CDFs, or of survival functions in the upper tail.
- **Quantiles:** a CDF that is linear between the predicted quantiles, with Gaussian tails matched to the outer two quantiles on each side.
- **Samples:** a Gaussian kernel on each sample, with Silverman's width but at least one bin. This is a kernel score, not the sampler's own distribution.
- **Bar distribution (TabPFN):** this uses the density TabPFN is trained with (FullSupportBarDistribution.forward), not its linear `cdf()`. Inner buckets are uniform. Each outer bucket's mass runs outward from its inner border as a half-normal with scale width / 0.6745, so half the mass lies within one bucket width of that border.

## Command line
```
python -m codelen preds.npz --verify [--out per_reading_bits.npz]
```
`preds.npz` has one entry per real reading:

| key | shape | content |
|---|---|---|
| `y` | (n,) float64 | the reading, °C |
| `mu`, `sd` | (n,) float64 | Gaussian forecast (sd in K); or one of the next two |
| `samples` | (n, S) float32 | generative forecast |
| `q_levels`, `q_values` | (k,) and (n, k) | quantile forecast |
| `bar_borders`, `bar_probs` | (B+1,) or (n, B+1), and (n, B) | bar distribution; bar_probs = softmax(logits) |
| `ref_bits` | (n,) float64, optional | reference bits (a fixture); the output adds `vs_ref_bits` |
| `block` | (n,) str | independent unit; bits are reported per block, with blocks won against each baseline |
| `hour` | (n,) int | with `block`, sets the coding order: block, then hour (stable) |
| `sensor` | (n,) int, optional | for the persistence baseline (the same sensor's last reading in the block) |
| `ctx` | (n, k) int, optional | row indices of context readings, −1 = none, for the context-mean baseline; they must come earlier in the coding order |
| `meta` | JSON string | model, set, units, delta, eps, lo, hi, post_hoc; delta, eps, lo and hi set the grid unless `--res`, `--eps`, `--lo` or `--hi` is given |

The output is JSON with four parts:
- the grid;
- for the model and each baseline: total bits, bits per reading and per-block bits, plus `roundtrip_ok` with `--verify`;
- the savings in bits per reading;
- the blocks won.

`--out` writes the per-reading bits, in the file's row order.

## Python interface
```python
import codelen
g = codelen.Grid(0.1, -60, 60, eps=1e-3)
codelen.score(y, codelen.Given(g, mu, sd), g, verify=True, per_reading=True)
codelen.score(y, codelen.GivenSamples(g, S), g)
codelen.score(y, codelen.GivenQuantiles(g, q_levels, q_values), g)
codelen.score(y, codelen.GivenBars(g, bar_borders, bar_probs), g)
codelen.score(y, codelen.ContextMean(g, ctx), g)       # baseline: mean of earlier context readings, sd fitted online
codelen.score(y, codelen.Persistence(g, group), g)     # baseline: the same group's last reading, sd fitted online
```
The array order is the coding order. A reading's forecast may use only earlier readings and side inputs that aren't scored. The scorer enforces this: a predictor's `y` argument holds NaN for the current reading and every later one. Any custom predictor needs only `reset()`, `table(i, y_seen)` and `update(i, k, y_seen)`.

## MDL for law-finding (`codelen.mdl`)
total = structure + parameters + residuals:
- **Structure:** the candidate's token list. An Elias-delta length, plus log2 |vocab| bits per token.
- **Parameters:** each one rounded to a step 2^-p. The code holds p, then round(value / step), both with signed Elias-gamma codes.
- **Residuals:** coded with codelen under N(f(x; rounded params), sd). The sd is sent as one more parameter.

```python
from codelen import mdl
r = mdl.fit_steps(tokens, vocab, f, params, x, y, res, lo, hi, verify=True)   # f(x, params) -> mean of y
```
`fit_steps` chooses each parameter's precision to minimise the total. The shortest total wins among candidates scored on the same data, grid and vocabulary.

## Tests
```
python -m codelen.test_codelen
```
The tests use synthetic data:
- The coder round trip.
- A 6-sensor field: the true model matches the entropy bound within 2%, and a biased model costs more.
- Both baselines, sample forecasts, and the no-peeking rule.
- The MDL helper picks the generating quadratic over a line and over a cubic.
- Per-reading agreement with `sparc_tf.bits`, for Gaussian and quantile forecasts.
- The pre-registered fixtures (`results/bits/fixture_{gaussian,quantile,bars}.npz`, made by `experiments/field/bits_fixture.py` and `bits_fixture_bars.py`; the bars reference is integrated from TabPFN 2.2.1's own `FullSupportBarDistribution.forward`) agree within 1e-4 bits per reading, with the coded length within overhead of the reference. These two parts are skipped if sparc-transformers isn't next to gridzip-lab.
- Bar distributions: bin masses sum to 1 and match the numerically integrated density.

Speed: about 2,000 readings per second with scipy, and about half that without.
