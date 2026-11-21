"""Self-tests: python -m codelen.test_codelen  (from the folder that holds codelen/). Synthetic data only.
Every scored set is also arithmetic-coded and decoded, and the round trip is asserted."""
import json, math, time
import numpy as np
from codelen import dists, scoring as sc, mdl, coder


def check(name, cond, info=""):
    print(f"{'ok  ' if cond else 'FAIL'} {name} {info}"); assert cond, name


def main():
    rng = np.random.default_rng(0); t0 = time.time()
    # 1. coder alone: random tables, random symbols
    enc = coder.Encoder(); syms = []
    for _ in range(3000):
        tot = int(rng.integers(2, 1 << 32)); lo = int(rng.integers(0, tot - 1)); hi = int(rng.integers(lo + 1, tot + 1))
        syms.append((lo, hi, tot)); enc.encode(lo, hi, tot)
    bits = enc.finish(); ideal = sum(math.log2(t / (h - l)) for l, h, t in syms); dec = coder.Decoder(bits); ok = True
    for l, h, t in syms:
        v = dec.target(t); ok &= l <= v < h; dec.consume(l, h, t)
    check("coder round trip", ok, f"{len(bits)} bits for {ideal:.1f} ideal")
    check("coder overhead < 4 bits", len(bits) - ideal < 4)

    # 2. a heat-field-like set: 6 sensors x 400 steps, a shared daily wave, per-sensor offsets, noise 0.15 C
    S, T, res, lo, hi = 6, 400, 0.01, -40.0, 60.0
    t = np.arange(T); wave = 20 + 5 * np.sin(2 * np.pi * t / 96); off = rng.normal(0, 1.0, S)
    truth = wave[:, None] + off[None, :]; Y = truth + rng.normal(0, 0.15, (T, S))
    y = Y.reshape(-1); group = np.tile(np.arange(S), T)            # coding order: time-major, sensors 0..5 per step
    idx = np.arange(len(y)).reshape(T, S)
    ctx = np.full((len(y), S - 1), -1)                              # context: the sensors already coded this step
    for s in range(1, S): ctx[idx[:, s], :s] = idx[:, :s]
    g = dists.Grid(res, lo, hi)
    good = sc.score(y, sc.Given(g, truth.reshape(-1), np.full(len(y), 0.15)), g, verify=True)
    ideal_ref = len(y) * (0.5 * math.log2(2 * math.pi * math.e * 0.15 ** 2) - math.log2(res))
    check("true model round trip", good["roundtrip_ok"], f"{good['bits_per_reading']:.3f} b/reading")
    check("true model ~ differential entropy", abs(good["total_bits"] - ideal_ref) / ideal_ref < 0.02,
          f"(theory {ideal_ref / len(y):.3f})")
    check("coded length ~ ideal", good["coded_bits"] - good["total_bits"] < 16, f"{good['coded_bits']} vs {good['total_bits']:.1f}")
    bad = sc.score(y, sc.Given(g, truth.reshape(-1) + 0.5, np.full(len(y), 0.15)), g, verify=True)
    check("biased model costs more", bad["total_bits"] > good["total_bits"], f"{bad['bits_per_reading']:.3f} b/reading")
    cm = sc.score(y, sc.ContextMean(g, ctx), g, verify=True)
    ps = sc.score(y, sc.Persistence(g, group), g, verify=True)
    check("baselines round trip", cm["roundtrip_ok"] and ps["roundtrip_ok"],
          f"context {cm['bits_per_reading']:.3f}, persistence {ps['bits_per_reading']:.3f} b/reading")
    check("true model beats baselines", good["total_bits"] < min(cm["total_bits"], ps["total_bits"]))
    # samples: 200 draws from the true forecast
    Sm = truth.reshape(-1)[:, None] + rng.normal(0, 0.15, (len(y), 200))
    sm = sc.score(y[:600], sc.GivenSamples(g, Sm[:600]), g, verify=True)
    gs = sc.score(y[:600], sc.Given(g, truth.reshape(-1)[:600], np.full(600, 0.15)), g)
    check("samples round trip", sm["roundtrip_ok"], f"{sm['bits_per_reading']:.3f} vs exact {gs['bits_per_reading']:.3f} b/reading")
    # a reading far outside the forecast still codes (uniform floor)
    far = sc.score(np.array([55.0]), sc.Given(g, [0.0], [0.01]), g, verify=True)
    check("ruled-out reading costs <= max bits", far["roundtrip_ok"] and far["total_bits"] <= g.max_bits() + 1e-6, f"{far['total_bits']:.1f} bits")
    # no peeking: a predictor that reads the current reading sees NaN
    class Peek(sc.Given):
        def table(self, i, y): assert np.isnan(y[i]); return super().table(i, y)
    sc.score(y[:50], Peek(g, truth.reshape(-1), np.full(len(y), 0.15)), g, verify=True); check("no peeking", True)

    # 3. MDL: data from y = 1.5 x^2 - 2 x + 3; the quadratic should beat a line and a cubic
    x = np.linspace(-3, 3, 300); yv = 1.5 * x ** 2 - 2 * x + 3 + rng.normal(0, 0.2, len(x))
    V = ["x", "c", "+", "*", "^2", "^3"]
    lin = mdl.fit_steps(["c", "x", "*", "c", "+"], V, lambda x, p: p[0] * x + p[1], [-2.0, 7.5], x, yv, 0.01, -50, 50)
    pq = np.polyfit(x, yv, 2); pc = np.polyfit(x, yv, 3)
    quad = mdl.fit_steps(["c", "x", "^2", "*", "c", "x", "*", "+", "c", "+"], V,
                         lambda x, p: p[0] * x ** 2 + p[1] * x + p[2], list(pq), x, yv, 0.01, -50, 50, verify=True)
    cub = mdl.fit_steps(["c", "x", "^3", "*", "c", "x", "^2", "*", "+", "c", "x", "*", "+", "c", "+"], V,
                        lambda x, p: p[0] * x ** 3 + p[1] * x ** 2 + p[2] * x + p[3], list(pc), x, yv, 0.01, -50, 50)
    print(f"     MDL bits: line {lin['total_bits']:.0f}, quadratic {quad['total_bits']:.0f}, cubic {cub['total_bits']:.0f}")
    check("MDL picks the quadratic", quad["total_bits"] < min(lin["total_bits"], cub["total_bits"]))
    check("MDL residuals round trip", quad["roundtrip_ok"])
    # 4. agreement with the analytic reference sparc_tf.bits (sparc-transformers), if it is on this machine
    import os, sys
    ref = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "sparc-transformers", "src")
    try:
        sys.path.insert(0, ref); from sparc_tf import bits as sb
    except Exception as e:
        print(f"skip sparc_tf.bits cross-check ({type(e).__name__})")
    else:
        g = dists.Grid(sb.DELTA, sb.LO, sb.HI, sb.EPS); n = 3000
        mu = rng.normal(15, 8, n); sd = np.exp(rng.normal(0, 1.2, n)); yv = mu + sd * rng.standard_t(4, n)
        yv[:20] = rng.uniform(-59, 59, 20); yv = np.clip(yv, -60, 60)          # some far misses
        mine = sc.score(yv, sc.Given(g, mu, sd), g, verify=True, per_reading=True); theirs = sb.gaussian_bits(yv, mu, sd)
        dmax = float(np.max(np.abs(mine["bits"] - theirs)))
        check("Gaussian bits = sparc_tf.bits per reading", dmax < 1e-4, f"max diff {dmax:.2e} bits")
        check("coded length = sparc_tf total + overhead", 0 <= mine["coded_bits"] - theirs.sum() < 16,
              f"{mine['coded_bits']} vs {theirs.sum():.1f}")
        qs = np.linspace(0.01, 0.99, 99); V = mu[:, None] + sd[:, None] * dists._ndtri(qs)[None] * 1.3
        mq = sc.score(yv, sc.GivenQuantiles(g, qs, V), g, verify=True, per_reading=True)
        tq = sb.cdf_bits(yv, sb.quantile_cdf(qs, V))
        dq = float(np.max(np.abs(mq["bits"] - tq)))
        check("quantile bits = sparc_tf.bits per reading", dq < 1e-4, f"max diff {dq:.2e} bits")
        # the pre-registered fixtures (sparc-transformers experiments/field/bits_fixture.py), if generated
        for name, mk in (("gaussian", lambda d, g: sc.Given(g, d["mu"], d["sd"])),
                         ("quantile", lambda d, g: sc.GivenQuantiles(g, d["q_levels"], d["q_values"])),
                         ("bars", lambda d, g: sc.GivenBars(g, d["bar_borders"], d["bar_probs"]))):   # ref from TabPFN's own forward()
            fp = os.path.join(os.path.dirname(ref), "results", "bits", f"fixture_{name}.npz")
            if not os.path.exists(fp): print(f"skip fixture_{name} (not generated)"); continue
            d = np.load(fp); m = json.loads(str(d["meta"])); g = dists.Grid(m["delta"], m["lo"], m["hi"], m["eps"])
            order = np.lexsort((d["hour"], d["block"].astype(str))); inv = np.argsort(order)
            r = sc.score(d["y"][order], mk({k: (d[k][order] if d[k].ndim and len(d[k]) == len(d["y"]) else d[k]) for k in d.files}, g),
                         g, verify=True, per_reading=True)
            dm = float(np.max(np.abs(r["bits"][inv] - d["ref_bits"])))
            check(f"fixture_{name}", r["roundtrip_ok"] and dm < 1e-4 and 0 <= r["coded_bits"] - d["ref_bits"].sum() < 16,
                  f"max diff {dm:.1e} bits, coded {r['coded_bits']} vs ref {d['ref_bits'].sum():.1f}")
    # 5. bar distribution: bin masses over a wide grid sum to 1, and match the density integrated numerically
    gw = dists.Grid(0.1, -200, 200, 1e-3)
    B = np.linspace(0, 30, 61); P = rng.dirichlet(np.ones(60)); t = dists.bars(gw, B, P)
    room = dists.SCALE - gw.K; cnt = np.diff(t.W)
    mass = ((cnt - 1) / room - 1e-3 * gw.u) / (1 - 1e-3); hn = dists.HN_MEDIAN
    def dens(v):
        w = 0.5; out = np.zeros_like(v); inner = (v >= B[1]) & (v <= B[-2])
        k = np.clip(np.searchsorted(B, v[inner], side="right") - 1, 1, 58); out[inner] = P[k] / w
        sl = w / hn; L = v < B[1]; out[L] = P[0] * 2 / (sl * math.sqrt(2 * math.pi)) * np.exp(-0.5 * ((B[1] - v[L]) / sl) ** 2)
        R = v > B[-2]; out[R] = P[-1] * 2 / (sl * math.sqrt(2 * math.pi)) * np.exp(-0.5 * ((v[R] - B[-2]) / sl) ** 2)
        return out
    E = gw.edges(t.a, t.b); fine = np.linspace(0, 1, 2001)
    num = np.array([np.trapezoid(dens(lo_ + fine * 0.1), dx=0.1 / 2000) for lo_ in E[:-1]])
    check("bar masses sum to 1", abs(mass.sum() - 1) < 1e-6, f"{mass.sum():.8f}")
    check("bar masses = integrated density", np.max(np.abs(mass - num)) < 1e-5, f"max diff {np.max(np.abs(mass - num)):.1e}")
    yb = rng.uniform(-3, 33, 300); rb = sc.score(yb, sc.GivenBars(gw, B, np.tile(P, (300, 1))), gw, verify=True)
    check("bars round trip", rb["roundtrip_ok"], f"{rb['bits_per_reading']:.3f} b/reading")
    print(f"all tests passed in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
