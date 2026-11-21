"""MDL for law-finding: total bits = bits for the equation (structure + parameters) + bits for the residuals.

The decoder is assumed to know the vocabulary (the operators and variables the search may use) and the inputs x;
it must learn the equation and then the readings y. Every part is a prefix code, so the parts add:
  structure   Elias-delta length of the token sequence + log2(|vocab|) per token (a prefix / RPN token list)
  parameters  each parameter rounded to a step 2**-p: p with the zigzag Elias-gamma code, then the signed integer
              round(value / step) with a zigzag Elias-gamma code
  residuals   readings coded with codelen under N(f(x; rounded params), sd), sd itself sent as one more parameter
`mdl` scores one candidate at given steps; `fit_steps` tries steps 2**-p and keeps the shortest total (the MDL
trade-off between precise parameters and small residuals).
"""
import math
import numpy as np
from . import scoring as sc


def gamma_bits(n):
    """Elias gamma, n >= 1."""
    assert n >= 1; return 2 * (n.bit_length() - 1) + 1


def delta_bits(n):
    """Elias delta, n >= 1."""
    assert n >= 1; L = n.bit_length(); return (L - 1) + gamma_bits(L)


def signed_bits(q):
    """Zigzag map to >= 0, then gamma of (value + 1)."""
    z = 2 * q if q >= 0 else -2 * q - 1
    return gamma_bits(z + 1)


def structure_bits(tokens, vocab):
    vocab = list(vocab); bad = [t for t in tokens if t not in vocab]
    assert not bad, f"tokens not in vocab: {bad}"
    return delta_bits(len(tokens) + 1) + len(tokens) * math.log2(max(len(vocab), 1))


def param_bits(value, p):
    """Bits for one parameter at step 2**-p (p may be negative), and the rounded value the decoder gets."""
    step = 2.0 ** -p; q = int(round(value / step))
    return signed_bits(p) + signed_bits(q), q * step


def mdl(tokens, vocab, f, params, steps, x, y, res, lo, hi, sd_step=None, verify=False):
    """f(x, params) -> predicted mean of y. steps: one p per parameter (step 2**-p). Returns a dict of bits."""
    s_bits = structure_bits(tokens, vocab)
    pb, rounded = zip(*[param_bits(v, p) for v, p in zip(params, steps)]) if len(params) else ((), ())
    mu = np.asarray(f(x, list(rounded)), np.float64) * np.ones(len(y))
    sd0 = float(np.sqrt(np.mean((np.asarray(y) - mu) ** 2)))
    p_sd = sd_step if sd_step is not None else max(0, -math.floor(math.log2(max(sd0, res)))) + 3
    sd_b, sd_r = param_bits(sd0, p_sd); sd_r = max(sd_r, res / 2)
    r = sc.score_gaussian(y, mu, np.full(len(y), sd_r), res, lo, hi, verify=verify)
    total = s_bits + sum(pb) + sd_b + r["total_bits"]
    return {"structure_bits": s_bits, "param_bits": float(sum(pb)), "sd_bits": sd_b, "residual_bits": r["total_bits"],
            "total_bits": total, "bits_per_reading": total / len(y), "rounded_params": list(rounded), "sd": sd_r,
            "steps": list(steps), **({"roundtrip_ok": r["roundtrip_ok"]} if verify else {})}


def fit_steps(tokens, vocab, f, params, x, y, res, lo, hi, p_range=range(-4, 25), verify=False):
    """Greedy coordinate search over each parameter's step; returns the best mdl() result."""
    steps = [10] * len(params); best = mdl(tokens, vocab, f, params, steps, x, y, res, lo, hi)
    for _ in range(3):
        changed = False
        for j in range(len(params)):
            for p in p_range:
                trial = steps[:j] + [p] + steps[j + 1:]
                r = mdl(tokens, vocab, f, params, trial, x, y, res, lo, hi)
                if r["total_bits"] < best["total_bits"] - 1e-9: best, steps, changed = r, trial, True
        if not changed: break
    return mdl(tokens, vocab, f, params, steps, x, y, res, lo, hi, verify=verify) if verify else best
