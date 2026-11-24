"""PROTOCOL-7: candidate terms the library does not have, built by a small grammar from a case's raw inputs.

Written fresh for SPARC / gridzip-lab. A candidate is either one input under one operation, or the product of two:

  single   op(x)                    op in OPS (14 choices, three of them with a constant at the 25/50/75 % quantile)
  product  a(x1) · b(x2), x1 != x2  a, b in BASIC (id, log, sqrt, sq, inv)         'basic product'
           gate(x1) · x2            gate in decay/hinge/step at a quantile         'gated product'

Its description length (bits) is what a decoder needs to rebuild it: 1 bit single or product, 1 bit which kind of
product, log2(ops) per operation, log2(inputs) per input, log2(3) per quantile constant. The engine adds one flag bit
(library or grammar) and the '+' token, so a grammar term always pays more to be named than a library term: a tie
goes to the named physics. The grammar is a search vocabulary, not physics: a grammar term that wins says "the
errors follow this shape", and naming the mechanism stays a human step.
"""
import math

import numpy as np

from .engine import Term

QS = (0.25, 0.5, 0.75)


def _pos(x): return np.all(x > 0)


def _nonneg(x): return np.all(x >= 0)


def _inv_ok(x):
    ax = np.abs(x); return (np.all(x > 0) or np.all(x < 0)) and ax.min() > 1e-3 * ax.max()


# name, formula pattern, applies(x), f(x, const), takes a quantile constant
BASIC = [
    ("id", "{x}", lambda x: True, lambda x, q: x, False),
    ("log", "ln {x}", _pos, lambda x, q: np.log(x), False),
    ("sqrt", "√{x}", _nonneg, lambda x, q: np.sqrt(x), False),
    ("sq", "{x}²", lambda x: True, lambda x, q: x ** 2, False),
    ("inv", "1/{x}", _inv_ok, lambda x, q: 1.0 / x, False),
]
GATES = [
    ("decay", "e^(−{x}/{q})", lambda x: _nonneg(x) and len(np.unique(x)) > 3, lambda x, q: np.exp(-x / q), True),
    ("hinge", "max(0, {x}−{q})", lambda x: len(np.unique(x)) > 3, lambda x, q: np.maximum(0.0, x - q), True),
    ("step", "[{x}>{q}]", lambda x: len(np.unique(x)) > 2, lambda x, q: (x > q).astype(float), True),
]
OPS = BASIC + GATES


def primitives(ctx, n, exclude=()):
    """Every numeric (n,) input in ctx that varies: the grammar's alphabet."""
    out = {}
    for k, v in ctx.items():
        if k in exclude or callable(v): continue
        a = np.asarray(v) if not np.isscalar(v) else None
        if a is None or a.ndim != 1 or a.shape[0] != n or not np.issubdtype(a.dtype, np.number): continue
        a = a.astype(np.float64)
        if np.all(np.isfinite(a)) and a.std() > 1e-12 * (1 + np.abs(a).max()): out[k] = a
    return out


def _consts(x, takes):
    if not takes: return [(None, "")]
    xs = x[x > 0] if np.all(x >= 0) and np.any(x > 0) else x
    out = []
    for q in QS:
        v = float(np.quantile(xs, q))
        if v == 0 and np.all(x >= 0): continue
        out.append((v, f"{v:.3g}"))
    return out


def _forms(P, ops):
    """(label, column, bits for op + input + constant) for each op on each input."""
    lp = math.log2(len(P)); lo = math.log2(len(ops)); out = []
    for k, x in P.items():
        for name, pat, ok, f, takes in ops:
            if not ok(x): continue
            for qv, qs in _consts(x, takes):
                v = f(x, qv)
                if not np.all(np.isfinite(v)) or v.std() < 1e-12: continue
                out.append((k, pat.format(x=k, q=qs), v, lo + lp + (math.log2(len(QS)) if takes else 0.0), name))
    return out


def candidates(ctx, n, exclude=(), products=True, max_products=None):
    """Grammar terms (engine.Term with cost_bits set) for a case. Duplicates of each other or of library terms are
    dropped later by engine.assemble (library terms come first in the pool)."""
    P = primitives(ctx, n, exclude)
    if not P: return []
    lp = math.log2(len(P))
    singles = _forms(P, OPS)
    terms = [Term(f"g:{lab}", "grammar", lab, "a shape the errors follow (grammar; name the mechanism by hand)", (k,), (lambda v: lambda c: v)(v),
                  physics=False, cost_bits=1.0 + b) for k, lab, v, b, _ in singles]
    if products:
        basic = [f for f in _forms(P, BASIC)]
        gates = [f for f in singles if f[4] in ("decay", "hinge", "step")]
        lb = math.log2(len(BASIC)); lg = math.log2(len(GATES))
        prods = []
        for i in range(len(basic)):
            for j in range(i + 1, len(basic)):
                a, b = basic[i], basic[j]
                if a[0] == b[0]: continue
                prods.append((f"{a[1]} · {b[1]}", a[2] * b[2], 2.0 + a[3] + b[3], (a[0], b[0])))
        for gk, glab, gv, gb, _ in gates:
            for k, x in P.items():
                if k == gk: continue
                prods.append((f"{glab} · {k}", gv * x, 2.0 + lg + lp + math.log2(len(QS)) + lp, (gk, k)))
        if max_products and len(prods) > max_products:
            prods = prods[:max_products]
        terms += [Term(f"g:{lab}", "grammar", lab, "a shape the errors follow (grammar; name the mechanism by hand)", ks,
                       (lambda v: lambda c: v)(v), physics=False, cost_bits=c) for lab, v, c, ks in prods if np.all(np.isfinite(v))]
    return terms
