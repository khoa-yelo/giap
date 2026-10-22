"""Exact oracles for tilted (p, complement) pairs (M8-04a review B1; spec M8-04 lead decision 11).

- ``cpb_pairs_exact(p_in, c_in, p_out, c_out, k, N)``: the conditional law of K given S = N when
  every unit carries its own (p, comp) weights, by the rational DP. The conditional law depends on
  the units only through the odds p / comp, so a pair whose float sum is not exactly 1 (p == 1.0
  with comp == 1e-21) is still a well-defined unit;
- ``cpb_grouped_exact(p_in, out_groups, k, N)``: the law of the RAW q for a large complement made
  of a few groups of equal q (binomial terms evaluated only at the counts needed), so that a
  universe of 1000 units stays cheap;
- ``tilted_reference(q, p, comp)``: the exact tilted (p, comp) of the raw q in 60-digit decimal
  arithmetic, expit(x + theta) and expit(-(x + theta)) with x = logit q computed from the exact
  rational q, never as 1 - p. theta is the median of log(p / comp) - x over the open entries the
  implementation returned with both values normal (a single theta; float rounding of theta moves
  each entry by ~1e-14 relative).

``Fraction(float(x))`` is the exact value the implementation receives (as in ``cpb_oracle``).
"""

from __future__ import annotations

import statistics
from collections.abc import Sequence
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import cache
from math import comb

import numpy as np

from helpers.cpb_oracle import ExactCpb

PREC = 60


def _finish(raw: list[Fraction], k: int) -> ExactCpb | None:
    total = sum(raw, Fraction(0))
    if total == 0:
        return None
    w = tuple(x / total for x in raw)
    nz = [j for j, x in enumerate(w) if x != 0]
    return ExactCpb(w, k, nz[0], nz[-1])


def pb_pairs_exact(p: Sequence[float], c: Sequence[float]) -> list[Fraction]:
    pmf = [Fraction(1)]
    for x, y in zip(p, c, strict=True):
        fx, fy = Fraction(float(x)), Fraction(float(y))
        nxt = [Fraction(0)] * (len(pmf) + 1)
        for j, v in enumerate(pmf):
            if v:
                nxt[j] += v * fy
                nxt[j + 1] += v * fx
        pmf = nxt
    return pmf


def cpb_pairs_exact(
    p_in: Sequence[float],
    c_in: Sequence[float],
    p_out: Sequence[float],
    c_out: Sequence[float],
    k: int,
    N: int,
) -> ExactCpb | None:
    a = pb_pairs_exact(p_in, c_in)
    b = pb_pairs_exact(p_out, c_out)
    raw = [a[j] * (b[N - j] if 0 <= N - j < len(b) else Fraction(0)) for j in range(len(a))]
    return _finish(raw, k)


def cpb_grouped_exact(
    p_in: Sequence[float], out_groups: Sequence[tuple[float, int]], k: int, N: int
) -> ExactCpb | None:
    """The law of the raw q: C = ``p_in`` (small), the complement = ``count`` units of each q."""
    a = pb_pairs_exact(p_in, [1.0 - Fraction(float(x)) for x in p_in])  # exact 1 - q
    groups = tuple((Fraction(float(q)), n) for q, n in out_groups)

    @cache
    def pb(g: int, m: int) -> Fraction:  # PB of groups g.. at count m
        if g == len(groups):
            return Fraction(1) if m == 0 else Fraction(0)
        q, n = groups[g]
        s = Fraction(0)
        for i in range(max(0, m - sum(c for _, c in groups[g + 1 :])), min(n, m) + 1):
            s += _binom_term(q, n, i) * pb(g + 1, m - i)
        return s

    raw = [a[j] * pb(0, N - j) if N - j >= 0 else Fraction(0) for j in range(len(a))]
    return _finish(raw, k)


def _binom_term(q: Fraction, n: int, i: int) -> Fraction:
    return comb(n, i) * q**i * (1 - q) ** (n - i)


def _dec(x: Fraction) -> Decimal:
    return Decimal(x.numerator) / Decimal(x.denominator)


def tilted_reference(
    q: np.ndarray, p: np.ndarray, comp: np.ndarray
) -> tuple[list[Decimal | None], list[Decimal | None]]:
    """Exact tilted (p, comp) of the open entries of q (None for exact 0 and 1 of q)."""
    with localcontext() as ctx:
        ctx.prec = PREC
        xs: list[Decimal | None] = []
        thetas = []
        for qi, pi, ci in zip(q, p, comp, strict=True):
            if qi in (0.0, 1.0):
                xs.append(None)
                continue
            f = Fraction(float(qi))
            x = _dec(f).ln() - _dec(1 - f).ln()
            xs.append(x)
            if pi > 1e-300 and ci > 1e-300:
                thetas.append(_dec(Fraction(float(pi))).ln() - _dec(Fraction(float(ci))).ln() - x)
        theta = statistics.median(thetas)
        ref_p: list[Decimal | None] = []
        ref_c: list[Decimal | None] = []
        for x in xs:
            if x is None:
                ref_p.append(None)
                ref_c.append(None)
                continue
            s = x + theta
            ref_p.append(1 / (1 + (-s).exp()))
            ref_c.append(1 / (1 + s.exp()))
        return ref_p, ref_c


def dec_rel_err(got: float, exact: Decimal) -> float:
    with localcontext() as ctx:
        ctx.prec = PREC
        return float(abs((Decimal(got) - exact) / exact))
