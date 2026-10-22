"""Exact-rational oracles for the conditional Poisson-binomial law (issue M8-04a; spec M8-04 §1).

Every probability is converted with ``Fraction(float(p))``, i.e. the exact value of the float the
implementation receives, so the oracle and the implementation see the same input and any
difference is the implementation's rounding.

- ``pb_pmf_exact(p)``: the Poisson-binomial pmf by the rational DP (j = 0..len(p));
- ``cpb_exact(p_in, p_out, k, N)``: the law of K = sum over C given S = N from the two rational
  pmfs, normalised by the support sum;
- ``cpb_enumerate(p_in, p_out, k, N)``: the same law by exhaustive enumeration of the 2^T
  configurations (T <= 14), the spec's test oracle (§1 "Not used, and why").

Both return ``ExactCpb`` with the rational ``p_over`` (P(K >= k | N)), ``p_under``
(P(K <= k | N)), ``p_eq`` (P(K = k | N)), ``mean`` (E[K | N]) and the support ``lo``/``hi`` of the
j with non-zero weight; ``None`` when the observation (S = N) has probability 0.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

TINY = Fraction(1, 10**300)  # the §1 precision contract applies to tails >= 1e-300


@dataclass(frozen=True)
class ExactCpb:
    weights: tuple[Fraction, ...]  # P(K = j | S = N), j = 0..len(p_in)
    k: int
    lo: int
    hi: int

    @property
    def p_over(self) -> Fraction:
        return sum(self.weights[max(self.k, 0) :], Fraction(0))

    @property
    def p_under(self) -> Fraction:
        return sum(self.weights[: self.k + 1], Fraction(0)) if self.k >= 0 else Fraction(0)

    @property
    def p_eq(self) -> Fraction:
        return self.weights[self.k] if 0 <= self.k < len(self.weights) else Fraction(0)

    @property
    def mean(self) -> Fraction:
        return sum((j * w for j, w in enumerate(self.weights)), Fraction(0))


def fr(p: Sequence[float]) -> list[Fraction]:
    return [Fraction(float(x)) for x in p]


def pb_pmf_exact(p: Sequence[float]) -> list[Fraction]:
    pmf = [Fraction(1)]
    for x in fr(p):
        nxt = [Fraction(0)] * (len(pmf) + 1)
        for j, v in enumerate(pmf):
            if v:
                nxt[j] += v * (1 - x)
                nxt[j + 1] += v * x
        pmf = nxt
    return pmf


def _finish(raw: list[Fraction], k: int) -> ExactCpb | None:
    total = sum(raw, Fraction(0))
    if total == 0:
        return None
    w = tuple(x / total for x in raw)
    nz = [j for j, x in enumerate(w) if x != 0]
    return ExactCpb(w, k, nz[0], nz[-1])


def cpb_exact(p_in: Sequence[float], p_out: Sequence[float], k: int, N: int) -> ExactCpb | None:
    a = pb_pmf_exact(p_in)
    b = pb_pmf_exact(p_out)
    raw = [a[j] * (b[N - j] if 0 <= N - j < len(b) else Fraction(0)) for j in range(len(a))]
    return _finish(raw, k)


def cpb_enumerate(p_in: Sequence[float], p_out: Sequence[float], k: int, N: int) -> ExactCpb | None:
    q = fr(p_in) + fr(p_out)
    n = len(p_in)
    if len(q) > 14:
        raise ValueError("enumeration oracle is for T <= 14")
    raw = [Fraction(0)] * (n + 1)
    for d in itertools.product((0, 1), repeat=len(q)):
        if sum(d) != N:
            continue
        pr = Fraction(1)
        for x, di in zip(q, d, strict=True):
            pr *= x if di else 1 - x
        raw[sum(d[:n])] += pr
    return _finish(raw, k)


def hypergeom_exact(M: int, n: int, N: int, j: int) -> Fraction:
    """P(K = j) for K ~ hypergeometric(M items, n in C, N drawn), exactly."""
    if j < 0 or j > n or N - j < 0 or N - j > M - n:
        return Fraction(0)
    return Fraction(math.comb(n, j) * math.comb(M - n, N - j), math.comb(M, N))


def rel_err(got: float, exact: Fraction) -> float:
    if exact == 0:
        return abs(got)
    return abs(float((Fraction(got) - exact) / exact))


def assert_tail_close(got: float, exact: Fraction, rtol: float = 1e-12) -> None:
    """The §1 precision contract: relative ``rtol`` when exact >= 1e-300, else a float >= 0."""
    assert math.isfinite(got), got
    assert got >= 0.0, got
    if exact >= TINY:
        assert rel_err(got, exact) <= rtol, (got, float(exact), rel_err(got, exact))
    else:
        assert got <= 1e-290, (got, float(exact))
