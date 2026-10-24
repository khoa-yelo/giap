#!/usr/bin/env python3
"""Write the exact Wallenius tails of the golden ORA inputs (M8-12b §3; ADR-0044 decision 1).

    python scripts/golden/make_wallenius_exact.py [--inputs DIR] [--out FILE] [--workers INT]
                                                  [--categories CAT_a,CAT_b,...] [--log FILE]

R1 is the urn recursion in 80-bit long double (tails >= 1e-3). R2 is the exact rational integral
in fixed-point integers with >= 200 surviving bits (the smaller tail when it is < 1e-3); the
larger tail is then 1 - small + pmf(k). Standard library and numpy (long double) only: no scipy,
no R, no goseq or BiasedUrn code, no giap.
"""

from __future__ import annotations

import argparse
import csv
import functools
import gzip
import io
import json
import math
import multiprocessing
import os
import platform
import sys
from collections.abc import Callable, Iterable, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DEFAULT_INPUTS = REPO / "tests" / "data" / "golden" / "r" / "inputs"
DEFAULT_OUT = REPO / "tests" / "data" / "golden" / "wallenius" / "exact_tails.tsv.gz"
SMALL = 1e-3  # tails below this are R2 (strict <)
SURVIVING_BITS = 200
TAILS = ("over", "under")
Cat = tuple[int, int, int, int, float]  # (M, n, N, k, omega)


def support(M: int, n: int, N: int) -> tuple[int, int]:
    return max(0, N - (M - n)), min(n, N)


def odds(p: Sequence[Fraction], members: Sequence[int]) -> float:
    """ω = mean_in / mean_out, computed exactly and correctly rounded once to float64."""
    M, n = len(p), len(members)
    if n == 0 or n >= M:
        raise ValueError(f"category size {n} leaves no genes on one side")
    s_in, s_all = sum((p[i] for i in members), Fraction(0)), sum(p, Fraction(0))
    if s_in == 0 or s_all == s_in:
        raise ValueError("ω is undefined (a zero mean p_detect)")
    return float((s_in / n) / ((s_all - s_in) / (M - n)))


def r1_pmf(M: int, n: int, N: int, omega: float, dtype: Any = np.longdouble) -> np.ndarray:
    """The urn recursion over N draws; ``pmf[x - lo]`` for x = lo..hi."""
    lo, hi = support(M, n, N)
    m2, w = M - n, dtype(omega)
    probs = np.ones(1, dtype=dtype)
    for t in range(N):
        r = np.arange(probs.size, dtype=dtype)
        red = w * (n - r)
        white = np.maximum(dtype(m2) - (t - r), dtype(0))
        den = red + white
        zero = np.zeros_like(den)
        p_red = np.divide(red, den, out=zero.copy(), where=den > 0)
        p_white = np.divide(white, den, out=zero, where=den > 0)
        new = np.zeros(min(probs.size + 1, n + 1), dtype=dtype)
        new[: probs.size] += probs * p_white
        new[1 : probs.size + 1] += (probs * p_red)[: new.size - 1]
        probs = new
    return probs[lo : hi + 1]


@functools.lru_cache(maxsize=8)
def _r2_products(x: int, M: int, n: int, N: int, omega: float) -> tuple[int, int, list[int]]:
    """(m!·2^{e(m+1)}, e, [∏ᵢ A_{j,i} for j = 0..x]) with ω = W/2^e (spec §3.3 step 3)."""
    W, den = omega.as_integer_ratio()
    e = den.bit_length() - 1
    m, m1, m2 = N - x, n, M - n
    base = math.factorial(m) << (e * (m + 1))
    prods = []
    for j in range(x + 1):
        c = W * j + W * (m1 - x) + den * (m2 - N + x)
        prods.append(math.prod(c + den * i for i in range(m + 1)))
    return base, e, prods


def r2_sum(x: int, M: int, n: int, N: int, omega: float, P: int) -> int:
    """S = Σⱼ (-1)ʲ ⌊C(x,j)·m!·2^{e(m+1)}·2ᴾ / ∏ᵢ A_{j,i}⌋ (absolute error ≤ x + 1)."""
    base, _, prods = _r2_products(x, M, n, N, omega)
    num = base << P
    return sum((-1) ** j * (math.comb(x, j) * num // q) for j, q in enumerate(prods))


def r2_pmf(x: int, M: int, n: int, N: int, omega: float, *, p_start: int | None = None) -> Fraction:
    """f(x) by the exact integral, with relative error < 2^-199 (self-verifying precision)."""
    m1, m2 = n, M - n
    W, den = omega.as_integer_ratio()
    d_num = W * (m1 - x) + den * (m2 - N + x)  # D = d_num / den
    if d_num == 0:  # N = M: the point mass at x = n
        return Fraction(1)
    P = 264 + 2 * x if p_start is None else p_start
    target = SURVIVING_BITS + (x + 1).bit_length() + 1
    while True:
        S = r2_sum(x, M, n, N, omega, P)
        if S > 0 and S.bit_length() >= target:
            break
        P = P * 2 if S <= 0 else P + target - S.bit_length() + 32
    e = den.bit_length() - 1
    return Fraction(math.comb(m1, x) * math.comb(m2, N - x) * d_num * S, 1 << (P + e))


def _ld_fraction(values: Iterable[Any]) -> Fraction:
    return sum((Fraction(*v.as_integer_ratio()) for v in values), Fraction(0))


def _r2_tail(points: Sequence[int], pmf: Callable[[int], Fraction]) -> Fraction:
    return sum((pmf(x) for x in points), Fraction(0))


def _regime(
    k: int, M: int, n: int, N: int, pmf1: np.ndarray
) -> tuple[Fraction, Fraction, str | None]:
    """R1 tails (exact sums of the long-double terms, k = lo/hi exact 1) and the R2 tail, if any."""
    lo, hi = support(M, n, N)
    over = Fraction(1) if k == lo else _ld_fraction(pmf1[k - lo :])
    under = Fraction(1) if k == hi else _ld_fraction(pmf1[: k - lo + 1])
    s_tail = "over" if over <= under else "under"
    small = float(min(over, under)) < SMALL
    return over, under, s_tail if small else None


def _tail_points(tail: str, k: int, lo: int, hi: int) -> range:
    return range(k, hi + 1) if tail == "over" else range(lo, k + 1)


def exact_tails(
    k: int, M: int, n: int, N: int, omega: float, *, pmf: Callable[[int], Fraction] | None = None
) -> tuple[Fraction, Fraction, dict[str, str]]:
    """(P(X ≥ k), P(X ≤ k), method): R1 when both tails are ≥ 1e-3, else R2 and its complement."""
    lo, hi = support(M, n, N)
    if not lo <= k <= hi:
        raise ValueError(f"k = {k} outside the support [{lo}, {hi}]")
    over, under, s_tail = _regime(k, M, n, N, r1_pmf(M, n, N, omega))
    if s_tail is None:
        return over, under, {"over": "R1", "under": "R1"}
    f = pmf or (lambda x: r2_pmf(x, M, n, N, omega))
    small = _r2_tail(_tail_points(s_tail, k, lo, hi), f)
    large = 1 - small + f(k)
    big = "under" if s_tail == "over" else "over"
    tails = {s_tail: small, big: large}
    return tails["over"], tails["under"], {s_tail: "R2", big: "R2 complement"}


# --- the CLI -----------------------------------------------------------------------------------


def _read_tsv(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))[1:]


def _load(inputs: Path, only: list[str] | None) -> dict[str, Cat]:
    genes = _read_tsv(inputs / "ora_genes.tsv.gz")
    index = {g: i for i, (g, _, _) in enumerate(genes)}
    p = [Fraction(q) for _, q, _ in genes]
    det = [int(d) for _, _, d in genes]
    members: dict[str, list[int]] = {}
    for cat, gene in _read_tsv(inputs / "ora_categories.tsv.gz"):
        if gene not in index:
            raise SystemExit(f"{cat}: gene {gene} is not in ora_genes")
        members.setdefault(cat, []).append(index[gene])
    if only is not None:
        missing = [c for c in only if c not in members]
        if missing:
            raise SystemExit(f"unknown categories: {','.join(missing)}")
        members = {c: v for c, v in members.items() if c in only}
    if not members:
        raise SystemExit("no categories")
    M, N = len(genes), sum(det)
    cats: dict[str, Cat] = {}
    for cat, idx in members.items():
        try:
            omega = odds(p, idx)
        except ValueError as exc:
            raise SystemExit(f"{cat}: {exc}") from None
        k = sum(det[i] for i in idx)
        lo, hi = support(M, len(idx), N)
        if not lo <= k <= hi:
            raise SystemExit(f"{cat}: k = {k} outside the support [{lo}, {hi}]")
        cats[cat] = (M, len(idx), N, k, omega)
    return cats


def _numpy_odds(inputs: Path, cats: dict[str, Cat]) -> dict[str, float]:
    genes = _read_tsv(inputs / "ora_genes.tsv.gz")
    index = {g: i for i, (g, _, _) in enumerate(genes)}
    p = np.array([float(q) for _, q, _ in genes])
    mask: dict[str, np.ndarray] = {c: np.zeros(p.size, dtype=bool) for c in cats}
    for cat, gene in _read_tsv(inputs / "ora_categories.tsv.gz"):
        if cat in mask:
            mask[cat][index[gene]] = True
    return {c: float(np.mean(p[m]) / np.mean(p[~m])) for c, m in mask.items()}


def _r1_task(c: Cat) -> tuple[np.ndarray, np.ndarray]:
    M, n, N, _, omega = c
    return r1_pmf(M, n, N, omega), r1_pmf(M, n, N, omega, dtype=np.float64)


def _r2_task(job: tuple[int, int, int, int, float]) -> Fraction:
    x, M, n, N, omega = job
    return r2_pmf(x, M, n, N, omega)


def _map(fn: Callable[[Any], Any], jobs: list[Any], workers: int) -> list[Any]:
    if workers <= 1 or len(jobs) <= 1:
        return [fn(j) for j in jobs]
    with multiprocessing.Pool(min(workers, len(jobs))) as pool:
        return list(pool.imap(fn, jobs, chunksize=1))


def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b)


def run(inputs: Path, cats: dict[str, Cat], workers: int, full: bool) -> tuple[list[Any], Any]:
    """The tails of every category and the §3.4 check summary (SystemExit on a failed check)."""
    r1 = dict(zip(cats, _map(_r1_task, list(cats.values()), workers), strict=True))
    regimes = {c: _regime(c_[3], *c_[:3], r1[c][0]) for c, c_ in cats.items()}
    need: dict[tuple[str, int], tuple[int, int, int, int, float]] = {}

    def want(cat: str, points: Iterable[int]) -> None:
        M, n, N, _, omega = cats[cat]
        for x in points:
            need[(cat, x)] = (x, M, n, N, omega)

    for cat, (_, _, s_tail) in regimes.items():
        if s_tail is not None:
            M, n, N, k, _ = cats[cat]
            want(cat, _tail_points(s_tail, k, *support(M, n, N)))
    spot: list[tuple[str, str]] = []
    sumf: list[str] = []
    if full:
        costs = []
        for cat, (M, n, N, k, _) in cats.items():
            for tail in TAILS if regimes[cat][2] is None else ():
                pts = _tail_points(tail, k, *support(M, n, N))
                costs.append((sum((x + 1) * (N - x + 1) for x in pts), cat, tail))
        spot = [(c, t) for _, c, t in sorted(costs)[:10]]
        sumf = [c for _, c in sorted((len(range(*_span(cats[c]))), c) for c in cats)[:3]]
        for cat, tail in spot:
            M, n, N, k, _ = cats[cat]
            want(cat, _tail_points(tail, k, *support(M, n, N)))
        for cat in sumf:
            want(cat, range(*_span(cats[cat])))
    keys = sorted(need)
    r2 = dict(zip(keys, _map(_r2_task, [need[k] for k in keys], workers), strict=True))

    rows: list[tuple[str, str, float]] = []
    log: dict[str, Any] = {"categories": len(cats), "r2_tails": [], "checks": {}}
    worst64 = worst12 = 0.0
    for cat, (M, n, N, k, omega) in cats.items():
        lo = support(M, n, N)[0]
        over1, under1, s_tail = regimes[cat]
        p64 = r1[cat][1]
        f64 = {"over": math.fsum(p64[k - lo :].tolist()), "under": math.fsum(p64[: k - lo + 1])}
        ld = {"over": float(over1), "under": float(under1)}
        for tail in TAILS:
            worst64 = max(worst64, _rel(f64[tail], ld[tail]))
        if s_tail is None:
            vals = {"over": over1, "under": under1}
        else:
            vals = dict(
                zip(
                    ("over", "under"),
                    exact_tails(k, M, n, N, omega, pmf=lambda x, c=cat: r2[(c, x)])[:2],
                    strict=True,
                )
            )
            worst12 = max(worst12, _rel(ld[s_tail], float(vals[s_tail])))
            log["r2_tails"].append(f"{cat} {s_tail}")
        for tail in TAILS:
            p = float(vals[tail])
            if not 1e-300 <= p <= 1.0:
                raise SystemExit(f"{cat} {tail}: p = {p!r} outside [1e-300, 1]")
            rows.append((cat, tail, p))
        if float(vals["over"]) + float(vals["under"]) < 1.0:
            raise SystemExit(f"{cat}: over + under < 1")
    checks = log["checks"]
    checks["r1_longdouble_vs_float64_max_rel"] = _check(worst64, 1e-13, "R1 80-bit vs float64")
    checks["r1_vs_r2_max_rel"] = _check(worst12, 1e-13, "R1 vs R2 on the R2 tails")
    if full:
        checks["spot_checks"] = {}
        for cat, tail in spot:
            M, n, N, k, _ = cats[cat]
            got = _r2_tail(_tail_points(tail, k, *support(M, n, N)), lambda x, c=cat: r2[(c, x)])
            ref = regimes[cat][0] if tail == "over" else regimes[cat][1]
            rel = _rel(float(got), float(ref))
            checks["spot_checks"][f"{cat} {tail}"] = _check(rel, 1e-13, f"spot {cat} {tail}")
        checks["sum_f"] = {}
        for cat in sumf:
            err = abs(sum((r2[(cat, x)] for x in range(*_span(cats[cat]))), Fraction(0)) - 1)
            if err > Fraction(1, 10**50):
                raise SystemExit(f"{cat}: |Σf - 1| = {float(err):.3e} > 1e-50")
            checks["sum_f"][cat] = f"{float(err):.3e}"
        values = [p for _, _, p in rows]
        deep, mid = sum(0 < p < 1e-20 for p in values), sum(1e-20 <= p < 1e-8 for p in values)
        if deep < 2 or mid < 1:
            raise SystemExit(
                f"vacuity floor: {deep} tails < 1e-20 (need 2), {mid} in [1e-20, 1e-8)"
            )
        checks["vacuity"] = {"below_1e-20": deep, "1e-20_to_1e-8": mid}
    np_odds = _numpy_odds(inputs, cats)
    ulps = [abs(np_odds[c] - cats[c][4]) / float(np.spacing(cats[c][4])) for c in cats]
    log["omega_max_ulp_vs_numpy"] = max(ulps)
    return rows, log


def _span(c: Cat) -> tuple[int, int]:
    lo, hi = support(c[0], c[1], c[2])
    return lo, hi + 1


def _check(value: float, limit: float, what: str) -> str:
    if value > limit:
        raise SystemExit(f"{what}: max rel {value:.3e} > {limit:g}")
    return f"{value:.3e}"


def _write(path: Path, rows: list[tuple[str, str, float]]) -> None:
    text = "category\ttail\tp_exact\n" + "".join(
        f"{c}\t{t}\t{format(p, '.17g')}\n" for c, t, p in rows
    )
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, compresslevel=9, mtime=0) as gz:
        gz.write(text.encode("utf-8"))
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        tmp.write_bytes(buf.getvalue())
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs", type=Path, default=DEFAULT_INPUTS)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--categories", default=None)
    ap.add_argument("--log", type=Path, default=None)
    args = ap.parse_args(argv)
    info = np.finfo(np.longdouble)
    if info.nmant < 63:
        raise SystemExit(f"np.longdouble has nmant {info.nmant} < 63: 80-bit x87 is required")
    only = None if args.categories is None else [c for c in args.categories.split(",") if c]
    cats = _load(args.inputs, only)
    rows, log = run(args.inputs, cats, args.workers, full=only is None)
    log["environment"] = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "longdouble": f"nmant {info.nmant}, bits {info.bits}",
        "machine": platform.machine(),
    }
    _write(args.out, rows)
    if args.log is not None:
        args.log.write_text(json.dumps(log, indent=1, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
