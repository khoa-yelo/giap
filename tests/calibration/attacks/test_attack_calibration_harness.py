"""Attack on the calibration harness itself (GIAP_SPEC.md §8.4, §8.5; statistician-owned).

A gate is only useful if it can fail. These tests run deliberately anti-conservative procedures
through ``_criteria`` at the configured replicate counts and require the harness to flag them, while
the matching exactly-at-level procedure (the worst case for false alarms) passes.

- FDR, stream ``gene_fdr.attack_harness_*``: m = 2,000 p-values, 1,600 U(0, 1) nulls and 400
  Beta(0.05, 1) alternatives; oracle BH (step-up at level t / pi0 with the true pi0 = 0.8). With
  independent p, FDR = pi0 * (t / pi0) = t exactly (Benjamini & Hochberg 1995). ``valid`` uses
  t = q (FDR = q, on the boundary), ``inflated`` uses t = 1.5 q and 2 q, gated through
  ``calibrate`` and required to be a confirmed ``fail`` (PR #69 decision 3). R =
  ``replicates_fast.gene_fdr`` / ``replicates_full.gene_fdr``.
- Level, streams ``regions.attack_harness_*``: ``valid`` p ~ U(0, 1) (rate = alpha exactly) over
  200 independent streams, whose failure share must match the exact binomial false-alarm rate;
  ``inflated`` p = U / 1.5 (rate = 1.5 alpha). R = ``replicates_full.regions``. Power of the level
  check against 1.5 alpha at R = 1,000 is 0.34 / 0.92 / 0.998 for alpha = 0.01 / 0.05 / 0.10
  (exact binomial, docs/methods/calibration.md), so the attack requires a flag at alpha = 0.10 and
  only records alpha = 0.01 and 0.05.
- NaN guard: a vector of NaN p-values must not pass ``check_level`` (fail closed).

- ADR-0015 (two-block confirmation) through ``calibrate``: exact-level p give a confirmed ``fail``
  at the rate f^2 over 200 streams, and BH at FDR = 1.5 q is still a confirmed ``fail``.

BH is written out here (not imported from ``giap.stats``) so the attack tests the harness alone.
"""

from __future__ import annotations

import math
from pathlib import Path

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels
from scipy.stats import binom

M, M0 = 2000, 1600
PI0 = M0 / M
INFLATION = 1.5
S_STREAMS = 200

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def _bh_reject(p: np.ndarray, t: float) -> np.ndarray:
    """Benjamini-Hochberg step-up rejection set at level t."""
    m = p.size
    order = np.argsort(p, kind="stable")
    below = p[order] <= t * np.arange(1, m + 1) / m
    out = np.zeros(m, dtype=bool)
    if below.any():
        out[order[: int(np.flatnonzero(below)[-1]) + 1]] = True
    return out


def _fdp(levels: tuple[float, ...], factor: float):  # type: ignore[no-untyped-def]
    def rep(rng: np.random.Generator) -> np.ndarray:
        p = np.concatenate([rng.uniform(size=M0), rng.beta(0.05, 1.0, size=M - M0)])
        null = np.arange(M) < M0
        out = np.zeros(len(levels))
        for j, q in enumerate(levels):
            called = _bh_reject(p, min(1.0, factor * q / PI0))
            n = int(called.sum())
            out[j] = (called & null).sum() / n if n else 0.0
        return out

    return rep


def _fdr_results(full: bool, kind: str, factor: float) -> list[C.CheckResult]:
    levels = C.load_tolerances().levels
    r = C.replicates("gene_fdr", full=full)
    fdp = np.vstack(C.run_replicates(f"gene_fdr.attack_harness_{kind}", _fdp(levels, factor), r))
    return [C.check_fdr(fdp[:, j], q) for j, q in enumerate(levels)]


@pytest.mark.parametrize("full", MODES)
def test_fdr_check_passes_exact_level_bh(full: bool) -> None:
    res = _fdr_results(full, "valid", 1.0)
    assert all(x.ok for x in res), [x.detail for x in res]


@pytest.mark.parametrize("factor", [INFLATION, 2.0])
@pytest.mark.parametrize("full", MODES)
def test_fdr_check_flags_inflated_bh(full: bool, factor: float, teeth_log: Path) -> None:
    """Teeth (PR #69 decision 3): oracle BH at FDR = 1.5 q and 2 q is a confirmed ``fail``
    (both blocks, ADR-0015) through ``calibrate`` at every q."""
    levels = C.load_tolerances().levels
    r = C.replicates("gene_fdr", full=full)
    gate = calibrate_levels(
        "gene_fdr", "fdr", levels, _fdp(levels, factor), r,
        stream=f"gene_fdr.attack_harness_inflated_{factor:g}", log_dir=teeth_log,
    )  # fmt: skip
    for rec, q in zip(gate.records, levels, strict=True):
        assert rec["observed"] > q  # the procedure really is anti-conservative
    assert [rec["status"] for rec in gate.records] == ["fail"] * len(levels), gate.detail()


def _level_results(stream: str, divisor: float) -> list[C.CheckResult]:
    r = C.replicates("regions", full=True)
    p = np.asarray(C.run_replicates(stream, lambda g: g.uniform() / divisor, r))
    return [C.check_level(p, a) for a in C.load_tolerances().levels]


def _false_alarm(alpha: float, r: int) -> float:
    """Exact P(level check fails) when the true rejection rate is exactly alpha."""
    c = math.floor(C.level_bound(alpha, r) * r + 1e-9)  # fails iff count > c
    return float(binom.sf(c, r, alpha))


@pytest.mark.calibration_full
def test_level_check_false_alarm_matches_binomial_theory() -> None:
    """Exactly-valid p (U(0, 1)) over S independent streams: the share of failing streams per alpha
    stays within the 99.9 % binomial band of the exact false-alarm rate (about 0.026-0.028 at
    R = 1,000). A single stream can and does fail: ``regions.attack_harness_valid`` has 122 of
    1,000 p <= 0.10 (P = 0.013), which is why this check is made over streams."""
    r = C.replicates("regions", full=True)
    levels = C.load_tolerances().levels
    fails = np.zeros(len(levels), dtype=int)
    for s in range(S_STREAMS):
        res = _level_results(f"regions.attack_harness_valid_{s}", 1.0)
        fails += [not x.ok for x in res]
    for a, f in zip(levels, fails, strict=True):
        assert f <= binom.ppf(0.999, S_STREAMS, _false_alarm(a, r)), (a, f)


@pytest.mark.calibration_full
def test_level_check_flags_one_and_a_half_alpha() -> None:
    levels = C.load_tolerances().levels
    res = dict(
        zip(levels, _level_results("regions.attack_harness_inflated", INFLATION), strict=True)
    )
    assert not res[0.10].ok, res[0.10].detail


def test_level_check_fails_closed_on_nan() -> None:
    """All-NaN p-values (a method that returned nothing) must not pass the level gate."""
    try:
        res = C.check_level([float("nan")] * 100, 0.05)
    except ValueError:
        return  # refusing the input is failing closed
    assert not res.ok, res.detail  # passing is the bug


def test_level_check_fails_closed_on_out_of_range_p() -> None:
    """p-values above 1 are a method bug; they must not count as 'not rejected' and pass."""
    try:
        res = C.check_level([2.0] * 100, 0.05)
    except ValueError:
        return
    assert not res.ok, res.detail


# --------------------------------------------------------------------------- ADR-0015 (calibrate)


@pytest.mark.calibration_full
def test_two_block_confirmed_false_alarm_matches_theory(tmp_path: Path) -> None:
    """Exact-level U(0, 1) p through ``calibrate`` over S streams: the share of confirmed ``fail``
    per alpha stays in the 99.9 % binomial band of f^2 (f = single-block false alarm, ≈ 0.027 at
    R = 1,000, so f^2 ≈ 0.0008), and ``unconfirmed`` occurs at about f."""
    r = C.replicates("regions", full=True)
    for a in C.load_tolerances().levels:
        status = [
            C.calibrate(
                "regions",
                "level",
                a,
                lambda g: g.uniform(),
                r,
                stream=f"regions.attack_harness_adr15_{s}",
                log_dir=tmp_path,
            )["status"]
            for s in range(S_STREAMS)
        ]
        f = _false_alarm(a, r)
        assert status.count("fail") <= binom.ppf(0.999, S_STREAMS, f * f), (a, status)
        n_unconfirmed = status.count("unconfirmed")
        assert n_unconfirmed <= binom.ppf(0.999, S_STREAMS, f), (a, n_unconfirmed)


@pytest.mark.parametrize("full", MODES)
def test_two_block_confirms_bh_at_one_and_a_half_q(full: bool, tmp_path: Path) -> None:
    """Confirmation must not rescue an anti-conservative procedure: oracle BH at FDR = 1.5 q is a
    confirmed ``fail`` at every q, and BH at FDR = q is never a confirmed ``fail``."""
    levels = C.load_tolerances().levels
    r = C.replicates("gene_fdr", full=full)
    for j, q in enumerate(levels):
        for kind, factor in (("valid", 1.0), ("inflated", INFLATION)):
            rep = _fdp(levels, factor)
            rec = C.calibrate(
                "gene_fdr",
                "fdr",
                q,
                lambda g, rep=rep, j=j: rep(g)[j],
                r,
                stream=f"gene_fdr.attack_harness_adr15_{kind}",
                log_dir=tmp_path,
            )
            if kind == "inflated":
                assert rec["status"] == "fail", rec
            else:
                assert rec["status"] != "fail", rec
