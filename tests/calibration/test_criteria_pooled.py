"""Pooled level criterion ``m=`` of ``_criteria`` (PR #69 lead decision 1; GIAP_SPEC.md §8.4).

A replicate that pools ``m`` independent p-values (a draft genome of ``m`` contigs) reports its
rejection share ``s_i = #{p <= alpha} / m``. With ``m=`` the level criterion is

    observed = sum_i s_i m / (R m)  <=  level_bound(alpha, R m)
                                     = alpha + k sqrt(alpha (1 - alpha) / (R m))

so the R m p-values count as independent trials, while the gate stays per replicate: two-block
confirmation (ADR-0015), sharding and ``merge_shards`` see one scalar per replicate. Shares that are
not multiples of 1 / m, or are outside [0, 1], are refused (fail closed); ``m=`` applies to the
level criterion only.

Unmarked tests run on every commit; the false-alarm test over streams is marked ``calibration``.
"""

from __future__ import annotations

import math
from pathlib import Path

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.stats import binom

LEVELS = C.load_tolerances().levels
FAST = settings(max_examples=100, deadline=None)


@pytest.fixture(autouse=True)
def _unsharded(monkeypatch: pytest.MonkeyPatch) -> None:
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV):
        monkeypatch.delenv(var, raising=False)


def _share(alpha: float, m: int):  # type: ignore[no-untyped-def]
    """Replicate: rejection share of m independent p-values whose rejection rate is ``alpha``."""

    def fn(rng: np.random.Generator) -> float:
        return float(rng.binomial(m, min(1.0, alpha))) / m

    return fn


# --------------------------------------------------------------------------- check(..., m=)


def test_pooled_check_small_example() -> None:
    res = C.check("level", [0.1, 0.0, 0.2], 0.05, m=10)  # 3 of 30 p-values rejected
    assert res.observed == pytest.approx(0.1)
    assert res.bound == pytest.approx(C.level_bound(0.05, 30))
    assert res.R == 3  # replicates, not p-values: merge and records key on the suite's R
    assert res.ok is (C.level_bound(0.05, 30) >= 0.1)
    assert "m=10" in res.detail


def test_pooled_bound_is_level_bound_of_r_times_m() -> None:
    shares = [0.02] * 100
    res = C.check("level", shares, 0.01, m=50)
    assert res.bound == pytest.approx(C.level_bound(0.01, 5000), abs=1e-15)
    assert res.bound < C.level_bound(0.01, 100)  # tighter than the per-replicate bound


def test_m_one_counts_indicators_like_check_level() -> None:
    p = np.linspace(0.001, 1.0, 200)
    ref = C.check_level(p, 0.05)
    res = C.check("level", (p <= 0.05).astype(float), 0.05, m=1)
    assert (res.ok, res.observed, res.bound, res.R) == (ref.ok, ref.observed, ref.bound, ref.R)


@FAST
@given(
    m=st.integers(1, 60),
    counts=st.lists(st.integers(0, 60), min_size=1, max_size=40),
    alpha=st.sampled_from(LEVELS),
)
def test_pooled_check_equals_check_level_on_the_pooled_p(m: int, counts: list[int], alpha: float):
    """Property: m= on shares == check_level on the R m pooled indicators (same rate and bound)."""
    counts = [min(c, m) for c in counts]
    res = C.check("level", [c / m for c in counts], alpha, m=m)
    pooled = np.concatenate([[0.0] * c + [1.0] * (m - c) for c in counts])  # 0 <= alpha rejects
    ref = C.check_level(pooled, alpha)
    assert res.observed == pytest.approx(ref.observed, abs=1e-12)
    assert res.bound == pytest.approx(ref.bound, abs=1e-12)
    assert res.ok == ref.ok or abs(ref.observed - ref.bound) < 1e-12
    assert len(counts) == res.R


@pytest.mark.parametrize("m", [0, -3])
def test_m_below_one_refused(m: int) -> None:
    with pytest.raises(ValueError, match="m"):
        C.check("level", [0.0, 0.1], 0.05, m=m)


@pytest.mark.parametrize("criterion", ["uniform", "fdr", "fwer"])
def test_m_only_for_level(criterion: str) -> None:
    with pytest.raises(ValueError, match="level"):
        C.check(criterion, [0.0, 0.1, 0.2], 0.05, m=10)


@pytest.mark.parametrize("bad", [float("nan"), -0.1, 1.1, 0.123])
def test_pooled_check_fails_closed(bad: float) -> None:
    """Non-finite, out-of-range or non-multiple-of-1/m shares are a method bug: refused."""
    with pytest.raises(ValueError):
        C.check("level", [0.1, bad, 0.0], 0.05, m=10)


# --------------------------------------------------------------------------- calibrate(..., m=)


def test_calibrate_records_m_and_uses_pooled_bound(tmp_path: Path) -> None:
    rec = C.calibrate(
        "regions", "level", 0.05, _share(0.05, 50), 100, stream="t.pooled.rec", m=50,
        log_dir=tmp_path,
    )  # fmt: skip
    assert rec["m"] == 50 and rec["R"] == 100 and len(rec["values"]) == 100
    assert rec["bound"] == pytest.approx(C.level_bound(0.05, 5000))
    logged = C.read_results(tmp_path / "regions.jsonl")
    assert all(r["m"] == 50 for r in logged)


def test_calibrate_without_m_records_none(tmp_path: Path) -> None:
    rec = C.calibrate(
        "regions", "level", 0.05, lambda g: g.uniform(), 20, stream="t.pooled.none",
        log_dir=tmp_path,
    )  # fmt: skip
    assert rec.get("m") is None


def test_exact_level_stream_passes(tmp_path: Path) -> None:
    """Exact-level pooled p (rate = alpha) through calibrate_levels: no confirmed fail."""
    m, r = 50, 100

    def fn(rng: np.random.Generator) -> np.ndarray:
        return np.array([_share(a, m)(rng) for a in LEVELS])

    gate = calibrate_levels(
        "regions", "level", LEVELS, fn, r, stream="t.pooled.exact", m=m, log_dir=tmp_path
    )
    assert all(rec["m"] == m for rec in gate.records)
    gate.assert_ok()


def test_one_and_a_half_alpha_stream_fails_confirmed(tmp_path: Path) -> None:
    """Rate 1.5 alpha with R m = 50,000 p-values: a confirmed fail (both blocks) at every alpha."""
    m, r = 500, 100

    def fn(rng: np.random.Generator) -> np.ndarray:
        return np.array([_share(1.5 * a, m)(rng) for a in LEVELS])

    gate = calibrate_levels(
        "regions", "level", LEVELS, fn, r, stream="t.pooled.inflated", m=m, log_dir=tmp_path
    )
    assert [rec["status"] for rec in gate.records] == ["fail"] * len(LEVELS), gate.detail()
    assert all(rec["block"] == 2 for rec in gate.records)


def test_pooled_check_is_not_the_per_replicate_check(tmp_path: Path) -> None:
    """Read as p-values (no m=), the shares of the 1.5 alpha attack pass; m= catches it."""
    m, r, a = 500, 100, 0.05
    rec = C.calibrate(
        "regions", "level", a, _share(1.5 * a, m), r, stream="t.pooled.contrast", m=m,
        log_dir=tmp_path,
    )  # fmt: skip
    assert rec["status"] == "fail"
    assert C.check("level", rec["values"], a).ok


def test_merge_shards_honours_m(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Sharded block-1 records merge to the unsharded verdict, with the pooled bound."""
    m, r, a, stream = 50, 40, 0.05, "t.pooled.shard"
    C.calibrate("regions", "level", a, _share(a, m), r, stream=stream, m=m,
                      log_dir=tmp_path / "ref")  # fmt: skip
    for start, stop in ((0, 15), (15, 40)):
        monkeypatch.setenv(C.SHARD_ENV, f"{start}:{stop}")
        rec = C.calibrate("regions", "level", a, _share(a, m), r, stream=stream, m=m,
                          log_dir=tmp_path / "sh")  # fmt: skip
        assert rec["status"] == C.SHARD_STATUS and rec["m"] == m
    monkeypatch.delenv(C.SHARD_ENV)
    (merged,) = C.merge_shards("regions", log_dir=tmp_path / "sh")
    assert merged["m"] == m
    assert merged["bound"] == pytest.approx(C.level_bound(a, r * m))
    ref1 = C.read_results(tmp_path / "ref" / "regions.jsonl")[0]
    assert ref1["block"] == 1 and merged["observed"] == pytest.approx(ref1["observed"])
    assert merged["status"] in ("pass", C.PENDING_STATUS)


def test_merge_refuses_mixed_m(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    stream = "t.pooled.mixed"
    monkeypatch.setenv(C.SHARD_ENV, "0:10")
    C.calibrate("regions", "level", 0.05, _share(0.05, 10), 20, stream=stream, m=10,
                log_dir=tmp_path)  # fmt: skip
    monkeypatch.setenv(C.SHARD_ENV, "10:20")
    C.calibrate("regions", "level", 0.05, _share(0.05, 20), 20, stream=stream, m=20,
                log_dir=tmp_path)  # fmt: skip
    monkeypatch.delenv(C.SHARD_ENV)
    with pytest.raises(ValueError):
        C.merge_shards("regions", log_dir=tmp_path)


# --------------------------------------------------------------------------- operating charact.

S_STREAMS = 200


def _false_alarm(alpha: float, n: int) -> float:
    """Exact P(one block fails) when n independent trials reject at rate exactly ``alpha``."""
    c = math.floor(C.level_bound(alpha, n) * n + 1e-9)
    return float(binom.sf(c, n, alpha))


@pytest.mark.calibration
def test_pooled_confirmed_false_alarm_matches_theory(tmp_path: Path) -> None:
    """Exact-level pooled p over S streams (R = regions fast, m = 50): confirmed ``fail`` per alpha
    within the 99.9 % binomial band of f^2, ``unconfirmed`` within that of f, where f is the exact
    single-block false alarm of R m binomial trials."""
    r, m = C.replicates("regions", full=False), 50
    for a in LEVELS:
        status = [
            C.calibrate(
                "regions", "level", a, _share(a, m), r, stream=f"t.pooled.oc_{s}", m=m,
                log_dir=tmp_path,
            )["status"]
            for s in range(S_STREAMS)
        ]  # fmt: skip
        f = _false_alarm(a, r * m)
        assert status.count("fail") <= binom.ppf(0.999, S_STREAMS, f * f), (a, status)
        assert status.count("unconfirmed") <= binom.ppf(0.999, S_STREAMS, f), (a, status)
