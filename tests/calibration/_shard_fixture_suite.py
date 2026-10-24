"""A tiny calibration suite for the shard-equivalence teeth (test_criteria_shard_equivalence.py).

Not a calibration claim: tiny fixed R, synthetic p-values, nothing read from ``tolerances.yaml``
but the levels. The leading underscore keeps it out of normal collection; the teeth test names this
file on the pytest command line and sets ``GIAP_SHARD_FIXTURE=1`` (without it the module skips).
It lives under tests/calibration so the real conftest applies (suite keyword, shard-0 rule for
direct tests and teeth, the owning-test record fields, the xfail outcome records).

Every check uses the patterns the real suites use: ``skip_without_verdict`` after the last
harness call, direct parts only where ``owns_direct_checks()``. Direct parts append one line to
``<log_dir>/fixture_direct.jsonl`` so the teeth can count where they ran.
"""

from __future__ import annotations

import json
import os
from collections.abc import Callable
from typing import Any

import _criteria as C
import pytest
from _gate import calibrate_levels

if os.environ.get("GIAP_SHARD_FIXTURE") != "1":
    pytest.skip("shard fixture suite: run by test_criteria_shard_equivalence.py only",
                allow_module_level=True)  # fmt: skip

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("discordance")]
S = "discordance"
FX = "discordance.shard_fixture"
STRICT = pytest.mark.xfail(strict=True, raises=AssertionError, reason="fixture: strict")
LOOSE = pytest.mark.xfail(strict=False, raises=AssertionError, reason="fixture: loose")


def direct(name: str) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "fixture_direct.jsonl").open("a") as fh:
        fh.write(json.dumps({"name": name}) + "\n")


def by_block(stream: str, R: int, b1: float, b2: float) -> Callable[[Any], float]:
    """p = ``b1`` on replicates 0..R-1 and ``b2`` on R..2R-1 (keyed by the replicate seed)."""
    block = {C.suite_seed(stream, i): 1 if i < R else 2 for i in range(2 * R)}
    return lambda rng: b1 if block[int(rng.bit_generator.seed_seq.entropy)] == 1 else b2


def level(name: str, R: int, b1: float, b2: float) -> dict[str, Any]:
    stream = f"{FX}.{name}"
    return C.calibrate(S, "level", 0.05, by_block(stream, R, b1, b2), R, stream=stream)


def verdict(rec: dict[str, Any]) -> None:
    C.skip_without_verdict(rec)
    assert rec["status"] in ("pass", "unconfirmed"), rec


def test_levels_at_every_alpha() -> None:
    gate = calibrate_levels(S, "level", C.load_tolerances().levels,
                            lambda rng: float(rng.uniform()), 7, stream=f"{FX}.levels")  # fmt: skip
    C.skip_without_verdict(gate.records)
    gate.assert_ok()


def test_unconfirmed() -> None:
    verdict(level("unconfirmed", 4, 0.0, 0.5))


def test_unmarked_fail() -> None:
    verdict(level("unmarked_fail", 4, 0.0, 0.0))


@STRICT
def test_strict_xfail_fail() -> None:
    verdict(level("strict_fail", 4, 0.0, 0.0))


@STRICT
def test_strict_xfail_pass() -> None:
    verdict(level("strict_pass", 4, 0.5, 0.5))


@LOOSE
def test_loose_xfail_fail() -> None:
    verdict(level("loose_fail", 5, 0.0, 0.0))


@STRICT
def test_strict_xfail_direct_fails() -> None:
    """The check passes; the direct part (shard 0 only) fails: pytest says XFAIL unsharded."""
    rec = level("strict_direct", 5, 0.5, 0.5)
    if not C.owns_direct_checks():
        C.skip_without_verdict(rec)
    direct("strict_direct")
    raise AssertionError("direct part of a strict-xfail test")


def _counts(i: int) -> dict[str, Any]:
    rng = C.replicate_rng(f"{FX}.counts", i)
    n = int(rng.integers(5, 12))
    f = int(rng.integers(0, 2))
    return {"f1": {"tp_p": n - f, "n_pred": n, "tp_r": n - f, "n_truth": n},
            "fdr": {"n_called": n, "n_false": f}}  # fmt: skip


def test_pooled_f1_with_r_below_n() -> None:
    rec = C.calibrate_pooled(S, f"{FX}.f1", "f1", 0.5, _counts, 2)
    verdict(rec)


def test_fdr_counts() -> None:
    verdict(C.calibrate_fdr_counts(S, f"{FX}.fdrc", "fdr", 0.3, _counts, 5))


def test_direct_only() -> None:
    direct("direct_only")


def test_mixed() -> None:
    rec = level("mixed", 5, 0.5, 0.5)
    if not C.owns_direct_checks():
        C.skip_without_verdict(rec)
    direct("mixed")
    verdict(rec)


def test_teeth(teeth_log: Any) -> None:
    stream = f"{FX}.teeth"
    rec = C.calibrate(S, "level", 0.05, by_block(stream, 4, 0.0, 0.0), 4, stream=stream,
                      log_dir=teeth_log)  # fmt: skip
    direct("teeth")
    assert rec["status"] == "fail", rec


@STRICT
def test_strict_xfail_assert_ok() -> None:
    """A strict xfail on ``gate.assert_ok()`` with no skip: unsharded XFAIL; a shard must not
    XPASS on its no-verdict records (conftest turns that XPASS into a skip, ADR-0057)."""
    gate = calibrate_levels(S, "level", [0.05], lambda _rng: 0.0, 4, stream=f"{FX}.assert_ok")
    gate.assert_ok()


R_PART = 5


def test_partitioned_direct() -> None:
    """A direct test over ``replicate_range`` (as G2 / A6a): every shard checks its own range."""
    for i in C.replicate_range(R_PART):
        direct(f"part{i}")
        assert C.replicate_rng(f"{FX}.part", i).uniform() < 1.0


@STRICT
def test_partitioned_xpass() -> None:
    """Strict xfail, partitioned direct check that is clean: unsharded a strict XPASS; sharded each
    shard skips its clean range, and ``--merge`` must call the test a strict XPASS."""
    bad = [i for i in C.replicate_range(R_PART) if C.replicate_rng(f"{FX}.px", i).uniform() > 1.0]
    if not bad and os.environ.get(C.SHARD_ENV):
        pytest.skip(C.SHARD_SKIP)
    assert not bad, bad
