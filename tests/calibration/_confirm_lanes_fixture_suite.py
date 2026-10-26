"""A tiny calibration suite for the block-2-on-lanes teeth (test_criteria_confirm_lanes.py).

M11-03c PR 2 (docs/specs/M11-03c-calibration-full-sharding.md, criterion 7). Not a calibration
claim: tiny fixed R per check (2, 4, 5, 7: uneven slices for L = 2, 3 and an empty slice for
R = 2, L = 3), synthetic p-values and count records keyed by the absolute replicate index, so
block 1 (0..R-1) and block 2 (R..2R-1) can be made to fail or pass on purpose. The leading
underscore keeps it out of normal collection; the teeth name this file on the pytest command line
and set ``GIAP_CONFIRM_LANES_FIXTURE=1`` (without it the module skips). It lives under
tests/calibration so the real conftest applies (suite keyword, owning-test records, xfail
outcomes, the confirm-run skips).

Every check fails block 1 unless its name says otherwise, so ``--merge`` must run its block 2:
``calibrate`` (level), ``calibrate_pooled`` (pooled F1) and ``calibrate_fdr_counts``. The
``vacuous`` check fails block 1 too but is vacuous (rule 6), so it must never reach block 2.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Mapping
from typing import Any

import _criteria as C
import pytest

if os.environ.get("GIAP_CONFIRM_LANES_FIXTURE") != "1":
    pytest.skip("confirm-lanes fixture suite: run by test_criteria_confirm_lanes.py only",
                allow_module_level=True)  # fmt: skip

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("discordance")]
S = "discordance"
FX = "discordance.confirm_lanes_fixture"
STRICT = pytest.mark.xfail(strict=True, raises=AssertionError, reason="fixture: strict")
F1_BOUND = 0.5
Q = 0.3


def _p_values(stream: str, R: int, bad1: bool, bad2: bool) -> Callable[[Any], float]:
    """p in [0, 0.04] (rejects at 0.05) on a bad block, in [0.06, 1] on a good one; every value
    is distinct (drawn from the replicate's own generator), so a wrong seed changes ``values``."""
    block = {C.suite_seed(stream, i): 1 if i < R else 2 for i in range(2 * R)}

    def fn(rng: Any) -> float:
        bad = bad1 if block[int(rng.bit_generator.seed_seq.entropy)] == 1 else bad2
        u = float(rng.uniform())
        return 0.04 * u if bad else 0.06 + 0.94 * u

    return fn


def _records(stream: str, R: int, bad1: bool, bad2: bool) -> Callable[[int], Mapping[str, Any]]:
    """Count records of absolute replicate ``i``: F1 below 0.5 and FDP above 0.5 on a bad block,
    F1 and FDP near 1 and 0 on a good one."""

    def rec(i: int) -> Mapping[str, Any]:
        rng = C.replicate_rng(stream, i)
        n = int(rng.integers(8, 15))
        bad = bad1 if i < R else bad2
        f = int(rng.integers(n // 2 + 2, n)) if bad else int(rng.integers(0, 2))
        tp = n - f
        return {"f1": {"tp_p": tp, "n_pred": n, "tp_r": tp, "n_truth": n},
                "fdr": {"n_called": n, "n_false": f, "n_targets": n}}  # fmt: skip

    return rec


def level(name: str, R: int, *, bad1: bool = True, bad2: bool) -> dict[str, Any]:
    stream = f"{FX}.{name}"
    return C.calibrate(S, "level", 0.05, _p_values(stream, R, bad1, bad2), R, stream=stream)


def pooled(name: str, R: int, *, bad2: bool) -> dict[str, Any]:
    stream = f"{FX}.{name}"
    return C.calibrate_pooled(S, stream, "f1", F1_BOUND, _records(stream, R, True, bad2), R)


def fdrc(name: str, R: int, *, bad2: bool, vacuous: Mapping[str, float] | None = None) -> Any:
    stream = f"{FX}.{name}"
    return C.calibrate_fdr_counts(S, stream, "fdr", Q, _records(stream, R, True, bad2), R,
                                  vacuous=vacuous)  # fmt: skip


def verdict(rec: Mapping[str, Any]) -> None:
    C.skip_without_verdict(rec)
    assert rec["status"] in ("pass", "unconfirmed"), rec


# --------------------------------------------------------------------------- scenario "ok"


def test_pass_b1() -> None:
    """Block 1 passes: never pending, never in a block-2 lane."""
    verdict(level("pass_b1", 5, bad1=False, bad2=True))


def test_cal_unconfirmed() -> None:
    verdict(level("cal_unconfirmed", 5, bad2=False))


def test_pooled_unconfirmed() -> None:
    verdict(pooled("pooled_unconfirmed", 5, bad2=False))


def test_fdrc_unconfirmed() -> None:
    verdict(fdrc("fdrc_unconfirmed", 7, bad2=False))


@STRICT
def test_strict_fail_cal() -> None:
    """Confirmed fail under a strict xfail: XFAIL, ``--merge`` exit 0 (ADR-0057)."""
    verdict(level("strict_fail_cal", 4, bad2=True))


@STRICT
def test_strict_fail_pooled() -> None:
    verdict(pooled("strict_fail_pooled", 4, bad2=True))


@STRICT
def test_strict_fail_fdrc() -> None:
    verdict(fdrc("strict_fail_fdrc", 5, bad2=True))


def test_vacuous() -> None:
    """Block 1 is bad but vacuous (mean n_called below the floor): no block 2, anywhere."""
    rec = fdrc("vacuous", 5, bad2=True, vacuous={"n_called": 1000.0})
    C.skip_without_verdict(rec)
    assert rec["status"] == C.VACUOUS_STATUS, rec


def test_r2_cal() -> None:
    """R = 2: with L = 3 block-2 lanes, lane 2's slice is empty (it must log nothing)."""
    verdict(level("r2_cal", 2, bad2=False))


def test_r2_pooled() -> None:
    verdict(pooled("r2_pooled", 2, bad2=False))


@STRICT
def test_r2_fdrc_strict_fail() -> None:
    verdict(fdrc("r2_fdrc_strict_fail", 2, bad2=True))


# --------------------------------------------------------------------------- scenario "fail"


def test_unmarked_fail() -> None:
    """Both blocks fail, no mark: FAILED unsharded, ``--merge`` MERGE_FAIL."""
    verdict(level("unmarked_fail", 4, bad2=True))


def test_unmarked_fdrc_fail() -> None:
    verdict(fdrc("unmarked_fdrc_fail", 5, bad2=True))


# --------------------------------------------------------------------------- scenario "xpass"


@STRICT
def test_xpass_strict_cal() -> None:
    """Block 1 fails, block 2 passes (``unconfirmed``) under a strict xfail: XPASS(strict), exit 1
    unsharded and at ``--merge``, however block 2 is split."""
    verdict(level("xpass_strict_cal", 5, bad2=False))


@STRICT
def test_xpass_strict_fdrc() -> None:
    verdict(fdrc("xpass_strict_fdrc", 7, bad2=False))
