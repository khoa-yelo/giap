# ruff: noqa: F811  (pytest fixture sb imported from test_m11_03b_launcher)
"""#349 review round of M11-03b, B3: the launcher's R14 job count and nightly window.

Spec: docs/specs/M11-03b-nightly-readiness-labels.md, rule 3 of ``calibration_full_launch.sh`` (as
amended in "Lead decisions, #349 review (2026-10-25)"): a launch queues 16 jobs (15 array tasks
and the aggregator), so it refuses when ``jobs + 16 > 20``; ``T = min(15, floor((128 - cpus - 1) /
4))``, refuse when ``T < 4``; the worst-case end ``now + ceil(15 / T) * 4 h + 20 min`` must be
before the next nightly start (``GIAP_NIGHTLY_AT``, default ``01:30`` local), else refuse. The clock
is ``GIAP_NOW`` (epoch seconds). Every refusal is exit 2, no ``sbatch`` call, and a message naming
the numbers (jobs, CPUs, T, end, nightly start).

Fixtures and helpers come from ``test_m11_03b_launcher`` (sandbox with a fake ``squeue`` printing
``<cpus> <name>`` lines, ``at(h, m)`` = epoch seconds of h:m local time on 2026-10-26).
"""

from __future__ import annotations

import re
import subprocess

import pytest
from helpers.slurm_sandbox import Sandbox, option
from unit.test_m11_03b_launcher import (  # noqa: F401  (fixture sb)
    _array,
    _launch,
    _out,
    _queue,
    at,
    sb,
)


def refused(sb: Sandbox, res: subprocess.CompletedProcess[str]) -> str:
    out = _out(res)
    assert res.returncode == 2, out
    assert sb.sbatch_calls() == [], out
    return out


def accepted(sb: Sandbox, res: subprocess.CompletedProcess[str], throttle: int) -> None:
    assert res.returncode == 0, _out(res)
    assert option(_array(sb)["argv"], "array", "a") == f"0-14%{throttle}"
    assert f"throttle={throttle}" in res.stdout, res.stdout


# ------------------------------------------------------------------ the nightly window


def test_empty_queue_at_22_00_reaches_the_nightly_and_refuses(sb: Sandbox) -> None:
    """T = 15, one wave: 22:00 + 4 h + 20 min = 02:20 >= 01:30."""
    _queue(sb, [])
    out = refused(sb, _launch(sb, GIAP_NOW=at(22)))
    assert "02:20" in out and "01:30" in out, out
    assert re.search(r"\b15\b", out), out  # the throttle


def test_empty_queue_at_10_00_ends_by_14_20(sb: Sandbox) -> None:
    _queue(sb, [])
    accepted(sb, _launch(sb, GIAP_NOW=at(10)), 15)


@pytest.mark.parametrize(
    ("hour", "minute", "ok"),
    [(21, 9, True),    # ends 01:29 < 01:30
     (21, 10, False),  # ends 01:30: "at or after" refuses
     (21, 11, False),
     (1, 0, False),    # the nightly starts 01:30 the same day
     (2, 0, True)],    # after today's start: the next one is tomorrow 01:30
)  # fmt: skip
def test_one_wave_boundary_at_the_default_nightly(
    sb: Sandbox, hour: int, minute: int, ok: bool
) -> None:
    _queue(sb, [])
    res = _launch(sb, GIAP_NOW=at(hour, minute))
    if ok:
        accepted(sb, res, 15)
    else:
        out = refused(sb, res)
        assert "01:30" in out, out


@pytest.mark.parametrize(
    ("lines", "throttle", "hour", "minute", "ok"),
    [(["68 giap-big"], 14, 17, 9, True),    # 2 waves: 17:09 + 8 h 20 = 01:29
     (["68 giap-big"], 14, 17, 10, False),  # 01:30
     (["111 giap-big"], 4, 9, 9, True),     # 4 waves: 09:09 + 16 h 20 = 01:29
     (["111 giap-big"], 4, 9, 10, False),   # 01:30
     (["64 giap-big"], 15, 21, 0, True),    # floor(63/4) = 15: one wave, 01:20
     (["72 giap-big"], 13, 10, 0, True)],   # floor(55/4) = 13: 2 waves, 18:20
)  # fmt: skip
def test_waves_follow_the_throttle(
    sb: Sandbox, lines: list[str], throttle: int, hour: int, minute: int, ok: bool
) -> None:
    _queue(sb, lines)
    res = _launch(sb, GIAP_NOW=at(hour, minute))
    if ok:
        accepted(sb, res, throttle)
    else:
        out = refused(sb, res)
        assert re.search(rf"\b{throttle}\b", out) and "01:30" in out, out


@pytest.mark.parametrize(
    ("nightly", "hour", "minute", "ok"),
    [("14:00", 10, 0, False),  # ends 14:20 >= 14:00
     ("14:20", 10, 0, False),  # at the start
     ("14:21", 10, 0, True),
     ("03:00", 22, 0, True),   # ends 02:20 < 03:00
     ("23:00", 22, 0, False),  # the next start is today 23:00, before the 02:20 end
     ("22:30", 22, 0, False)],
)  # fmt: skip
def test_giap_nightly_at_overrides_the_start(
    sb: Sandbox, nightly: str, hour: int, minute: int, ok: bool
) -> None:
    _queue(sb, [])
    res = _launch(sb, GIAP_NOW=at(hour, minute), GIAP_NIGHTLY_AT=nightly)
    if ok:
        accepted(sb, res, 15)
    else:
        out = refused(sb, res)
        assert nightly in out, out


def test_the_clock_is_giap_now(sb: Sandbox) -> None:
    """Same queue, same wall clock: only GIAP_NOW differs, and so does the outcome."""
    _queue(sb, [])
    assert _launch(sb, GIAP_NOW=at(23)).returncode == 2
    assert sb.sbatch_calls() == []
    accepted(sb, _launch(sb, GIAP_NOW=at(11)), 15)


# ------------------------------------------------------------------ the 16-job count


@pytest.mark.parametrize(
    ("lines", "ok"),
    [(["1 giap-j"] * 4, True),                       # 4 + 16 = 20
     (["1 giap-j"] * 4 + ["2 other"] * 30, True),    # only giap-* jobs count
     (["1 giap-j"] * 5, False),                      # 21 > 20
     (["0 giap-pending"] * 5, False),                # pending tasks count, CPUs or not
     (["1 giap-j"] * 10, False)],                    # the review's 10 + 16 = 26 case
)  # fmt: skip
def test_a_launch_counts_its_16_jobs(sb: Sandbox, lines: list[str], ok: bool) -> None:
    _queue(sb, lines)
    res = _launch(sb, GIAP_NOW=at(10))
    if ok:
        accepted(sb, res, 15)
    else:
        out = refused(sb, res)
        n = sum(ln.split()[1].startswith("giap-") for ln in lines)
        assert re.search(rf"\b{n}\b", out) and re.search(r"\b16\b", out), out
        assert re.search(r"\b20\b", out), out


def test_the_throttle_does_not_reserve_the_nightly_cpus(sb: Sandbox) -> None:
    """No standing reserve (B3 replaces Q3): 60 CPUs queued still allow T = 15 at 10:00."""
    _queue(sb, ["60 giap-other"])
    accepted(sb, _launch(sb, GIAP_NOW=at(10)), 15)


def test_cpu_refusal_names_the_numbers(sb: Sandbox) -> None:
    _queue(sb, ["112 giap-big"])  # floor(15/4) = 3 < 4
    out = refused(sb, _launch(sb, GIAP_NOW=at(10)))
    assert re.search(r"\b112\b", out) and re.search(r"\b3\b", out), out
