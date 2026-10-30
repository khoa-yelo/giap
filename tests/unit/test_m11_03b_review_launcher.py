# ruff: noqa: F811  (pytest fixture sb imported from test_m11_03b_launcher)
"""#349 review round of M11-03b, B3: the launcher's R14 job count and nightly window.

Spec: docs/specs/M11-03b-nightly-readiness-labels.md, rule 3 of ``calibration_full_launch.sh``, as
replaced by M11-03c PR 4 (docs/specs/M11-03c-calibration-full-sharding.md): a launch queues 16 jobs
(15 per-suite jobs and the aggregator), so it refuses when ``jobs + 16 > 20``; it refuses when
``cpus + Σ + 2 > 144`` (Σ = 138 allocated CPUs; 128 is a soft cap that only warns, Q-K59); the
worst-case end ``now + max(time) + 20 min`` (8 h 20 m) must be
before the next nightly start (``GIPPA_NIGHTLY_AT``, default ``01:30`` local), else refuse. The
clock
is ``GIPPA_NOW`` (epoch seconds). Every refusal is exit 2, no ``sbatch`` call, and a message naming
the numbers (jobs, CPUs, end, nightly start).

Fixtures and helpers come from ``test_m11_03b_launcher`` (sandbox with a fake ``squeue`` printing
``<cpus> <name>`` lines, ``at(h, m)`` = epoch seconds of h:m local time on 2026-10-26).
"""

from __future__ import annotations

import re
import subprocess

import pytest
from helpers.slurm_sandbox import Sandbox
from unit.test_m11_03b_launcher import (  # noqa: F401  (fixture sb)
    _launch,
    _out,
    _queue,
    _suites,
    at,
    sb,
)


def refused(sb: Sandbox, res: subprocess.CompletedProcess[str]) -> str:
    out = _out(res)
    assert res.returncode == 2, out
    assert sb.sbatch_calls() == [], out
    return out


def accepted(sb: Sandbox, res: subprocess.CompletedProcess[str]) -> None:
    assert res.returncode == 0, _out(res)
    assert len(_suites(sb)) == 15
    assert "cpus=138" in res.stdout, res.stdout


# ------------------------------------------------------------------ the nightly window


def test_empty_queue_at_22_00_reaches_the_nightly_and_refuses(sb: Sandbox) -> None:
    """22:00 + 8 h + 20 min = 06:20 >= 01:30."""
    _queue(sb, [])
    out = refused(sb, _launch(sb, GIPPA_NOW=at(22)))
    assert "06:20" in out and "01:30" in out, out


def test_empty_queue_at_10_00_ends_by_18_20(sb: Sandbox) -> None:
    _queue(sb, [])
    accepted(sb, _launch(sb, GIPPA_NOW=at(10)))


@pytest.mark.parametrize(
    ("hour", "minute", "ok"),
    [(17, 9, True),    # ends 01:29 < 01:30
     (17, 10, False),  # ends 01:30: "at or after" refuses
     (17, 11, False),
     (1, 0, False),    # the nightly starts 01:30 the same day
     (2, 0, True)],    # after today's start: the next one is tomorrow 01:30
)  # fmt: skip
def test_window_boundary_at_the_default_nightly(
    sb: Sandbox, hour: int, minute: int, ok: bool
) -> None:
    _queue(sb, [])
    res = _launch(sb, GIPPA_NOW=at(hour, minute))
    if ok:
        accepted(sb, res)
    else:
        out = refused(sb, res)
        assert "01:30" in out, out


@pytest.mark.parametrize(
    ("nightly", "hour", "minute", "ok"),
    [("18:00", 10, 0, False),  # ends 18:20 >= 18:00
     ("18:20", 10, 0, False),  # at the start
     ("18:21", 10, 0, True),
     ("06:21", 22, 0, True),   # ends 06:20 < 06:21
     ("23:00", 22, 0, False),  # the next start is today 23:00, before the 06:20 end
     ("22:30", 22, 0, False)],
)  # fmt: skip
def test_gippa_nightly_at_overrides_the_start(
    sb: Sandbox, nightly: str, hour: int, minute: int, ok: bool
) -> None:
    _queue(sb, [])
    res = _launch(sb, GIPPA_NOW=at(hour, minute), GIPPA_NIGHTLY_AT=nightly)
    if ok:
        accepted(sb, res)
    else:
        out = refused(sb, res)
        assert nightly in out, out


def test_the_clock_is_gippa_now(sb: Sandbox) -> None:
    """Same queue, same wall clock: only GIPPA_NOW differs, and so does the outcome."""
    _queue(sb, [])
    assert _launch(sb, GIPPA_NOW=at(23)).returncode == 2
    assert sb.sbatch_calls() == []
    accepted(sb, _launch(sb, GIPPA_NOW=at(11)))


# ------------------------------------------------------------------ the 16-job count


@pytest.mark.parametrize(
    ("lines", "ok"),
    [(["0 gippa-j"] * 4, True),                       # 4 + 16 = 20 (and 0 + 140 = 140 CPUs)
     (["0 gippa-j"] * 4 + ["2 other"] * 30, True),    # only gippa-* jobs count
     (["1 gippa-j"] * 5, False),                      # 21 > 20
     (["0 gippa-pending"] * 5, False),                # pending tasks count, CPUs or not
     (["1 gippa-j"] * 10, False)],                    # the review's 10 + 16 = 26 case
)  # fmt: skip
def test_a_launch_counts_its_16_jobs(sb: Sandbox, lines: list[str], ok: bool) -> None:
    _queue(sb, lines)
    res = _launch(sb, GIPPA_NOW=at(10))
    if ok:
        accepted(sb, res)
    else:
        out = refused(sb, res)
        n = sum(ln.split()[1].startswith("gippa-") for ln in lines)
        assert re.search(rf"\b{n}\b", out) and re.search(r"\b16\b", out), out
        assert re.search(r"\b20\b", out), out


def test_no_standing_reserve_for_the_nightly_cpus(sb: Sandbox) -> None:
    """No standing reserve (B3 replaces Q3): an empty queue fits (0 + 138 + 2 = 140 <= 144); five
    pending CPUs (three whole cores, 6) do not (146 > 144)."""
    _queue(sb, [])
    accepted(sb, _launch(sb, GIPPA_NOW=at(10)))
    _queue(sb, ["5 gippa-other"])
    refused_res = _launch(sb, GIPPA_NOW=at(10))
    assert refused_res.returncode == 2, _out(refused_res)


def test_cpu_refusal_names_the_numbers(sb: Sandbox) -> None:
    _queue(sb, ["112 gippa-big"])  # 112 + 138 + 2 = 252 > 144
    out = refused(sb, _launch(sb, GIPPA_NOW=at(10)))
    assert all(re.search(rf"\b{n}\b", out) for n in (112, 144, 252)), out
