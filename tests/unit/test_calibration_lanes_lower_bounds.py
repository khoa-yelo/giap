"""Probe-log lower bounds across several logs (backlog 2026-10-27, #445 review nits N1, N3).

- N1: ``durations --probe-log`` with several logs keeps each test's largest lower bound, whatever
  the order of the logs; under ``--max`` a measured value replaces a ``:lower-bound`` row even when
  it is smaller, and a lower bound replaces a row only when it is larger.
- N3 (mutant M8): a test finished in log A and unfinished in log B gets A's measured value, not
  B's (larger) lower bound.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "ci" / "calibration_lanes.py"
X = "tests/calibration/test_t.py::test_x"
Y = "tests/calibration/test_t.py::test_y"
Z = "tests/calibration/test_t.py::test_z"

# A stamped line's true time is the previous stamp (xdist ends a line at the next event).
LOG_A = (  # X starts at 1000 and passes at 1003: measured 3 s
    "start 999\n"
    "1000 ============================= test session starts ==============================\n"
    f"1001 {X} \n"
    "1003 scheduling\n"
    f"1004 [gw0] [100%] PASSED {X} \n"
    "1004 ============================== 1 passed ===============================\n"
)
LOG_B = (  # X starts at 2000, Y at 2001; neither finishes; last stamp 2500
    "2000 ============================= test session starts ==============================\n"
    f"2001 {X} \n"
    f"2002 {Y} \n"
    "2500 slurmstepd: *** JOB CANCELLED DUE TO TIME LIMIT ***\n"
)
LOG_C = (  # Y starts at 3000, never finishes; last stamp 3030
    "3000 ============================= test session starts ==============================\n"
    f"3001 {Y} \n"
    "3030 slurmstepd: *** JOB CANCELLED DUE TO TIME LIMIT ***\n"
)


@pytest.fixture(scope="module")
def mod() -> ModuleType:
    name = "calibration_lanes_lower_bounds"
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    assert spec is not None and spec.loader is not None
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


@pytest.fixture
def logs(tmp_path: Path) -> dict[str, Path]:
    out = {}
    for key, text in {"a": LOG_A, "b": LOG_B, "c": LOG_C}.items():
        out[key] = tmp_path / f"gippa-probe-{key}.1.out"
        out[key].write_text(text)
    return out


def _cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def _rows(path: Path) -> dict[str, list[str]]:
    return {ln.split("\t")[0]: ln.split("\t")[1:] for ln in path.read_text().splitlines()[1:]}


def test_fixture_logs_parse_as_documented(mod: ModuleType, logs: dict[str, Path]) -> None:
    assert mod.durations_from_probe_log(logs["a"]) == {X: 3.0}
    assert mod.probe_log_lower_bounds(logs["a"]) == {}
    assert mod.durations_from_probe_log(logs["b"]) == {}
    assert mod.probe_log_lower_bounds(logs["b"]) == {X: 500.0, Y: 499.0}
    assert mod.probe_log_lower_bounds(logs["c"]) == {Y: 30.0}


def test_merge_lower_bounds_keeps_the_largest_in_any_order(mod: ModuleType) -> None:
    b, c = {X: 500.0, Y: 499.0}, {Y: 30.0}
    assert mod.merge_lower_bounds([b, c]) == {X: 500.0, Y: 499.0}
    assert mod.merge_lower_bounds([c, b]) == {X: 500.0, Y: 499.0}
    assert mod.merge_lower_bounds([{Y: 1.0}, {Y: 0.5}, {Y: 2.0}]) == {Y: 2.0}
    assert mod.merge_lower_bounds([]) == {}


@pytest.mark.parametrize("order", [("a", "b", "c"), ("c", "b", "a"), ("b", "a", "c")])
def test_cli_measured_beats_bound_and_bounds_keep_the_max(
    logs: dict[str, Path], tmp_path: Path, order: tuple[str, ...]
) -> None:
    out = tmp_path / "d.tsv"
    res = _cli("durations", "--source", "p", "--probe-log", "--out", str(out),
               *(str(logs[k]) for k in order))  # fmt: skip
    assert res.returncode == 0, res.stderr
    rows = _rows(out)
    assert rows[X] == ["3.0", "p"]  # finished in A: measured, not B's 500 s bound (M8)
    assert rows[Y] == ["499.0", "p:lower-bound"]  # B's 499 s, not C's 30 s (N1)


def test_cli_max_measured_replaces_a_larger_lower_bound_row(logs: dict[str, Path],
                                                            tmp_path: Path) -> None:  # fmt: skip
    out = tmp_path / "d.tsv"
    out.write_text("nodeid\tseconds\tsource\n"
                   f"{X}\t500.0\told:lower-bound\n"
                   f"{Y}\t900.0\told:lower-bound\n"
                   f"{Z}\t10.0\told\n")  # fmt: skip
    res = _cli("durations", "--source", "p", "--probe-log", "--max", "--out", str(out),
               str(logs["a"]), str(logs["c"]))  # fmt: skip
    assert res.returncode == 0, res.stderr
    rows = _rows(out)
    assert rows[X] == ["3.0", "p"]  # a measurement beats a bound, even a larger one
    assert rows[Y] == ["900.0", "old:lower-bound"]  # a smaller bound never lowers a bound
    assert rows[Z] == ["10.0", "old"]  # untouched


def test_cli_max_measured_replaces_a_smaller_lower_bound_row(logs: dict[str, Path],
                                                             tmp_path: Path) -> None:  # fmt: skip
    out = tmp_path / "d.tsv"
    out.write_text(f"nodeid\tseconds\tsource\n{X}\t1.0\told:lower-bound\n")
    res = _cli("durations", "--source", "p", "--probe-log", "--max", "--out", str(out),
               str(logs["a"]))  # fmt: skip
    assert res.returncode == 0, res.stderr
    assert _rows(out)[X] == ["3.0", "p"]  # a larger measurement replaces a bound too


def test_cli_max_bounds_only_raise(logs: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "d.tsv"
    out.write_text("nodeid\tseconds\tsource\n"
                   f"{X}\t600.0\told\n"
                   f"{Y}\t100.0\told:lower-bound\n")  # fmt: skip
    res = _cli("durations", "--source", "p", "--probe-log", "--max", "--out", str(out),
               str(logs["b"]))  # fmt: skip
    assert res.returncode == 0, res.stderr
    rows = _rows(out)
    assert rows[X] == ["600.0", "old"]  # a smaller bound never replaces a measurement
    assert rows[Y] == ["499.0", "p:lower-bound"]  # a larger bound raises a bound


def test_cli_append_never_replaces(logs: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "d.tsv"
    out.write_text(f"nodeid\tseconds\tsource\n{X}\t500.0\told:lower-bound\n")
    res = _cli("durations", "--source", "p", "--probe-log", "--append", "--out", str(out),
               str(logs["a"]))  # fmt: skip
    assert res.returncode == 0, res.stderr
    assert _rows(out)[X] == ["500.0", "old:lower-bound"]  # --append only adds new tests
