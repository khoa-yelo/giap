"""``scripts/ci/calibration_rows.py`` as a command (issue M11-03b, examples 4 and 9; acceptance
criterion 7; CLI exit codes).

The log directory is written by ``helpers.calibration_rows_fixture.World.write`` (records
``task-<i>.json`` plus the JUnit files their ``parts`` list); the TSV lives in ``tmp_path``.

Convention settled by the test author (the spec leaves the file layout open): ``--sacct-json FILE``
holds a JSON object mapping the array task index (a string) to its sacct state, e.g.
``{"7": "OUT_OF_MEMORY"}``; the state is added to the ``no record`` reason of that task.
"""

from __future__ import annotations

import fcntl
import json
import re
import subprocess
import sys
import time
from pathlib import Path

from helpers.calibration_rows_fixture import (
    JOB,
    SCRIPT,
    C,
    World,
    by_label,
    record,
    run_cli,
    split_row,
)
from helpers.release_readiness import load_module


def _args(log: Path, tsv: Path, *extra: str) -> list[str]:
    return ["--log-dir", str(log), "--job", JOB, "--commit", C, "--expect", "15",
            "--tsv", str(tsv), *extra]  # fmt: skip


def _lines(tsv: Path) -> list[str]:
    return tsv.read_text().splitlines(keepends=True) if tsv.is_file() else []


def test_cli_writes_two_passing_rows(tmp_path: Path) -> None:
    log, tsv = World.default().write(tmp_path / "log"), tmp_path / "nightly.tsv"
    proc = run_cli(*_args(log, tsv))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    lines = _lines(tsv)
    assert len(lines) == 2, lines
    for line in lines:
        date, job, commit, ex, _s = split_row(line)
        assert re.fullmatch(r"\d{4}-\d\d-\d\d", date) and job == f"{JOB}_calibration"
        assert commit == C and ex == 0
    got = by_label(lines)
    assert got["calibration_full"][1].startswith("calibration_full: 15 suites, 0 failed suites;")
    rows, bad = load_module().parse_log(tsv.read_text())
    assert bad == 0 and all(r.ok() for r in rows)


def test_cli_appends_to_an_existing_tsv(tmp_path: Path) -> None:
    log, tsv = World.default().write(tmp_path / "log"), tmp_path / "nightly.tsv"
    old = "2026-11-01\t1800000_2\tabcdef12\t0\tcalibration_nightly: 1 passed in 1.00s\n"
    tsv.write_text(old)
    assert run_cli(*_args(log, tsv)).returncode == 0
    lines = _lines(tsv)
    assert lines[0] == old and len(lines) == 3, lines


def test_example9_second_run_on_the_same_log_dir_writes_nothing(tmp_path: Path) -> None:
    log, tsv = World.default().write(tmp_path / "log"), tmp_path / "nightly.tsv"
    assert run_cli(*_args(log, tsv)).returncode == 0
    assert (log / "rows.done").exists()
    first = tsv.read_text()
    proc = run_cli(*_args(log, tsv))
    assert proc.returncode == 0, proc.stderr
    assert tsv.read_text() == first


def test_append_waits_for_the_tsv_lock(tmp_path: Path) -> None:
    """Both rows are appended under ``flock`` on ``<tsv>.lock``."""
    log, tsv = World.default().write(tmp_path / "log"), tmp_path / "nightly.tsv"
    lock = tmp_path / "nightly.tsv.lock"
    with open(lock, "a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        cmd = [sys.executable, str(SCRIPT), *_args(log, tsv)]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            time.sleep(2.0)
            assert proc.poll() is None, "the writer did not wait for the lock"
            assert _lines(tsv) == [], "rows appended while another process held the lock"
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)
        out, err = proc.communicate(timeout=60)
    assert proc.returncode == 0, out + err
    assert len(_lines(tsv)) == 2


def test_cli_exit_1_when_one_row_fails(tmp_path: Path) -> None:
    w = World.default()
    w.records[3] = record(3, exit_=1)
    w.set_outcome(3, attack=False, outcome="failed")
    log, tsv = w.write(tmp_path / "log"), tmp_path / "nightly.tsv"
    proc = run_cli(*_args(log, tsv))
    assert proc.returncode == 1, proc.stdout + proc.stderr
    got = by_label(_lines(tsv))
    assert got["calibration_full"][0] == 1 and got["calibration_attacks"][0] == 0


def test_example4_missing_record_names_the_suite(tmp_path: Path) -> None:
    w = World.default()
    w.records[6] = None  # discordance (task 6; the spec table says 7)
    log, tsv = w.write(tmp_path / "log"), tmp_path / "nightly.tsv"
    proc = run_cli(*_args(log, tsv))
    assert proc.returncode == 1, proc.stdout + proc.stderr
    for ex, s in by_label(_lines(tsv)).values():
        assert ex == 1 and "1 failed suites (discordance no record" in s, s


def test_example4_sacct_state_is_named(tmp_path: Path) -> None:
    w = World.default()
    w.records[6] = None  # discordance (task 6; the spec table says 7)
    log, tsv = w.write(tmp_path / "log"), tmp_path / "nightly.tsv"
    sacct = tmp_path / "sacct.json"
    sacct.write_text(json.dumps({str(t): "COMPLETED" for t in range(15)} | {"6": "OUT_OF_MEMORY"}))
    proc = run_cli(*_args(log, tsv, "--sacct-json", str(sacct)))
    assert proc.returncode == 1, proc.stdout + proc.stderr
    for ex, s in by_label(_lines(tsv)).values():
        assert ex == 1 and "(discordance no record, OUT_OF_MEMORY)" in s, s


def test_missing_junit_file_fails_both_rows(tmp_path: Path) -> None:
    log, tsv = World.default().write(tmp_path / "log"), tmp_path / "nightly.tsv"
    (log / "task-4.run.junit.xml").unlink()
    assert run_cli(*_args(log, tsv)).returncode == 1
    assert all(ex == 1 for ex, _s in by_label(_lines(tsv)).values())


def test_unparsable_junit_or_record_fails_both_rows(tmp_path: Path) -> None:
    log, tsv = World.default().write(tmp_path / "log"), tmp_path / "nightly.tsv"
    (log / "task-4.run.junit.xml").write_text("<testsuites><testsuite>")
    (log / "task-9.json").write_text('{"task": 9, "exit"')
    assert run_cli(*_args(log, tsv)).returncode == 1
    for ex, s in by_label(_lines(tsv)).values():
        assert ex == 1 and "2 failed suites (" in s, s


def test_empty_log_dir_gives_15_failed_suites(tmp_path: Path) -> None:
    log, tsv = tmp_path / "log", tmp_path / "nightly.tsv"
    log.mkdir()
    proc = run_cli(*_args(log, tsv))
    assert proc.returncode == 1, proc.stdout + proc.stderr
    ex, s = by_label(_lines(tsv))["calibration_full"]
    assert ex == 1 and "15 failed suites (" in s, s


def test_unreadable_log_dir_exits_2_and_still_writes_failing_rows(tmp_path: Path) -> None:
    tsv = tmp_path / "nightly.tsv"
    proc = run_cli(*_args(tmp_path / "no-such-dir", tsv))
    assert proc.returncode == 2, proc.stdout + proc.stderr
    got = by_label(_lines(tsv))
    assert all(ex == 1 for ex, _s in got.values()), got


def test_expect_other_than_the_suite_count_is_a_usage_error(tmp_path: Path) -> None:
    log, tsv = World.default().write(tmp_path / "log"), tmp_path / "nightly.tsv"
    args = _args(log, tsv)
    args[args.index("--expect") + 1] = "14"
    assert run_cli(*args).returncode == 2


def test_bad_commit_is_a_usage_error(tmp_path: Path) -> None:
    log, tsv = World.default().write(tmp_path / "log"), tmp_path / "nightly.tsv"
    args = _args(log, tsv)
    args[args.index("--commit") + 1] = "abc1234"  # not 40-hex
    assert run_cli(*args).returncode == 2
