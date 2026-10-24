"""Teeth of sharding: a small suite run unsharded and as 3 shards + ``--merge`` gives the same
verdicts, the same pooled rates and the same exit (R3: sharding moves replicates between
processes, it never changes a record). The suite is ``_shard_fixture_suite.py`` (tiny R, every
pattern of the real suites: multi-level gates, an unconfirmed check that needs ``--merge``'s
block-2 confirm run, pooled F1 with R below N, FDR counts, direct-only and mixed tests, a teeth
test, and strict / non-strict xfails of ADR-0057).

Driven through ``_criteria.main`` exactly as ``calibration_full.sbatch`` drives it
(``--suite discordance --fast --shard K/3``, then ``--merge``), with ``pytest_args`` pointed at the
fixture file. Scenarios (``-k``) give the four exit cases: ``ok`` (strict-xfail fail, loose xfail:
exit 0 both ways), ``xpass`` (strict-xfail pass: non-zero both ways), ``fail`` (unmarked confirmed
fail: non-zero both ways). Unmarked: runs on every commit (about 30 s).
"""

from __future__ import annotations

import os
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import _criteria as C
import pytest

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "_shard_fixture_suite.py"
SUITE = "discordance"
N = 3
SCENARIOS = {
    "ok": "not unmarked_fail and not strict_xfail_pass",
    "xpass": "strict_xfail_pass or levels",
    "fail": "unmarked_fail or levels",
}


def _args(selection: str) -> list[str]:
    return [str(FIXTURE), "-p", "no:cacheprovider", "-q", "-m", "calibration", "-k", selection]


def _env(log: Path) -> dict[str, str]:
    keep = {k: v for k, v in os.environ.items() if k not in (C.SHARD_ENV, C.BLOCK_ENV,
                                                            C.CONFIRM_ENV, "PYTEST_XDIST_WORKER")}  # fmt: skip
    return {**keep, C.LOG_DIR_ENV: str(log), "GIAP_SHARD_FIXTURE": "1"}


def _unsharded(log: Path, selection: str) -> int:
    cmd = [sys.executable, "-m", "pytest", *_args(selection)]
    return subprocess.run(cmd, env=_env(log), cwd=C.REPO, check=False).returncode


def _sharded(log: Path, selection: str, monkeypatch: pytest.MonkeyPatch) -> tuple[list[int], int]:
    monkeypatch.setattr(C, "pytest_args", lambda _suite, _full: _args(selection))
    for k, v in _env(log).items():
        monkeypatch.setenv(k, v)
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.chdir(C.REPO)
    shards = [C.main(["--suite", SUITE, "--fast", "--shard", f"{k}/{N}"]) for k in range(N)]
    return shards, C.main(["--suite", SUITE, "--fast", "--merge"])


def _verdicts(records: list[dict[str, Any]]) -> dict[tuple[Any, ...], tuple[Any, ...]]:
    return {(r["stream"], r["criterion"], r["level"], r["R"], str(r.get("key"))):
            (r["status"], r["block"], r["observed"], r["bound"]) for r in records}  # fmt: skip


def _direct(log: Path) -> Counter[str]:
    path = log / "fixture_direct.jsonl"
    return Counter(r["name"] for r in C.read_results(path)) if path.exists() else Counter()


@pytest.mark.parametrize("scenario", sorted(SCENARIOS))
def test_three_shards_and_merge_equal_the_unsharded_run(
    scenario: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sel = SCENARIOS[scenario]
    one, many = tmp_path / "one", tmp_path / "many"
    rc_one = _unsharded(one, sel)
    shards, rc_merge = _sharded(many, sel, monkeypatch)

    ref = _verdicts(C.merge_shards(SUITE, one))
    got = _verdicts(C.read_results(many / f"{SUITE}.merged.jsonl"))
    assert got == ref  # same checks, verdicts, deciding block and pooled rates
    assert ref, "the fixture logged no check"
    assert _direct(many) == _direct(one)  # every direct part ran exactly once, as unsharded

    want = {"ok": 0, "xpass": 1, "fail": 1}[scenario]
    assert rc_one == want, rc_one
    assert all(rc == 0 for rc in shards), shards  # shards have no verdict of their own
    assert rc_merge == (0 if want == 0 else C.MERGE_FAIL), rc_merge


def test_ok_scenario_reports_the_expected_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _sharded(tmp_path, SCENARIOS["ok"], monkeypatch)
    merged = {r["stream"]: r for r in C.read_results(tmp_path / f"{SUITE}.merged.jsonl")}
    fx = "discordance.shard_fixture"
    assert merged[f"{fx}.strict_fail"]["outcome"] == "xfail"
    assert merged[f"{fx}.loose_fail"]["outcome"] == "xfail"
    assert merged[f"{fx}.strict_direct"]["outcome"] == "xfail"  # shard 0's direct failure
    assert merged[f"{fx}.unconfirmed"]["status"] == "unconfirmed"
    assert merged[f"{fx}.unconfirmed"].get("outcome") is None
