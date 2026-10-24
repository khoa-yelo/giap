"""Teeth of sharding: a small suite run unsharded and as 3 shards + ``--merge`` gives the same
verdicts, the same pooled rates and the same exit (R3: sharding moves replicates between
processes, it never changes a record). The suite is ``_shard_fixture_suite.py`` (tiny R, every
pattern of the real suites: multi-level gates, an unconfirmed check that needs ``--merge``'s
block-2 confirm run, pooled F1 with R below N, FDR counts, direct-only and mixed tests, a teeth
test, and strict / non-strict xfails of ADR-0057).

Driven through ``_criteria.main`` exactly as ``calibration_full.sbatch`` drives it
(``--suite discordance --fast --shard K/3``, then ``--merge``), with ``pytest_args`` pointed at the
fixture file. Scenarios (``-k``) give the exit cases: ``ok`` (strict-xfail fail, loose xfail, a
strict xfail on ``assert_ok``: exit 0 both ways), ``xpass`` (strict-xfail check passes: non-zero
both ways), ``direct_xpass`` (strict-xfail partitioned direct test is clean: non-zero both ways),
``fail`` (unmarked confirmed fail: non-zero both ways). Review of PR #327 adds ``shard0_fail`` (an
unmarked direct assertion fails in shard 0: shard 0 exits 1 and so must the merge, B2),
``shard0_xpass`` (a strict-xfail test without harness checks passes in shard 0: strict XPASS
there, never turned into a skip, and the merge exits 1), ``raises_other`` / ``loose_xpass``
(``raises`` and ``strict=False`` reach the records) and ``mixed_xfail`` / ``mixed_xpass`` (a test
whose shard-0 statuses are [pass, shard]: the shard defers, the merge decides). The fail-closed
shard-exit cases (usage error, collection error, ``-k`` selecting nothing) are at the end.
Unmarked: runs on every commit (about 90 s).
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
    "ok": "not unmarked_fail and not strict_xfail_pass and not partitioned_xpass and not rr_",
    "xpass": "strict_xfail_pass or levels",
    "direct_xpass": "partitioned_xpass or levels",
    "fail": "unmarked_fail or levels",
    # review of PR #327
    "shard0_fail": "rr_shard0_assert_fails or levels",  # B2: verdict only in shard 0's exit
    "shard0_xpass": "rr_shard0_strict_xpass or levels",  # B2 + mutant (a)
    "raises_other": "rr_raises_other or levels",  # mutant (c): raises= is honoured
    "loose_xpass": "rr_loose_xpass or levels",  # mutant (d): strict=False is honoured
    "mixed_xfail": "rr_mixed_xfail or levels",  # mutant (b): one check without a verdict defers
    "mixed_xpass": "rr_mixed_xpass or levels",
}
#: unsharded (= merged) exit of each scenario: 0, or non-zero (merge: MERGE_FAIL)
WANT = {"ok": 0, "xpass": 1, "direct_xpass": 1, "fail": 1, "shard0_fail": 1, "shard0_xpass": 1,
        "raises_other": 1, "loose_xpass": 0, "mixed_xfail": 0, "mixed_xpass": 1}  # fmt: skip
#: shard exits per scenario (default: all 0, a shard has no verdict of its own; None: any). A
#: direct assertion (or a strict XPASS of a test without harness checks) is judged in shard 0
#: only, so shard 0 exits 1 there, and ``--merge`` must fold that exit in (B2).
SHARD_RCS: dict[str, list[int] | None] = {"shard0_fail": [1, 0, 0], "shard0_xpass": [1, 0, 0],
                                         "mixed_xpass": None}  # fmt: skip


def _args(selection: str) -> list[str]:
    return [str(FIXTURE), "-p", "no:cacheprovider", "-q", "-m", "calibration", "-k", selection]


def _env(log: Path) -> dict[str, str]:
    drop = (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV, "PYTEST_XDIST_WORKER")
    keep = {k: v for k, v in os.environ.items() if k not in drop}
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

    want = WANT[scenario]
    assert rc_one == want, rc_one
    rcs = SHARD_RCS.get(scenario, [0] * N)
    if rcs is not None:
        assert shards == rcs, shards
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
    assert merged[f"{fx}.assert_ok"]["outcome"] == "xfail"  # no shard XPASSed it


def test_partitioned_direct_test_runs_every_replicate_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A direct test over ``replicate_range`` runs in every shard (each on its own range)."""
    _sharded(tmp_path, "partitioned_direct", monkeypatch)
    assert _direct(tmp_path) == Counter(f"part{i}" for i in range(5))


# --------------------------------------------------------------------------- B1 (review of #327)
# Shards whose pytest exits non-zero without logging a calibrate record must never merge to 0
# (origin/main: MERGE_ERROR "no calibrate records"). Exit 1 (a test failed) is MERGE_FAIL; any
# other non-zero exit (usage error 4, collection error 2, nothing selected 5) is MERGE_ERROR.


def test_usage_error_in_every_shard_is_a_merge_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PYTEST_ADDOPTS", "--no-such-option")
    shards, rc_merge = _sharded(tmp_path, SCENARIOS["ok"], monkeypatch)
    assert shards == [4] * N, shards
    assert rc_merge == C.MERGE_ERROR, rc_merge


def test_nothing_selected_is_a_merge_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    shards, rc_merge = _sharded(tmp_path, "rr_no_such_test_selected", monkeypatch)
    assert shards == [5] * N, shards
    assert rc_merge == C.MERGE_ERROR, rc_merge


def test_collection_error_is_a_merge_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bad = tmp_path / "suite" / "test_broken_import.py"
    bad.parent.mkdir()
    bad.write_text("import giap_no_such_module_xyz  # noqa: F401\n\n\ndef test_x() -> None:\n"
                   "    pass\n")  # fmt: skip
    for k, v in _env(tmp_path / "log").items():
        monkeypatch.setenv(k, v)
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.chdir(C.REPO)
    monkeypatch.setattr(C, "pytest_args", lambda _s, _f: [str(bad), "-p", "no:cacheprovider", "-q"])
    shards = [C.main(["--suite", SUITE, "--fast", "--shard", f"{k}/{N}"]) for k in range(N)]
    assert all(rc not in (0, 1) for rc in shards), shards
    assert C.main(["--suite", SUITE, "--fast", "--merge"]) == C.MERGE_ERROR


def test_direct_failure_of_a_suite_without_records_is_a_merge_fail(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A denovo_precision-like suite (no harness check at all): its only verdict is shard 0's exit
    (1); the "no calibrate records; every shard finished" branch must not return 0."""
    shards, rc_merge = _sharded(tmp_path, "rr_shard0_assert_fails", monkeypatch)
    assert not (tmp_path / f"{SUITE}.jsonl").exists()
    assert shards == [1, 0, 0], shards
    assert rc_merge == C.MERGE_FAIL, rc_merge
