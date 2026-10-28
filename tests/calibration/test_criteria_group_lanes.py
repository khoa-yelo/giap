"""M11-03c PR 2: ``shard0_only`` and ``--calibration-lane-report`` (conftest.py), and group-lane
equivalence (docs/specs/M11-03c-calibration-full-sharding.md, "Group selection uses the harness's
own rule", acceptance criterion 8, statistician S6).

- ``shard0_only(item)`` = ``calibration_unshardable`` marker, or the ``teeth_log`` fixture, or no
  harness entry point / ``replicate_range`` reached; ``pytest_runtest_setup`` calls it, so its
  shard-0-only skips in shard K > 0 are exactly the items the helper flags.
- ``--calibration-lane-report=<path>`` writes ``{nodeid: {"shard0_only": bool}}`` for the
  collected items and changes nothing else.
- Group lanes (unsharded runs of disjoint test subsets via ``--deselect``, each with its own
  log directory) give the same per-check records and JUnit outcomes for G = 1, 2, 3; and the
  ``direct`` layout (shard lanes without the shard-0-only tests, those in group lanes, block 2 on
  lanes) gives the unsharded verdicts.

Runs on ``_shard_fixture_suite.py`` (every pattern of the real suites: direct-only, teeth,
unshardable, mixed, partitioned direct, strict / loose xfails). Unmarked: every commit.
"""

from __future__ import annotations

import importlib.util
import inspect
import json
import os
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any

import _criteria as C
import pytest

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "_shard_fixture_suite.py"
SUITE = "discordance"
#: the ``ok`` scenario of test_criteria_shard_equivalence.py (exit 0 unsharded and merged)
OK = "not unmarked_fail and not strict_xfail_pass and not partitioned_xpass and not rr_"
#: shard-0-only tests of the whole fixture file, by today's three rules (read off the file)
SHARD0_ONLY = {
    "test_direct_only",  # no harness call, no replicate_range
    "test_teeth",  # teeth_log
    "test_unshardable",  # calibration_unshardable
    "test_rr_shard0_assert_fails",  # no harness call
    "test_rr_shard0_strict_xpass",  # no harness call
}
#: the shard-0-only skip reason (conftest constant) of each flagged test in a shard K > 0
SHARD0_REASON = {"test_teeth": "TEETH_SKIP", "test_unshardable": "UNSHARDABLE_SKIP"}
LANE_VARS = ("GIPPA_CALIBRATION_CONFIRM_SLICE", "GIPPA_CALIBRATION_CONFIRM_LANES",
             "GIPPA_CALIBRATION_CONFIRM_JUNIT")  # fmt: skip


def _env(log: Path, extra: dict[str, str] | None = None) -> dict[str, str]:
    drop = (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV, *LANE_VARS, "PYTEST_ADDOPTS",
            "PYTEST_XDIST_WORKER")  # fmt: skip
    keep = {k: v for k, v in os.environ.items() if k not in drop}
    return {**keep, C.LOG_DIR_ENV: str(log), "GIPPA_SHARD_FIXTURE": "1", **(extra or {})}


def _args(selection: str | None) -> list[str]:
    out = [str(FIXTURE), "-p", "no:cacheprovider", "-q", "-m", "calibration"]
    return out if selection is None else [*out, "-k", selection]


def _pytest(args: list[str], env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    cmd = [sys.executable, "-m", "pytest", *args]
    return subprocess.run(cmd, env=env, cwd=C.REPO, check=False, capture_output=True, text=True)


def _lane_report(tmp: Path, selection: str | None = None) -> dict[str, bool]:
    path = tmp / "lane_report.json"
    res = _pytest([*_args(selection), "--collect-only", f"--calibration-lane-report={path}"],
                  _env(tmp / "log"))  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr
    raw = json.loads(path.read_text())
    assert all(set(v) == {"shard0_only"} and isinstance(v["shard0_only"], bool)
               for v in raw.values()), raw  # fmt: skip
    return {k: v["shard0_only"] for k, v in raw.items()}


def _collected(
    tmp: Path, selection: str | None = None, extra: list[str] | None = None
) -> list[str]:
    res = _pytest([*_args(selection), "--collect-only", *(extra or [])], _env(tmp / "log"))
    assert res.returncode == 0, res.stdout + res.stderr
    return sorted(line.strip() for line in res.stdout.splitlines() if "::" in line)


def _junit(path: Path) -> dict[str, tuple[str, str]]:
    """``{nodeid-ish: (outcome, message)}``; outcome passed / skipped / failure / error."""
    out = {}
    for tc in ET.parse(path).getroot().iter("testcase"):
        kinds = [c for c in tc if c.tag in ("failure", "error", "skipped")]
        outcome = kinds[0].tag if kinds else "passed"
        msg = kinds[0].get("message", "") if kinds else ""
        out[f"{tc.get('classname')}::{tc.get('name')}"] = (outcome, msg)
    return out


def _name(nodeid: str) -> str:
    return nodeid.rsplit("::", 1)[1]


# --------------------------------------------------------------------------- shard0_only


def _load_conftest() -> Any:
    spec = importlib.util.spec_from_file_location("_m11_03c_conftest_probe", HERE / "conftest.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_shard0_only_exists_and_pytest_runtest_setup_calls_it() -> None:
    """One rule, two users (lane report, runtest setup): they cannot drift."""
    mod = _load_conftest()
    assert callable(getattr(mod, "shard0_only", None))
    assert list(inspect.signature(mod.shard0_only).parameters) == ["item"]
    assert "shard0_only(" in inspect.getsource(mod.pytest_runtest_setup)


def test_lane_report_flags_exactly_todays_shard0_only_tests(tmp_path: Path) -> None:
    report = _lane_report(tmp_path)
    assert sorted(report) == _collected(tmp_path)  # every collected item, nothing else
    assert {_name(n) for n, flag in report.items() if flag} == SHARD0_ONLY


def test_lane_report_changes_nothing_else(tmp_path: Path) -> None:
    """Same collection with and without the option; a normal run with it gives the same
    outcomes as without."""
    path = tmp_path / "r.json"
    assert _collected(tmp_path, OK, [f"--calibration-lane-report={path}"]) == _collected(
        tmp_path, OK
    )
    a, b = tmp_path / "a", tmp_path / "b"
    ra = _pytest([*_args(OK), f"--junitxml={a}.xml"], _env(a))
    rb = _pytest([*_args(OK), f"--junitxml={b}.xml",
                  f"--calibration-lane-report={tmp_path / 'r2.json'}"], _env(b))  # fmt: skip
    assert (ra.returncode, _junit(Path(f"{a}.xml"))) == (rb.returncode, _junit(Path(f"{b}.xml")))


@pytest.mark.parametrize("shard", ["1/3", "2/3", "0/3"])
def test_shard0_only_agrees_with_runtest_setup_skips(shard: str, tmp_path: Path) -> None:
    """In a block-1 shard K > 0, the tests ``pytest_runtest_setup`` skips as shard-0-only are
    exactly the ``shard0_only`` items; in shard 0 none is skipped for that reason."""
    report = _lane_report(tmp_path)
    k, n = (int(x) for x in shard.split("/"))
    rng = C.shard_replicates(10, n, k)
    junit = tmp_path / "shard.xml"
    env = _env(tmp_path / "log", {C.SHARD_ENV: f"{rng.start}:{rng.stop}@{k}/{n}",
                                  C.BLOCK_ENV: "1"})  # fmt: skip
    _pytest([*_args(None), f"--junitxml={junit}"], env)
    outcomes = _junit(junit)
    assert {_name(t) for t in outcomes} == {_name(t) for t in report}
    skipped0 = {_name(t) for t, (o, msg) in outcomes.items()
                if o == "skipped" and "shard 0 only" in msg}  # fmt: skip
    flagged = {_name(t) for t, flag in report.items() if flag}
    assert skipped0 == (flagged if k > 0 else set()), (skipped0, flagged)
    # and each for its own reason (review S3): swapping two skip messages is caught
    mod = _load_conftest()
    want = {name: SHARD0_REASON.get(name, "SHARD_DIRECT_SKIP") for name in flagged}
    got = {_name(t): msg for t, (o, msg) in outcomes.items() if _name(t) in skipped0}
    assert got.keys() == (want.keys() if k > 0 else set())
    for name, const in want.items() if k > 0 else ():
        reason = getattr(mod, const)
        assert got[name] == reason or got[name].endswith(f": {reason}"), (name, const, got[name])


# --------------------------------------------------------------------------- group lanes


def _addopts(deselect: list[str], junit: Path | None = None) -> str:
    opts = [f"--deselect={n}" for n in deselect]
    if junit is not None:
        opts.append(f"--junitxml={junit}")
    return " ".join(shlex.quote(o) for o in opts)


def _groups(nodeids: list[str], G: int) -> list[list[str]]:
    return [nodeids[g::G] for g in range(G)]


def _no_prefix_pairs(nodeids: list[str]) -> None:
    """``--deselect`` is a prefix match: the fixture must hold no proper-prefix pair."""
    bad = [(a, b) for a in nodeids for b in nodeids if a != b and b.startswith(a)]
    assert not bad, bad


def _records(logs: list[Path]) -> dict[tuple[Any, ...], tuple[Any, ...]]:
    out: dict[tuple[Any, ...], tuple[Any, ...]] = {}
    for log in logs:
        if not (log / f"{SUITE}.jsonl").exists():
            continue
        for r in C.merge_shards(SUITE, log):
            key = (r["stream"], r["criterion"], r["level"], r["R"], str(r.get("key")))
            assert key not in out, key  # every check in exactly one lane
            out[key] = (r["status"], r["block"], r["observed"], r["bound"], r["seed"],
                        json.dumps(r["values"]), r["test"])  # fmt: skip
    return out


def _direct(logs: list[Path]) -> Counter[str]:
    out: Counter[str] = Counter()
    for log in logs:
        path = log / "fixture_direct.jsonl"
        if path.exists():
            out += Counter(r["name"] for r in C.read_results(path))
    return out


Outcomes = dict[str, tuple[str, str]]


def _run_groups(
    tmp: Path, nodeids: list[str], G: int, universe: list[str] | None = None
) -> tuple[list[int], list[Path], Outcomes]:
    """Group lane g deselects every collected test (``universe``, default ``nodeids``) not in g."""
    rcs: list[int] = []
    logs: list[Path] = []
    outcomes: Outcomes = {}
    for g, mine in enumerate(_groups(nodeids, G)):
        log = tmp / "lanes" / f"{SUITE}-group{g}"
        junit = tmp / f"task-5.group{g}.junit.xml"
        others = [n for n in (universe or nodeids) if n not in mine]
        env = _env(log, {"PYTEST_ADDOPTS": _addopts(others, junit)})
        res = _pytest(_args(OK), env)
        rcs.append(res.returncode)
        logs.append(log)
        got = _junit(junit)
        assert {_name(t) for t in got} == {_name(t) for t in mine}, (g, got)  # expected set ran
        outcomes |= got
    return rcs, logs, outcomes


@pytest.mark.parametrize("G", [2, 3])
def test_group_lanes_give_the_single_process_records(G: int, tmp_path: Path) -> None:
    nodeids = _collected(tmp_path, OK)
    _no_prefix_pairs(nodeids)
    rcs1, logs1, out1 = _run_groups(tmp_path / "g1", nodeids, 1)
    rcs, logs, out = _run_groups(tmp_path / f"g{G}", nodeids, G)
    assert rcs1 == [0] and rcs == [0] * G, (rcs1, rcs)
    assert _records(logs) == _records(logs1) and _records(logs1)  # values, status, order aside
    assert _direct(logs) == _direct(logs1)
    assert out == out1  # same JUnit outcome per test


def test_direct_layout_shard_lanes_plus_group_lanes_equal_unsharded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Plan row ``2 shards, 2 groups, direct``: the shard-0-only tests (lane report) go to group
    lanes, the shard lanes and the merge deselect them, block 2 runs on 2 lanes. Verdicts, direct
    parts and exits equal the unsharded run's."""
    report = _lane_report(tmp_path, OK)
    group_tests = sorted(n for n, flag in report.items() if flag)
    assert {_name(n) for n in group_tests} == SHARD0_ONLY & {"test_direct_only", "test_teeth",
                                                             "test_unshardable"}  # fmt: skip
    one = tmp_path / "one"
    rc_one = _pytest(_args(OK), _env(one)).returncode

    shards = tmp_path / "log"
    monkeypatch.setattr(C, "pytest_args", lambda _s, _f: _args(OK))
    for k, v in _env(shards, {"PYTEST_ADDOPTS": _addopts(group_tests)}).items():
        monkeypatch.setenv(k, v)
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV, *LANE_VARS):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.chdir(C.REPO)
    rcs = [C.main(["--suite", SUITE, "--fast", "--shard", f"{k}/2"]) for k in range(2)]
    monkeypatch.setenv("GIPPA_CALIBRATION_CONFIRM_LANES", "2")
    monkeypatch.setenv("GIPPA_CALIBRATION_CONFIRM_JUNIT", str(shards / "task-5.merge"))
    rc_merge = C.main(["--suite", SUITE, "--fast", "--merge"])

    g_rcs, g_logs, _out = _run_groups(tmp_path, group_tests, 2, universe=sorted(report))
    assert (rc_one, rcs, rc_merge, g_rcs) == (0, [0, 0], 0, [0, 0])

    ref = _records([one])
    merged = {}
    for r in C.read_results(shards / f"{SUITE}.merged.jsonl"):
        key = (r["stream"], r["criterion"], r["level"], r["R"], str(r.get("key")))
        merged[key] = (r["status"], r["block"], r["observed"], r["bound"], r["seed"],
                       json.dumps(r["values"]), r["test"])  # fmt: skip
    groups = _records(g_logs)
    assert not merged.keys() & groups.keys()  # a check runs in the shard lanes or a group lane
    assert merged | groups == ref
    assert _direct([shards, *g_logs]) == _direct([one])  # each direct part ran once
    parts = json.loads((shards / f"{SUITE}.confirm_parts.json").read_text())
    assert sorted(parts) == ["task-5.merge.lane0.junit.xml", "task-5.merge.lane1.junit.xml"]
