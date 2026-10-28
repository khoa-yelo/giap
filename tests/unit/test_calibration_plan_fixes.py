# ruff: noqa: F811  (pytest fixture sb imported from test_m11_03b_launcher)
"""Full-calibration plan fixes after dry run 45bc3599 (jobs 1900336-1900351, 2026-10-27).

- The plan is sized for wheat (Skylake 6148, 2.8-3.2x slower than barley): placement, regions,
  autocorrelation, determinants get 08:00:00 (``MAX_TIME_S``), gene_fdr 3 groups, autocorrelation
  3 x 2 CPUs, determinants 12 shards + 3 ``direct`` groups. Every row passes ``read_plan``.
- R14 counts the CPUs Slurm allocates: ``max(lanes x cpus_per_lane, ceil(lanes x mem / 4000))``
  (partition MaxMemPerCPU 4000 MB); discordance asks 2 x 16000 MB and gets 8 CPUs, not 4. The plan
  allocates 126 CPUs, 127 with the rows job; the launcher sums the same.
- The ⅔ rule (docs/specs/M11-03c-calibration-full-sharding.md, "Initial plan"): each suite job's
  LPT estimate from ``calibration_durations.tsv`` is at most ⅔ of its ``time``.
- A probe-log test that started and never finished gets a lower bound (last stamp minus start).
- A nameless ``<testcase time="0.001"/>`` (an interrupted lane) is never counted as a test.
"""

from __future__ import annotations

import importlib.util
import json
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from types import ModuleType

import pytest
import yaml
from helpers.slurm_sandbox import REPO, Sandbox, option
from unit.test_m11_03b_launcher import (  # noqa: F401  (fixture sb)
    _launch,
    _out,
    _queue,
    sb,
)
from unit.test_m11_03c_lanes import PROBE_LOG

PLAN = REPO / "scripts" / "slurm" / "calibration_plan.tsv"
DURATIONS = REPO / "scripts" / "slurm" / "calibration_durations.tsv"
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"
SUITE_SBATCH = REPO / "scripts" / "slurm" / "calibration_suite.sbatch"
R14_CPUS = 128

# suite: (shards, groups, group_select, cpus_per_lane, mem_per_lane_mb, time)
CHANGED = {
    "placement": (2, 8, "direct", 2, 4000, "08:00:00"),
    "gene_fdr": (0, 3, "*", 2, 4000, "06:00:00"),
    "regions": (0, 5, "*", 2, 4000, "08:00:00"),
    "autocorrelation": (0, 3, "*", 2, 4000, "08:00:00"),
    "determinants": (12, 3, "direct", 2, 4000, "08:00:00"),
}

# The 5-group regions row is the lead's plan (cal_timeouts.md: its longest wheat lane ran 5.8 h).
# With the recorded durations (max of the barley/probe rows and run 1900336) its LPT estimate is
# about 5.5 h, over ⅔ of 8 h (5.33 h); 6 groups (4.6 h) would need 2 CPUs the R14 budget lacks.
OVER_TWO_THIRDS = {"regions"}


def _load(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def lanes_mod() -> ModuleType:
    return _load("calibration_lanes_fixes", REPO / "scripts" / "ci" / "calibration_lanes.py")


@pytest.fixture(scope="module")
def rows_mod() -> ModuleType:
    return _load("calibration_rows_fixes", REPO / "scripts" / "ci" / "calibration_rows.py")


def _suites() -> list[str]:
    with TOLERANCES.open() as fh:
        return [str(k) for k in yaml.safe_load(fh)["replicates_full"]]


# ------------------------------------------------------------------ plan rows and R14


def test_every_plan_row_passes_read_plan(lanes_mod: ModuleType) -> None:
    rows = lanes_mod.read_plan(PLAN, _suites())
    assert [r.suite for r in rows] == _suites() and len(rows) == 15
    for r in rows:
        assert 0 < r.time_s <= lanes_mod.MAX_TIME_S == 8 * 3600, r
        if r.suite in CHANGED:
            got = (r.shards, r.groups, r.group_select, r.cpus_per_lane, r.mem_per_lane_mb, r.time)
            assert got == CHANGED[r.suite], r


def test_discordance_gets_the_cpus_its_memory_needs(lanes_mod: ModuleType) -> None:
    rows = {r.suite: r for r in lanes_mod.read_plan(PLAN, _suites())}
    disc = rows["discordance"]
    assert (disc.shards + disc.groups, disc.mem_per_lane_mb) == (2, 16000)
    assert disc.requested_cpus == 4
    assert lanes_mod.allocated_cpus(disc) == 8  # 32000 MB / MaxMemPerCPU 4000
    assert lanes_mod.allocated_cpus(disc, 16000) == 4
    assert lanes_mod.allocated_cpus(rows["missed_operons"]) == 2  # 8000 MB: exactly 2 CPUs
    for r in rows.values():
        lanes = r.shards + r.groups
        want = max(lanes * r.cpus_per_lane, -(-lanes * r.mem_per_lane_mb // 4000))
        assert lanes_mod.allocated_cpus(r) == want, r


def test_planned_allocated_cpus_with_the_rows_job_fit_r14(lanes_mod: ModuleType) -> None:
    rows = lanes_mod.read_plan(PLAN, _suites())
    requested = sum(r.requested_cpus for r in rows)
    allocated = sum(lanes_mod.allocated_cpus(r) for r in rows)
    assert (requested, allocated) == (122, 126)
    assert allocated + 1 == 127 <= R14_CPUS
    assert len(rows) + 1 <= 20  # jobs


def test_launcher_sums_allocated_cpus(sb: Sandbox) -> None:
    """1 queued CPU + 126 allocated + 1 rows job = 128 launches; 2 queued (129) refuses."""
    _queue(sb, ["1 gippa-other"])
    res = _launch(sb)
    assert res.returncode == 0, _out(res)
    assert re.search(r"\bcpus=126\b", res.stdout), res.stdout
    jobs = {option(c["argv"], "job-name", "J"): c["argv"] for c in sb.sbatch_calls()}
    disc = jobs["gippa-cal-discordance"]
    assert option(disc, "cpus-per-task", "c") == "4", disc  # the request itself is unchanged
    assert re.fullmatch(r"32000M?", option(disc, "mem") or ""), disc


def test_launcher_refuses_when_allocated_cpus_pass_r14(sb: Sandbox) -> None:
    _queue(sb, ["2 gippa-other"])
    res = _launch(sb)
    assert res.returncode == 2 and sb.sbatch_calls() == [], _out(res)
    assert re.search(r"\b129\b", _out(res)) and "allocated" in _out(res), _out(res)


# ------------------------------------------------------------------ the ⅔ rule

_SUITE_PLUGIN = """
import json, os

def pytest_collection_finish(session):
    out = {}
    for item in session.items:
        names = [k for k in getattr(item, "extra_keyword_matches", ()) if k.startswith("calsuite_")]
        if names:
            out[item.nodeid] = names[0][len("calsuite_"):]
    with open(os.environ["GIPPA_TEST_SUITE_MAP"], "w") as fh:
        json.dump(out, fh)
"""


@pytest.fixture(scope="module")
def collection(tmp_path_factory: pytest.TempPathFactory) -> dict[str, dict[str, bool]]:
    """``suite -> {nodeid: shard0_only}`` of the full tier, from one ``--collect-only``."""
    d = tmp_path_factory.mktemp("collect")
    (d / "gippa_suite_map_plugin.py").write_text(_SUITE_PLUGIN)
    report, suites = d / "report.json", d / "suites.json"
    env_path = [str(d), str(REPO / "src"), str(REPO / "tests" / "sim"), str(REPO / "tests")]
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(env_path),
           "GIPPA_TEST_SUITE_MAP": str(suites)}  # fmt: skip
    env.pop("PYTEST_ADDOPTS", None)
    res = subprocess.run(
        [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "-p", "gippa_suite_map_plugin",
         "-p", "no:xdist", "tests/calibration", "-m", "calibration_full", "--collect-only", "-q",
         f"--calibration-lane-report={report}"],
        cwd=REPO, env=env, capture_output=True, text=True, timeout=600,
    )  # fmt: skip
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    flags = {n: bool(v["shard0_only"]) for n, v in json.loads(report.read_text()).items()}
    out: dict[str, dict[str, bool]] = {}
    for n, s in json.loads(suites.read_text()).items():
        out.setdefault(s, {})[n] = flags[n]
    return out


def _estimate_s(mod: ModuleType, row: object, report: dict[str, bool],
                durations: dict[str, float]) -> float:  # fmt: skip
    """The job's estimate: its longest group lane (LPT) or a shard lane's share of the shard tests
    (their recorded single-process time over ``shards``); a missing test weighs ``unknown_min``."""
    lanes = mod.lanes(row, report, durations)

    def w(n: str) -> float:
        return float(durations.get(n, row.unknown_min * 60))  # type: ignore[attr-defined]

    loads = [sum(w(n) for n in ln.expected) for ln in lanes if ln.kind == "group"]
    shard = [ln for ln in lanes if ln.kind == "shard"]
    if shard:
        loads.append(sum(w(n) for n in shard[0].expected) / len(shard))
    return max(loads)


@pytest.mark.parametrize(
    "suite",
    [
        pytest.param(s, marks=pytest.mark.xfail(strict=True, reason="regions 5 x 8 h: see "
                                                 "OVER_TWO_THIRDS")) if s in OVER_TWO_THIRDS
        else s
        for s in _suites()
    ],
)  # fmt: skip
def test_lpt_estimate_within_two_thirds_of_time(lanes_mod: ModuleType,
                                                collection: dict[str, dict[str, bool]],
                                                suite: str) -> None:  # fmt: skip
    rows = {r.suite: r for r in lanes_mod.read_plan(PLAN, _suites())}
    durations = lanes_mod.read_durations(DURATIONS)
    row = rows[suite]
    est = _estimate_s(lanes_mod, row, collection[suite], durations)
    assert est <= 2 / 3 * row.time_s, (suite, est / 3600, row.time)


def test_durations_rows_are_the_max_of_old_and_new(lanes_mod: ModuleType) -> None:
    """Re-recorded with ``--max`` from run 1900336: tests it measured slower than the old rows
    were raised (the 7 ``arm_b`` autocorrelation tests never finished there and stay unknown)."""
    d = lanes_mod.read_durations(DURATIONS)
    lines = DURATIONS.read_text().splitlines()[1:]
    sources = {ln.split("\t")[0]: ln.split("\t")[2] for ln in lines}
    assert "1900336" in set(sources.values())
    raised = "tests/calibration/attacks/test_attack_determinants_null.py::" \
             "test_null_level_under_attack[full-detect-elastic_net]"  # fmt: skip
    assert d[raised] > 2398.0 and sources[raised] == "1900336"  # probe-1895353 had 2398 s


# ------------------------------------------------------------------ probe-log lower bounds


def test_unfinished_probe_tests_get_a_lower_bound(lanes_mod: ModuleType, tmp_path: Path) -> None:
    log = tmp_path / "gippa-probe-x.1.out"
    log.write_text(PROBE_LOG.split("1792998939 [gw1]")[0])  # [6] and [3] never finish
    f = "tests/calibration/test_s.py::test_s"
    lb = lanes_mod.probe_log_lower_bounds(log)
    assert lb == {f"{f}[6]": 4.0, f"{f}[3]": 3.0}  # last stamp 1792998937 minus the starts
    assert lb[f"{f}[6]"] <= 6.0 and lb[f"{f}[3]"] <= 3.0  # the true durations (sleep 6 / 3)
    log.write_text(PROBE_LOG)
    assert lanes_mod.probe_log_lower_bounds(log) == {}


def _cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(REPO / "scripts" / "ci" / "calibration_lanes.py"),
                           *args], capture_output=True, text=True)  # fmt: skip


def test_durations_cli_writes_lower_bounds_marked(lanes_mod: ModuleType, tmp_path: Path) -> None:
    log = tmp_path / "gippa-probe-x.1.out"
    log.write_text(PROBE_LOG.split("1792998939 [gw1]")[0])
    out = tmp_path / "d.tsv"
    res = _cli("durations", "--source", "probe-9", "--probe-log", "--out", str(out), str(log))
    assert res.returncode == 0, res.stderr
    f = "tests/calibration/test_s.py::test_s"
    rows = {ln.split("\t")[0]: ln.split("\t")[1:] for ln in out.read_text().splitlines()[1:]}
    assert rows[f"{f}[6]"] == ["4.0", "probe-9:lower-bound"]
    assert rows[f"{f}[1]"] == ["1.0", "probe-9"]
    assert lanes_mod.read_durations(out)[f"{f}[3]"] == 3.0  # the planner weighs it, not unknown


def test_durations_cli_max_keeps_the_larger(lanes_mod: ModuleType, tmp_path: Path) -> None:
    out = tmp_path / "d.tsv"
    f = "tests/calibration/test_s.py::test_s"
    out.write_text(f"nodeid\tseconds\tsource\n{f}[1]\t0.5\told\n{f}[2]\t9.0\told\n"
                   f"{f}[6]\t2.0\told\n")  # fmt: skip
    log = tmp_path / "gippa-probe-x.1.out"
    log.write_text(PROBE_LOG.split("1792998939 [gw1]")[0])
    res = _cli("durations", "--source", "p", "--probe-log", "--max", "--out", str(out), str(log))
    assert res.returncode == 0, res.stderr
    rows = {ln.split("\t")[0]: ln.split("\t")[1:] for ln in out.read_text().splitlines()[1:]}
    assert rows[f"{f}[1]"] == ["1.0", "p"]  # raised
    assert rows[f"{f}[2]"] == ["9.0", "old"]  # larger old row kept
    assert rows[f"{f}[6]"] == ["4.0", "p:lower-bound"]  # a lower bound above the old row
    assert rows[f"{f}[3]"] == ["3.0", "p:lower-bound"]  # added
    both = _cli("durations", "--source", "p", "--append", "--max", "--out", str(out), str(log))
    assert both.returncode == 2, both.stderr  # one mode only


# ------------------------------------------------------------------ nameless testcases

PHANTOM = ('<testsuites><testsuite name="pytest" tests="3">'
           '<testcase classname="tests.calibration.test_a" name="test_x[full]" time="5.0"/>'
           '<testcase classname="" name="tests/calibration/test_b.py" time="0.1">'
           '<error message="collection failure">boom</error></testcase>'
           '<testcase time="0.001"/>'
           '</testsuite></testsuites>')  # fmt: skip


def test_nameless_testcase_is_not_counted_by_the_rows(rows_mod: ModuleType) -> None:
    cases = rows_mod._testcases([ET.fromstring(PHANTOM)])
    assert ("", "") not in cases and len(cases) == 2, cases
    assert cases[("tests.calibration.test_a", "test_x[full]")] == "passed"
    assert cases[("", "tests/calibration/test_b.py")] == "errors"  # a collection error still counts
    only = rows_mod._testcases([ET.fromstring('<testsuite><testcase time="0.001"/></testsuite>')])
    assert only == {}


def test_nameless_testcase_is_not_counted_by_the_suite_job(tmp_path: Path) -> None:
    """``calibration_suite.sbatch`` counts JUnit testcases for ``verify_counts``: 2, not 3."""
    text = SUITE_SBATCH.read_text()
    m = re.search(r"cat > \"\$JOBPY\" <<'PY'\n(.*?)\nPY\n", text, re.S)
    assert m, "suite_job.py heredoc not found"
    body = m.group(1).split("\nCMDS = ")[0]  # the helpers, without the command dispatch
    (tmp_path / "suite_job_fixes.py").write_text(body)
    job = _load("suite_job_fixes", tmp_path / "suite_job_fixes.py")
    junit = tmp_path / "task-1.group1.junit.xml"
    junit.write_text(PHANTOM)
    assert job._junit_count(junit) == 2
