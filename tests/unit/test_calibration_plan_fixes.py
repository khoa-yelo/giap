# ruff: noqa: F811  (pytest fixture sb imported from test_m11_03b_launcher)
"""Full-calibration plan fixes after dry run 45bc3599 (jobs 1900336-1900351, 2026-10-27).

- The plan is sized for wheat (Skylake 6148, 2.8-3.2x slower than barley): placement, regions,
  autocorrelation, determinants get 08:00:00 (``MAX_TIME_S``), gene_fdr 3 groups, autocorrelation
  3 x 2 CPUs, determinants 12 shards + 3 ``direct`` groups. Every row passes ``read_plan``.
- R14 counts the CPUs Slurm allocates: ``max(lanes x cpus_per_lane, ceil(lanes x mem / 4000))``
  (partition MaxMemPerCPU 4000 MB); discordance asks 2 x 16000 MB and gets 8 CPUs, not 4. The plan
  allocates 126 CPUs, 127 with the rows job; the launcher sums the same.
- The ⅔ rule (docs/specs/M11-03c-calibration-full-sharding.md, "Initial plan"): each suite job's
  LPT estimate from ``calibration_durations.tsv`` is at most ⅔ of its ``time``. A shard lane's
  estimate includes the merge, which runs in the same job time (backlog 2026-10-27).
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
        assert lanes_mod.allocated_cpus(r) == want + want % 2, r


@pytest.mark.parametrize(("lanes", "cpl", "mpl", "want"),
                         [(2, 2, 4001, 4),   # ceil(8002 / 4000) = 3, a whole core more: 4
                          (1, 1, 1000, 2),   # an odd request holds a whole core (#445 S1)
                          (3, 1, 2000, 4),   # 3 CPUs requested: 4
                          (2, 2, 4000, 4)])  # exact: no raise  # fmt: skip
def test_allocated_cpus_rounds_memory_up_then_to_whole_cores(
    lanes_mod: ModuleType, lanes: int, cpl: int, mpl: int, want: int
) -> None:
    row = lanes_mod.PlanRow("x", 0, lanes, "*", cpl, mpl, "01:00:00", 0.0)
    assert lanes_mod.allocated_cpus(row) == want


def test_planned_allocated_cpus_with_the_rows_job_fit_r14(lanes_mod: ModuleType) -> None:
    rows = lanes_mod.read_plan(PLAN, _suites())
    requested = sum(r.requested_cpus for r in rows)
    allocated = sum(lanes_mod.allocated_cpus(r) for r in rows)
    assert (requested, allocated) == (122, 126)
    assert allocated + 2 == 128 <= R14_CPUS  # the 1-CPU rows job holds a whole core
    assert len(rows) + 1 <= 20  # jobs


def test_launcher_sums_allocated_cpus(sb: Sandbox) -> None:
    """0 queued + 126 allocated + 2 (rows job, a whole core) = 128 launches; a pending 1-CPU job
    (a whole core: 2) refuses."""
    _queue(sb, ["0 gippa-other"])
    res = _launch(sb)
    assert res.returncode == 0, _out(res)
    assert re.search(r"\bcpus=126\b", res.stdout), res.stdout
    jobs = {option(c["argv"], "job-name", "J"): c["argv"] for c in sb.sbatch_calls()}
    disc = jobs["gippa-cal-discordance"]
    assert option(disc, "cpus-per-task", "c") == "4", disc  # the request itself is unchanged
    assert re.fullmatch(r"32000M?", option(disc, "mem") or ""), disc


def test_launcher_refuses_when_allocated_cpus_pass_r14(sb: Sandbox) -> None:
    _queue(sb, ["1 gippa-other"])
    res = _launch(sb)
    assert res.returncode == 2 and sb.sbatch_calls() == [], _out(res)
    assert re.search(r"\b130\b", _out(res)) and "allocated" in _out(res), _out(res)


def test_launcher_max_mem_per_cpu_override_changes_the_sum(sb: Sandbox) -> None:
    """MaxMemPerCPU 16000: discordance needs no raise (4 CPUs), so 122 suite CPUs."""
    _queue(sb, ["4 gippa-other"])  # 4 + 122 + 2 = 128 fits; at 4000 it would be 132
    res = _launch(sb, GIPPA_MAX_MEM_PER_CPU_MB="16000")
    assert res.returncode == 0, _out(res)
    assert re.search(r"\bcpus=122\b", res.stdout), res.stdout


@pytest.mark.parametrize("bad", ["0", "00", "abc", "-1"])
def test_launcher_refuses_a_bad_max_mem_per_cpu(sb: Sandbox, bad: str) -> None:
    _queue(sb, [])
    res = _launch(sb, GIPPA_MAX_MEM_PER_CPU_MB=bad)
    assert res.returncode == 2 and sb.sbatch_calls() == [], _out(res)
    assert "GIPPA_MAX_MEM_PER_CPU_MB" in _out(res), _out(res)


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


# The merge phase (block 2 of the pending checks, ``--merge`` on the shard lanes after the shards)
# runs inside the same job time. Its time is recorded in ``calibration_durations.tsv`` only through
# the JUnit sum: a row is ``durations_from_junit`` over every part of the run, shard lanes *and*
# ``task-<i>.merge.lane<k>.junit.xml``, so a shard test's row is block 1 + block 2 and
# ``row / shards`` is one lane's shard work plus its merge work. That holds where the merge
# finished: placement and variants (run 45bc3599, jobs 1900336 / 1900349, merge exit 0; e.g.
# variants g1 = 7730.8 s shards + 8604.4 s merge = its 16335.2 s row).
#
# novel_orfs' merge (job 1900347) was killed at 0.95 h: only merge lane 0 of 6 wrote its JUnit, so
# its rows hold 1/6 of the merge. The block-2 records in that run's ``novel_orfs.jsonl`` hold the
# rest (``duration`` of each lane's block-2 slice); the seeds are fixed, so the same checks are
# pending again (R7). Seconds of merge work missing from the rows, per nodeid:
_GC = ("tests/calibration/attacks/test_attack_m7_gate.py::"
       "test_attack_gc_poor_intergenic_fdp_at_level_full")  # fmt: skip
_MULTI = ("tests/calibration/attacks/test_attack_novel_class_fdr.py::"
          "test_denovo_orf_fdr_with_multi_peptide_false_orfs[full]")  # fmt: skip
_G4K3 = ("tests/calibration/test_novel_orfs_gate.py::"
         "test_g4k3_clustered_false_orfs_full[intergenic-0.01]")  # fmt: skip
# Each correction holds only while its row still has the source it was derived from; a row
# re-recorded from a run whose merge finished already includes the merge (``--max``).
MERGE_NOT_IN_DURATIONS: dict[str, dict[str, tuple[float, str]]] = {
    "novel_orfs": {
        # 6 block-2 slices (3328-3534 s) minus lane 0's JUnit time, which the row holds
        _GC: (20487.4 - 3406.041, "1900347"),
        # the row is the probe's (97.0 s), not 1900347's, so it holds none of lane 0's 9.672 s;
        # the probe ran unsharded, so its row may hold block 2 already (59 s counted twice at most)
        _MULTI: (59.2, "probe-1895355-n8"),
    },
}
# novel_orfs ``g4k3.intergenic`` (``fdr_counts``, R = 300) was in the merge's confirm list but the
# kill came before it logged block 2. Block 2 runs the same R replicates as block 1, so it weighs
# its block-1 row (which holds no merge part: lane 0's JUnit has no g4k3 testcase).
MERGE_AS_BLOCK1: dict[str, dict[str, str]] = {"novel_orfs": {_G4K3: "1900347"}}
# Where each sharded suite's merge time comes from. determinants has no recorded merge: it ran
# unsharded in 45bc3599 and its per-model rows (``1901050:measured``) are block 1 only (#445 N6).
# Block 2 runs only for checks whose block 1 is not decided; if all its shard checks went to block
# 2, its lane would be about twice its block-1 estimate.
MERGE_SOURCE = {
    "placement": "rows (1900336 merge finished)",
    "variants": "rows (1900349 merge finished)",
    "novel_orfs": "rows + MERGE_NOT_IN_DURATIONS + MERGE_AS_BLOCK1 (1900347 merge killed)",
    "determinants": "none (block 1 only, 1901050:measured)",
}


def _row_sources() -> dict[str, str]:
    lines = DURATIONS.read_text().splitlines()[1:]
    return {c[0]: c[2] for ln in lines if len(c := ln.split("\t")) == 3}


def _merge_extra_s(suite: str, durations: dict[str, float]) -> float:
    """Merge seconds (summed over lanes) that the durations rows of ``suite`` do not hold. Fails
    once a corrected row was re-recorded: the new row includes the merge, the correction must go."""
    sources = _row_sources()
    named = {n: src for n, (_, src) in MERGE_NOT_IN_DURATIONS.get(suite, {}).items()}
    named.update(MERGE_AS_BLOCK1.get(suite, {}))
    for n, src in named.items():
        if sources.get(n) != src:
            pytest.fail(
                f"{suite}: the durations row of {n} is now from {sources.get(n)!r}, not {src!r}; "
                "it was re-recorded and includes the merge, so delete its correction from "
                "MERGE_NOT_IN_DURATIONS / MERGE_AS_BLOCK1 (test_calibration_plan_fixes.py)"
            )
    extra = sum(sec for sec, _ in MERGE_NOT_IN_DURATIONS.get(suite, {}).values())
    return extra + sum(durations[n] for n in MERGE_AS_BLOCK1.get(suite, {}))


def _estimate_s(mod: ModuleType, row: object, report: dict[str, bool],
                durations: dict[str, float], *, merge: bool = True) -> float:  # fmt: skip
    """The job's estimate: its longest group lane (LPT) or a shard lane's share of the shard tests
    and of the merge (their recorded single-process time, plus the merge work the rows miss, over
    ``shards``); a missing test weighs ``unknown_min``."""
    lanes = mod.lanes(row, report, durations)

    def w(n: str) -> float:
        return float(durations.get(n, row.unknown_min * 60))  # type: ignore[attr-defined]

    loads = [sum(w(n) for n in ln.expected) for ln in lanes if ln.kind == "group"]
    shard = [ln for ln in lanes if ln.kind == "shard"]
    if shard:
        work = sum(w(n) for n in shard[0].expected)
        if merge:
            work += _merge_extra_s(row.suite, durations)  # type: ignore[attr-defined]
        loads.append(work / len(shard))
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


def test_every_sharded_suite_states_its_merge_source(
    lanes_mod: ModuleType, collection: dict[str, dict[str, bool]]
) -> None:
    """A new sharded suite must say where its merge time comes from; the merge corrections name
    shard tests of their suite (a group test never reaches the merge)."""
    rows = lanes_mod.read_plan(PLAN, _suites())
    assert {r.suite for r in rows if r.shards > 0} == set(MERGE_SOURCE)
    durations = lanes_mod.read_durations(DURATIONS)
    for r in rows:
        named = [*MERGE_NOT_IN_DURATIONS.get(r.suite, {}), *MERGE_AS_BLOCK1.get(r.suite, {})]
        if not named:
            continue
        shard = [ln for ln in lanes_mod.lanes(r, collection[r.suite], durations)
                 if ln.kind == "shard"]  # fmt: skip
        for n in named:
            assert n in shard[0].expected, (r.suite, n)
            assert n in durations, n


def test_durations_rows_hold_the_merge_parts(lanes_mod: ModuleType, tmp_path: Path) -> None:
    """Why ``row / shards`` covers a finished merge: the rows sum shard and merge JUnit parts."""
    case = '<testcase classname="tests.calibration.test_v" name="test_g1_full" time="{}"/>'
    parts = {"task-13.shard0": 3865.0, "task-13.shard1": 3866.0,
             "task-13.merge.lane0": 4305.0, "task-13.merge.lane1": 4299.0}  # fmt: skip
    paths = []
    for name, t in parts.items():
        p = tmp_path / f"{name}.junit.xml"
        p.write_text(f"<testsuites><testsuite>{case.format(t)}</testsuite></testsuites>")
        paths.append(p)
    got = lanes_mod.durations_from_junit(paths)
    assert got == {"tests/calibration/test_v.py::test_g1_full": pytest.approx(sum(parts.values()))}


@pytest.mark.parametrize("nodeid", [_GC, _G4K3])
def test_merge_correction_expires_when_its_row_is_re_recorded(
    lanes_mod: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nodeid: str
) -> None:
    """A run whose merge finishes raises the row (``--max``); keeping the correction would count
    the merge twice, so the estimate refuses and says to delete it."""
    here = sys.modules[__name__]
    assert here._merge_extra_s("novel_orfs", lanes_mod.read_durations(DURATIONS)) > 0
    lines = DURATIONS.read_text().splitlines()
    out = [f"{ln.split(chr(9))[0]}\t99999.0\t1909999" if ln.startswith(nodeid + "\t") else ln
           for ln in lines]  # fmt: skip
    assert out != lines
    d = tmp_path / "calibration_durations.tsv"
    d.write_text("\n".join(out) + "\n")
    monkeypatch.setattr(here, "DURATIONS", d)
    with pytest.raises(pytest.fail.Exception, match=r"re-recorded .* delete its correction"):
        here._merge_extra_s("novel_orfs", lanes_mod.read_durations(d))


def test_novel_orfs_estimate_counts_its_killed_merge(
    lanes_mod: ModuleType, collection: dict[str, dict[str, bool]]
) -> None:
    """Dry run 1900347 killed the merge after 0.95 h; the rows hold lane 0's part only. The
    estimate adds the other lanes' block-2 gc time and g4k3's pending block 2 over 8 shards."""
    row = {r.suite: r for r in lanes_mod.read_plan(PLAN, _suites())}["novel_orfs"]
    d = lanes_mod.read_durations(DURATIONS)
    lanes_only = _estimate_s(lanes_mod, row, collection["novel_orfs"], d, merge=False)
    est = _estimate_s(lanes_mod, row, collection["novel_orfs"], d)
    assert est - lanes_only == pytest.approx((17081.359 + 59.2 + d[_G4K3]) / row.shards)
    assert est > 4 * 3600  # 1900347: six shard lanes of 4.7 h, then the merge


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
    failed = '<testsuite><testcase time="0.1"><failure message="x"/></testcase></testsuite>'
    assert rows_mod._testcases([ET.fromstring(failed)]) == {("", ""): "failed"}  # still counts


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
    junit.write_text('<testsuite><testcase time="0.1"><error message="x"/></testcase></testsuite>')
    assert job._junit_count(junit) == 1  # a nameless testcase with an error still counts
