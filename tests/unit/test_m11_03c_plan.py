"""``scripts/ci/calibration_lanes.py`` part 1: the plan (issue M11-03c, PR 1a).

Spec: docs/specs/M11-03c-calibration-full-sharding.md, "The plan", "Plan validation", "Exit
folding", "Record"; example 1 (the committed plan), example 4 (``-`` with shard-0-only tests),
example 5 (``fold_exit``); acceptance criteria 1 and the ``fold_exit`` part of 3.

Conventions settled by the test author where the spec leaves them open:

- ``PlanRow`` exposes the plan columns as attributes of the same names; ``shards``, ``groups``,
  ``cpus_per_lane`` and ``mem_per_lane_mb`` are ``int``. ``time`` may be kept as the ``HH:MM:SS``
  string, as seconds or as a ``timedelta`` (``_seconds`` below accepts all three).
- ``Lane`` is a dataclass built with keywords ``kind`` (``"shard"`` / ``"group"``), ``index``,
  ``n``, ``deselect``, ``part`` and ``expected``. An empty group lane has ``expected == ()`` and
  ``part is None``. ``record_parts`` lists shard parts, then non-empty group parts, in the order
  given, then ``merge_parts``.
- CLI: ``calibration_lanes.py check [--plan PATH]``; the default plan is
  ``scripts/slurm/calibration_plan.tsv`` next to the script's repository and the suites are the
  ``replicates_full`` keys of ``tests/calibration/tolerances.yaml``. Exit 0 = valid, 2 = PlanError
  (message on stderr).

Dependency: the ``-`` check needs the collection, so it is tested through ``lanes(row, report,
durations)`` (spec name; the full LPT assignment is PR 1b). PR 1a must at least raise
``PlanError`` there; the positive ``lanes`` cases live in ``tests/unit/test_m11_03c_lanes.py`` (1b).

No statistical claim is made here (scheduling only, R3), so no calibration test applies.
"""

from __future__ import annotations

import datetime
import importlib.util
import itertools
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import yaml
from hypothesis import given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "ci" / "calibration_lanes.py"
PLAN = REPO / "scripts" / "slurm" / "calibration_plan.tsv"
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"
HEADER = ("suite", "shards", "groups", "group_select", "cpus_per_lane", "mem_per_lane_mb", "time",
          "unknown_min")  # fmt: skip

# The spec's "Initial plan" table as re-sized for wheat after dry run 45bc3599 (jobs
# 1900336-1900351, local/handoff/cal_timeouts.md): placement, regions, autocorrelation and
# determinants 08:00:00; gene_fdr 3 groups; autocorrelation 3 x 2 CPUs (was 1 x 4); determinants
# 12 shards + 3 direct groups (was 0 + 5). Columns: suite, shards, groups, group_select,
# cpus_per_lane, mem, time. Full run b3d84613 (2026-10-29): the 8 h guard interrupted placement
# group lanes (7.0-7.8 h) and novel_orfs shards (6.9-7.8 h); Khoa's decision "More lanes, ~140
# CPUs": placement 8 -> 10 direct groups, novel_orfs 8 -> 12 shards.
INITIAL = [
    ("placement", 2, 10, "direct", 2, 4000, "08:00:00"),
    ("gene_fdr", 0, 3, "*", 2, 4000, "06:00:00"),
    ("detectability", 0, 1, "*", 2, 4000, "02:00:00"),
    ("regions", 0, 5, "*", 2, 4000, "08:00:00"),
    ("autocorrelation", 0, 3, "*", 2, 4000, "08:00:00"),
    ("islands", 0, 1, "*", 2, 4000, "06:00:00"),
    ("discordance", 0, 2, "*", 2, 16000, "06:00:00"),
    ("missed_operons", 0, 1, "*", 2, 8000, "02:00:00"),
    ("ora", 0, 2, "*", 4, 4000, "06:00:00"),
    ("gsea", 0, 2, "*", 2, 4000, "06:00:00"),
    ("determinants", 12, 3, "direct", 2, 4000, "08:00:00"),
    ("novel_orfs", 12, 1, "direct", 2, 4000, "08:00:00"),  # b3d84613: 8 shards 6.9-7.8 h
    ("denovo_precision", 0, 1, "*", 2, 4000, "02:00:00"),
    ("variants", 2, 1, "direct", 2, 4000, "06:00:00"),
    ("compare", 0, 1, "*", 2, 4000, "02:00:00"),
]
SUITES = [r[0] for r in INITIAL]


@pytest.fixture(scope="module")
def mod() -> ModuleType:
    assert SCRIPT.is_file(), f"{SCRIPT} does not exist (M11-03c PR 1a)"
    spec = importlib.util.spec_from_file_location("gippa_calibration_lanes", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _seconds(t: Any) -> int:
    if isinstance(t, datetime.timedelta):
        return int(t.total_seconds())
    if isinstance(t, int | float):
        return int(t)
    h, m, s = (int(x) for x in str(t).split(":"))
    return h * 3600 + m * 60 + s


def _rows(**override: dict[str, str]) -> list[dict[str, str]]:
    """The initial plan as string cells; ``override[suite]`` replaces cells of that row."""
    rows = []
    for suite, shards, groups, sel, cpus, mem, time in INITIAL:
        row = {"suite": suite, "shards": str(shards), "groups": str(groups), "group_select": sel,
               "cpus_per_lane": str(cpus), "mem_per_lane_mb": str(mem), "time": time,
               "unknown_min": "30"}  # fmt: skip
        row.update(override.get(suite, {}))
        rows.append(row)
    return rows


def _write(path: Path, rows: list[dict[str, str]], header: tuple[str, ...] = HEADER) -> Path:
    lines = ["\t".join(header)] + ["\t".join(r.get(c, "") for c in header) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


# ------------------------------------------------------------------ the committed plan (example 1)


def test_suites_are_the_replicates_full_keys() -> None:
    with TOLERANCES.open() as fh:
        assert [str(k) for k in yaml.safe_load(fh)["replicates_full"]] == SUITES


def test_committed_plan_is_the_initial_plan(mod: ModuleType) -> None:
    assert PLAN.is_file(), f"{PLAN} missing"
    assert PLAN.read_text().splitlines()[0].split("\t") == list(HEADER)
    rows = mod.read_plan(PLAN, SUITES)
    assert len(rows) == 15
    got = [(r.suite, r.shards, r.groups, r.group_select, r.cpus_per_lane, r.mem_per_lane_mb,
            _seconds(r.time)) for r in rows]  # fmt: skip
    assert got == [(*r[:6], _seconds(r[6])) for r in INITIAL]
    assert all(float(r.unknown_min) > 0 for r in rows)


def test_example_1_totals(mod: ModuleType) -> None:
    rows = mod.read_plan(PLAN, SUITES)
    assert sum(r.shards + r.groups for r in rows) == 65
    # requested (lanes x cpus_per_lane) 134; allocated 138: discordance's 2 x 16000 MB gets
    # ceil(32000 / MaxMemPerCPU 4000) = 8 CPUs, not 4
    assert sum((r.shards + r.groups) * r.cpus_per_lane for r in rows) == 134
    assert sum(mod.allocated_cpus(r) for r in rows) == 138
    assert max(_seconds(r.time) for r in rows) == 8 * 3600 == mod.MAX_TIME_S
    assert [r.suite for r in rows if r.cpus_per_lane != 2] == ["ora"]  # Q6; autocorrelation 3 x 2


def test_plan_lives_outside_calib_paths() -> None:
    rel = PLAN.relative_to(REPO).as_posix()
    assert not rel.startswith(("src/", "tests/calibration/", "ci/", "tests/sim/", "tests/helpers/"))


# ------------------------------------------------------------------ plan validation (criterion 1)

S = "placement"
BAD = {  # rule -> (overrides, columns of which the message names one)
    "shards_negative": ({S: {"shards": "-1"}}, ("shards",)),
    "groups_negative": ({S: {"groups": "-2"}}, ("groups",)),
    "shards_not_integer": ({S: {"shards": "1.5"}}, ("shards",)),
    "shards_word": ({S: {"shards": "two"}}, ("shards",)),
    "groups_not_integer": ({S: {"groups": "x"}}, ("groups",)),
    "groups_empty": ({S: {"groups": ""}}, ("groups",)),
    "no_lanes": ({"compare": {"shards": "0", "groups": "0"}}, ("shards", "groups")),
    "star_with_shards": ({"gene_fdr": {"shards": "2"}}, ("group_select", "shards")),
    "dash_with_groups": ({"compare": {"group_select": "-"}}, ("group_select", "groups")),
    "direct_without_shards": ({"compare": {"group_select": "direct"}}, ("group_select", "shards")),
    "direct_list_without_shards": (
        {"compare": {"group_select": "direct,tests/calibration/test_compare.py::"}},
        ("group_select", "shards"),
    ),
    "cpus_odd": ({S: {"cpus_per_lane": "3"}}, ("cpus_per_lane",)),
    "cpus_one": ({"ora": {"cpus_per_lane": "1"}}, ("cpus_per_lane",)),
    "time_words": ({S: {"time": "6h"}}, ("time",)),
    "time_two_fields": ({S: {"time": "06:00"}}, ("time",)),
    "time_bad_minutes": ({S: {"time": "05:60:00"}}, ("time",)),
    "time_empty": ({S: {"time": ""}}, ("time",)),
    "time_over_8h": ({S: {"time": "08:00:01"}}, ("time",)),
    "time_9h": ({S: {"time": "09:00:00"}}, ("time",)),
    "prefix_without_groups": (
        {S: {"groups": "0", "group_select": "tests/calibration/test_placement_fdr_sim.py::"}},
        ("group_select", "groups"),
    ),
    "direct_without_groups": ({S: {"groups": "0"}}, ("group_select", "groups")),
    "cpus_zero": ({S: {"cpus_per_lane": "0"}}, ("cpus_per_lane",)),
    "mem_zero": ({S: {"mem_per_lane_mb": "0"}}, ("mem_per_lane_mb",)),
    "time_zero": ({S: {"time": "00:00:00"}}, ("time",)),
    "unknown_zero": ({S: {"unknown_min": "0"}}, ("unknown_min",)),
    "unknown_negative": ({S: {"unknown_min": "-5"}}, ("unknown_min",)),
    "unknown_nan": ({S: {"unknown_min": "nan"}}, ("unknown_min",)),
    "unknown_inf": ({S: {"unknown_min": "inf"}}, ("unknown_min",)),
    "unknown_word": ({S: {"unknown_min": "x"}}, ("unknown_min",)),
    "select_trailing_comma": ({S: {"group_select": "direct,"}}, ("group_select",)),
    "select_leading_comma": (
        {S: {"group_select": ",tests/calibration/test_placement_fdr_sim.py::"}},
        ("group_select",),
    ),
}


@pytest.mark.parametrize("rule", sorted(BAD))
def test_read_plan_rejects_rule_naming_row_and_column(mod: ModuleType, tmp_path: Path,
                                                      rule: str) -> None:  # fmt: skip
    override, columns = BAD[rule]
    (suite,) = override
    path = _write(tmp_path / "plan.tsv", _rows(**override))
    with pytest.raises(mod.PlanError) as exc:
        mod.read_plan(path, SUITES)
    msg = str(exc.value)
    row = SUITES.index(suite)
    assert suite in msg or f"row {row}" in msg or f"row {row + 1}" in msg, msg
    assert any(c in msg for c in columns), msg


SUITE_ORDER = {
    "swapped": lambda r: [r[1], r[0], *r[2:]],
    "missing_last": lambda r: r[:-1],
    "extra_row": lambda r: [*r, {**r[0], "suite": "bogus"}],
    "duplicate": lambda r: [*r[:-1], r[0]],
    "renamed": lambda r: [*r[:-1], {**r[-1], "suite": "Compare"}],
}


@pytest.mark.parametrize("case", sorted(SUITE_ORDER))
def test_read_plan_rejects_suites_not_the_keys_in_order(mod: ModuleType, tmp_path: Path,
                                                        case: str) -> None:  # fmt: skip
    path = _write(tmp_path / "plan.tsv", SUITE_ORDER[case](_rows()))
    with pytest.raises(mod.PlanError) as exc:
        mod.read_plan(path, SUITES)
    assert "suite" in str(exc.value)


def test_read_plan_rejects_empty_plan_and_bad_header(mod: ModuleType, tmp_path: Path) -> None:
    for text in ("", "\t".join(HEADER) + "\n"):
        (tmp_path / "p.tsv").write_text(text)
        with pytest.raises(mod.PlanError):
            mod.read_plan(tmp_path / "p.tsv", SUITES)
    with pytest.raises(mod.PlanError, match="unknown_min"):
        mod.read_plan(_write(tmp_path / "q.tsv", _rows(), HEADER[:-1]), SUITES)
    short = tmp_path / "r.tsv"
    _write(short, _rows())
    lines = short.read_text().splitlines()
    lines[3] = lines[3].rsplit("\t", 1)[0]  # one cell missing in the detectability row
    short.write_text("\n".join(lines) + "\n")
    with pytest.raises(mod.PlanError):
        mod.read_plan(short, SUITES)


@pytest.mark.parametrize("override", [
    {S: {"time": "08:00:00"}},  # the bound itself
    {S: {"group_select": "direct,tests/calibration/test_placement_fdr_sim.py::"}},
    {S: {"group_select": "tests/calibration/test_placement_fdr_sim.py::,tests/calibration/x.py"}},
    {"compare": {"shards": "3", "groups": "0", "group_select": "-"}},
    {"compare": {"cpus_per_lane": "4"}},
])  # fmt: skip
def test_read_plan_accepts_valid_variants(mod: ModuleType, tmp_path: Path,
                                          override: dict[str, dict[str, str]]) -> None:  # fmt: skip
    rows = mod.read_plan(_write(tmp_path / "plan.tsv", _rows(**override)), SUITES)
    assert [r.suite for r in rows] == SUITES


ROW = st.fixed_dictionaries({
    "shards": st.integers(0, 12), "groups": st.integers(0, 12),
    "cpus": st.sampled_from([2, 4, 6, 8]), "mem": st.integers(500, 64000),
    "secs": st.integers(60, 8 * 3600), "unknown": st.integers(1, 240),
}).filter(lambda d: d["shards"] + d["groups"] > 0)  # fmt: skip


@settings(max_examples=60, deadline=None)
@given(st.lists(ROW, min_size=15, max_size=15), st.data())
def test_valid_random_plans_round_trip(
    mod: ModuleType,
    tmp_path_factory: pytest.TempPathFactory,
    spec: list[dict[str, int]],
    data: st.DataObject,
) -> None:
    override = {}
    for suite, d in zip(SUITES, spec, strict=True):
        if d["shards"] == 0:
            sel = "*"
        elif d["groups"] == 0:
            sel = "-"
        else:
            sel = data.draw(st.sampled_from(["direct", "direct,tests/calibration/test_a.py::t_x"]))
        h, rem = divmod(d["secs"], 3600)
        override[suite] = {"shards": str(d["shards"]), "groups": str(d["groups"]),
                           "group_select": sel, "cpus_per_lane": str(d["cpus"]),
                           "mem_per_lane_mb": str(d["mem"]),
                           "time": f"{h:02d}:{rem // 60:02d}:{rem % 60:02d}",
                           "unknown_min": str(d["unknown"])}  # fmt: skip
    path = _write(tmp_path_factory.mktemp("plan") / "plan.tsv", _rows(**override))
    rows = mod.read_plan(path, SUITES)
    for r, d in zip(rows, spec, strict=True):
        assert (r.shards, r.groups, r.cpus_per_lane, r.mem_per_lane_mb) == (
            d["shards"], d["groups"], d["cpus"], d["mem"])  # fmt: skip
        assert _seconds(r.time) == d["secs"]
        assert float(r.unknown_min) == d["unknown"]


# ------------------------------------------------------------------ `-` with shard-0-only (ex. 4)


def test_dash_with_shard0_only_tests_is_a_plan_error(mod: ModuleType, tmp_path: Path) -> None:
    path = _write(
        tmp_path / "plan.tsv", _rows(novel_orfs={"shards": "6", "groups": "0", "group_select": "-"})
    )
    (row,) = [r for r in mod.read_plan(path, SUITES) if r.suite == "novel_orfs"]
    f = "tests/calibration/test_novel_orfs_full.py::"
    report = {f + "test_planted_tir_calibration_300k_full[a]": True,
              f + "test_planted_tir_calibration_300k_full[b]": True,
              f + "test_sixframe_decoys_null_and_coin_full": True,
              f + "test_attack_gc_poor_intergenic_fdp_at_level_full": False,
              f + "test_harness_check_full": False}  # fmt: skip
    with pytest.raises(mod.PlanError) as exc:
        mod.lanes(row, report, {})
    msg = str(exc.value)
    assert "novel_orfs" in msg and "test_planted_tir_calibration_300k_full" in msg, msg


# ------------------------------------------------------------------ fold_exit (example 5)


@pytest.mark.parametrize(("lanes", "merge", "want"), [
    ([0, 0, 1], 0, 1), ([0, 143, 1], 0, 143), ([0, 0], 2, 2), ([0, 0], None, 0),
])  # fmt: skip
def test_fold_exit_example_5(mod: ModuleType, lanes: list[int], merge: int | None,
                             want: int) -> None:  # fmt: skip
    assert mod.fold_exit(lanes, merge) == want


def _fold_rule(lanes: tuple[int, ...], merge: int | None) -> int:
    seq = [*lanes, *([] if merge is None else [merge])]
    crash = [e for e in seq if e not in (0, 1)]
    return crash[0] if crash else int(1 in seq)


def test_fold_exit_rule_over_all_combinations(mod: ModuleType) -> None:
    codes = (0, 1, 2, 143)
    for n in range(1, 5):
        for lanes in itertools.product(codes, repeat=n):
            for merge in (None, *codes):
                assert mod.fold_exit(list(lanes), merge) == _fold_rule(lanes, merge), (lanes, merge)


def test_fold_exit_first_crash_wins_in_order(mod: ModuleType) -> None:
    assert mod.fold_exit([143, 2], 137) == 143
    assert mod.fold_exit([1, 137, 2], 143) == 137
    assert mod.fold_exit([0, 1], 2) == 2  # the merge's crash beats a lane's 1
    assert mod.fold_exit([1, 0], 1) == 1


# ------------------------------------------------------------------ record_parts (Record)


def _lane(mod: ModuleType, kind: str, index: int, n: int, part: str | None,
          expected: tuple[str, ...]) -> Any:  # fmt: skip
    return mod.Lane(kind=kind, index=index, n=n, deselect=(), part=part, expected=expected)


def test_record_parts_lists_parts_launched(mod: ModuleType) -> None:
    lanes = [_lane(mod, "shard", 0, 2, "task-0.shard0.junit.xml", ("a::t",)),
             _lane(mod, "shard", 1, 2, "task-0.shard1.junit.xml", ("a::t",)),
             _lane(mod, "group", 0, 3, "task-0.group0.junit.xml", ("a::u",)),
             _lane(mod, "group", 1, 3, None, ()),  # empty group lane: dropped, no part
             _lane(mod, "group", 2, 3, "task-0.group2.junit.xml", ("a::v", "a::w"))]  # fmt: skip
    merge = ["task-0.merge.lane0.junit.xml", "task-0.merge.lane1.junit.xml"]
    assert mod.record_parts(lanes, merge) == [
        "task-0.shard0.junit.xml", "task-0.shard1.junit.xml", "task-0.group0.junit.xml",
        "task-0.group2.junit.xml", *merge]  # fmt: skip


def test_record_parts_example_2_unsharded(mod: ModuleType) -> None:
    lanes = [_lane(mod, "group", g, 5, f"task-10.group{g}.junit.xml", (f"d::t{g}",))
             for g in range(5)]  # fmt: skip
    assert mod.record_parts(lanes, []) == [f"task-10.group{g}.junit.xml" for g in range(5)]


# ------------------------------------------------------------------ check CLI


def _check(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPT), "check", *args], capture_output=True,
                          text=True, timeout=60, check=False)  # fmt: skip


def test_check_cli_accepts_committed_plan() -> None:
    assert SCRIPT.is_file()
    res = _check()
    assert res.returncode == 0, res.stderr
    res = _check("--plan", str(PLAN))
    assert res.returncode == 0, res.stderr


@pytest.mark.parametrize("override", [
    {S: {"cpus_per_lane": "3"}}, {S: {"time": "09:00:00"}}, {"gene_fdr": {"shards": "1"}},
])  # fmt: skip
def test_check_cli_exits_2_on_plan_error(tmp_path: Path,
                                         override: dict[str, dict[str, str]]) -> None:  # fmt: skip
    assert SCRIPT.is_file()
    res = _check("--plan", str(_write(tmp_path / "plan.tsv", _rows(**override))))
    assert res.returncode == 2
    (suite,) = override
    assert suite in res.stderr


def test_check_cli_exits_2_when_suites_differ_from_tolerances(tmp_path: Path) -> None:
    assert SCRIPT.is_file()
    res = _check("--plan", str(_write(tmp_path / "plan.tsv", _rows()[:-1])))
    assert res.returncode == 2 and "suite" in res.stderr
