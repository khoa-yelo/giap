"""``scripts/ci/calibration_lanes.py`` part 2: lane assignment (issue M11-03c, PR 1b).

Spec: docs/specs/M11-03c-calibration-full-sharding.md, "The plan" (group selection), "Group
assignment", "``--deselect`` is a prefix match", "Plan validation" (the checks that need the
collection), "In-job partition check", "Edge cases"; examples 2 and 3; acceptance criterion 2 and
the ``verify_partition`` / ``verify_counts`` part of criterion 3.

Conventions settled by the test author where the spec leaves them open:

- ``lanes(row, report, durations)``: ``report`` is ``{nodeid: shard0_only}`` (the values of the
  conftest lane report, PR 2); ``durations`` is ``{nodeid: seconds}``; ``unknown_min`` is minutes,
  so a missing test weighs ``unknown_min * 60`` seconds. The result lists the ``shards`` shard
  lanes (index 0..N-1), then all ``groups`` group lanes (index 0..G-1), empty ones included with
  ``expected == ()`` and ``part is None`` (the record names them with ``tests: 0``). A shard lane
  has ``n == shards``, a group lane ``n == groups``. ``part`` is ``task-<i>.shard<k>.junit.xml`` /
  ``task-<i>.group<g>.junit.xml``, ``i`` the row's index among the ``replicates_full`` keys. Rows
  in these tests come from ``read_plan`` of the committed plan (``dataclasses.replace`` for
  variants), so ``i`` may come from the suite name or from a ``PlanRow`` field.
- Deselect lists are exact: a shard lane deselects exactly the group tests, group lane ``g``
  exactly the collected tests not in ``g``. ``deselect`` and ``expected`` are compared as sets in
  the examples; the order-independence property compares whole ``Lane`` objects.
- A ``group_select`` prefix list with ``shards == 0`` (review point N2 of PR 1a): the spec makes a
  test a group test iff it matches a prefix, and with no shard lane an unmatched test would run
  nowhere ("No test is ever dropped"). Such a row is a ``PlanError``, either in ``read_plan`` or,
  when a collected test matches no prefix, in ``lanes``.
- An empty collection with ``shards == 0`` may raise ``PlanError`` or give only empty lanes (the
  job treats an empty collection as a collection error either way); with ``shards > 0`` it is the
  "no test left for the shard lanes" ``PlanError``.
- ``verify_partition(row, lanes, collected)``: ``collected`` maps each non-empty lane's ``part`` to
  the nodeids its ``--collect-only`` printed and the key ``"all"`` to the suite's full collection.
  A missing key is a problem, not an exception.
- ``verify_counts(lanes, junit_counts)``: ``junit_counts`` maps ``part`` to the JUnit testcase
  count; a non-empty lane whose part is missing is a problem (fail closed).
- ``read_durations(path)`` reads ``nodeid seconds source`` (tab-separated, header); a negative or
  non-numeric ``seconds`` raises ``ValueError`` (``PlanError`` qualifies).
- ``durations_from_junit(paths)`` maps ``classname`` + ``name`` of each pytest testcase back to
  ``<path>.py::<name>`` with its ``time`` (calibration tests are module-level functions).
- CLI ``lanes --suite S --report PATH [--durations PATH] [--plan PATH]`` prints the lanes as a
  JSON list of objects with the ``Lane`` field names (exit 2 + stderr on ``PlanError``);
  ``durations --source ID --out PATH JUNIT...`` writes the durations TSV.

No statistical claim is made here (scheduling only, R3), so no calibration test applies.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "ci" / "calibration_lanes.py"
PLAN = REPO / "scripts" / "slurm" / "calibration_plan.tsv"
DURATIONS = REPO / "scripts" / "slurm" / "calibration_durations.tsv"
SUITES = ["placement", "gene_fdr", "detectability", "regions", "autocorrelation", "islands",
          "discordance", "missed_operons", "ora", "gsea", "determinants", "novel_orfs",
          "denovo_precision", "variants", "compare"]  # fmt: skip
DET = "tests/calibration/attacks/test_attack_determinants_null.py::test_null_level_under_attack"
ATTACK5 = ("tests/calibration/attacks/test_attack5_reversed_stretches.py::"
           "test_attack5_reversed_stretches_shuffle_decoys_full")  # fmt: skip


@pytest.fixture(scope="module")
def mod() -> ModuleType:
    assert SCRIPT.is_file(), f"{SCRIPT} does not exist (M11-03c PR 1a)"
    spec = importlib.util.spec_from_file_location("giap_calibration_lanes_1b", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def plan(mod: ModuleType) -> dict[str, Any]:
    return {r.suite: r for r in mod.read_plan(PLAN, SUITES)}


def _row(plan: Mapping[str, Any], suite: str, **change: Any) -> Any:
    return dataclasses.replace(plan[suite], **change) if change else plan[suite]


def _part(suite: str, kind: str, k: int) -> str:
    return f"task-{SUITES.index(suite)}.{kind}{k}.junit.xml"


def _split(out: list[Any]) -> tuple[list[Any], list[Any]]:
    shard = [ln for ln in out if ln.kind == "shard"]
    group = [ln for ln in out if ln.kind == "group"]
    assert len(shard) + len(group) == len(out), [ln.kind for ln in out]
    return shard, group


def _kept(collection: list[str], deselect: tuple[str, ...]) -> set[str]:
    """What pytest runs: ``--deselect`` drops every nodeid that starts with an entry."""
    prefixes = tuple(deselect)
    return {n for n in collection if not (prefixes and n.startswith(prefixes))}


# ------------------------------------------------------------------ reference model (spec text)


class _Err(Exception):
    pass


def _model(row: Any, report: Mapping[str, bool], durations: Mapping[str, float]
           ) -> tuple[set[str], list[list[str]]]:  # fmt: skip
    """(shard set, group lanes in assignment order); ``_Err`` where the spec asks for PlanError."""
    ids = sorted(report)
    sel = row.group_select
    if sel == "*":
        group = set(ids)
    elif sel == "-":
        if any(report.values()):
            raise _Err("dash with shard-0-only")
        group = set()
    else:
        prefixes = tuple(t.strip() for t in sel.split(",") if t.strip() != "direct")
        group = {n for n in ids if (row.shards > 0 and report[n])
                 or (prefixes and n.startswith(prefixes))}  # fmt: skip
        if row.shards == 0 and set(ids) - group:
            raise _Err("prefix list with shards == 0 drops a test")
    shard_set = set(ids) - group
    if row.shards > 0 and not shard_set:
        raise _Err("no test left for the shard lanes")

    def weight(n: str) -> float:
        return durations[n] if n in durations else row.unknown_min * 60

    loads = [0.0] * row.groups
    lanes_: list[list[str]] = [[] for _ in range(row.groups)]
    for n in sorted(group, key=lambda n: (-weight(n), n)):
        g = min(range(row.groups), key=lambda j: (loads[j], j))
        loads[g] += weight(n)
        lanes_[g].append(n)
    checks = [(group, shard_set)] if row.shards > 0 else []
    checks += [(set(ids) - set(ln), set(ln)) for ln in lanes_ if ln]
    for deselect, keep in checks:
        for k in keep:
            if any(k != d and k.startswith(d) for d in deselect):
                raise _Err(f"prefix collision at {k}")
    return shard_set, lanes_


# ------------------------------------------------------------------ hand examples


def _lpt_example(plan: Mapping[str, Any]) -> tuple[Any, dict[str, bool], dict[str, float]]:
    row = _row(plan, "regions", groups=3, unknown_min=5.0)
    report = {f"tests/calibration/test_r.py::t_{c}": False for c in "gfedcba"}
    minutes = {"a": 10, "b": 9, "c": 8, "d": 7, "e": 6}  # f, g: unknown_min = 5 min each
    durations = {f"tests/calibration/test_r.py::t_{c}": m * 60.0 for c, m in minutes.items()}
    return row, report, durations


def test_lpt_by_hand_with_ties_and_unknown_durations(mod: ModuleType, plan: dict[str, Any]) -> None:
    # sorted (10, 9, 8, 7, 6, 5 f, 5 g): a->0, b->1, c->2, d->2 (8), e->1 (9), f->0 (10),
    # g: loads 15/15/15 tie -> lowest lane 0.
    row, report, durations = _lpt_example(plan)
    out = mod.lanes(row, report, durations)
    shard, group = _split(out)
    assert shard == []
    assert [ln.index for ln in group] == [0, 1, 2] and {ln.n for ln in group} == {3}
    f = "tests/calibration/test_r.py::t_"
    assert [set(ln.expected) for ln in group] == [
        {f + "a", f + "f", f + "g"},
        {f + "b", f + "e"},
        {f + "c", f + "d"},
    ]
    assert [ln.part for ln in group] == [_part("regions", "group", g) for g in range(3)]
    for ln in group:
        assert set(ln.deselect) == set(report) - set(ln.expected)


def test_empty_durations_means_unknown_min_for_all(mod: ModuleType, plan: dict[str, Any]) -> None:
    row = _row(plan, "gene_fdr")  # 0 shards, 2 groups, '*'
    report = {f"tests/calibration/test_g.py::t_{c}": False for c in "dcba"}
    out = mod.lanes(row, report, {})
    _, group = _split(out)
    f = "tests/calibration/test_g.py::t_"
    assert [set(ln.expected) for ln in group] == [{f + "a", f + "c"}, {f + "b", f + "d"}]


def _example2(plan: Mapping[str, Any]) -> tuple[Any, dict[str, bool], dict[str, float]]:
    params = [
        f"{p}-{m}"
        for p in (
            "autocorr_pos",
            "detect",
            "size",
            "clustered",
            "gc",
            "len",
            "codon",
            "phase",
            "tss",
        )
        for m in ("lightgbm", "elastic_net")
    ]
    power = "tests/calibration/test_determinants_null_calibration.py::test_power_planted[lightgbm]"
    ids = [*[f"{DET}[full-{p}]" for p in params][:18], power]
    assert len(ids) == 19
    durations = {n: 600.0 + 37.0 * j for j, n in enumerate(sorted(ids))}
    durations[f"{DET}[full-autocorr_pos-lightgbm]"] = 8364.0
    durations[f"{DET}[full-detect-lightgbm]"] = 7782.0
    return _row(plan, "determinants"), dict.fromkeys(ids, True), durations


def test_example_2_determinants(mod: ModuleType, plan: dict[str, Any]) -> None:
    row, report, durations = _example2(plan)
    assert (row.shards, row.groups, row.group_select) == (0, 5, "*")
    out = mod.lanes(row, report, durations)
    shard, group = _split(out)
    assert shard == [] and len(group) == 5
    assert [ln.index for ln in group] == list(range(5)) and {ln.n for ln in group} == {5}
    sets = [set(ln.expected) for ln in group]
    assert sorted(n for s in sets for n in s) == sorted(report)  # each nodeid in exactly one lane
    assert f"{DET}[full-autocorr_pos-lightgbm]" in sets[0]
    assert f"{DET}[full-detect-lightgbm]" in sets[1]
    for g, ln in enumerate(group):
        others = set().union(*(s for j, s in enumerate(sets) if j != g))
        assert set(ln.deselect) == others
        assert ln.part == f"task-10.group{g}.junit.xml"
    assert mod.record_parts(out, []) == [f"task-10.group{g}.junit.xml" for g in range(5)]
    _, want = _model(row, report, durations)
    assert sets == [set(w) for w in want]


def _example3(plan: Mapping[str, Any]) -> tuple[Any, dict[str, bool], dict[str, float]]:
    direct = [f"tests/calibration/test_placement_fdr_sim.py::test_direct_{j:02d}"
              for j in range(28)]  # fmt: skip
    harness = [f"tests/calibration/test_placement_harness.py::test_check_{j:02d}_full"
               for j in range(29)] + [ATTACK5]  # fmt: skip
    report = {**dict.fromkeys(direct, True), **dict.fromkeys(harness, False)}
    durations = {n: 4140.0 + 60.0 * (j % 11) for j, n in enumerate(direct)}
    durations[ATTACK5] = 12096.0
    return _row(plan, "placement"), report, durations


def test_example_3_placement(mod: ModuleType, plan: dict[str, Any]) -> None:
    row, report, durations = _example3(plan)
    assert (row.shards, row.groups, row.group_select) == (2, 8, "direct")
    direct = {n for n, s in report.items() if s}
    harness = set(report) - direct
    assert (len(direct), len(harness)) == (28, 30)
    out = mod.lanes(row, report, durations)
    shard, group = _split(out)
    assert [(ln.index, ln.n) for ln in shard] == [(0, 2), (1, 2)]
    for k, ln in enumerate(shard):
        assert set(ln.deselect) == direct
        assert set(ln.expected) == harness and ATTACK5 in ln.expected
        assert ln.part == f"task-0.shard{k}.junit.xml"
    assert shard[0].deselect == shard[1].deselect and shard[0].expected == shard[1].expected
    assert [(ln.index, ln.n) for ln in group] == [(g, 8) for g in range(8)]
    sets = [set(ln.expected) for ln in group]
    assert sorted(n for s in sets for n in s) == sorted(direct)
    assert all(sets) and [ln.part for ln in group] == [f"task-0.group{g}.junit.xml"
                                                        for g in range(8)]  # fmt: skip
    _, want = _model(row, report, durations)
    assert sets == [set(w) for w in want]
    for ln in out:
        assert _kept(sorted(report), ln.deselect) == set(ln.expected)


def test_direct_makes_new_tests_land_by_flag(mod: ModuleType, plan: dict[str, Any]) -> None:
    # Edge case "new test added after the plan": shard-0-only -> group lane, else shard lanes.
    row = _row(plan, "variants")  # 2 shards, 1 group, 'direct'
    report = {"tests/calibration/test_v.py::t_new_direct": True,
              "tests/calibration/test_v.py::t_new_harness": False}  # fmt: skip
    shard, group = _split(mod.lanes(row, report, {}))
    assert [ln.expected for ln in shard] == [("tests/calibration/test_v.py::t_new_harness",)] * 2
    assert [ln.expected for ln in group] == [("tests/calibration/test_v.py::t_new_direct",)]
    assert [ln.part for ln in shard + group if ln.part is not None] == [
        "task-13.shard0.junit.xml",
        "task-13.shard1.junit.xml",
        "task-13.group0.junit.xml",
    ]


def test_hand_prefix_adds_harness_test_and_shard0_only_still_grouped(
    mod: ModuleType, plan: dict[str, Any]
) -> None:
    f = "tests/calibration/test_p.py::"
    report = {f + "t_heavy_full": False, f + "t_light_full": False, f + "t_direct": True}
    for sel in (f + "t_heavy", f"direct,{f}t_heavy"):  # shards > 0: shard-0-only is grouped anyway
        row = _row(plan, "placement", groups=2, group_select=sel)
        shard, group = _split(mod.lanes(row, report, {}))
        assert {ln.expected for ln in shard} == {(f + "t_light_full",)}
        assert set().union(*(set(ln.expected) for ln in group)) == {
            f + "t_heavy_full",
            f + "t_direct",
        }


def test_dash_without_shard0_only_runs_everything_in_every_shard(
    mod: ModuleType, plan: dict[str, Any]
) -> None:
    row = _row(plan, "compare", shards=3, groups=0, group_select="-")
    report = {"tests/calibration/test_c.py::t_a": False, "tests/calibration/test_c.py::t_b": False}
    out = mod.lanes(row, report, {})
    shard, group = _split(out)
    assert group == [] and [(ln.index, ln.n) for ln in shard] == [(0, 3), (1, 3), (2, 3)]
    assert all(set(ln.expected) == set(report) and not ln.deselect for ln in shard)
    assert [ln.part for ln in shard] == [f"task-14.shard{k}.junit.xml" for k in range(3)]


def test_fewer_group_tests_than_lanes_leaves_empty_lanes(mod: ModuleType,
                                                         plan: dict[str, Any]) -> None:  # fmt: skip
    row = _row(plan, "regions")  # 5 groups
    report = {"tests/calibration/test_r.py::t_a": False, "tests/calibration/test_r.py::t_b": False}
    out = mod.lanes(row, report, {"tests/calibration/test_r.py::t_a": 50.0})
    _, group = _split(out)
    assert [ln.index for ln in group] == list(range(5))
    assert [ln.expected for ln in group[:2]] == [
        ("tests/calibration/test_r.py::t_b",),
        ("tests/calibration/test_r.py::t_a",),
    ]
    assert all(ln.expected == () and ln.part is None for ln in group[2:])
    assert mod.record_parts(out, []) == [_part("regions", "group", 0), _part("regions", "group", 1)]


# ------------------------------------------------------------------ PlanErrors at collection


def test_prefix_collision_is_a_plan_error(mod: ModuleType, plan: dict[str, Any]) -> None:
    f = "tests/calibration/test_x.py::"
    row = _row(plan, "placement")
    report = {f + "test_x": True, f + "test_x_full": False}  # shard lanes deselect test_x
    with pytest.raises(mod.PlanError) as exc:
        mod.lanes(row, report, {})
    assert "test_x" in str(exc.value)
    # group lane collision: t_x and t_x_full on different lanes under '*'
    row = _row(plan, "gene_fdr")
    with pytest.raises(mod.PlanError):
        mod.lanes(row, {f + "t_x": False, f + "t_x_full": False}, {})
    # same lane: no collision
    assert mod.lanes(_row(plan, "detectability"), {f + "t_x": False, f + "t_x_full": False}, {})
    # the group lane keeping t_y_full deselects the shard test t_y: also a collision
    with pytest.raises(mod.PlanError):
        mod.lanes(_row(plan, "placement", groups=1), {f + "t_y": False, f + "t_y_full": True}, {})
    # no deselect at all ('-', no group lane): prefix pairs are fine
    out = mod.lanes(_row(plan, "compare", shards=2, groups=0, group_select="-"),
                    {f + "t_y": False, f + "t_y_full": False}, {})  # fmt: skip
    assert [set(ln.expected) for ln in out] == [{f + "t_y", f + "t_y_full"}] * 2


def test_no_test_left_for_shard_lanes_is_a_plan_error(mod: ModuleType,
                                                      plan: dict[str, Any]) -> None:  # fmt: skip
    row = _row(plan, "novel_orfs")
    report = {"tests/calibration/test_n.py::t_a": True, "tests/calibration/test_n.py::t_b": True}
    with pytest.raises(mod.PlanError):
        mod.lanes(row, report, {})
    with pytest.raises(mod.PlanError):
        mod.lanes(row, {}, {})  # empty collection, shards > 0


def test_prefix_list_with_zero_shards_never_drops_a_test(mod: ModuleType, plan: dict[str, Any],
                                                         tmp_path: Path) -> None:  # fmt: skip
    """Review point N2 of PR 1a."""
    sel = "tests/calibration/test_compare.py::"
    lines = PLAN.read_text().splitlines()
    lines[-1] = "\t".join(["compare", "0", "1", sel, "2", "4000", "02:00:00", "30"])
    path = tmp_path / "plan.tsv"
    path.write_text("\n".join(lines) + "\n")
    try:
        (row,) = [r for r in mod.read_plan(path, SUITES) if r.suite == "compare"]
    except mod.PlanError:
        return  # rejected in validation: also fine
    report = {sel + "t_a": False, "tests/calibration/test_other.py::t_b": False}
    with pytest.raises(mod.PlanError) as exc:
        mod.lanes(row, report, {})
    assert "test_other.py::t_b" in str(exc.value)


# ------------------------------------------------------------------ property (criterion 2)

FILES = ("tests/calibration/test_p.py::", "tests/calibration/attacks/test_q.py::")
PREFIXES = ("tests/calibration/test_p.py::", "tests/calibration/test_p.py::t_a",
            "tests/calibration/attacks/test_q.py::t_b")  # fmt: skip
_PAIRED = st.builds(lambda f, n, s: f"{f}t_{n}{s}", st.sampled_from(FILES),
                    st.text("ab", min_size=1, max_size=3),
                    st.sampled_from(["", "_full", "[p1]", "[p10]"]))  # fmt: skip
_FLAT = st.builds(lambda f, n, k: f"{f}t_{n}_{k}", st.sampled_from(FILES),
                  st.sampled_from(["aa", "ab", "ba", "bb"]), st.integers(10, 19))  # fmt: skip


@st.composite
def _case(draw: st.DrawFn) -> tuple[str, int, int, str, int, dict[str, bool], dict[str, float]]:
    ids = sorted(draw(st.sets(draw(st.sampled_from([_PAIRED, _FLAT])), max_size=80)))
    shards = draw(st.integers(0, 6))
    groups = draw(st.integers(0 if shards else 1, 8))
    if shards == 0:
        sel = "*"
    elif groups == 0:
        sel = "-"
    else:
        picks = draw(st.lists(st.sampled_from(PREFIXES), max_size=2, unique=True))
        sel = ",".join([*(["direct"] if draw(st.booleans()) or not picks else []), *picks])
    flags = {n: draw(st.booleans()) for n in ids}
    if sel == "-" and draw(st.booleans()):
        flags = dict.fromkeys(ids, False)
    durations = {n: float(draw(st.integers(1, 20000))) for n in ids if draw(st.booleans())}
    durations["tests/calibration/test_gone.py::t_removed"] = 99999.0  # not collected: ignored
    suite = draw(st.sampled_from(SUITES))
    return suite, shards, groups, sel, draw(st.integers(1, 120)), flags, durations


def _canon(out: list[Any]) -> list[tuple[Any, ...]]:
    return [(ln.kind, ln.index, ln.n, frozenset(ln.deselect), ln.part, frozenset(ln.expected))
            for ln in out]  # fmt: skip


@settings(max_examples=400, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(_case(), st.data())
def test_lanes_property(mod: ModuleType, plan: dict[str, Any], case: Any,
                        data: st.DataObject) -> None:  # fmt: skip
    suite, shards, groups, sel, unknown, flags, durations = case
    row = _row(plan, suite, shards=shards, groups=groups, group_select=sel,
               unknown_min=float(unknown))  # fmt: skip
    ids = sorted(flags)
    try:
        want_shard, want_groups = _model(row, flags, durations)
    except _Err:
        with pytest.raises(mod.PlanError):
            mod.lanes(row, flags, durations)
        return
    if not ids and shards == 0:
        try:
            out = mod.lanes(row, flags, durations)
        except mod.PlanError:
            return
        assert all(ln.expected == () for ln in out)
        return
    out = mod.lanes(row, flags, durations)
    shard, group = _split(out)
    assert [(ln.index, ln.n) for ln in shard] == [(k, shards) for k in range(shards)]
    assert [(ln.index, ln.n) for ln in group] == [(g, groups) for g in range(groups)]
    i = SUITES.index(suite)
    assert [ln.part for ln in shard] == [f"task-{i}.shard{k}.junit.xml" for k in range(shards)]
    for g, ln in enumerate(group):
        assert ln.part == (f"task-{i}.group{g}.junit.xml" if ln.expected else None)

    # pytest's startswith deselection gives exactly `expected` in every non-empty lane
    for ln in out:
        if ln.expected:
            assert _kept(ids, ln.deselect) == set(ln.expected)
    # every nodeid in exactly one group lane or in every shard lane, never both
    group_sets = [set(ln.expected) for ln in group]
    for n in ids:
        in_groups = sum(n in s for s in group_sets)
        in_shards = sum(n in ln.expected for ln in shard)
        assert (in_groups, in_shards) in ((1, 0), (0, shards)) and in_groups + in_shards > 0, n
    if shards > 0:
        assert all(n in set().union(*group_sets) for n in ids if flags[n])
        assert {frozenset(ln.expected) for ln in shard} == {frozenset(want_shard)}
    assert group_sets == [set(w) for w in want_groups]  # the LPT assignment, exactly

    # LPT bound (list scheduling): max load <= mean load + longest test
    def w(n: str) -> float:
        return durations.get(n, unknown * 60.0)

    group_tests = [n for s in group_sets for n in s]
    if group_tests:
        loads = [sum(w(n) for n in s) for s in group_sets]
        assert max(loads) <= sum(loads) / groups + max(w(n) for n in group_tests) + 1e-6

    # output does not depend on input order
    perm = data.draw(st.permutations(ids))
    dperm = data.draw(st.permutations(sorted(durations)))
    again = mod.lanes(row, {n: flags[n] for n in perm}, {n: durations[n] for n in dperm})
    assert again == out
    assert _canon(again) == _canon(out)


# ------------------------------------------------------------------ verify_partition (crit. 3)


def _collected(out: list[Any], report: Mapping[str, bool]) -> dict[str, list[str]]:
    got: dict[str, list[str]] = {"all": sorted(report)}
    for ln in out:
        if ln.part is not None:
            got[ln.part] = sorted(_kept(sorted(report), ln.deselect))
    return got


def test_verify_partition_accepts_examples(mod: ModuleType, plan: dict[str, Any]) -> None:
    for row, report, durations in (_example2(plan), _example3(plan), _lpt_example(plan)):
        out = mod.lanes(row, report, durations)
        assert mod.verify_partition(row, out, _collected(out, report)) == []


def _break(case: str, got: dict[str, list[str]], out: list[Any]) -> None:
    s0, s1 = "task-0.shard0.junit.xml", "task-0.shard1.junit.xml"
    g0, g1 = "task-0.group0.junit.xml", "task-0.group1.junit.xml"
    if case == "shard_lanes_differ":
        got[s1] = got[s1][1:]
    elif case == "shard_lane_extra":  # both shard lanes run a group test too
        got[s0] = sorted([*got[s0], got[g0][0]])
        got[s1] = list(got[s0])
    elif case == "group_lane_short":
        got[g0] = got[g0][1:]
    elif case == "group_lane_extra":
        got[g0] = sorted([*got[g0], got[s0][0]])
    elif case == "groups_overlap":
        got[g1] = sorted([*got[g1], got[g0][0]])
    elif case == "union_short_of_collection":
        got["all"] = sorted([*got["all"], "tests/calibration/test_new.py::t_unplanned"])
    elif case == "lane_not_collected":
        del got[g1]
    elif case == "all_missing":
        del got["all"]


@pytest.mark.parametrize("case", [
    "shard_lanes_differ", "shard_lane_extra", "group_lane_short", "group_lane_extra",
    "groups_overlap", "union_short_of_collection", "lane_not_collected", "all_missing",
])  # fmt: skip
def test_verify_partition_reports_each_problem(mod: ModuleType, plan: dict[str, Any],
                                               case: str) -> None:  # fmt: skip
    row, report, durations = _example3(plan)
    out = mod.lanes(row, report, durations)
    got = _collected(out, report)
    _break(case, got, out)
    problems = mod.verify_partition(row, out, got)
    assert isinstance(problems, list) and problems, case
    assert all(isinstance(p, str) and p for p in problems)


# ------------------------------------------------------------------ verify_counts (crit. 3)


def test_verify_counts(mod: ModuleType, plan: dict[str, Any]) -> None:
    out = mod.lanes(*_fewer(plan))
    counts = {ln.part: len(ln.expected) for ln in out if ln.part is not None}
    assert mod.verify_counts(out, counts) == []
    bad = dict(counts)
    part = next(iter(bad))
    bad[part] -= 1
    problems = mod.verify_counts(out, bad)
    assert problems and any(part in p for p in problems)
    missing = dict(counts)
    del missing[part]
    assert mod.verify_counts(out, missing)
    assert mod.verify_counts(out, {**counts, part: counts[part] + 1})


def _fewer(plan: Mapping[str, Any]) -> tuple[Any, dict[str, bool], dict[str, float]]:
    f = "tests/calibration/test_v.py::"
    report = {f + "t_h1": False, f + "t_h2": False, f + "t_d": True}
    return _row(plan, "variants", groups=3), report, {}  # 2 shard lanes, 1 used + 2 empty groups


# ------------------------------------------------------------------ durations


def _junit(path: Path, cases: list[tuple[str, str, str]]) -> Path:
    body = "".join(f'<testcase classname="{c}" name="{n}" time="{t}" />' for c, n, t in cases)
    path.write_text('<?xml version="1.0" encoding="utf-8"?><testsuites name="pytest tests">'
                    f'<testsuite name="pytest" tests="{len(cases)}">{body}</testsuite>'
                    "</testsuites>")  # fmt: skip
    return path


def test_durations_from_junit(mod: ModuleType, tmp_path: Path) -> None:
    a = _junit(tmp_path / "task-10.run.junit.xml", [
        ("tests.calibration.attacks.test_attack_determinants_null",
         "test_null_level_under_attack[full-autocorr_pos-lightgbm]", "8364.120"),
        ("tests.calibration.test_regions_calibration", "test_c1_level_full[rho0.0]", "4140.5"),
    ])  # fmt: skip
    b = tmp_path / "task-2.run.junit.xml"
    b.write_text(
        '<?xml version="1.0" encoding="utf-8"?><testsuites><testsuite name="pytest">'
        '<testcase classname="tests.calibration.test_d" name="test_fails" time="12.0">'
        '<failure message="x">boom</failure></testcase></testsuite></testsuites>'
    )
    got = mod.durations_from_junit([a, b])
    assert got == {
        f"{DET}[full-autocorr_pos-lightgbm]": pytest.approx(8364.12),
        "tests/calibration/test_regions_calibration.py::test_c1_level_full[rho0.0]": 4140.5,
        "tests/calibration/test_d.py::test_fails": 12.0,
    }
    assert mod.durations_from_junit([]) == {}


def test_read_durations(mod: ModuleType, tmp_path: Path) -> None:
    p = tmp_path / "d.tsv"
    p.write_text("nodeid\tseconds\tsource\n"
                 f"{DET}[full-detect-lightgbm]\t7782\t1894441\n"
                 "tests/calibration/test_r.py::t_a[x-1.5]\t0.25\tprobe-1895352\n")  # fmt: skip
    assert mod.read_durations(p) == {
        f"{DET}[full-detect-lightgbm]": 7782.0,
        "tests/calibration/test_r.py::t_a[x-1.5]": 0.25,
    }
    p.write_text("nodeid\tseconds\tsource\n")
    assert mod.read_durations(p) == {}
    for bad in ("-3", "x", ""):
        p.write_text(f"nodeid\tseconds\tsource\ntests/calibration/test_r.py::t\t{bad}\tj\n")
        with pytest.raises(ValueError):
            mod.read_durations(p)


def test_committed_durations_file(mod: ModuleType) -> None:
    assert DURATIONS.is_file(), f"{DURATIONS} missing (M11-03c PR 1b data)"
    lines = DURATIONS.read_text().splitlines()
    assert lines[0].split("\t") == ["nodeid", "seconds", "source"]
    ids = [ln.split("\t")[0] for ln in lines[1:]]
    assert len(ids) == len(set(ids)) and all(n.startswith("tests/calibration/") for n in ids)
    assert all(ln.split("\t")[2] for ln in lines[1:])
    got = mod.read_durations(DURATIONS)
    assert round(got[f"{DET}[full-autocorr_pos-lightgbm]"]) == 8364  # example 2
    assert round(got[f"{DET}[full-detect-lightgbm]"]) == 7782
    rel = DURATIONS.relative_to(REPO).as_posix()  # outside CALIB_PATHS (Q6)
    assert not rel.startswith(("src/", "tests/calibration/", "ci/", "tests/sim/", "tests/helpers/"))


# ------------------------------------------------------------------ CLI


def _cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True,
                          timeout=60, check=False)  # fmt: skip


def test_lanes_cli_reads_the_conftest_lane_report(mod: ModuleType, plan: dict[str, Any],
                                                  tmp_path: Path) -> None:  # fmt: skip
    row, report, durations = _example2(plan)
    rep = tmp_path / "report.json"
    rep.write_text(json.dumps({n: {"shard0_only": s} for n, s in report.items()}))
    dur = tmp_path / "d.tsv"
    dur.write_text("nodeid\tseconds\tsource\n"
                   + "".join(f"{n}\t{s}\tt\n" for n, s in durations.items()))  # fmt: skip
    res = _cli("lanes", "--suite", "determinants", "--report", str(rep), "--durations", str(dur))
    assert res.returncode == 0, res.stderr
    got = json.loads(res.stdout)
    want = mod.lanes(row, report, durations)
    assert [(d["kind"], d["index"], d["n"], d["part"], set(d["expected"]), set(d["deselect"]))
            for d in got] == [(ln.kind, ln.index, ln.n, ln.part, set(ln.expected),
                               set(ln.deselect)) for ln in want]  # fmt: skip


def test_lanes_cli_exits_2_on_plan_error(tmp_path: Path) -> None:
    rep = tmp_path / "report.json"
    rep.write_text(json.dumps({"tests/calibration/test_n.py::t_a": {"shard0_only": True}}))
    res = _cli("lanes", "--suite", "novel_orfs", "--report", str(rep), "--durations",
               str(DURATIONS))  # fmt: skip
    assert res.returncode == 2 and "novel_orfs" in res.stderr


def test_durations_cli_round_trip(mod: ModuleType, tmp_path: Path) -> None:
    a = _junit(
        tmp_path / "a.junit.xml",
        [
            ("tests.calibration.test_a", "test_x[1]", "3.5"),
            ("tests.calibration.test_a", "test_y", "70.0"),
        ],
    )
    out = tmp_path / "durations.tsv"
    res = _cli("durations", "--source", "1894441", "--out", str(out), str(a))
    assert res.returncode == 0, res.stderr
    assert out.read_text().splitlines()[0].split("\t") == ["nodeid", "seconds", "source"]
    assert mod.read_durations(out) == mod.durations_from_junit([a])
    assert all(ln.endswith("\t1894441") for ln in out.read_text().splitlines()[1:])
