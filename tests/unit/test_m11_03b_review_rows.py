"""#349 review round of M11-03b: ``suite_result`` and ``full_rows`` (B1, B2, S1, S2).

Spec: docs/specs/M11-03b-nightly-readiness-labels.md, "Lead decisions, #349 review (2026-10-25)".

- B1: a crash exit (unsharded exit not in {0, 1}; sharded exit >= 2) fails **both** rows whatever
  the suite's attack-test count. The collection-error JUnit below is what pytest writes (checked
  with the venv's pytest): one ``<error message="collection failure">`` testcase with
  ``classname=""`` and the module path as ``name``.
- B2: ``full_rows(..., expect=15)``; too few suites, 0 attack suites or 0 attack tests, or an
  attack suite whose attack tests are all skipped fail the row; ``full_rows([])`` never passes.
  Suites without attack tests are data (today ``denovo_precision``): any other suite that holds no
  attack testcase fails the attack row.
- S1: testcases are de-duplicated by (classname, name) across a task's JUnit parts. Task 13 writes
  the same testcases in ``shard0`` (direct checks, verdicts) and ``shard1`` (``SHARD_SKIP``); the
  merge's block-2 confirm run, when it happens, writes them again. A verdict beats a no-verdict
  skip and a failure is never hidden by a skip in another part.
- S2: direct ``suite_result`` tests for the fail-closed paths (one per surviving mutant of the
  review's mutation run), so the ``full_rows`` fallback cannot mask a rule bug.
"""

from __future__ import annotations

import copy
import re
from types import ModuleType
from typing import Any

import pytest
from helpers.calibration_rows_fixture import (
    DATE,
    JOB,
    NO_ATTACKS,
    SHARDED,
    SUITES,
    C,
    Case,
    World,
    by_label,
    junit_xml,
    load_rows_module,
    plain_class,
    record,
    split_row,
)
from helpers.release_readiness import load_module as load_readiness

EXPECT = 15
EX1_FULL = ("calibration_full: 15 suites, 0 failed suites; 422 tests: 0 failed, 0 errors, "
            "400 passed, 20 skipped, 2 xfailed, 0 xpassed; array 1900000")  # fmt: skip
EX1_ATTACKS_COUNTS = ("230 attack tests: 0 failed, 0 errors, 220 passed, 10 skipped, 0 xfailed, "
                      "0 xpassed; array 1900000")  # fmt: skip
COUNT = re.compile(r"(\d+) (?:failed|errors?)\b")
ISLANDS = 5


@pytest.fixture(scope="module")
def mod() -> ModuleType:
    return load_rows_module()


def collection_error_xml(module: str) -> str:
    """pytest's JUnit after a collection error (exit 2): one error testcase, ``classname=""``."""
    return ('<?xml version="1.0" encoding="utf-8"?><testsuites name="pytest tests">'
            '<testsuite name="pytest" errors="1" failures="0" skipped="0" tests="1" time="0.111" '
            'timestamp="2026-11-02T03:00:00.000000-08:00" hostname="barley-01">'
            f'<testcase classname="" name="{module}" time="0.000">'
            '<error message="collection failure">ImportError while importing test module '
            f"'/farmshare/ci/calibration/{module.replace('.', '/')}.py'.\n"
            "E   ModuleNotFoundError: No module named 'nosuchmod'</error></testcase>"
            "</testsuite></testsuites>")  # fmt: skip


ATTACK_MODULE = "tests.calibration.attacks.test_attack_islands"


def rows(mod: ModuleType, results: list[Any], expect: int = EXPECT) -> list[str]:
    return list(mod.full_rows(results, JOB, C, DATE, expect=expect))


def wrows(mod: ModuleType, w: World) -> dict[str, tuple[int, str]]:
    return by_label(rows(mod, w.results(mod)))


def readiness_ok(lines: list[str]) -> dict[str, bool]:
    parsed, bad = load_readiness().parse_log("".join(lines))
    assert bad == 0, lines
    return {r.label: r.ok() for r in parsed}


def failed_suites(summary: str) -> int:
    m = re.search(r", (\d+) failed suites\b", summary)
    assert m, summary
    return int(m.group(1))


def assert_red(ex: int, summary: str) -> None:
    """A failing row: exit 1 **and** ``failed suites`` >= 1 (readiness false by exit and count)."""
    assert ex == 1, summary
    assert failed_suites(summary) >= 1, summary


# ------------------------------------------------------------------ B1: crash exits


def collection_error_world(task: int, module: str, exit_: int = 2) -> World:
    w = World.default()
    w.records[task] = record(task, exit_=exit_)
    texts: dict[str, str | None] = {p: None for p in w.records[task]["parts"]}  # type: ignore[index]
    first = w.records[task]["parts"][0]  # type: ignore[index]
    texts[first] = collection_error_xml(module)
    w.records[task]["parts"] = [first]  # type: ignore[index]
    w.junit_override[task] = texts
    return w


def test_b1_suite_result_collection_error_fails_the_attacks(mod: ModuleType) -> None:
    w = collection_error_world(ISLANDS, ATTACK_MODULE)
    res = mod.suite_result(w.records[ISLANDS], w.junits(ISLANDS), ISLANDS, C)
    assert res.problems, res
    assert res.attacks_failed is True, res
    assert res.attack_suite is True, res  # counted in the attack row, not dropped from it


def test_b1_one_collection_error_fails_both_rows(mod: ModuleType) -> None:
    lines = rows(mod, collection_error_world(ISLANDS, ATTACK_MODULE).results(mod))
    got = by_label(lines)
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert_red(ex, s)
        assert "islands" in s, s
    assert got["calibration_attacks"][1].startswith("calibration_attacks: 14 suites, "), got
    assert readiness_ok(lines) == {"calibration_full": False, "calibration_attacks": False}


def test_b1_all_fifteen_collection_errors_fail_both_rows(mod: ModuleType) -> None:
    """The review's reproduction: one ImportError under tests/calibration stops every suite."""
    w = World.default()
    for t in range(EXPECT):
        w.records[t] = record(t, exit_=2)
        w.records[t]["parts"] = [w.records[t]["parts"][0]]  # type: ignore[index]
        w.junit_override[t] = {p: None for p in record(t)["parts"]}
        w.junit_override[t][w.records[t]["parts"][0]] = collection_error_xml(ATTACK_MODULE)  # type: ignore[index]
    lines = rows(mod, w.results(mod))
    got = by_label(lines)
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert_red(ex, s)
    ex, s = got["calibration_attacks"]
    assert failed_suites(s) >= 14, s
    assert readiness_ok(lines) == {"calibration_full": False, "calibration_attacks": False}


@pytest.mark.parametrize("exit_", [2, 3, 4, 5, 143])
def test_b1_unsharded_crash_of_the_attack_free_suite_fails_the_attack_row(
    mod: ModuleType, exit_: int
) -> None:
    """``denovo_precision`` holds no attack test, but a crash leaves its JUnit partial."""
    w = collection_error_world(NO_ATTACKS, "tests.calibration.test_denovo_precision_calibration",
                               exit_)  # fmt: skip
    res = mod.suite_result(w.records[NO_ATTACKS], w.junits(NO_ATTACKS), NO_ATTACKS, C)
    assert res.attacks_failed is True and res.attack_suite is True, res
    got = by_label(rows(mod, w.results(mod)))
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert_red(ex, s)
        assert "denovo_precision" in s, s


@pytest.mark.parametrize("exit_", [2, 3])
def test_b1_sharded_exit_2_without_attack_testcases_fails_both_rows(
    mod: ModuleType, exit_: int
) -> None:
    """Task 13 merge error (exit 2) whose JUnit has no attack testcase: still fail closed."""
    w = World.default()
    w.records[SHARDED] = record(SHARDED, exit_=exit_)
    w.cases[SHARDED] = [c for c in w.cases[SHARDED] if not c.attack]
    res = mod.suite_result(w.records[SHARDED], w.junits(SHARDED), SHARDED, C)
    assert res.attacks_failed is True and res.attack_suite is True, res
    got = by_label(rows(mod, w.results(mod)))
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert_red(ex, s)
        assert "variants" in s, s


def test_b1_sharded_collection_error_fails_both_rows(mod: ModuleType) -> None:
    w = World.default()
    w.records[SHARDED] = record(SHARDED, exit_=2)
    shards = ["task-13.shard0.junit.xml", "task-13.shard1.junit.xml"]
    w.records[SHARDED]["parts"] = shards  # type: ignore[index]
    err = collection_error_xml("tests.calibration.attacks.test_attack_variants")
    w.junit_override[SHARDED] = {shards[0]: err, shards[1]: err, "task-13.merge.junit.xml": None}
    got = by_label(rows(mod, w.results(mod)))
    for label in ("calibration_full", "calibration_attacks"):
        assert_red(*got[label])


def test_b1_sharded_exit_1_without_attack_testcases_still_spares_the_attack_row(
    mod: ModuleType,
) -> None:
    """Example 6 is unchanged by B1 (merge exit 1 is a confirmed fail, not a crash)."""
    w = World.default()
    w.records[SHARDED] = record(SHARDED, exit_=1)
    w.cases[SHARDED] = [c for c in w.cases[SHARDED] if not c.attack]
    got = by_label(rows(mod, w.results(mod)))
    assert got["calibration_full"][0] == 1
    assert got["calibration_attacks"][0] == 0, got


# ------------------------------------------------------------------ B2: expect and vacuity


def test_b2_example1_passes_with_expect(mod: ModuleType) -> None:
    got = wrows(mod, World.default())
    assert got["calibration_full"] == (0, EX1_FULL)
    ex, s = got["calibration_attacks"]
    assert ex == 0, s
    assert s.startswith("calibration_attacks: 14 suites, 0 failed suites"), s
    assert s.endswith(EX1_ATTACKS_COUNTS), s


def test_b2_empty_results_never_pass(mod: ModuleType) -> None:
    lines = rows(mod, [])
    got = by_label(lines)
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert ex != 0, s
    assert readiness_ok(lines) == {"calibration_full": False, "calibration_attacks": False}


@pytest.mark.parametrize("drop", [0, ISLANDS, NO_ATTACKS, SHARDED, 14])
def test_b2_too_few_suites_fail_the_full_row(mod: ModuleType, drop: int) -> None:
    results = [r for r in World.default().results(mod) if r.task != drop]
    lines = rows(mod, results)
    ex, s = by_label(lines)["calibration_full"]
    assert ex == 1, s
    assert readiness_ok(lines)["calibration_full"] is False, s


def test_b2_expect_is_honoured_not_hard_coded(mod: ModuleType) -> None:
    results = World.default().results(mod)
    assert by_label(rows(mod, results, expect=16))["calibration_full"][0] == 1
    assert by_label(rows(mod, results, expect=15))["calibration_full"][0] == 0


def test_b2_misclassified_attacks_fail_the_attack_row(mod: ModuleType) -> None:
    """0 attack suites and 0 attack tests (e.g. a classname prefix change): never a pass."""
    w = World.default()
    for t, cases in w.cases.items():
        for c in cases:
            if c.attack:
                c.classname = plain_class(t)
    lines = rows(mod, w.results(mod))
    got = by_label(lines)
    ex, s = got["calibration_attacks"]
    assert ex == 1, s
    assert " 0 attack tests" in s, s
    assert readiness_ok(lines)["calibration_attacks"] is False, s


def test_b2_attack_suite_without_attack_testcases_fails_the_attack_row(mod: ModuleType) -> None:
    """islands is not attack-free data: losing every attack testcase fails, never drops out."""
    w = World.default()
    w.cases[ISLANDS] = [c for c in w.cases[ISLANDS] if not c.attack]
    got = by_label(rows(mod, w.results(mod)))
    ex, s = got["calibration_attacks"]
    assert_red(ex, s)
    assert "islands" in s, s


def test_b2_the_attack_free_suite_is_data_and_passes(mod: ModuleType) -> None:
    """denovo_precision holds no attack testcase in example 1 and the attack row passes."""
    w = World.default()
    assert not any(c.attack for c in w.cases[NO_ATTACKS])
    res = mod.suite_result(w.records[NO_ATTACKS], w.junits(NO_ATTACKS), NO_ATTACKS, C)
    assert res.attacks_failed is False and res.problems == [], res
    assert wrows(mod, w)["calibration_attacks"][0] == 0


def test_b2_all_skipped_attack_tests_fail_the_attack_row(mod: ModuleType) -> None:
    w = World.default()
    for c in w.cases[ISLANDS]:
        if c.attack:
            c.outcome = "skipped"
    got = by_label(rows(mod, w.results(mod)))
    assert got["calibration_full"][0] == 0, got  # a skip is not a failure of the full tier
    ex, s = got["calibration_attacks"]
    assert_red(ex, s)
    assert "islands" in s, s


def test_b2_one_attack_verdict_is_enough(mod: ModuleType) -> None:
    """Not vacuous in the other direction: all but one attack test of a suite skipped passes."""
    w = World.default()
    att = [c for c in w.cases[ISLANDS] if c.attack]
    for c in att[1:]:
        c.outcome = "skipped"
    att[0].outcome = "passed"
    assert wrows(mod, w)["calibration_attacks"][0] == 0


# ------------------------------------------------------------------ S1: task 13 de-duplication

SHARD0, SHARD1, MERGE = (f"task-{SHARDED}.{p}.junit.xml" for p in ("shard0", "shard1", "merge"))


def skipped_copies(cases: list[Case]) -> list[Case]:
    return [Case(c.classname, c.name, "skipped") for c in cases]


def clean_sharded_world() -> World:
    """A passing task 13 as pytest writes it: shard0 holds the verdicts of the direct checks,
    shard1 the same testcases skipped (``SHARD_SKIP``); no merge JUnit (block 2 did not run)."""
    w = World.default()
    w.records[SHARDED]["parts"] = [SHARD0, SHARD1]  # type: ignore[index]
    w.junit_override[SHARDED] = {SHARD1: junit_xml(skipped_copies(w.cases[SHARDED])),
                                 MERGE: None}  # fmt: skip
    return w


def test_s1_shard_files_are_not_double_counted(mod: ModuleType) -> None:
    w = clean_sharded_world()
    cases = w.cases[SHARDED]
    res = mod.suite_result(w.records[SHARDED], w.junits(SHARDED), SHARDED, C)
    assert res.problems == [], res
    assert res.all.total == len(cases), res.all
    assert res.attacks.total == sum(c.attack for c in cases), res.attacks


def test_s1_a_verdict_beats_a_shard_skip(mod: ModuleType) -> None:
    w = clean_sharded_world()
    cases = w.cases[SHARDED]
    res = mod.suite_result(w.records[SHARDED], w.junits(SHARDED), SHARDED, C)
    assert res.all.passed == sum(c.outcome == "passed" for c in cases), res.all
    assert res.attacks.passed == sum(c.attack and c.outcome == "passed" for c in cases)


def test_s1_rows_of_a_clean_sharded_run_equal_example1(mod: ModuleType) -> None:
    got = wrows(mod, clean_sharded_world())
    assert got["calibration_full"] == (0, EX1_FULL), got
    ex, s = got["calibration_attacks"]
    assert ex == 0 and s.endswith(EX1_ATTACKS_COUNTS), s


def test_s1_merge_confirm_run_rewrites_the_testcases(mod: ModuleType) -> None:
    """Block 2 ran: shard0 and shard1 skip, the merge JUnit holds the verdicts (one failed)."""
    w = World.default()
    cases = w.cases[SHARDED]
    w.records[SHARDED] = record(SHARDED, exit_=1)
    verdicts = copy.deepcopy(cases)
    plain = next(c for c in verdicts if not c.attack and c.outcome == "passed")
    plain.outcome = "failed"
    w.junit_override[SHARDED] = {SHARD0: junit_xml(skipped_copies(cases)),
                                 SHARD1: junit_xml(skipped_copies(cases)),
                                 MERGE: junit_xml(verdicts)}  # fmt: skip
    res = mod.suite_result(w.records[SHARDED], w.junits(SHARDED), SHARDED, C)
    assert res.all.total == len(cases), res.all
    assert res.all.failed == 1, res.all
    assert res.attacks.total == sum(c.attack for c in cases), res.attacks
    ex, s = by_label(rows(mod, w.results(mod)))["calibration_full"]
    assert ex == 1 and " 422 tests: 1 failed, 0 errors," in s, s


def test_s1_a_failure_in_a_later_part_is_not_hidden_by_an_earlier_skip(mod: ModuleType) -> None:
    w = clean_sharded_world()
    cases = w.cases[SHARDED]
    w.records[SHARDED] = record(SHARDED, exit_=1)
    w.records[SHARDED]["parts"] = [SHARD0, SHARD1]  # type: ignore[index]
    s0 = skipped_copies(cases)
    s1 = skipped_copies(cases)
    s1[0].outcome = "failed"
    s1[1].outcome = "error"
    w.junit_override[SHARDED] = {SHARD0: junit_xml(s0), SHARD1: junit_xml(s1), MERGE: None}
    res = mod.suite_result(w.records[SHARDED], w.junits(SHARDED), SHARDED, C)
    assert res.all.total == len(cases), res.all
    assert (res.all.failed, res.all.errors) == (1, 1), res.all


def test_s1_same_name_in_another_class_is_a_different_testcase(mod: ModuleType) -> None:
    w = World.default()
    cases = w.cases[3]
    twin = Case(cases[0].classname + "_twin", cases[0].name, "passed")
    w.cases[3] = [*cases, twin]
    res = mod.suite_result(record(3), w.junits(3), 3, C)
    assert res.all.total == len(cases) + 1, res.all


# ------------------------------------------------------------------ S2: fail-closed paths


def base(task: int = 3) -> tuple[dict[str, Any], list[Any]]:
    w = World.default()
    return record(task), w.junits(task)


def test_s2_task_mismatch(mod: ModuleType) -> None:  # mutant 3
    rec, junits = base(3)
    rec["task"] = 4
    res = mod.suite_result(rec, junits, 3, C)
    assert any("task mismatch" in p for p in res.problems), res.problems
    assert res.attacks_failed is True


@pytest.mark.parametrize(
    ("key", "value"),
    [("args", None), ("commit", None), ("exit", "0"), ("dirty", "false"), ("parts", None),
     ("args", "--fast"), ("dirty", 0), ("exit", 0.0)],
)  # fmt: skip
def test_s2_malformed_record(mod: ModuleType, key: str, value: object) -> None:  # mutant 4
    rec, junits = base(3)
    rec[key] = value
    res = mod.suite_result(rec, junits, 3, C)
    assert "record malformed" in res.problems, res.problems
    assert res.attacks_failed is True


def test_s2_malformed_record_fails_both_rows(mod: ModuleType) -> None:
    w = World.default()
    w.records[3]["args"] = None  # type: ignore[index]
    got = wrows(mod, w)
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert_red(ex, s)
        assert "regions record malformed" in s, s


@pytest.mark.parametrize("value", [False, True])
def test_s2_exit_bool_is_not_an_int(mod: ModuleType, value: bool) -> None:  # mutant 28
    rec, junits = base(3)
    rec["exit"] = value
    res = mod.suite_result(rec, junits, 3, C)
    assert res.problems, res
    assert res.attacks_failed is True


def test_s2_task_true_is_not_task_1(mod: ModuleType) -> None:
    rec, junits = record(1), World.default().junits(1)
    rec["task"] = True
    res = mod.suite_result(rec, junits, 1, C)
    assert res.problems, res
    assert res.attacks_failed is True


def test_s2_suite_mismatch(mod: ModuleType) -> None:  # mutant 27
    rec, junits = base(3)
    rec["suite"] = SUITES[ISLANDS]
    res = mod.suite_result(rec, junits, 3, C, suite=SUITES[3])
    assert any("suite mismatch" in p for p in res.problems), res.problems
    assert res.attacks_failed is True


@pytest.mark.parametrize(
    "part", ["../task-3.run.junit.xml", "sub/task-3.run.junit.xml", "/tmp/task-3.run.junit.xml",
             "..", ""],
)  # fmt: skip
def test_s2_path_traversal_in_parts(mod: ModuleType, part: str) -> None:  # mutant 26
    rec, junits = base(3)
    rec["parts"] = [part]
    res = mod.suite_result(rec, junits, 3, C)
    assert res.problems, (part, res)
    assert res.attacks_failed is True


def test_s2_unsharded_exit_0_with_a_failed_testcase(mod: ModuleType) -> None:  # mutants 8, 15
    w = World.default()
    w.set_outcome(3, attack=False, outcome="failed")
    res = mod.suite_result(record(3), w.junits(3), 3, C)
    assert "failed testcases" in res.problems, res.problems


def test_s2_unsharded_exit_0_with_a_failed_attack_testcase(mod: ModuleType) -> None:
    w = World.default()
    w.set_outcome(3, attack=True, outcome="failed")
    res = mod.suite_result(record(3), w.junits(3), 3, C)
    assert "failed testcases" in res.problems, res.problems
    assert res.attacks_failed is True


def test_s2_sharded_exit_1_non_attack_failure_with_attack_tests(mod: ModuleType) -> None:  # 11
    w = World.default()
    w.set_outcome(SHARDED, attack=False, outcome="failed")
    assert any(c.attack for c in w.cases[SHARDED])
    res = mod.suite_result(record(SHARDED, exit_=1), w.junits(SHARDED), SHARDED, C)
    assert res.problems, res
    assert res.attacks_failed is True, res
    assert res.attacks.bad == 0  # the attack row fails closed, not by an attack count


def test_s2_unsharded_exit_2_with_a_non_attack_failure(mod: ModuleType) -> None:  # mutant 12
    w = World.default()
    w.set_outcome(3, attack=False, outcome="failed")
    res = mod.suite_result(record(3, exit_=2), w.junits(3), 3, C)
    assert res.problems, res
    assert res.attacks_failed is True, res


def test_s2_attack_error_only_with_a_non_attack_failure(mod: ModuleType) -> None:  # mutant 10
    w = World.default()
    w.set_outcome(3, attack=False, outcome="failed")
    w.set_outcome(3, attack=True, outcome="error")
    res = mod.suite_result(record(3, exit_=1), w.junits(3), 3, C)
    assert (res.attacks.failed, res.attacks.errors) == (0, 1), res.attacks
    assert res.attacks_failed is True, res
    assert res.problems, res


def test_s2_non_attack_failure_exit_1_spares_the_attacks(mod: ModuleType) -> None:
    """The counterpart of the previous tests: example 2 stays attack-clean."""
    w = World.default()
    w.set_outcome(3, attack=False, outcome="failed")
    res = mod.suite_result(record(3, exit_=1), w.junits(3), 3, C)
    assert res.problems and res.attacks_failed is False, res


@pytest.mark.parametrize(
    "args",
    [["3 failed"], ["--k", "1 error"], ["2 errors"], ["7\tfailed"], ["5\nfailed"], ["9  errors"]],
)  # fmt: skip
def test_s2_count_like_args_do_not_read_as_counts(mod: ModuleType, args: list[str]) -> None:
    w = World.default()
    w.records[3] = record(3, args=args)
    _check_counts_honest(mod, rows(mod, w.results(mod)))


@pytest.mark.parametrize("state", ["1 failed", "OUT_OF_MEMORY 3 errors", "0 failed\t4 error"])
def test_s2_count_like_missing_reason_does_not_read_as_counts(mod: ModuleType, state: str) -> None:
    results = World.default().results(mod)
    results[6] = mod.suite_result(None, [], 6, C, suite=SUITES[6], missing=f"no record, {state}")
    _check_counts_honest(mod, rows(mod, results))


def _check_counts_honest(mod: ModuleType, lines: list[str]) -> None:
    """Every ``<n> failed|errors`` match is a structured count (``<f> failed suites`` and the
    two test counts); the row keeps five fields and one line; it is red by exit and by count."""
    for line in lines:
        _d, _j, _c, ex, summary = split_row(line)
        tail = summary.rsplit("; ", 2)
        assert len(tail) == 3 and tail[2] == f"array {JOB}", summary
        tests = COUNT.findall(tail[1])
        assert len(tests) == 2, summary
        assert COUNT.findall(summary) == [str(failed_suites(summary)), *tests], summary
        assert_red(ex, summary)
    assert readiness_ok(lines) == {"calibration_full": False, "calibration_attacks": False}
