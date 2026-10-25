"""``scripts/ci/calibration_rows.py``: ``suite_result`` and ``full_rows`` on fixture records and
hand-written JUnit XML (issue M11-03b, examples 1-8; acceptance criteria 2 and 3).

Spec: docs/specs/M11-03b-nightly-readiness-labels.md. The world of example 1 is
``helpers.calibration_rows_fixture.World.default()``; every test changes one task. The functions
are pure: they get parsed records (``dict | None``) and parsed JUnit roots (``ET.Element``).
"""

from __future__ import annotations

import re
from types import ModuleType

import pytest
from helpers.calibration_rows_fixture import (
    DATE,
    JOB,
    NO_ATTACKS,
    OTHER,
    PREFIX,
    SHARDED,
    SUITES,
    C,
    World,
    by_label,
    load_rows_module,
    record,
    split_row,
)
from hypothesis import given, settings
from hypothesis import strategies as st

EX1_FULL = ("calibration_full: 15 suites, 0 failed suites; 422 tests: 0 failed, 0 errors, "
            "400 passed, 20 skipped, 2 xfailed, 0 xpassed; array 1900000")  # fmt: skip
EX1_ATTACKS = ("calibration_attacks: 14 suites, 0 failed suites; 230 attack tests: 0 failed, "
               "0 errors, 220 passed, 10 skipped, 0 xfailed, 0 xpassed; array 1900000")  # fmt: skip
COUNT = re.compile(r"(\d+) (?:failed|errors?)\b")


@pytest.fixture(scope="module")
def mod() -> ModuleType:
    return load_rows_module()


def test_suite_names_of_the_fixture_match_the_spec() -> None:
    assert len(SUITES) == 15
    assert (SUITES[3], SUITES[5], SUITES[6], SUITES[12], SUITES[13]) == (
        "regions", "islands", "discordance", "denovo_precision", "variants")  # fmt: skip


# ------------------------------------------------------------------ suite_result


def test_suite_result_counts_one_passing_suite(mod: ModuleType) -> None:
    w = World.default()
    res = mod.suite_result(record(3), w.junits(3), 3, C)
    assert (res.task, res.suite, res.exit, res.sharded) == (3, "regions", 0, False)
    assert res.problems == []
    plain = [c for c in w.cases[3] if not c.attack]
    att = [c for c in w.cases[3] if c.attack]
    assert res.all.passed == sum(c.outcome == "passed" for c in w.cases[3])
    assert res.all.skipped == sum(c.outcome == "skipped" for c in w.cases[3])
    assert res.attacks.passed == sum(c.outcome == "passed" for c in att)
    assert res.attacks.skipped == sum(c.outcome == "skipped" for c in att)
    total = res.all
    n = total.passed + total.failed + total.errors + total.skipped + total.xfailed + total.xpassed
    assert n == len(plain) + len(att)
    assert (total.failed, total.errors, total.xpassed) == (0, 0, 0)


def test_suite_result_reads_xfail_failure_and_error(mod: ModuleType) -> None:
    w = World.default()
    w.set_outcome(10, attack=False, outcome="failed")
    w.set_outcome(10, attack=True, outcome="error")
    res = mod.suite_result(record(10, exit_=1), w.junits(10), 10, C)
    assert res.all.xfailed == 1  # <skipped type="pytest.xfail"> (task 10 holds one xfail)
    assert (res.all.failed, res.all.errors) == (1, 1)
    assert (res.attacks.failed, res.attacks.errors) == (0, 1)
    assert res.problems


def test_suite_result_sums_the_shard_files_of_task_13(mod: ModuleType) -> None:
    w = World.default()
    res = mod.suite_result(record(SHARDED), w.junits(SHARDED), SHARDED, C)
    assert res.sharded is True and res.suite == "variants" and res.problems == []
    assert res.all.passed + res.all.skipped + res.all.xfailed == len(w.cases[SHARDED])


def test_suite_result_attack_is_by_classname_prefix(mod: ModuleType) -> None:
    w = World.default()
    w.cases[3][0].classname = "tests.calibration.attacks_like.test_x"  # not the package
    res = mod.suite_result(record(3), w.junits(3), 3, C)
    assert res.attacks.passed == sum(1 for c in w.cases[3] if c.attack and c.outcome == "passed")


@pytest.mark.parametrize(
    ("change", "needle"),
    [(lambda r: r.update(commit=OTHER), "commit mismatch"),
     (lambda r: r.update(dirty=True), "dirty"),
     (lambda r: r.update(args=["--fast"]), "args --fast"),
     (lambda r: r.update(exit=2), "exit 2")],
)  # fmt: skip
def test_suite_result_names_the_problem(mod: ModuleType, change, needle: str) -> None:  # type: ignore[no-untyped-def]
    w = World.default()
    rec = record(4)
    change(rec)
    res = mod.suite_result(rec, w.junits(4), 4, C)
    assert any(needle in p for p in res.problems), res.problems


def test_suite_result_without_record(mod: ModuleType) -> None:
    res = mod.suite_result(None, [], 6, C)
    assert res.task == 6
    assert any("no record" in p for p in res.problems), res.problems


def test_suite_result_with_zero_testcases_is_a_problem(mod: ModuleType) -> None:
    w = World.default()
    w.cases[4] = []
    res = mod.suite_result(record(4), w.junits(4), 4, C)
    assert res.problems


# ------------------------------------------------------------------ full_rows: examples 1-8


def test_example1_all_pass(mod: ModuleType) -> None:
    rows = World.default().rows(mod)
    assert rows == [f"{PREFIX}0\t{EX1_FULL}\n", f"{PREFIX}0\t{EX1_ATTACKS}\n"]


def test_example2_non_attack_failure_spares_the_attack_row(mod: ModuleType) -> None:
    w = World.default()
    w.records[3] = record(3, exit_=1)
    w.set_outcome(3, attack=False, outcome="failed")
    got = by_label(w.rows(mod))
    assert got["calibration_full"] == (1, (
        "calibration_full: 15 suites, 1 failed suites (regions exit 1); 422 tests: 1 failed, "
        "0 errors, 399 passed, 20 skipped, 2 xfailed, 0 xpassed; array 1900000"))  # fmt: skip
    assert got["calibration_attacks"] == (0, EX1_ATTACKS)


def test_example3_attack_failure_fails_both_rows(mod: ModuleType) -> None:
    w = World.default()
    w.records[5] = record(5, exit_=1)
    w.set_outcome(5, attack=True, outcome="failed")
    got = by_label(w.rows(mod))
    assert got["calibration_full"] == (1, (
        "calibration_full: 15 suites, 1 failed suites (islands exit 1); 422 tests: 1 failed, "
        "0 errors, 399 passed, 20 skipped, 2 xfailed, 0 xpassed; array 1900000"))  # fmt: skip
    assert got["calibration_attacks"] == (1, (
        "calibration_attacks: 14 suites, 1 failed suites (islands exit 1); 230 attack tests: "
        "1 failed, 0 errors, 219 passed, 10 skipped, 0 xfailed, 0 xpassed; array 1900000"
    ))  # fmt: skip


def test_example4_missing_record_fails_both_rows(mod: ModuleType) -> None:
    w = World.default()
    w.records[6] = None  # discordance is task 6 of replicates_full (the spec's table says 7)
    got = by_label(w.rows(mod))
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert ex == 1 and "1 failed suites (" in s and "no record" in s, s


def test_example5_sharded_merge_exit_fails_closed_for_attacks(mod: ModuleType) -> None:
    w = World.default()
    w.records[SHARDED] = record(SHARDED, exit_=1)
    att = [c for c in w.cases[SHARDED] if c.attack]
    for c in w.cases[SHARDED]:
        c.outcome = "skipped"  # shard tests skip without a verdict
    from helpers.calibration_rows_fixture import junit_xml

    w.junit_override[SHARDED] = {"task-13.shard1.junit.xml": junit_xml(att)}
    got = by_label(w.rows(mod))
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert ex == 1 and "1 failed suites (" in s, s
        assert re.search(r"variants [^,)]*exit 1", s), s


def test_example6_sharded_failure_without_attacks_spares_the_attack_row(mod: ModuleType) -> None:
    w = World.default()
    w.records[SHARDED] = record(SHARDED, exit_=1)
    n_att = sum(c.attack for c in w.cases[SHARDED])
    w.cases[SHARDED] = [c for c in w.cases[SHARDED] if not c.attack]
    got = by_label(w.rows(mod))
    ex, s = got["calibration_full"]
    assert ex == 1 and re.search(r"1 failed suites \(variants [^,)]*exit 1\)", s), s
    ex, s = got["calibration_attacks"]
    assert ex == 0, s
    assert s.startswith("calibration_attacks: 13 suites, 0 failed suites; "
                        f"{230 - n_att} attack tests: 0 failed, 0 errors, "), s  # fmt: skip


def test_example7_fast_arguments_never_pass(mod: ModuleType) -> None:
    w = World.default()
    w.records[0] = record(0, args=["--fast"])
    got = by_label(w.rows(mod))
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert ex == 1 and "1 failed suites (placement args --fast)" in s, s


def test_example8_commit_mismatch_fails_both_rows(mod: ModuleType) -> None:
    w = World.default()
    w.records[0] = record(0, commit=OTHER)
    got = by_label(w.rows(mod))
    for label in ("calibration_full", "calibration_attacks"):
        ex, s = got[label]
        assert ex == 1 and "1 failed suites (placement commit mismatch)" in s, s


# ------------------------------------------------------------------ the other rules


def test_dirty_checkout_fails_both_rows(mod: ModuleType) -> None:
    w = World.default()
    w.records[2] = record(2, dirty=True)
    got = by_label(w.rows(mod))
    assert all(ex == 1 and "dirty" in s for ex, s in got.values()), got


def test_zero_testcases_fails_both_rows(mod: ModuleType) -> None:
    w = World.default()
    w.cases[4] = []
    got = by_label(w.rows(mod))
    assert all(ex == 1 and "1 failed suites (" in s for ex, s in got.values()), got


def test_crash_exit_fails_the_attack_row_of_an_attack_suite(mod: ModuleType) -> None:
    """Unsharded exit not in {0, 1} (crash, interrupt, collection error): JUnit may be partial."""
    w = World.default()
    w.records[8] = record(8, exit_=2)
    got = by_label(w.rows(mod))
    assert got["calibration_full"][0] == 1 and got["calibration_attacks"][0] == 1, got


def test_term_exit_fails_the_attack_row(mod: ModuleType) -> None:
    w = World.default()
    w.records[8] = record(8, exit_=143)
    got = by_label(w.rows(mod))
    assert got["calibration_attacks"][0] == 1 and "exit 143" in got["calibration_attacks"][1]


def test_unattributable_exit_1_fails_the_attack_row(mod: ModuleType) -> None:
    """Unsharded exit 1 with no failed testcase at all: the failure cannot be attributed."""
    w = World.default()
    w.records[8] = record(8, exit_=1)
    got = by_label(w.rows(mod))
    assert got["calibration_full"][0] == 1 and got["calibration_attacks"][0] == 1, got


def test_attack_error_testcase_fails_the_attack_row(mod: ModuleType) -> None:
    w = World.default()
    w.records[8] = record(8, exit_=1)
    w.set_outcome(8, attack=True, outcome="error")
    got = by_label(w.rows(mod))
    ex, s = got["calibration_attacks"]
    assert ex == 1 and re.search(r"attack tests: 0 failed, 1 errors?,", s), s


def test_suite_without_attacks_is_not_counted_in_the_attack_row(mod: ModuleType) -> None:
    w = World.default()
    w.records[NO_ATTACKS] = record(NO_ATTACKS, exit_=1)
    w.set_outcome(NO_ATTACKS, attack=False, outcome="failed")
    got = by_label(w.rows(mod))
    assert got["calibration_full"][0] == 1
    assert got["calibration_attacks"] == (0, EX1_ATTACKS)


def test_missing_junit_part_fails_both_rows(mod: ModuleType) -> None:
    """A JUnit file listed in ``parts`` that is absent: ``junits`` has fewer roots than parts."""
    w = World.default()
    w.junit_override[SHARDED] = {"task-13.shard1.junit.xml": None}
    got = by_label(w.rows(mod))
    assert got["calibration_full"][0] == 1 and got["calibration_attacks"][0] == 1, got


def test_two_failed_suites_are_both_named(mod: ModuleType) -> None:
    w = World.default()
    w.records[0] = record(0, args=["--fast"])
    w.records[6] = record(6, commit=OTHER)
    ex, s = by_label(w.rows(mod))["calibration_full"]
    assert ex == 1 and "2 failed suites (" in s, s
    assert "placement args --fast" in s and "discordance commit mismatch" in s, s


def test_full_rows_are_pure(mod: ModuleType) -> None:
    w = World.default()
    results = w.results(mod)
    assert mod.full_rows(results, JOB, C, DATE) == mod.full_rows(results, JOB, C, DATE)
    other = mod.full_rows(results, "1900500", C, "2026-11-03")
    assert other[0].startswith(f"2026-11-03\t1900500_calibration\t{C}\t0\t"), other
    assert other[0].rstrip("\n").endswith("; array 1900500"), other


# ------------------------------------------------------------------ AC3: property over outcomes

KINDS = ("ok", "missing", "plain_fail", "attack_fail", "attack_error", "crash", "unattributed",
         "dirty", "args", "mismatch", "empty")  # fmt: skip


def _apply(w: World, task: int, kind: str) -> None:
    has_att = any(c.attack for c in w.cases[task])
    if kind == "missing":
        w.records[task] = None
    elif kind == "plain_fail":
        w.records[task] = record(task, exit_=1)
        w.set_outcome(task, attack=False, outcome="failed")
    elif kind in ("attack_fail", "attack_error") and has_att:
        w.records[task] = record(task, exit_=1)
        w.set_outcome(task, attack=True, outcome=kind.split("_")[1].replace("fail", "failed"))
    elif kind == "crash":
        w.records[task] = record(task, exit_=2)
    elif kind == "unattributed":
        w.records[task] = record(task, exit_=1)
    elif kind == "dirty":
        w.records[task] = record(task, dirty=True)
    elif kind == "args":
        w.records[task] = record(task, args=["--fast"])
    elif kind == "mismatch":
        w.records[task] = record(task, commit=OTHER)
    elif kind == "empty":
        w.cases[task] = []


def _expect(task: int, kind: str, has_att: bool) -> tuple[bool, bool]:
    """(full ok, attack ok) of one suite under the spec's rules."""
    if kind == "ok" or (kind in ("attack_fail", "attack_error") and not has_att):
        return True, True
    if kind in ("missing", "dirty", "args", "mismatch", "empty"):
        return False, False  # fail closed
    if not has_att:
        return False, True  # not an attack suite
    if task == SHARDED:
        return False, False  # sharded, exit != 0, holds an attack testcase
    return False, kind == "plain_fail"


@settings(max_examples=60, deadline=None)
@given(st.dictionaries(st.integers(0, 14), st.sampled_from(KINDS), max_size=5))
def test_rows_are_well_formed_and_exit_matches_counts(kinds: dict[int, str]) -> None:
    mod = load_rows_module()
    w = World.default()
    want_full = want_att = True
    for t, kind in kinds.items():
        if kind in ("crash", "unattributed") and t == NO_ATTACKS:
            kind = "plain_fail"  # unspecified for a suite without attacks
        has_att = any(c.attack for c in w.cases[t])
        _apply(w, t, kind)
        f, a = _expect(t, kind, has_att)
        want_full, want_att = want_full and f, want_att and a
    rows = w.rows(mod)
    assert len(rows) == 2
    labels = []
    for line in rows:
        date, job, commit, ex, summary = split_row(line)
        assert (date, job) == (DATE, f"{JOB}_calibration")
        assert re.fullmatch(r"[0-9a-f]{40}", commit) and commit == C
        label = summary.split(":", 1)[0]
        labels.append(label)
        failed = re.search(r", (\d+) failed suites\b", summary)
        assert failed, summary
        zero = int(failed.group(1)) == 0 and all(int(n) == 0 for n in COUNT.findall(summary))
        assert (ex == 0) == zero, summary
        assert ex in (0, 1)
        if ex == 1:
            assert int(failed.group(1)) >= 1, summary
        assert summary.endswith(f"; array {JOB}"), summary
        assert (ex == 0) == (want_full if label == "calibration_full" else want_att), (kinds, line)
    assert labels == ["calibration_full", "calibration_attacks"]
