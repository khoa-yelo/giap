"""Rows of ``calibration_rows.full_rows`` read back by the merged ``readiness.py`` (issue M11-03b,
acceptance criterion 1; lead decision Q6).

Examples 1-8 of docs/specs/M11-03b-nightly-readiness-labels.md: the two rows, written for the
fixture release commit R (``helpers.release_readiness.Fixture``) next to a passing
``calibration_nightly`` row, give exactly the R3/R4 statuses of the spec's table through
``parse_log`` and ``_calib`` (``r3``/``r4``). Q6: ``tests/sim`` and ``tests/helpers`` are inputs of
the calibration suites, so a change there makes an older calibration row stale (pending).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from types import ModuleType

import pytest
from helpers.calibration_rows_fixture import (
    JOB,
    OTHER,
    SHARDED,
    C,
    World,
    junit_xml,
    load_rows_module,
    record,
)
from helpers.release_readiness import VERSION, Fixture, load_module, nightly_row

NIGHTLY_OK = "calibration_nightly: 17290 passed, 108 skipped, 13 xfailed in 20556.28s"


def ex1(w: World) -> None:
    pass


def ex2(w: World) -> None:
    w.records[3] = record(3, exit_=1)
    w.set_outcome(3, attack=False, outcome="failed")


def ex3(w: World) -> None:
    w.records[5] = record(5, exit_=1)
    w.set_outcome(5, attack=True, outcome="failed")


def ex4(w: World) -> None:
    w.records[6] = None  # discordance


def ex5(w: World) -> None:
    w.records[SHARDED] = record(SHARDED, exit_=1)
    for c in w.cases[SHARDED]:
        c.outcome = "skipped"
    att = [c for c in w.cases[SHARDED] if c.attack]
    w.junit_override[SHARDED] = {"task-13.shard1.junit.xml": junit_xml(att)}


def ex6(w: World) -> None:
    w.records[SHARDED] = record(SHARDED, exit_=1)
    w.cases[SHARDED] = [c for c in w.cases[SHARDED] if not c.attack]


def ex7(w: World) -> None:
    w.records[0] = record(0, args=["--fast"])


def ex8(w: World) -> None:
    w.records[0] = record(0, commit=OTHER)


EXAMPLES: list[tuple[str, Callable[[World], None], str, str]] = [
    ("1", ex1, "pass", "pass"),
    ("2", ex2, "fail", "pass"),
    ("3", ex3, "fail", "fail"),
    ("4", ex4, "fail", "fail"),
    ("5", ex5, "fail", "fail"),
    ("6", ex6, "fail", "pass"),
    ("7", ex7, "fail", "fail"),
    ("8", ex8, "fail", "fail"),
]


@pytest.fixture(scope="module")
def rows_mod() -> ModuleType:
    return load_rows_module()


def _rows_for(rows_mod: ModuleType, w: World, commit: str, date: str) -> list[str]:
    return list(rows_mod.full_rows(w.results(rows_mod, commit), JOB, commit, date))


def _pin(w: World, commit: str) -> None:
    """Records of the fixture commit ``C`` become records of the release commit."""
    for rec in w.records.values():
        if rec is not None and rec["commit"] == C:
            rec["commit"] = commit


def _statuses(fx: Fixture) -> tuple[str, str]:
    rd = load_module()
    rd.code_equal.cache_clear()  # functools.cache: HEAD moves between calls in one test
    c = rd.Ctx(repo=fx.repo, version=VERSION, logs=fx.logs)
    return rd.r3(c)["status"], rd.r4(c)["status"]


@pytest.mark.parametrize(("ex", "change", "r3", "r4"), EXAMPLES, ids=[e[0] for e in EXAMPLES])
def test_rows_round_trip_through_readiness(
    tmp_path: Path, rows_mod: ModuleType, ex: str, change: Callable[[World], None], r3: str,
    r4: str,
) -> None:  # fmt: skip
    fx = Fixture.build(tmp_path)
    w = World.default()
    change(w)
    _pin(w, fx.r)
    rows = _rows_for(rows_mod, w, fx.r, "2026-11-02")
    fx.set_nightly(nightly_row("2026-11-02", "2000_2", fx.short(fx.r), 0, NIGHTLY_OK)
                   + "".join(rows))  # fmt: skip
    parsed, bad = load_module().parse_log((fx.logs / "nightly.tsv").read_text())
    assert bad == 0, rows
    assert [r.label for r in parsed] == ["calibration_nightly", "calibration_full",
                                         "calibration_attacks"]  # fmt: skip
    assert _statuses(fx) == (r3, r4), (ex, rows)


def test_a_later_failing_rerun_beats_an_earlier_pass(tmp_path: Path, rows_mod: ModuleType) -> None:
    """Edge case "Re-run": same date, the larger array ID decides (newest)."""
    fx = Fixture.build(tmp_path)
    good = World.default()
    bad = good.copy()
    ex3(bad)
    _pin(good, fx.r)
    _pin(bad, fx.r)
    first = list(rows_mod.full_rows(good.results(rows_mod, fx.r), "1900000", fx.r, "2026-11-02"))
    second = list(rows_mod.full_rows(bad.results(rows_mod, fx.r), "1900500", fx.r, "2026-11-02"))
    nightly = nightly_row("2026-11-02", "2000_2", fx.short(fx.r), 0, NIGHTLY_OK)
    fx.set_nightly(nightly + "".join(second) + "".join(first))
    assert _statuses(fx) == ("fail", "fail")
    fx.set_nightly(nightly + "".join(first + second).replace("1900500", "1899999"))
    assert _statuses(fx) == ("pass", "pass")


# ------------------------------------------------------------------ Q6: CALIB_PATHS


def test_calib_paths_cover_the_simulator_and_the_test_helpers() -> None:
    paths = load_module().CALIB_PATHS
    assert {"src", "tests/calibration", "ci", "tests/sim", "tests/helpers"} <= set(paths), paths


@pytest.mark.parametrize("rel", ["tests/sim/giapsim/nulls.py", "tests/helpers/region_oracle.py"])
def test_a_change_in_a_calibration_input_makes_the_rows_stale(tmp_path: Path, rel: str) -> None:
    fx = Fixture.build(tmp_path)  # its nightly rows (all labels) are for R
    assert _statuses(fx) == ("pass", "pass")
    fx.write_text(rel, '"""changed null generator (fixture)."""\n')
    fx.commit(f"change {rel}")
    r3, r4 = _statuses(fx)
    assert (r3, r4) == ("pending", "pending"), (rel, r3, r4)


def test_a_docs_change_keeps_the_rows(tmp_path: Path) -> None:
    fx = Fixture.build(tmp_path)
    fx.write_text("docs/user-guide/install.md", "# Install (edited)\n")
    fx.commit("docs")
    assert _statuses(fx) == ("pass", "pass")
