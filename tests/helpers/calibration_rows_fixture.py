"""Fixture world for ``scripts/ci/calibration_rows.py`` (issue M11-03b, examples 1-9).

``World.default()`` is example 1 of docs/specs/M11-03b-nightly-readiness-labels.md: 15 task records
(``task-<i>.json``, exit 0, commit ``C``) and hand-written JUnit XML with 422 testcases (400
passed, 20 skipped, 2 xfailed), 230 of them attack testcases (classname
``tests.calibration.attacks.…``; 220 passed, 10 skipped) spread over the 14 suites other than
``denovo_precision`` (task 12). Task 13 (``variants``) is sharded: its testcases sit in the
``shard0`` JUnit, ``shard1`` and ``merge`` hold an empty ``<testsuite>``.

A test changes one task (``set_outcome``, ``records[i] = …``) and either calls the pure
``suite_result``/``full_rows`` on parsed inputs (``World.results``) or writes the log directory
(``World.write``) and runs the CLI.
"""

from __future__ import annotations

import copy
import functools
import importlib.util
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any
from xml.sax.saxutils import quoteattr

import yaml

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "ci" / "calibration_rows.py"
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"
SUITES: tuple[str, ...] = tuple(yaml.safe_load(TOLERANCES.read_text())["replicates_full"])
C = "c0ffee00" * 5  # 40-hex commit of the run
OTHER = "deadbeef" * 5
JOB = "1900000"
DATE = "2026-11-02"
SHARDED = 13  # variants
NO_ATTACKS = 12  # denovo_precision
PREFIX = f"{DATE}\t{JOB}_calibration\t{C}\t"

OUTCOMES = ("passed", "failed", "error", "skipped", "xfailed")


@functools.cache
def load_rows_module() -> ModuleType:
    """``scripts/ci/calibration_rows.py`` as a module (its CLI sits behind ``__main__``)."""
    assert SCRIPT.is_file(), f"{SCRIPT} does not exist (M11-03b row writer)"
    spec = importlib.util.spec_from_file_location("giap_calibration_rows", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@dataclass
class Case:
    classname: str
    name: str
    outcome: str = "passed"

    @property
    def attack(self) -> bool:
        return self.classname.startswith("tests.calibration.attacks.")


def attack_class(task: int) -> str:
    return f"tests.calibration.attacks.test_attack_{SUITES[task]}"


def plain_class(task: int) -> str:
    return f"tests.calibration.test_{SUITES[task]}_calibration"


def default_cases() -> dict[int, list[Case]]:
    cases: dict[int, list[Case]] = {t: [] for t in range(len(SUITES))}
    plain = ["passed"] * 180 + ["skipped"] * 10 + ["xfailed"] * 2
    for k, o in enumerate(plain):
        t = k % 15
        cases[t].append(Case(plain_class(t), f"test_p{k}", o))
    attack_tasks = [t for t in range(15) if t != NO_ATTACKS]
    for k, o in enumerate(["passed"] * 220 + ["skipped"] * 10):
        t = attack_tasks[k % 14]
        cases[t].append(Case(attack_class(t), f"test_a{k}", o))
    return cases


def parts_of(task: int) -> list[str]:
    if task == SHARDED:
        return [f"task-{task}.shard0.junit.xml", f"task-{task}.shard1.junit.xml",
                f"task-{task}.merge.junit.xml"]  # fmt: skip
    return [f"task-{task}.run.junit.xml"]


def record(task: int, *, exit_: int = 0, commit: str = C, dirty: bool = False,
           args: list[str] | None = None) -> dict[str, Any]:  # fmt: skip
    return {"task": task, "suite": SUITES[task], "job": JOB, "commit": commit, "dirty": dirty,
            "exit": exit_, "args": list(args or []), "sharded": task == SHARDED,
            "parts": parts_of(task), "date": DATE}  # fmt: skip


def junit_xml(cases: list[Case]) -> str:
    body = []
    for c in cases:
        attrs = f'classname={quoteattr(c.classname)} name={quoteattr(c.name)} time="0.010"'
        inner = {
            "passed": "",
            "failed": '<failure message="assert 0.07 &lt;= 0.05">AssertionError</failure>',
            "error": '<error message="failed on setup">fixture error</error>',
            "skipped": '<skipped type="pytest.skip" message="not in this shard">skip</skipped>',
            "xfailed": '<skipped type="pytest.xfail" message="known power gap" />',
        }[c.outcome]
        body.append(f"<testcase {attrs}>{inner}</testcase>" if inner else f"<testcase {attrs} />")
    n = len(cases)
    f = sum(c.outcome == "failed" for c in cases)
    e = sum(c.outcome == "error" for c in cases)
    s = sum(c.outcome in ("skipped", "xfailed") for c in cases)
    return ('<?xml version="1.0" encoding="utf-8"?><testsuites name="pytest tests">'
            f'<testsuite name="pytest" errors="{e}" failures="{f}" skipped="{s}" tests="{n}" '
            f'time="1.0" timestamp="2026-11-02T03:00:00" hostname="barley-01">'
            + "".join(body) + "</testsuite></testsuites>")  # fmt: skip


@dataclass
class World:
    records: dict[int, dict[str, Any] | None] = field(default_factory=dict)
    cases: dict[int, list[Case]] = field(default_factory=dict)
    # per-task override of the JUnit texts by part name (None = file missing)
    junit_override: dict[int, dict[str, str | None]] = field(default_factory=dict)

    @classmethod
    def default(cls) -> World:
        return cls({t: record(t) for t in range(15)}, default_cases())

    def copy(self) -> World:
        return copy.deepcopy(self)

    def set_outcome(self, task: int, *, attack: bool, outcome: str) -> Case:
        """Turn the first passed (attack or non-attack) testcase of ``task`` into ``outcome``."""
        for c in self.cases[task]:
            if c.attack == attack and c.outcome == "passed":
                c.outcome = outcome
                return c
        raise AssertionError(f"task {task} has no passed {'attack' if attack else 'plain'} case")

    def junit_texts(self, task: int) -> dict[str, str | None]:
        parts = parts_of(task)
        texts: dict[str, str | None] = {p: junit_xml([]) for p in parts}
        texts[parts[0]] = junit_xml(self.cases[task])
        texts.update(self.junit_override.get(task, {}))
        return texts

    def junits(self, task: int) -> list[ET.Element]:
        return [ET.fromstring(t) for t in self.junit_texts(task).values() if t is not None]

    def results(self, mod: ModuleType, expect_commit: str = C) -> list[Any]:
        out = []
        for t in range(15):
            rec = self.records.get(t)
            junits = self.junits(t) if rec is not None else []
            out.append(mod.suite_result(copy.deepcopy(rec), junits, t, expect_commit))
        return out

    def rows(self, mod: ModuleType) -> list[str]:
        return list(mod.full_rows(self.results(mod), JOB, C, DATE))

    def write(self, log_dir: Path) -> Path:
        log_dir.mkdir(parents=True, exist_ok=True)
        for t in range(15):
            rec = self.records.get(t)
            if rec is None:
                continue
            (log_dir / f"task-{t}.json").write_text(json.dumps(rec))
            for name, text in self.junit_texts(t).items():
                if text is not None:
                    (log_dir / name).write_text(text)
        return log_dir


def split_row(line: str) -> tuple[str, str, str, int, str]:
    assert line.endswith("\n") and line.count("\n") == 1, repr(line)
    f = line.rstrip("\n").split("\t")
    assert len(f) == 5, f
    return f[0], f[1], f[2], int(f[3]), f[4]


def by_label(rows: list[str]) -> dict[str, tuple[int, str]]:
    out = {}
    for line in rows:
        _d, _j, _c, ex, summary = split_row(line)
        out[summary.split(":", 1)[0]] = (ex, summary)
    assert set(out) == {"calibration_full", "calibration_attacks"}, rows
    return out


def run_cli(*args: str, timeout: float = 60) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True,
                          timeout=timeout, check=False)  # fmt: skip
