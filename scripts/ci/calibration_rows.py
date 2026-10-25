"""``calibration_full`` and ``calibration_attacks`` rows of ``nightly.tsv`` (issue M11-03b).

Reads the task records (``task-<i>.json``) and JUnit files that ``calibration_full.sbatch`` leaves
in one log directory and appends two rows for ``scripts/release/readiness.py`` (R3, R4):

    <date>\\t<array>_calibration\\t<40-hex commit>\\t<exit>\\tcalibration_full: …
    <date>\\t<array>_calibration\\t<40-hex commit>\\t<exit>\\tcalibration_attacks: …

``suite_result`` and ``full_rows`` are pure (no Slurm, git, clock or file I/O); the command line
part reads the log directory, takes ``flock`` on ``<tsv>.lock``, appends both rows in one write and
leaves ``<log-dir>/rows.done`` so a second run on the same directory writes nothing.
Spec: docs/specs/M11-03b-nightly-readiness-labels.md. Exit 0 rows written and both pass, 1 rows
written and one fails, 2 usage error or unreadable log directory (failing rows are still written
for the latter).
"""

from __future__ import annotations

import argparse
import contextlib
import datetime
import fcntl
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"
ATTACK_PREFIX = "tests.calibration.attacks."
REPLICATES_FULL = 15  # suites of tolerances.yaml replicates_full = calibration_full.sbatch tasks
# suites without any attack testcase (data, #349 B2); any other suite without one fails R4
ATTACK_FREE = frozenset({"denovo_precision"})
HEX40 = re.compile(r"[0-9a-f]{40}")
_COUNTLIKE = re.compile(r"(\d) (failed|error)")  # keep reasons from reading as test counts


@dataclass(frozen=True)
class Counts:
    passed: int = 0
    failed: int = 0
    errors: int = 0
    skipped: int = 0
    xfailed: int = 0
    xpassed: int = 0

    def __add__(self, other: Counts) -> Counts:
        return Counts(*(a + b for a, b in zip(self.astuple(), other.astuple(), strict=True)))

    def astuple(self) -> tuple[int, int, int, int, int, int]:
        return (self.passed, self.failed, self.errors, self.skipped, self.xfailed, self.xpassed)

    @property
    def total(self) -> int:
        return sum(self.astuple())

    @property
    def bad(self) -> int:
        return self.failed + self.errors


@dataclass(frozen=True)
class SuiteResult:
    task: int
    suite: str
    exit: int | None
    sharded: bool
    all: Counts
    attacks: Counts
    problems: list[str] = field(default_factory=list)
    fail_closed: bool = False  # the outcome is unknown (record, JUnit, commit, args, crash, …)
    attacks_failed: bool = False
    attack_problems: list[str] = field(default_factory=list)  # R4-only: no or only-skipped attacks

    @property
    def full_failed(self) -> bool:
        return bool(self.problems)

    @property
    def attack_suite(self) -> bool:
        """Counted in the attack row: holds an attack testcase, should hold one (not in
        ``ATTACK_FREE``), or its outcome is unknown."""
        return self.fail_closed or self.attacks.total > 0 or bool(self.attack_problems)


def _outcome(case: ET.Element) -> str:
    """pytest's JUnit: <failure> (also strict XPASS), <error>, <skipped type="pytest.xfail">."""
    tags = {child.tag: child for child in case}
    if "failure" in tags:
        return "failed"
    if "error" in tags:
        return "errors"
    if "skipped" in tags:
        return "xfailed" if tags["skipped"].get("type") == "pytest.xfail" else "skipped"
    return "passed"


# a failure is never hidden, a verdict beats a no-verdict skip (SHARD_SKIP); ties: the later part
_RANK = {"skipped": 0, "passed": 1, "xfailed": 1, "xpassed": 1, "errors": 2, "failed": 3}


def _is_attack(classname: str, name: str) -> bool:
    """An attack testcase; a collection error has ``classname=""`` and the module as ``name``."""
    return classname.startswith(ATTACK_PREFIX) or (not classname and name.startswith(ATTACK_PREFIX))


def _testcases(junits: list[ET.Element]) -> dict[tuple[str, str], str]:
    """Outcome per (classname, name) over every part (#349 S1: task 13 shards and merge write the
    same testcases), keeping the worst-ranked one."""
    out: dict[tuple[str, str], str] = {}
    for root in junits:
        for case in root.iter("testcase"):
            key, outcome = (case.get("classname", ""), case.get("name", "")), _outcome(case)
            if key not in out or _RANK[outcome] >= _RANK[out[key]]:
                out[key] = outcome
    return out


def _count(outcomes: list[str]) -> Counts:
    n = dict.fromkeys(("passed", "failed", "errors", "skipped", "xfailed", "xpassed"), 0)
    for outcome in outcomes:
        n[outcome] += 1
    return Counts(**n)


def _safe_part(part: object) -> bool:
    """A plain file name in the log directory (no path, no ``.``/``..``)."""
    return isinstance(part, str) and part not in ("", ".", "..") and Path(part).name == part


def _record_problem(record: dict[str, Any], task: int, expect_commit: str) -> list[str]:
    """Fail-closed problems of a parsed record (types first: a malformed record says nothing)."""
    types = {"task": int, "exit": int, "commit": str, "dirty": bool, "args": list, "parts": list}
    if any(not isinstance(record.get(k), t) or isinstance(record.get(k), bool) is not (t is bool)
           for k, t in types.items()) or not all(map(_safe_part, record["parts"])):  # fmt: skip
        return ["record malformed"]
    if record["task"] != task:
        return ["task mismatch"]
    out = []
    if record["commit"] != expect_commit:
        out.append("commit mismatch")
    if record["dirty"]:
        out.append("dirty")
    if record["args"]:
        out.append("args " + " ".join(str(a) for a in record["args"]))
    return out


def suite_result(
    record: dict[str, Any] | None,
    junits: list[ET.Element],
    task: int,
    expect_commit: str,
    *,
    suite: str | None = None,
    missing: str = "no record",
) -> SuiteResult:
    """One array task's outcome. ``junits`` are the parsed files of ``record["parts"]`` that could
    be read; ``suite`` names the task when the record cannot (``replicates_full`` key); ``missing``
    is the problem of an absent record (e.g. ``"no record, OUT_OF_MEMORY"``).

    Fail closed (both rows): record problems, a missing JUnit part, no testcases, and a crash
    (exit not in {0, 1}, sharded or not: the JUnit may be partial, #349 B1). R4 only: an attack
    testcase failed or errored; exit 1 that no testcase explains, or exit 1 of the sharded task
    with attack testcases; a suite outside ``ATTACK_FREE`` without attack testcases (except the
    sharded task's exit 1, example 6) or with only skipped ones (#349 B2)."""
    name = suite or (record or {}).get("suite") or f"task {task}"
    if not isinstance(name, str):
        name = f"task {task}"
    if record is None:
        return SuiteResult(task, name, None, False, Counts(), Counts(), [missing], True, True)
    closed = _record_problem(record, task, expect_commit)
    if suite is not None and record.get("suite") != suite and not closed:
        closed.append("suite mismatch")
    rc = record.get("exit") if isinstance(record.get("exit"), int) else None
    sharded = record.get("sharded") is True
    raw_parts = record.get("parts")
    parts: list[Any] = raw_parts if isinstance(raw_parts, list) else []
    cases = _testcases(junits)
    every = _count(list(cases.values()))
    att = _count([o for (cls, n), o in cases.items() if _is_attack(cls, n)])
    if "record malformed" not in closed:
        if len(junits) < len(parts):
            closed.append("junit missing")
        elif not cases:
            closed.append("no testcases")
        if rc not in (0, 1):
            closed.append(f"exit {rc}")  # crash, interrupt, collection or merge error (B1)
    problems = list(closed)
    if not closed and rc != 0:
        problems.append(f"exit {rc}")
    elif not closed and every.bad:
        problems.append("failed testcases")
    extra = []
    if not closed and not att.total and name not in ATTACK_FREE and not (sharded and rc == 1):
        extra.append("no attack tests")
    elif not closed and att.total and att.skipped == att.total:
        extra.append("attack tests all skipped")
    attacks_failed = bool(closed) or att.bad > 0 or bool(extra)
    if not closed and att.total and rc == 1 and (sharded or not every.bad):
        attacks_failed = True  # unattributable exit 1 or merge exit 1
    return SuiteResult(task, name, rc, sharded, every, att, problems, bool(closed), attacks_failed,
                       extra)  # fmt: skip


def _clean(text: str) -> str:
    return _COUNTLIKE.sub(r"\1_\2", " ".join(text.split()))


def _summary(label: str, unit: str, results: list[SuiteResult], failed: list[str],
             counts: Counts, job: str) -> str:  # fmt: skip
    named = ", ".join(_clean(f) for f in failed)
    c = counts
    return (f"{label}: {len(results)} suites, {len(failed)} failed suites"
            + (f" ({named})" if named else "")
            + f"; {c.total} {unit}: {c.failed} failed, {c.errors} errors, {c.passed} passed, "
            f"{c.skipped} skipped, {c.xfailed} xfailed, {c.xpassed} xpassed; array {job}"
            )  # fmt: skip


def full_rows(results: list[SuiteResult], job: str, commit: str, date: str, *,
              expect: int = REPLICATES_FULL) -> list[str]:  # fmt: skip
    """The ``calibration_full`` and ``calibration_attacks`` rows (each ends with one newline).

    Both rows fail unless ``expect`` suites are present; the attack row also fails with 0 attack
    suites or 0 attack tests (#349 B2). Such a row-level reason counts as one failed suite, so a
    failing row always has ``failed suites`` >= 1."""
    att = [r for r in results if r.attack_suite]
    whole = [] if len(results) == expect else [f"expected {expect} suites, got {len(results)}"]
    full_counts = sum((r.all for r in results), Counts())
    att_counts = sum((r.attacks for r in att), Counts())
    vacuous = ["no attack tests in the row"] if not att or not att_counts.total else []
    tiers = [
        ("calibration_full", "tests", results, full_counts,
         [f"{r.suite} {', '.join(r.problems)}" for r in results if r.full_failed] + whole),
        ("calibration_attacks", "attack tests", att, att_counts,
         [f"{r.suite} {', '.join(r.problems + r.attack_problems) or 'attack tests'}"
          for r in att if r.attacks_failed] + whole + vacuous),
    ]  # fmt: skip
    rows = []
    for label, unit, members, counts, failed in tiers:
        ok = not failed and counts.bad == 0
        if not ok and not failed:  # cannot happen (a bad count fails its suite); fail visibly
            failed = [f"{r.suite} failed testcases" for r in members if r.all.bad or r.attacks.bad]
        summary = _summary(label, unit, members, failed, counts, job)
        rows.append(f"{date}\t{job}_calibration\t{commit}\t{0 if ok else 1}\t{summary}\n")
    return rows


# ------------------------------------------------------------------ command line (I/O)


def _suites() -> list[str]:
    import yaml

    with TOLERANCES.open() as fh:
        return [str(k) for k in yaml.safe_load(fh)["replicates_full"]]


def _read_task(log: Path, task: int, commit: str, suite: str, sacct: dict[str, str]
               ) -> SuiteResult:  # fmt: skip
    path = log / f"task-{task}.json"
    if not path.is_file():
        state = sacct.get(str(task))
        missing = f"no record, {state}" if state else "no record"
        return suite_result(None, [], task, commit, suite=suite, missing=missing)
    try:
        record = json.loads(path.read_text())
    except (OSError, ValueError):
        record = None
    if not isinstance(record, dict):
        return suite_result(None, [], task, commit, suite=suite, missing="record unreadable")
    junits = []
    parts = record.get("parts")
    for part in parts if isinstance(parts, list) else []:
        if not _safe_part(part):
            return suite_result(None, [], task, commit, suite=suite, missing="record malformed")
        with contextlib.suppress(OSError, ET.ParseError):  # fewer roots than parts: "junit missing"
            junits.append(ET.parse(log / part).getroot())
    return suite_result(record, junits, task, commit, suite=suite)


def _append(tsv: Path, rows: list[str], done: Path | None) -> bool:
    """Append ``rows`` in one write under ``flock`` on ``<tsv>.lock``; False if already done."""
    with open(f"{tsv}.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if done is not None and done.exists():
            return False
        fd = os.open(tsv, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
        try:
            os.write(fd, "".join(rows).encode())
        finally:
            os.close(fd)
        if done is not None:
            done.write_text(f"{rows[0].split(chr(9))[1]}\n")
    return True


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--log-dir", type=Path, required=True)
    ap.add_argument("--job", required=True, help="array job ID of calibration_full.sbatch")
    ap.add_argument("--commit", required=True, help="40-hex commit every record must carry")
    ap.add_argument("--expect", type=int, required=True, help="number of replicates_full suites")
    ap.add_argument("--tsv", type=Path, required=True)
    ap.add_argument("--sacct-json", type=Path, help='{"<task>": "<sacct state>", …}')
    ns = ap.parse_args(argv)
    if not HEX40.fullmatch(ns.commit) or not re.fullmatch(r"[0-9]+", ns.job):
        print("calibration_rows: --commit must be 40-hex and --job an array ID", file=sys.stderr)
        return 2
    try:
        suites = _suites()
        sacct = json.loads(ns.sacct_json.read_text()) if ns.sacct_json else {}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"calibration_rows: {exc}", file=sys.stderr)
        return 2
    if ns.expect != len(suites) or not isinstance(sacct, dict):
        print(f"calibration_rows: --expect {ns.expect} but tolerances.yaml has {len(suites)} "
              "replicates_full suites (or --sacct-json is not an object)",
              file=sys.stderr)  # fmt: skip
        return 2
    sacct = {str(k): _clean(str(v)) for k, v in sacct.items()}
    readable = ns.log_dir.is_dir() and os.access(ns.log_dir, os.R_OK | os.X_OK)
    if readable:
        results = [_read_task(ns.log_dir, t, ns.commit, s, sacct) for t, s in enumerate(suites)]
    else:
        results = [suite_result(None, [], t, ns.commit, suite=s, missing="no log dir")
                   for t, s in enumerate(suites)]  # fmt: skip
    rows = full_rows(results, ns.job, ns.commit, datetime.date.today().isoformat(),
                     expect=ns.expect)  # fmt: skip
    try:
        wrote = _append(ns.tsv, rows, ns.log_dir / "rows.done" if readable else None)
    except OSError as exc:
        print(f"calibration_rows: cannot append to {ns.tsv}: {exc}", file=sys.stderr)
        return 2
    if not wrote:
        print(f"calibration_rows: {ns.log_dir}/rows.done exists; nothing written")
        return 0
    sys.stdout.write("".join(rows))
    if not readable:
        print(f"calibration_rows: log directory {ns.log_dir} unreadable", file=sys.stderr)
        return 2
    return 0 if all(r.split("\t")[3] == "0" for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
