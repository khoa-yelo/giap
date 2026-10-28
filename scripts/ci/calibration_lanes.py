"""Per-suite lane jobs of the full calibration tier (issue M11-03c).

Reads ``scripts/slurm/calibration_plan.tsv`` (one row per ``replicates_full`` suite of
``tests/calibration/tolerances.yaml``, in that order) and turns a row into the lanes one Slurm job
runs side by side: ``shards`` replicate-shard lanes and ``groups`` group lanes.

Pure functions apart from the readers (``read_plan``, ``read_durations``, ``durations_from_junit``,
``durations_from_probe_log``, ``probe_log_lower_bounds``). Spec:
docs/specs/M11-03c-calibration-full-sharding.md.

CLI (exit 2 and a message on stderr on a ``PlanError``):
- ``check [--plan PATH]``: exit 0 when the plan is valid;
- ``lanes --suite S --report PATH [--durations PATH] [--plan PATH]``: the suite's lanes as JSON
  (``report`` is the conftest lane report ``{nodeid: {"shard0_only": bool}}``);
- ``durations --source ID --out PATH [--append | --max] [--probe-log] FILE...``: write the
  durations TSV from JUnit files (the parts of one run; a nodeid in several parts gets the sum) or,
  with ``--probe-log``, from timestamped ``pytest -v`` xdist logs. ``--append`` keeps the rows
  already in ``--out`` and adds only nodeids it lacks; ``--max`` also keeps them, but a nodeid in
  both gets the larger seconds with that row's source (weights that hold on the slowest node class
  measured so far). A probe-log test that started but never finished gets a lower bound (the log's
  last stamp minus its start) with source ``ID:lower-bound``, so the planner weighs it at least
  that much instead of ``unknown_min``. Only nodeids under ``tests/calibration/`` are kept.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
PLAN = REPO / "scripts" / "slurm" / "calibration_plan.tsv"
DURATIONS = REPO / "scripts" / "slurm" / "calibration_durations.tsv"
DURATIONS_HEADER = ("nodeid", "seconds", "source")
CALIBRATION = "tests/calibration/"
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"
HEADER = ("suite", "shards", "groups", "group_select", "cpus_per_lane", "mem_per_lane_mb", "time",
          "unknown_min")  # fmt: skip
MAX_TIME_S = 8 * 3600
#: the partition's MaxMemPerCPU (MB): Slurm raises a job's CPUs to ceil(mem / this)
MAX_MEM_PER_CPU_MB = 4000
LOWER_BOUND = ":lower-bound"  # source suffix of a duration that is only a lower bound
_UINT = re.compile(r"[0-9]+")
_TIME = re.compile(r"([0-9]{2}):([0-5][0-9]):([0-5][0-9])")


class PlanError(ValueError):
    """The plan (or a row against its collection) breaks a rule of "Plan validation"."""


@dataclass(frozen=True)
class PlanRow:
    """One plan row; ``time`` stays the ``HH:MM:SS`` string handed to ``sbatch --time``; ``task``
    is the row index (the task index ``i`` of ``task-<i>.*``), set by ``read_plan``."""

    suite: str
    shards: int
    groups: int
    group_select: str
    cpus_per_lane: int
    mem_per_lane_mb: int
    time: str
    unknown_min: float
    task: int = -1

    @property
    def time_s(self) -> int:
        h, m, s = (int(x) for x in self.time.split(":"))
        return h * 3600 + m * 60 + s

    @property
    def requested_cpus(self) -> int:
        """``--cpus-per-task`` of the suite job: lanes x cpus_per_lane."""
        return (self.shards + self.groups) * self.cpus_per_lane


def allocated_cpus(row: PlanRow, max_mem_per_cpu_mb: int = MAX_MEM_PER_CPU_MB) -> int:
    """The CPUs Slurm allocates (and squeue shows) for ``row``'s job: the requested CPUs, raised
    to ``ceil(lanes x mem_per_lane_mb / MaxMemPerCPU)`` when the memory needs more, rounded up to
    whole cores of 2 threads (CR_CORE_MEMORY; #445 S1). R14 sums this."""
    mem = (row.shards + row.groups) * row.mem_per_lane_mb
    cpus = max(row.requested_cpus, -(-mem // max_mem_per_cpu_mb))
    return -(-cpus // 2) * 2


@dataclass(frozen=True)
class Lane:
    """One lane: ``part`` is its JUnit file name (``None`` for an empty group lane) and
    ``expected`` the nodeids it must run."""

    kind: str  # "shard" or "group"
    index: int
    n: int
    deselect: tuple[str, ...]
    part: str | None
    expected: tuple[str, ...]


def _uint(cells: Mapping[str, str], col: str, where: str) -> int:
    value = cells[col]
    if not _UINT.fullmatch(value):
        raise PlanError(f"{where}: column {col} must be a non-negative integer, got {value!r}")
    return int(value)


def _row(cells: Mapping[str, str], where: str) -> PlanRow:
    shards = _uint(cells, "shards", where)
    groups = _uint(cells, "groups", where)
    if shards + groups == 0:
        raise PlanError(f"{where}: shards + groups must be >= 1")
    sel = cells["group_select"]
    if sel == "*" and shards > 0:
        raise PlanError(f"{where}: group_select '*' needs shards == 0")
    if sel == "-" and groups > 0:
        raise PlanError(f"{where}: group_select '-' needs groups == 0")
    if sel not in ("*", "-"):
        tokens = [t.strip() for t in sel.split(",")]
        if "direct" in tokens and shards == 0:
            raise PlanError(f"{where}: group_select 'direct' needs shards > 0 (use '*')")
        if groups == 0:
            raise PlanError(f"{where}: group_select {sel!r} needs groups > 0")
        if any(not t for t in tokens):
            raise PlanError(f"{where}: group_select {sel!r} has an empty entry")
    cpus = _uint(cells, "cpus_per_lane", where)
    if cpus == 0 or cpus % 2:
        raise PlanError(f"{where}: cpus_per_lane must be even and positive, got {cpus}")
    mem = _uint(cells, "mem_per_lane_mb", where)
    if mem == 0:
        raise PlanError(f"{where}: mem_per_lane_mb must be positive")
    time = cells["time"]
    if not _TIME.fullmatch(time):
        raise PlanError(f"{where}: time must be HH:MM:SS, got {time!r}")
    row = PlanRow(cells["suite"], shards, groups, sel, cpus, mem, time, 0.0)
    if not 0 < row.time_s <= MAX_TIME_S:  # Slurm reads --time=0 as no limit
        raise PlanError(f"{where}: time must be in (00:00:00, 08:00:00], got {time}")
    try:
        unknown = float(cells["unknown_min"])
    except ValueError:
        unknown = -1.0
    if not (math.isfinite(unknown) and unknown > 0):
        raise PlanError(f"{where}: unknown_min must be a positive number")
    return PlanRow(row.suite, shards, groups, sel, cpus, mem, time, unknown)


def read_plan(path: Path, suites: Sequence[str]) -> list[PlanRow]:
    """Parse and validate the plan; its suites must be ``suites``, in order."""
    lines = Path(path).read_text().splitlines()
    if not lines or tuple(lines[0].split("\t")) != HEADER:
        raise PlanError(f"{path}: header must be {' '.join(HEADER)}")
    rows: list[PlanRow] = []
    for n, line in enumerate(lines[1:]):
        cells = line.split("\t")
        where = f"row {n} ({cells[0]})"
        if len(cells) != len(HEADER):
            raise PlanError(f"{where}: {len(cells)} cells, header has {len(HEADER)}")
        row = _row(dict(zip(HEADER, cells, strict=True)), where)
        rows.append(dataclasses.replace(row, task=n))
    got = [r.suite for r in rows]
    if got != list(suites):
        raise PlanError(f"{path}: suite column {got} is not the replicates_full keys {suites}")
    return rows


def _group_tests(row: PlanRow, ids: Sequence[str], report: Mapping[str, bool]) -> set[str]:
    """The group tests (rule "group selection"); ``ids`` is the sorted collection."""
    sel = row.group_select
    if sel == "*":
        return set(ids)
    if sel == "-":
        only = [n for n in ids if report[n]]
        if only:
            raise PlanError(f"{row.suite}: group_select '-' but shard-0-only tests exist: "
                            + ", ".join(only))  # fmt: skip
        return set()
    prefixes = tuple(t.strip() for t in sel.split(",") if t.strip() != "direct")
    group = {
        n for n in ids if (row.shards > 0 and report[n]) or (prefixes and n.startswith(prefixes))
    }
    dropped = [n for n in ids if n not in group]
    if row.shards == 0 and dropped:  # no shard lane would run them ("No test is ever dropped")
        raise PlanError(f"{row.suite}: group_select {sel!r} with shards == 0 matches no prefix "
                        "for: " + ", ".join(dropped))  # fmt: skip
    return group


def _collisions(suite: str, what: str, deselect: Sequence[str], keep: Sequence[str]) -> None:
    """``--deselect`` is a prefix match: no deselected nodeid may be a prefix of a kept one."""
    hits = [f"{d} -> {k}" for k in keep for d in deselect if k.startswith(d)]
    if hits:
        raise PlanError(f"{suite}: {what} deselects a prefix of a test it keeps: "
                        + "; ".join(hits))  # fmt: skip


def lanes(row: PlanRow, report: Mapping[str, bool], durations: Mapping[str, float]) -> list[Lane]:
    """The lanes of ``row``: ``shards`` shard lanes, then ``groups`` group lanes (empty included).

    ``report`` maps each collected nodeid to ``shard0_only``; ``durations`` maps nodeids to
    seconds, a missing test weighs ``unknown_min * 60``. Group tests go by LPT (duration
    descending, nodeid ascending; smallest load, lowest lane on ties). A pure function of its
    inputs, independent of their order.
    """
    if row.task < 0:
        raise PlanError(f"{row.suite}: plan row has no task index (use read_plan)")
    ids = sorted(report)
    group = _group_tests(row, ids, report)
    shard_set = [n for n in ids if n not in group]
    if row.shards > 0 and not shard_set:
        raise PlanError(f"{row.suite}: shards = {row.shards} but no test is left for the shard "
                        "lanes (use '*')")  # fmt: skip

    def weight(n: str) -> float:
        return float(durations[n]) if n in durations else row.unknown_min * 60

    loads = [0.0] * row.groups
    members: list[list[str]] = [[] for _ in range(row.groups)]
    for n in sorted(group, key=lambda n: (-weight(n), n)):
        g = min(range(row.groups), key=lambda j: (loads[j], j))
        loads[g] += weight(n)
        members[g].append(n)

    group_ids = tuple(sorted(group))
    out: list[Lane] = []
    if row.shards > 0:
        _collisions(row.suite, "the shard lanes", group_ids, shard_set)
    for k in range(row.shards):
        out.append(Lane("shard", k, row.shards, group_ids, f"task-{row.task}.shard{k}.junit.xml",
                        tuple(shard_set)))  # fmt: skip
    for g, tests in enumerate(members):
        if not tests:
            out.append(Lane("group", g, row.groups, (), None, ()))
            continue
        keep = tuple(sorted(tests))
        drop = tuple(n for n in ids if n not in set(keep))
        _collisions(row.suite, f"group lane {g}", drop, keep)
        out.append(Lane("group", g, row.groups, drop, f"task-{row.task}.group{g}.junit.xml", keep))
    return out


def verify_partition(row: PlanRow, lane_list: Sequence[Lane],
                     collected: Mapping[str, Sequence[str]]) -> list[str]:  # fmt: skip
    """Problems of the per-lane ``--collect-only`` sets; ``[]`` = ok.

    ``collected`` maps each non-empty lane's ``part`` to its collected nodeids and ``"all"`` to
    the suite's full collection. ``lanes()`` builds disjoint group sets and identical shard sets,
    so each lane's own ``expected`` check covers overlaps and differing shard lanes.
    """
    problems: list[str] = []
    union: set[str] = set()
    for ln in lane_list:
        if ln.part is None:
            continue
        if ln.part not in collected:
            problems.append(f"{row.suite}: {ln.kind} lane {ln.index}: no collection for {ln.part}")
            continue
        got = list(collected[ln.part])
        have = set(got)
        if len(have) != len(got):
            problems.append(f"{row.suite}: {ln.part}: a nodeid is collected twice")
        if have != set(ln.expected):
            extra = sorted(have - set(ln.expected))
            short = sorted(set(ln.expected) - have)
            problems.append(f"{row.suite}: {ln.part} collected {len(have)} tests, planned "
                            f"{len(ln.expected)}; extra {extra}, missing {short}")  # fmt: skip
        union |= have
    if "all" not in collected:
        problems.append(f"{row.suite}: no full collection ('all')")
    else:
        everything = set(collected["all"])
        if union != everything:
            problems.append(f"{row.suite}: the lanes run {len(union)} tests, the collection has "
                            f"{len(everything)}; unplanned {sorted(everything - union)}, "
                            f"uncollected {sorted(union - everything)}")  # fmt: skip
    return problems


def verify_counts(lane_list: Sequence[Lane], junit_counts: Mapping[str, int]) -> list[str]:
    """Problems of the JUnit testcase counts (``part -> count``) against ``len(expected)``."""
    problems: list[str] = []
    for ln in lane_list:
        if ln.part is None:
            continue
        if ln.part not in junit_counts:
            problems.append(f"{ln.part}: no JUnit count ({len(ln.expected)} tests planned)")
        elif junit_counts[ln.part] != len(ln.expected):
            problems.append(f"{ln.part}: JUnit has {junit_counts[ln.part]} testcases, "
                            f"planned {len(ln.expected)}")  # fmt: skip
    return problems


def fold_exit(lane_exits: Sequence[int], merge_exit: int | None) -> int:
    """First exit not in {0, 1} (lanes, then merge); else 1 if any returned 1, else 0."""
    seq = [*lane_exits, *([] if merge_exit is None else [merge_exit])]
    for code in seq:
        if code not in (0, 1):
            return code
    return int(1 in seq)


def record_parts(lane_list: Sequence[Lane], merge_parts: Sequence[str]) -> list[str]:
    """The record's ``parts``: shard parts, non-empty group parts (in order), then the merge's."""
    shard = [ln.part for ln in lane_list if ln.kind == "shard" and ln.part is not None]
    group = [ln.part for ln in lane_list if ln.kind == "group" and ln.part is not None]
    return [*shard, *group, *merge_parts]


def read_durations(path: Path) -> dict[str, float]:
    """``nodeid -> seconds`` of a durations TSV; a bad row raises ``PlanError`` (a ValueError)."""
    lines = Path(path).read_text().splitlines()
    if not lines or tuple(lines[0].split("\t")) != DURATIONS_HEADER:
        raise PlanError(f"{path}: header must be {' '.join(DURATIONS_HEADER)}")
    out: dict[str, float] = {}
    for n, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        cells = line.split("\t")
        try:
            seconds = float(cells[1]) if len(cells) == 3 else math.nan
        except ValueError:
            seconds = math.nan
        if not (math.isfinite(seconds) and seconds >= 0) or not cells[0] or not cells[2]:
            raise PlanError(f"{path}: line {n}: need a nodeid, seconds >= 0 and a source")
        if cells[0] in out:
            raise PlanError(f"{path}: line {n}: {cells[0]} is listed twice")
        out[cells[0]] = seconds
    return out


def _junit_nodeid(classname: str, name: str) -> str:
    """``tests.calibration.test_x`` + ``name`` -> ``tests/calibration/test_x.py::name``
    (classname segments starting with an upper-case letter are test classes)."""
    parts = classname.split(".")
    cut = len(parts)
    while cut > 1 and parts[cut - 1][:1].isupper():
        cut -= 1
    return "::".join(["/".join(parts[:cut]) + ".py", *parts[cut:], name])


def durations_from_junit(paths: Iterable[Path]) -> dict[str, float]:
    """``nodeid -> seconds`` from pytest JUnit files; a nodeid in several files gets the sum."""
    out: dict[str, float] = {}
    for path in paths:
        for case in ET.parse(path).getroot().iter("testcase"):
            classname, name = case.get("classname"), case.get("name")
            if not classname or not name:
                continue
            nodeid = _junit_nodeid(classname, name)
            out[nodeid] = out.get(nodeid, 0.0) + float(case.get("time") or 0.0)
    return out


_STAMP = re.compile(r"([0-9]+)(?: (.*))?")
_PROBE_START = re.compile(r"(tests/\S+::\S+) ?")
_PROBE_END = re.compile(r"\[gw[0-9]+\] \[ *[0-9]+%\] (?:PASSED|FAILED|XFAIL|XPASS|SKIPPED|ERROR) "
                        r"(tests/\S+::\S+) ?")  # fmt: skip


def durations_from_probe_log(path: Path) -> dict[str, float]:
    """``nodeid -> seconds`` from a ``pytest -n N -v`` log whose lines carry an epoch-second stamp:
    from a test's start line to its outcome line. xdist ends each line only when the next event
    arrives, so a line's true time is the stamp of the previous stamped line. Unfinished tests are
    left out here; ``probe_log_lower_bounds`` gives them a lower bound."""
    start: dict[str, int] = {}
    out: dict[str, float] = {}
    prev: int | None = None
    for line in Path(path).read_text().splitlines():
        if not (stamp := _STAMP.fullmatch(line)):
            continue
        text, at = stamp.group(2) or "", prev
        prev = int(stamp.group(1))
        if at is None:
            continue
        if m := _PROBE_END.fullmatch(text):
            if m.group(1) in start:
                out[m.group(1)] = float(at - start[m.group(1)])
        elif m := _PROBE_START.fullmatch(text):
            start.setdefault(m.group(1), at)
    return out


def probe_log_lower_bounds(path: Path) -> dict[str, float]:
    """``nodeid -> seconds`` of the tests of a stamped probe log (as ``durations_from_probe_log``)
    that started but have no outcome line: the log's last stamp minus the start, a lower bound of
    the test's duration (it was still running when the log ended or the job was killed).

    Caveat: the bound assumes the test ran until the log's last stamp. A worker that died
    silently (no outcome line, no crash report) while the other workers kept logging gets a bound
    that overstates its duration; ``merge_lower_bounds`` and ``--max`` then keep that inflated
    value until a real measurement of the test replaces it."""
    start: dict[str, int] = {}
    done: set[str] = set()
    prev: int | None = None
    for line in Path(path).read_text().splitlines():
        if not (stamp := _STAMP.fullmatch(line)):
            continue
        text, at = stamp.group(2) or "", prev
        prev = int(stamp.group(1))
        if at is None:
            continue
        if m := _PROBE_END.fullmatch(text):
            done.add(m.group(1))
        elif m := _PROBE_START.fullmatch(text):
            start.setdefault(m.group(1), at)
    if prev is None:
        return {}
    return {n: float(prev - t) for n, t in start.items() if n not in done}


def merge_lower_bounds(bounds: Iterable[Mapping[str, float]]) -> dict[str, float]:
    """Merge the ``probe_log_lower_bounds`` of several probe logs: each test keeps its largest
    bound across the logs (every bound holds, so the largest is the tightest), whatever the
    order of the logs."""
    out: dict[str, float] = {}
    for b in bounds:
        for n, s in b.items():
            out[n] = max(float(s), out.get(n, float(s)))
    return out


def _write_durations(out: Path, found: Mapping[str, float], source: str, append: bool,
                     keep_max: bool = False,
                     lower: Mapping[str, float] | None = None) -> int:  # fmt: skip
    """Write ``out``; ``lower`` holds lower bounds (source ``<source>:lower-bound``), never
    replacing a measured value of ``found``. Returns the number of rows added or replaced.

    A row of ``out`` (with ``append`` or ``keep_max``) is replaced as follows; without
    ``keep_max`` an existing row is never replaced, only new tests are added. With ``keep_max``:

    - a measured value replaces a ``:lower-bound`` row, even when it is smaller (a measurement
      beats a bound);
    - otherwise (measured over measured, bound over bound, bound over measured) the new value
      replaces the row only when it is strictly larger."""
    rows: dict[str, tuple[str, str]] = {}
    if (append or keep_max) and out.exists():
        old = read_durations(out)
        sources = {ln.split("\t")[0]: ln.split("\t")[2] for ln in out.read_text().splitlines()[1:]
                   if ln.strip()}  # fmt: skip
        rows = {n: (repr(s), sources[n]) for n, s in old.items()}
    new = {n: (float(s), source) for n, s in found.items()}
    for n, s in (lower or {}).items():
        if n not in new:
            new[n] = (float(s), source + LOWER_BOUND)
    added = 0
    for nodeid, (seconds, src) in sorted(new.items()):
        if not nodeid.startswith(CALIBRATION):
            continue
        if nodeid not in rows or (keep_max and _replaces(rows[nodeid], seconds, src)):
            rows[nodeid] = (repr(seconds), src)
            added += 1
    body = "".join(f"{n}\t{s}\t{src}\n" for n, (s, src) in sorted(rows.items()))
    out.write_text("\t".join(DURATIONS_HEADER) + "\n" + body)
    return added


def _replaces(row: tuple[str, str], seconds: float, src: str) -> bool:
    """``keep_max`` rule of ``_write_durations``: does ``(seconds, src)`` replace ``row``?"""
    if row[1].endswith(LOWER_BOUND) and not src.endswith(LOWER_BOUND):
        return True
    return seconds > float(row[0])


def _suites() -> list[str]:
    with TOLERANCES.open() as fh:
        return [str(k) for k in yaml.safe_load(fh)["replicates_full"]]


def _lanes_cli(ns: argparse.Namespace) -> list[dict[str, object]]:
    rows = {r.suite: r for r in read_plan(ns.plan, _suites())}
    if ns.suite not in rows:
        raise PlanError(f"{ns.suite}: not a suite of the plan")
    raw = json.loads(Path(ns.report).read_text())
    if not isinstance(raw, dict):
        raise PlanError(f"{ns.report}: the lane report must be a JSON object")
    report: dict[str, bool] = {}
    for n, v in raw.items():
        flag = v.get("shard0_only") if isinstance(v, dict) else None
        if not isinstance(flag, bool):
            raise PlanError(f"{ns.report}: {n}: shard0_only must be a JSON bool, got {v!r}")
        report[str(n)] = flag
    path = ns.durations if ns.durations is not None else DURATIONS
    durations = {} if ns.durations is None and not path.exists() else read_durations(path)
    return [dataclasses.asdict(ln) for ln in lanes(rows[ns.suite], report, durations)]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    check = sub.add_parser("check", help="validate the plan")
    check.add_argument("--plan", type=Path, default=PLAN)
    lane = sub.add_parser("lanes", help="print a suite's lanes as JSON")
    lane.add_argument("--suite", required=True)
    lane.add_argument("--report", type=Path, required=True)
    lane.add_argument("--durations", type=Path, default=None)  # default: DURATIONS if present
    lane.add_argument("--plan", type=Path, default=PLAN)
    dur = sub.add_parser("durations", help="write the durations TSV")
    dur.add_argument("--source", required=True)
    dur.add_argument("--out", type=Path, required=True)
    mode = dur.add_mutually_exclusive_group()
    mode.add_argument("--append", action="store_true")
    mode.add_argument("--max", action="store_true", dest="keep_max")
    dur.add_argument("--probe-log", action="store_true")
    dur.add_argument("files", type=Path, nargs="+")
    ns = parser.parse_args(argv)
    try:
        if ns.cmd == "check":
            read_plan(ns.plan, _suites())
        elif ns.cmd == "lanes":
            print(json.dumps(_lanes_cli(ns), indent=1))
        else:
            found: dict[str, float] = {}
            lower: dict[str, float] = {}
            if ns.probe_log:
                for path in ns.files:
                    found.update(durations_from_probe_log(path))
                lower = merge_lower_bounds(probe_log_lower_bounds(p) for p in ns.files)
            else:
                found = durations_from_junit(ns.files)
            added = _write_durations(ns.out, found, ns.source, ns.append, ns.keep_max, lower)
            print(f"calibration_lanes: {added} rows from {ns.source} -> {ns.out}", file=sys.stderr)
    except (ValueError, OSError, ET.ParseError, KeyError, TypeError) as exc:  # all exit 2
        print(f"calibration_lanes: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
