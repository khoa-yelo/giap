"""Per-suite lane jobs of the full calibration tier (issue M11-03c).

Reads ``scripts/slurm/calibration_plan.tsv`` (one row per ``replicates_full`` suite of
``tests/calibration/tolerances.yaml``, in that order) and turns a row into the lanes one Slurm job
runs side by side: ``shards`` replicate-shard lanes and ``groups`` group lanes.

``read_plan``, ``fold_exit`` and ``record_parts`` are pure apart from ``read_plan`` reading the
plan file. ``lanes`` (the LPT assignment) is PR 1b; PR 1a only applies its ``-`` check.
Spec: docs/specs/M11-03c-calibration-full-sharding.md. ``check [--plan PATH]`` exits 0 when the
plan is valid and 2 on a ``PlanError`` (message on stderr).
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
PLAN = REPO / "scripts" / "slurm" / "calibration_plan.tsv"
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"
HEADER = ("suite", "shards", "groups", "group_select", "cpus_per_lane", "mem_per_lane_mb", "time",
          "unknown_min")  # fmt: skip
MAX_TIME_S = 8 * 3600
_UINT = re.compile(r"[0-9]+")
_TIME = re.compile(r"([0-9]{2}):([0-5][0-9]):([0-5][0-9])")


class PlanError(ValueError):
    """The plan (or a row against its collection) breaks a rule of "Plan validation"."""


@dataclass(frozen=True)
class PlanRow:
    """One plan row; ``time`` stays the ``HH:MM:SS`` string handed to ``sbatch --time``."""

    suite: str
    shards: int
    groups: int
    group_select: str
    cpus_per_lane: int
    mem_per_lane_mb: int
    time: str
    unknown_min: float

    @property
    def time_s(self) -> int:
        h, m, s = (int(x) for x in self.time.split(":"))
        return h * 3600 + m * 60 + s


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
        rows.append(_row(dict(zip(HEADER, cells, strict=True)), where))
    got = [r.suite for r in rows]
    if got != list(suites):
        raise PlanError(f"{path}: suite column {got} is not the replicates_full keys {suites}")
    return rows


def lanes(row: PlanRow, report: Mapping[str, bool], durations: Mapping[str, float]) -> list[Lane]:
    """The lanes of ``row``; ``report`` maps each collected nodeid to ``shard0_only``.

    PR 1a applies only the ``-`` check; the LPT assignment is PR 1b.
    """
    if row.group_select == "-" and row.shards > 0:
        only = sorted(nodeid for nodeid, shard0 in report.items() if shard0)
        if only:
            raise PlanError(f"{row.suite}: group_select '-' but shard-0-only tests exist: "
                            + ", ".join(only))  # fmt: skip
    raise NotImplementedError("lane assignment is M11-03c PR 1b")


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


def _suites() -> list[str]:
    with TOLERANCES.open() as fh:
        return [str(k) for k in yaml.safe_load(fh)["replicates_full"]]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    check = sub.add_parser("check", help="validate the plan")
    check.add_argument("--plan", type=Path, default=PLAN)
    ns = parser.parse_args(argv)
    try:
        read_plan(ns.plan, _suites())
    except (PlanError, OSError) as exc:
        print(f"calibration_lanes: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
