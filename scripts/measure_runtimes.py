#!/usr/bin/env python3
"""Measured runtimes from a Nextflow ``trace.tsv`` (M2-21; GIPPA_SPEC.md §3.8, §13 M2 gate).

Usage: ``measure_runtimes.py TRACE.tsv [--sacct JOBID] [--genome NAME]``

Prints markdown rows ``step | genome | wall | peak RSS | node class (job ID)`` for
``PROGRESS.md`` > Measured numbers (lead decision 6): one row per process (the last ``:``
component of the trace ``name``, lowercased), the largest peak RSS and the longest wall time when a
process ran several tasks (``promoter_calculator (16 tasks)``), ``CACHED``/``FAILED`` tasks left
out. ``wall`` is the task ``realtime`` (else ``duration``) as ``H:MM:SS``; peak RSS keeps the
trace's unit. The node class (``wheat``/``barley``/``rye``/``oat``) is the prefix of the node
each task ran on, from the trace ``hostname`` field (lead decisions 16 and 20; under ``-profile
farmshare`` the tasks run on other nodes than the head); tasks without it (older traces, or
``-``) fall back to the node that ``sacct -j JOBID`` reports for the head job.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

_DUR = re.compile(r"([0-9.]+)\s*(ms|d|h|m|s)")
_DUR_S = {"ms": 1e-3, "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0}
_RSS = re.compile(r"([0-9.]+)\s*([KMGT]?B)")
_RSS_B = {"B": 1, "KB": 1024, "MB": 1024**2, "GB": 1024**3, "TB": 1024**4}
_SKIP = {"CACHED", "FAILED", "ABORTED"}


@dataclass(frozen=True)
class Task:
    """One trace row: process name, wall seconds, peak RSS (text and bytes), Bakta DB type."""

    process: str
    wall_s: float
    rss_text: str
    rss_bytes: float
    variant: str = ""
    node: str = ""  # node class of the task's own node ("" when the trace has no hostname)


def seconds(text: str) -> float | None:
    """Nextflow duration (``13m 15s``, ``50.4s``, ``746ms``) → seconds; None if absent."""
    parts = _DUR.findall(text or "")
    return sum(float(n) * _DUR_S[u] for n, u in parts) if parts else None


def hms(secs: float) -> str:
    """Seconds → ``H:MM:SS`` (fraction truncated)."""
    s = int(secs)
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


def rss_bytes(text: str) -> float | None:
    """``11 GB`` / ``601.4 MB`` → bytes; None if absent."""
    m = _RSS.fullmatch((text or "").strip())
    return float(m.group(1)) * _RSS_B[m.group(2)] if m else None


def process_name(name: str) -> str:
    """``ANNOTATE:PROMOTER_CALCULATOR (TU_12)`` → ``promoter_calculator``."""
    return name.split(" (")[0].split(":")[-1].strip().lower()


def node_class(node: str) -> str:
    """``barley-01`` / ``wheat-[03-04]`` → ``barley`` / ``wheat``."""
    m = re.match(r"[A-Za-z]+", node.strip())
    return m.group(0).lower() if m else "?"


def bakta_variant(workdir: str) -> str:
    """``full`` or ``light`` from the task's ``.command.sh`` when it is still on disk."""
    try:
        cmd = (Path(workdir) / ".command.sh").read_text()
    except OSError:
        return ""
    if "db-light" in cmd:
        return "light"
    return "full" if re.search(r"--db\s+\S*/db\b", cmd) else ""


def parse_trace(rows: list[dict[str, str]]) -> list[Task]:
    """Trace rows → measured tasks (CACHED/FAILED rows and rows without a time dropped)."""
    tasks = []
    for r in rows:
        if r.get("status", "").strip().upper() in _SKIP:
            continue
        wall = seconds(r.get("realtime", ""))
        if wall is None:
            wall = seconds(r.get("duration", ""))
        if wall is None:
            continue
        proc = process_name(r.get("name", ""))
        rss = r.get("peak_rss", "").strip()
        variant = bakta_variant(r.get("workdir", "")) if proc == "bakta" else ""
        host = (r.get("hostname") or "").strip()
        node = node_class(host) if host and host != "-" else ""
        tasks.append(Task(proc, wall, rss or "-", rss_bytes(rss) or 0.0, variant, node))
    return tasks


def summarise(tasks: list[Task], genome: str, node: str, job: str | None = None) -> list[str]:
    """One markdown row per process, in order of first appearance.

    ``node`` is the head job's class, used for tasks without their own; ``job`` labels the cell.
    """
    by: dict[str, list[Task]] = {}
    for t in tasks:
        by.setdefault(t.process, []).append(t)
    out = ["| step | genome | wall | peak RSS | node class (job ID) |", "|---|---|---|---|---|"]
    for proc, ts in by.items():
        top = max(ts, key=lambda t: t.rss_bytes)
        notes = [v for v in sorted({t.variant for t in ts}) if v]
        if len(ts) > 1:
            notes.append(f"{len(ts)} tasks")
        step = f"{proc} ({', '.join(notes)})" if notes else proc
        wall = hms(max(t.wall_s for t in ts))
        klass = "/".join(sorted({t.node or node for t in ts}))
        cell = f"{klass} ({job})" if job else klass
        out.append(f"| {step} | {genome} | {wall} | {top.rss_text} | {cell} |")
    return out


def sacct_node(job: str) -> str:
    """Node name of the job's first allocation row from ``sacct``."""
    res = subprocess.run(
        ["sacct", "-j", job, "--format=JobID,NodeList", "-P", "-n"],
        capture_output=True,
        text=True,
        check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(f"sacct -j {job} failed: {res.stderr.strip()}")
    for line in res.stdout.splitlines():
        cells = line.split("|")
        if len(cells) >= 2 and cells[1].strip() and cells[1].strip() != "None assigned":
            return cells[1].strip()
    return "?"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("trace", type=Path, help="Nextflow trace.tsv")
    ap.add_argument("--sacct", metavar="JOBID", help="Slurm job that ran the head (node class)")
    ap.add_argument("--genome", default="-", help="genome label for the rows")
    args = ap.parse_args(argv)
    try:
        with args.trace.open(newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    except OSError as exc:
        print(f"measure_runtimes: cannot read {args.trace}: {exc.strerror}", file=sys.stderr)
        return 2
    tasks = parse_trace(rows)
    node = "?"
    if args.sacct and any(not t.node for t in tasks):
        try:
            node = node_class(sacct_node(args.sacct))
        except (OSError, RuntimeError) as exc:
            print(f"measure_runtimes: {exc}", file=sys.stderr)
            return 2
    print("\n".join(summarise(tasks, args.genome, node, args.sacct)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
