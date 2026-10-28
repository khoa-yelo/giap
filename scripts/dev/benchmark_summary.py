"""Standard-run benchmark summary (docs/specs/M10-12-benchmark-standard.md, "Interface").

Reads every ``r<k>[.<attempt>]/<genome>_<variant>/bench_run.json`` written by
``scripts/slurm/benchmark_standard.sbatch`` under one benchmark directory and writes
``summary.json`` and ``summary.md`` (the PROGRESS > Measured numbers rows and the user-guide
table) next to them. Per (genome, variant) cell the latest valid attempt of each repeat counts;
a cell with 3 valid repeats and equal input sha256 values is ``complete`` (median of 3), 1-2 is
``incomplete`` (no median; the mean of the valid ones, labelled), 0 is ``no_data``. With
``--write`` the complete cells go into ``standard_run.<node class>.<genome>_<variant>`` of
``baseline.json``; every other key is kept byte-for-byte, and an entry measured on another
commit is replaced only with ``--replace`` (else exit 1, file unchanged). The noisy threshold
is ``benchmark_slowdown_max`` of ``ci/gates.yaml``. Exit 2 on an unreadable or schema-invalid
record (a minimal ``record_error`` record is accepted: a failed repeat) or a bad ``--adr``.
The numbers are reported against GIPPA_SPEC.md §13.2, never gated.

    python scripts/dev/benchmark_summary.py BENCH_DIR/<bench_id> [--baseline B --write
        [--replace] [--adr ADR-NNNN]]
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

REPO = Path(__file__).resolve().parents[2]
NOISY_MAX: float = float(
    yaml.safe_load((REPO / "ci" / "gates.yaml").read_text())["benchmark_slowdown_max"]
)
SCHEMA = "gippa-bench-run/1"
RUN_KEYS = ("schema", "bench_id", "commit", "genome", "variant", "repeat", "attempt", "argv",
            "inputs", "sacct", "wall_s", "status", "trace", "steps", "node", "cpu_sku",
            "job_id")  # fmt: skip
MIN_KEYS = ("schema", "genome", "variant", "repeat", "attempt", "status")  # a record_error record
SAME_FLAGS = ("--profile", "--memory", "--modules", "--threads", "--seed", "--set")
GENOMES = ("k12", "sim5m")
VARIANTS = ("standard", "core")
ROUNDS = 3
THREADS = 16
NODE_CLASSES = {"CPU_SKU:6148": ("wheat", "Intel Xeon Gold 6148 (CPU_SKU:6148)")}
Run = Mapping[str, Any]


@dataclass(frozen=True)
class Estimate:
    """GIPPA_SPEC.md §13.2 runtime estimate of one variant (seconds; range when one is given)."""

    seconds: float
    low: float | None = None
    high: float | None = None


ESTIMATES = {"standard": Estimate(2100, 960, 3300), "core": Estimate(1500)}


def run_name(run: Run) -> str:
    """Run-directory name of a record's repeat: ``r2`` (attempt ``a``), ``r2.b``."""
    att = str(run.get("attempt") or "a")
    return f"r{run['repeat']}" + ("" if att == "a" else f".{att}")


def select_attempts(runs: Sequence[Run]) -> tuple[dict[int, Run], list[str], list[Run]]:
    """Per repeat, the latest ``ok`` attempt; the other attempts of those repeats (superseded);
    and the latest attempt of each repeat that has no valid one."""
    by: dict[int, list[Run]] = {}
    for r in runs:
        by.setdefault(int(r["repeat"]), []).append(r)
    chosen: dict[int, Run] = {}
    superseded: list[str] = []
    unresolved: list[Run] = []
    for k in sorted(by):
        atts = sorted(by[k], key=lambda r: str(r.get("attempt") or "a"))
        ok = [r for r in atts if r["status"] == "ok"]
        if ok:
            chosen[k] = ok[-1]
            superseded += [run_name(r) for r in atts if r is not ok[-1]]
        else:
            unresolved.append(atts[-1])
    return chosen, superseded, unresolved


def _median(xs: Sequence[float]) -> float | None:
    return statistics.median(xs) if xs else None


def _max(xs: Sequence[float | None]) -> float | None:
    vals = [x for x in xs if x is not None]
    return max(vals) if vals else None


def _steps(valid: Sequence[Run], median_run: Run) -> dict[str, dict[str, Any]]:
    """Per step: median of realtime over repeats, max peak RSS; median repeat's start order."""
    names: list[str] = []
    for r in (median_run, *valid):
        order = sorted(r["steps"].items(), key=lambda kv: kv[1].get("order", 0))
        names += [n for n, _ in order if n not in names]
    out = {}
    for n in names:
        rows = [r["steps"][n] for r in valid if n in r["steps"]]
        out[n] = {
            "median_s": _median([s["realtime_s"] for s in rows]),
            "peak_rss_max_bytes": _max([s.get("peak_rss_bytes") for s in rows]),
        }
    return out


def _flag(argv: Sequence[str], name: str) -> str | None:
    return next((argv[i + 1] for i, a in enumerate(argv[:-1]) if a == name), None)


def cell_summary(runs: Sequence[Run], estimate: Estimate) -> dict[str, Any]:
    """Every record of one cell (all repeats and attempts) → the cell's numbers and status."""
    cell: dict[str, Any] = {
        "status": "no_data",
        "median_s": None,
        "mean_of_valid_s": None,
        "spread": None,
        "noisy": None,
        "ratio_to_estimate": None,
        "within_estimate_range": None,
        "superseded": [],
        "reasons": [],
        "not_valid": [],
        "estimate_s": estimate.seconds,
    }
    if runs:
        cell |= {"genome": runs[0]["genome"], "variant": runs[0]["variant"]}
    chosen, cell["superseded"], unresolved = select_attempts(runs)
    cell["not_valid"] = [
        {
            "run": run_name(r),
            "status": r["status"],
            "exit_code": r.get("exit_code"),
            "invalid_reasons": r.get("invalid_reasons", []),
            "failed_processes": r.get("trace", {}).get("failed_processes", []),
        }
        for r in unresolved
    ]
    valid = [chosen[k] for k in sorted(chosen)]
    if not valid:
        return cell
    walls = [float(r["wall_s"]) for r in valid]
    same = {
        json.dumps([{n: v["sha256"] for n, v in r["inputs"].items()}, r["commit"],
                    [_flag(r["argv"], f) for f in SAME_FLAGS]], sort_keys=True)
        for r in valid
    }  # fmt: skip
    classes = {r["cpu_sku"] for r in valid}
    cell["reasons"] = (
        (["inputs_differ"] if len(same) > 1 else [])
        + (["mixed_node_class"] if len(classes) > 1 else [])
        + [f"node_class:{c}" for c in sorted(classes) if c not in NODE_CLASSES]  # never summarised
    )
    cell["valid_runs"] = len(valid)
    if len(valid) != ROUNDS or cell["reasons"]:
        cell |= {"status": "incomplete", "mean_of_valid_s": statistics.fmean(walls)}
        return cell
    med = statistics.median(walls)
    median_run = valid[walls.index(med)]
    sku = valid[0]["cpu_sku"]
    node_class, cpu = NODE_CLASSES[sku]
    within = None
    if estimate.low is not None and estimate.high is not None:
        within = estimate.low <= med <= estimate.high
    spread = (max(walls) - min(walls)) / med
    eff = [
        r["sacct"]["total_cpu_s"] / (r["sacct"]["elapsed_s"] * THREADS)
        for r in valid
        if r["sacct"].get("total_cpu_s") is not None and r["sacct"].get("elapsed_s")
    ]
    argv = list(valid[0]["argv"])
    day = str(valid[0]["bench_id"])[:10]
    cell |= {
        "status": "complete",
        "median_s": med,
        "spread": spread,
        "noisy": spread > NOISY_MAX,
        "ratio_to_estimate": med / estimate.seconds,
        "within_estimate_range": within,
        "node_class": node_class,
        "runs_s": [r["wall_s"] for r in valid],
        "nodes": [r["node"] for r in valid],
        "cpu": cpu,
        "profile": _flag(argv, "--profile"),
        "modules": _flag(argv, "--modules") or "core,mobile,bgc,function",
        "bakta": "recorded JSON" if _flag(argv, "--annotation") else "live, full DB",
        "peak_rss_sacct_max_bytes": _max([r["sacct"].get("max_rss_bytes") for r in valid]),
        "task_peak_rss_max_bytes": _max([r["trace"].get("task_peak_rss_max_bytes") for r in valid]),
        "cpu_efficiency_median": _median(eff),
        "steps": _steps(valid, median_run),
        "slurm_jobs": [r["job_id"] for r in valid],
        "commit": valid[0]["commit"],
        "date": day if re.fullmatch(r"\d{4}-\d{2}-\d{2}", day) else None,
        "inputs_sha256": {n: v["sha256"] for n, v in valid[0]["inputs"].items()},
    }
    return cell


def baseline_entry(cell: Mapping[str, Any]) -> dict[str, Any]:
    """The ``standard_run`` entry of a complete cell (spec "Baseline entry")."""
    if cell["status"] != "complete":
        raise ValueError(f"only a complete cell has a baseline entry, not {cell['status']}")
    keys = ("median_s", "runs_s", "spread", "noisy", "nodes", "cpu", "profile", "modules", "bakta",
            "peak_rss_sacct_max_bytes", "task_peak_rss_max_bytes", "cpu_efficiency_median", "steps",
            "estimate_s", "ratio_to_estimate", "within_estimate_range", "slurm_jobs", "commit",
            "date", "inputs_sha256")  # fmt: skip
    entry = {k: cell[k] for k in keys}
    return entry | {
        "threads": THREADS,
        "mem_mb": 64000,
        "task_memory": "56.GB",
        "rounds": ROUNDS,
        "adr": cell.get("adr"),
    }


def _hms(s: float) -> str:
    s = round(s)
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


def _gib(b: float | None) -> str:
    return "-" if b is None else f"{b / 1024**3:.1f} GiB"


def markdown_rows(cells: Sequence[Mapping[str, Any]]) -> list[str]:
    """PROGRESS > Measured numbers rows, one per cell (``step = gippa run <variant>``)."""
    rows = []
    for c in cells:
        notes = list(c.get("reasons", []))
        for nv in c.get("not_valid", []):
            why = ", ".join(nv["invalid_reasons"] or nv["failed_processes"])
            notes.append(
                f"{nv['run']} {nv['status']} (exit {nv['exit_code']}"
                + (f"; {why})" if why else ")")
            )
        if c["status"] == "complete":
            wall = (
                f"{_hms(c['median_s'])} (median of 3; {_hms(min(c['runs_s']))}-"
                f"{_hms(max(c['runs_s']))})"
            )
            if c["noisy"]:
                notes.append(f"noisy (spread {c['spread'] * 100:.1f} % > {NOISY_MAX * 100:g} %)")
            notes.append(f"{c['ratio_to_estimate']:.2f} x the §13.2 estimate")
            node = f"{c['node_class']} ({', '.join(map(str, c['slurm_jobs']))})"
            rss = _gib(c["peak_rss_sacct_max_bytes"])
        elif c["status"] == "incomplete":
            wall = f"{_hms(c['mean_of_valid_s'])} (mean of {c['valid_runs']} of 3 repeats)"
            node, rss = "-", "-"
        else:
            wall, node, rss = "no data", "-", "-"
        cell = f"{c.get('genome', '?')}"
        note = f" ({'; '.join(notes)})" if notes else ""
        rows.append(
            f"| gippa run {c.get('variant', '?')} | {cell} | {wall}{note} | {rss} | {node} |"
        )
    return rows


def guide_rows(cells: Sequence[Mapping[str, Any]]) -> list[str]:
    """User-guide runtime table: the complete ``k12`` cells only (lead decisions L8, L9)."""
    rows = [
        "| variant | median (h:mm) | range of 3 runs | peak memory | node class | date | commit |",
        "|---|---|---|---|---|---|---|",
    ]
    for c in cells:
        if c.get("genome") == "k12" and c["status"] == "complete":
            m = round(c["median_s"] / 60)
            rows.append(
                f"| {c['variant']} | {m // 60}:{m % 60:02d} | {_hms(min(c['runs_s']))}-"
                f"{_hms(max(c['runs_s']))} | {_gib(c['peak_rss_sacct_max_bytes'])} | "
                f"{c['node_class']} | {c['date']} | {str(c['commit'])[:7]} |"
            )
    return rows


def _load(bench: Path) -> list[dict[str, Any]]:
    runs = []
    for p in sorted(bench.glob("r*/*/bench_run.json")):
        try:
            rec = json.loads(p.read_text())
        except (OSError, ValueError) as exc:
            raise ValueError(f"{p}: unreadable ({exc})") from exc
        err = isinstance(rec, dict) and any(
            str(x).startswith("record_error") for x in rec.get("invalid_reasons") or []
        )
        missing = [
            k for k in (MIN_KEYS if err else RUN_KEYS) if not isinstance(rec, dict) or k not in rec
        ]
        if missing or rec.get("schema") != SCHEMA:
            raise ValueError(f"{p}: not a {SCHEMA} record (missing {missing})")
        runs.append(rec)
    return runs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("bench", type=Path, help="BENCH_DIR/<bench_id>")
    ap.add_argument("--baseline", type=Path, default=REPO / "benchmarks" / "baseline.json")
    ap.add_argument("--write", action="store_true", help="add the complete cells to --baseline")
    ap.add_argument("--replace", action="store_true", help="replace entries of another commit")
    ap.add_argument("--adr", help="ADR-NNNN of the entries (kept on a same-commit rewrite)")
    args = ap.parse_args(argv)
    if args.adr is not None and not re.fullmatch(r"ADR-\d{4}", args.adr):
        print(f"benchmark_summary: --adr {args.adr!r} is not ADR-NNNN", file=sys.stderr)
        return 2
    try:
        runs = _load(args.bench)
    except ValueError as exc:
        print(f"benchmark_summary: {exc}", file=sys.stderr)
        return 2
    cells = {}
    for g in GENOMES:
        for v in VARIANTS:
            group = [r for r in runs if (r["genome"], r["variant"]) == (g, v)]
            cells[f"{g}_{v}"] = {"genome": g, "variant": v} | cell_summary(group, ESTIMATES[v])
    summary = {"bench_id": args.bench.name, "noisy_max": NOISY_MAX, "cells": cells}
    (args.bench / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    md = [
        "| step | genome | wall | peak RSS | node class (job ID) |",
        "|---|---|---|---|---|",
        *markdown_rows(list(cells.values())),
        "",
        *guide_rows(list(cells.values())),
    ]
    (args.bench / "summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    entries = {k: c for k, c in cells.items() if c["status"] == "complete"}
    if not args.write or not entries:
        return 0
    data = json.loads(args.baseline.read_text())
    sr = data.get("standard_run", {})
    for key, c in entries.items():
        old = sr.get(c["node_class"], {}).get(key)
        if old and old.get("commit") != c["commit"] and not args.replace:
            print(
                f"benchmark_summary: standard_run.{c['node_class']}.{key} is from commit "
                f"{old.get('commit')}; pass --replace (an ADR's update)",
                file=sys.stderr,
            )
            return 1
    for key, c in entries.items():
        old = sr.get(c["node_class"], {}).get(key) or {}
        keep = old.get("adr") if old.get("commit") == c["commit"] else None
        sr.setdefault(c["node_class"], {})[key] = baseline_entry(c | {"adr": args.adr or keep})
    data["standard_run"] = sr
    args.baseline.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
