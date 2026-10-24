"""One run's ``bench_run.json`` (docs/specs/M10-12-benchmark-standard.md, "Run validity").

Pure functions (``parse_sacct``, ``trace_stats``, ``invalid_reasons``, ``run_status``) that never
raise on malformed input, and a thin ``main`` called by ``scripts/slurm/benchmark_standard.sbatch``
from the submitting checkout (``$GIAP_REPO``), so an older ``BENCH_COMMIT`` is still recorded:

    bench_record.py RUN_DIR GENOME VARIANT ORDER STEP RC T0 T1 FORCED SACCT_ROW BENCH_ID COMMIT
                    REPEAT ATTEMPT JOB_ID TASK NODE SKU CMD...

It writes ``RUN_DIR/bench_run.json`` and prints the status.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import statistics
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from measure_runtimes import parse_trace, process_name, rss_bytes

SCHEMA = "giap-bench-run/1"
SKU = "CPU_SKU:6148"
CPUS = 16
SACCT_KEYS = ("elapsed_s", "max_rss_bytes", "total_cpu_s", "alloc_cpus", "state", "exit")
_INT = re.compile(r"[0-9]+")
_RSS = re.compile(r"([0-9]*\.?[0-9]+)([KMGT]?)")


def _cpu_seconds(text: str) -> float | None:
    """sacct TotalCPU ``[D-][HH:]MM:SS[.mmm]`` → seconds; None if unparseable."""
    try:
        d, _, t = text.strip().rpartition("-")
        parts = [float(x) for x in t.split(":")]
        if not t or len(parts) > 3:
            return None
        h, m, s = [0.0, 0.0, *parts][-3:]
        return (int(d) if d else 0) * 86400 + h * 3600 + m * 60 + s
    except ValueError:
        return None


def _max_rss(text: str) -> int | None:
    """sacct MaxRSS (binary K/M/G/T) → bytes; empty or 0 → None, never 0."""
    m = _RSS.fullmatch(text.strip())
    b = int(float(m.group(1)) * 1024 ** "_KMGT".find(m.group(2) or "_")) if m else 0
    return b or None


def parse_sacct(row: str | None) -> dict[str, Any] | None:
    """One ``JobIDRaw|ElapsedRaw|MaxRSS|TotalCPU|AllocCPUS|State|ExitCode`` row; None for no row
    or an incomplete one (ElapsedRaw, AllocCPUS or State still empty), so the caller retries."""
    f = (row or "").strip().split("|")
    if len(f) < 7 or not _INT.fullmatch(f[1]) or not _INT.fullmatch(f[4]) or not f[5].split():
        return None
    return {"elapsed_s": int(f[1]), "max_rss_bytes": _max_rss(f[2]),
            "total_cpu_s": _cpu_seconds(f[3]), "alloc_cpus": int(f[4]),
            "state": f[5].split()[0], "exit": f[6]}  # fmt: skip


def _attempt(r: Mapping[str, str]) -> int:
    a = str(r.get("attempt") or "")
    return int(a) if _INT.fullmatch(a) else 1


def trace_stats(rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    """``trace.tsv`` rows → the record's ``trace`` and ``steps`` (ordered by first ``submit``)."""
    st = [str(r.get("status") or "").strip().upper() for r in rows]
    att = [_attempt(r) for r in rows]
    names = [str(r.get("name") or "") for r in rows]
    last: dict[str, int] = {}
    for n, a in zip(names, att, strict=True):
        last[n] = max(a, last.get(n, 0))
    dead = {n for n, s, a in zip(names, st, att, strict=True) if s == "FAILED" and a == last[n]}
    failed = [process_name(n) for n in dict.fromkeys(names) if n in dead]
    hosts = {str(r.get("hostname") or "").strip() for r in rows}
    peaks = [rss_bytes(str(r.get("peak_rss") or "")) for r, s in zip(rows, st, strict=True)
             if s != "CACHED"]  # fmt: skip
    steps: dict[str, dict[str, Any]] = {}
    first: dict[str, tuple[str, int]] = {}
    for i, r in enumerate(rows):
        for t in parse_trace([{k: str(v) for k, v in r.items()}]):
            s = steps.setdefault(t.process, {"realtime_s": 0.0, "peak_rss_bytes": 0})
            s["realtime_s"] += t.wall_s
            s["peak_rss_bytes"] = max(s["peak_rss_bytes"], int(t.rss_bytes))
            sub = str(r.get("submit") or "").strip()
            key = (sub if sub[:1].isdigit() else "~", i)  # no submit time: last
            first[t.process] = min(first.get(t.process, key), key)
    for k, n in enumerate(sorted(steps, key=lambda n: first[n])):
        steps[n]["order"] = k
    return {"trace": {"n_tasks": len(rows), "n_cached": st.count("CACHED"), "n_failed": len(failed),
                      "n_retried": sum(a > 1 for a in att), "hosts": sorted(hosts - {"", "-"}),
                      "failed_processes": failed,
                      "task_peak_rss_max_bytes": max([p for p in peaks if p] or [0]) or None},
            "steps": steps}  # fmt: skip


def _short(host: str) -> str:
    return host.strip().split(".")[0]


def invalid_reasons(*, variant: str, commit: str, node: str, sku: str,
                    trace_rows: Sequence[Mapping[str, str]], manifest: Mapping[str, Any] | None,
                    sacct: Mapping[str, Any] | None, giap_other: Iterable[Mapping[str, Any]],
                    caches_before: Iterable[str], caches_after: Iterable[str],
                    ) -> tuple[list[str], list[str]]:  # fmt: skip
    """The run's invalid reasons (all of them) and its warnings (``contended_other_nodes``)."""
    reasons: list[str] = []
    warnings: list[str] = []
    n_cached = sum(str(r.get("status") or "").strip().upper() == "CACHED" for r in trace_rows)
    if n_cached:
        reasons.append(f"cached_tasks:{n_cached}")
    if manifest is None:
        reasons.append("manifest_missing")
    else:
        if manifest.get("resume") is not False:
            reasons.append("manifest_resume")
        steps = manifest.get("steps")
        n_cs = sum(isinstance(s, Mapping) and s.get("status") == "cached" for s in steps or [])
        if n_cs:
            reasons.append(f"cached_steps:{n_cs}")
        gc = str(manifest.get("git_commit") or "")
        if not gc:
            reasons.append("commit_missing")
        other = gc and not (len(gc) >= 7 and commit.startswith(gc))
        if manifest.get("git_dirty") or other:
            reasons.append("commit_mismatch")
    hosts = {str(r.get("hostname") or "").strip() for r in trace_rows} - {"", "-"}
    if any(_short(h) != _short(node) for h in hosts):
        reasons.append("host_mismatch")
    if sku != SKU:
        reasons.append(f"node_class:{sku or 'unknown'}")
    procs = {process_name(str(r.get("name") or "")) for r in trace_rows}
    if variant == "core" and "bakta" in procs:
        reasons.append("bakta_ran_in_core")
    for o in giap_other:
        if _short(str(o.get("node") or "")) == _short(node):
            reasons.append(f"same_node_giap:{o.get('job')}")
        elif int(o.get("cpus") or 0) >= 8 and "contended_other_nodes" not in warnings:
            warnings.append("contended_other_nodes")
    if set(caches_after) - set(caches_before) - {""}:
        reasons.append("image_or_env_created")
    if sacct is None:
        reasons.append("sacct_missing")
    elif sacct.get("alloc_cpus") != CPUS:
        reasons.append(f"alloc_cpus:{sacct.get('alloc_cpus')}")
    return reasons, warnings


def run_status(rc: int, sacct_state: str | None, reasons: Sequence[str],
               forced: str | None = None) -> str:  # fmt: skip
    """forced > failed (rc not 0/124/137, or OUT_OF_MEMORY) > timeout (124, 137) > invalid > ok."""
    if forced:
        return forced
    if str(sacct_state or "").startswith("OUT_OF_MEMORY") or rc not in (0, 124, 137):
        return "failed"
    return "timeout" if rc else "invalid" if reasons else "ok"


def _lines(p: Path) -> list[str]:
    return p.read_text().split("\n") if p.is_file() else []


def main(argv: list[str]) -> int:
    (rd_s, genome, variant, order, step, rc_s, t0, t1, forced, row, bench_id, commit, repeat,
     attempt, job_id, task, node, sku, *cmd) = argv  # fmt: skip
    rd, rc = Path(rd_s), int(rc_s)
    gargv = cmd[cmd.index("giap") :] if "giap" in cmd else cmd
    out = rd / "out"
    tp = out / "logs" / "trace.tsv"
    rows = list(csv.DictReader(tp.open(), delimiter="\t")) if tp.is_file() else []
    try:
        manifest = json.loads((out / "run_manifest.json").read_text())
    except (OSError, ValueError):
        manifest = None
    samples = [ln.split("\t") for ln in _lines(rd / "cotenancy.tsv") if ln.strip()]
    allocs = [int(s[1]) for s in samples if len(s) > 1 and s[1].isdigit()]
    others: dict[str, dict[str, Any]] = {}
    for s in samples:
        for item in filter(None, (s[4] if len(s) > 4 else "").split(";")):
            j, name, cpus, where = [*item.split("|"), "", "", ""][:4]
            n = int(cpus) if cpus.isdigit() else 0
            others.setdefault(j, {"job": j, "name": name, "cpus": n, "node": where})
    sa = parse_sacct(row)
    reasons, warnings = invalid_reasons(
        variant=variant, commit=commit, node=node, sku=sku, trace_rows=rows,
        manifest=manifest if isinstance(manifest, dict) else None, sacct=sa,
        giap_other=list(others.values()), caches_before=_lines(rd / "caches.before"),
        caches_after=_lines(rd / "caches.after"))  # fmt: skip
    script_wall = round(float(t1) - float(t0), 3)
    wall: float = script_wall if sa is None else sa["elapsed_s"]
    if sa is None:
        warnings.append("wall_from_script: sacct gave no complete row for the step")
    elif abs(wall - script_wall) > 5:
        warnings.append(f"wall_mismatch: sacct {wall} s, script {script_wall} s")
    status = run_status(rc, None if sa is None else str(sa["state"]), reasons, forced or None)
    inputs = {}
    for name in ("genome", "peptides", "config", "annotation"):
        p = Path(gargv[gargv.index("--" + name) + 1]) if "--" + name in gargv[:-1] else None
        if p and p.is_file():
            inputs[name] = {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
    ts = trace_stats(rows)
    rec = {"schema": SCHEMA, "bench_id": bench_id, "commit": commit, "genome": genome,
           "variant": variant, "repeat": int(repeat), "attempt": attempt, "order_index": int(order),
           "job_id": job_id, "array_task": int(task), "step_id": int(step), "node": node,
           "cpu_sku": sku, "argv": gargv, "inputs": inputs,
           "sacct": sa or dict.fromkeys(SACCT_KEYS),
           "script_wall_s": script_wall, "wall_s": wall, "exit_code": rc, "status": status,
           "invalid_reasons": reasons, "warnings": warnings, **ts,
           "cotenancy": {"samples": len(samples), "giap_other": list(others.values()),
                         "node_cpualloc_median": statistics.median_low(allocs) if allocs else None,
                         "node_cpualloc_max": max(allocs) if allocs else None}}  # fmt: skip
    (rd / "bench_run.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
    print(status)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
