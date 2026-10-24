"""Evidence for the pass-1 pilot runs: rule sets A (one arm) and C (the compare) of M10-04.

Spec: docs/specs/M10-04-pilot-runs.md ("run_evidence.py", lead decisions of 2026-10-24). The
evidence holds provenance only (paths, sha256s, byte sizes, exit codes, states, resources), never
a value from the tables (Khoa's 2026-10-23 rule).

Usage (``PYTHONPATH=$GIAP_REPO/src``; exit 0 when ``failures`` is empty, else 1)::

    python validation/pilot/run_evidence.py arm OUT --arm A1|A2 --job JOBID --expect-commit SHA \
        --expect-input genome=SHA --expect-input annotation=SHA --expect-input peptides=SHA \
        --out FILE
    python validation/pilot/run_evidence.py compare CMP --second CMP2 --labels A1,A2 --out FILE
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections.abc import Mapping, Sequence
from importlib import resources
from pathlib import Path
from typing import Any

MAX_CPUS = 14  # rule A6: one arm at --threads 14 --memory 54.GB (lead decision 1, Q1)
CHECKED_INPUTS = ("genome", "peptides")  # rule A5 (lead decision 3: annotation recorded only)
RECORDED_BAKTA = {"bakta_version": "1.12.2", "db_type": "full"}  # rule A4
COMPARE_TABLES = ("compare_detection", "compare_abundance", "compare_stability")
MANIFEST_KEYS = (
    "status", "exit_code", "failed_step", "failed_log", "git_commit", "git_dirty", "profile",
    "resume", "started_at", "finished_at", "seeds", "modules_off", "detection_tier",
    "annotation", "databases", "warnings",
)  # fmt: skip
HEAD_FIELDS = "JobID,State,ExitCode,Elapsed,MaxRSS,NodeList"
TASK_FIELDS = "JobID,AllocCPUS,MaxRSS,State"
_RSS_UNIT = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}


# ------------------------------------------------------------------ small pure helpers


def file_entry(path: Path) -> dict[str, Any]:
    """Path, sha256 and byte size of one file (its content is never copied)."""
    data = path.read_bytes()
    return {"path": str(path), "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}


def rss_bytes(text: str | None) -> int | None:
    """``sacct`` MaxRSS (``2104832K``, ``1.5G``, ``''``) in bytes."""
    m = re.fullmatch(r"([\d.]+)([KMGT]?)", (text or "").strip())
    return None if m is None else int(float(m.group(1)) * _RSS_UNIT[m.group(2)])


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _job_rows(rows: Sequence[Mapping[str, str]], job: str) -> list[Mapping[str, str]]:
    """The rows of ``job`` and its steps (``job.batch``, ``job.extern``)."""
    return [r for r in rows if r.get("JobID") == job or r.get("JobID", "").startswith(job + ".")]


def _max_rss(rows: Sequence[Mapping[str, str]]) -> str:
    """The largest MaxRSS of the rows, as ``sacct`` printed it ('' when none)."""
    best = max(rows, key=lambda r: rss_bytes(r.get("MaxRSS")) or -1, default=None)
    return "" if best is None else best.get("MaxRSS", "")


def _process(name: str) -> str:
    """``GIAP:ANNOTATION:BAKTA (genome)`` → ``BAKTA``."""
    return name.split(" ", 1)[0].rsplit(":", 1)[-1]


def manifest_errors(doc: Any) -> list[str]:
    """JSON pointers (with a missing property's name) of the schema violations of ``doc``."""
    import jsonschema  # type: ignore[import-untyped]

    res = resources.files("giap.schemas").joinpath("run_manifest.schema.json")
    schema = json.loads(res.read_text(encoding="utf-8"))
    out = []
    for e in jsonschema.Draft202012Validator(schema).iter_errors(doc):
        ptr = "/" + "/".join(str(p) for p in e.absolute_path)
        if e.validator == "required" and isinstance(e.instance, dict):
            missing = [p for p in e.validator_value if p not in e.instance]
            ptr = ptr.rstrip("/") + "/" + ",".join(missing)
        out.append(ptr)
    return sorted(set(out))


# ------------------------------------------------------------------ rule set A


def _task_summary(
    trace: Sequence[Mapping[str, str]], sacct: Sequence[Mapping[str, str]]
) -> dict[str, Any]:
    ids = sorted({r["native_id"] for r in trace if r.get("native_id", "-") not in ("", "-")})
    rows = [row for jid in ids for row in _job_rows(sacct, jid)]
    cpus = [int(r["AllocCPUS"]) for r in rows if str(r.get("AllocCPUS", "")).isdigit()]
    rss = [b for b in (rss_bytes(r.get("MaxRSS")) for r in rows) if b is not None]
    return {
        "count": len(trace),
        "failed": sum(r.get("status") in ("FAILED", "ABORTED") for r in trace),
        "retried": sum(int(r.get("attempt") or 1) > 1 for r in trace),
        "max_alloc_cpus": max(cpus, default=0),
        "peak_max_rss_bytes": max(rss, default=None),
    }


def _arm_failures(
    out: Path, man: Mapping[str, Any] | None, job: Mapping[str, Any],
    trace: Sequence[Mapping[str, str]], tasks: Mapping[str, Any], expect: Mapping[str, Any],
) -> list[str]:  # fmt: skip
    fail: list[str] = []
    if man is None:  # A1: the head job died before the tail wrote the manifest
        return ["no_manifest"]
    fail += [f"manifest_invalid:{p}" for p in manifest_errors(man)]
    failed_steps = [s.get("name") for s in man.get("steps") or [] if s.get("status") == "failed"]
    ok_head = job.get("state") == "COMPLETED" and job.get("exit_code") == "0:0"
    if man.get("status") != "completed" or man.get("exit_code") != 0 or failed_steps or not ok_head:
        why = man.get("failed_step") or (failed_steps[0] if failed_steps else None)
        if why is None:
            code = man.get("exit_code")
            why = code if code not in (0, None) else f"{job.get('state')} {job.get('exit_code')}"
        fail.append(f"run_failed:{why}")  # A2
    commit, want = str(man.get("git_commit") or ""), str(expect.get("commit") or "")
    if not commit or not want or not (commit.startswith(want) or want.startswith(commit)):
        fail.append(f"commit_mismatch:{commit}")  # A3
    if man.get("git_dirty") is not False:
        fail.append("dirty_checkout")
    ann = man.get("annotation") or {}
    bakta_task = any(_process(r.get("name", "")) == "BAKTA" for r in trace)
    if bakta_task or any(ann.get(k) != v for k, v in RECORDED_BAKTA.items()):
        fail.append("live_bakta")  # A4
    seen = {r.get("role"): r.get("sha256") for r in man.get("inputs") or []}
    want_in = expect.get("inputs") or {}
    fail += [f"input_mismatch:{r}" for r in CHECKED_INPUTS if seen.get(r) != want_in.get(r)]
    if tasks["max_alloc_cpus"] > MAX_CPUS:
        fail.append(f"cpu_cap:{tasks['max_alloc_cpus']}")  # A6
    if not _nonempty(out / "report" / "report.html"):
        fail.append("no_report")  # A7
    fail += _genes_failure(out / "tables" / "genes.tsv")  # A8
    fail += [f"module_off:{m}" for m in man.get("modules_off") or []]  # A9
    return fail


def _nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _genes_failure(path: Path) -> list[str]:
    from giap.io.output_checks import check_tsv

    if not path.is_file():
        return [f"genes_table:missing {path}"]
    try:
        check_tsv("genes", path, "genes")
    except ValueError as exc:
        return [f"genes_table:{exc}"]
    return []


def arm_evidence(
    outdir: Path, sacct_rows: Sequence[Mapping[str, str]], expect: Mapping[str, Any]
) -> dict[str, Any]:
    """Rule set A on one ``giap run`` output directory; ``failures`` empty when the arm passes.

    ``sacct_rows``: the head job's rows (``JobID == expect["job"]`` and its steps) and the task
    jobs' rows (``JobID`` in the trace ``native_id`` column). ``expect``: ``arm``, ``job``,
    ``commit`` and ``inputs`` (role → sha256).
    """
    out = Path(outdir)
    head_id = str(expect["job"])
    head_rows = _job_rows(sacct_rows, head_id)
    head = next((r for r in head_rows if r.get("JobID") == head_id), {})
    job = {"id": head_id, "state": head.get("State"), "exit_code": head.get("ExitCode"),
           "elapsed": head.get("Elapsed"), "max_rss": _max_rss(head_rows),
           "node": head.get("NodeList")}  # fmt: skip
    trace_path = out / "logs" / "trace.tsv"
    trace = read_rows(trace_path) if trace_path.is_file() else []
    tasks = _task_summary(trace, sacct_rows)
    mpath = out / "run_manifest.json"
    man: dict[str, Any] | None = None
    if mpath.is_file():
        loaded = json.loads(mpath.read_text(encoding="utf-8"))
        man = loaded if isinstance(loaded, dict) else {}
    report = out / "report" / "report.html"
    tables = sorted((out / "tables").glob("*.tsv")) if (out / "tables").is_dir() else []
    return {
        "arm": expect["arm"],
        "outdir": str(out),
        "job": job,
        "manifest": {k: (man or {}).get(k) for k in MANIFEST_KEYS},
        "inputs": {r["role"]: r.get("sha256") for r in (man or {}).get("inputs") or []},
        "report": file_entry(report) if report.is_file() else None,
        "tables": {p.name: file_entry(p) for p in tables},
        "tasks": tasks,
        "failures": _arm_failures(out, man, job, trace, tasks, expect),
    }


# ------------------------------------------------------------------ rule set C


def _stability_families() -> list[str]:
    res = resources.files("giap.schemas").joinpath("tables", "compare_stability.schema.json")
    fields = json.loads(res.read_text(encoding="utf-8"))["fields"]
    return list(next(f for f in fields if f["name"] == "family")["constraints"]["enum"])


def _compare_failures(cmp: Path, cmp2: Path, labels: Sequence[str]) -> list[str]:
    from giap.io.output_checks import check_tsv

    a, b = labels[0], labels[1]
    fail: list[str] = []
    rows: dict[str, list[dict[str, str]]] = {}
    for name in COMPARE_TABLES:  # C1
        path = cmp / "tables" / f"{name}.tsv"
        if not path.is_file():
            fail.append(f"table:{name}:missing")
            continue
        try:
            check_tsv(name, path, name)
        except ValueError as exc:
            fail.append(f"table:{name}:{exc}")
        rows[name] = read_rows(path)
    if "compare_abundance" in rows:  # C2
        ab = rows["compare_abundance"]
        if len(ab) != 1 or (ab[0].get("run_a"), ab[0].get("run_b")) != (a, b):
            fail.append(f"abundance_rows:{len(ab)}")
    if "compare_detection" in rows:  # C3: G = the intersection rows' n_genes (lead decision 1)
        det = rows["compare_detection"]
        inter = [r for r in det if r.get("section") == "intersection"]
        g = sum(int(r["n_genes"]) for r in inter if str(r.get("n_genes")).isdigit())
        matrix = [(r.get("locus_tag") or "", r.get("run") or "") for r in det
                  if r.get("section") == "matrix"]  # fmt: skip
        tags = {t for t, _ in matrix}
        per_run = all(sorted(t for t, r in matrix if r == lab) == sorted(tags) for lab in (a, b))
        if len(matrix) != 2 * g or len(set(matrix)) != len(matrix) or not per_run:
            fail.append(f"matrix_rows:{len(matrix)}")
    if "compare_stability" in rows:  # C4
        have = {r.get("family") for r in rows["compare_stability"]
                if (r.get("run_a"), r.get("run_b")) == (a, b)}  # fmt: skip
        fail += [f"stability_family:{f}" for f in _stability_families() if f not in have]
    if not all(_nonempty(cmp / p) for p in ("report/compare.html", "plots/V16_compare.svg")):
        fail.append("no_compare_report")  # C5
    for name in COMPARE_TABLES:  # C6
        p1, p2 = cmp / "tables" / f"{name}.tsv", cmp2 / "tables" / f"{name}.tsv"
        if p1.is_file() and (not p2.is_file() or p1.read_bytes() != p2.read_bytes()):
            fail.append(f"nondeterministic:{name}")
    return fail


def compare_evidence(cmp: Path, cmp2: Path, labels: Sequence[str]) -> dict[str, Any]:
    """Rule set C on a ``giap compare`` directory and its same-seed repeat ``cmp2``."""
    cmp, cmp2 = Path(cmp), Path(cmp2)
    tables = {n: file_entry(cmp / "tables" / f"{n}.tsv") for n in COMPARE_TABLES
              if (cmp / "tables" / f"{n}.tsv").is_file()}  # fmt: skip
    html = cmp / "report" / "compare.html"
    return {
        "cmp": str(cmp),
        "second": str(cmp2),
        "labels": list(labels),
        "tables": tables,
        "compare_html": file_entry(html) if html.is_file() else None,
        "failures": _compare_failures(cmp, cmp2, labels),
    }


# ------------------------------------------------------------------ CLI (I/O only here)


def _sacct(jobs: Sequence[str], fields: str) -> list[dict[str, str]]:
    names = fields.split(",")
    out: list[dict[str, str]] = []
    for i in range(0, len(jobs), 200):
        cmd = ["sacct", "-P", "-n", "-j", ",".join(jobs[i : i + 200]), "-o", fields]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        out += [dict(zip(names, ln.split("|"), strict=False))
                for ln in res.stdout.splitlines() if ln.strip()]  # fmt: skip
    return out


def _write(doc: Mapping[str, Any], dest: Path, row: str) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(row)
    for f in doc["failures"]:
        print(f"FAIL {f}", file=sys.stderr)
    return 0 if not doc["failures"] else 1


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    arm = sub.add_parser("arm")
    arm.add_argument("outdir", type=Path)
    arm.add_argument("--arm", required=True, choices=("A1", "A2"))
    arm.add_argument("--job", required=True)
    arm.add_argument("--expect-commit", required=True)
    arm.add_argument("--expect-input", action="append", default=[], metavar="ROLE=SHA")
    arm.add_argument("--out", required=True, type=Path)
    cmp = sub.add_parser("compare")
    cmp.add_argument("cmp", type=Path)
    cmp.add_argument("--second", required=True, type=Path)
    cmp.add_argument("--labels", default="A1,A2")
    cmp.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)
    try:
        if args.cmd == "arm":
            inputs = dict(s.split("=", 1) for s in args.expect_input)
            trace_path = args.outdir / "logs" / "trace.tsv"
            trace = read_rows(trace_path) if trace_path.is_file() else []
            ids = sorted({r["native_id"] for r in trace if r.get("native_id", "").isdigit()})
            rows = _sacct([args.job], HEAD_FIELDS) + (_sacct(ids, TASK_FIELDS) if ids else [])
            exp = {"arm": args.arm, "job": args.job, "commit": args.expect_commit,
                   "inputs": inputs}  # fmt: skip
            doc = arm_evidence(args.outdir, rows, exp)
            j, m, t = doc["job"], doc["manifest"], doc["tasks"]
            row = (
                f"| {doc['arm']} | {j['id']} | {j['state']} | {j['exit_code']} | "
                f"{j['elapsed']} | {j['max_rss']} | {j['node']} | `{m['git_commit']}` | "
                f"{t['max_alloc_cpus']} | {'pass' if not doc['failures'] else 'FAIL'} |"
            )
            return _write(doc, args.out, row)
        labels = tuple(args.labels.split(","))
        if len(labels) != 2:
            ap.error("--labels takes exactly two labels")
        doc = compare_evidence(args.cmp, args.second, labels)
        html = doc["compare_html"] or {}
        row = (
            f"| compare | {args.cmp.name} | `{html.get('sha256', 'NA')}` | "
            f"{html.get('bytes', 'NA')} bytes | {'pass' if not doc['failures'] else 'FAIL'} |"
        )
        return _write(doc, args.out, row)
    except (OSError, subprocess.CalledProcessError) as exc:
        print(f"run_evidence: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
