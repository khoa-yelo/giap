"""Evidence for the pass-1 pilot runs: rule sets A (one arm) and C (the compare) of M10-04.

Spec: docs/specs/M10-04-pilot-runs.md ("run_evidence.py", lead decisions of 2026-10-24). The
evidence holds provenance only (paths, sha256s, byte sizes, exit codes, states, resources), never
a value from the tables (Khoa's 2026-10-23 rule).

Usage (``PYTHONPATH=$GIAP_REPO/src``; exit 0 when ``failures`` is empty, else 1)::

    python validation/pilot/run_evidence.py arm OUT --arm A1|A2 --job JOBID --expect-commit SHA \
        --expect-input genome=SHA --expect-input annotation=SHA --expect-input peptides=SHA \
        --out FILE
    python validation/pilot/run_evidence.py compare CMP --second CMP2 --labels A1,A2 \
        --genes OUT_A1/tables/genes.tsv --out FILE
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
# rule A5: every --expect-input role is checked; genome and peptides always (review S1)
CHECKED_INPUTS = ("genome", "peptides")
RECORDED_BAKTA = {"bakta_version": "1.12.2", "db_type": "full"}  # rule A4
COMPARE_TABLES = ("compare_detection", "compare_abundance", "compare_stability")
MANIFEST_KEYS = (
    "status", "exit_code", "failed_step", "failed_log", "git_commit", "git_dirty", "profile",
    "resume", "started_at", "finished_at", "seeds", "modules_off", "detection_tier",
    "annotation", "databases", "warnings", "skipped", "layers_missing",
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


def _process(name: str) -> str:
    """``GIAP:ANNOTATION:BAKTA (genome)`` → ``BAKTA``."""
    return name.split(" ", 1)[0].rsplit(":", 1)[-1]


def native_ids(trace: Sequence[Mapping[str, str]]) -> list[str]:
    """The trace's Slurm job IDs (one filter for the CLI's ``sacct`` query and rule A6)."""
    return sorted({j for r in trace if (j := r.get("native_id") or "-") != "-"})


def manifest_errors(doc: Any) -> list[str]:
    """JSON pointers of the schema violations of ``doc``."""
    import jsonschema  # type: ignore[import-untyped]

    res = resources.files("giap.schemas").joinpath("run_manifest.schema.json")
    schema = json.loads(res.read_text(encoding="utf-8"))
    errors = jsonschema.Draft202012Validator(schema).iter_errors(doc)
    return sorted({"/" + "/".join(str(p) for p in e.absolute_path) for e in errors})


# ------------------------------------------------------------------ rule set A


def _task_summary(
    trace: Sequence[Mapping[str, str]], sacct: Sequence[Mapping[str, str]]
) -> dict[str, Any]:
    rows = {jid: _job_rows(sacct, jid) for jid in native_ids(trace)}
    cpus = {jid: [int(r["AllocCPUS"]) for r in rs if str(r.get("AllocCPUS", "")).isdigit()]
            for jid, rs in rows.items()}  # fmt: skip
    rss = [b for rs in rows.values() for r in rs if (b := rss_bytes(r.get("MaxRSS"))) is not None]
    return {
        "count": len(trace),
        "failed": sum(r.get("status") in ("FAILED", "ABORTED") for r in trace),
        "retried": sum(int(a) > 1 for r in trace if (a := r.get("attempt") or "").isdigit()),
        "max_alloc_cpus": max((c for cs in cpus.values() for c in cs), default=0),
        # A6 cannot be judged for tasks without a sacct allocation row (review B1)
        "unknown_cpus": sum(not cpus.get(r.get("native_id") or "-") for r in trace),
        "peak_max_rss_bytes": max(rss, default=None),
    }


def _arm_failures(
    out: Path, man: Mapping[str, Any] | str | None, job: Mapping[str, Any],
    trace: Sequence[Mapping[str, str]], tasks: Mapping[str, Any], expect: Mapping[str, Any],
) -> list[str]:  # fmt: skip
    fail: list[str] = []
    if man is None:  # A1: the head job died before the tail wrote the manifest
        return ["no_manifest"]
    if isinstance(man, str):  # A1: a truncated or unreadable manifest (review S2)
        return ["manifest_invalid:/"]
    fail += [f"manifest_invalid:{p}" for p in manifest_errors(man)]
    failed_steps = [s.get("name") for s in man.get("steps") or [] if s.get("status") == "failed"]
    ok_head = job.get("state") == "COMPLETED" and job.get("exit_code") == "0:0"
    if man.get("status") != "completed" or man.get("exit_code") != 0 or failed_steps or not ok_head:
        why = man.get("failed_step") or (failed_steps[0] if failed_steps else None)
        if why is None:
            code, state = man.get("exit_code"), job.get("state")
            why = code if code not in (0, None) else "no_sacct" if state is None else state
        fail.append(f"run_failed:{why}")  # A2
    commit, want = str(man.get("git_commit") or ""), str(expect.get("commit") or "")
    if not want or not commit.startswith(want):
        fail.append(f"commit_mismatch:{commit}")  # A3
    if man.get("git_dirty") is not False:
        fail.append("dirty_checkout")
    ann = man.get("annotation") or {}
    bakta_task = any(_process(r.get("name", "")) == "BAKTA" for r in trace)
    if bakta_task or any(ann.get(k) != v for k, v in RECORDED_BAKTA.items()):
        fail.append("live_bakta")  # A4
    seen = {r.get("role"): r.get("sha256") for r in man.get("inputs") or []}
    want_in = expect.get("inputs") or {}
    for role in dict.fromkeys([*CHECKED_INPUTS, *want_in]):  # A5
        if role not in seen:
            fail.append(f"input_missing:{role}")
        elif seen[role] != want_in.get(role):
            fail.append(f"input_mismatch:{role}")
    if tasks["unknown_cpus"]:
        fail.append(f"cpu_cap:unknown:{tasks['unknown_cpus']}")  # A6 (review B1)
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
           "elapsed": head.get("Elapsed"),
           "max_rss": max((r.get("MaxRSS", "") for r in head_rows),
                          key=lambda s: rss_bytes(s) or -1, default=""),
           "node": head.get("NodeList")}  # fmt: skip
    trace_path = out / "logs" / "trace.tsv"
    trace = read_rows(trace_path) if trace_path.is_file() else []
    tasks = _task_summary(trace, sacct_rows)
    mpath = out / "run_manifest.json"
    man: dict[str, Any] | str | None = None
    if mpath.is_file():
        try:
            loaded = json.loads(mpath.read_text(encoding="utf-8"))
            man = loaded if isinstance(loaded, dict) else {}
        except ValueError:  # JSONDecodeError, UnicodeDecodeError: a truncated manifest
            man = "unparseable"
    md = man if isinstance(man, dict) else {}
    report = out / "report" / "report.html"
    tables = sorted((out / "tables").glob("*.tsv")) if (out / "tables").is_dir() else []
    return {
        "arm": expect["arm"],
        "outdir": str(out),
        "job": job,
        "manifest": {k: md.get(k) for k in MANIFEST_KEYS},
        "inputs": {r["role"]: r.get("sha256") for r in md.get("inputs") or []},
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


def _compare_failures(cmp: Path, cmp2: Path, labels: Sequence[str], genes: Path) -> list[str]:
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
    tags = {r.get("locus_tag") for r in read_rows(genes)} if genes.is_file() else None
    if tags is None:
        fail.append(f"genes_table:missing {genes}")
    elif "compare_detection" in rows:  # C3: one row per (genes.tsv locus tag, run) (review B4)
        matrix = [(r.get("locus_tag"), r.get("run")) for r in rows["compare_detection"]
                  if r.get("section") == "matrix"]  # fmt: skip
        want = {(t, lab) for t in tags for lab in (a, b)}
        if set(matrix) != want or len(matrix) != len(want):
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


def compare_evidence(
    cmp: Path, cmp2: Path, labels: Sequence[str], *, genes: Path
) -> dict[str, Any]:
    """Rule set C on a ``giap compare`` directory and its same-seed repeat ``cmp2``; ``genes`` is
    ``out_a1/tables/genes.tsv`` (G of rule C3)."""
    cmp, cmp2, genes = Path(cmp), Path(cmp2), Path(genes)
    tables = {n: file_entry(cmp / "tables" / f"{n}.tsv") for n in COMPARE_TABLES
              if (cmp / "tables" / f"{n}.tsv").is_file()}  # fmt: skip
    html = cmp / "report" / "compare.html"
    return {
        "cmp": str(cmp),
        "second": str(cmp2),
        "labels": list(labels),
        "tables": tables,
        "compare_html": file_entry(html) if html.is_file() else None,
        "genes": file_entry(genes) if genes.is_file() else None,
        "failures": _compare_failures(cmp, cmp2, labels, genes),
    }


# ------------------------------------------------------------------ CLI (I/O only here)


def _sacct(jobs: Sequence[str], fields: str) -> list[dict[str, str]]:
    names = fields.split(",")
    out: list[dict[str, str]] = []
    for i in range(0, len(jobs), 200):
        cmd = ["sacct", "-P", "-n", "-j", ",".join(jobs[i : i + 200]), "-o", fields]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for vals in (ln.split("|") for ln in res.stdout.splitlines() if ln.strip()):
            if len(vals) != len(names):  # exit 2, no evidence (review N3)
                raise ValueError(f"sacct printed {len(vals)} fields, not {len(names)}: {vals}")
            out.append(dict(zip(names, vals, strict=True)))
    return out


def _write(doc: Mapping[str, Any], dest: Path, row: str) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(row)
    for f in doc["failures"]:
        print(f"FAIL {f}", file=sys.stderr)
    return 0 if not doc["failures"] else 1


def _commit(text: str) -> str:
    if not re.fullmatch(r"[0-9a-f]{7,40}", text):
        raise argparse.ArgumentTypeError(f"not a 7-40 character hex sha: {text!r}")
    return text


def _role_sha(text: str) -> tuple[str, str]:
    role, eq, sha = text.partition("=")
    if not eq or not role:
        raise argparse.ArgumentTypeError(f"expected ROLE=SHA: {text!r}")
    return role, sha


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    arm = sub.add_parser("arm")
    arm.add_argument("outdir", type=Path)
    arm.add_argument("--arm", required=True, choices=("A1", "A2"))
    arm.add_argument("--job", required=True)
    arm.add_argument("--expect-commit", required=True, type=_commit)
    arm.add_argument("--expect-input", action="append", default=[], metavar="ROLE=SHA",
                     type=_role_sha)  # fmt: skip
    arm.add_argument("--out", required=True, type=Path)
    cmp = sub.add_parser("compare")
    cmp.add_argument("cmp", type=Path)
    cmp.add_argument("--second", required=True, type=Path)
    cmp.add_argument("--labels", default="A1,A2")
    cmp.add_argument("--genes", required=True, type=Path, help="OUT_A1/tables/genes.tsv")
    cmp.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)
    try:
        if args.cmd == "arm":
            inputs = dict(args.expect_input)
            trace_path = args.outdir / "logs" / "trace.tsv"
            trace = read_rows(trace_path) if trace_path.is_file() else []
            ids = native_ids(trace)
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
        doc = compare_evidence(args.cmp, args.second, labels, genes=args.genes)
        html = doc["compare_html"] or {}
        row = (
            f"| compare | {args.cmp.name} | `{html.get('sha256', 'NA')}` | "
            f"{html.get('bytes', 'NA')} bytes | {'pass' if not doc['failures'] else 'FAIL'} |"
        )
        return _write(doc, args.out, row)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"run_evidence: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
