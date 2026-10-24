"""DIA-NN pilot helpers for the M10-02 Slurm jobs (docs/specs/M10-02-diann-searches.md).

``python -I validation/pilot/diann_report.py <command> …``:

* ``runs ACQUISITION_JSON``: the run labels of M10-01's acquisition file, one per line; exit 2 if
  the file is missing or malformed, 3 if its mode forbids the DIA-NN settings (§17, ADR-0051);
* ``provenance --arm A --out OUT …``: ``OUT/provenance.json`` of one search attempt (stdlib only);
* ``check REPORT --acquisition … --sums … --out check.json``: the report's five check rules;
* ``sample REPORT --out DIR [--n 1000] [--decoys 100]``: the stratified, hash-keyed sample.

Each command is a pure function plus a thin CLI; nothing is written into the repository. pyarrow
is imported only by ``check`` and ``sample``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import subprocess
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

PIN_NAME = "diann_2.3.1.sif"
REQUIRED = ("Stripped.Sequence", "Modified.Sequence", "Protein.Ids", "Protein.Group", "Run",
            "Precursor.Charge", "Decoy", "PG.MaxLFQ")  # fmt: skip
Q_COLUMNS = ("Global.Q.Value", "Q.Value")
VARIANT = re.compile(r"\b(synchro|slice|midia)[-_ ]?pasef\b", re.I)
SUMS_LINE = re.compile(r"^([0-9a-f]{64}) [ *](.+)$")


class AcquisitionError(Exception):
    """The acquisition file is missing or malformed (exit 2)."""


# ------------------------------------------------------------------ runs


def load_acquisition(path: Path) -> dict[str, Any]:
    try:
        doc = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise AcquisitionError(f"{path}: unreadable acquisition file ({exc})") from exc
    if not isinstance(doc, dict) or not isinstance(doc.get("runs"), list):
        raise AcquisitionError(f"{path}: not an object with a 'runs' list (M10-01 shape)")
    for r in doc["runs"]:
        if not isinstance(r, dict) or not isinstance(r.get("run"), str) or not r["run"]:
            raise AcquisitionError(f"{path}: a run entry has no 'run' label: {r!r}")
    return doc


def run_labels(doc: Mapping[str, Any]) -> list[str]:
    return [str(r["run"]) for r in doc["runs"]]


def acquisition_refusal(doc: Mapping[str, Any]) -> str | None:
    """Why these DIA-NN settings must not run on this acquisition (exit 3), or ``None``."""
    mode = doc.get("dataset_mode")
    if mode != "diaPASEF":
        return (f"§17: the acquisition mode is {mode}; the DIA-NN settings must change "
                "(ADR-0051 amendment, Khoa re-checks the flags)")  # fmt: skip
    for r in doc["runs"]:
        m = VARIANT.search(str(r.get("notes") or ""))
        if m:
            return (f"§17: run {r['run']} notes name a diaPASEF variant ({m.group(0)}); the "
                    "DIA-NN settings must be confirmed for its window scheme first (question "
                    "Q5, ADR-0051)")  # fmt: skip
    return None


# ------------------------------------------------------------------ sample selection


def row_key(run: object, precursor_id: object) -> str:
    return hashlib.sha256(f"{run}\t{precursor_id}".encode()).hexdigest()


def select_rows(
    run: Sequence[Any], decoy: Sequence[Any], precursor_id: Sequence[Any], n: int = 1000,
    d: int = 100,
) -> list[int]:  # fmt: skip
    """File-order indices of the stratified sample (spec "Row selection"), ascending.

    Strata are (run, decoy); ``n - d`` targets and ``d`` decoys are split evenly over the runs
    sorted by name, the remainder one each to the first runs. Each stratum takes its rows with
    the smallest ``(sha256(run\\tprecursor_id), index)``; a shortfall is filled by the smallest
    keys not yet taken over all strata.
    """
    keyed = sorted((row_key(run[i], precursor_id[i]), i) for i in range(len(run)))
    runs = sorted(set(run))
    buckets: dict[tuple[Any, int], list[int]] = {}
    for _, i in keyed:
        buckets.setdefault((run[i], int(bool(decoy[i]))), []).append(i)
    taken: set[int] = set()
    for kind, total in ((0, n - d), (1, d)):
        base, rem = divmod(total, len(runs)) if runs else (0, 0)
        for j, name in enumerate(runs):
            taken.update(buckets.get((name, kind), [])[: base + (1 if j < rem else 0)])
    for _, i in keyed:
        if len(taken) >= n:
            break
        taken.add(i)
    return sorted(taken)


# ------------------------------------------------------------------ check


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run_matches(value: str, label: str) -> bool:
    base = re.split(r"[/\\]", value)[-1]
    return value == label or base.removesuffix(".d") == label


def structural_failures(report: Path) -> tuple[list[str], Any]:
    """Rules 1 and 2: ``(failures, parquet metadata or None)``."""
    import pyarrow.parquet as pq

    if not report.is_file() or report.stat().st_size == 0:
        return [f"rule 1: {report} is missing or empty"], None
    try:
        meta = pq.read_metadata(report)
    except Exception as exc:
        return [f"rule 1: parquet footer unreadable ({type(exc).__name__}: {exc})"], None
    if meta.num_rows == 0:
        return ["rule 1: the report has no rows"], meta
    names = set(meta.schema.to_arrow_schema().names)
    fails = [f"rule 2: missing column {c}" for c in REQUIRED if c not in names]
    if not names & set(Q_COLUMNS):
        fails.append("rule 2: missing column Global.Q.Value or Q.Value")
    return fails, meta


def sums_failure(report: Path, sums: Path) -> str | None:
    """Rule 5: the report's sha256 equals its line in ``outputs.sha256``."""
    try:
        lines = sums.read_text().splitlines()
    except OSError:
        return f"rule 5: {sums} is unreadable, so the report's sha256 cannot be confirmed"
    listed = [m.group(1) for m in map(SUMS_LINE.match, lines)
              if m and Path(m.group(2)).name == report.name]  # fmt: skip
    if not listed:
        return f"rule 5: {report.name} has no sha256 line in {sums}"
    have = sha256_file(report) if report.is_file() else "missing"
    if listed[-1] != have:
        return f"rule 5: sha256 of {report.name} is {have}, {sums.name} records {listed[-1]}"
    return None


def max_q(table: Any) -> float | None:
    import pyarrow.compute as pc

    vals = [pc.max(table.column(c)).as_py() for c in Q_COLUMNS if c in table.column_names]
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None


def check_report(report: Path, acquisition: Path, sums: Path) -> dict[str, Any]:
    import pyarrow as pa
    import pyarrow.parquet as pq

    fails, meta = structural_failures(report)
    doc: dict[str, Any] = {"ok": False, "failures": fails, "rows": None, "columns": None,
                           "runs": None, "decoy_rows": None, "max_q": None,
                           "pyarrow": pa.__version__}  # fmt: skip
    if meta is not None:
        names = meta.schema.to_arrow_schema().names
        doc["rows"], doc["columns"] = meta.num_rows, names
        want = [c for c in ("Run", "Decoy", *Q_COLUMNS) if c in names]
        table = pq.read_table(report, columns=want)
        doc["max_q"] = max_q(table)
        if "Decoy" in names:
            doc["decoy_rows"] = sum(1 for v in table.column("Decoy").to_pylist() if v)
            if doc["decoy_rows"] == 0 and meta.num_rows:
                fails.append("rule 4: no decoy rows (Decoy = 1): --report-decoys did not "
                             "reach the output")  # fmt: skip
        if "Run" in names:
            counts: dict[str, int] = {}
            for v in table.column("Run").to_pylist():
                counts[str(v)] = counts.get(str(v), 0) + 1
            doc["runs"] = counts
            fails.extend(run_failures(counts, acquisition))
    fail5 = sums_failure(report, sums)
    if fail5:
        fails.append(fail5)
    doc["ok"] = not fails
    return doc


def run_failures(counts: Mapping[str, int], acquisition: Path) -> list[str]:
    """Rule 3: every acquisition run matches some ``Run`` value and no value is left over."""
    try:
        labels = run_labels(load_acquisition(acquisition))
    except AcquisitionError as exc:
        return [f"rule 3: {exc}"]
    fails = [f"rule 3: run {lab} of the acquisition file has no Run value in the report"
             for lab in labels if not any(run_matches(v, lab) for v in counts)]  # fmt: skip
    fails += [f"rule 3: Run value {v} matches no run of the acquisition file"
              for v in counts if not any(run_matches(v, lab) for lab in labels)]  # fmt: skip
    return fails


# ------------------------------------------------------------------ sample


def sample_report(report: Path, out: Path, n: int, d: int) -> tuple[list[str], dict[str, Any]]:
    """Write ``report.sample.parquet``, ``schema.tsv`` and ``summary.json`` under ``out``."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    fails, _ = structural_failures(report)
    if fails:
        return fails, {}
    table = pq.read_table(report)
    names = table.column_names
    run = table.column("Run").to_pylist()
    decoy = table.column("Decoy").to_pylist()
    if "Precursor.Id" in names:
        pid = table.column("Precursor.Id").to_pylist()
    else:
        pid = [f"{m}{c}" for m, c in zip(table.column("Modified.Sequence").to_pylist(),
                                          table.column("Precursor.Charge").to_pylist(),
                                          strict=True)]  # fmt: skip
    rows = select_rows(run, decoy, pid, n=n, d=d)
    keys = {row_key(r, p) for r, p in zip(run, pid, strict=True)}
    qcol = next(c for c in Q_COLUMNS if c in names)
    stripped = table.column("Stripped.Sequence").to_pylist()
    q = table.column(qcol).to_pylist()
    per_run: dict[str, int] = {}
    for r in run:
        per_run[str(r)] = per_run.get(str(r), 0) + 1
    summary = {
        "rows": table.num_rows,
        "columns": len(names),
        "rows_per_run": per_run,
        "decoy_rows": sum(1 for v in decoy if v),
        "distinct_precursors": len(set(pid)),
        "distinct_stripped_q01": len({s for s, dv, qv in zip(stripped, decoy, q, strict=True)
                                      if not dv and qv is not None and qv <= 0.01}),
        "max_q": max_q(table),
        "report_sha256": sha256_file(report),
        "report_bytes": report.stat().st_size,
        "sample_rows": len(rows),
        "sample_rule": f"strata(run,decoy) sha256(Run\tPrecursor.Id) n={n} d={d}",
        "duplicate_keys": table.num_rows - len(keys),
        "pyarrow": pa.__version__,
    }  # fmt: skip
    out.mkdir(parents=True, exist_ok=True)
    sample = table.take(pa.array(rows, pa.int64()))
    pq.write_table(sample, out / "report.sample.parquet", compression="zstd")
    tsv = ["column\tarrow_type\tnullable", *(f"{f.name}\t{f.type}\t{f.nullable}"
                                             for f in table.schema)]  # fmt: skip
    (out / "schema.tsv").write_text("\n".join(tsv) + "\n")
    write_json(out / "summary.json", summary)
    return [], summary


# ------------------------------------------------------------------ provenance


def file_record(path: Path, *, digest: bool = True) -> dict[str, Any]:
    rec: dict[str, Any] = {"path": str(path), "sha256": None}
    if path.is_file():
        st = path.stat()
        rec.update(size=st.st_size, mode=format(stat.S_IMODE(st.st_mode), "o"))
        if digest:
            rec["sha256"] = sha256_file(path)
    return rec


def git_state(repo: Path) -> dict[str, Any]:
    def git(*args: str) -> str | None:
        res = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True,
                             check=False)  # fmt: skip
        return res.stdout.strip() if res.returncode == 0 else None

    status = git("status", "--porcelain")
    return {"path": str(repo), "commit": git("rev-parse", "HEAD"),
            "dirty": None if status is None else bool(status)}  # fmt: skip


def cpu_model() -> str | None:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return None


def read_pin(pin_file: Path) -> str | None:
    try:
        lines = pin_file.read_text().splitlines()
    except OSError:
        return None
    m = re.fullmatch(r"([0-9a-f]{64})  " + re.escape(PIN_NAME), lines[0]) if lines else None
    return m.group(1) if m and len(lines) == 1 else None


def iso(t: float) -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t))


def provenance(a: argparse.Namespace, env: Mapping[str, str], now: float) -> dict[str, Any]:
    out = Path(a.out)
    steps = {s: {"status": getattr(a, s), "exit_code": _int(getattr(a, f"{s}_rc"))}
             for s in ("step1", "step2")}  # fmt: skip
    commands: dict[str, list[str] | None] = {}
    for s, n in (("step1", 1), ("step2", 2)):
        f = out / f"command_step{n}.txt"
        ran = steps[s]["status"] != "not_reached" and f.is_file()
        commands[s] = f.read_text().splitlines() if ran else None
    banner_file = out / "diann_version.txt"
    image: dict[str, Any] = file_record(Path(a.image)) if a.image else {"path": None}
    image.update(pin=read_pin(Path(a.pin_file)), pin_file=a.pin_file,
                 banner=[ln for ln in banner_file.read_text().splitlines() if ln.strip()][:2]
                 if banner_file.is_file() else None)  # fmt: skip
    image.setdefault("mode", None)
    acq: dict[str, Any] = file_record(Path(a.acquisition))
    runs: list[dict[str, Any]] = []
    try:
        doc = load_acquisition(Path(a.acquisition))
        acq["dataset_mode"] = doc.get("dataset_mode")
        for label in run_labels(doc):
            d = Path(a.raw_dir) / f"{label}.d"
            tdf_bin = d / "analysis.tdf_bin"
            runs.append({"label": label, "path": str(d),
                         "analysis_tdf_sha256": file_record(d / "analysis.tdf")["sha256"],
                         "analysis_tdf_bin_size": tdf_bin.stat().st_size
                         if tdf_bin.is_file() else None})  # fmt: skip
    except AcquisitionError:
        acq["dataset_mode"] = None
    fastas = [file_record(Path(p)) for p in dict.fromkeys(a.fastas.splitlines()) if p]
    sums = out / "outputs.sha256"
    return {
        "arm": a.arm,
        "status": a.status,
        "exit_code": int(a.exit_code),
        **steps,
        "job": {
            "id": env.get("SLURM_JOB_ID"),
            "name": env.get("SLURM_JOB_NAME"),
            "partition": env.get("SLURM_JOB_PARTITION"),
            "cpus": env.get("SLURM_CPUS_PER_TASK"),
            "mem": env.get("SLURM_MEM_PER_NODE"),
            "node": a.node,
            "cpu_model": cpu_model(),
            "start_utc": iso(a.start),
            "end_utc": iso(now),
            "elapsed_s": round(now - a.start, 3),
        },
        "image": image,
        "apptainer": {"version": a.apptainer or None},
        "cfg": file_record(Path(a.cfg)),
        "repo": git_state(Path(a.repo)),
        "inputs": {"fasta": fastas, "acquisition_json": acq, "runs": runs},
        "commands": commands,
        "outputs": str(sums) if a.status == "ok" and sums.is_file() else None,
    }  # fmt: skip


def _int(value: str | None) -> int | None:
    return int(value) if value not in (None, "") else None


# ------------------------------------------------------------------ CLI


def write_json(path: Path, doc: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n")
    tmp.replace(path)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="diann_report.py", description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("runs").add_argument("acquisition")
    c = sub.add_parser("check")
    c.add_argument("report")
    for flag in ("--acquisition", "--sums", "--out"):
        c.add_argument(flag, required=True)
    s = sub.add_parser("sample")
    s.add_argument("report")
    s.add_argument("--out", required=True)
    s.add_argument("--n", type=int, default=1000)
    s.add_argument("--decoys", type=int, default=100)
    v = sub.add_parser("provenance")
    for flag in ("--arm", "--out", "--status", "--exit-code"):
        v.add_argument(flag, required=True)
    for flag in ("--step1", "--step2"):
        v.add_argument(flag, default="not_reached", choices=("ran", "skipped", "not_reached"))
    for flag in ("--step1-rc", "--step2-rc", "--image", "--apptainer", "--node"):
        v.add_argument(flag, default="")
    for flag in ("--pin-file", "--cfg", "--repo", "--acquisition", "--raw-dir", "--fastas"):
        v.add_argument(flag, default="")
    v.add_argument("--start", type=float, default=time.time())
    return p


def main(argv: Sequence[str] | None = None) -> int:
    a = parser().parse_args(argv)
    if a.command == "runs":
        try:
            doc = load_acquisition(Path(a.acquisition))
        except AcquisitionError as exc:
            print(exc, file=sys.stderr)
            return 2
        refusal = acquisition_refusal(doc)
        if refusal:
            print(refusal, file=sys.stderr)
            return 3
        print("\n".join(run_labels(doc)))
        return 0
    if a.command == "check":
        result = check_report(Path(a.report), Path(a.acquisition), Path(a.sums))
        write_json(Path(a.out), result)
        for f in result["failures"]:
            print(f"check: {f}", file=sys.stderr)
        return 0 if result["ok"] else 1
    if a.command == "sample":
        if not 0 <= a.decoys <= a.n:
            print("sample: need 0 <= --decoys <= --n", file=sys.stderr)
            return 2
        fails, summary = sample_report(Path(a.report), Path(a.out), a.n, a.decoys)
        for f in fails:
            print(f"sample: {f}", file=sys.stderr)
        if not fails:
            print(f"sample: {summary['sample_rows']} of {summary['rows']} rows -> {a.out}")
        return 1 if fails else 0
    write_json(Path(a.out) / "provenance.json", provenance(a, os.environ, time.time()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
