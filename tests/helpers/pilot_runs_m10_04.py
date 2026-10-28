"""Hand-made ``gippa run`` output trees, ``sacct`` rows and a job sandbox for the M10-04 tests.

Spec: docs/specs/M10-04-pilot-runs.md (rule sets A and C, acceptance criteria 1-7, lead decisions
of 2026-10-24). Nothing here runs a pipeline, touches the network or reads the real
``$GIPPA_ROOT``: every value is invented (Khoa's 2026-10-23 rule: nothing derived from the DIA-NN
reports is committed).

* ``write_arm(out, ...)`` builds an ``OUT/`` tree as ``gippa run`` leaves it: a
  ``run_manifest.json`` produced by the real ``gippa.manifest.build_manifest`` (Nextflow executor,
  the trace rows below, the recorded-annotation metadata), ``logs/trace.tsv`` with the
  ``nextflow.config`` trace fields and Slurm ``native_id``s, ``report/report.html`` and a
  ``tables/genes.tsv`` written through the ``genes`` schema (``helpers.compare_data``). The same
  tree is a valid ``gippa compare`` input (the genome digest is the manifest's ``genome`` row).
* ``sacct_rows(...)`` is what ``sacct -P -n`` prints for the head job (``JobID``, ``State``,
  ``ExitCode``, ``Elapsed``, ``MaxRSS``, ``NodeList``; plus its ``.batch`` step) and for each
  task job (``JobID`` = the trace ``native_id``, ``AllocCPUS``, ``MaxRSS``, ``State``).
* ``arm_doc``, ``run_evidence_cli`` and ``arm_cli_args`` call ``arm_evidence`` and the CLI (the
  CLI with a fake ``sacct`` serving ``sacct_rows``; ``sacct_script`` replaces the fake).
* ``compare_job_sandbox`` builds on ``helpers.slurm_sandbox.make_sandbox`` and adds ``src/`` and
  ``validation/`` to the sandbox checkout (committed), so ``PYTHONPATH=$GIPPA_REPO/src`` runs this
  branch's ``gippa compare`` and ``repo_is_clean`` means "the job wrote nothing into the checkout".

Interface assumed (test-author reading of the spec; deliverable 1):

``arm_evidence(outdir, sacct_rows, expect) -> dict``
    ``outdir``: the run directory (``Path``). ``sacct_rows``: a list of dicts keyed by ``sacct``
    field names, holding the head job's row (``JobID == expect["job"]``), possibly its steps
    (``<job>.batch``, ignored) and the task jobs' rows (``JobID`` in the trace ``native_id``
    column). ``expect``: ``{"arm": "A1"|"A2", "job": str, "commit": str, "inputs": {role: sha}}``
    (the CLI's ``--arm``, ``--job``, ``--expect-commit`` and ``--expect-input role=sha``).
``compare_evidence(cmp, cmp2, labels, genes=PATH) -> dict``
    ``cmp`` and ``cmp2``: the two ``gippa compare`` output directories; ``labels``: ``("A1",
    "A2")``;
    ``genes``: ``out_a1/tables/genes.tsv`` (the CLI's required ``--genes``; lead decision on B4).

Both return a JSON-serialisable dict with a ``failures`` list of the spec's exact rule strings.
The ``sacct`` rows' head ``MaxRSS``/``Elapsed`` and the task rows' ``AllocCPUS`` are strings, as
``sacct -P`` prints them.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

from helpers.compare_data import gene, write_genes, write_table
from helpers.container_guard import guarded_path
from helpers.slurm_sandbox import REPO, Sandbox, make_sandbox, write_exe

RUN_EVIDENCE = REPO / "validation" / "pilot" / "run_evidence.py"
COMPARE_JOB = REPO / "scripts" / "slurm" / "pilot_compare.sbatch"
README = REPO / "validation" / "pilot" / "README.md"

# the minimal PATH of the `gippa compare` subprocess (container_guard.guarded_path keeps the
# container pull guard on it)
MIN_PATH = "/usr/bin:/bin"
HEAD_JOB = "1900001"
TASK_JOBS = ("1900002", "1900003")
COMMIT = "0123abc"  # the short sha in ATTEMPTS.tsv and --expect-commit
FULL_COMMIT = COMMIT + "4567890abcdef0123456789abcdef0123"  # the manifest's `git rev-parse HEAD`
SHA_GENOME = "1" * 64
SHA_ANNOTATION = "2" * 64
SHA_PEPTIDES = "3" * 64
EXPECT_INPUTS = {"genome": SHA_GENOME, "annotation": SHA_ANNOTATION, "peptides": SHA_PEPTIDES}
SEED = 20261014
N_GENES = 12
TRACE_FIELDS = [
    "task_id", "hash", "native_id", "name", "tag", "status", "exit", "attempt", "submit",
    "duration", "realtime", "%cpu", "peak_rss", "peak_vmem", "rchar", "wchar", "workdir",
    "hostname", "container",
]  # fmt: skip
ANNOTATION_META = {"annotation_source": "bakta", "bakta_version": "1.12.2", "db_version": "6.0",
                   "db_type": "full"}  # fmt: skip
TABLES = ("compare_detection", "compare_abundance", "compare_stability")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------ run directories


def trace_task(
    name: str, native_id: str, *, status: str = "COMPLETED", attempt: int = 1, task_id: int = 1
) -> dict[str, str]:
    row = dict.fromkeys(TRACE_FIELDS, "-")
    row.update({"task_id": str(task_id), "hash": f"ab/{task_id:06d}", "native_id": native_id,
                "name": name, "status": status, "exit": "0" if status == "COMPLETED" else "1",
                "attempt": str(attempt), "realtime": "2m 3s", "%cpu": "95.0%",
                "peak_rss": "1.2 GB", "workdir": f"/w/ab/{task_id:06d}"})  # fmt: skip
    return row


def default_trace() -> list[dict[str, str]]:
    """Two task jobs (criterion 1): gene table and detection, with Slurm ``native_id``s."""
    return [
        trace_task("GIPPA:CORE:GENE_TABLE (gene_table)", TASK_JOBS[0], task_id=1),
        trace_task("GIPPA:CORE:DETECTION (detection)", TASK_JOBS[1], task_id=2),
    ]


def bakta_task(native_id: str = "1900009") -> dict[str, str]:
    """The task a live Bakta annotation leaves in the trace (rule A4)."""
    return trace_task("GIPPA:ANNOTATION:BAKTA (genome)", native_id, task_id=9)


def arm_genes(*, shift: float = 0.0, detected: bool = True) -> list[dict[str, object]]:
    """``N_GENES`` invented genes; ``detected=False`` gives the empty-result arm (criterion 3)."""
    rows = []
    for i in range(N_GENES):
        log2 = (10.0 + ((i * 7) % 5) + shift * (i % 3)) if detected else None
        rows.append(gene(f"TST_{i:05d}", start=1000 * i, detected=detected and i < N_GENES - 1,
                         log2=log2, tu=f"TU_{i // 2}"))  # fmt: skip
    return rows


def manifest_doc(
    *,
    outdir: Path,
    trace: Sequence[Mapping[str, str]],
    commit: str | None = FULL_COMMIT,
    dirty: bool | None = False,
    exit_code: int = 0,
    failed_step: str | None = None,
    modules: Sequence[str] = ("core", "mobile", "bgc", "function"),
    annotation_meta: Mapping[str, Any] | None = None,
    inputs: Mapping[str, str] | None = None,
    events: Sequence[Mapping[str, Any]] = (),
    tier: str = "peptides",
) -> dict[str, Any]:
    """A manifest exactly as ``gippa.manifest.build_manifest`` writes it for a Nextflow run."""
    from gippa.config import GippaConfig
    from gippa.manifest import ManifestSources, build_manifest
    from helpers.manifest_data import preflight_doc

    log = None if failed_step is None else str(outdir / "work" / "ab" / "x" / ".command.log")
    shas = {**EXPECT_INPUTS, **(inputs or {})}
    files = {"genome": "nissle1917.fna.gz", "annotation": "nissle1917.bakta.json",
             "peptides": "a2.tsv.gz"}  # fmt: skip
    rows = [{"role": r, "file": files[r], "sha256": shas[r], "size": 1000 + i}
            for i, r in enumerate(("genome", "annotation", "peptides"))]  # fmt: skip
    src = ManifestSources(
        cfg=GippaConfig(seed=SEED, modules=list(modules)),
        preflight=preflight_doc(inputs=rows, tier=tier),
        registered=[],
        executor="nextflow",
        profile="farmshare,apptainer",
        resume=False,
        outdir=str(outdir),
        command_line=["gippa", "run", "--outdir", str(outdir)],
        started_at="2026-10-24T08:00:00+00:00",
        finished_at="2026-10-24T14:00:00+00:00",
        exit_code=exit_code,
        failed_step=failed_step,
        failed_log=log,
        git_commit=commit,
        git_dirty=dirty,
        input_paths={r["role"]: f"/pilot/inputs/{r['file']}" for r in rows},
        trace=list(trace),
        events=list(events),
        annotation_meta=dict(ANNOTATION_META if annotation_meta is None else annotation_meta),
    )
    doc = build_manifest(src)
    assert isinstance(doc, dict)
    return doc


def write_trace(out: Path, rows: Sequence[Mapping[str, str]]) -> Path:
    path = out / "logs" / "trace.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=TRACE_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "-") for k in TRACE_FIELDS})
    return path


def write_arm(
    out: Path,
    *,
    genes: Sequence[Mapping[str, object]] | None = None,
    trace: Sequence[Mapping[str, str]] | None = None,
    report: str | None = "<!doctype html><html><body>GIPPA report</body></html>\n",
    manifest: bool = True,
    novel_orfs_empty: bool = False,
    **manifest_kw: Any,
) -> Path:
    """An ``OUT/`` tree of one arm; the defaults pass every rule of set A."""
    out.mkdir(parents=True, exist_ok=True)
    rows = list(default_trace() if trace is None else trace)
    write_trace(out, rows)
    write_genes(out, arm_genes() if genes is None else genes)
    if novel_orfs_empty:
        write_table(out, "novel_orfs", [])
    if report is not None:
        (out / "report").mkdir(parents=True, exist_ok=True)
        (out / "report" / "report.html").write_text(report, encoding="utf-8")
    if manifest:
        doc = manifest_doc(outdir=out, trace=rows, **manifest_kw)
        (out / "run_manifest.json").write_text(
            json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8"
        )
    return out


def edit_manifest(out: Path, **changes: Any) -> None:
    """Overwrite top-level manifest keys (``None`` deletes the key) after the fact."""
    path = out / "run_manifest.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    for k, v in changes.items():
        if v is None:
            doc.pop(k, None)
        else:
            doc[k] = v
    path.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def sacct_rows(
    *,
    job: str = HEAD_JOB,
    state: str = "COMPLETED",
    exit_code: str = "0:0",
    tasks: Mapping[str, int] | None = None,
) -> list[dict[str, str]]:
    """Head job (+ ``.batch`` step) and task-job rows as ``sacct -P`` prints them."""
    alloc = dict.fromkeys(TASK_JOBS, 14) if tasks is None else dict(tasks)
    rows = [
        {"JobID": job, "State": state, "ExitCode": exit_code, "Elapsed": "05:59:12",
         "MaxRSS": "", "NodeList": "wheat-07", "AllocCPUS": "2"},
        {"JobID": f"{job}.batch", "State": state, "ExitCode": exit_code, "Elapsed": "05:59:12",
         "MaxRSS": "2104832K", "NodeList": "wheat-07", "AllocCPUS": "2"},
    ]  # fmt: skip
    for jid, cpus in alloc.items():
        rows.append({"JobID": jid, "State": "COMPLETED", "ExitCode": "0:0",
                     "Elapsed": "00:02:03", "MaxRSS": "1258291K", "NodeList": "oat-03",
                     "AllocCPUS": str(cpus)})  # fmt: skip
    return rows


def expect(**over: Any) -> dict[str, Any]:
    base: dict[str, Any] = {"arm": "A2", "job": HEAD_JOB, "commit": COMMIT,
                            "inputs": dict(EXPECT_INPUTS)}  # fmt: skip
    base.update(over)
    return base


def load_run_evidence() -> ModuleType:
    """Import ``validation/pilot/run_evidence.py`` (deliverable 1)."""
    assert RUN_EVIDENCE.is_file(), f"{RUN_EVIDENCE.relative_to(REPO)} does not exist (M10-04)"
    spec = importlib.util.spec_from_file_location("gippa_pilot_run_evidence", RUN_EVIDENCE)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def arm_doc(
    out: Path, sacct: Sequence[Mapping[str, str]] | None = None, **exp: Any
) -> dict[str, Any]:
    """``arm_evidence`` on ``out`` (default ``sacct_rows()``, ``expect(**exp)``), JSON-checked."""
    mod = load_run_evidence()
    doc = mod.arm_evidence(out, list(sacct_rows() if sacct is None else sacct), expect(**exp))
    assert isinstance(doc, dict)
    json.dumps(doc)  # JSON-serialisable
    return doc


def fake_sacct_env(
    tmp: Path, rows: Sequence[Mapping[str, str]], *, sacct_script: str | None = None
) -> dict[str, str]:
    """A minimal environment whose ``sacct`` is the fake (or ``sacct_script``) serving ``rows``."""
    from helpers.slurm_sandbox import FAKE_SACCT

    bindir, state = tmp / "fakebin", tmp / "fakestate"
    bindir.mkdir(exist_ok=True)
    state.mkdir(exist_ok=True)
    text = FAKE_SACCT if sacct_script is None else sacct_script
    write_exe(bindir / "sacct", text.replace("{python}", sys.executable))
    (state / "sacct_rows.json").write_text(json.dumps(list(rows)))
    return {"PATH": f"{bindir}:/usr/bin:/bin", "HOME": str(tmp), "LANG": "C.UTF-8",
            "FAKE_SLURM_DIR": str(state), "PYTHONPATH": str(REPO / "src"),
            "GIPPA_ROOT": str(tmp / "no_root"), "PYTHONDONTWRITEBYTECODE": "1"}  # fmt: skip


def run_evidence_cli(
    tmp: Path,
    *args: str,
    rows: Sequence[Mapping[str, str]] | None = None,
    sacct_script: str | None = None,
) -> subprocess.CompletedProcess[str]:
    """``python validation/pilot/run_evidence.py ARGS`` with the fake ``sacct``."""
    assert RUN_EVIDENCE.is_file(), "validation/pilot/run_evidence.py does not exist (M10-04)"
    env = fake_sacct_env(tmp, sacct_rows() if rows is None else rows, sacct_script=sacct_script)
    return subprocess.run([sys.executable, str(RUN_EVIDENCE), *args], capture_output=True,
                          text=True, env=env, cwd=str(tmp), timeout=300,
                          check=False)  # fmt: skip


def arm_cli_args(out: Path, dest: Path, **kw: str) -> list[str]:
    """``arm OUT --arm --job --expect-commit --expect-input (3) --out`` (``kw`` overrides)."""
    return ["arm", str(out), "--arm", kw.get("arm", "A2"), "--job", kw.get("job", HEAD_JOB),
            "--expect-commit", kw.get("commit", COMMIT),
            "--expect-input", f"genome={kw.get('genome', SHA_GENOME)}",
            "--expect-input", f"annotation={kw.get('annotation', SHA_ANNOTATION)}",
            "--expect-input", f"peptides={kw.get('peptides', SHA_PEPTIDES)}",
            "--out", str(dest)]  # fmt: skip


# ------------------------------------------------------------------ gippa compare


def gippa_compare(
    a: Path, b: Path, cmp: Path, *, env_root: Path, bootstrap: int = 200
) -> subprocess.CompletedProcess[str]:
    """The real ``gippa compare a b --outdir cmp --labels A1,A2`` of this checkout."""
    env_root.mkdir(parents=True, exist_ok=True)
    env = {"PATH": guarded_path(MIN_PATH), "HOME": str(env_root), "LANG": "C.UTF-8",
           "PYTHONPATH": str(REPO / "src"), "GIPPA_DB_DIR": str(env_root),
           "GIPPA_ROOT": str(env_root), "PYTHONDONTWRITEBYTECODE": "1",
           "MPLCONFIGDIR": str(env_root / "mpl")}  # fmt: skip
    cmd = [sys.executable, "-c", "from gippa.cli import cli; cli()", "compare", str(a), str(b),
           "--outdir", str(cmp), "--labels", "A1,A2", "--set",
           f"stats.bootstrap={bootstrap}"]  # fmt: skip
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300, check=False)


def two_arms(tmp: Path, *, detected: bool = True, sha_b: str | None = None) -> tuple[Path, Path]:
    """``out_a1`` and ``out_a2`` under ``tmp`` (``sha_b``: another genome digest for A2)."""
    a1 = write_arm(tmp / "out_a1", genes=arm_genes(shift=0.0, detected=detected))
    a2 = write_arm(
        tmp / "out_a2",
        genes=arm_genes(shift=0.5, detected=detected),
        inputs=None if sha_b is None else {"genome": sha_b},
    )
    return a1, a2


# ------------------------------------------------------------------ the compare job sandbox

_RECORDING_PYTHON = r"""#!/bin/bash
# Records every call (argv as one JSON line), then runs the real interpreter.
{python} - "$@" <<'PY' >> "$FAKE_SLURM_DIR/gippa_python.jsonl"
import json, sys
print(json.dumps(sys.argv[1:]))
PY
exec {python} "$@"
"""


@dataclass
class CompareJobSandbox:
    sb: Sandbox
    runs: Path

    @property
    def script(self) -> Path:
        return self.sb.repo / "scripts" / "slurm" / "pilot_compare.sbatch"

    def run(
        self, *args: str, env: dict[str, str] | None = None
    ) -> subprocess.CompletedProcess[str]:
        return self.sb.run(["bash", str(self.script), *args], env=env, timeout=600)

    def python_calls(self) -> list[list[str]]:
        return self.sb.jsonl("gippa_python.jsonl")

    def recording_python(self) -> Path:
        return self.sb.bin / "gippa_python"

    def repo_is_clean(self) -> bool:
        res = subprocess.run(
            ["git", "-C", str(self.sb.repo), "status", "--porcelain", "--ignored"],
            capture_output=True, text=True, check=True, env=self.sb.env,
        )  # fmt: skip
        return res.stdout.strip() == ""


def compare_job_sandbox(tmp: Path) -> CompareJobSandbox:
    sb = make_sandbox(tmp)
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    shutil.copytree(REPO / "src", sb.repo / "src", ignore=ignore)
    shutil.copytree(REPO / "validation", sb.repo / "validation", ignore=ignore)
    for args in (["add", "-A"], ["commit", "-q", "-m", "M10-04 sandbox: src, validation"]):
        subprocess.run(
            ["git", "-C", str(sb.repo), *args], check=True, capture_output=True, env=sb.env
        )
    write_exe(sb.bin / "gippa_python", _RECORDING_PYTHON.replace("{python}", sys.executable))
    sb.env.update({"PYTHONDONTWRITEBYTECODE": "1", "MPLCONFIGDIR": str(tmp / "mpl"),
                   "SLURM_JOB_ID": "1900100"})  # fmt: skip
    sb.env.pop("GIPPA_PYTHON", None)
    sb.env.pop("PYTHONPATH", None)
    runs = tmp / "runs"
    runs.mkdir()
    return CompareJobSandbox(sb=sb, runs=runs)
