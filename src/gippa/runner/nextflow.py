"""Nextflow launcher (GIPPA_SPEC.md §6.3 "Launch", "Error strategy"; M2-06).

Python validates everything and writes ``OUT/work/params.json`` plus ``OUT/resolved_config.yaml``;
Nextflow only reads the params file. A failed run maps to exit code 4 naming the failed step and the
``.command.log`` of its task directory (read from ``OUT/logs/trace.tsv``). All runner messages go to
stderr. After Nextflow ends, the step events each ``GIPPA_STEP`` task wrote (published as
``OUT/logs/steps/<step>.events.jsonl``) are appended to ``OUT/logs/events.jsonl`` in trace order
(M3-17 lead decision 15).

Resources: the user's ``--threads``/``--memory`` (``params.threads``/``params.memory``) and the
profile's limits combine as the minimum of the two, never a silent override: e.g. under
``farmshare`` a task gets at most ``min(16, threads)`` CPUs and ``min(64 GB, memory)``, under
``farmshare_local`` at most ``min(SLURM_CPUS_PER_TASK, threads)`` CPUs and
``min(SLURM_MEM_PER_NODE - 2 GB, memory)``.
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import yaml

from gippa.config import GippaConfig
from gippa.exitcodes import ExitCode
from gippa.runner.workflow_path import workflow_main as workflow_main  # re-export
from gippa.steps import load_step_modules
from gippa.steps import registry as reg

DEFAULT_MEMORY = "64 GB"  # §5.1 `gippa run --memory 64.GB`, written in a Nextflow-parseable form
_FAILED = {"FAILED", "ABORTED"}
_EMPTY = {"", "-"}
_NAME = re.compile(r"(?P<proc>[^ ]+)(?: \((?P<tag>.*)\))?")
_DONE = {"COMPLETED", "CACHED", "FAILED"}


def write_params(
    outdir: Path,
    cfg: GippaConfig,
    inputs: Mapping[str, Path],
    memory: str = DEFAULT_MEMORY,
    db_dir: Path | None = None,
) -> tuple[Path, Path]:
    """Write ``OUT/work/params.json`` and ``OUT/resolved_config.yaml`` deterministically;
    ``params.db_dir`` (absolute) only when ``db_dir`` is given. ``params.run_uniop`` is true
    when ``cfg.annotation.uniop.enabled`` (already resolved by preflight, M11-14d) is ``true``
    and no precomputed operon layer replaces UniOP."""
    outdir = Path(outdir).absolute()
    resolved = outdir / "resolved_config.yaml"
    params = outdir / "work" / "params.json"
    params.parent.mkdir(parents=True, exist_ok=True)
    dumped = cfg.model_dump(mode="json")
    # precomputed layers are read by py steps inside Nextflow task directories, so a relative
    # path would resolve there, not at launch (code review of #146)
    pre = dumped["annotation"]["precomputed"]
    for key, value in pre.items():
        if value is not None:
            pre[key] = str(Path(value).absolute())
    resolved.write_text(yaml.safe_dump(dumped, sort_keys=True), encoding="utf-8")
    load_step_modules()
    names = [spec.name for spec in reg.all_steps()]
    data = {
        "annotation": dumped["annotation"],
        "config": str(resolved),
        "config_hashes": reg.config_hashes(dumped, names),  # task hashes (M7-15 decision 7)
        "inputs": {k: str(Path(v).absolute()) for k, v in inputs.items()},
        "mapping": dumped["mapping"],
        "memory": memory,
        "modules": list(cfg.modules),
        "outdir": str(outdir),
        "run_uniop": cfg.annotation.uniop.enabled is True
        and cfg.annotation.precomputed.operon_pairs is None,
        "threads": cfg.threads,
    }
    if db_dir is not None:
        data["db_dir"] = str(Path(db_dir).absolute())
    params.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return params, resolved


def nextflow_memory(memory: str) -> str:
    """``--memory 64.GB`` (§5.1 form) → ``"64 GB"``; other forms pass through unchanged."""
    return re.sub(r"^\s*(\d+(?:\.\d+)?)\s*\.?\s*([KMGT]B)\s*$", r"\1 \2", memory)


def launch_command(
    workflow: Path,
    params: Path,
    profile: str,
    outdir: Path,
    resume: bool,
    local_tools: Path | None = None,
) -> list[str]:
    """The ``nextflow run`` command line of §6.3 (no plugins, no ``--param`` flags), with
    ``-c <local_tools>`` when given (a local UniOP build; M11-14 item 1b)."""
    logs = outdir / "logs"
    cmd = ["nextflow", "run", str(workflow), "-params-file", str(params), "-profile", profile]
    cmd += ["-c", str(local_tools)] if local_tools is not None else []
    cmd += ["-work-dir", str(outdir / "work"), "-with-trace", str(logs / "trace.tsv")]
    cmd += ["-with-timeline", str(logs / "timeline.html")]
    cmd += ["-with-report", str(logs / "nf_report.html")]
    cmd += ["-with-dag", str(logs / "dag.html")]  # ADR-0069
    return [*cmd, "-resume"] if resume else cmd


# ``GIPPA_STEP`` aliases whose name is not ``*_STEP`` nor the step
# (``REPORT:REPORT_HTML (report)``), and the sharded TerminatorNet processes (M2-15b).
_ALIAS_STEPS = {
    "report_html": "report",
    "terminatornet_plan": "terminatornet",
    "terminatornet_chunk": "terminatornet",
    "terminatornet_merge": "terminatornet",
}


def _step_name(name: str) -> str:
    """Step of a trace task name: the process name, or the tag for ``GIPPA_STEP`` and an alias
    named ``*_STEP`` (``OPERONS:OPERONS_STEP (operons)`` runs ``operons``); ``REPORT_HTML`` runs
    ``report``."""
    m = _NAME.fullmatch(name.strip())
    proc, tag = (m["proc"], m["tag"]) if m else (name, None)
    step = proc.split(":")[-1].lower()
    return tag if step.endswith("_step") and tag else _ALIAS_STEPS.get(step, step)


def _task_dir(row: Mapping[str, str | None], work: Path) -> str:
    """The ``workdir`` column, else the task directory matching the abbreviated ``hash``."""
    workdir = (row.get("workdir") or "").strip()
    if workdir not in _EMPTY:
        return workdir
    hash_ = (row.get("hash") or "").strip()
    if "/" not in hash_:
        return str(work)
    prefix, rest = hash_.split("/", 1)
    hits = sorted((work / prefix).glob(rest + "*")) if (work / prefix).is_dir() else []
    return str(hits[0] if hits else work / hash_)


def failed_step(trace_tsv: Path) -> tuple[str, str] | None:
    """``(step, workdir)`` of the last FAILED task whose name has no later COMPLETED row.

    A retried-then-succeeded attempt never counts; tasks ABORTED because another task failed (killed
    by the terminate strategy) count only when no FAILED row is left. Without a ``workdir`` column
    (Nextflow's default trace fields) the task directory is rebuilt from ``hash`` under
    ``OUT/work`` (the trace lives in ``OUT/logs``).
    """
    trace_tsv = Path(trace_tsv)
    with trace_tsv.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    completed: set[str] = set()
    aborted: Mapping[str, str | None] | None = None
    for row in reversed(rows):
        name = (row.get("name") or "").strip()
        status = (row.get("status") or "").strip()
        if status == "COMPLETED":
            completed.add(name)
        elif status in _FAILED and name not in completed:
            if status == "FAILED":
                return _step_name(name), _task_dir(row, trace_tsv.parent.parent / "work")
            aborted = aborted or row
    if aborted is None:
        return None
    return _step_name(aborted.get("name") or ""), _task_dir(
        aborted, trace_tsv.parent.parent / "work"
    )


def _task_events(row: Mapping[str, str | None], step: str, outdir: Path) -> Path:
    """The task's published events file, else ``logs/events.jsonl`` in its task directory (a
    failed task publishes nothing)."""
    published = outdir / "logs" / "steps" / f"{step}.events.jsonl"
    if published.is_file():
        return published
    return Path(_task_dir(row, outdir / "work")) / "logs" / "events.jsonl"


def append_step_events(outdir: Path) -> int:
    """Append each step task's events to ``OUT/logs/events.jsonl`` in ``trace.tsv`` order (task
    completion); returns the number of lines appended.

    Lines are copied unchanged, so every event keeps its own ``step``; a step's events are taken
    once (its first finished row) and a missing trace or events file appends nothing.
    """
    trace = outdir / "logs" / "trace.tsv"
    if not trace.is_file():
        return 0
    with trace.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    lines: list[str] = []
    seen: set[str] = set()
    for row in rows:
        m = _NAME.fullmatch((row.get("name") or "").strip())
        step = (m["tag"] if m else None) or (row.get("tag") or "").strip()
        if (row.get("status") or "").strip() not in _DONE or step in _EMPTY or step in seen:
            continue
        events = _task_events(row, step, outdir)
        if events.is_file():
            seen.add(step)
            text = events.read_text(encoding="utf-8")
            lines += [ln for ln in text.splitlines() if ln.strip()]
    if lines:
        with (outdir / "logs" / "events.jsonl").open("a", encoding="utf-8") as fh:
            fh.write("".join(f"{ln}\n" for ln in lines))
    return len(lines)


RUNTIME_ENV = "GIPPA_CONTAINER_RUNTIME"


CACHE_ENV = "NXF_APPTAINER_CACHEDIR"


def launch_env(
    environ: Mapping[str, str], runtime: str | None, db_dir: Path | None
) -> dict[str, str]:
    """``environ`` for ``nextflow run``: ``GIPPA_CONTAINER_RUNTIME`` set to ``runtime`` (an
    inherited
    value is dropped when ``runtime`` is ``None``) and ``GIPPA_DB_DIR`` to the absolute ``db_dir``
    when given, so the function step finds the fgsea image where preflight checked it (M11-14
    14c). With ``db_dir`` and no non-empty inherited ``NXF_APPTAINER_CACHEDIR``, that variable is
    set to :func:`gippa.io.image_key.image_cache` of ``environ`` and the absolute ``db_dir``, so
    Nextflow's ``apptainer.cacheDir`` is the cache setup pulled into and preflight checked.
    ``environ`` itself is not modified."""
    from gippa.io.image_key import image_cache

    env = {k: v for k, v in environ.items() if k != RUNTIME_ENV}
    if runtime is not None:
        env[RUNTIME_ENV] = runtime
    if db_dir is not None:
        db = Path(db_dir).absolute()
        env["GIPPA_DB_DIR"] = str(db)
        if not env.get(CACHE_ENV):
            env[CACHE_ENV] = str(image_cache(environ, db))
    return env


def run_nextflow(
    cmd: Sequence[str], outdir: Path, env: Mapping[str, str] | None = None
) -> ExitCode:
    """Run Nextflow in ``outdir`` (with ``env``, else the inherited environment) and append the
    step events (:func:`append_step_events`); 0 on success, else 4 with a message naming step
    and log."""
    proc = subprocess.run(
        list(cmd), cwd=outdir, check=False, env=None if env is None else dict(env)
    )
    append_step_events(outdir)
    if proc.returncode == 0:
        return ExitCode.SUCCESS
    trace = outdir / "logs" / "trace.tsv"
    hit = failed_step(trace) if trace.is_file() else None
    if hit is None:
        msg = f"nextflow exited {proc.returncode}; log: {outdir / '.nextflow.log'}"
    else:
        msg = f"step {hit[0]} failed; log: {hit[1]}/.command.log"
    print(f"gippa run: {msg}", file=sys.stderr)
    return ExitCode.STEP_FAILED
