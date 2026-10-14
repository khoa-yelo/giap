"""Nextflow launcher (GIAP_SPEC.md §6.3 "Launch", "Error strategy"; M2-06).

Python validates everything and writes ``OUT/work/params.json`` plus ``OUT/resolved_config.yaml``;
Nextflow only reads the params file. A failed run maps to exit code 4 naming the failed step and the
``.command.log`` of its task directory (read from ``OUT/logs/trace.tsv``). All runner messages go to
stderr.

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

from giap.config import GiapConfig
from giap.exitcodes import ExitCode

DEFAULT_MEMORY = "64 GB"  # §5.1 `giap run --memory 64.GB`, written in a Nextflow-parseable form
_FAILED = {"FAILED", "ABORTED"}
_EMPTY = {"", "-"}
_NAME = re.compile(r"(?P<proc>[^ ]+)(?: \((?P<tag>.*)\))?")


def write_params(
    outdir: Path, cfg: GiapConfig, inputs: Mapping[str, Path], memory: str = DEFAULT_MEMORY
) -> tuple[Path, Path]:
    """Write ``OUT/work/params.json`` and ``OUT/resolved_config.yaml`` deterministically."""
    outdir = Path(outdir).absolute()
    resolved = outdir / "resolved_config.yaml"
    params = outdir / "work" / "params.json"
    params.parent.mkdir(parents=True, exist_ok=True)
    dumped = cfg.model_dump(mode="json")
    resolved.write_text(yaml.safe_dump(dumped, sort_keys=True), encoding="utf-8")
    data = {
        "annotation": dumped["annotation"],
        "config": str(resolved),
        "inputs": {k: str(Path(v).absolute()) for k, v in inputs.items()},
        "memory": memory,
        "modules": list(cfg.modules),
        "outdir": str(outdir),
        "threads": cfg.threads,
    }
    params.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return params, resolved


def launch_command(
    workflow: Path, params: Path, profile: str, outdir: Path, resume: bool
) -> list[str]:
    """The ``nextflow run`` command line of §6.3 (no plugins, no ``--param`` flags)."""
    logs = outdir / "logs"
    cmd = ["nextflow", "run", str(workflow), "-params-file", str(params), "-profile", profile]
    cmd += ["-work-dir", str(outdir / "work"), "-with-trace", str(logs / "trace.tsv")]
    cmd += ["-with-timeline", str(logs / "timeline.html")]
    cmd += ["-with-report", str(logs / "nf_report.html")]
    return [*cmd, "-resume"] if resume else cmd


def _step_name(name: str) -> str:
    """Step of a trace task name: the process name, or the tag for ``GIAP_STEP``."""
    m = _NAME.fullmatch(name.strip())
    proc, tag = (m["proc"], m["tag"]) if m else (name, None)
    step = proc.split(":")[-1].lower()
    return tag if step == "giap_step" and tag else step


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
    """``(step, workdir)`` of the last FAILED/ABORTED task whose name has no later COMPLETED row.

    A retried-then-succeeded attempt never counts. Without a ``workdir`` column (Nextflow's default
    trace fields) the task directory is rebuilt from ``hash`` under ``OUT/work`` (the trace lives in
    ``OUT/logs``).
    """
    trace_tsv = Path(trace_tsv)
    with trace_tsv.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    completed: set[str] = set()
    for row in reversed(rows):
        name = (row.get("name") or "").strip()
        status = (row.get("status") or "").strip()
        if status == "COMPLETED":
            completed.add(name)
        elif status in _FAILED and name not in completed:
            return _step_name(name), _task_dir(row, trace_tsv.parent.parent / "work")
    return None


def run_nextflow(cmd: Sequence[str], outdir: Path) -> ExitCode:
    """Run Nextflow in ``outdir``; 0 on success, else 4 with a message naming step and log."""
    proc = subprocess.run(list(cmd), cwd=outdir, check=False)
    if proc.returncode == 0:
        return ExitCode.SUCCESS
    trace = outdir / "logs" / "trace.tsv"
    hit = failed_step(trace) if trace.is_file() else None
    if hit is None:
        msg = f"nextflow exited {proc.returncode}; log: {outdir / '.nextflow.log'}"
    else:
        msg = f"step {hit[0]} failed; log: {hit[1]}/.command.log"
    print(f"giap run: {msg}", file=sys.stderr)
    return ExitCode.STEP_FAILED
