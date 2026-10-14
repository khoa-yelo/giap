"""Nextflow launcher (GIAP_SPEC.md §6.3 "Launch", "Error strategy"; M2-06).

Python validates everything and writes ``OUT/work/params.json`` plus ``OUT/resolved_config.yaml``;
Nextflow only reads the params file. A failed run maps to exit code 4 naming the failed step and the
``.command.log`` of its task directory (read from ``OUT/logs/trace.tsv``).
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


def failed_step(trace_tsv: Path) -> tuple[str, str] | None:
    """``(step, workdir)`` of the first FAILED or ABORTED task in a Nextflow trace, else None."""
    with Path(trace_tsv).open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status") not in _FAILED:
                continue
            m = _NAME.fullmatch(row["name"].strip())
            proc, tag = (m["proc"], m["tag"]) if m else (row["name"], None)
            step = proc.split(":")[-1].lower()
            if step == "giap_step" and tag:
                step = tag
            return step, row["workdir"]
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
