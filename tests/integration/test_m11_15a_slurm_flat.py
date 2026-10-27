"""M11-15a: the generic ``slurm`` profile merged by Nextflow (spec M11-15 item 3, example 3.1).

``nextflow config workflow -flat -profile slurm,apptainer`` with no ``slurm_*`` params set:
``process.executor = 'slurm'``, no ``process.queue``, an empty ``clusterOptions`` and
``executor.queueSize = 20``. The flag-driven examples (3.2-3.9) belong to 15e. Nothing is
submitted.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.needs_nextflow

WF = Path(__file__).resolve().parents[2] / "workflow"


def _flat(profile: str) -> dict[str, str]:
    drop = {"SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU"}
    env = {k: v for k, v in os.environ.items() if k not in drop} | {"NXF_ANSI_LOG": "false"}
    cmd = ["nextflow", "config", str(WF), "-flat", "-profile", profile]
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300, check=False)
    assert res.returncode == 0, f"{' '.join(cmd)} failed:\n{res.stdout}\n{res.stderr}"
    out: dict[str, str] = {}
    for line in res.stdout.splitlines():
        m = re.match(r"^([\w.$'\":-]+) = (.*)$", line)
        if m:
            out[m.group(1)] = m.group(2).strip()
    return out


def test_3_1_slurm_apptainer_defaults() -> None:
    flat = _flat("slurm,apptainer")
    assert flat.get("process.executor", "").strip("'\"") == "slurm"
    assert flat.get("process.queue", "null").strip("'\"") in ("", "null"), "no -p by default"
    assert flat.get("process.clusterOptions", "''").strip("'\"") == ""
    assert flat.get("executor.queueSize", "").strip("'\"") == "20"
    assert "executor.submitRateLimit" in flat
    assert "process.scratch" not in flat
