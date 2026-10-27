"""Nextflow configuration and profiles (M2-06 examples 5-7, edge cases; GIAP_SPEC.md §6.3, §3.8).

Examples 5 and 6 parse ``nextflow config workflow -flat``. Example 7 needs evaluated closures
(memory x attempt, error strategy), so it runs a throw-away pipeline in ``tmp_path`` whose
processes carry the four labels, with the repository's ``workflow/nextflow.config`` as config and
a small override (local executor with 64 slots, a trace with cpus/memory/attempt). Only ``true``
and ``exit N`` run; nothing is submitted to Slurm.
"""

from __future__ import annotations

import csv
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

from giap.config import load_defaults

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"

_UNITS = {"B": 1, "KB": 2**10, "MB": 2**20, "GB": 2**30, "TB": 2**40}


def _env(**extra: str) -> dict[str, str]:
    drop = {"SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU"}
    env = {k: v for k, v in os.environ.items() if k not in drop}
    env["NXF_ANSI_LOG"] = "false"
    env.update(extra)
    return env


def _config(profile: str | None, **env: str) -> tuple[dict[str, str], str]:
    # without a nextflow.config Nextflow silently ignores unknown profiles
    assert (WF / "nextflow.config").is_file(), "workflow/nextflow.config missing"
    cmd = ["nextflow", "config", str(WF), "-flat"]
    if profile:
        cmd += ["-profile", profile]
    res = subprocess.run(cmd, capture_output=True, text=True, env=_env(**env), timeout=300)
    assert res.returncode == 0, f"{' '.join(cmd)} failed:\n{res.stdout}\n{res.stderr}"
    flat: dict[str, str] = {}
    for line in res.stdout.splitlines():
        m = re.match(r"^([\w.$'\":-]+) = (.*)$", line)
        if m:
            flat[m.group(1)] = m.group(2).strip()
    return flat, res.stdout + res.stderr


def _unq(value: str) -> str:
    return value.strip().strip("'\"")


def _mem(value: str) -> int:
    m = re.fullmatch(r"\s*'?\"?([\d.]+)\s*\.?\s*([KMGT]?B)'?\"?\s*", value, flags=re.I)
    assert m, f"not a memory value: {value!r}"
    return int(float(m.group(1)) * _UNITS[m.group(2).upper()])


def _hours(value: str) -> float:
    v = _unq(value)
    m = re.fullmatch(r"([\d.]+)\s*\.?\s*(d|h|day|days|hour|hours)", v)
    assert m, f"not a duration: {value!r}"
    n = float(m.group(1))
    return n * 24 if m.group(2).startswith("d") else n


def _get(flat: dict[str, str], *keys: str) -> str:
    for k in keys:
        if k in flat:
            return flat[k]
    raise AssertionError(f"none of {keys} in config:\n" + "\n".join(sorted(flat)))


def _executor(flat: dict[str, str]) -> str:
    return _unq(flat.get("process.executor", flat.get("executor.name", "'local'")))


# ---------------------------------------------------------------- example 5: farmshare


def test_farmshare_profile_example_5() -> None:
    flat, _ = _config("farmshare")
    assert _executor(flat) == "slurm"
    assert _unq(flat["process.queue"]) == "normal"
    # slow-node stopgap (#416 D3): Python/JVM start-up stalls on these nodes
    slow = "--exclude=barley-01,barley-03,rye-02"
    assert _unq(flat["process.clusterOptions"]) == f"--qos=normal {slow}"
    # review decision 3 (R14): queueSize 8 x resourceLimits.cpus 16 <= 128 (was 20 / 64 / 240 GB)
    assert _get(flat, "executor.queueSize", "executor.$slurm.queueSize") == "8"
    rate = _get(flat, "executor.submitRateLimit", "executor.$slurm.submitRateLimit")
    assert _unq(rate) == "10/1min"
    assert _unq(flat["process.scratch"]) == "/tmp"
    assert flat["process.resourceLimits.cpus"] == "16"
    assert _mem(flat["process.resourceLimits.memory"]) == 64 * 2**30
    assert _hours(flat["process.resourceLimits.time"]) == 48


def test_farmshare_job_names_start_with_giap() -> None:
    """R14: every GIAP job name starts with ``giap-`` (the closure is printed as source)."""
    flat, _ = _config("farmshare")
    job = _get(flat, "executor.jobName", "executor.$slurm.jobName")
    assert "giap-" in job
    assert "toLowerCase" in job and "task.process" in job


def test_farmshare_with_apptainer_resolves() -> None:
    flat, _ = _config("farmshare,apptainer")
    assert _executor(flat) == "slurm"
    assert flat["apptainer.enabled"] == "true"


# ---------------------------------------------------------------- example 6: farmshare_local


# Review decision 4: farmshare_local follows SLURM_CPUS_PER_TASK / SLURM_MEM_PER_NODE - 2 GB via
# executor.cpus/memory (was SLURM_CPUS_ON_NODE via resourceLimits); the limits never exceed them.


def test_farmshare_local_profile_example_6() -> None:
    flat, _ = _config("farmshare_local", SLURM_CPUS_PER_TASK="8", SLURM_MEM_PER_NODE="32768")
    assert _executor(flat) == "local"
    assert _get(flat, "executor.cpus", "executor.$local.cpus") == "8"
    assert _mem(_get(flat, "executor.memory", "executor.$local.memory")) == 30 * 2**30
    assert int(flat["process.resourceLimits.cpus"]) <= 8
    assert _mem(flat["process.resourceLimits.memory"]) <= 30 * 2**30


def test_farmshare_local_without_slurm_vars_falls_back_to_4_cpus_30_gb() -> None:
    flat, _ = _config("farmshare_local")
    assert _executor(flat) == "local"
    assert _get(flat, "executor.cpus", "executor.$local.cpus") == "4"
    assert _mem(_get(flat, "executor.memory", "executor.$local.memory")) == 30 * 2**30
    assert int(flat["process.resourceLimits.cpus"]) <= 4
    assert _mem(flat["process.resourceLimits.memory"]) <= 30 * 2**30


# ---------------------------------------------------------------- other profiles


def test_apptainer_profile() -> None:
    flat, _ = _config("apptainer")
    assert flat["apptainer.enabled"] == "true"
    assert flat["apptainer.autoMounts"] == "true"


def test_conda_profile_uses_micromamba() -> None:
    flat, _ = _config("conda")
    assert flat["conda.useMicromamba"] == "true"


def test_docker_profile_parses() -> None:
    flat, _ = _config("docker")
    assert flat["docker.enabled"] == "true"


def test_standard_profile_is_local() -> None:
    flat, _ = _config("standard")
    assert _executor(flat) == "local"


@pytest.mark.parametrize("profile", ["test", "test,stub"])
def test_test_profile_resolves(profile: str) -> None:
    """``ci.yml`` runs ``nf-test test workflow/tests --profile test,stub``."""
    flat, _ = _config(profile)
    assert _executor(flat) == "local"


def test_resource_limits_follow_params() -> None:
    """§6.3: process.resourceLimits = [cpus: params.threads, memory: params.memory]."""
    flat, _ = _config(None)
    assert flat["process.resourceLimits.cpus"] == flat["params.threads"]
    assert _mem(flat["process.resourceLimits.memory"]) == _mem(flat["params.memory"])


def test_params_defaults_match_defaults_yaml() -> None:
    """Edge case: no params.* default differs from ``defaults.yaml``."""
    flat, _ = _config(None)
    d = load_defaults()
    assert int(flat["params.threads"]) == d["threads"]
    assert re.findall(r"[a-z_]+", flat["params.modules"]) == d["modules"]
    for key in ("params.memory", "params.outdir"):
        assert key in flat, f"{key} has no default"
    assert any(k == "params.annotation" or k.startswith("params.annotation.") for k in flat)


# ---------------------------------------------------------------- example 7: labels and retry

_MINI = """
process SINGLE {
    label 'process_single'
    script:
    "true"
}
process LOW {
    label 'process_low'
    script:
    "true"
}
process MEDIUM {
    label 'process_medium'
    script:
    "true"
}
process HIGH {
    label 'process_high'
    script:
    "true"
}
process OOM137 {
    label 'process_medium'
    script:
    "if [ ${task.attempt} -eq 1 ]; then exit 137; fi"
}
process TERM143 {
    label 'process_low'
    script:
    "if [ ${task.attempt} -eq 1 ]; then exit 143; fi"
}
workflow {
    SINGLE()
    LOW()
    MEDIUM()
    HIGH()
    OOM137()
    TERM143()
}
"""

_FAIL = """
process FAIL1 {
    label 'process_medium'
    script:
    "exit 1"
}
workflow {
    FAIL1()
}
"""

_OVERRIDE = """
executor {
    name = 'local'
    cpus = 64
    memory = '512 GB'
}
trace {
    enabled = true
    overwrite = true
    file = 'trace.tsv'
    fields = 'name,status,exit,cpus,memory,attempt'
}
"""


def _mini_run(tmp_path: Path, script: str) -> tuple[int, list[dict[str, str]]]:
    (tmp_path / "main.nf").write_text(script)
    (tmp_path / "override.config").write_text(_OVERRIDE)
    (tmp_path / "params.json").write_text(json.dumps({"threads": 64, "memory": "256 GB"}))
    cmd = ["nextflow", "-q", "run", "main.nf", "-c", str(WF / "nextflow.config")]
    cmd += ["-c", "override.config", "-params-file", "params.json", "-work-dir", "work"]
    res = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True, env=_env(), timeout=600)
    with (tmp_path / "trace.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return res.returncode, rows


def _gb(value: str) -> float:
    return _mem(value) / 2**30


def test_labels_and_retry_on_137_143_example_7(tmp_path: Path) -> None:
    code, rows = _mini_run(tmp_path, _MINI)
    assert code == 0, rows
    first = {r["name"]: r for r in rows if r["attempt"] == "1" and r["status"] == "COMPLETED"}
    expect = {"SINGLE": (1, 2), "LOW": (2, 6), "MEDIUM": (8, 16), "HIGH": (16, 48)}
    for name, (cpus, gb) in expect.items():
        assert int(first[name]["cpus"]) == cpus, name
        assert _gb(first[name]["memory"]) == gb, name
    for name, code_, gb in (("OOM137", "137", 16), ("TERM143", "143", 6)):
        mine = sorted((r for r in rows if r["name"] == name), key=lambda r: int(r["attempt"]))
        assert [r["attempt"] for r in mine] == ["1", "2"], name
        assert (mine[0]["status"], mine[0]["exit"]) == ("FAILED", code_)
        assert mine[1]["status"] == "COMPLETED"
        assert _gb(mine[0]["memory"]) == gb
        assert _gb(mine[1]["memory"]) == 2 * gb, f"{name}: memory not doubled on retry"


def test_exit_1_terminates_without_retry_example_7(tmp_path: Path) -> None:
    code, rows = _mini_run(tmp_path, _FAIL)
    assert code != 0
    mine = [r for r in rows if r["name"] == "FAIL1"]
    assert len(mine) == 1
    assert (mine[0]["status"], mine[0]["exit"], mine[0]["attempt"]) == ("FAILED", "1", "1")


def test_retry_twice_on_137(tmp_path: Path) -> None:
    """maxRetries = 2 (M7-15, Khoa): exit 137 is retried twice, then the run fails."""
    script = """
process ALWAYS137 {
    label 'process_single'
    script:
    "exit 137"
}
workflow {
    ALWAYS137()
}
"""
    code, rows = _mini_run(tmp_path, script)
    assert code != 0
    mine = [r for r in rows if r["name"] == "ALWAYS137"]
    assert sorted(r["attempt"] for r in mine) == ["1", "2", "3"]
