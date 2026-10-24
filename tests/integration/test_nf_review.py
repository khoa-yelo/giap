"""Nextflow-side review fixes (M2-06 "Lead decisions after code review", items 1, 3, 4, 5, 7).

- (1) ``trace.fields`` in the config includes ``workdir``; a real failing stub run (exit 1, and
  exit 137 which is retried once and fails again) driven through ``run_nextflow`` exits 4 and names
  the step and its existing ``.command.log``.
- (3) Under ``farmshare``: ``executor.queueSize`` x ``process.resourceLimits.cpus`` <= 128 (R14),
  resourceLimits cpus <= 16 and memory <= 64 GB, and every label's evaluated task has
  ``cpus >= ceil(memory / 4 GB)`` and stays within the limits, also after the memory x 2 retry.
  The closures are evaluated by a throw-away pipeline run with ``-profile farmshare`` and the local
  executor forced by a later ``-c`` file (nothing reaches Slurm).
- (4) ``farmshare_local``: ``executor.cpus`` = ``SLURM_CPUS_PER_TASK`` (default 4) and
  ``executor.memory`` = ``SLURM_MEM_PER_NODE`` (MB) - 2 GB (default 30 GB).
- (5) User ``threads``/``memory`` (params.json) and the profile limits combine as the minimum.
- (7) ``apptainer.cacheDir`` defaults to ``$GIAP_ROOT/containers`` when ``NXF_APPTAINER_CACHEDIR``
  is unset.
"""

from __future__ import annotations

import csv
import importlib
import json
import math
import os
import re
import subprocess
from pathlib import Path

import pytest

from giap.config import resolve_config
from giap.exitcodes import ExitCode

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
GB = 2**30
_UNITS = {"B": 1, "KB": 2**10, "MB": 2**20, "GB": 2**30, "TB": 2**40}
_SLURM_VARS = (
    "SLURM_CPUS_ON_NODE",
    "SLURM_CPUS_PER_TASK",
    "SLURM_MEM_PER_NODE",
    "SLURM_MEM_PER_CPU",
)


def _env(**extra: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if k not in _SLURM_VARS}
    env.pop("NXF_APPTAINER_CACHEDIR", None)
    env["NXF_ANSI_LOG"] = "false"
    env.update(extra)
    return env


def _config(profile: str, **env: str) -> dict[str, str]:
    cmd = ["nextflow", "config", str(WF), "-flat", "-profile", profile]
    res = subprocess.run(cmd, capture_output=True, text=True, env=_env(**env), timeout=300)
    assert res.returncode == 0, f"{' '.join(cmd)} failed:\n{res.stdout}\n{res.stderr}"
    flat: dict[str, str] = {}
    for line in res.stdout.splitlines():
        m = re.match(r"^([\w.$'\":-]+) = (.*)$", line)
        if m:
            flat[m.group(1)] = m.group(2).strip()
    return flat


def _unq(value: str) -> str:
    return value.strip().strip("'\"")


def _mem(value: str) -> int:
    m = re.fullmatch(r"\s*'?\"?([\d.]+)\s*\.?\s*([KMGT]?B)'?\"?\s*", value, flags=re.I)
    assert m, f"not a memory value: {value!r}"
    return int(float(m.group(1)) * _UNITS[m.group(2).upper()])


def _get(flat: dict[str, str], *keys: str) -> str:
    for k in keys:
        if k in flat:
            return flat[k]
    raise AssertionError(f"none of {keys} in config:\n" + "\n".join(sorted(flat)))


# ---------------------------------------------------------------- decision 1


def test_trace_fields_include_workdir() -> None:
    flat = _config("standard")
    fields = [f.strip() for f in _unq(_get(flat, "trace.fields")).split(",")]
    for f in ("task_id", "hash", "name", "tag", "status", "exit", "workdir"):
        assert f in fields, f"trace.fields lacks {f}: {fields}"


_BROKEN = """
process BROKEN {
    tag 'chr1'
    label 'process_single'
    script:
    "exit %(code)d"
    stub:
    "exit %(code)d"
}
workflow {
    BROKEN()
}
"""


@pytest.mark.parametrize("code", [1, 137])
def test_real_failing_stub_run_through_run_nextflow(
    tmp_path: Path, code: int, capfd: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    nf = importlib.import_module("giap.runner.nextflow")
    proj = tmp_path / "proj"
    proj.mkdir()
    (proj / "main.nf").write_text(_BROKEN % {"code": code})
    (proj / "nextflow.config").write_text(f"includeConfig '{WF / 'nextflow.config'}'\n")
    out = tmp_path / "OUT"
    params, _ = nf.write_params(out, resolve_config(), {})
    cmd = nf.launch_command(proj / "main.nf", params, "test,stub", out, resume=False)
    for k in list(os.environ):
        if k in _SLURM_VARS:
            monkeypatch.delenv(k)
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    rc = nf.run_nextflow([*cmd, "-stub-run"], out)
    assert rc == ExitCode.STEP_FAILED
    err = capfd.readouterr().err
    m = re.search(r"step (\S+) failed; log: (\S+)", err)
    assert m, f"no step message on stderr:\n{err}"
    assert m.group(1) == "broken"
    log = Path(m.group(2))
    assert log.name == ".command.log"
    assert log.parent.is_dir() and out / "work" in log.parents, log
    assert log.is_file(), f"{log} does not exist"


# ---------------------------------------------------------------- decision 3


def test_farmshare_r14_concurrent_cpus_at_most_128() -> None:
    flat = _config("farmshare")
    queue = int(_unq(_get(flat, "executor.queueSize", "executor.$slurm.queueSize")))
    cpus = int(_unq(flat["process.resourceLimits.cpus"]))
    assert queue == 8
    assert cpus <= 16
    assert queue * cpus <= 128
    assert _mem(flat["process.resourceLimits.memory"]) <= 64 * GB


_LABELS = {
    "SINGLE": "process_single",
    "LOW": "process_low",
    "MEDIUM": "process_medium",
    "HIGH": "process_high",
}


def _labels_script(fail_first: bool = True) -> str:
    """One process per label; ``fail_first`` makes attempt 1 exit 137 (the memory retry)."""
    body = "if [ ${task.attempt} -eq 1 ]; then exit 137; fi" if fail_first else "true"
    procs = [
        f"process {name} {{\n    label '{label}'\n    script:\n    \"{body}\"\n}}\n"
        for name, label in _LABELS.items()
    ]
    return "".join(procs) + "workflow {\n" + "".join(f"    {n}()\n" for n in _LABELS) + "}\n"


def _override(executor_cpus: int | None) -> str:
    cpus = f"    cpus = {executor_cpus}\n    memory = '512 GB'\n" if executor_cpus else ""
    return (
        "process.executor = 'local'\nprocess.scratch = false\n"
        f"executor {{\n    name = 'local'\n{cpus}}}\n"
        "trace {\n    enabled = true\n    overwrite = true\n    file = 'trace.tsv'\n"
        "    fields = 'name,status,exit,cpus,memory,attempt'\n}\n"
    )


def _profile_run(
    tmp_path: Path,
    profile: str,
    threads: int,
    memory: str,
    executor_cpus: int | None,
    *,
    fail_first: bool = True,
    **env: str,
) -> list[dict[str, str]]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / "main.nf").write_text(_labels_script(fail_first))
    (tmp_path / "override.config").write_text(_override(executor_cpus))
    (tmp_path / "params.json").write_text(json.dumps({"threads": threads, "memory": memory}))
    cmd = ["nextflow", "-q", "run", "main.nf", "-c", str(WF / "nextflow.config")]
    cmd += ["-c", "override.config", "-profile", profile]
    cmd += ["-params-file", "params.json", "-work-dir", "work"]
    res = subprocess.run(
        cmd, cwd=tmp_path, capture_output=True, text=True, env=_env(**env), timeout=600
    )
    assert res.returncode == 0, f"{' '.join(cmd)}\n{res.stdout}\n{res.stderr}"
    with (tmp_path / "trace.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_farmshare_labels_cpus_cover_memory_incl_retry(tmp_path: Path) -> None:
    rows = _profile_run(tmp_path, "farmshare", 64, "256 GB", executor_cpus=64)
    seen = {(r["name"], r["attempt"]) for r in rows}
    for name in _LABELS:
        assert {(name, "1"), (name, "2")} <= seen, f"{name}: no first attempt and retry"
    for r in rows:
        cpus, mem = int(r["cpus"]), _mem(r["memory"])
        where = f"{r['name']} attempt {r['attempt']}: cpus={cpus} memory={r['memory']}"
        assert cpus <= 16, where
        assert mem <= 64 * GB, where
        assert cpus >= math.ceil(mem / (4 * GB)), where
    by = {(r["name"], r["attempt"]): r for r in rows}
    # the retry still doubles memory up to the 64 GB cap
    for name in _LABELS:
        m1, m2 = _mem(by[(name, "1")]["memory"]), _mem(by[(name, "2")]["memory"])
        assert m2 == min(2 * m1, 64 * GB), name


# ---------------------------------------------------------------- decision 4


def test_farmshare_local_executor_follows_slurm_allocation() -> None:
    flat = _config("farmshare_local", SLURM_CPUS_PER_TASK="8", SLURM_MEM_PER_NODE="65536")
    assert _unq(flat.get("process.executor", "'local'")) == "local"
    assert int(_unq(_get(flat, "executor.cpus", "executor.$local.cpus"))) == 8
    assert _mem(_get(flat, "executor.memory", "executor.$local.memory")) == 62 * GB


def test_farmshare_local_executor_defaults_4_cpus_30_gb() -> None:
    flat = _config("farmshare_local")
    assert int(_unq(_get(flat, "executor.cpus", "executor.$local.cpus"))) == 4
    assert _mem(_get(flat, "executor.memory", "executor.$local.memory")) == 30 * GB


# ---------------------------------------------------------------- decision 5


# The limits tests run a script that succeeds: under the caps here the first request is already at
# the cap, so a 137 terminates without a retry (M7-15 decision 11 B1; test_m7_15_retry_growth.py).


def test_farmshare_user_limits_below_profile_win(tmp_path: Path) -> None:
    rows = _profile_run(tmp_path, "farmshare", 4, "8 GB", executor_cpus=64, fail_first=False)
    for r in rows:
        assert int(r["cpus"]) <= 4, r
        assert _mem(r["memory"]) <= 8 * GB, r
    high = next(r for r in rows if r["name"] == "HIGH" and r["attempt"] == "1")
    assert int(high["cpus"]) == 4
    assert _mem(high["memory"]) == 8 * GB


def test_farmshare_local_user_limits_vs_allocation(tmp_path: Path) -> None:
    env = {"SLURM_CPUS_PER_TASK": "8", "SLURM_MEM_PER_NODE": "34816"}  # 34 GB -> 32 GB
    big = _profile_run(
        tmp_path / "big", "farmshare_local", 64, "256 GB", None, fail_first=False, **env
    )
    high = [r for r in big if r["name"] == "HIGH"]
    assert high
    for r in high:
        assert int(r["cpus"]) == 8, r
        assert _mem(r["memory"]) <= 32 * GB, r
    small = _profile_run(
        tmp_path / "small", "farmshare_local", 2, "4 GB", None, fail_first=False, **env
    )
    for r in small:
        assert int(r["cpus"]) <= 2, r
        assert _mem(r["memory"]) <= 4 * GB, r


# ---------------------------------------------------------------- decision 7


def test_apptainer_cache_dir_defaults_to_giap_root() -> None:
    flat = _config("apptainer", GIAP_ROOT="/scratch/giap_root")
    assert _unq(flat["apptainer.cacheDir"]) == "/scratch/giap_root/containers"


def test_apptainer_cache_dir_env_wins() -> None:
    flat = _config("apptainer", GIAP_ROOT="/scratch/giap_root", NXF_APPTAINER_CACHEDIR="/c/x")
    assert _unq(flat["apptainer.cacheDir"]) == "/c/x"
