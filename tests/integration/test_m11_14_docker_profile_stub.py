"""The docker profile as Nextflow resolves and launches it (M11-14 item 3b, acceptance 11; 14p).

Two layers, no Docker daemon and no network:

1. ``nextflow config workflow -flat -profile <p>``: the resolved docker settings match the
   apptainer ones (TMPDIR ``beforeScript``, ps shim on ``BAKTA|ANTISMASH``, the same conda-off
   selector, also under ``docker,conda``), ``docker.runOptions`` always carries
   ``-u $(id -u):$(id -g)`` and no DB bind (and no ``null``) when ``params.db_dir`` is unset.
2. ``-stub-run`` of ``main.nf`` (modules core + bgc) under ``-profile test,docker,stub`` with a
   logging fake ``docker`` on ``PATH``: it records each ``docker run`` argv and then runs the task
   command on the host. The BAKTA call carries ``-v <db_dir>:<db_dir>:ro``, ``-u <uid>:<gid>``
   and the ps shim; ANTISMASH the ps shim; the task's ``TMPDIR`` is its own ``tmp/`` and is passed
   into the container (a ``beforeScript`` runs on the host, so Docker needs ``-e TMPDIR``;
   Apptainer passes the environment through by default).
"""

from __future__ import annotations

import csv
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

from gippa.runner.nextflow import launch_command

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
MAIN = WF / "main.nf"
DATA = WF / "tests" / "data" / "annotate"
TMPDIR_SCRIPT = "mkdir -p tmp; export TMPDIR=$PWD/tmp"
PS_TARGET = "/assets/ps/ps:/usr/bin/ps:ro"
CONDA_OFF_KEY = re.compile(r"^process\.'withName:([A-Za-z0-9_|]+)'\.conda$")

FAKE_DOCKER = """#!/usr/bin/env bash
# logging fake docker (M11-14 14p tests): record each `docker run` argv, run the task on the host
if [ "$1" = run ]; then
  { printf '%s\\0' "$@"; printf '\\n'; } >> "$FAKE_DOCKER_LOG"
  shift
  while [ $# -gt 0 ] && [ "$1" != /bin/bash ] && [ "$1" != /bin/sh ]; do shift; done
  exec "$@"
fi
exit 0
"""


def _env(**extra: str) -> dict[str, str]:
    drop = {"SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU"}
    env = {k: v for k, v in os.environ.items() if k not in drop}
    env.pop("NXF_APPTAINER_CACHEDIR", None)
    env["NXF_ANSI_LOG"] = "false"
    env.update(extra)
    return env


def _unq(value: str) -> str:
    return value.strip().strip("'\"")


def _config(profile: str) -> dict[str, str]:
    cmd = ["nextflow", "config", str(WF), "-flat", "-profile", profile]
    res = subprocess.run(cmd, capture_output=True, text=True, env=_env(), timeout=300)
    assert res.returncode == 0, f"{' '.join(cmd)} failed:\n{res.stdout}\n{res.stderr}"
    flat: dict[str, str] = {}
    for line in res.stdout.splitlines():
        m = re.match(r"^(\S+) = (.*)$", line)
        if m:
            flat[m.group(1)] = m.group(2).strip()
    return flat


def _conda_off(flat: dict[str, str]) -> set[str]:
    names: set[str] = set()
    for key, value in flat.items():
        m = CONDA_OFF_KEY.match(key)
        if m and _unq(value) == "":
            names |= set(m.group(1).split("|"))
    return names


def _shim(flat: dict[str, str]) -> str:
    for key, value in flat.items():
        m = re.match(r"^process\.'withName:([A-Za-z0-9_|]+)'\.containerOptions$", key)
        if m and set(m.group(1).split("|")) == {"BAKTA", "ANTISMASH"}:
            return _unq(value)
    raise AssertionError("no containerOptions for BAKTA|ANTISMASH:\n" + "\n".join(sorted(flat)))


# --- 1. resolved configuration --------------------------------------------------------------------


@pytest.fixture(scope="module")
def docker_cfg() -> dict[str, str]:
    return _config("docker")


@pytest.fixture(scope="module")
def apptainer_cfg() -> dict[str, str]:
    return _config("apptainer")


def test_docker_before_script_sets_task_tmpdir(docker_cfg: dict[str, str]) -> None:
    assert _unq(docker_cfg.get("process.beforeScript", "")) == TMPDIR_SCRIPT


def test_docker_ps_shim_matches_apptainer(
    docker_cfg: dict[str, str], apptainer_cfg: dict[str, str]
) -> None:
    assert _shim(docker_cfg) == f"-v {WF}{PS_TARGET}"
    assert _shim(apptainer_cfg) == f"-B {WF}{PS_TARGET}"


def test_docker_runs_as_user_without_db_bind_when_unset(docker_cfg: dict[str, str]) -> None:
    opts = _unq(docker_cfg.get("docker.runOptions", ""))
    assert "-u $(id -u):$(id -g)" in opts, opts
    assert "-v" not in opts.split() and "null" not in opts, opts


def test_docker_passes_tmpdir_into_the_container(docker_cfg: dict[str, str]) -> None:
    allowed = _unq(docker_cfg.get("docker.envWhitelist", ""))
    opts = _unq(docker_cfg.get("docker.runOptions", ""))
    assert "TMPDIR" in re.split(r"[\s,\[\]'\"]+", allowed) or re.search(
        r"(-e|--env)\s+TMPDIR\b", opts
    ), "a host beforeScript's TMPDIR does not reach a docker container without -e TMPDIR"


@pytest.mark.parametrize("profile", ["docker", "docker,conda"])
def test_same_conda_off_selector_as_apptainer(profile: str, apptainer_cfg: dict[str, str]) -> None:
    flat = _config(profile)
    expected = _conda_off(apptainer_cfg)
    assert expected, "apptainer profile lost its conda-off selector"
    assert _conda_off(flat) == expected
    if "conda" in profile:
        assert _unq(flat.get("conda.enabled", "")) == "true"


# --- 2. stub run with a logging fake docker -------------------------------------------------------


def _argv_by_workdir(log: Path) -> dict[str, list[str]]:
    calls: dict[str, list[str]] = {}
    for record in log.read_bytes().split(b"\0\n"):
        argv = [a.decode() for a in record.split(b"\0") if a]
        if argv and "-w" in argv:
            calls[argv[argv.index("-w") + 1]] = argv
    return calls


def _pairs(argv: list[str], flag: str) -> list[str]:
    return [argv[i + 1] for i, a in enumerate(argv[:-1]) if a == flag]


@pytest.fixture(scope="module")
def stub_run(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    tmp = tmp_path_factory.mktemp("docker_stub")
    fake = tmp / "bin"
    fake.mkdir()
    (fake / "docker").write_text(FAKE_DOCKER)
    (fake / "docker").chmod(0o755)
    db = tmp / "db"
    db.mkdir()
    (db / "version.json").write_text('{"type": "light"}\n')
    out = tmp / "OUT"
    out.mkdir()
    params = tmp / "params.json"
    params.write_text(
        json.dumps(
            {
                "outdir": str(out),
                "modules": ["core", "bgc"],
                "db_dir": str(db),
                "inputs": {
                    "genome": str(DATA / "genome.fa"),
                    "replicons": str(DATA / "replicons.tsv"),
                },
                "annotation": json.loads((DATA / "annotation.json").read_text()),
            }
        )
    )
    log = tmp / "docker.log"
    log.touch()
    cmd = [*launch_command(MAIN, params, "test,docker,stub", out, resume=False), "-stub-run"]
    env = _env(PATH=f"{fake}:{os.environ['PATH']}", FAKE_DOCKER_LOG=str(log))
    res = subprocess.run(cmd, cwd=out, capture_output=True, text=True, env=env, timeout=1200)
    assert res.returncode == 0, res.stdout + res.stderr
    with (out / "logs" / "trace.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    calls = _argv_by_workdir(log)
    by_proc = {}
    for r in rows:
        proc = r["name"].split(" (")[0].split(":")[-1]
        if proc in {"BAKTA", "ANTISMASH"}:
            by_proc[proc] = (Path(r["workdir"]), calls.get(r["workdir"]))
    return {"db": db, "calls": calls, "procs": by_proc}


def _call(stub_run: dict[str, object], proc: str) -> tuple[Path, list[str]]:
    procs = stub_run["procs"]
    assert isinstance(procs, dict) and proc in procs, f"{proc} did not run"
    workdir, argv = procs[proc]
    assert argv, f"{proc} was not launched through docker run"
    return workdir, argv


def test_bakta_sees_the_db_read_only(stub_run: dict[str, object]) -> None:
    _, argv = _call(stub_run, "BAKTA")
    db = stub_run["db"]
    assert f"{db}:{db}:ro" in _pairs(argv, "-v"), argv


def test_every_docker_call_runs_as_the_user(stub_run: dict[str, object]) -> None:
    calls = stub_run["calls"]
    assert isinstance(calls, dict) and calls, "no docker run was logged"
    uid_gid = f"{os.getuid()}:{os.getgid()}"
    for workdir, argv in calls.items():
        assert uid_gid in _pairs(argv, "-u") + _pairs(argv, "--user"), (workdir, argv)


@pytest.mark.parametrize("proc", ["BAKTA", "ANTISMASH"])
def test_ps_shim_bound(stub_run: dict[str, object], proc: str) -> None:
    _, argv = _call(stub_run, proc)
    assert any(v.endswith(PS_TARGET) for v in _pairs(argv, "-v")), argv


def test_bakta_task_tmpdir_is_its_own_and_passed_in(stub_run: dict[str, object]) -> None:
    workdir, argv = _call(stub_run, "BAKTA")
    assert (workdir / "tmp").is_dir(), "beforeScript did not create the per-task tmp/"
    assert TMPDIR_SCRIPT in (workdir / ".command.run").read_text()
    env_args = _pairs(argv, "-e") + _pairs(argv, "--env")
    assert any(a == "TMPDIR" or a.startswith("TMPDIR=") for a in env_args), argv
