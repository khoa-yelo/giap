"""M11-19: Nextflow's docker tasks run with ``--pull=never`` (examples 18-21, AC 7 and 8).

docs/specs/M11-19-nextflow-docker-pull-never.md item 1. Two layers, no Docker daemon, no
network (helpers of test_m11_14_docker_profile_stub.py):

1. ``nextflow config <wf> -flat -profile <p>``: ``docker.runOptions`` is
   ``--pull=never -u $(id -u):$(id -g)`` under ``docker`` and ``docker,conda``, and
   ``--pull=never -v /db:/db:ro -u $(id -u):$(id -g)`` when ``params.db_dir`` is ``/db``
   (``nextflow config`` takes no params, so a copy of ``nextflow.config`` with that default is
   read, every other entry of ``workflow/`` symlinked); no ``--pull`` under ``apptainer``, and
   no ``docker.runOptions`` under ``apptainer``, ``conda`` and ``test,stub``.
2. ``-stub-run`` of ``main.nf`` under ``-profile test,docker,stub`` with the logging fake
   ``docker``: every recorded ``docker run`` argv carries ``--pull=never`` exactly once, before
   the image ref.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path

import pytest
from test_m11_14_docker_profile_stub import (
    DATA,
    FAKE_DOCKER,
    MAIN,
    WF,
    _env,
    _pairs,
    _unq,
)
from test_m11_14_docker_profile_stub import (
    _config as _flat,
)

from gippa.runner.nextflow import launch_command

pytestmark = pytest.mark.needs_nextflow

USER = "-u $(id -u):$(id -g)"


def _config_at(wf: Path, profile: str) -> dict[str, str]:
    cmd = ["nextflow", "config", str(wf), "-flat", "-profile", profile]
    res = subprocess.run(cmd, capture_output=True, text=True, env=_env(), timeout=300)
    assert res.returncode == 0, f"{' '.join(cmd)} failed:\n{res.stdout}\n{res.stderr}"
    flat: dict[str, str] = {}
    for line in res.stdout.splitlines():
        m = re.match(r"^(\S+) = (.*)$", line)
        if m:
            flat[m.group(1)] = m.group(2).strip()
    return flat


# --- 1. resolved configuration --------------------------------------------------------------------


@pytest.mark.parametrize("profile", ["docker", "docker,conda"])
def test_examples_18_20_run_options_start_with_pull_never(profile: str) -> None:
    opts = _unq(_flat(profile).get("docker.runOptions", ""))
    assert opts == f"--pull=never {USER}", opts


def test_example_19_with_db_dir(tmp_path: Path) -> None:
    wf = tmp_path / "wf"
    wf.mkdir()
    for entry in WF.iterdir():
        if entry.name != "nextflow.config":
            (wf / entry.name).symlink_to(entry)
    text = (WF / "nextflow.config").read_text(encoding="utf-8")
    patched, n = re.subn(r"(?m)^(\s*db_dir\s*=\s*)null\b", r"\1'/db'", text)
    assert n == 1, "params.db_dir default not found in nextflow.config"
    (wf / "nextflow.config").write_text(patched, encoding="utf-8")
    opts = _unq(_config_at(wf, "docker").get("docker.runOptions", ""))
    assert opts == f"--pull=never -v /db:/db:ro {USER}", opts


def test_example_20_apptainer_has_no_pull_and_no_docker_options() -> None:
    flat = _flat("apptainer")
    assert "--pull" not in flat.get("apptainer.runOptions", "")
    assert "docker.runOptions" not in flat, flat.get("docker.runOptions")


@pytest.mark.parametrize("profile", ["test,stub", "conda"])
def test_example_20_no_docker_options_without_docker(profile: str) -> None:
    assert "docker.runOptions" not in _flat(profile)


# --- 2. stub run with a logging fake docker -------------------------------------------------------


def _runs(log: Path) -> list[list[str]]:
    out = []
    for record in log.read_bytes().split(b"\0\n"):
        argv = [a.decode() for a in record.split(b"\0") if a]
        if argv:
            out.append(argv)
    return out


def _image_index(argv: list[str]) -> int:
    """The image ref precedes the task's shell (the fake runs from ``/bin/bash``/``/bin/sh``)."""
    for i, a in enumerate(argv):
        if a in ("/bin/bash", "/bin/sh"):
            return i - 1
    raise AssertionError(f"no task shell in {argv}")


@pytest.fixture(scope="module")
def stub_run(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    tmp = tmp_path_factory.mktemp("docker_pull_never_stub")
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
    params.write_text(json.dumps({
        "outdir": str(out),
        "modules": ["core", "bgc"],
        "db_dir": str(db),
        "inputs": {"genome": str(DATA / "genome.fa"), "replicons": str(DATA / "replicons.tsv")},
        "annotation": json.loads((DATA / "annotation.json").read_text()),
    }))  # fmt: skip
    log = tmp / "docker.log"
    log.touch()
    cmd = [*launch_command(MAIN, params, "test,docker,stub", out, resume=False), "-stub-run"]
    env = _env(PATH=f"{fake}:{os.environ['PATH']}", FAKE_DOCKER_LOG=str(log))
    res = subprocess.run(cmd, cwd=out, capture_output=True, text=True, env=env, timeout=1200)
    assert res.returncode == 0, res.stdout + res.stderr
    return {"db": db, "runs": _runs(log)}


def test_example_21_every_docker_run_never_pulls(stub_run: dict[str, object]) -> None:
    runs = stub_run["runs"]
    assert isinstance(runs, list) and runs, "no docker run was logged"
    for argv in runs:
        assert argv.count("--pull=never") == 1, argv
        assert not any(a == "--pull" or (a.startswith("--pull=") and a != "--pull=never")
                       for a in argv), argv  # fmt: skip
        assert argv.index("--pull=never") < _image_index(argv), argv


@pytest.mark.parametrize("image", ["oschwengers/bakta@", "antismash/standalone@"])
def test_example_21_pinned_processes_are_covered(stub_run: dict[str, object], image: str) -> None:
    runs = stub_run["runs"]
    assert isinstance(runs, list)
    hits = [a for a in runs if a[_image_index(a)].startswith(image)]
    assert hits, f"no docker run of {image}"
    assert all("--pull=never" in a for a in hits)


def test_example_19_pull_never_heads_run_options_in_the_task(stub_run: dict[str, object]) -> None:
    """``--pull=never`` is the first token of runOptions: the DB bind and ``-u`` follow it."""
    db = stub_run["db"]
    runs = stub_run["runs"]
    assert isinstance(runs, list)
    bakta = [a for a in runs if a[_image_index(a)].startswith("oschwengers/bakta@")]
    assert bakta
    argv = bakta[0]
    i = argv.index("--pull=never")
    assert argv[i + 1 : i + 3] == ["-v", f"{db}:{db}:ro"], argv
    assert argv[i + 3] == "-u" or argv[i + 3].startswith("-u"), argv
    assert f"{db}:{db}:ro" in _pairs(argv, "-v")
