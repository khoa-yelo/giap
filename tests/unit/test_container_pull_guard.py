"""The suite-wide guard against real container pulls (tests/conftest.py ``_no_container_pulls``).

Unit tests of ``gippa setup`` once ran the host's ``apptainer pull`` (about 6 GB per image) and
filled ``/tmp`` (#437, #441). The guard puts fake ``apptainer``/``singularity``/``docker`` first
on ``PATH`` that refuse network pulls and builds and hand everything else to the real binary.

Nothing here can reach a registry: the "real" runtimes are recorder scripts that log their argv
and exit 0. Two layers are tested:

- the guard scripts themselves (``tests/helpers/container_guard.py``), per refused and allowed
  command line;
- the autouse fixture, in a child pytest whose ``PATH`` starts with the recorders: the real
  ``gippa.setupdb.images`` pull and ``gippa.setupdb.uniop_build`` build code get the refusal and
  the recorders never run; a test's own fake first on ``PATH`` wins; ``real_container_pull``
  opts out.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest
from helpers.container_guard import (
    GUARD_ENV,
    NAMES,
    REFUSED,
    REFUSED_EXIT,
    real_binaries,
    write_guard,
)

REPO = Path(__file__).resolve().parents[2]
REF = "quay.io/biocontainers/toy@sha256:" + "0" * 64

RECORDER = """#!/bin/sh
echo "$(basename "$0") $*" >> "{log}"
exit 0
"""


def _recorders(bindir: Path, log: Path) -> dict[str, str]:
    """Recorder ``apptainer``/``singularity``/``docker`` in ``bindir``: log argv, exit 0."""
    bindir.mkdir(parents=True, exist_ok=True)
    out = {}
    for name in NAMES:
        exe = bindir / name
        exe.write_text(RECORDER.format(log=log), encoding="utf-8")
        exe.chmod(0o755)
        out[name] = str(exe)
    return out


@pytest.fixture
def guarded(tmp_path: Path) -> tuple[Path, Path]:
    log = tmp_path / "real.log"
    reals = _recorders(tmp_path / "real", log)
    return write_guard(tmp_path / "guard", reals), log


def _run(guard: Path, argv: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run([str(guard / argv[0]), *argv[1:]], cwd=cwd, capture_output=True,
                          text=True, check=False, timeout=60)  # fmt: skip


def _log(log: Path) -> str:
    return log.read_text(encoding="utf-8") if log.exists() else ""


BOOT_DOCKER = "Bootstrap: docker\nFrom: mambaorg/micromamba@sha256:" + "f" * 64 + "\n"
BOOT_LOCAL = "Bootstrap: localimage\nFrom: base.sif\n"

REFUSED_ARGVS = [
    ["apptainer", "pull", "--disable-cache", "x.img", f"docker://{REF}"],
    ["apptainer", "--debug", "pull", "x.img", f"docker://{REF}"],
    ["apptainer", "pull", "x.sif", "library://alpine"],
    ["singularity", "pull", "x.sif", f"docker://{REF}"],
    ["apptainer", "build", "x.sif", f"docker://{REF}"],
    ["apptainer", "build", "--fakeroot", "x.sif", "remote.def"],
    ["singularity", "build", "x.sif", "remote.def"],
    ["apptainer", "exec", f"docker://{REF}", "true"],
    ["apptainer", "run", "oras://ghcr.io/x/y:1"],
    ["docker", "pull", REF],
    ["docker", "build", "-f", "ctx/Dockerfile", "-t", "localhost/x:1", "ctx"],
    ["docker", "buildx", "build", "ctx"],
    ["docker", "image", "pull", REF],
    ["docker", "--context", "default", "pull", REF],
    ["docker", "run", "--rm", REF, "true"],
    ["docker", "container", "create", REF],
    ["docker", "compose", "up"],
]

ALLOWED_ARGVS = [
    ["apptainer", "--version"],
    ["singularity", "version"],
    ["apptainer", "exec", "--cleanenv", "local.sif", "true"],
    ["apptainer", "build", "x.sif", "local.def"],
    ["docker", "info"],
    ["docker", "image", "inspect", "--format", "{{.Id}}", "localhost/x:1"],
    ["docker", "run", "--rm", "--pull=never", "localhost/x:1", "true"],
    ["docker", "run", "--pull", "never", "localhost/x:1", "true"],
]


def _defs(cwd: Path) -> None:
    (cwd / "remote.def").write_text(BOOT_DOCKER, encoding="utf-8")
    (cwd / "local.def").write_text(BOOT_LOCAL, encoding="utf-8")


@pytest.mark.parametrize("argv", REFUSED_ARGVS, ids=lambda a: " ".join(a)[:60])
def test_guard_script_refuses_network_pulls(
    guarded: tuple[Path, Path], tmp_path: Path, argv: list[str]
) -> None:
    guard, log = guarded
    _defs(tmp_path)
    got = _run(guard, argv, tmp_path)
    assert got.returncode == REFUSED_EXIT, (got.stdout, got.stderr)
    assert REFUSED in got.stderr
    assert "real_container_pull" in got.stderr  # says how to opt out
    assert _log(log) == ""  # the real binary never ran


@pytest.mark.parametrize("argv", ALLOWED_ARGVS, ids=lambda a: " ".join(a)[:60])
def test_guard_script_passes_other_commands_to_the_real_binary(
    guarded: tuple[Path, Path], tmp_path: Path, argv: list[str]
) -> None:
    guard, log = guarded
    _defs(tmp_path)
    got = _run(guard, argv, tmp_path)
    assert got.returncode == 0, got.stderr
    assert _log(log).splitlines() == [" ".join(argv)]


def test_guard_script_without_a_real_binary_exits_127(tmp_path: Path) -> None:
    guard = write_guard(tmp_path / "guard", {"apptainer": None})
    got = _run(guard, ["apptainer", "--version"], tmp_path)
    assert got.returncode == 127
    got = _run(guard, ["apptainer", "pull", "x.img", f"docker://{REF}"], tmp_path)
    assert got.returncode == REFUSED_EXIT and REFUSED in got.stderr


def test_real_binaries_skips_an_outer_guard(guarded: tuple[Path, Path], tmp_path: Path) -> None:
    guard, _log_path = guarded
    path = os.pathsep.join([str(guard), str(tmp_path / "real")])
    assert real_binaries(path)["apptainer"] == str(guard / "apptainer")
    assert real_binaries(path, (str(guard),))["apptainer"] == str(tmp_path / "real" / "apptainer")


# ---------------------------------------------------------------- the autouse fixture


TOY = """
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from gippa.setupdb import images, uniop_build

REF = {ref!r}
REFUSED = {refused!r}
GUARD_ENV = {guard_env!r}


def test_setup_apptainer_pull_is_refused(tmp_path):
    rec = images._pull("apptainer", "toy", REF, "1", tmp_path / "cache")
    assert rec.status == "absent", rec
    assert REFUSED in rec.message, rec.message
    assert not any((tmp_path / "cache").iterdir())


def test_setup_docker_pull_is_refused(tmp_path):
    rec = images._pull("docker", "toy", REF, "1", tmp_path)
    assert rec.status == "absent", rec
    assert REFUSED in rec.message, rec.message


def test_uniop_apptainer_build_is_refused(tmp_path):
    (tmp_path / "uniop.def").write_text({boot!r})
    digest, why = uniop_build._apptainer_build(tmp_path, tmp_path / "uniop.sif")
    assert digest == ""
    assert REFUSED in why, why


def test_which_resolves_to_the_guard():
    guard = Path(os.environ[GUARD_ENV])
    for name in ("apptainer", "singularity", "docker"):
        assert Path(shutil.which(name)).parent == guard


def test_own_fake_first_on_path_wins(tmp_path, monkeypatch):
    own = tmp_path / "own"
    own.mkdir()
    exe = own / "apptainer"
    exe.write_text("#!/bin/sh\\necho own-fake\\nexit 0\\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{{own}}{{os.pathsep}}{{os.environ['PATH']}}")
    got = subprocess.run(["apptainer", "pull", "x.img", "docker://" + REF],
                         capture_output=True, text=True, check=False)
    assert got.returncode == 0 and got.stdout.strip() == "own-fake"


@pytest.mark.real_container_pull
def test_opt_out_runs_the_real_binary(tmp_path):
    assert GUARD_ENV not in os.environ
    got = subprocess.run(["apptainer", "pull", "x.img", "docker://" + REF],
                         capture_output=True, text=True, check=False)
    assert got.returncode == 0  # the recorder ran
"""


def test_autouse_fixture_guards_setup_in_every_test(tmp_path: Path) -> None:
    """A child pytest under the repository's conftests, with recorder runtimes first on PATH."""
    log = tmp_path / "real.log"
    real_dir = tmp_path / "real"
    _recorders(real_dir, log)
    outer = os.environ.get(GUARD_ENV)
    rest = [d for d in os.environ.get("PATH", "").split(os.pathsep) if d and d != outer]
    env = {k: v for k, v in os.environ.items() if k != GUARD_ENV}
    env["PATH"] = os.pathsep.join([str(real_dir), *rest])
    src = TOY.format(ref=REF, refused=REFUSED, guard_env=GUARD_ENV, boot=BOOT_DOCKER)
    toy = REPO / "tests" / "unit" / f"test_zz_container_guard_toy_{os.getpid()}_{tmp_path.name}.py"
    toy.write_text(src, encoding="utf-8")
    try:
        run = subprocess.run(
            [sys.executable, "-m", "pytest", str(toy), "-q", "-p", "no:cacheprovider",
             "-p", "no:xdist", "-p", "no:randomly", "-rA",
             "--basetemp", str(tmp_path / "child")],
            cwd=REPO, env=env, capture_output=True, text=True, timeout=600,
        )  # fmt: skip
    finally:
        toy.unlink(missing_ok=True)
    out = run.stdout + run.stderr
    assert run.returncode == 0, out
    assert "6 passed" in out, out
    # only the opt-out test reached the recorder, once
    assert _log(log).splitlines() == [f"apptainer pull x.img docker://{REF}"], out
