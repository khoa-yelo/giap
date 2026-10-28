"""The suite-wide guard against real container pulls (tests/conftest.py ``_no_container_pulls``).

Unit tests of ``gippa setup`` once ran the host's ``apptainer pull`` (about 6 GB per image) and
filled ``/tmp`` (#437, #441). The guard puts fake ``apptainer``/``singularity``/``docker`` first
on ``PATH`` that refuse network pulls and builds and hand everything else to the real binary.

Nothing here can reach a registry: the "real" runtimes are recorder scripts that log their argv
and exit 0. Two layers are tested:

- the guard scripts themselves (``tests/helpers/container_guard.py``), per refused and allowed
  command line;
- ``guarded_path``, which keeps the guard on a subprocess's own minimal ``PATH`` (the
  ``gippa compare`` helpers of M10-04 and M10-10);
- the autouse fixture, in a child pytest whose ``PATH`` starts with the recorders: the real
  ``gippa.setupdb.images`` pull and ``gippa.setupdb.uniop_build`` build code get the refusal and
  the recorders never run; a test's own fake first on ``PATH`` wins; ``real_container_pull``
  opts out.
"""

from __future__ import annotations

import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from helpers.container_guard import (
    GUARD_ENV,
    NAMES,
    REFUSED,
    REFUSED_EXIT,
    guarded_path,
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
                          stdin=subprocess.DEVNULL, text=True, check=False,
                          timeout=60)  # fmt: skip


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
    ["docker", "compose", "-f", "c.yml", "up"],
    ["docker", "compose", "-p", "p", "up"],
    ["docker", "compose", "--project-name", "p", "pull"],
    ["apptainer", "--config", "x.conf", "pull", "x.sif", "ubuntu"],
    ["apptainer", "-c", "x.conf", "pull", "x.sif", "ubuntu"],
    ["singularity", "exec", f"docker://{REF}", "true"],
    ["apptainer", "build", "x.sif", "spaced.def"],
    # N2: remote schemes in any case; a definition on stdin; the image after value flags
    ["apptainer", "pull", "x.sif", f"DOCKER://{REF}"],
    ["apptainer", "exec", f"DOCKER://{REF}", "true"],
    ["apptainer", "exec", f"Docker://{REF}", "true"],
    ["apptainer", "run", "ORAS://ghcr.io/x/y:1"],
    ["apptainer", "exec", "LIBRARY://alpine", "true"],
    ["apptainer", "exec", "HTTPS://example.org/x.sif", "true"],
    ["apptainer", "build", "x.sif", f"DOCKER://{REF}"],
    ["apptainer", "build", "out.sif", "-"],
    ["apptainer", "build", "--fakeroot", "out.sif", "-"],
    ["singularity", "build", "out.sif", "-"],
    ["apptainer", "exec", "--env", "URL=https://example.org", f"docker://{REF}", "true"],
    ["apptainer", "exec", "--bind", "bound", f"docker://{REF}", "true"],  # bound/ exists
    ["apptainer", "exec", "-B", "bound", f"docker://{REF}", "true"],
    ["apptainer", "exec", "--home", "bound", f"docker://{REF}", "true"],
    ["apptainer", "exec", "--cleanenv", f"docker://{REF}", "local.sif"],
    ["apptainer", "instance", "start", f"docker://{REF}", "inst"],
    # no word is a local path (the image is missing): every word is scanned (a false refusal)
    ["apptainer", "exec", "missing.sif", "curl", "https://example.org"],
    # N3: docker import of a URL, plugin/service/stack, --pull=never only before the image
    ["docker", "import", "https://example.org/rootfs.tar.gz"],
    ["docker", "import", "HTTPS://example.org/rootfs.tar.gz", "localhost/x:1"],
    ["docker", "import", "-c", "ENV A=1", "https://example.org/rootfs.tar.gz"],
    ["docker", "image", "import", "https://example.org/rootfs.tar.gz"],
    ["docker", "plugin", "install", "vieux/sshfs"],
    ["docker", "plugin", "install", "--grant-all-permissions", "vieux/sshfs"],
    ["docker", "plugin", "upgrade", "vieux/sshfs"],
    ["docker", "service", "create", "--name", "s", REF],
    ["docker", "service", "update", "--image", REF, "s"],
    ["docker", "stack", "deploy", "-c", "stack.yml", "s"],
    ["docker", "stack", "up", "-c", "stack.yml", "s"],
    ["docker", "--context", "x", "service", "create", REF],
    ["docker", "run", "a", "sh", "-c", "--pull=never"],
    ["docker", "run", "--rm", "localhost/x:1", "--pull", "never"],
    ["docker", "container", "run", "localhost/x:1", "--pull=never"],
    ["docker", "create", "localhost/x:1", "--pull=never"],
    ["docker", "run", "--pull=always", "localhost/x:1"],
    ["docker", "run", "--pull", "missing", "localhost/x:1", "never"],
    # an unknown value flag: its value is taken for the image, so --pull=never after it is
    # missed (a refusal, never a pull)
    ["docker", "run", "--made-up-flag", "v", "--pull=never", "localhost/x:1"],
    # #470 review S1: an unknown flag or its value never ends the scan (fail-safe)
    ["apptainer", "exec", "--docker-host", "unix:///d.sock", "docker://alpine", "true"],
    ["apptainer", "shell", "-s", "/bin/sh", "docker://alpine"],
    ["apptainer", "exec", "--cwd", "/tmp", "docker://alpine", "true"],
    ["apptainer", "exec", "--no-mount", "/etc/hosts", "docker://alpine", "true"],
    ["apptainer", "instance", "start", "--dmtcp-launch", "/tmp", "docker://alpine", "i1"],
    ["apptainer", "instance", "start", "docker://alpine", "i1"],  # cwd holds a file "start"
    ["apptainer", "exec", "bound", "docker://alpine", "true"],  # a plain directory, no sandbox
    ["apptainer", "exec", "local.def", "docker://alpine", "true"],  # a file, not *.sif/*.img
    # S2: docker takes the last --pull value
    ["docker", "run", "--pull=never", "--pull=always", "alpine"],
    ["docker", "run", "--pull", "never", "--pull", "missing", "alpine"],
    ["docker", "run", "--pull=never", "--pull", "missing", "alpine"],
    ["docker", "container", "create", "--pull", "never", "--pull=always", "alpine"],
    # N1: every Bootstrap: line of a multi-stage definition
    ["apptainer", "build", "x.sif", "multi.def"],
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
    ["docker", "container", "run", "--pull=never", "localhost/x:1", "true"],
    ["docker", "compose", "version"],
    ["apptainer", "--config", "x.conf", "exec", "local.sif", "true"],
    # N2: a flag's value and the container's own command are not pulls
    ["apptainer", "exec", "--env", "URL=https://example.org", "local.sif", "true"],
    ["apptainer", "exec", "--env", "https://example.org", "local.sif", "true"],
    ["apptainer", "exec", "--env=URL=https://example.org", "local.sif", "true"],
    ["apptainer", "exec", "local.sif", "curl", "https://example.org"],
    ["apptainer", "exec", "local.sif", "sh", "-c", "apptainer pull x docker://y"],
    ["apptainer", "run", "--bind", "bound:/b", "local.sif", "wget", "HTTPS://example.org"],
    ["apptainer", "exec", "sandbox", "curl", "https://example.org"],  # a sandbox directory
    ["apptainer", "exec", "local.img", "curl", "https://example.org"],
    ["apptainer", "exec", "--bind", "bound", "local.sif", "true"],
    ["apptainer", "exec", "instance://inst", "curl", "https://example.org"],
    ["apptainer", "instance", "start", "local.sif", "inst"],
    ["apptainer", "exec", "docker-daemon://localhost/x:1", "true"],
    ["apptainer", "build", "out.sif", "local.def"],
    # N3: local import, read-only plugin/service/stack commands, --pull=never before the image
    ["docker", "import", "rootfs.tar", "localhost/x:1"],
    ["docker", "import", "-", "localhost/x:1"],
    ["docker", "image", "import", "rootfs.tar"],
    ["docker", "plugin", "ls"],
    ["docker", "service", "ls"],
    ["docker", "stack", "ls"],
    ["docker", "run", "--pull=never", "localhost/x:1", "sh", "-c", "--pull=always"],
    ["docker", "run", "-e", "A=1", "--name", "n", "--pull", "never", "localhost/x:1"],
    ["docker", "run", "-v", "/a:/b", "-w", "/b", "--pull=never", "localhost/x:1", "true"],
    ["docker", "run", "-it", "--rm", "--pull=never", "localhost/x:1"],
    ["docker", "--context", "x", "run", "--pull=never", "localhost/x:1"],
    ["docker", "-c", "x", "container", "create", "--pull", "never", "localhost/x:1"],
    ["docker", "create", "--entrypoint", "sh", "--pull=never", "localhost/x:1"],
    ["docker", "run", "--pull=always", "--pull=never", "alpine"],
    ["docker", "run", "--pull", "missing", "--pull", "never", "alpine"],
    ["apptainer", "build", "x.sif", "multi_local.def"],
]


def _defs(cwd: Path) -> None:
    (cwd / "remote.def").write_text(BOOT_DOCKER, encoding="utf-8")
    (cwd / "local.def").write_text(BOOT_LOCAL, encoding="utf-8")
    (cwd / "spaced.def").write_text("  bootstrap : library\nFrom: alpine\n", encoding="utf-8")
    (cwd / "local.sif").write_bytes(b"")
    (cwd / "rootfs.tar").write_bytes(b"")
    (cwd / "bound").mkdir(exist_ok=True)
    (cwd / "sandbox" / ".singularity.d").mkdir(parents=True, exist_ok=True)
    (cwd / "local.img").write_bytes(b"")
    (cwd / "start").write_bytes(b"")  # a stray cwd file named like a positional word
    (cwd / "multi.def").write_text(BOOT_LOCAL + "\n" + BOOT_DOCKER, encoding="utf-8")
    (cwd / "multi_local.def").write_text(
        "Bootstrap: localimage\nFrom: a.sif\nStage: one\n\nBootStrap: scratch\n",
        encoding="utf-8",
    )


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


@pytest.mark.needs_apptainer
def test_needs_apptainer_keeps_the_real_runtime():
    assert GUARD_ENV not in os.environ
    assert Path(shutil.which("apptainer")).parent == Path({real_dir!r})


@pytest.mark.needs_docker
def test_needs_docker_keeps_the_real_runtime():
    assert GUARD_ENV not in os.environ
    assert Path(shutil.which("docker")).parent == Path({real_dir!r})


@pytest.mark.e2e
def test_e2e_keeps_the_real_runtime():
    assert GUARD_ENV not in os.environ
    assert Path(shutil.which("apptainer")).parent == Path({real_dir!r})


@pytest.mark.real_container_pull
def test_opt_out_runs_the_real_binary(tmp_path):
    assert GUARD_ENV not in os.environ
    got = subprocess.run(["apptainer", "pull", "x.img", "docker://" + REF],
                         capture_output=True, text=True, check=False)
    assert got.returncode == 0  # the recorder ran
"""

TOY_E2E_DIR = """
import os
import shutil
from pathlib import Path


def test_tests_e2e_keeps_the_real_runtime():
    assert {guard_env!r} not in os.environ
    assert Path(shutil.which("apptainer")).parent == Path({real_dir!r})
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
    env["GIPPA_TEST_CONTAINER_GUARD_TOY"] = "1"  # tests/conftest.py ignores the toys otherwise
    fmt = {"ref": REF, "refused": REFUSED, "guard_env": GUARD_ENV, "boot": BOOT_DOCKER,
           "real_dir": str(real_dir)}  # fmt: skip
    name = f"test_zz_container_guard_toy_{os.getpid()}_{tmp_path.name}.py"
    e2e_name = name.replace("_toy_", "_toy_e2e_")  # unique basenames (rootdir-less imports)
    toys = [REPO / "tests" / "unit" / name, REPO / "tests" / "e2e" / e2e_name]
    toys[0].write_text(TOY.format(**fmt), encoding="utf-8")
    toys[1].write_text(TOY_E2E_DIR.format(**fmt), encoding="utf-8")
    try:
        run = subprocess.run(
            [sys.executable, "-m", "pytest", *map(str, toys), "-q", "-p", "no:cacheprovider",
             "-p", "no:xdist", "-p", "no:randomly", "-rA",
             "--basetemp", str(tmp_path / "child")],
            cwd=REPO, env=env, capture_output=True, text=True, timeout=600,
        )  # fmt: skip
    finally:
        for toy in toys:
            toy.unlink(missing_ok=True)
    out = run.stdout + run.stderr
    assert run.returncode == 0, out
    assert "10 passed" in out, out
    # only the opt-out test's pull reached the recorder, once; besides it only the capability
    # probes of the needs_apptainer/needs_docker tests ran (tests/conftest.py CHECKS)
    lines = _log(log).splitlines()
    pull = f"apptainer pull x.img docker://{REF}"
    assert lines.count(pull) == 1, lines
    assert set(lines) <= {pull, "apptainer --version", "docker info"}, lines


def test_left_over_toys_are_never_collected(tmp_path: Path) -> None:
    """A toy left by a killed run is ignored without the child's env var (and gitignored)."""
    # a small directory: an explicit file argument would bypass collect_ignore_glob
    toy = REPO / "tests" / "unit" / "compare" / f"test_zz_container_guard_toy_left_{os.getpid()}.py"
    toy.write_text("def test_should_not_run():\n    raise AssertionError\n", encoding="utf-8")
    env = {k: v for k, v in os.environ.items() if k != "GIPPA_TEST_CONTAINER_GUARD_TOY"}
    try:
        run = subprocess.run(
            [sys.executable, "-m", "pytest", str(toy.parent), "-q", "-p", "no:cacheprovider",
             "-p", "no:xdist", "--collect-only", "-k", "should_not_run",
             "--basetemp", str(tmp_path / "child")],
            cwd=REPO, env=env, capture_output=True, text=True, timeout=600,
        )  # fmt: skip
        ignored = subprocess.run(["git", "check-ignore", "-q", str(toy)], cwd=REPO, check=False)
    finally:
        toy.unlink(missing_ok=True)
    assert run.returncode == 5, run.stdout + run.stderr  # 5: nothing collected (-k)
    assert "should_not_run" not in run.stdout, run.stdout
    assert ignored.returncode == 0  # .gitignore covers the pattern


# ---------------------------------------------------------------- N4: a subprocess's own PATH


def test_guarded_path_is_unchanged_when_the_path_holds_no_runtime(tmp_path: Path) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    assert guarded_path(str(empty)) == str(empty)


def test_guarded_path_guards_exactly_the_runtimes_of_that_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The scripts forward to the subprocess PATH's own binaries, so ``which`` finds the same
    runtimes as without the guard, and refuse pulls."""
    monkeypatch.delenv(GUARD_ENV, raising=False)  # independent of the per-test guard
    own_log = tmp_path / "own.log"
    own = tmp_path / "own"
    _recorders(own, own_log)
    (own / "singularity").unlink()
    (own / "docker").unlink()
    path = os.pathsep.join([str(own), "/usr/bin", "/bin"])
    got = guarded_path(path)
    first, rest = got.split(os.pathsep, 1)
    assert rest == path
    assert sorted(p.name for p in Path(first).iterdir()) == ["apptainer"]
    assert guarded_path(path) == got  # cached per path
    which = shutil.which("apptainer", path=got)
    assert which is not None and Path(which).parent == Path(first)
    pull = _run(Path(first), ["apptainer", "pull", "x.img", f"docker://{REF}"], tmp_path)
    assert pull.returncode == REFUSED_EXIT and REFUSED in pull.stderr
    ok = _run(Path(first), ["apptainer", "--version"], tmp_path)
    assert ok.returncode == 0, ok.stderr
    assert _log(own_log).splitlines() == ["apptainer --version"]


def test_guarded_path_writes_under_the_session_directory(tmp_path: Path) -> None:
    """tests/conftest.py publishes a session directory before any module-scoped fixture."""
    from helpers import container_guard

    base = container_guard._session_dir
    assert base is not None and base.is_dir()
    own = tmp_path / "own"
    _recorders(own, tmp_path / "own.log")
    first = Path(guarded_path(str(own)).split(os.pathsep)[0])
    assert first.parent == base


def test_write_guard_replaces_scripts_atomically(tmp_path: Path) -> None:
    bindir = write_guard(tmp_path / "g", {"apptainer": "/a", "docker": "/d"})
    old = (bindir / "apptainer").stat().st_ino
    write_guard(bindir, {"apptainer": "/b"})
    assert (bindir / "apptainer").stat().st_ino != old  # a new file renamed over the old one
    assert "real=/b" in (bindir / "apptainer").read_text(encoding="utf-8")
    assert os.access(bindir / "apptainer", os.X_OK)
    assert sorted(p.name for p in bindir.iterdir()) == ["apptainer", "docker"]  # no temp left


def _helper_compare_env(module: str, tmp: Path, mp: pytest.MonkeyPatch) -> dict[str, str]:
    """The env a helper's ``gippa compare`` subprocess gets, with recorders as its MIN_PATH."""
    mod = importlib.import_module(module)
    own = tmp / "own"
    _recorders(own, tmp / "own.log")
    mp.setattr(mod, "MIN_PATH", str(own))
    seen: dict[str, dict[str, str]] = {}

    def fake_run(cmd: list[str], **kw: object) -> subprocess.CompletedProcess[str]:
        seen["env"] = kw["env"]  # type: ignore[assignment]
        return subprocess.CompletedProcess(cmd, 0, "", "")

    runs = [tmp / "a", tmp / "b"]
    with mp.context() as m:  # subprocess.run is the real one again afterwards
        m.setattr(mod.subprocess, "run", fake_run)
        if module.endswith("m10_04"):
            mod.gippa_compare(*runs, tmp / "cmp", env_root=tmp / "env")
        else:
            mod.gippa_compare_n(runs, ["A1", "A2"], tmp / "cmp", env_root=tmp / "env")
    return seen["env"]


HELPERS = ["helpers.pilot_runs_m10_04", "helpers.a3_run_m10_10"]


@pytest.fixture(scope="module")
def module_scoped_envs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, tuple]:
    """Like ``compare_pair``/``compare3``: a module-scoped fixture runs before the function-scoped
    autouse guard, so ``GUARD_ENV`` is not set yet."""
    out = {}
    with pytest.MonkeyPatch.context() as mp:
        for module in HELPERS:
            tmp = tmp_path_factory.mktemp("helper_env")
            out[module] = (os.environ.get(GUARD_ENV), tmp, _helper_compare_env(module, tmp, mp))
    return out


@pytest.mark.parametrize("module", HELPERS)
def test_helpers_called_from_a_module_scoped_fixture_keep_the_guard(
    module: str, module_scoped_envs: dict[str, tuple]
) -> None:
    guard_env, tmp, env = module_scoped_envs[module]
    assert guard_env is None  # the per-test guard was not active when the helper ran
    first, rest = env["PATH"].split(os.pathsep, 1)
    assert rest == str(tmp / "own")
    pull = _run(Path(first), ["apptainer", "pull", "x.img", f"docker://{REF}"], tmp)
    assert pull.returncode == REFUSED_EXIT and REFUSED in pull.stderr
    assert _log(tmp / "own.log") == ""


@pytest.mark.parametrize("module", HELPERS)
def test_helpers_with_a_minimal_path_keep_the_guard(
    module: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    env = _helper_compare_env(module, tmp_path, monkeypatch)
    first, rest = env["PATH"].split(os.pathsep, 1)
    assert rest == str(tmp_path / "own")
    pull = _run(Path(first), ["apptainer", "pull", "x.img", f"docker://{REF}"], tmp_path)
    assert pull.returncode == REFUSED_EXIT and REFUSED in pull.stderr
    assert _log(tmp_path / "own.log") == ""


@pytest.mark.parametrize("module", HELPERS)
def test_helpers_minimal_path_is_the_system_path(module: str) -> None:
    assert importlib.import_module(module).MIN_PATH == "/usr/bin:/bin"
