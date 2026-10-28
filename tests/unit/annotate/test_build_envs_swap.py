"""build_envs.sh keeps the old env until the new one verifies, and CUDA solves on CPU nodes.

Backlog (build_envs): job 1888853 lost the cascadia env and its locks because the old env was
removed before the new one solved; job 1901404 needed ``CONDA_OVERRIDE_CUDA`` by hand. No conda is
needed: ``micromamba``, ``git`` and ``nvidia-smi`` are shell stubs in ``$GIPPA_ROOT/bin`` (first on
PATH via ``farmshare_env.sh``). The stub ``create`` writes a half-built prefix before it fails, so a
restored env is told apart from a new one by its contents. Nothing is ever solved.
"""

from __future__ import annotations

import contextlib
import os
import shutil
import signal
import subprocess
import sys
import time
from collections.abc import Iterator
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]

FAKE_MICROMAMBA = r"""#!/bin/bash
echo "micromamba $* CUDA=${CONDA_OVERRIDE_CUDA-unset}" >> "$FAKE_DIR/calls.log"
envs="$MAMBA_ROOT_PREFIX/envs"
name_of() { while [ "$#" -gt 0 ]; do if [ "$1" = "-n" ]; then echo "$2"; return; fi; shift; done; }
case "$1" in
  env)
    case "$2" in
      export)
        n="$(name_of "$@")"
        [ -e "$FAKE_DIR/fail_export_$n" ] && exit 1
        echo "@EXPLICIT"
        echo "https://conda.anaconda.org/conda-forge/linux-64/new-2.0-0.conda#$n"
        exit 0 ;;
    esac
    exit 0 ;;
  create)
    n="$(name_of "$@")"
    mkdir -p "$envs/$n/conda-meta" "$envs/$n/bin"
    echo new > "$envs/$n/bin/payload"
    if [ -e "$FAKE_DIR/hold_create_$n" ] && mkdir "$FAKE_DIR/held" 2>/dev/null; then
      touch "$FAKE_DIR/in_create"
      # Bounded: a test that fails before it writes the release file must not leave this
      # loop running (two such children once ran 7.7 h after their test had ended).
      end=$((SECONDS + 120))
      while [ ! -e "$FAKE_DIR/release" ]; do
        [ "$SECONDS" -lt "$end" ] || { echo "fake create: never released" >&2; exit 1; }
        sleep 0.05
      done
      exit 0
    fi
    if [ -e "$FAKE_DIR/hang_create_$n" ]; then
      touch "$FAKE_DIR/in_create"
      sleep 60 &
      wait
      exit 1
    fi
    [ -e "$FAKE_DIR/fail_create_$n" ] && exit 1
    exit 0 ;;
  run)
    shift
    n=""
    while [ "$#" -gt 0 ]; do
      case "$1" in -n) n="$2"; shift 2 ;; -*) shift ;; *) break ;; esac
    done
    [ -d "$envs/$n" ] || { echo "no env $n" >&2; exit 1; }
    case " $* " in
      *"pip freeze"*)
        [ -e "$FAKE_DIR/fail_pip_$n" ] && exit 1
        echo "newpkg==2.0"
        exit 0 ;;
    esac
    [ -e "$FAKE_DIR/fail_version_$n" ] && exit 1
    cat "$FAKE_DIR/version_$n" ;;
  clean) exit 0 ;;
esac
exit 0
"""

FAKE_GIT = r"""#!/bin/bash
for a in "$@"; do [ -e "$FAKE_DIR/fail_git_$a" ] && exit 1; done
if [ "$1" = "clone" ]; then for a in "$@"; do last="$a"; done; mkdir -p "$last"; fi
exit 0
"""

FAKE_NVIDIA_SMI = r"""#!/bin/bash
[ -e "$FAKE_DIR/gpu" ] || { echo "NVIDIA-SMI has failed" >&2; exit 9; }
echo "GPU 0: NVIDIA L40S (UUID: GPU-00000000-0000-0000-0000-000000000000)"
"""

OLD_LOCK = "@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/old-1.0-0.conda\n"
OLD_PIP = "oldpkg==1.0\n"


class Box:
    def __init__(self, tmp: Path) -> None:
        self.repo, self.root, self.fake = tmp / "repo", tmp / "root", tmp / "fake"
        self.locks = self.repo / "workflow" / "envs" / "locks"
        for d in (self.locks, self.repo / "workflow" / "conf", self.fake, self.root / "bin"):
            d.mkdir(parents=True)
        (self.repo / "src").symlink_to(REPO / "src")
        (self.repo / "scripts").symlink_to(REPO / "scripts")
        for yml in (REPO / "workflow" / "envs").glob("*.yml"):
            shutil.copy(yml, self.repo / "workflow" / "envs" / yml.name)
        shutil.copy(
            REPO / "workflow" / "conf" / "containers.config",
            self.repo / "workflow" / "conf" / "containers.config",
        )
        for name, body in (
            ("micromamba", FAKE_MICROMAMBA),
            ("git", FAKE_GIT),
            ("nvidia-smi", FAKE_NVIDIA_SMI),
        ):
            path = self.root / "bin" / name
            path.write_text(body)
            path.chmod(0o755)
        (self.fake / "calls.log").write_text("")
        self.envs = self.root / "micromamba" / "envs"

    def prefix(self, tool: str) -> Path:
        return self.envs / f"gippa-{tool}"

    def version(self, tool: str, text: str) -> None:
        (self.fake / f"version_gippa-{tool}").write_text(text + "\n")

    def flag(self, what: str, tool: str = "") -> None:
        (self.fake / (f"{what}_gippa-{tool}" if tool else what)).write_text("")

    def old_env(self, tool: str, marker: str = "0" * 64) -> None:
        """A built env with a stale marker (a rebuild is due) and its committed locks."""
        prefix = self.prefix(tool)
        (prefix / "bin").mkdir(parents=True)
        (prefix / "bin" / "payload").write_text("old\n")
        (prefix / ".gippa-built").write_text(marker + "\n")
        (self.locks / f"{tool}.explicit.txt").write_text(OLD_LOCK)
        (self.locks / f"{tool}.pip.txt").write_text(OLD_PIP)

    def env(self, **extra: str) -> dict[str, str]:
        env = {k: v for k, v in os.environ.items() if k in {"PATH", "LANG", "TMPDIR"}}
        env.update(
            GIPPA_REPO=str(self.repo),
            GIPPA_ROOT=str(self.root),
            GIPPA_PYTHON=sys.executable,
            FAKE_DIR=str(self.fake),
            HOME=str(self.root / "home"),
        )
        env.pop("CONDA_OVERRIDE_CUDA", None)
        env.update(extra)
        return env

    def cmd(self, *tools: str) -> list[str]:
        return ["bash", str(self.repo / "scripts" / "build_envs.sh"), *tools]

    def run(self, *tools: str, **extra: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            self.cmd(*tools),
            capture_output=True,
            text=True,
            cwd=self.repo,
            env=self.env(**extra),
            timeout=300,
            check=False,
        )

    def creates(self) -> list[str]:
        lines = (self.fake / "calls.log").read_text().splitlines()
        return [c for c in lines if c.startswith("micromamba create")]

    def yml_sha(self, tool: str) -> str:
        import hashlib

        data = (self.repo / "workflow" / "envs" / f"{tool}.yml").read_bytes()
        return hashlib.sha256(data).hexdigest()


def _survivors(fake: Path) -> list[int]:
    """PIDs of live processes started for this box (they inherit ``FAKE_DIR=<fake>``)."""
    needle = f"FAKE_DIR={fake}".encode()
    pids = []
    for proc in Path("/proc").glob("[0-9]*"):
        try:
            environ = (proc / "environ").read_bytes()
        except OSError:  # gone, a zombie or not ours
            continue
        if needle in environ.split(b"\0"):
            pids.append(int(proc.name))
    return pids


@pytest.fixture
def box(tmp_path: Path) -> Iterator[Box]:
    b = Box(tmp_path)
    b.version("bakta", "bakta 1.12.2")
    yield b
    # No child of the test may outlive it (a held fake create once ran 7.7 h).
    deadline = time.monotonic() + 10
    while (left := _survivors(b.fake)) and time.monotonic() < deadline:
        time.sleep(0.1)
    for pid in left:
        with contextlib.suppress(ProcessLookupError):
            os.kill(pid, signal.SIGKILL)
    assert not left, f"processes outlived the test: {left}"


def _assert_old_kept(box: Box, tool: str) -> None:
    prefix = box.prefix(tool)
    assert (prefix / "bin" / "payload").read_text() == "old\n", "the old env was not restored"
    assert (prefix / ".gippa-built").read_text() == "0" * 64 + "\n"
    assert not (prefix / "conda-meta").exists(), "half-built new env left in place"
    assert not Path(f"{prefix}.gippa-old").exists()
    assert (box.locks / f"{tool}.explicit.txt").read_text() == OLD_LOCK
    assert (box.locks / f"{tool}.pip.txt").read_text() == OLD_PIP
    assert not list(box.locks.glob(f"{tool}.*.tmp"))


# ---------------------------------------------------------------- the old env survives a failure


@pytest.mark.parametrize(
    "step",
    ["create", "version", "pip", "export"],
)
def test_failed_rebuild_keeps_old_env_and_locks(box: Box, step: str) -> None:
    box.old_env("bakta")
    box.flag(f"fail_{step}", "bakta")
    proc = box.run("bakta")
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "FAILED: bakta" in proc.stderr
    _assert_old_kept(box, "bakta")


def test_version_mismatch_keeps_old_env_and_locks(box: Box) -> None:
    box.old_env("bakta")
    box.version("bakta", "bakta 1.11.0")
    assert box.run("bakta").returncode == 1
    _assert_old_kept(box, "bakta")


def test_failed_tool_does_not_touch_other_tool(box: Box) -> None:
    box.old_env("bakta")
    box.version("miniprot", "0.18")
    box.flag("fail_create", "bakta")
    proc = box.run("bakta", "miniprot")
    assert proc.returncode == 1
    _assert_old_kept(box, "bakta")
    assert (box.prefix("miniprot") / ".gippa-built").read_text().strip() == box.yml_sha("miniprot")


# ---------------------------------------------------------------- a verified rebuild swaps


def test_successful_rebuild_replaces_old_env_and_locks(box: Box) -> None:
    box.old_env("bakta")
    proc = box.run("bakta")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    prefix = box.prefix("bakta")
    assert (prefix / "bin" / "payload").read_text() == "new\n"
    assert (prefix / ".gippa-built").read_text().strip() == box.yml_sha("bakta")
    assert not Path(f"{prefix}.gippa-old").exists()
    assert "new-2.0-0.conda" in (box.locks / "bakta.explicit.txt").read_text()
    assert (box.locks / "bakta.pip.txt").read_text() == "newpkg==2.0\n"


def test_create_runs_at_the_real_prefix(box: Box) -> None:
    """Envs are not relocatable: the new env is built under its final name, never moved in."""
    box.old_env("bakta")
    assert box.run("bakta").returncode == 0
    (create,) = box.creates()
    assert " -n gippa-bakta " in create
    assert " -p " not in create and "--prefix" not in create


def test_first_build_failure_leaves_no_env_and_keeps_committed_locks(box: Box) -> None:
    (box.locks / "bakta.explicit.txt").write_text(OLD_LOCK)
    box.flag("fail_version", "bakta")
    assert box.run("bakta").returncode == 1
    assert not (box.prefix("bakta") / ".gippa-built").exists(), "unverified env marked built"
    assert not Path(f"{box.prefix('bakta')}.gippa-old").exists()
    assert (box.locks / "bakta.explicit.txt").read_text() == OLD_LOCK


# ---------------------------------------------------------------- a killed run


def test_leftover_old_env_is_restored_before_the_next_build(box: Box) -> None:
    """A run killed mid-build left the old env aside and a half-built new prefix."""
    box.old_env("bakta")
    prefix = box.prefix("bakta")
    prefix.rename(f"{prefix}.gippa-old")
    (prefix / "conda-meta").mkdir(parents=True)
    box.flag("fail_create", "bakta")
    assert box.run("bakta").returncode == 1
    _assert_old_kept(box, "bakta")


def test_leftover_old_env_is_dropped_when_the_new_env_finished(box: Box) -> None:
    box.old_env("bakta")
    prefix = box.prefix("bakta")
    prefix.rename(f"{prefix}.gippa-old")
    (prefix / "bin").mkdir(parents=True)
    (prefix / "bin" / "payload").write_text("new\n")
    (prefix / ".gippa-built").write_text(box.yml_sha("bakta") + "\n")
    proc = box.run("bakta")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "gippa-bakta exists" in proc.stdout
    assert box.creates() == []
    assert (prefix / "bin" / "payload").read_text() == "new\n"
    assert not Path(f"{prefix}.gippa-old").exists()


def test_sigterm_during_create_restores_old_env(box: Box) -> None:
    box.old_env("bakta")
    box.flag("hang_create", "bakta")
    proc = subprocess.Popen(
        box.cmd("bakta"),
        cwd=box.repo,
        env=box.env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    try:
        deadline = time.monotonic() + 120
        while not (box.fake / "in_create").exists():
            assert proc.poll() is None, proc.communicate()
            assert time.monotonic() < deadline, "create never started"
            time.sleep(0.05)
        os.killpg(proc.pid, signal.SIGTERM)  # Slurm signals every process of the job
        proc.wait(timeout=60)
    finally:
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
    assert proc.returncode != 0
    _assert_old_kept(box, "bakta")


# ---------------------------------------------------------------- CONDA_OVERRIDE_CUDA


def _cascadia(box: Box) -> None:
    from gippa.annotate.tools import PILOT_TOOLS

    box.version("cascadia", PILOT_TOOLS["cascadia"].version.removeprefix("git:"))


def test_cuda_build_on_cpu_node_sets_override_11_8(box: Box) -> None:
    _cascadia(box)
    proc = box.run("cascadia")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    (create,) = box.creates()
    assert create.endswith("CUDA=11.8")


def test_cuda_override_value_is_overridable(box: Box) -> None:
    _cascadia(box)
    assert box.run("cascadia", CONDA_OVERRIDE_CUDA="12.4").returncode == 0
    (create,) = box.creates()
    assert create.endswith("CUDA=12.4")


def test_cuda_build_on_gpu_node_sets_no_override(box: Box) -> None:
    _cascadia(box)
    box.flag("gpu")
    assert box.run("cascadia").returncode == 0
    (create,) = box.creates()
    assert create.endswith("CUDA=unset")


def test_cpu_env_sets_no_override(box: Box) -> None:
    assert box.run("bakta").returncode == 0
    (create,) = box.creates()
    assert create.endswith("CUDA=unset")


@pytest.mark.parametrize(
    ("line", "wants"),
    [
        ("  - pytorch=2.0.*=cuda*", True),
        ("  - pytorch 2.0.* cuda*", True),
        ("  - pytorch[build=cuda*]", True),
        ("  - pytorch=2.0.*=cpu*", False),
        ("  - cuda-version=11.8", False),
        ("  # - pytorch=2.0.*=cuda*", False),
    ],
)
def test_cuda_build_detection(box: Box, line: str, wants: bool) -> None:
    yml = box.repo / "workflow" / "envs" / "bakta.yml"
    yml.write_text(yml.read_text() + line + "\n")
    proc = box.run("bakta")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    (create,) = box.creates()
    assert create.endswith("CUDA=11.8" if wants else "CUDA=unset"), create


# ---------------------------------------------------------------- #457 review: B1, S1-S4


def _shim(box: Box, name: str, match: str, action: str) -> None:
    """A ``name`` wrapper first on PATH: runs ``action`` once when an argument matches ``match``."""
    real = shutil.which(name)
    assert real, name
    body = f"""#!/bin/bash
for a in "$@"; do
  case "$a" in
    {match})
      if mkdir "$FAKE_DIR/shim_{name}" 2>/dev/null; then {action}; fi ;;
  esac
done
exec {real} "$@"
"""
    path = box.root / "bin" / name
    path.write_text(body)
    path.chmod(0o755)


_HANG = 'touch "$FAKE_DIR/in_shim"; sleep 60 & wait; exit 1'


def _term_at_shim(box: Box, *tools: str) -> subprocess.Popen[bytes]:
    proc = subprocess.Popen(
        box.cmd(*tools),
        cwd=box.repo,
        env=box.env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    try:
        _wait_for(box.fake / "in_shim", proc)
        os.killpg(proc.pid, signal.SIGTERM)
        proc.wait(timeout=60)
    finally:
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
    return proc


def _wait_for(path: Path, proc: subprocess.Popen[bytes]) -> None:
    deadline = time.monotonic() + 120
    while not path.exists():
        assert proc.poll() is None, proc.communicate()
        assert time.monotonic() < deadline, f"never reached {path.name}"
        time.sleep(0.05)


def _assert_new(box: Box, tool: str = "bakta", *, locks: bool = True) -> None:
    prefix = box.prefix(tool)
    assert (prefix / "bin" / "payload").read_text() == "new\n"
    assert (prefix / ".gippa-built").read_text().strip() == box.yml_sha(tool)
    assert not Path(f"{prefix}.gippa-old").exists()
    if locks:
        assert "new-2.0-0.conda" in (box.locks / f"{tool}.explicit.txt").read_text()
        assert (box.locks / f"{tool}.pip.txt").read_text() == "newpkg==2.0\n"
    assert not list(box.locks.glob(f"{tool}.*.prev")) and not list(box.locks.glob(f"{tool}.*.tmp"))


def _leftovers(box: Box, tool: str = "bakta") -> list[str]:
    return sorted(p.name for p in box.envs.glob(f"gippa-{tool}.gippa-*"))


def test_term_while_old_env_is_deleted_keeps_the_new_env(box: Box) -> None:
    """B1: the old env is renamed to trash before rm -rf; a TERM then never restores it."""
    box.old_env("bakta")
    (box.prefix("bakta") / "lib").mkdir()
    (box.prefix("bakta") / "lib" / "libold.so").write_text("x")
    _shim(box, "rm", "*.gippa-trash.*", '/bin/rm -rf "$a/lib"; ' + _HANG)
    proc = _term_at_shim(box, "bakta")
    assert proc.returncode == 143
    _assert_new(box)
    assert not (box.prefix("bakta") / "lib").exists()
    (box.root / "bin" / "rm").unlink()
    assert box.run("bakta").returncode == 0  # the next run sweeps the trash
    assert _leftovers(box) == []
    _assert_new(box)


def test_term_while_locks_move_keeps_new_env_and_heals_locks(box: Box) -> None:
    """B1: the marker is written before the locks move; the next run re-exports the locks."""
    box.old_env("bakta")
    _shim(box, "mv", "*.explicit.txt.tmp", _HANG)
    proc = _term_at_shim(box, "bakta")
    assert proc.returncode == 143
    prefix = box.prefix("bakta")
    assert (prefix / "bin" / "payload").read_text() == "new\n"
    assert (prefix / ".gippa-built").read_text().strip() == box.yml_sha("bakta")
    (box.root / "bin" / "mv").unlink()
    assert box.run("bakta").returncode == 0
    _assert_new(box)
    assert _leftovers(box) == []


@pytest.mark.parametrize("impl", ["flock", "mkdir"])
def test_second_build_of_one_env_waits_out_the_first(box: Box, impl: str) -> None:
    """S1: a second build of an env that is being built fails at once and touches nothing."""
    extra = {"GIPPA_BUILD_LOCK": "mkdir"} if impl == "mkdir" else {}
    box.old_env("bakta")
    box.flag("hold_create", "bakta")
    first = subprocess.Popen(
        box.cmd("bakta"),
        cwd=box.repo,
        env=box.env(**extra),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,  # its fake create child is killed with it in the finally
    )
    try:
        _wait_for(box.fake / "in_create", first)
        box.flag("fail_create", "bakta")  # the second build's solve would fail
        second = box.run("bakta", **extra)
        assert second.returncode == 1
        assert "another build holds the lock" in second.stderr
        assert "restored" not in second.stderr
        old = Path(f"{box.prefix('bakta')}.gippa-old")
        assert (old / "bin" / "payload").read_text() == "old\n", "the first build was disturbed"
        assert (box.prefix("bakta") / "conda-meta").is_dir()
        (box.fake / "fail_create_gippa-bakta").unlink()
        (box.fake / "release").write_text("")
        out, err = first.communicate(timeout=120)
    finally:
        (box.fake / "release").write_text("")  # frees a held create if an assert failed above
        if first.poll() is None:
            with contextlib.suppress(subprocess.TimeoutExpired):
                first.communicate(timeout=30)
        with contextlib.suppress(ProcessLookupError):
            os.killpg(first.pid, signal.SIGKILL)  # the whole group, not only bash
        first.wait()
    assert first.returncode == 0, (out + err).decode()
    _assert_new(box)
    assert box.run("bakta", **extra).returncode == 0, "the lock was not released"


def test_stale_mkdir_lock_of_a_dead_build_is_taken_over(box: Box) -> None:
    lock = box.envs / ".gippa-bakta.gippa-lock.d"
    lock.mkdir(parents=True)
    dead = subprocess.Popen(["true"])
    dead.wait()
    (lock / "pid").write_text(f"{dead.pid}\n")
    proc = box.run("bakta", GIPPA_BUILD_LOCK="mkdir")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert not lock.exists()


def test_leftover_old_env_is_reconciled_without_nesting(box: Box) -> None:
    """S2: with a leftover .gippa-old, a successful rebuild must not nest the env inside it."""
    box.old_env("bakta")
    prefix = box.prefix("bakta")
    prefix.rename(f"{prefix}.gippa-old")
    (prefix / "conda-meta").mkdir(parents=True)  # half-built, no marker
    proc = box.run("bakta")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    _assert_new(box)
    assert not (prefix / "gippa-bakta").exists()
    assert _leftovers(box) == []


def _fgsea(box: Box) -> str:
    lock = (
        "@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/new-2.0-0.conda#gippa-fgsea\n"
    )
    box.version("fgsea", "[1] \u20181.36.2\u2019")  # R prints curly quotes
    prefix = box.prefix("fgsea")
    (prefix / "bin").mkdir(parents=True)
    (prefix / "bin" / "payload").write_text("old\n")
    (prefix / ".gippa-built").write_text("0" * 64 + "\n")
    (box.locks / "fgsea.explicit.txt").write_text(lock)
    return lock


def test_frozen_rebuild_removes_the_old_env(box: Box) -> None:
    """S3: the frozen path retires the old env too, and never rewrites the committed lock."""
    lock = _fgsea(box)
    proc = box.run("fgsea")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    _assert_new(box, "fgsea", locks=False)
    assert _leftovers(box, "fgsea") == []
    assert (box.locks / "fgsea.explicit.txt").read_text() == lock
    assert not (box.locks / "fgsea.pip.txt").exists()


def test_failed_frozen_rebuild_keeps_old_env_and_lock(box: Box) -> None:
    lock = _fgsea(box)
    box.flag("fail_version", "fgsea")
    assert box.run("fgsea").returncode == 1
    prefix = box.prefix("fgsea")
    assert (prefix / "bin" / "payload").read_text() == "old\n"
    assert _leftovers(box, "fgsea") == []
    assert (box.locks / "fgsea.explicit.txt").read_text() == lock


@pytest.mark.parametrize("step", ["clone", "checkout"])
def test_failed_git_step_keeps_old_env(box: Box, step: str) -> None:
    """S3: a failed clone or checkout of a git-pinned tool's source keeps the old env."""
    from gippa.annotate.tools import TOOLS

    box.version("uniop", TOOLS["uniop"].version.removeprefix("git:"))
    box.old_env("uniop")
    (box.fake / f"fail_git_{step}").write_text("")
    assert box.run("uniop").returncode == 1
    _assert_old_kept(box, "uniop")


@pytest.mark.parametrize("which", ["explicit", "pip"])
def test_failed_lock_move_keeps_old_env_and_locks(box: Box, which: str) -> None:
    """S4: a lock that cannot be moved fails the tool before the old env is deleted."""
    box.old_env("bakta")
    _shim(box, "mv", f"*.{which}.txt.tmp", "exit 1")
    proc = box.run("bakta")
    assert proc.returncode == 1, proc.stdout + proc.stderr
    _assert_old_kept(box, "bakta")
    assert not list(box.locks.glob("bakta.*.prev"))


def test_next_run_sweeps_trash_left_by_a_killed_delete(box: Box) -> None:
    """B1: a SIGKILL mid-delete leaves ``.gippa-trash.<pid>``; the next build removes it."""
    trash = box.envs / "gippa-bakta.gippa-trash.12345"
    (trash / "lib").mkdir(parents=True)
    assert box.run("bakta").returncode == 0
    assert _leftovers(box) == []


def test_lock_left_aside_by_a_killed_swap_is_put_back(box: Box) -> None:
    """A kill mid-swap leaves ``<lock>.prev``; a later failed run puts it back, never loses it."""
    box.old_env("bakta")
    _shim(box, "mv", "*.explicit.txt.tmp", _HANG)
    assert _term_at_shim(box, "bakta").returncode == 143
    (box.root / "bin" / "mv").unlink()
    box.flag("fail_version", "bakta")
    assert box.run("bakta").returncode == 1
    assert (box.locks / "bakta.explicit.txt").read_text() == OLD_LOCK
    assert not list(box.locks.glob("bakta.*.prev"))
