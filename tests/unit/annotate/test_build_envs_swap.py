"""build_envs.sh keeps the old env until the new one verifies, and CUDA solves on CPU nodes.

Backlog (build_envs): job 1888853 lost the cascadia env and its locks because the old env was
removed before the new one solved; job 1901404 needed ``CONDA_OVERRIDE_CUDA`` by hand. No conda is
needed: ``micromamba``, ``git`` and ``nvidia-smi`` are shell stubs in ``$GIPPA_ROOT/bin`` (first on
PATH via ``farmshare_env.sh``). The stub ``create`` writes a half-built prefix before it fails, so a
restored env is told apart from a new one by its contents. Nothing is ever solved.
"""

from __future__ import annotations

import os
import shutil
import signal
import subprocess
import sys
import time
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


@pytest.fixture
def box(tmp_path: Path) -> Box:
    b = Box(tmp_path)
    b.version("bakta", "bakta 1.12.2")
    return b


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
