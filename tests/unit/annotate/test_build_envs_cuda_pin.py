"""build_envs.sh takes the CPU-node CONDA_OVERRIDE_CUDA from the yml's ``cuda-version`` pin.

Backlog 2026-10-28 (#457/#458 review N2): the override default was a hard-coded 11.8. Order now:
an exported ``CONDA_OVERRIDE_CUDA``, then the yml's exact ``cuda-version`` pin, then 11.8.
``micromamba`` and ``nvidia-smi`` are shell stubs in ``$GIPPA_ROOT/bin`` (first on PATH via
``farmshare_env.sh``); the stub ``create`` logs the CONDA_OVERRIDE_CUDA it sees and never solves.
The harness is a small copy of the one in ``test_build_envs_swap.py``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]

FAKE_MICROMAMBA = r"""#!/bin/bash
echo "micromamba $* CUDA=${CONDA_OVERRIDE_CUDA-unset}" >> "$FAKE_DIR/calls.log"
envs="$MAMBA_ROOT_PREFIX/envs"
name_of() { while [ "$#" -gt 0 ]; do if [ "$1" = "-n" ]; then echo "$2"; return; fi; shift; done; }
case "$1" in
  env)
    [ "$2" = "export" ] && { echo "@EXPLICIT"; echo "https://example.invalid/x-1.0-0.conda"; }
    exit 0 ;;
  create)
    n="$(name_of "$@")"
    mkdir -p "$envs/$n/conda-meta" "$envs/$n/bin"
    exit 0 ;;
  run)
    shift
    n=""
    while [ "$#" -gt 0 ]; do
      case "$1" in -n) n="$2"; shift 2 ;; -*) shift ;; *) break ;; esac
    done
    case " $* " in *"pip freeze"*) echo "pkg==1.0"; exit 0 ;; esac
    cat "$FAKE_DIR/version_$n" ;;
esac
exit 0
"""

FAKE_NVIDIA_SMI = r"""#!/bin/bash
[ -e "$FAKE_DIR/gpu" ] || { echo "NVIDIA-SMI has failed" >&2; exit 9; }
echo "GPU 0: NVIDIA L40S (UUID: GPU-00000000-0000-0000-0000-000000000000)"
"""

CUDA_BUILD = "  - pytorch=2.0.*=cuda*\n"


class Box:
    def __init__(self, tmp: Path) -> None:
        self.repo, self.root, self.fake = tmp / "repo", tmp / "root", tmp / "fake"
        for d in (
            self.repo / "workflow" / "envs" / "locks",
            self.repo / "workflow" / "conf",
            self.fake,
            self.root / "bin",
        ):
            d.mkdir(parents=True)
        (self.repo / "src").symlink_to(REPO / "src")
        (self.repo / "scripts").symlink_to(REPO / "scripts")
        for yml in (REPO / "workflow" / "envs").glob("*.yml"):
            shutil.copy(yml, self.repo / "workflow" / "envs" / yml.name)
        shutil.copy(
            REPO / "workflow" / "conf" / "containers.config",
            self.repo / "workflow" / "conf" / "containers.config",
        )
        for name, body in (("micromamba", FAKE_MICROMAMBA), ("nvidia-smi", FAKE_NVIDIA_SMI)):
            path = self.root / "bin" / name
            path.write_text(body)
            path.chmod(0o755)
        (self.fake / "calls.log").write_text("")

    def yml(self, tool: str) -> Path:
        return self.repo / "workflow" / "envs" / f"{tool}.yml"

    def append(self, tool: str, text: str) -> None:
        self.yml(tool).write_text(self.yml(tool).read_text() + text)

    def run(self, tool: str, **extra: str) -> str:
        """Build ``tool``; return the CONDA_OVERRIDE_CUDA its create saw (``unset`` for none)."""
        from gippa.annotate.tools import PILOT_TOOLS, TOOLS

        t = {**PILOT_TOOLS, **TOOLS}[tool]
        (self.fake / f"version_{t.env_name}").write_text(t.version.removeprefix("git:") + "\n")
        env = {k: v for k, v in os.environ.items() if k in {"PATH", "LANG", "TMPDIR"}}
        env.update(
            GIPPA_REPO=str(self.repo),
            GIPPA_ROOT=str(self.root),
            GIPPA_PYTHON=sys.executable,
            FAKE_DIR=str(self.fake),
            HOME=str(self.root / "home"),
            **extra,
        )
        proc = subprocess.run(
            ["bash", str(self.repo / "scripts" / "build_envs.sh"), tool],
            capture_output=True,
            text=True,
            cwd=self.repo,
            env=env,
            timeout=300,
            check=False,
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
        lines = (self.fake / "calls.log").read_text().splitlines()
        (create,) = [c for c in lines if c.startswith("micromamba create")]
        return create.rsplit("CUDA=", 1)[1]


@pytest.fixture
def box(tmp_path: Path) -> Box:
    return Box(tmp_path)


def test_cascadia_pins_cuda_11_8_and_gets_11_8(box: Box) -> None:
    assert "  - cuda-version=11.8\n" in box.yml("cascadia").read_text()
    assert box.run("cascadia") == "11.8"


def test_cascadia_override_follows_its_pin(box: Box) -> None:
    yml = box.yml("cascadia")
    yml.write_text(yml.read_text().replace("cuda-version=11.8", "cuda-version=12.1"))
    assert box.run("cascadia") == "12.1"


@pytest.mark.parametrize(
    ("pin", "want"),
    [
        ("  - cuda-version=12.1\n", "12.1"),
        ("  - cuda-version==12.4\n", "12.4"),
        ("  - cuda-version 11.8.*\n", "11.8"),
        ("  - cuda-version=12.1.*\n", "12.1"),
        ("  - cuda-version=12.2=h1234_0\n", "12.2"),
        ('  - "cuda-version=12.3"\n', "12.3"),
        ("  - cuda-version=12.1  # L40S\n", "12.1"),
        ("  - cuda-version>=11.8,<12\n", "11.8"),  # a range is not a pin: the default
        ("  # - cuda-version=12.1\n", "11.8"),  # commented out
        ("  - cuda-versions=12.1\n", "11.8"),  # another package
        ("", "11.8"),  # no pin
    ],
)
def test_override_comes_from_the_pin(box: Box, pin: str, want: str) -> None:
    box.append("bakta", pin + CUDA_BUILD)
    assert box.run("bakta") == want


def test_exported_value_wins_over_the_pin(box: Box) -> None:
    box.append("bakta", "  - cuda-version=12.1\n" + CUDA_BUILD)
    assert box.run("bakta", CONDA_OVERRIDE_CUDA="12.4") == "12.4"


def test_pin_without_cuda_build_sets_no_override(box: Box) -> None:
    box.append("bakta", "  - cuda-version=12.1\n")
    assert box.run("bakta") == "unset"


def test_pin_on_gpu_node_sets_no_override(box: Box) -> None:
    box.append("bakta", "  - cuda-version=12.1\n" + CUDA_BUILD)
    (box.fake / "gpu").write_text("")
    assert box.run("bakta") == "unset"
