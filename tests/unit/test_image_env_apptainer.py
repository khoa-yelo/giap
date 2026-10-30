"""The published images must work under `apptainer exec`, which skips the image's entrypoint.

Found by the M11-14 wheel smoke: no `python` on PATH and no CONDA_PREFIX in the conda-only images,
TERMINATORNET_PLAN with no container, and no libGL for Promotech's PyQt5 import.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[2]
CONDA_ONLY = ["promotech", "promoter_calculator", "terminatornet", "bactermfinder", "ostir"]
REF = "ghcr.io/khoa-yelo/gippa-terminatornet@sha256:" + "a" * 64


def _update_digests() -> ModuleType:
    path = REPO / "scripts" / "release" / "update_digests.py"
    spec = importlib.util.spec_from_file_location("update_digests_env", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("tool", CONDA_ONLY)
def test_image_sets_conda_env_without_the_entrypoint(tool: str) -> None:
    text = (REPO / "workflow" / "containers" / tool / "Dockerfile").read_text()
    env = "\n".join(line for line in text.splitlines() if line.startswith("ENV "))
    assert '"/opt/conda/bin:$PATH"' in env
    assert "CONDA_PREFIX=/opt/conda" in env


def test_promotech_image_has_libgl_and_stays_unprivileged() -> None:
    text = (REPO / "workflow" / "containers" / "promotech" / "Dockerfile").read_text()
    assert "libgl1" in text
    users = re.findall(r"^USER (.+)$", text, flags=re.M)
    assert users[-1] == "$MAMBA_USER"


def test_every_process_on_a_published_tool_env_is_pinnable() -> None:
    """A process on a published tool's env needs that tool's image under a container profile."""
    tools = _update_digests().PROCESS_TOOL
    envs = "|".join(CONDA_ONLY)
    found = []
    for module in (REPO / "workflow" / "modules").rglob("main.nf"):
        text = module.read_text()
        env = re.search(rf'envs/({envs})\.yml"', text)
        process = re.search(r"^process (\w+)", text, flags=re.M)
        if env and process:
            found.append(process.group(1))
            assert tools.get(process.group(1)) == env.group(1), process.group(1)
    assert "TERMINATORNET_CHUNK" in found


def test_a_terminatornet_digest_pins_the_sharded_processes() -> None:
    pins = _update_digests().validate({"TERMINATORNET": REF})
    assert pins == dict.fromkeys(
        ["TERMINATORNET", "TERMINATORNET_PLAN", "TERMINATORNET_CHUNK", "TERMINATORNET_MERGE"], REF
    )
