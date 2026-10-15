"""farmshare,apptainer: conda only for the conda-only tools (M2-07 decision 4, PR #65 review).

The five tools without an image (uniop, promotech, promoter_calculator, terminatornet, ostir) run
in the prebuilt micromamba envs under ``$GIAP_ROOT`` (``conda.enabled``, ``conda.useMicromamba``);
the pinned-container tools (bakta, transtermhp, genomad, antiSMASH) run in their containers with
conda off. Parses ``nextflow config workflow -flat -profile farmshare,apptainer`` with a fake
``GIAP_ROOT``; the per-process ``conda`` is the last matching ``withName`` setting (the selector
is a regex on the process name, ``ANNOTATE:<TOOL>`` in the pipeline), else the module's env file.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
FAKE_ROOT = "/fake/giap_root"

CONDA_TOOLS = ["UNIOP", "PROMOTECH", "PROMOTER_CALCULATOR", "TERMINATORNET", "OSTIR"]
CONTAINER_TOOLS = ["BAKTA", "TRANSTERMHP", "GENOMAD", "ANTISMASH"]
_OFF = {"null", "''", '""', "false"}


def _flat(profile: str) -> list[str]:
    env = {**os.environ, "NXF_ANSI_LOG": "false", "GIAP_ROOT": FAKE_ROOT}
    env.pop("NXF_APPTAINER_CACHEDIR", None)
    cmd = ["nextflow", "config", str(WF), "-flat", "-profile", profile]
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)
    assert res.returncode == 0, f"{' '.join(cmd)} failed:\n{res.stdout}\n{res.stderr}"
    return res.stdout.splitlines()


def _value(lines: list[str], key: str) -> str | None:
    for line in lines:
        if line.startswith(f"{key} = "):
            return line.split(" = ", 1)[1].strip()
    return None


def _matches(selector: str, name: str) -> bool:
    neg = selector.startswith("!")
    pat = selector[1:] if neg else selector
    hit = any(re.fullmatch(pat, n) for n in (name, name.split(":")[-1]))
    return hit != neg


def _setting(lines: list[str], name: str, directive: str) -> str | None:
    """Last ``withName`` value of ``directive`` that applies to process ``name`` (None: unset)."""
    got = None
    for line in lines:
        m = re.match(r"^process\.'withName:(.+)'\.(\w+) = (.*)$", line)
        if m and m.group(2) == directive and _matches(m.group(1), name):
            got = m.group(3).strip()
    return got


@pytest.fixture(scope="module")
def flat() -> list[str]:
    return _flat("farmshare,apptainer")


def test_conda_is_enabled_with_micromamba(flat: list[str]) -> None:
    assert _value(flat, "conda.enabled") == "true"
    assert _value(flat, "conda.useMicromamba") == "true"
    assert _value(flat, "apptainer.enabled") == "true"


@pytest.mark.parametrize("tool", CONDA_TOOLS)
def test_conda_only_tool_uses_the_prebuilt_env_under_giap_root(flat: list[str], tool: str) -> None:
    value = _setting(flat, f"ANNOTATE:{tool}", "conda")
    assert value is not None, f"{tool}: no withName conda setting (would build the yml env)"
    path = value.strip("'\"")
    assert path.startswith(FAKE_ROOT + "/"), f"{tool}: conda {path} is not under $GIAP_ROOT"
    assert tool.lower() in Path(path).name, f"{tool}: env {path} is not the tool's env"
    container = _setting(flat, f"ANNOTATE:{tool}", "container")
    assert container is None or container in _OFF, f"{tool} is conda-only but has {container}"


@pytest.mark.parametrize("tool", CONTAINER_TOOLS)
def test_container_tool_has_conda_off(flat: list[str], tool: str) -> None:
    container = _setting(flat, f"ANNOTATE:{tool}", "container")
    assert container and container.strip("'\"").startswith("docker://"), f"{tool}: no container pin"
    value = _setting(flat, f"ANNOTATE:{tool}", "conda")
    assert value in _OFF, f"{tool}: conda {value!r} with conda.enabled (must be off for this tool)"


def test_exactly_the_five_tools_get_conda(flat: list[str]) -> None:
    on = sorted(
        t
        for t in CONDA_TOOLS + CONTAINER_TOOLS
        if (_setting(flat, f"ANNOTATE:{t}", "conda") or "module-env") not in _OFF
    )
    assert on == sorted(CONDA_TOOLS)


def test_stub_profile_keeps_conda_off() -> None:
    assert _value(_flat("test,stub"), "conda.enabled") in {"false", None}
