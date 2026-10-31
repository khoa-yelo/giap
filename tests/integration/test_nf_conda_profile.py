"""farmshare,apptainer: conda only for the conda-only tools (M2-07 decisions 4, 12 and 13).

The one process tool without a published image (uniop) runs in its prebuilt micromamba env under
``$GIPPA_ROOT`` (``conda.enabled``, ``conda.useMicromamba``); the pinned-container tools (bakta,
transtermhp, genomad, antiSMASH and the GHCR images of promotech, promoter_calculator,
terminatornet, ostir and the opt-in bactermfinder) run in their containers with conda off; no
conda env is ever built on the fly
under farmshare or farmshare,apptainer. The ``ps`` shim bound into the bakta/antiSMASH containers
is read-only (decision 13). Parses ``nextflow config workflow -flat -profile <profile>`` with a
fake ``GIPPA_ROOT``; the per-process ``conda`` is the last matching ``withName`` setting (the
selector is a regex on the process name, ``ANNOTATE:<TOOL>`` in the pipeline), else the module's
env file.
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
FAKE_ROOT = "/fake/gippa_root"

CONDA_TOOLS = ["UNIOP"]
PROFILES = ["farmshare", "farmshare,apptainer"]
CONTAINER_TOOLS = [
    "BAKTA",
    "TRANSTERMHP",
    "GENOMAD",
    "ANTISMASH",
    "PROMOTECH",
    "PROMOTER_CALCULATOR",
    "TERMINATORNET",
    "OSTIR",
    "BACTERMFINDER",
]
_OFF = {"null", "''", '""', "false"}


def _flat(profile: str) -> list[str]:
    env = {**os.environ, "NXF_ANSI_LOG": "false", "GIPPA_ROOT": FAKE_ROOT}
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


_FLAT_CACHE: dict[str, list[str]] = {}


def _flat_cached(profile: str) -> list[str]:
    if profile not in _FLAT_CACHE:
        _FLAT_CACHE[profile] = _flat(profile)
    return _FLAT_CACHE[profile]


@pytest.fixture(scope="module")
def flat() -> list[str]:
    return _flat_cached("farmshare,apptainer")


def test_conda_is_enabled_with_micromamba(flat: list[str]) -> None:
    assert _value(flat, "conda.enabled") == "true"
    assert _value(flat, "conda.useMicromamba") == "true"
    assert _value(flat, "apptainer.enabled") == "true"


@pytest.mark.parametrize("profile", PROFILES)
@pytest.mark.parametrize("tool", CONDA_TOOLS)
def test_conda_only_tool_uses_the_prebuilt_env_under_gippa_root(profile: str, tool: str) -> None:
    """Decision 12: a conda-only tool uses a prebuilt env."""
    flat = _flat_cached(profile)
    value = _setting(flat, f"ANNOTATE:{tool}", "conda")
    assert value is not None, f"{tool}: no withName conda setting (would build the yml env)"
    path = value.strip("'\"")
    assert path.startswith(FAKE_ROOT + "/"), f"{tool}: conda {path} is not under $GIPPA_ROOT"
    assert tool.lower() in Path(path).name, f"{tool}: env {path} is not the tool's env"
    container = _setting(flat, f"ANNOTATE:{tool}", "container")
    assert container is None or container in _OFF, f"{tool} is conda-only but has {container}"


@pytest.mark.parametrize("tool", CONTAINER_TOOLS)
def test_container_tool_has_conda_off(flat: list[str], tool: str) -> None:
    container = _setting(flat, f"ANNOTATE:{tool}", "container")
    assert container and re.fullmatch(
        r"[a-z0-9][a-z0-9._/-]*@sha256:[0-9a-f]{64}", container.strip("'\"")
    ), f"{tool}: no container pin"
    value = _setting(flat, f"ANNOTATE:{tool}", "conda")
    assert value in _OFF, f"{tool}: conda {value!r} with conda.enabled (must be off for this tool)"


def test_exactly_the_conda_only_tools_get_conda(flat: list[str]) -> None:
    on = sorted(
        t
        for t in CONDA_TOOLS + CONTAINER_TOOLS
        if (_setting(flat, f"ANNOTATE:{t}", "conda") or "module-env") not in _OFF
    )
    assert on == sorted(CONDA_TOOLS)


def test_stub_profile_keeps_conda_off() -> None:
    assert _value(_flat("test,stub"), "conda.enabled") in {"false", None}


@pytest.mark.parametrize("tool", ["BAKTA", "ANTISMASH"])
def test_ps_shim_bind_is_read_only(flat: list[str], tool: str) -> None:
    """Decision 13: the /proc-based ps shim is bound read-only over /usr/bin/ps."""
    opts = _setting(flat, f"ANNOTATE:{tool}", "containerOptions")
    assert opts, f"{tool}: no containerOptions (ps shim not bound)"
    m = re.search(r"-B\s+(\S*assets/ps/ps):/usr/bin/ps(:\S+)?", opts.strip("'\""))
    assert m, f"{tool}: no ps shim bind in {opts}"
    assert m.group(2) == ":ro", f"{tool}: ps shim bind is not read-only: {m.group(0)}"
