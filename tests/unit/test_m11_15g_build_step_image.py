"""M11-15g: the dev build script ``scripts/dev/build_step_image.sh`` and the image names/tags.

Spec docs/specs/M11-15-portable-clusters.md item 8a "Dev builds": ``scripts/dev/build_step_image.sh
[--runtime docker|apptainer|singularity]`` "builds the wheel and both images from the checkout,
tags ``<pkg>-step:dev-<sha12 of the wheel>`` and ``<pkg>-step-r:dev-<sha12>``, and prints the refs
for ``--step-image``"; "Names and tags (amends ADR-0060 decision 1, in the Q-K39 ADR)". The
slice table lists the Q-K39 ADR with 15g. Static, plus one argument-validation run with stub
binaries (nothing is built).

Readings chosen where the spec leaves room:

* The script is bash, executable, passes ``bash -n``, and builds from both Dockerfiles
  (``workflow/containers/step/Dockerfile``, ``workflow/containers/step-r/Dockerfile``).
* An unknown ``--runtime`` is a usage error (exit 2) detected before any side effect: no
  ``uv``/``docker``/``apptainer``/``singularity``/``python`` stub on PATH is called.
* The tag is ``dev-`` plus the first 12 hex of the wheel's sha256 (``sha256sum`` or
  ``shasum -a 256``; ``cut -c1-12``, ``${x:0:12}`` or ``head -c 12``), with image names
  ``<pkg>-step`` and ``<pkg>-step-r`` (``<pkg>`` literal or read from pyproject).
* The ADR: a ``docs/adr/ADR-*.md`` naming ``Q-K39`` and ADR-0060, the ``-step`` and ``-step-r``
  names, and that the tag is the version alone.
"""

from __future__ import annotations

import os
import re
import stat
import subprocess
import tomllib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "dev" / "build_step_image.sh"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
STUBS = ("uv", "docker", "apptainer", "singularity", "python", "python3", "pip")


def _text() -> str:
    assert SCRIPT.is_file(), f"{SCRIPT.relative_to(REPO)} is missing (M11-15g)"
    return SCRIPT.read_text(encoding="utf-8")


def _code() -> str:
    return "\n".join(ln for ln in _text().splitlines() if not ln.lstrip().startswith("#"))


def test_script_exists_and_is_executable_bash() -> None:
    text = _text()
    assert text.startswith("#!") and "bash" in text.splitlines()[0], text.splitlines()[0]
    assert SCRIPT.stat().st_mode & stat.S_IXUSR, "not executable"
    done = subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr


def test_script_is_strict() -> None:
    assert re.search(r"set -[a-z]*e[a-z]*u|set -euo pipefail|set -eu", _code()), "no set -eu"


@pytest.mark.parametrize("runtime", ["docker", "apptainer", "singularity"])
def test_script_knows_each_runtime(runtime: str) -> None:
    assert re.search(rf"\b{runtime}\b", _code()), runtime


def test_script_builds_both_dockerfiles() -> None:
    code = _code()
    assert "workflow/containers/step/Dockerfile" in code
    assert "workflow/containers/step-r/Dockerfile" in code


def test_script_builds_the_wheel() -> None:
    assert re.search(r"\buv\s+build\b|\bpython3?\s+-m\s+build\b", _code()), "no wheel build"


def test_tags_are_dev_sha12_of_the_wheel() -> None:
    code = _code()
    assert re.search(r"sha256sum|shasum\s+-a\s+256|hashlib\.sha256", code), "no sha256 of the wheel"
    assert re.search(r"cut\s+-c\s*1-12|:0:12\}|head\s+-c\s*12|\[:12\]", code), "no 12-hex cut"
    assert re.search(r"-step:dev-", code), "no <pkg>-step:dev-<sha12> tag"
    assert re.search(r"-step-r:dev-", code), "no <pkg>-step-r:dev-<sha12> tag"
    assert PKG in code or "pyproject" in code, f"image name not derived from {PKG}"


def test_unknown_runtime_is_a_usage_error_before_any_side_effect(tmp_path: Path) -> None:
    _text()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    calls = tmp_path / "calls.log"
    for name in STUBS:
        stub = bin_dir / name
        stub.write_text(f'#!/bin/sh\necho "{name} $*" >> "{calls}"\nexit 1\n', encoding="utf-8")
        stub.chmod(0o755)
    env = {
        "PATH": f"{bin_dir}:/usr/bin:/bin",
        "HOME": str(tmp_path),
        "TMPDIR": str(tmp_path),
        "LC_ALL": "C",
    }
    env.update({k: v for k, v in os.environ.items() if k == "SYSTEMROOT"})
    done = subprocess.run(
        ["bash", str(SCRIPT), "--runtime", "podmanx"],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert done.returncode == 2, (done.returncode, done.stdout, done.stderr)
    assert "podmanx" in done.stdout + done.stderr or "runtime" in done.stderr.lower()
    assert not calls.exists(), calls.read_text()


# ---------------------------------------------------------------- the Q-K39 ADR


def test_q_k39_adr_records_names_and_version_tag() -> None:
    adrs = sorted((REPO / "docs" / "adr").glob("ADR-*.md"))
    hits = [p for p in adrs if "Q-K39" in p.read_text(encoding="utf-8")]
    assert hits, "no ADR records Q-K39 (step images; amends ADR-0060 decision 1)"
    text = hits[0].read_text(encoding="utf-8")
    assert "ADR-0060" in text
    assert re.search(r"-step\b", text) and "-step-r" in text, "the step image names are missing"
    assert "<version>" in text or "__version__" in text, "the version tag is missing"
