"""#465 review S1 (backlog N6): under docker, preflight's fgsea row checks the image is local.

Runs never pull (``docker run --pull=never``, N6), so the M11-14c rule "under docker the
``images.json`` entry suffices (``docker run`` pulls it)" no longer holds: a pinned but absent
image would pass preflight and then fail inside the function step. Interface pinned here
(implementer, ``gippa.preflight.environment.missing_tool_reasons``):

- under ``docker`` without a built ``gippa-fgsea`` env, fgsea is present iff ``docker image
  inspect <images.json["fgsea"]["ref"]>`` succeeds (the existing ``_docker_image_exists``); when
  it fails (or docker is not on ``PATH``) the row has a reason naming the ref and ``<pkg> setup
  --images --profile docker``;
- a built env still satisfies fgsea under docker without any docker call;
- plain ``apptainer`` is unchanged: the ``.img`` in the cache counts, no docker call, and a
  missing ``.img`` names ``<pkg> setup --images --profile apptainer``.

Fakes only: a fake workflow tree, temporary cache/env directories and a logging fake ``docker``
first on ``PATH`` whose ``image inspect`` succeeds only for the refs in ``FAKE_PRESENT``. No
network, no real docker, never the real ``$GIPPA_ROOT``/``$GIPPA_DB_DIR``.
"""

from __future__ import annotations

import hashlib
import json
import os
import stat
import sys
import tomllib
from pathlib import Path
from typing import Any

import pytest

from gippa.annotate.tools import TOOLS

REPO = Path(__file__).resolve().parents[3]
WORKFLOW = REPO / "workflow"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
REF = f"ghcr.io/khoa-yelo/{PKG}-fgsea@sha256:" + "c3" * 32
SETUP_DOCKER = f"{PKG} setup --images --profile docker"
SETUP_APPTAINER = f"{PKG} setup --images --profile apptainer"

FAKE_DOCKER = r'''#!{python} -I
"""Fake docker (tests only): ``image inspect`` succeeds only for refs in FAKE_PRESENT."""
import json, os, sys
argv = sys.argv[1:]
with open(os.environ["FAKE_RT_LOG"], "a") as fh:
    fh.write(json.dumps(argv) + "\n")
present = [r for r in os.environ.get("FAKE_PRESENT", "").split(",") if r]
if argv[:2] == ["image", "inspect"]:
    refs = argv[2:]
    if refs and all(r in present for r in refs):
        sys.stdout.write(json.dumps([{{"Id": "sha256:" + "ab" * 32}}]) + "\n")
        sys.exit(0)
    sys.stdout.write("[]\n")
    sys.stderr.write("Error: No such image: " + refs[-1] + "\n")
    sys.exit(1)
sys.stderr.write("fake docker: unexpected call\n")
sys.exit(3)
'''


def _simple_name(ref: str) -> str:
    return ref.removeprefix("docker://").replace("/", "-").replace(":", "-") + ".img"


def _env_mod() -> Any:
    from gippa.preflight import environment

    return environment


@pytest.fixture(autouse=True)
def log(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    for var in (
        "GIPPA_DB_DIR",
        "GIPPA_ROOT",
        "NXF_APPTAINER_CACHEDIR",
        "GIPPA_FGSEA_PREFIX",
        "GIPPA_CONTAINER_RUNTIME",
        "FAKE_PRESENT",
    ):
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "docker"
    exe.write_text(FAKE_DOCKER.format(python=sys.executable))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    path = tmp_path / "docker.log"
    monkeypatch.setenv("FAKE_RT_LOG", str(path))
    _fake_tree(tmp_path, monkeypatch)
    return path


def _fake_tree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    wf = tmp_path / "wf"
    (wf / "conf").mkdir(parents=True)
    (wf / "main.nf").write_text("// fake\n")
    (wf / "envs").symlink_to(WORKFLOW / "envs")
    (wf / "bin").symlink_to(WORKFLOW / "bin")
    (wf / "conf" / "containers.config").write_text(
        (WORKFLOW / "conf" / "containers.config").read_text(encoding="utf-8")
    )
    images = json.loads((WORKFLOW / "conf" / "images.json").read_text())
    images["fgsea"] = {"key": "0123456789ab", "tag": "1.36.2-0123456789ab", "ref": REF,
                       "processes": [], "built_from": "0" * 40, "licence": "MIT"}  # fmt: skip
    (wf / "conf" / "images.json").write_text(json.dumps(images, indent=2) + "\n")
    monkeypatch.setattr(_env_mod(), "workflow_main", lambda: wf / "main.nf")
    return wf


def _cache(tmp_path: Path, img: bool) -> Path:
    c = tmp_path / "cache"
    c.mkdir(exist_ok=True)
    if img:
        (c / _simple_name(REF)).write_bytes(b"fake sif\n")
    return c


def _envs(tmp_path: Path, built: bool) -> Path:
    root = tmp_path / "envs"
    root.mkdir(exist_ok=True)
    if built:
        prefix = root / "gippa-fgsea"
        (prefix / "bin").mkdir(parents=True)
        (prefix / "bin" / ("R" + "script")).write_text("#!/bin/sh\n")
        sha = hashlib.sha256((WORKFLOW / "envs" / "fgsea.yml").read_bytes()).hexdigest()
        (prefix / ".gippa-built").write_text(sha + "\n")
    return root


def _reasons(tmp_path: Path, profile: str, *, img: bool = False, built: bool = False) -> Any:
    return _env_mod().missing_tool_reasons(
        [TOOLS["fgsea"]], profile, _cache(tmp_path, img), _envs(tmp_path, built), None
    )


def _inspected(log: Path) -> list[list[str]]:
    calls = [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []
    return [c for c in calls if c[:2] == ["image", "inspect"]]


# ---------------------------------------------------------------- docker


def test_docker_missing_local_image_fails_the_fgsea_row(tmp_path: Path, log: Path) -> None:
    reasons = _reasons(tmp_path, "docker")
    assert set(reasons) == {"fgsea"}, (
        "under docker a pinned but absent fgsea image must fail preflight: runs never pull (N6)"
    )
    assert SETUP_DOCKER in reasons["fgsea"], reasons
    assert REF in reasons["fgsea"], reasons
    assert any(REF in c for c in _inspected(log)), "the ref was never checked with image inspect"


def test_docker_missing_image_is_reported_by_missing_tools(tmp_path: Path) -> None:
    fn = _env_mod().missing_tools
    assert fn(
        [TOOLS["fgsea"]], "docker", _cache(tmp_path, False), _envs(tmp_path, False), None
    ) == ["fgsea"]


def test_docker_present_local_image_passes(
    tmp_path: Path, log: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_PRESENT", REF)
    assert _reasons(tmp_path, "docker") == {}
    assert any(REF in c for c in _inspected(log)), "the ref was never checked with image inspect"


def test_docker_without_docker_on_path_fails_the_fgsea_row(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "nobin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    reasons = _reasons(tmp_path, "docker")
    assert set(reasons) == {"fgsea"}, reasons
    assert SETUP_DOCKER in reasons["fgsea"], reasons


def test_docker_built_env_still_satisfies_fgsea(tmp_path: Path) -> None:
    """The env comes first (prefix -> env -> image): no image needed."""
    assert _reasons(tmp_path, "docker", built=True) == {}


# ---------------------------------------------------------------- apptainer unchanged


def test_apptainer_cached_img_passes_without_docker(tmp_path: Path, log: Path) -> None:
    assert _reasons(tmp_path, "apptainer", img=True) == {}
    assert _inspected(log) == [], "apptainer must not consult docker"


def test_apptainer_missing_img_names_setup_apptainer(tmp_path: Path, log: Path) -> None:
    reasons = _reasons(tmp_path, "apptainer", img=False)
    assert set(reasons) == {"fgsea"}, reasons
    assert SETUP_APPTAINER in reasons["fgsea"], reasons
    assert SETUP_DOCKER not in reasons["fgsea"], reasons
    assert _inspected(log) == [], "apptainer must not consult docker"
