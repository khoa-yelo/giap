"""Backlog 2026-10-28 (lead) N6: the UniOP preflight probe never lets ``docker run`` pull.

``probe_uniop`` under docker runs the image that ``gippa setup --uniop`` built locally and wrote
into ``<db_dir>/local_tools.config``. No registry holds that tag, so an implicit pull can never
succeed; it only turns a missing local image into a slow, network-dependent failure. The argv
must carry ``--pull=never`` after ``run`` and before the image, so a missing image fails fast
(docker exit 125, "Unable to find image ... locally").

Interface pinned here (implementer, ``src/gippa/preflight/environment.py``):

- ``probe_uniop("docker", image_cache=..., env_root=..., local_tools=<file>)`` runs
  ``docker run --rm --pull=never <tag> sh -c <help>`` (``argv[0]`` the bare name or a path);
  the rest of the argv is unchanged.
- The apptainer and micromamba argvs are unchanged (``--pull`` is a docker flag).

Fakes only: a fake ``docker`` / ``apptainer`` first on ``PATH`` that logs its argv and behaves like
docker for a missing image (with ``--pull=never``: exit 125 at once; without: a refused "pull").
No network, no real ``$GIPPA_DB_DIR``.
"""

from __future__ import annotations

import json
import os
import stat
import sys
from pathlib import Path
from typing import Any

import pytest

TAG = "localhost/gippa-uniop:0123456789ab"

FAKE = r'''#!{python} -I
"""Fake docker/apptainer (tests only): log argv; docker reports a missing image like docker."""
import json, os, sys
from pathlib import Path
exe = Path(sys.argv[0]).name
argv = sys.argv[1:]
with open(os.environ["FAKE_RT_LOG"], "a") as fh:
    fh.write(json.dumps({{"exe": exe, "argv": argv}}) + "\n")
if exe == "docker" and os.environ.get("FAKE_IMAGE_MISSING") == "1":
    image = os.environ["FAKE_IMAGE"]
    if "--pull=never" in argv:
        sys.stderr.write("Unable to find image '" + image + "' locally\n"
                         "docker: Error response from daemon: No such image: " + image + ".\n")
        sys.exit(125)
    sys.stderr.write("fake docker: network pull of " + image + " attempted\n")
    sys.exit(97)
sys.stderr.write("usage: UniOP [-h]\n")
sys.exit(0)
'''


def _env_mod() -> Any:
    from gippa.preflight import environment

    return environment


@pytest.fixture
def fake_rt(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    for exe in ("docker", "apptainer"):
        p = bindir / exe
        p.write_text(FAKE.format(python=sys.executable))
        p.chmod(p.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    log = tmp_path / "runtime.log"
    monkeypatch.setenv("FAKE_RT_LOG", str(log))
    monkeypatch.setenv("FAKE_IMAGE", TAG)
    monkeypatch.delenv("FAKE_IMAGE_MISSING", raising=False)
    return log


def _calls(log: Path) -> list[dict[str, Any]]:
    return [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []


def _local_tools(tmp_path: Path, runtime: str, container: str) -> Path:
    """``local_tools.config`` as ``gippa setup --uniop`` writes it (one ``runtime`` block)."""
    from gippa.setupdb.uniop_build import local_tools_text

    path = tmp_path / "db" / "local_tools.config"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(local_tools_text("", runtime, container), encoding="utf-8")
    return path


def _probe(profile: str, tmp_path: Path, local_tools: Path) -> tuple[int, str]:
    code, err = _env_mod().probe_uniop(
        profile, image_cache=tmp_path / "cache", env_root=tmp_path / "envs", local_tools=local_tools
    )
    return int(code), str(err)


@pytest.mark.real_uniop_probe
def test_docker_probe_passes_pull_never_before_the_image(tmp_path: Path, fake_rt: Path) -> None:
    code, err = _probe("docker", tmp_path, _local_tools(tmp_path, "docker", TAG))
    assert code == 0, err
    (call,) = [c for c in _calls(fake_rt) if c["exe"] == "docker"]
    argv = call["argv"]
    assert argv[0] == "run", argv
    assert TAG in argv, argv
    assert "--pull=never" in argv[: argv.index(TAG)], f"no --pull=never before the image: {argv}"


@pytest.mark.real_uniop_probe
def test_docker_probe_argv_is_otherwise_unchanged(tmp_path: Path, fake_rt: Path) -> None:
    """Only ``--pull=never`` is added: ``run --rm <tag> sh -c <help>`` around it."""
    _probe("docker", tmp_path, _local_tools(tmp_path, "docker", TAG))
    (call,) = [c for c in _calls(fake_rt) if c["exe"] == "docker"]
    argv = call["argv"]
    assert argv.count("--pull=never") == 1, argv
    rest = [a for a in argv if a != "--pull=never"]
    assert rest == ["run", "--rm", TAG, "sh", "-c", _env_mod()._UNIOP_HELP], argv


@pytest.mark.real_uniop_probe
def test_missing_local_image_fails_fast_without_a_pull(
    tmp_path: Path, fake_rt: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A missing local UniOP image gives docker's own exit 125, never a registry pull."""
    monkeypatch.setenv("FAKE_IMAGE_MISSING", "1")
    code, err = _probe("docker", tmp_path, _local_tools(tmp_path, "docker", TAG))
    assert "network pull" not in err, err
    assert code == 125, (code, err)
    assert "Unable to find image" in err, err


@pytest.mark.real_uniop_probe
def test_apptainer_probe_has_no_pull_flag(tmp_path: Path, fake_rt: Path) -> None:
    """Guard: ``--pull`` is docker-only; the apptainer argv is unchanged."""
    sif = tmp_path / "uniop.sif"
    sif.write_bytes(b"fake sif\n")
    code, err = _probe("apptainer", tmp_path, _local_tools(tmp_path, "apptainer", str(sif)))
    assert code == 0, err
    (call,) = [c for c in _calls(fake_rt) if c["exe"] == "apptainer"]
    assert not any(a.startswith("--pull") for a in call["argv"]), call["argv"]
    assert call["argv"][:3] == ["exec", "--cleanenv", str(sif)], call["argv"]
