"""Backlog 2026-10-28 (lead) N6: the R fgsea docker launcher never pulls at run time.

Lead decision (2026-10-28): ``gippa setup`` is the only step that pulls images (R6, offline
runs). The docker argv of ``gippa.io.fgsea_r.command`` therefore carries ``--pull=never`` after
``run`` and before the pinned ``images.json["fgsea"]["ref"]``; the rest of the M11-14c argv is
unchanged. When the image is not present locally the call fails at once, and the
``FgseaEngineError`` it raises tells the user to run ``gippa setup``.

Interface pinned here (implementer, ``src/gippa/io/fgsea_r.py``):

- docker: ``docker run --rm --pull=never --cidfile <d>.cid -u <uid>:<gid> -v <d>:<d> -w <d>
  -e LC_ALL=C <ref> Rscript --vanilla /opt/<pkg>/gippa_fgsea.R <options>`` (where exactly
  ``--pull=never`` sits among the options before ``<ref>`` is free);
- apptainer: unchanged, no ``--pull``;
- a missing local image (docker: exit 125, "Unable to find image '<ref>' locally"), whether
  detected by the failing ``docker run`` or by a ``docker image inspect`` beforehand, raises
  ``FgseaEngineError`` (or a subclass) from ``resolve_launcher`` + ``RFgseaEngine`` whose message
  names ``gippa setup``; nothing tries a registry.

Fakes only: a fake workflow tree (``conf/images.json``) and a fake ``docker`` first on ``PATH``
that logs its argv and acts as if the image were absent. No network, no R.
"""

from __future__ import annotations

import json
import os
import stat
import sys
import tomllib
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[3]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
REF = f"ghcr.io/khoa-yelo/{PKG}-fgsea@sha256:" + "c3" * 32
BAKED = f"/opt/{PKG}/gippa_fgsea.R"
UIDGID = f"{os.getuid()}:{os.getgid()}"

FAKE_DOCKER = r'''#!{python} -I
"""Fake docker (tests only): the image is absent. Logs argv; never reaches a registry."""
import json, os, sys
argv = sys.argv[1:]
with open(os.environ["FAKE_RT_LOG"], "a") as fh:
    fh.write(json.dumps(argv) + "\n")
image = os.environ["FAKE_IMAGE"]
if argv[:2] == ["image", "inspect"]:
    sys.stderr.write("Error: No such image: " + image + "\n")
    sys.exit(1)
if argv[:1] in (["kill"], ["rm"]):
    sys.exit(0)
if argv[:1] == ["run"] and "--pull=never" in argv:
    sys.stderr.write("Unable to find image '" + image + "' locally\n"
                     "docker: Error response from daemon: No such image: " + image + ".\n")
    sys.exit(125)
if argv[:1] == ["run"]:
    sys.stderr.write("fake docker: network pull of " + image + " attempted\n")
    sys.exit(97)
sys.stderr.write("fake docker: unexpected call\n")
sys.exit(3)
'''


def _fr() -> Any:
    from gippa.io import fgsea_r

    return fgsea_r


def _call(seed: int = 424242) -> Any:
    from gippa.config.model import GseaConfig
    from gippa.function.gsea import EngineCall, GseaJob

    ids = tuple(f"TU_{i}" for i in (5, 2, 7, 1, 4, 6, 3, 8))
    return EngineCall(
        tu_ids=ids,
        scores=np.array([3.0, 2.5, 1.0, 0.5, -0.25, -1.0, -2.0, -3.5], dtype=np.float64),
        jobs=(
            GseaJob("main:cog", {"C1": ("TU_1", "TU_2", "TU_3"), "C2": ("TU_4", "TU_5")}, 2, 7),
            GseaJob("main:kegg", {"K1": ("TU_6", "TU_7", "TU_8")}, 2, 7),
        ),
        params=GseaConfig(),
        seed=seed,
    )


def _options(call: Any) -> list[str]:
    p = call.params
    return [
        "--stats", "stats.tsv",
        "--sets", "sets.tsv",
        "--jobs", "jobs.tsv",
        "--sample-size", str(int(p.sample_size)),
        "--eps", repr(float(p.eps)),
        "--n-perm-simple", str(int(p.n_perm_simple)),
        "--seed", str(int(call.seed)),
        "--out", "results.tsv",
        "--versions", "versions.txt",
    ]  # fmt: skip


def _old_docker_argv(d: Path, call: Any, ref: str = REF) -> list[str]:
    """The M11-14c docker argv without ``--pull=never``."""
    return ["docker", "run", "--rm", "--cidfile", f"{d}.cid", "-u", UIDGID,
            "-v", f"{d}:{d}", "-w", str(d), "-e", "LC_ALL=C", ref,
            "Rscript", "--vanilla", BAKED, *_options(call)]  # fmt: skip


def _norm(argv: list[str]) -> list[str]:
    return [Path(argv[0]).name, *argv[1:]]


def _assert_pull_never(argv: list[str], ref: str) -> None:
    assert argv[1] == "run", argv
    assert argv.count(ref) == 1, argv
    assert argv.count("--pull=never") == 1, f"want exactly one --pull=never: {argv}"
    assert argv.index("--pull=never") < argv.index(ref), f"--pull=never after the image: {argv}"


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    for var in (
        "GIPPA_DB_DIR",
        "GIPPA_ROOT",
        "NXF_APPTAINER_CACHEDIR",
        "GIPPA_FGSEA_PREFIX",
        "GIPPA_CONTAINER_RUNTIME",
        "GIPPA_FGSEA_IMAGE",
    ):
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "docker"
    exe.write_text(FAKE_DOCKER.format(python=sys.executable))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    log = tmp_path / "docker.log"
    monkeypatch.setenv("FAKE_RT_LOG", str(log))
    monkeypatch.setenv("FAKE_IMAGE", REF)
    return log


def _wf(root: Path, ref: str = REF) -> Path:
    """A fake ``workflow/`` tree: ``main.nf`` and ``conf/images.json`` pinning ``ref``."""
    wf = root / "wf"
    (wf / "conf").mkdir(parents=True, exist_ok=True)
    (wf / "main.nf").write_text("// fake\n")
    entry = {"key": "0123456789ab", "tag": "1.36.2-0123456789ab", "ref": ref, "processes": [],
             "built_from": "0" * 40, "licence": "MIT"}  # fmt: skip
    (wf / "conf" / "images.json").write_text(json.dumps({"fgsea": entry}, indent=2) + "\n")
    return wf


@pytest.fixture
def wf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    tree = _wf(tmp_path)
    monkeypatch.setattr(_fr(), "workflow_main", lambda: tree / "main.nf")
    return tree


def _docker_launcher() -> Any:
    return _fr().resolve_launcher({"GIPPA_CONTAINER_RUNTIME": "docker"})


# ---------------------------------------------------------------- the argv


def test_docker_argv_has_pull_never_before_the_image(wf: Path, tmp_path: Path) -> None:
    d = tmp_path / "call"
    argv = _fr().command(_call(), launcher=_docker_launcher(), call_dir=d)
    _assert_pull_never(argv, REF)


def test_docker_argv_is_otherwise_the_m11_14c_argv(wf: Path, tmp_path: Path) -> None:
    d = tmp_path / "call"
    call = _call()
    argv = _norm(_fr().command(call, launcher=_docker_launcher(), call_dir=d))
    assert [a for a in argv if a != "--pull=never"] == _old_docker_argv(d, call), argv


def test_apptainer_argv_has_no_pull_flag(tmp_path: Path) -> None:
    """Guard: ``--pull`` is docker-only."""
    fr = _fr()
    img = tmp_path / "fgsea.img"
    launcher = fr.Launcher("apptainer", img, REF)
    argv = fr.command(_call(), launcher=launcher, call_dir=tmp_path / "call")
    assert not any(a.startswith("--pull") for a in argv), argv


_REF_PART = st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789-", min_size=1, max_size=12)


@settings(
    max_examples=25, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    host=_REF_PART,
    name=_REF_PART,
    digest=st.text(alphabet="0123456789abcdef", min_size=64, max_size=64),
    seed=st.integers(min_value=0, max_value=2**31 - 1),
)
def test_pull_never_precedes_any_pinned_ref(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, host: str, name: str, digest: str, seed: int
) -> None:
    """Property: for any pinned ref and seed, exactly one ``--pull=never`` precedes the ref."""
    ref = f"{host}.io/{name}@sha256:{digest}"
    tree = _wf(tmp_path / f"t{abs(hash((ref, seed)))}", ref)
    monkeypatch.setattr(_fr(), "workflow_main", lambda: tree / "main.nf")
    argv = _fr().command(_call(seed), launcher=_docker_launcher(), call_dir=tmp_path / "call")
    _assert_pull_never(argv, ref)


# ---------------------------------------------------------------- a missing image


def test_missing_image_fails_without_a_pull_and_names_gippa_setup(
    wf: Path, tmp_path: Path, _isolated: Path
) -> None:
    fr = _fr()
    root = tmp_path / "tmp_root"
    root.mkdir()
    with pytest.raises(fr.FgseaEngineError) as info:
        engine = fr.RFgseaEngine(_docker_launcher(), timeout_s=60.0, tmp_root=root)
        engine(_call())
    msg = str(info.value)
    assert "network pull" not in msg, msg
    assert "gippa setup" in msg, f"the missing-image error should say to run gippa setup: {msg}"
    runs = [json.loads(x) for x in _isolated.read_text().splitlines()] if _isolated.exists() else []
    assert all("--pull=never" in a for a in runs if a[:1] == ["run"]), runs
    assert list(root.iterdir()) == [], "the call directory must be removed"
