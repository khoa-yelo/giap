"""M11-14c: ``gippa run`` exports ``GIPPA_CONTAINER_RUNTIME``; preflight's fgsea row.

Spec docs/specs/M11-14-container-install.md item 1 (resolution step 3: "the runtime that ``gippa
run`` passes in ``GIPPA_CONTAINER_RUNTIME`` (``docker`` | ``apptainer``; unset →
``FgseaEngineUnavailable``). The cache directory comes from the same helper as preflight"), item 3
(the table's ``fgsea`` row: "prefix → env → image" under ``docker`` and plain ``apptainer``, "as
today" under ``conda``/``farmshare*``; "plus fgsea through ``images.json``").

Interface pinned here:

- ``gippa run`` (Nextflow executor) starts ``nextflow`` with an environment in which
  ``GIPPA_CONTAINER_RUNTIME`` equals 14f1's ``container_runtime(profile)``; when that is ``None``
  the variable is absent (an inherited value is dropped). With ``--db-dir D`` the child also has
  ``GIPPA_DB_DIR=D`` (absolute), so the function step's ``image_cache(environ, GIPPA_DB_DIR)`` is
  the directory preflight checked. ``os.environ`` of the caller is not modified.
- ``gippa.preflight.environment.missing_tools(tools, profile, image_cache, env_root, local_tools)``
  (14f1 signature) for ``fgsea``: under ``docker``, present iff the prebuilt env is built or
  ``images.json`` has an ``fgsea`` entry; under plain ``apptainer``, present iff the env is built
  or ``<image_cache>/<simpleName(ref)>.img`` exists; under ``conda``/``farmshare*`` profiles only
  the env counts (as today). ``images.json`` is read from ``workflow_main().parent / "conf"``
  through the module attribute ``gippa.preflight.environment.workflow_main`` (as 14f1 does for
  ``containers.config``). A missing fgsea image under apptainer names ``<pkg> setup --images
  --profile apptainer`` in ``missing_tool_reasons``.

Fakes only: a fake ``nextflow`` that dumps its environment (preflight stubbed for those tests),
a fake workflow tree, temporary DB
and cache directories. Never the real ``$GIPPA_ROOT``/``$GIPPA_DB_DIR``; no network.
"""

from __future__ import annotations

import hashlib
import json
import os
import random
import stat
import tomllib
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import write_setup_manifest

from gippa.annotate.tools import TOOLS

REPO = Path(__file__).resolve().parents[2]
WORKFLOW = REPO / "workflow"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
REF = f"ghcr.io/khoa-yelo/{PKG}-fgsea@sha256:" + "c3" * 32
SETUP_IMAGES = f"{PKG} setup --images --profile apptainer"
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


def simple_name(ref: str) -> str:
    return ref.removeprefix("docker://").replace("/", "-").replace(":", "-") + ".img"


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for var in (
        "GIPPA_DB_DIR",
        "GIPPA_ROOT",
        "NXF_APPTAINER_CACHEDIR",
        "GIPPA_FGSEA_PREFIX",
        "GIPPA_CONTAINER_RUNTIME",
    ):
        monkeypatch.delenv(var, raising=False)


# ---------------------------------------------------------------- gippa run → nextflow env

FAKE_NF = """#!/bin/sh
env -0 > "$FAKE_NF_ENV"
exit 0
"""


@pytest.fixture
def fake_nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A fake ``nextflow``; preflight is stubbed out (``gippa.cli._run_preflight`` returns no
    normalised inputs) because ``gippa run`` accepts only real profiles, under which the session
    has no tools. The runtime is a launch concern, not a preflight one."""
    import gippa.cli

    monkeypatch.setattr(gippa.cli, "_run_preflight", lambda cfg, *a, **k: ({}, cfg))
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    out = tmp_path / "nf_env.txt"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ENV", str(out))
    return out


def _inputs(tmp_path: Path) -> dict[str, Path]:
    rng = random.Random(0)
    d = tmp_path / "in"
    d.mkdir(parents=True, exist_ok=True)
    g = d / "g.fa"
    g.write_text(
        f">chr1 [topology=circular]\n{''.join(rng.choice('ACGT') for _ in range(120_000))}\n"
    )
    p = d / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    return {"genome": g, "proteins": p}


def _gippa_run(tmp_path: Path, profile: str, *extra: str) -> tuple[int, str]:
    from gippa.cli import cli

    f = _inputs(tmp_path)
    argv = ["run", "--genome", str(f["genome"]), "--proteins", str(f["proteins"]),
            "--outdir", str(tmp_path / "OUT"), "--profile", profile, *extra]  # fmt: skip
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output


def _nf_env(path: Path) -> dict[str, str]:
    assert path.is_file(), "nextflow was not launched"
    out = {}
    for line in path.read_text().split("\0"):
        k, sep, v = line.partition("=")
        if sep:
            out[k] = v
    return out


@pytest.mark.parametrize(
    ("profile", "want"),
    [("docker", "docker"), ("apptainer", "apptainer"), ("farmshare,apptainer", "apptainer"),
     ("conda", None), ("standard", None)],
)  # fmt: skip
def test_gippa_run_exports_the_runtime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_nf: Path, profile: str, want: str | None
) -> None:
    db = write_setup_manifest(tmp_path / "db")
    code, text = _gippa_run(tmp_path, profile, "--db-dir", str(db))
    assert fake_nf.is_file(), (code, text)  # launched; the tail may fail on a fake run
    env = _nf_env(fake_nf)
    assert env.get("GIPPA_CONTAINER_RUNTIME") == want
    assert "GIPPA_CONTAINER_RUNTIME" not in os.environ, "gippa run must not mutate os.environ"


@pytest.mark.parametrize(
    ("inherited", "profile", "want"),
    [("apptainer", "docker", "docker"), ("docker", "conda", None)],
)
def test_gippa_run_overrides_an_inherited_runtime(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fake_nf: Path,
    inherited: str,
    profile: str,
    want: str | None,
) -> None:
    monkeypatch.setenv("GIPPA_CONTAINER_RUNTIME", inherited)
    db = write_setup_manifest(tmp_path / "db")
    code, text = _gippa_run(tmp_path, profile, "--db-dir", str(db))
    assert fake_nf.is_file(), (code, text)  # launched; the tail may fail on a fake run
    assert _nf_env(fake_nf).get("GIPPA_CONTAINER_RUNTIME") == want
    assert os.environ["GIPPA_CONTAINER_RUNTIME"] == inherited


def test_gippa_run_passes_the_db_dir_for_the_image_cache(tmp_path: Path, fake_nf: Path) -> None:
    """No NXF cache, no GIPPA_ROOT: the step's cache is ``<db_dir>/containers``, as preflight's."""
    from gippa.io.image_key import image_cache

    db = write_setup_manifest(tmp_path / "db")
    code, text = _gippa_run(tmp_path, "apptainer", "--db-dir", str(db))
    assert fake_nf.is_file(), (code, text)  # launched; the tail may fail on a fake run
    env = _nf_env(fake_nf)
    assert Path(env["GIPPA_DB_DIR"]) == db.absolute()
    assert image_cache(env, Path(env["GIPPA_DB_DIR"])) == image_cache({}, db)


def test_gippa_run_keeps_the_nxf_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_nf: Path
) -> None:
    from gippa.io.image_key import image_cache

    monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", str(tmp_path / "nxf"))
    db = write_setup_manifest(tmp_path / "db")
    code, text = _gippa_run(tmp_path, "apptainer", "--db-dir", str(db))
    assert fake_nf.is_file(), (code, text)  # launched; the tail may fail on a fake run
    env = _nf_env(fake_nf)
    assert image_cache(env, Path(env["GIPPA_DB_DIR"])) == tmp_path / "nxf"


# ---------------------------------------------------------------- preflight: the fgsea row


def _env_mod() -> Any:
    from gippa.preflight import environment

    return environment


def _fn(name: str) -> Any:
    fn = getattr(_env_mod(), name, None)
    assert fn is not None, f"gippa.preflight.environment.{name} is not implemented"
    return fn


def _fake_tree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, with_fgsea: bool) -> Path:
    wf = tmp_path / "wf"
    (wf / "conf").mkdir(parents=True)
    (wf / "main.nf").write_text("// fake\n")
    (wf / "envs").symlink_to(WORKFLOW / "envs")
    (wf / "bin").symlink_to(WORKFLOW / "bin")
    (wf / "conf" / "containers.config").write_text(
        (WORKFLOW / "conf" / "containers.config").read_text(encoding="utf-8")
    )
    images = json.loads((WORKFLOW / "conf" / "images.json").read_text())
    images.pop("fgsea", None)
    if with_fgsea:
        images["fgsea"] = {"key": "0123456789ab", "tag": "1.36.2-0123456789ab", "ref": REF,
                           "processes": [], "built_from": "0" * 40, "licence": "MIT"}  # fmt: skip
    (wf / "conf" / "images.json").write_text(json.dumps(images, indent=2) + "\n")
    monkeypatch.setattr(_env_mod(), "workflow_main", lambda: wf / "main.nf")
    return wf


def _cache(tmp_path: Path, img: bool) -> Path:
    c = tmp_path / "cache"
    c.mkdir(exist_ok=True)
    if img:
        (c / simple_name(REF)).write_bytes(b"fake sif\n")
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


def _missing(tmp_path: Path, profile: str, *, img: bool, built: bool) -> list[str]:
    fn = _fn("missing_tools")
    return list(fn([TOOLS["fgsea"]], profile, _cache(tmp_path, img), _envs(tmp_path, built), None))


@pytest.mark.parametrize(
    ("profile", "entry", "img", "built", "missing"),
    [
        # docker: env, else the images.json entry (docker run pulls it)
        ("docker", True, False, False, False),
        ("docker", False, False, False, True),
        ("docker", False, False, True, False),
        # plain apptainer: env, else the .img in the cache
        ("apptainer", True, True, False, False),
        ("apptainer", True, False, False, True),
        ("apptainer", False, True, False, True),
        ("apptainer", True, False, True, False),
        # conda / farmshare*: as today (the env only)
        ("farmshare,apptainer", True, True, False, True),
        ("farmshare,apptainer", True, True, True, False),
        ("farmshare_local,apptainer", True, True, False, True),
        ("conda", True, True, False, True),
        ("conda", False, False, True, False),
        # test/stub need nothing
        ("test", False, False, False, False),
        ("test,stub", False, False, False, False),
    ],
)
def test_fgsea_row_of_the_item_3_table(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    profile: str,
    entry: bool,
    img: bool,
    built: bool,
    missing: bool,
) -> None:
    _fake_tree(tmp_path, monkeypatch, entry)
    assert _missing(tmp_path, profile, img=img, built=built) == (["fgsea"] if missing else [])


def test_fgsea_missing_image_under_apptainer_names_setup_images(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _fake_tree(tmp_path, monkeypatch, True)
    reasons = _fn("missing_tool_reasons")(
        [TOOLS["fgsea"]], "apptainer", _cache(tmp_path, False), _envs(tmp_path, False), None
    )
    assert set(reasons) == {"fgsea"}
    assert SETUP_IMAGES in reasons["fgsea"], reasons


def test_fgsea_under_docker_is_not_reported_as_conda_only_once_pinned(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _fake_tree(tmp_path, monkeypatch, True)
    reasons = _fn("missing_tool_reasons")(
        [TOOLS["fgsea"]], "docker", _cache(tmp_path, False), _envs(tmp_path, False), None
    )
    assert reasons == {}


def test_fgsea_has_no_containers_config_line() -> None:
    """The ref lives only in images.json: no ``withName: 'FGSEA'`` selector (item 1)."""
    text = (WORKFLOW / "conf" / "containers.config").read_text(encoding="utf-8")
    assert "FGSEA" not in text
    data = json.loads((WORKFLOW / "conf" / "images.json").read_text())
    if "fgsea" in data:
        assert data["fgsea"]["processes"] == []


# ---------------------------------------------------------------- review nits (m10, local executor)


def test_launch_env_makes_a_relative_db_dir_absolute(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """m10: Nextflow runs with ``cwd=outdir``; a relative ``GIPPA_DB_DIR`` would point elsewhere."""
    from gippa.runner.nextflow import launch_env

    monkeypatch.chdir(tmp_path)
    env = launch_env({"GIPPA_DB_DIR": "stale"}, "apptainer", Path("rel") / "db")
    assert env["GIPPA_DB_DIR"] == str(tmp_path / "rel" / "db")
    assert Path(env["GIPPA_DB_DIR"]).is_absolute()


@pytest.mark.parametrize("how", ["flag", "inherited"])
def test_gippa_run_exports_an_absolute_db_dir_for_a_relative_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_nf: Path, how: str
) -> None:
    write_setup_manifest(tmp_path / "db")
    monkeypatch.chdir(tmp_path)
    extra: tuple[str, ...] = ()
    if how == "flag":
        extra = ("--db-dir", "db")
    else:
        monkeypatch.setenv("GIPPA_DB_DIR", "db")
    code, text = _gippa_run(tmp_path, "apptainer", *extra)
    assert fake_nf.is_file(), (code, text)
    assert _nf_env(fake_nf)["GIPPA_DB_DIR"] == str(tmp_path / "db")


def test_local_executor_does_not_honour_an_inherited_runtime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``--executor local`` has no profile, so its runtime is ``None``: the function step must not
    see a ``GIPPA_CONTAINER_RUNTIME`` left in the user's shell. The spec says only that ``gippa
    run``
    passes the runtime (item 1, step 3); this pins the reading that ``gippa run`` always sets or
    clears it, under either executor. The caller's value is back after the run."""
    import gippa.cli
    from gippa.runner import local

    seen: dict[str, str | None] = {}

    def run_local(*a: Any, **k: Any) -> None:
        seen["runtime"] = os.environ.get("GIPPA_CONTAINER_RUNTIME")

    monkeypatch.setattr(gippa.cli, "_local_roots", lambda *a, **k: {})
    monkeypatch.setattr(local, "plan_local", lambda *a, **k: None)
    monkeypatch.setattr(local, "run_local", run_local)
    monkeypatch.setenv("GIPPA_CONTAINER_RUNTIME", "docker")
    _gippa_run(tmp_path, "docker", "--executor", "local")
    assert "runtime" in seen, "run_local was not reached"
    assert seen["runtime"] is None
    assert os.environ["GIPPA_CONTAINER_RUNTIME"] == "docker"
