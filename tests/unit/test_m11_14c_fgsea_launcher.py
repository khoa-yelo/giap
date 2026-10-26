"""M11-14c: ``giap.io.fgsea_r`` runs R fgsea from the image when no prefix or env is found.

Spec docs/specs/M11-14-container-install.md item 1 ("fgsea (14c)": resolution order, the image
argv), item 3 (``image_cache``, shared with preflight; 14a), examples 15 and 15b, AC 7.

Interface pinned here (implementer, ``src/giap/io/fgsea_r.py``):

- ``resolve_launcher(environ)`` keeps its signature and tries, in order:
  1. ``GIAP_FGSEA_PREFIX``: an explicit prefix never falls through (bad prefix raises
     ``FgseaEngineUnavailable`` even when an image is available; example 15b);
  2. the prebuilt ``$GIAP_ROOT/micromamba/envs/giap-fgsea`` (marker = sha256 of ``fgsea.yml``);
  3. the image, through ``environ["GIAP_CONTAINER_RUNTIME"]`` (``docker`` | ``apptainer``; unset or
     any other value → ``FgseaEngineUnavailable``). The ref is ``images.json["fgsea"]["ref"]`` of
     ``workflow_main().parent / "conf" / "images.json"`` (looked up as the module attribute
     ``giap.io.fgsea_r.workflow_main``, which the tests patch with a fake workflow tree); no
     ``fgsea`` entry → ``FgseaEngineUnavailable``. Under apptainer the image is
     ``image_cache(environ, Path(environ["GIAP_DB_DIR"])) / <simpleName>.img`` (14a's helper:
     ``NXF_APPTAINER_CACHEDIR`` → ``$GIAP_ROOT/containers`` → ``<db_dir>/containers``;
     ``simpleName``: ``/`` and ``:`` → ``-``, ADR-0059); a missing ``.img`` →
     ``FgseaEngineUnavailable``. Under docker nothing local is checked (``docker run`` pulls).
  The returned ``Launcher.kind`` is ``"prefix"``, ``"docker"`` or ``"apptainer"``. The shebang
  front end is not read on the image path.
- ``command(call, launcher=None, call_dir=None)``: with no launcher or a prefix launcher, today's
  list (script + nine options). For an image launcher, the exact argv of item 1 with ``call_dir``
  (the call's temporary directory) bound:
  - docker: ``docker run --rm -u <uid>:<gid> -v <d>:<d> -w <d> -e LC_ALL=C <ref> Rscript
    /opt/<pkg>/giap_fgsea.R <options>``;
  - apptainer: ``apptainer exec --cleanenv -B <d> --pwd <d> <cache>/<simpleName>.img Rscript
    /opt/<pkg>/giap_fgsea.R <options>``.
  ``argv[0]`` may be the bare name or an absolute path to it.
- ``RFgseaEngine`` runs that argv (cwd = the call dir), so a fake ``docker``/``apptainer`` on
  ``PATH`` sees it, and the engine parses what the "image" wrote as before.

Fakes only: a fake workflow tree, fake prefixes (``helpers.fgsea_r_data``), and fake ``docker`` /
``apptainer`` binaries that hand the call to the fake R front end. No network, no real
``$GIAP_ROOT``/``$GIAP_DB_DIR``.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import stat
import sys
import tempfile
import tomllib
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers.fgsea_r_data import fake_log, fake_prefix
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
WORKFLOW = REPO / "workflow"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
REF = f"ghcr.io/khoa-yelo/{PKG}-fgsea@sha256:" + "c3" * 32
BAKED = f"/opt/{PKG}/giap_fgsea.R"
UIDGID = f"{os.getuid()}:{os.getgid()}"


def simple_name(ref: str) -> str:
    """Nextflow's apptainer cache file name for ``ref`` (ADR-0059)."""
    return ref.removeprefix("docker://").replace("/", "-").replace(":", "-") + ".img"


def _fr() -> Any:
    from giap.io import fgsea_r

    return fgsea_r


def _call(seed: int = 424242) -> Any:
    from giap.config.model import GseaConfig
    from giap.function.gsea import EngineCall, GseaJob

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
    """The nine options of ``command()`` (file names relative to the call directory)."""
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


def docker_argv(d: Path, call: Any, ref: str = REF) -> list[str]:
    return ["docker", "run", "--rm", "-u", UIDGID, "-v", f"{d}:{d}", "-w", str(d),
            "-e", "LC_ALL=C", ref, "Rscript", BAKED, *_options(call)]  # fmt: skip


def apptainer_argv(d: Path, call: Any, img: Path) -> list[str]:
    return ["apptainer", "exec", "--cleanenv", "-B", str(d), "--pwd", str(d), str(img),
            "Rscript", BAKED, *_options(call)]  # fmt: skip


def _norm(argv: list[str]) -> list[str]:
    """``argv[0]`` may be absolute; compare its name."""
    return [Path(argv[0]).name, *argv[1:]]


# ---------------------------------------------------------------- fakes

FAKE_RUNTIME = r'''#!{python} -I
"""Fake docker/apptainer (tests only): log the argv, then run the fake R front end."""
import json, os, sys
from pathlib import Path
exe = Path(sys.argv[0]).name
argv = sys.argv[1:]
with open(os.environ["FAKE_RT_LOG"], "a") as fh:
    fh.write(json.dumps({{"exe": exe, "argv": argv, "cwd": os.getcwd()}}) + "\n")
flag = "-w" if exe == "docker" else "--pwd"
if "Rscript" not in argv or flag not in argv:
    sys.stderr.write("fake " + exe + ": unexpected call\n")
    sys.exit(3)
os.chdir(argv[argv.index(flag) + 1])
fake_r = os.environ["FAKE_R"]
rest = argv[argv.index("Rscript") + 1:]
os.execv(fake_r, [fake_r, "--vanilla", *rest])
'''


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for var in (
        "GIAP_DB_DIR",
        "GIAP_ROOT",
        "NXF_APPTAINER_CACHEDIR",
        "GIAP_FGSEA_PREFIX",
        "GIAP_CONTAINER_RUNTIME",
        "GIAP_FGSEA_IMAGE",
    ):
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    for exe in ("docker", "apptainer"):
        p = bindir / exe
        p.write_text(FAKE_RUNTIME.format(python=sys.executable))
        p.chmod(p.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_RT_LOG", str(tmp_path / "runtime.log"))
    r = fake_prefix(tmp_path / "image_r")
    monkeypatch.setenv("FAKE_R", str(r / "bin" / ("R" + "script")))


def runtime_calls(tmp_path: Path) -> list[dict[str, Any]]:
    log = tmp_path / "runtime.log"
    return [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []


def fake_workflow(root: Path, images: dict[str, Any] | None, *, shebang: str | None = None) -> Path:
    """A fake ``workflow/`` tree: main.nf, bin/giap_fgsea.R, envs (repo), conf/images.json."""
    wf = root / "wf"
    (wf / "bin").mkdir(parents=True, exist_ok=True)
    (wf / "conf").mkdir(parents=True, exist_ok=True)
    (wf / "main.nf").write_text("// fake\n")
    script = wf / "bin" / "giap_fgsea.R"
    shutil.copy2(WORKFLOW / "bin" / "giap_fgsea.R", script)
    if shebang is not None:
        body = script.read_text(encoding="utf-8").split("\n", 1)[1]
        script.write_text(shebang + "\n" + body, encoding="utf-8")
    if not (wf / "envs").exists():
        (wf / "envs").symlink_to(WORKFLOW / "envs")
    write_images(wf, images)
    return wf


def write_images(wf: Path, images: dict[str, Any] | None) -> None:
    data = {} if images is None else images
    (wf / "conf" / "images.json").write_text(json.dumps(data, indent=2) + "\n")


def fgsea_entry(ref: str = REF) -> dict[str, Any]:
    return {
        "fgsea": {
            "key": "0123456789ab",
            "tag": "1.36.2-0123456789ab",
            "ref": ref,
            "processes": [],
            "built_from": "0" * 40,
            "licence": "MIT",
        }
    }


@pytest.fixture
def wf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    tree = fake_workflow(tmp_path, fgsea_entry())
    monkeypatch.setattr(_fr(), "workflow_main", lambda: tree / "main.nf")
    return tree


def built_env(root: Path) -> Path:
    env = root / "micromamba" / "envs" / "giap-fgsea"
    fake_prefix(env)
    sha = hashlib.sha256((WORKFLOW / "envs" / "fgsea.yml").read_bytes()).hexdigest()
    (env / ".giap-built").write_text(sha + "\n")
    return env


def cached_img(cache: Path, ref: str = REF) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    img = cache / simple_name(ref)
    img.write_bytes(b"fake sif\n")
    return img


# ---------------------------------------------------------------- example 15: docker


def test_example_15_docker_no_prefix_no_env(wf: Path, tmp_path: Path) -> None:
    fr = _fr()
    launcher = fr.resolve_launcher({"GIAP_CONTAINER_RUNTIME": "docker"})
    assert launcher.kind == "docker"
    d = tmp_path / "call"
    d.mkdir()
    call = _call(seed=1999999999)
    got = fr.command(call, launcher=launcher, call_dir=d)
    assert _norm(got) == docker_argv(d, call)


def test_docker_argv_runs_as_the_calling_user_and_binds_only_the_call_dir(
    wf: Path, tmp_path: Path
) -> None:
    fr = _fr()
    d = tmp_path / "x" / "call"
    d.mkdir(parents=True)
    got = fr.command(_call(), launcher=fr.resolve_launcher({"GIAP_CONTAINER_RUNTIME": "docker"}),
                     call_dir=d)  # fmt: skip
    assert got[got.index("-u") + 1] == UIDGID
    assert [got[i + 1] for i, w in enumerate(got) if w == "-v"] == [f"{d}:{d}"]
    assert got[got.index("Rscript") + 1] == BAKED


def test_docker_needs_no_local_image_and_no_root(wf: Path, tmp_path: Path) -> None:
    """Nothing local is checked under docker: no cache, no GIAP_ROOT, no db dir."""
    launcher = _fr().resolve_launcher({"GIAP_CONTAINER_RUNTIME": "docker"})
    assert launcher.kind == "docker"
    assert runtime_calls(tmp_path) == [], "resolution must not call docker"


def test_docker_ref_comes_from_images_json(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    other = f"ghcr.io/khoa-yelo/{PKG}-fgsea@sha256:" + "9f" * 32
    tree = fake_workflow(tmp_path, fgsea_entry(other))
    monkeypatch.setattr(_fr(), "workflow_main", lambda: tree / "main.nf")
    fr = _fr()
    d = tmp_path / "call"
    d.mkdir()
    got = fr.command(_call(), launcher=fr.resolve_launcher({"GIAP_CONTAINER_RUNTIME": "docker"}),
                     call_dir=d)  # fmt: skip
    assert _norm(got) == docker_argv(d, _call(), other)


# ---------------------------------------------------------------- apptainer


@pytest.mark.parametrize("where", ["nxf", "giap_root", "db_dir"])
def test_apptainer_argv_uses_the_shared_cache_helper(wf: Path, tmp_path: Path, where: str) -> None:
    fr = _fr()
    environ = {"GIAP_CONTAINER_RUNTIME": "apptainer", "GIAP_DB_DIR": str(tmp_path / "db")}
    if where == "nxf":
        cache = tmp_path / "nxf"
        environ["NXF_APPTAINER_CACHEDIR"] = str(cache)
        environ["GIAP_ROOT"] = str(tmp_path / "root")  # ignored for the cache: NXF wins
    elif where == "giap_root":
        cache = tmp_path / "root" / "containers"
        environ["GIAP_ROOT"] = str(tmp_path / "root")
    else:
        cache = tmp_path / "db" / "containers"
    img = cached_img(cache)
    from giap.io.image_key import image_cache

    assert image_cache(environ, tmp_path / "db") == cache
    launcher = fr.resolve_launcher(environ)
    assert launcher.kind == "apptainer"
    d = tmp_path / "call"
    d.mkdir()
    got = fr.command(_call(), launcher=launcher, call_dir=d)
    assert _norm(got) == apptainer_argv(d, _call(), img)


def test_apptainer_without_the_img_is_unavailable(wf: Path, tmp_path: Path) -> None:
    fr = _fr()
    (tmp_path / "nxf").mkdir()
    environ = {
        "GIAP_CONTAINER_RUNTIME": "apptainer",
        "NXF_APPTAINER_CACHEDIR": str(tmp_path / "nxf"),
    }
    with pytest.raises(fr.FgseaEngineUnavailable):
        fr.resolve_launcher(environ)


def test_apptainer_other_img_in_cache_does_not_count(wf: Path, tmp_path: Path) -> None:
    fr = _fr()
    cached_img(tmp_path / "nxf", f"ghcr.io/khoa-yelo/{PKG}-fgsea@sha256:" + "00" * 32)
    environ = {
        "GIAP_CONTAINER_RUNTIME": "apptainer",
        "NXF_APPTAINER_CACHEDIR": str(tmp_path / "nxf"),
    }
    with pytest.raises(fr.FgseaEngineUnavailable):
        fr.resolve_launcher(environ)


# ---------------------------------------------------------------- order and example 15b


def test_runtime_unset_is_unavailable_even_with_an_image(wf: Path, tmp_path: Path) -> None:
    fr = _fr()
    cached_img(tmp_path / "nxf")
    with pytest.raises(fr.FgseaEngineUnavailable):
        fr.resolve_launcher({"NXF_APPTAINER_CACHEDIR": str(tmp_path / "nxf")})


@pytest.mark.parametrize("value", ["", "podman", "singularity", "Docker"])
def test_other_runtime_values_are_unavailable(wf: Path, value: str) -> None:
    fr = _fr()
    with pytest.raises(fr.FgseaEngineUnavailable):
        fr.resolve_launcher({"GIAP_CONTAINER_RUNTIME": value})


@pytest.mark.parametrize("runtime", ["docker", "apptainer"])
def test_no_fgsea_entry_in_images_json_is_unavailable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runtime: str
) -> None:
    tree = fake_workflow(tmp_path, {})
    monkeypatch.setattr(_fr(), "workflow_main", lambda: tree / "main.nf")
    cached_img(tmp_path / "nxf")
    fr = _fr()
    with pytest.raises(fr.FgseaEngineUnavailable):
        fr.resolve_launcher(
            {"GIAP_CONTAINER_RUNTIME": runtime, "NXF_APPTAINER_CACHEDIR": str(tmp_path / "nxf")}
        )


@pytest.mark.parametrize("runtime", ["docker", "apptainer"])
def test_the_built_env_wins_over_the_image(wf: Path, tmp_path: Path, runtime: str) -> None:
    env = built_env(tmp_path / "root")
    cached_img(tmp_path / "root" / "containers")
    launcher = _fr().resolve_launcher(
        {"GIAP_CONTAINER_RUNTIME": runtime, "GIAP_ROOT": str(tmp_path / "root")}
    )
    assert launcher.kind == "prefix" and Path(launcher.path) == env


def test_the_explicit_prefix_wins_over_env_and_image(wf: Path, tmp_path: Path) -> None:
    prefix = fake_prefix(tmp_path / "explicit")
    built_env(tmp_path / "root")
    launcher = _fr().resolve_launcher(
        {
            "GIAP_FGSEA_PREFIX": str(prefix),
            "GIAP_ROOT": str(tmp_path / "root"),
            "GIAP_CONTAINER_RUNTIME": "docker",
        }
    )
    assert launcher.kind == "prefix" and Path(launcher.path) == prefix


@pytest.mark.parametrize("runtime", ["docker", "apptainer"])
@pytest.mark.parametrize("bad", ["no_dir", "no_r"])
def test_example_15b_a_bad_prefix_never_falls_through(
    wf: Path, tmp_path: Path, runtime: str, bad: str
) -> None:
    fr = _fr()
    prefix = tmp_path / "bad_prefix"
    if bad == "no_r":
        (prefix / "bin").mkdir(parents=True)
    cached_img(tmp_path / "root" / "containers")
    built_env(tmp_path / "root")
    with pytest.raises(fr.FgseaEngineUnavailable) as info:
        fr.resolve_launcher(
            {
                "GIAP_FGSEA_PREFIX": str(prefix),
                "GIAP_ROOT": str(tmp_path / "root"),
                "GIAP_CONTAINER_RUNTIME": runtime,
            }
        )
    assert "GIAP_FGSEA_PREFIX" in str(info.value)
    assert runtime_calls(tmp_path) == []


def test_unbuilt_env_message_is_unchanged_without_a_runtime(wf: Path, tmp_path: Path) -> None:
    """FarmShare (no runtime): the M8-17 message still names build_envs."""
    fr = _fr()
    with pytest.raises(fr.FgseaEngineUnavailable) as info:
        fr.resolve_launcher({"GIAP_ROOT": str(tmp_path / "empty")})
    assert "sbatch scripts/slurm/build_envs.sbatch fgsea" in str(info.value)


@pytest.mark.parametrize("runtime", ["docker", "apptainer"])
def test_image_path_does_not_read_the_shebang(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runtime: str
) -> None:
    tree = fake_workflow(tmp_path, fgsea_entry(), shebang="#!/bin/sh")
    monkeypatch.setattr(_fr(), "workflow_main", lambda: tree / "main.nf")
    cached_img(tmp_path / "nxf")
    fr = _fr()
    launcher = fr.resolve_launcher(
        {"GIAP_CONTAINER_RUNTIME": runtime, "NXF_APPTAINER_CACHEDIR": str(tmp_path / "nxf")}
    )
    assert launcher.kind == runtime
    d = tmp_path / "call"
    d.mkdir()
    assert "Rscript" in fr.command(_call(), launcher=launcher, call_dir=d)


# ---------------------------------------------------------------- command() for prefixes


def test_command_without_a_launcher_is_unchanged(wf: Path, tmp_path: Path) -> None:
    fr = _fr()
    call = _call()
    assert fr.command(call) == [str(fr.script_path()), *_options(call)]
    prefix = fr.resolve_launcher({"GIAP_FGSEA_PREFIX": str(fake_prefix(tmp_path / "p"))})
    d = tmp_path / "call"
    d.mkdir()
    assert fr.command(call, launcher=prefix, call_dir=d) == fr.command(call)


# ---------------------------------------------------------------- the engine through the image


def _engine(launcher: Any, work: Path) -> Any:
    return _fr().r_fgsea_engine(launcher, timeout_s=60, tmp_root=work)


@pytest.mark.parametrize("runtime", ["docker", "apptainer"])
def test_engine_runs_the_image_argv_and_parses_its_output(
    wf: Path, tmp_path: Path, runtime: str
) -> None:
    img = cached_img(tmp_path / "nxf")
    fr = _fr()
    launcher = fr.resolve_launcher(
        {"GIAP_CONTAINER_RUNTIME": runtime, "NXF_APPTAINER_CACHEDIR": str(tmp_path / "nxf")}
    )
    work = tmp_path / "work"
    work.mkdir()
    engine = _engine(launcher, work)
    call = _call(seed=77)
    result = engine(call)
    assert [r.term for r in result["main:cog"]] == ["C1", "C2"]
    assert [r.term for r in result["main:kegg"]] == ["K1"]
    assert engine.provenance is not None and engine.provenance.seed == 77
    (logged,) = runtime_calls(tmp_path)
    assert logged["exe"] == runtime
    d = Path(logged["cwd"])
    assert d.parent == work
    want = docker_argv(d, call) if runtime == "docker" else apptainer_argv(d, call, img)
    assert [runtime, *logged["argv"]] == want
    assert list(work.iterdir()) == [], "the call directory is removed"
    # the fake R in the "image" saw the baked-in path, not the checkout's script
    r_log = fake_log(Path(os.environ["FAKE_R"]).parents[1])
    assert r_log["argv"][1] == BAKED


def test_engine_image_failure_is_an_engine_error(wf: Path, tmp_path: Path) -> None:
    from helpers.fgsea_r_data import set_behaviour

    set_behaviour(Path(os.environ["FAKE_R"]).parents[1], {"mode": "fail", "code": 1,
                                                           "stderr": "boom\n"})  # fmt: skip
    fr = _fr()
    work = tmp_path / "work"
    work.mkdir()
    engine = _engine(fr.resolve_launcher({"GIAP_CONTAINER_RUNTIME": "docker"}), work)
    with pytest.raises(fr.FgseaEngineError) as info:
        engine(_call())
    assert "boom" in str(info.value)
    assert list(work.iterdir()) == []


def test_engine_docker_missing_on_path_is_unavailable(
    wf: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "fakebin" / "docker").unlink()
    empty = tmp_path / "emptybin"
    empty.mkdir()
    monkeypatch.setenv("PATH", f"{tmp_path / 'fakebin'}{os.pathsep}{empty}")
    fr = _fr()
    work = tmp_path / "work"
    work.mkdir()
    with pytest.raises(fr.FgseaEngineUnavailable):
        launcher = fr.resolve_launcher({"GIAP_CONTAINER_RUNTIME": "docker"})
        _engine(launcher, work)(_call())


# ---------------------------------------------------------------- the order as a property


PREFIX_STATES = st.sampled_from(["none", "good", "no_dir", "no_r"])
RUNTIMES = st.sampled_from([None, "docker", "apptainer", "podman"])


@settings(
    max_examples=60, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    prefix_state=PREFIX_STATES,
    env=st.booleans(),
    runtime=RUNTIMES,
    entry=st.booleans(),
    img=st.booleans(),
)
def test_property_resolution_order(
    wf: Path, prefix_state: str, env: bool, runtime: str | None, entry: bool, img: bool
) -> None:
    fr = _fr()
    write_images(wf, fgsea_entry() if entry else {})
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        environ: dict[str, str] = {"GIAP_ROOT": str(t / "root"), "GIAP_DB_DIR": str(t / "db")}
        (t / "root").mkdir()
        if prefix_state == "good":
            prefix = fake_prefix(t / "explicit")
        elif prefix_state == "no_r":
            prefix = t / "explicit"
            (prefix / "bin").mkdir(parents=True)
        else:
            prefix = t / "explicit"
        if prefix_state != "none":
            environ["GIAP_FGSEA_PREFIX"] = str(prefix)
        env_prefix = built_env(t / "root") if env else None
        if img:
            cached_img(t / "root" / "containers")
        if runtime is not None:
            environ["GIAP_CONTAINER_RUNTIME"] = runtime

        if prefix_state == "good":
            want: tuple[str, Path | None] | None = ("prefix", prefix)
        elif prefix_state != "none":
            want = None
        elif env_prefix is not None:
            want = ("prefix", env_prefix)
        elif runtime == "docker" and entry:
            want = ("docker", None)
        elif runtime == "apptainer" and entry and img:
            want = ("apptainer", None)
        else:
            want = None

        if want is None:
            with pytest.raises(fr.FgseaEngineUnavailable):
                fr.resolve_launcher(environ)
        else:
            got = fr.resolve_launcher(environ)
            assert got.kind == want[0]
            if want[1] is not None:
                assert Path(got.path) == want[1]
