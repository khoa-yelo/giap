"""#425 review gaps in the five-step install path (docs/specs/M11-14-container-install.md, "The
target 1.0 user path" and item 3; review B1 and mismatch claims 1, 4, 5; nit N8).

G1, one image cache everywhere. Setup's image pull, setup's own ``ToolRunner`` (the Bakta and
antiSMASH downloaders that run inside images), run-mode preflight and the Nextflow launch all use
one directory: ``setup --cache-dir`` when given (setup only; the spec defines no such option for
run/preflight), else ``NXF_APPTAINER_CACHEDIR``, else ``$GIPPA_ROOT/containers``, else
``<db_dir>/containers`` (``gippa.io.image_key.image_cache``, Q-L11). "``gippa run`` passes the same
directory to Nextflow (``apptainer.cacheDir``)".

G2, images before databases. ``setup --images`` pulls the images before it runs a downloader that
needs one, so on an Apptainer-only host (no docker) the first run already records Bakta (and
antiSMASH with ``--antismash``) as ``ok``.

G3, ``gippa preflight --profile P``: the same choices and validation as ``gippa run --profile``,
then the tool check exactly as ``gippa run`` does it (step 4 of the user path). Without
``--profile`` nothing changes.

G4, ``gippa run`` without Nextflow on PATH points a pip user at the install docs, not AGENTS.md.

Interface pinned here (readings chosen where the spec leaves room):

- ``ToolRunner(*, db_dir=None, cache_dir=None)``: with ``cache_dir`` the apptainer image is looked
  up there (and nowhere else); without it, as today (``image_cache(os.environ, db_dir)``).
  ``gippa setup`` passes ``cache_dir`` = ``--cache-dir`` or ``image_cache(os.environ, db_dir)``.
- ``gippa.runner.nextflow.launch_env(environ, runtime, db_dir)``: when ``db_dir`` is given and
  ``environ`` has no non-empty ``NXF_APPTAINER_CACHEDIR``, the returned env sets it to
  ``image_cache(environ, absolute db_dir)`` (so ``$GIPPA_ROOT/containers`` when ``$GIPPA_ROOT`` is
  set, else ``<abs db_dir>/containers``). A non-empty inherited value is kept unchanged.
  ``workflow/nextflow.config`` reads that variable, so the launcher's env is what is tested; no
  real Nextflow runs.
- ``gippa preflight --profile P`` reuses ``gippa run``'s profile validation (``unknown profile(s)``
  usage error) and calls ``run_preflight(..., profile=P)`` with the DB directory from
  ``$GIPPA_DB_DIR`` and ``$GIPPA_ROOT`` as ``gippa run`` does; its ``environment.tools`` result
  equals ``run_preflight(..., "run", ..., db_dir, profile=P)``'s.
- The "nextflow not found" message of ``gippa run`` contains ``get-started/install`` and
  ``<pkg> setup --images``, and not ``AGENTS.md``.

Fakes only: a logging fake ``apptainer`` (no docker anywhere on PATH), a fake ``nextflow`` that
dumps its environment, temporary DB and cache directories. Never the real ``$GIPPA_ROOT`` or
``$GIPPA_DB_DIR``; no network.
"""

from __future__ import annotations

import json
import os
import random
import shutil
import stat
import tomllib
from collections import namedtuple
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import write_setup_manifest
from helpers.prereqs import link_prereqs
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.annotate.tools import TOOLS, read_containers_config, required_tools
from gippa.config import resolve_config

REPO = Path(__file__).resolve().parents[2]
CONTAINERS_CONFIG = REPO / "workflow" / "conf" / "containers.config"
NEXTFLOW_CONFIG = REPO / "workflow" / "nextflow.config"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"
LOW = ["--set", "preflight.min_n50_warn=1", "--set", "preflight.min_largest_contig_warn=1"]
DBS_BAKTA = ["--bakta", "light", "--no-genomad", "--no-promotech-models", "--no-uniop"]

# Shell builtins only (plus a linked mkdir): PATH holds nothing but the fakes, so no docker.
FAKE_APPTAINER = r"""#!/bin/sh
printf '%s\n' "${0##*/} $*" >> "$FAKE_RUNTIME_LOG"
last=""; prev=""; out=""; grab=0
for a in "$@"; do
    prev=$last; last=$a
    if [ "$grab" = 1 ]; then out=$a; grab=0; fi
    [ "$a" = "--output" ] && grab=1
done
case "$1" in
pull)
    printf 'fake sif %s\n' "$last" > "$prev"
    exit 0;;
exec)
    shift
    while [ "$1" = "--bind" ] || [ "$1" = "-B" ]; do shift 2; done
    if [ ! -f "$1" ]; then echo "FATAL: could not open image $1" >&2; exit 255; fi
    if [ -n "$out" ]; then
        mkdir -p "$out"
        printf '{"major": 6, "minor": 0, "type": "light"}\n' > "$out/version.json"
    fi
    exit 0;;
esac
echo "fake apptainer: unexpected call: $*" >&2
exit 3
"""

FAKE_NF = """#!/bin/sh
env -0 > "$FAKE_NF_ENV"
exit 0
"""


def simple_name(ref: str) -> str:
    """Nextflow's apptainer cache file name for ``ref`` (ADR-0059)."""
    return ref.removeprefix("docker://").replace("/", "-").replace(":", "-") + ".img"


def pins() -> dict[str, str]:
    return {
        p.lower(): r
        for p, r in read_containers_config(CONTAINERS_CONFIG).items()
        if p.lower() in TOOLS
    }


def _exe(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for var in (
        "GIPPA_DB_DIR",
        "GIPPA_ROOT",
        "GIPPA_REPO",
        "NXF_APPTAINER_CACHEDIR",
        "GIPPA_CONTAINER_RUNTIME",
        "GIPPA_FGSEA_PREFIX",
    ):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(tmp_path / "no_mamba"))
    monkeypatch.setenv("FAKE_RUNTIME_LOG", str(tmp_path / "runtime.log"))


@pytest.fixture
def apptainer_only(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """PATH = a fake ``apptainer`` and ``mkdir``; no docker; disk space is plenty."""
    from gippa.setupdb import run

    mkdir = shutil.which("mkdir")
    assert mkdir is not None
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    _exe(bindir / "apptainer", FAKE_APPTAINER)
    (bindir / "mkdir").symlink_to(mkdir)
    link_prereqs(bindir)  # M11-15d: setup checks java and finds nextflow (no launcher download)
    monkeypatch.setenv("PATH", str(bindir))
    assert shutil.which("docker") is None
    usage = namedtuple("usage", "total used free")
    monkeypatch.setattr(run.shutil, "disk_usage", lambda _p: usage(1 << 50, 0, 1 << 50))
    return bindir


def calls(tmp_path: Path, prefix: str = "") -> list[list[str]]:
    log = tmp_path / "runtime.log"
    lines = log.read_text().splitlines() if log.exists() else []
    return [ln.split() for ln in lines if ln.startswith(prefix)]


def setup_cli(args: list[str]) -> Any:
    from gippa.cli import cli

    return CliRunner().invoke(cli, ["setup", *args], catch_exceptions=False)


def manifest_items(db: Path) -> dict[str, dict[str, Any]]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    return {item["name"]: item for item in data["items"]}


def nextflow_cache(env: dict[str, str]) -> Path | None:
    """What ``workflow/nextflow.config``'s apptainer block resolves ``apptainer.cacheDir`` to
    (``NXF_APPTAINER_CACHEDIR ?: GIPPA_ROOT + '/containers'``; Groovy treats "" as false)."""
    if env.get("NXF_APPTAINER_CACHEDIR"):
        return Path(env["NXF_APPTAINER_CACHEDIR"])
    if env.get("GIPPA_ROOT"):
        return Path(env["GIPPA_ROOT"] + "/containers")
    return None


def test_nextflow_config_reads_the_env_var() -> None:
    """The emulation above matches the config line it stands for."""
    text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
    line = next(ln for ln in text.splitlines() if "apptainer.cacheDir" in ln)
    assert "System.getenv('NXF_APPTAINER_CACHEDIR')" in line, line


# ---------------------------------------------------------------- G1a: setup's ToolRunner


def test_toolrunner_cache_dir_keyword_wins(tmp_path: Path, apptainer_only: Path) -> None:
    from gippa.setupdb.run import ToolRunner

    db, cache = tmp_path / "db", tmp_path / "cache"
    img = cache / simple_name(pins()["bakta"])
    img.parent.mkdir(parents=True)
    img.write_bytes(b"sif")
    cmd = ToolRunner(db_dir=db, cache_dir=cache).command("bakta", ["bakta_db", "--help"])
    assert cmd[0] == "apptainer" and str(img) in cmd, cmd


def test_toolrunner_cache_dir_is_the_only_place_looked(
    tmp_path: Path, apptainer_only: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An image only in the env-var cache does not count once ``cache_dir`` is explicit."""
    from gippa.setupdb.run import ToolRunner, ToolUnavailable

    other = tmp_path / "nxf"
    (other / simple_name(pins()["bakta"])).parent.mkdir(parents=True)
    (other / simple_name(pins()["bakta"])).write_bytes(b"sif")
    monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", str(other))
    with pytest.raises(ToolUnavailable):
        ToolRunner(db_dir=tmp_path / "db", cache_dir=tmp_path / "cache").command(
            "bakta", ["bakta_db", "--help"]
        )


def test_setup_cache_dir_reaches_the_bakta_downloader(tmp_path: Path, apptainer_only: Path) -> None:
    """B1: the image is already in ``--cache-dir`` (an earlier ``--images`` run); a plain setup
    rerun with the same ``--cache-dir`` downloads Bakta through it instead of reporting absent."""
    db, cache = tmp_path / "db", tmp_path / "cache"
    img = cache / simple_name(pins()["bakta"])
    img.parent.mkdir(parents=True)
    img.write_bytes(b"sif")
    res = setup_cli(["--db-dir", str(db), "--cache-dir", str(cache), *DBS_BAKTA])
    items = manifest_items(db)
    assert items["bakta_light"]["status"] == "ok", (res.output, items["bakta_light"])
    execs = calls(tmp_path, "apptainer exec")
    assert execs and str(img) in execs[0], calls(tmp_path)
    assert not (db / "containers").exists()


# ---------------------------------------------------------------- G2: images before databases


def _first(log: list[list[str]], verb: str, needle: str) -> int:
    for i, c in enumerate(log):
        if c[1] == verb and any(needle in a for a in c):
            return i
    raise AssertionError(f"no apptainer {verb} with {needle}: {log}")


def test_images_then_bakta_on_an_apptainer_only_host(tmp_path: Path, apptainer_only: Path) -> None:
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *DBS_BAKTA, "--images", "--profile", "apptainer"])
    items = manifest_items(db)
    assert items["bakta_light"]["status"] == "ok", (res.output, items["bakta_light"])
    assert items["image:bakta"]["status"] == "ok"
    log = calls(tmp_path, "apptainer")
    img = str(db / "containers" / simple_name(pins()["bakta"]))
    assert _first(log, "pull", img) < _first(log, "exec", img), log
    assert res.exit_code == 0, res.output


def test_images_then_bakta_with_cache_dir(tmp_path: Path, apptainer_only: Path) -> None:
    """G1 and G2 together: pull and downloader both use ``--cache-dir``."""
    db, cache = tmp_path / "db", tmp_path / "cache"
    res = setup_cli(["--db-dir", str(db), *DBS_BAKTA, "--images", "--profile", "apptainer",
                     "--cache-dir", str(cache)])  # fmt: skip
    assert manifest_items(db)["bakta_light"]["status"] == "ok", res.output
    img = str(cache / simple_name(pins()["bakta"]))
    log = calls(tmp_path, "apptainer")
    assert _first(log, "pull", img) < _first(log, "exec", img), log
    assert not (db / "containers").exists()


def test_images_then_bakta_with_the_env_cache(
    tmp_path: Path, apptainer_only: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db, cache = tmp_path / "db", tmp_path / "nxf"
    monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", str(cache))
    res = setup_cli(["--db-dir", str(db), *DBS_BAKTA, "--images", "--profile", "apptainer"])
    assert manifest_items(db)["bakta_light"]["status"] == "ok", res.output
    assert (cache / simple_name(pins()["bakta"])).is_file()


def test_images_then_antismash_on_an_apptainer_only_host(
    tmp_path: Path, apptainer_only: Path
) -> None:
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), "--bakta", "skip", "--no-genomad",
                     "--no-promotech-models", "--no-uniop", "--antismash",
                     "--images", "all", "--profile", "apptainer"])  # fmt: skip
    items = manifest_items(db)
    assert items["antismash"]["status"] == "ok", (res.output, items["antismash"])
    log = calls(tmp_path, "apptainer")
    img = str(db / "containers" / simple_name(pins()["antismash"]))
    assert _first(log, "pull", img) < _first(log, "exec", img), log


def test_without_images_bakta_stays_absent_on_an_apptainer_only_host(
    tmp_path: Path, apptainer_only: Path
) -> None:
    """Control: no image, no docker, no env: the downloader cannot run (unchanged)."""
    db = tmp_path / "db"
    setup_cli(["--db-dir", str(db), *DBS_BAKTA])
    rec = manifest_items(db)["bakta_light"]
    assert rec["status"] == "absent" and "setup --images" in rec["message"], rec


# ---------------------------------------------------------------- G1c: the Nextflow launch env


def _launch_env() -> Any:
    from gippa.runner.nextflow import launch_env

    return launch_env


@pytest.mark.parametrize("runtime", ["apptainer", "docker", None])
def test_launch_env_exports_the_db_dir_cache(tmp_path: Path, runtime: str | None) -> None:
    env = _launch_env()({}, runtime, tmp_path / "db")
    assert env.get("NXF_APPTAINER_CACHEDIR") == str(tmp_path / "db" / "containers"), env


def test_launch_env_relative_db_dir_gives_an_absolute_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    env = _launch_env()({}, "apptainer", Path("rel") / "db")
    assert env.get("NXF_APPTAINER_CACHEDIR") == str(tmp_path / "rel" / "db" / "containers")


def test_launch_env_keeps_an_inherited_cache() -> None:
    env = _launch_env()({"NXF_APPTAINER_CACHEDIR": "/fake/nxf"}, "apptainer", Path("/fake/db"))
    assert env["NXF_APPTAINER_CACHEDIR"] == "/fake/nxf"


def test_launch_env_empty_inherited_cache_counts_as_unset() -> None:
    env = _launch_env()({"NXF_APPTAINER_CACHEDIR": ""}, "apptainer", Path("/fake/db"))
    assert env.get("NXF_APPTAINER_CACHEDIR") == "/fake/db/containers"


def test_launch_env_gippa_root_cache() -> None:
    env = _launch_env()({"GIPPA_ROOT": "/fake/root"}, "apptainer", Path("/fake/db"))
    assert nextflow_cache(env) == Path("/fake/root/containers")


def test_launch_env_does_not_modify_environ() -> None:
    environ = {"PATH": "/x"}
    _launch_env()(environ, "apptainer", Path("/fake/db"))
    assert environ == {"PATH": "/x"}


_val = st.one_of(st.none(), st.just(""), st.sampled_from(["/a/nxf", "/b/root"]))


@settings(max_examples=60, deadline=None)
@given(nxf=_val, root=_val, runtime=st.sampled_from(["apptainer", "docker", None]))
def test_property_nextflow_sees_preflights_cache(
    nxf: str | None, root: str | None, runtime: str | None
) -> None:
    """For every environment, the directory Nextflow's config resolves from the launch env is
    ``image_cache`` of the caller's environment: the one preflight checked."""
    from gippa.io.image_key import image_cache

    environ: dict[str, str] = {}
    if nxf is not None:
        environ["NXF_APPTAINER_CACHEDIR"] = nxf
    if root is not None:
        environ["GIPPA_ROOT"] = root
    db = Path("/fake/db")
    env = _launch_env()(environ, runtime, db)
    assert nextflow_cache(env) == image_cache(environ, db)


# ---------------------------------------------------------------- G1c through `gippa run`


@pytest.fixture
def fake_nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A fake ``nextflow`` (preflight stubbed: a launch concern, as in the 14c tests)."""
    import gippa.cli

    monkeypatch.setattr(gippa.cli, "_run_preflight", lambda cfg, *a, **k: ({}, cfg))
    bindir = tmp_path / "nfbin"
    bindir.mkdir()
    _exe(bindir / "nextflow", FAKE_NF)
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


@pytest.mark.parametrize("how", ["flag", "inherited"])
def test_gippa_run_exports_the_db_dir_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_nf: Path, how: str
) -> None:
    db = write_setup_manifest(tmp_path / "db")
    extra: tuple[str, ...] = ("--db-dir", str(db)) if how == "flag" else ()
    if how == "inherited":
        monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    code, text = _gippa_run(tmp_path, "apptainer", *extra)
    env = _nf_env(fake_nf)
    assert env.get("NXF_APPTAINER_CACHEDIR") == str(db.absolute() / "containers"), (code, text)
    assert "NXF_APPTAINER_CACHEDIR" not in os.environ


def test_gippa_run_farmshare_with_gippa_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_nf: Path
) -> None:
    """FarmShare is unchanged: ``$GIPPA_ROOT/containers``."""
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "root"))
    db = write_setup_manifest(tmp_path / "db")
    _gippa_run(tmp_path, "farmshare,apptainer", "--db-dir", str(db))
    assert nextflow_cache(_nf_env(fake_nf)) == tmp_path / "root" / "containers"


def test_setup_preflight_and_launch_agree(
    tmp_path: Path, apptainer_only: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The user path without ``$GIPPA_ROOT`` or ``NXF_APPTAINER_CACHEDIR``: the images setup
    pulled are where run-mode preflight looks and where Nextflow will look."""
    from gippa.io.image_key import image_cache
    from gippa.preflight.environment import missing_tools

    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), "--bakta", "skip", "--no-genomad",
                     "--no-promotech-models", "--no-uniop", "--images", "all",
                     "--profile", "apptainer"])  # fmt: skip
    assert res.exit_code == 0, res.output
    pinned = [TOOLS[t] for t in pins()]
    assert missing_tools(pinned, "apptainer", image_cache(os.environ, db), None, None) == []
    env = _launch_env()(dict(os.environ), "apptainer", db)
    got = nextflow_cache(env)
    assert got is not None and (got / simple_name(pins()["bakta"])).is_file(), env


# ---------------------------------------------------------------- G3: gippa preflight --profile


def _preflight_cli(tmp_path: Path, *extra: str) -> tuple[int, str, Path]:
    from gippa.cli import cli

    f = _inputs(tmp_path)
    out = tmp_path / "PF"
    argv = ["preflight", "--genome", str(f["genome"]), "--proteins", str(f["proteins"]),
            "--outdir", str(out), *LOW, *extra]  # fmt: skip
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output, out


def _tools(out: Path) -> list[dict[str, Any]]:
    data = json.loads((out / "preflight" / "preflight.json").read_text(encoding="utf-8"))
    return [c for c in data["checks"] if c["rule_id"] == "environment.tools"]


@pytest.fixture
def pf_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    db = write_setup_manifest(tmp_path / "db")
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    return db


def test_preflight_profile_reports_missing_images(tmp_path: Path, pf_db: Path) -> None:
    code, text, out = _preflight_cli(tmp_path, "--profile", "apptainer")
    assert "No such option" not in text, text
    assert code == 2, text
    checks = _tools(out)
    assert len(checks) == 1 and checks[0]["level"] == "error", checks
    assert f"{PKG} setup --images --profile apptainer" in checks[0]["message"]
    assert "bakta" in checks[0]["value"]["missing"]


def test_preflight_profile_finds_the_db_dir_cache(tmp_path: Path, pf_db: Path) -> None:
    cache = pf_db / "containers"
    cache.mkdir(parents=True, exist_ok=True)
    for ref in set(pins().values()):
        (cache / simple_name(ref)).write_bytes(b"sif")
    _code, text, out = _preflight_cli(tmp_path, "--profile", "apptainer")
    checks = _tools(out)
    assert len(checks) == 1, (checks, text)
    missing = set((checks[0].get("value") or {}).get("missing") or [])
    assert not set(pins()) & missing, missing


@pytest.mark.parametrize("db_dir", ["unset", "no-manifest"])
def test_preflight_profile_without_a_manifest_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, db_dir: str
) -> None:
    """With ``--profile`` a missing ``db_manifest.json`` is the run-mode error, missing tools
    included, not the plain-preflight info."""
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "root"))
    if db_dir == "unset":
        monkeypatch.delenv("GIPPA_DB_DIR", raising=False)
    else:
        (tmp_path / "db").mkdir()
        monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db"))
    code, text, out = _preflight_cli(tmp_path, "--profile", "apptainer")
    assert code == 2, text
    checks = _tools(out)
    assert len(checks) == 1 and checks[0]["level"] == "error", checks
    assert "bakta" in checks[0]["value"]["missing"]


@pytest.mark.parametrize("profile", ["apptainer", "docker", "farmshare,apptainer"])
def test_preflight_profile_matches_run_mode(
    tmp_path: Path, pf_db: Path, monkeypatch: pytest.MonkeyPatch, profile: str
) -> None:
    """The tool check of ``gippa preflight --profile P`` is the one ``gippa run --profile P``
    runs (run mode, same DB directory and ``$GIPPA_ROOT``)."""
    from gippa.preflight.run import PreflightInputs, run_preflight

    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "root"))
    _code, text, out = _preflight_cli(tmp_path, "--profile", profile)
    got = _tools(out)
    assert len(got) == 1, (got, text)
    f = _inputs(tmp_path)
    cfg = resolve_config(None, [LOW[1], LOW[3]])
    report, _ = run_preflight(PreflightInputs(genome=f["genome"], proteins=f["proteins"]), cfg,
                              None, "run", tmp_path / "RUN", pf_db,
                              gippa_root=tmp_path / "root", profile=profile)  # fmt: skip
    want = [c for c in report.checks if c.rule_id == "environment.tools"]
    assert len(want) == 1
    assert got[0]["level"] == want[0].level
    assert sorted((got[0].get("value") or {}).get("missing") or []) == sorted(
        (want[0].value or {}).get("missing") or []
    )


def test_preflight_profile_validation_is_runs(tmp_path: Path, pf_db: Path) -> None:
    code, text, _out = _preflight_cli(tmp_path, "--profile", "bogus")
    assert code == 2 and "unknown profile" in text, text
    assert "No such option" not in text


def test_preflight_profile_in_help() -> None:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["preflight", "--help"])
    assert "--profile" in res.output


def test_preflight_profile_none_runs_no_tool_check(tmp_path: Path, pf_db: Path) -> None:
    """An empty cache is not an error under --profile none (inputs only; a bare preflight is
    --profile auto since M11-15a)."""
    code, text, out = _preflight_cli(tmp_path, "--profile", "none")
    assert all(c["level"] != "error" for c in _tools(out)), text
    assert code == 0, text


def _needed_pinned() -> list[str]:
    cfg = resolve_config(None, [])
    return sorted(t.name for t in required_tools(cfg) if t.name in pins())


@settings(
    max_examples=8, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(present=st.sets(st.sampled_from(sorted(set(pins())))))
def test_property_preflight_profile_missing_images(
    tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch, present: set[str]
) -> None:
    """A pinned tool the default run needs is missing iff its image is not in the cache."""
    tmp = tmp_path_factory.mktemp("prop")
    db = write_setup_manifest(tmp / "db")
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    cache = tmp / "nxf"
    cache.mkdir()
    monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", str(cache))
    refs = {pins()[t] for t in present}
    for ref in refs:
        (cache / simple_name(ref)).write_bytes(b"sif")
    _code, text, out = _preflight_cli(tmp, "--profile", "apptainer")
    checks = _tools(out)
    assert len(checks) == 1, text
    missing = set((checks[0].get("value") or {}).get("missing") or [])
    for tool in _needed_pinned():
        assert (tool in missing) == (pins()[tool] not in refs), (tool, missing)


# ---------------------------------------------------------------- G4: no Nextflow on PATH


def test_run_without_nextflow_points_at_the_install_docs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import gippa.cli

    monkeypatch.setattr(gippa.cli, "_run_preflight", lambda cfg, *a, **k: ({}, cfg))
    bindir = tmp_path / "javaonly"
    link_prereqs(bindir, ("java",))  # Java 21, no nextflow anywhere (M11-15d)
    monkeypatch.setenv("PATH", str(bindir))
    db = write_setup_manifest(tmp_path / "db")
    code, text = _gippa_run(tmp_path, "docker", "--db-dir", str(db))
    assert code == 2, text
    assert "nextflow" in text.lower()
    assert "AGENTS.md" not in text, text
    assert "get-started/install" in text, text
    # M11-15d: the lookup is PATH, then <db_dir>/bin/nextflow, which `gippa setup` installs
    assert str(db / "bin" / "nextflow") in text, text
    assert f"{PKG} setup" in text, text
