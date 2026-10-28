"""M11-14u1: the docker build of UniOP, ``local_tools.config``, the record, ``setup --images``.

docs/specs/M11-14-container-install.md item 1b ("Included by default", "Fetch" staging context,
"Build per runtime: docker (14u1)", "``local_tools.config``", "Exit codes", "What is recorded"),
item 3 (the table row "uniop with a local build"), examples 16, 18/18b (the no-build part only),
19a, AC 9 (14u1 list).

Interface pinned here (readings chosen where the spec leaves room; listed in the PR):

- docker build argv, exactly:
  ``docker build -f <context>/Dockerfile -t localhost/<pkg>-uniop:git-3ed14dfc-<key> <context>``
  with ``<context>`` = ``<cache>/uniop/context`` and ``<cache>`` = ``--cache-dir`` else
  ``image_cache(os.environ, db_dir)``. The image ID is read afterwards with
  ``docker image inspect`` (``--format`` with ``.Id``, or the JSON). Nothing is ever pushed.
- The staging context holds ``Dockerfile``, ``uniop.def`` (root or its repo path),
  ``uniop-src.tar.gz`` (the verified archive bytes) and every ``COPY`` source at its repo path.
- ``<db_dir>/local_tools.config`` is ``profiles { docker { process { withName: 'UNIOP' {
  container = '<tag>'; conda = '' } } } }``; it parses with 14f1's reader; no ``conda.enabled``.
- ``db_manifest.json`` item ``image:uniop``: ``version`` ``git:<commit>``, ``url`` the archive URL,
  ``checksum`` ``sha256:<archive>``, ``path`` ``""`` (docker), ``status``, and the object
  ``local_build`` = ``{"runtime", "ref", "digest", "key", "terms"}`` (``digest`` = the image ID
  as docker prints it, ``sha256:<hex>``). The item survives later manifest writes.
- A build failure: ``image:uniop`` not ``ok`` with the last stderr line in its message, exit 2,
  no ``docker`` block written; 14d's ``uniop_setup_record`` turns it into ``(last setup: …)``.
- Terms (the gate only; the notice and the prompt are 14u2): without ``--accept-uniop-terms`` and
  without a TTY nothing is fetched or built; ``--uniop`` then exits 2, ``--images`` exits 0.
  ``--accept-uniop-terms`` records ``terms: "flag"``.
- ``setup --images --profile docker`` runs the UniOP step when ``uniop`` is in
  ``required_tools(cfg)``; ``--no-uniop`` and precomputed operon pairs skip it.
- ``gippa run`` adds ``-c <db_dir>/local_tools.config`` only when the file exists and R(P) is
  docker, apptainer, or micromamba under plain ``conda``; never under ``farmshare*``.
"""

from __future__ import annotations

import itertools
import os
import random
import re
import shutil
import stat
import subprocess
from pathlib import Path
from typing import Any

import pytest
from helpers.uniop_14u1 import (
    ACCEPT,
    ARCHIVE_URL,
    DB_PATH,
    FAKE_IMAGE_ID,
    PKG,
    REPO,
    SKIP_DBS,
    VERSION,
    FakeNet,
    docker_calls,
    fixture_archive,
    install_fake_docker,
    install_fake_net,
    isolate,
    manifest_items,
    pin_uniop_source,
    setup_cli,
    sha256_hex,
    uniop_key,
    uniop_tag,
)

pytestmark = pytest.mark.usefixtures("no_network")

LOCAL_TOOLS = "local_tools.config"


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    isolate(tmp_path, monkeypatch)
    old = os.umask(0o022)
    yield
    os.umask(old)


@pytest.fixture
def state(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    return install_fake_docker(tmp_path, monkeypatch)


@pytest.fixture
def net(monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    return install_fake_net(monkeypatch)


@pytest.fixture
def archive(net: FakeNet, monkeypatch: pytest.MonkeyPatch) -> bytes:
    data = fixture_archive()
    pin_uniop_source(monkeypatch, sha256_hex(data))
    net.serve_github(data)
    return data


def uniop_args(db: Path, *extra: str, profile: str = "docker") -> list[str]:
    return ["--db-dir", str(db), *SKIP_DBS, "--uniop", "--profile", profile, *extra]


def built(tmp_path: Path, db: Path) -> None:
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 0, out


# ---------------------------------------------------------------- example 16: docker build


def test_example_16_docker_build_argv(state: Path, archive: bytes, tmp_path: Path) -> None:
    db = tmp_path / "db"
    built(tmp_path, db)
    builds = docker_calls(tmp_path, "build")
    assert len(builds) == 1, builds
    ctx = builds[0][-1]
    assert builds[0] == ["build", "-f", f"{ctx}/Dockerfile", "-t", uniop_tag(), ctx]
    assert Path(ctx).resolve() == (db / "containers" / "uniop" / "context").resolve()
    assert docker_calls(tmp_path, "push") == [], "nothing is ever pushed"


def test_tag_is_localhost_pkg_uniop_version_key(
    state: Path, archive: bytes, tmp_path: Path
) -> None:
    from gippa.io.image_key import tag

    built(tmp_path, tmp_path / "db")
    (b,) = docker_calls(tmp_path, "build")
    assert b[4] == f"localhost/{PKG}-uniop:{tag(VERSION, uniop_key())}"
    assert re.fullmatch(rf"localhost/{PKG}-uniop:git-3ed14dfc-[0-9a-f]{{12}}", b[4]), b[4]


def test_staging_context_holds_the_bundled_files_and_the_archive(
    state: Path, archive: bytes, tmp_path: Path
) -> None:
    built(tmp_path, tmp_path / "db")
    ctx = state / "context"  # the fake copied it at build time
    assert (ctx / "uniop-src.tar.gz").read_bytes() == archive
    dockerfile = ctx / "Dockerfile"
    assert dockerfile.read_bytes() == (REPO / "workflow/containers/uniop/Dockerfile").read_bytes()
    assert (ctx / "uniop.def").is_file() or (ctx / "workflow/containers/uniop/uniop.def").is_file()
    assert (ctx / "workflow/envs/uniop.yml").read_bytes() == (
        REPO / "workflow/envs/uniop.yml"
    ).read_bytes()
    from gippa.io.image_key import _copy_sources

    for src in _copy_sources(dockerfile.read_text(encoding="utf-8")):
        assert (ctx / src.removeprefix("./")).exists(), f"COPY source {src} not staged"


def test_no_builder_fetches_anything_itself(state: Path, archive: bytes, net: FakeNet,
                                            tmp_path: Path) -> None:  # fmt: skip
    built(tmp_path, tmp_path / "db")
    assert len(net.contacted) == 2  # the one GET and its redirect, nothing during the build
    assert {c[0] for c in docker_calls(tmp_path)} <= {"build", "image", "inspect", "version",
                                                      "info"}  # fmt: skip


# ---------------------------------------------------------------- local_tools.config


def test_local_tools_config_has_only_the_docker_block(
    state: Path, archive: bytes, tmp_path: Path
) -> None:
    from gippa.preflight.environment import _local_containers

    db = tmp_path / "db"
    built(tmp_path, db)
    path = db / LOCAL_TOOLS
    text = path.read_text(encoding="utf-8")
    assert _local_containers(path, "docker") == {"uniop": uniop_tag()}
    assert _local_containers(path, "apptainer") == {}
    assert "conda.enabled" not in text
    assert re.search(r"^\s*profiles\s*\{", text, re.MULTILINE), text
    assert re.search(r"\bdocker\s*\{", text) and not re.search(r"\b(apptainer|conda)\s*\{", text)
    for line in text.splitlines():
        if re.search(r"\bcontainer\s*=", line):
            assert re.search(r"\bconda\s*=\s*''", line), f"no conda = '' on: {line}"
    assert re.search(r"withName\s*:\s*'UNIOP'", text), text


def test_preflight_finds_the_local_image_under_docker_only(
    state: Path, archive: bytes, tmp_path: Path
) -> None:
    """14f1's row: present iff the local tag exists (``docker image inspect``, the fake)."""
    from gippa.annotate.tools import TOOLS
    from gippa.preflight.environment import missing_tool_reasons

    db = tmp_path / "db"
    built(tmp_path, db)
    lt, cache = db / LOCAL_TOOLS, db / "containers"
    assert missing_tool_reasons([TOOLS["uniop"]], "docker", cache, None, lt) == {}
    assert "uniop" in missing_tool_reasons([TOOLS["uniop"]], "apptainer", cache, None, lt)
    (state / "tags").write_text("")  # the local tag removed
    assert "uniop" in missing_tool_reasons([TOOLS["uniop"]], "docker", cache, None, lt)


def test_local_tools_config_keeps_other_runtime_blocks(
    state: Path, archive: bytes, tmp_path: Path
) -> None:
    """Each builder adds its block: an existing apptainer block survives a docker build."""
    from gippa.preflight.environment import _local_containers

    db = tmp_path / "db"
    db.mkdir()
    sif = str(tmp_path / "x.sif")
    (db / LOCAL_TOOLS).write_text(
        "profiles {\n"
        f"    apptainer {{ process {{ withName: 'UNIOP' {{ container = '{sif}'; "
        "conda = '' } } }\n"
        "}\n"
    )
    built(tmp_path, db)
    path = db / LOCAL_TOOLS
    assert _local_containers(path, "apptainer") == {"uniop": sif}
    assert _local_containers(path, "docker") == {"uniop": uniop_tag()}
    assert "conda.enabled" not in path.read_text()


@pytest.mark.needs_nextflow
def test_local_tools_config_parses_under_nextflow(
    state: Path, archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = tmp_path / "db"
    built(tmp_path, db)
    proj = tmp_path / "proj"
    proj.mkdir()
    (proj / "main.nf").write_text("workflow {}\n")
    (proj / "nextflow.config").write_text("")
    env = {**os.environ, "NXF_OFFLINE": "true", "NXF_HOME": str(tmp_path / "nxf_home")}
    nf = shutil.which("nextflow")
    assert nf is not None
    proc = subprocess.run(
        [nf, "-c", str(db / LOCAL_TOOLS), "config", "-profile", "docker", str(proj)],
        capture_output=True, text=True, env=env, timeout=300, check=False,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    assert uniop_tag() in proc.stdout


# ---------------------------------------------------------------- the record


def test_example_16_image_uniop_record(state: Path, archive: bytes, tmp_path: Path) -> None:
    db = tmp_path / "db"
    built(tmp_path, db)
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "ok", item
    assert item["version"] == VERSION
    assert item["url"] == ARCHIVE_URL
    assert item["checksum"] == f"sha256:{sha256_hex(archive)}"
    assert item["path"] == ""  # docker: the image lives in docker's store
    assert item["local_build"] == {
        "runtime": "docker",
        "ref": uniop_tag(),
        "digest": FAKE_IMAGE_ID,
        "key": uniop_key(),
        "terms": "flag",
    }


def test_the_record_survives_later_manifest_writes(
    state: Path, archive: bytes, tmp_path: Path
) -> None:
    from gippa.setupdb.manifest import read_records, write_manifest

    db = tmp_path / "db"
    built(tmp_path, db)
    assert "image:uniop" in read_records(db)
    write_manifest(db, [])  # any later setup merges into the manifest
    item = manifest_items(db)["image:uniop"]
    assert item["local_build"]["runtime"] == "docker"
    assert item["status"] == "ok"


def test_uniop_source_record_when_present(state: Path, archive: bytes, tmp_path: Path) -> None:
    db = tmp_path / "db"
    built(tmp_path, db)
    src = manifest_items(db).get("uniop_source")
    if src is not None:  # reading: recorded like any other item when the implementation does
        assert src["status"] == "ok" and src["path"] == DB_PATH
        assert src["checksum"] == f"sha256:{sha256_hex(archive)}"


# ---------------------------------------------------------------- build failure


def test_build_failure_is_recorded_and_exits_2(
    state: Path, archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_BUILD_FAIL", "ERROR: failed to solve: fixture boom")
    db = tmp_path / "db"
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 2, out
    item = manifest_items(db)["image:uniop"]
    assert item["status"] != "ok", item
    assert "ERROR: failed to solve: fixture boom" in item["message"], item
    lt = db / LOCAL_TOOLS
    if lt.exists():
        from gippa.preflight.environment import _local_containers

        assert _local_containers(lt, "docker") == {}, "no docker block for a failed build"
    from gippa.preflight import environment

    record = getattr(environment, "uniop_setup_record", None)  # 14d (#387)
    if record is not None:
        assert record(db) == f"(last setup: {item['message']})"


def test_docker_missing_from_path_is_recorded(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No fake docker and an empty PATH for it: a recorded failure, not a crash."""
    bindir = tmp_path / "onlysh"
    bindir.mkdir()
    for exe in ("sh", "tar", "gzip"):
        found = shutil.which(exe)
        if found:
            (bindir / exe).symlink_to(found)
    monkeypatch.setenv("PATH", str(bindir))
    db = tmp_path / "db"
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 2, out
    item = manifest_items(db)["image:uniop"]
    assert item["status"] != "ok" and "docker" in item["message"], item


def test_os_error_in_the_build_is_recorded_and_exits_2(
    state: Path, archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An OSError (e.g. a shared cache dir owned by another user) is a record, not a traceback."""
    from gippa.setupdb import uniop_build

    def denied(context: Path, archive: Path, root: Path) -> None:
        raise PermissionError(13, "Permission denied", str(context))

    monkeypatch.setattr(uniop_build, "stage_context", denied)
    db = tmp_path / "db"
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 2, out
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "absent" and "Permission denied" in item["message"], item


# ---------------------------------------------------------------- terms gate


def test_uniop_without_acceptance_and_no_tty_builds_nothing(
    state: Path, archive: bytes, net: FakeNet, tmp_path: Path
) -> None:
    db = tmp_path / "db"
    code, out = setup_cli(uniop_args(db))
    assert code == 2, out
    assert net.contacted == [], "nothing is fetched before the terms are accepted"
    assert docker_calls(tmp_path, "build") == []
    assert not (db / DB_PATH).exists()
    assert "--accept-uniop-terms" in out


def test_images_without_acceptance_builds_nothing_and_exits_0(
    state: Path, archive: bytes, net: FakeNet, tmp_path: Path
) -> None:
    db = tmp_path / "db"
    code, out = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "docker"])
    assert code == 0, out
    assert net.contacted == []
    assert docker_calls(tmp_path, "build") == []
    assert docker_calls(tmp_path, "pull"), "the other images are still pulled"
    assert manifest_items(db)["image:uniop"]["status"] == "absent"


# ---------------------------------------------------------------- included in setup --images


def test_images_includes_the_uniop_build_by_default(
    state: Path, archive: bytes, tmp_path: Path
) -> None:
    db = tmp_path / "db"
    code, out = setup_cli(
        ["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "docker", ACCEPT]
    )
    assert code == 0, out
    assert len(docker_calls(tmp_path, "build")) == 1
    items = manifest_items(db)
    assert items["image:uniop"]["status"] == "ok"
    assert items["image:bakta"]["status"] == "ok", "the pinned images are pulled as before"
    assert docker_calls(tmp_path, "push") == []


@pytest.mark.parametrize(
    "extra",
    [["--no-uniop"], ["--set", 'annotation.precomputed.operon_pairs="/x/operon_pairs.tsv"']],
    ids=["no-uniop", "precomputed-pairs"],
)
def test_images_skips_uniop_when_not_needed_or_disabled(
    state: Path, archive: bytes, net: FakeNet, tmp_path: Path, extra: list[str]
) -> None:
    db = tmp_path / "db"
    code, out = setup_cli(
        ["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "docker", ACCEPT, *extra]
    )
    assert code == 0, out
    assert net.contacted == []
    assert docker_calls(tmp_path, "build") == []
    assert "image:uniop" not in manifest_items(db)


def test_plain_setup_never_builds_uniop(
    state: Path, archive: bytes, net: FakeNet, tmp_path: Path
) -> None:
    db = tmp_path / "db"
    # M11-15c2 pulls images by default (and then builds UniOP after consent, spec example 6.1):
    # the pre-15c2 "plain" setup, without images, is --no-images now
    code, out = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--no-images", ACCEPT])
    assert code == 0, out
    assert net.contacted == [] and docker_calls(tmp_path) == []
    assert not (db / LOCAL_TOOLS).exists()


def test_offline_under_images_installs_the_others_then_exits_2(
    state: Path, net: FakeNet, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from helpers.uniop_14u1 import ARCHIVE_PATH, Route

    pin_uniop_source(monkeypatch, sha256_hex(fixture_archive()))
    net.routes[("github.com", ARCHIVE_PATH)] = Route(error="[Errno 101] Network is unreachable")
    db = tmp_path / "db"
    code, out = setup_cli(
        ["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "docker", ACCEPT]
    )
    assert code == 2, out
    items = manifest_items(db)
    assert items["image:uniop"]["status"] == "absent"
    assert items["image:bakta"]["status"] == "ok"


# ---------------------------------------------------------------- example 19a: gippa run -c


FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""


@pytest.fixture
def run_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """A genome, proteins, a DB dir with ``local_tools.config``, a fake ``nextflow``; preflight
    stubbed (it is not under test here)."""
    from helpers.fake_env import fake_gippa_root, write_setup_manifest

    import gippa.cli

    rng = random.Random(5)
    genome = tmp_path / "in" / "genome.fa"
    genome.parent.mkdir()
    seq = "".join(rng.choice("ACGT") for _ in range(3000))
    genome.write_text(">chr1 [topology=circular]\n" + seq + "\n")
    proteins = tmp_path / "in" / "proteins.faa"
    proteins.write_text(">p1\nM" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(80)))
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    (db / LOCAL_TOOLS).write_text(
        "profiles {\n    docker { process { withName: 'UNIOP' { container = "
        "'localhost/x-uniop:t'; conda = '' } } }\n}\n"
    )
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
    bindir = tmp_path / "nfbin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))
    monkeypatch.setattr(gippa.cli, "_run_preflight", lambda cfg, *a, **k: ({}, cfg))
    return {"genome": genome, "proteins": proteins, "db": db, "args": tmp_path / "nf_args.txt"}


def _nf_argv(env: dict[str, Path], profile: str, tmp_path: Path) -> list[str]:
    from click.testing import CliRunner

    from gippa.cli import cli

    CliRunner().invoke(
        cli,
        ["run", "--genome", str(env["genome"]), "--proteins", str(env["proteins"]),
         "--outdir", str(tmp_path / "OUT"), "--db-dir", str(env["db"]), "--profile", profile],
    )  # fmt: skip
    assert env["args"].is_file(), "gippa run did not launch nextflow"
    return env["args"].read_text().splitlines()


def _has_local_tools(argv: list[str], db: Path) -> bool:
    want = str(db / LOCAL_TOOLS)
    pairs = itertools.pairwise(argv)
    return any(a in ("-c", "-config") and Path(b) == Path(want) for a, b in pairs)


@pytest.mark.parametrize(
    ("profile", "adds"),
    [
        ("docker", True),
        ("apptainer", True),
        ("conda", True),
        ("farmshare,apptainer", False),  # example 19a
        ("farmshare_local", False),
        ("farmshare", False),
        ("conda,docker", False),  # ambiguous
    ],
)
def test_run_adds_local_tools_config_only_for_its_runtimes(
    run_env: dict[str, Path], tmp_path: Path, profile: str, adds: bool
) -> None:
    argv = _nf_argv(run_env, profile, tmp_path)
    assert _has_local_tools(argv, run_env["db"]) is adds, argv


def test_run_passes_an_absolute_local_tools_path(
    run_env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Nextflow runs in outdir, so a relative --db-dir must reach -c as an absolute path."""
    from click.testing import CliRunner

    from gippa.cli import cli

    monkeypatch.chdir(tmp_path)
    env = run_env
    CliRunner().invoke(
        cli,
        ["run", "--genome", str(env["genome"]), "--proteins", str(env["proteins"]),
         "--outdir", "OUT", "--db-dir", "db", "--profile", "docker"],
    )  # fmt: skip
    argv = env["args"].read_text().splitlines()
    given = [b for a, b in itertools.pairwise(argv) if a in ("-c", "-config")]
    assert str(env["db"].absolute() / LOCAL_TOOLS) in given, argv


def test_run_without_the_file_adds_nothing(run_env: dict[str, Path], tmp_path: Path) -> None:
    (run_env["db"] / LOCAL_TOOLS).unlink()
    argv = _nf_argv(run_env, "docker", tmp_path)
    assert not any(a in ("-c", "-config") for a in argv), argv
