"""M11-14u2 example 21: a stub run after a docker UniOP build records the local build.

docs/specs/M11-14-container-install.md item 1b ("What is recorded (14u2)": the run manifest's
``uniop`` entry), examples 19 and 21, AC 9a.

``giap run --profile docker`` goes through the CLI, preflight and the run tail. The state "after
example 16" is written directly (no fetch, no build): ``db_manifest.json`` with an ``ok``
``image:uniop`` whose ``local_build`` is the docker build of the current ``key("uniop")``, and
``local_tools.config`` naming its tag. A ``nextflow`` wrapper forwards the runner's command line
to the real Nextflow with the profile swapped for ``test,docker,stub`` and ``-stub-run`` added
(the ``-c <db_dir>/local_tools.config`` that ``giap run`` adds is kept). A logging fake
``docker`` answers ``image inspect`` for the tag with the fake image ID, exits 0 for the UniOP
``--help`` probe, and runs each task's command on the host. Preflight's tool rows other than
``uniop`` are dropped (the conda-only promoter tools are missing under docker by design, item 3;
they are not under test here). No network, no real docker, never the real ``$GIAP_DB_DIR``.
"""

from __future__ import annotations

import itertools
import json
import os
import shutil
import stat
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_giap_root, setup_manifest
from helpers.uniop_14u1 import ARCHIVE_URL, FAKE_IMAGE_ID, VERSION, uniop_key, uniop_tag

from giap.exitcodes import ExitCode

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "workflow" / "tests" / "data" / "annotate"
ARCHIVE_SHA = "ab" * 32

WRAPPER = """#!/bin/sh
set -e
n=$#; i=0; prev=""
for a in "$@"; do
  i=$((i + 1))
  if [ "$prev" = "-profile" ]; then a="test,docker,stub"; fi
  set -- "$@" "$a"
  prev="$a"
done
shift "$n"
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exec "{real}" "$@" -stub-run
"""

FAKE_DOCKER = """#!/usr/bin/env bash
# logging fake docker (M11-14u2 example 21)
printf '%s ' "$@" >> "$FAKE_DOCKER_LOG"; printf '\\n' >> "$FAKE_DOCKER_LOG"
last=""; for a in "$@"; do last=$a; done
case "$1" in
image|inspect)
  [ "$last" = "$FAKE_TAG" ] || { echo "Error: No such image: $last" >&2; exit 1; }
  case "$*" in
  *"{{"*) echo "$FAKE_IMAGE_ID";;
  *) printf '[{"Id": "%s", "RepoTags": ["%s"]}]\\n' "$FAKE_IMAGE_ID" "$last";;
  esac
  exit 0;;
run)
  shift
  while [ $# -gt 0 ] && [ "$1" != /bin/bash ] && [ "$1" != /bin/sh ]; do shift; done
  [ $# -gt 0 ] || exit 0   # the UniOP --help probe: `docker run --rm <tag> sh -c ...`
  exec "$@";;
push) echo "fake docker: push is never allowed" >&2; exit 7;;
esac
exit 0
"""


def _exe(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _db(db: Path) -> Path:
    """The state after example 16 (docker): the manifest item and ``local_tools.config``."""
    db.mkdir()
    (db / "version.json").write_text('{"type": "light"}\n')
    data = setup_manifest()
    items = [i for i in data["items"] if i["name"] != "image:uniop"]  # type: ignore[attr-defined]
    items.append({"name": "image:uniop", "version": VERSION, "url": ARCHIVE_URL,
                  "checksum": f"sha256:{ARCHIVE_SHA}", "size_bytes": 0,
                  "retrieved": "2026-10-26", "path": "", "status": "ok", "message": "",
                  "local_build": {"runtime": "docker", "ref": uniop_tag(),
                                  "digest": FAKE_IMAGE_ID, "key": uniop_key(),
                                  "terms": "flag"}})  # fmt: skip
    data["items"] = items
    (db / "db_manifest.json").write_text(json.dumps(data, indent=2) + "\n")
    (db / "local_tools.config").write_text(
        "profiles {\n    docker { process { withName: 'UNIOP' { container = "
        f"'{uniop_tag()}'; conda = '' }} }} }}\n}}\n"
    )
    return db


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    import socket

    import giap.preflight.run as pre

    tmp = tmp_path_factory.mktemp("uniop_stub")
    mp = pytest.MonkeyPatch()
    try:

        def deny(*args: object, **_kw: object) -> None:
            raise OSError(f"no_network: network access is blocked in this test ({args!r})")

        for name in ("create_connection", "getaddrinfo"):
            mp.setattr(socket, name, deny)
        mp.setattr(socket.socket, "connect", deny)
        real_nf = shutil.which("nextflow")
        assert real_nf is not None
        bindir = tmp / "bin"
        bindir.mkdir()
        _exe(bindir / "nextflow", WRAPPER.replace("{real}", real_nf))
        _exe(bindir / "docker", FAKE_DOCKER)
        for var in ("GIAP_DB_DIR", "NXF_APPTAINER_CACHEDIR", "GIAP_CONTAINER_RUNTIME",
                    "SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE",
                    "SLURM_MEM_PER_CPU", "UNIOP_HOME"):  # fmt: skip
            mp.delenv(var, raising=False)
        mp.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
        mp.setenv("NXF_ANSI_LOG", "false")
        mp.setenv("NXF_OFFLINE", "true")
        mp.setenv("FAKE_NF_ARGS", str(tmp / "nf_args.txt"))
        mp.setenv("FAKE_DOCKER_LOG", str(tmp / "docker.log"))
        mp.setenv("FAKE_TAG", uniop_tag())
        mp.setenv("FAKE_IMAGE_ID", FAKE_IMAGE_ID)
        mp.setenv("GIAP_ROOT", str(fake_giap_root(tmp / "giap_root")))
        real_reasons = pre.missing_tool_reasons

        def uniop_only(*a: Any, **k: Any) -> dict[str, str]:
            return {k_: v for k_, v in real_reasons(*a, **k).items() if k_ == "uniop"}

        mp.setattr(pre, "missing_tool_reasons", uniop_only)
        db = _db(tmp / "db")
        proteins = tmp / "proteins.faa"
        proteins.write_text(">prot1\nM" + "ACDEFGHIKLMNPQRSTVWY" * 8 + "\n")
        out = tmp / "OUT"
        from giap.cli import cli

        argv = ["run", "--genome", str(DATA / "genome.fa"), "--proteins", str(proteins),
                "--outdir", str(out), "--db-dir", str(db), "--modules", "core",
                "--profile", "docker"]  # fmt: skip
        res = CliRunner().invoke(cli, argv)
        if res.exception is not None and not isinstance(res.exception, SystemExit):
            raise res.exception
        log = res.output
    finally:
        mp.undo()
    return {"code": res.exit_code, "log": log, "out": out, "db": db,
            "nf_args": tmp / "nf_args.txt", "docker_log": tmp / "docker.log"}  # fmt: skip


def _manifest(run: dict[str, Any]) -> dict[str, Any]:
    assert run["code"] == ExitCode.SUCCESS, run["log"]
    path = run["out"] / "run_manifest.json"
    assert path.is_file(), run["log"]
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def test_example_19_run_adds_the_local_tools_config(run: dict[str, Any]) -> None:
    assert run["code"] == ExitCode.SUCCESS, run["log"]
    args = run["nf_args"].read_text().splitlines()
    want = str(run["db"] / "local_tools.config")
    assert any(a in ("-c", "-config") and b == want for a, b in itertools.pairwise(args)), args


def test_example_21_operon_source_is_uniop(run: dict[str, Any]) -> None:
    assert _manifest(run)["operon_source"] == "uniop"


def test_example_21_uniop_entry_records_the_local_build(run: dict[str, Any]) -> None:
    m = _manifest(run)
    entries = [t for t in m["tools"] if t["name"] == "uniop"]
    assert len(entries) == 1, m["tools"]
    (entry,) = entries
    assert entry["local_build"] is True, entry
    assert entry["runtime"] == "docker"
    assert entry["source_sha256"] in (ARCHIVE_SHA, f"sha256:{ARCHIVE_SHA}"), entry
    assert entry["container_digest"] == FAKE_IMAGE_ID, entry
    assert "version" in entry


def test_example_21_uniop_ran_in_the_local_image(run: dict[str, Any]) -> None:
    assert run["code"] == ExitCode.SUCCESS, run["log"]
    lines = run["docker_log"].read_text().splitlines()
    assert any(ln.startswith("run ") and uniop_tag() in ln.split() for ln in lines), lines
    assert not any(ln.startswith("push") for ln in lines)


def test_example_21_manifest_validates(run: dict[str, Any]) -> None:
    from helpers.manifest_data import validate

    validate(_manifest(run))
