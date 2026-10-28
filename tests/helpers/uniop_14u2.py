"""Fakes for the M11-14u2 tests (UniOP local build: apptainer, micromamba, prebuilt, terms).

docs/specs/M11-14-container-install.md item 1b ("Terms notice", "Build per runtime: apptainer
(14u2)", "micromamba (14u2)", "What is recorded (14u2)"), item 3 (the uniop rows), examples
16b-16d, 18-18c, 19-19d, 20, 21; AC 9 (14u2 list). Builds on ``helpers.uniop_14u1``.

- ``FAKE_APPTAINER``: a logging ``apptainer`` for ``PATH``. ``build [--fakeroot
  [--ignore-fakeroot-command]] <sif> <def>`` needs the staged context as its working directory
  (``%files`` sources are relative to it, as with the real apptainer): ``uniop-src.tar.gz`` and
  ``workflow/envs/uniop.yml`` in ``$PWD``. ``$FAKE_APPTAINER_FAIL`` = ``fakeroot`` makes a plain
  ``--fakeroot`` build fail with FarmShare's ``GLIBC_2.38 not found`` (the retry with
  ``--ignore-fakeroot-command`` succeeds); ``all`` makes every build fail with
  ``$FAKE_APPTAINER_FAIL_LINE`` last on stderr. A good build writes a few bytes to ``<sif>`` and
  copies ``<def>`` to ``$FAKE_APPTAINER_STATE/def_seen``. ``pull <img> docker://<ref>`` writes
  ``<img>`` (the 14f2 path); ``exec`` exits 0; ``push`` is only logged.
- ``FAKE_MICROMAMBA``: a logging ``micromamba``. ``create … -p|--prefix <prefix> … -f|--file
  <lock>`` makes ``<prefix>/conda-meta`` and ``<prefix>/bin`` and copies ``<lock>`` to
  ``$FAKE_MAMBA_STATE/lock_seen``; ``$FAKE_MAMBA_FAIL`` makes it print that line and exit 1.
  ``install``/``run``/``env``/``clean`` exit 0.
"""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path
from typing import Any

import pytest

from helpers.uniop_14u1 import COMMIT, PKG, REPO, install_fake_docker, uniop_key

NOTICE = (
    "UniOP (github.com/hongsua/UniOP @3ed14dfc) is fetched from its authors and built on this "
    f"machine; it has no licence file, so {PKG} does not distribute it and whether your use is "
    "permitted is yours to decide."
)
PROMPT = "Build UniOP for your own use? [y/N]"
TERMS_NOT_ACCEPTED = "terms not accepted; rerun with --accept-uniop-terms"
GLIBC = "/.singularity.d/libs/faked: /lib64/libc.so.6: version `GLIBC_2.38' not found"
APPTAINER_FAILED = (
    "apptainer build failed (fakeroot unavailable or build error): {line}; use --profile docker "
    "or conda, or set annotation.uniop.enabled=auto"
)
PREBUILT_HINT = "sbatch scripts/slurm/build_envs.sbatch uniop"
UNIOP_DEF = REPO / "workflow" / "containers" / "uniop" / "uniop.def"
UNIOP_LOCK = REPO / "workflow" / "envs" / "locks" / "uniop.explicit.txt"
DB_SCHEMA = REPO / "src" / "gippa" / "schemas" / "db_manifest.schema.json"
RUN_SCHEMA = REPO / "src" / "gippa" / "schemas" / "run_manifest.schema.json"

FAKE_APPTAINER = r"""#!/bin/sh
printf '%s\n' "apptainer $*" >> "$FAKE_APPTAINER_LOG"
state="$FAKE_APPTAINER_STATE"
mkdir -p "$state"
n=$#; i=0; prev=""; last=""
for a in "$@"; do i=$((i + 1)); prev=$last; last=$a; done
case "$1" in
build)
    fakeroot=""; ignore=""
    for a in "$@"; do
        case "$a" in --fakeroot) fakeroot=1;; --ignore-fakeroot-command) ignore=1;; esac
    done
    if [ "$FAKE_APPTAINER_FAIL" = all ]; then
        echo "INFO:    Starting build..." >&2
        echo "$FAKE_APPTAINER_FAIL_LINE" >&2
        exit 1
    fi
    if [ "$FAKE_APPTAINER_FAIL" = fakeroot ] && [ -n "$fakeroot" ] && [ -z "$ignore" ]; then
        echo "INFO:    Using fakeroot command combined with root mapping" >&2
        echo "$FAKE_GLIBC_LINE" >&2
        exit 1
    fi
    [ -f "$last" ] || { echo "fake apptainer: no definition file $last" >&2; exit 9; }
    [ -f "$PWD/uniop-src.tar.gz" ] || {
        echo "fake apptainer: %files uniop-src.tar.gz not in the working directory $PWD" >&2
        exit 9; }
    [ -f "$PWD/workflow/envs/uniop.yml" ] || {
        echo "fake apptainer: %files workflow/envs/uniop.yml not in $PWD" >&2; exit 9; }
    cp "$last" "$state/def_seen"
    printf 'fake sif built from %s\n' "$(cat "$PWD/uniop-src.tar.gz" | cksum)" > "$prev"
    echo "INFO:    Build complete: $prev"
    exit 0;;
pull)
    printf 'fake sif %s\n' "$last" > "$prev"
    exit 0;;
exec|run|inspect|version|--version) exit 0;;
push) exit 0;;
esac
echo "fake apptainer: unexpected call: $*" >&2
exit 3
"""

FAKE_MICROMAMBA = r"""#!/bin/sh
printf '%s\n' "micromamba $*" >> "$FAKE_MAMBA_LOG"
state="$FAKE_MAMBA_STATE"
mkdir -p "$state"
case "$1" in
create)
    if [ -n "$FAKE_MAMBA_FAIL" ]; then
        echo "error    libmamba Could not solve for environment specs" >&2
        echo "$FAKE_MAMBA_FAIL" >&2
        exit 1
    fi
    p=""; f=""; prev=""
    for a in "$@"; do
        case "$prev" in -p|--prefix) p=$a;; -f|--file) f=$a;; esac
        case "$a" in --prefix=*) p=${a#--prefix=};; --file=*) f=${a#--file=};; esac
        prev=$a
    done
    [ -n "$p" ] || { echo "fake micromamba: create without -p" >&2; exit 9; }
    [ -f "$f" ] || { echo "fake micromamba: create without a lock file" >&2; exit 9; }
    mkdir -p "$p/conda-meta" "$p/bin"
    cp "$f" "$state/lock_seen"
    exit 0;;
install|run|env|clean|info|--version) exit 0;;
esac
echo "fake micromamba: unexpected call: $*" >&2
exit 3
"""


def _exe(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def install_fake_runtimes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Fake ``docker`` (14u1), ``apptainer`` and ``micromamba`` first on ``PATH``; returns the
    docker state directory."""
    # the real apptainer is shadowed here, so the conftest FarmShare probe pin stays
    state = install_fake_docker(tmp_path, monkeypatch, docker_host=False)
    bindir = tmp_path / "fakebin"
    _exe(bindir / "apptainer", FAKE_APPTAINER)
    _exe(bindir / "micromamba", FAKE_MICROMAMBA)
    monkeypatch.setenv("FAKE_APPTAINER_LOG", str(tmp_path / "apptainer.log"))
    monkeypatch.setenv("FAKE_APPTAINER_STATE", str(tmp_path / "apptainer_state"))
    monkeypatch.setenv("FAKE_GLIBC_LINE", GLIBC)
    monkeypatch.setenv("FAKE_MAMBA_LOG", str(tmp_path / "micromamba.log"))
    monkeypatch.setenv("FAKE_MAMBA_STATE", str(tmp_path / "micromamba_state"))
    for var in ("FAKE_APPTAINER_FAIL", "FAKE_APPTAINER_FAIL_LINE", "FAKE_MAMBA_FAIL",
                "MAMBA_EXE", "CONDA_PREFIX", "CONDA_EXE"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    return state


def calls(tmp_path: Path, prog: str, sub: str = "") -> list[list[str]]:
    """The logged argv (without the program name) of ``prog``, optionally only ``sub`` calls."""
    log = tmp_path / f"{prog}.log"
    lines = log.read_text().splitlines() if log.exists() else []
    out = [ln.split()[1:] for ln in lines]
    return [c for c in out if not sub or c[:1] == [sub]]


def sif_path(cache: Path) -> Path:
    """``<cache>/localhost-<pkg>-uniop-<version>-<key>.sif`` (spec item 1b)."""
    return cache / f"localhost-{PKG}-uniop-git-{COMMIT[:8]}-{uniop_key()}.sif"


def prebuilt_root(root: Path, *, share: bool = True) -> Path:
    """A fake ``$GIPPA_ROOT`` with ``micromamba/envs/gippa-uniop`` (and its ``share/uniop``)."""
    prefix = root / "micromamba" / "envs" / "gippa-uniop"
    (prefix / "bin").mkdir(parents=True)
    if share:
        (prefix / "share" / "uniop").mkdir(parents=True)
        (prefix / "share" / "uniop" / "GIPPA_COMMIT").write_text(COMMIT + "\n")
    return root


def validate_db_manifest(db: Path) -> dict[str, Any]:
    """``db_manifest.json`` validated against ``db_manifest.schema.json``; returns it."""
    import jsonschema

    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    schema = json.loads(DB_SCHEMA.read_text(encoding="utf-8"))
    jsonschema.validate(data, schema)
    return dict(data)


def set_local_build(db: Path, **changes: Any) -> None:
    """Edit the ``local_build`` object of the recorded ``image:uniop`` item in place."""
    path = db / "db_manifest.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    for item in data["items"]:
        if item["name"] == "image:uniop":
            item["local_build"].update(changes)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def as_root(monkeypatch: pytest.MonkeyPatch) -> None:
    """``os.geteuid() == 0`` for the code under test (the fake still runs as the user)."""
    monkeypatch.setattr(os, "geteuid", lambda: 0)
