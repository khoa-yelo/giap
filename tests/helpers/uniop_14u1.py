"""Fakes for the M11-14u1 tests (UniOP local build: fetch and docker).

docs/specs/M11-14-container-install.md item 1b, examples 16, 16a, 16e, 17, 19a, AC 9.

- ``FakeNet``: a fake ``http.client`` for the real ``HttpsFetcher``. Routes are
  ``(host, path) -> Route``; a route may answer without ``Content-Length`` (codeload streams
  that way) or raise ``OSError`` (offline). No socket is ever opened.
- ``fixture_archive()``: a few-byte ``UniOP-<commit>/`` tarball made here, never UniOP code (R5).
- ``pin_uniop_source(monkeypatch, sha_hex)``: the ``uniop_source`` row with its checksum replaced,
  patched into every loaded ``gippa`` module that holds the original ``SOURCES`` mapping (the pin
  of the real archive is computed on FarmShare by the implementer, so the fixture cannot match
  it).
- ``FAKE_DOCKER``: a logging ``docker`` for ``PATH``. ``build`` copies its context to
  ``$FAKE_DOCKER_STATE/context`` and records the tag; ``image inspect``/``inspect`` answer for a
  recorded tag with ``$FAKE_IMAGE_ID`` (the ``.Id`` of ``--format`` or the JSON form);
  ``pull`` prints ``Digest: <pin>``; ``push`` is only logged. ``$FAKE_BUILD_FAIL`` makes
  ``build`` print that line last on stderr and exit 1.
"""

from __future__ import annotations

import contextlib
import dataclasses
import gzip
import hashlib
import http.client
import io
import os
import stat
import sys
import tarfile
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
COMMIT = "3ed14dfc2a5d69b993d99a26e50c1f8e9723c087"
VERSION = f"git:{COMMIT}"
ARCHIVE_URL = f"https://github.com/hongsua/UniOP/archive/{COMMIT}.tar.gz"
ARCHIVE_PATH = f"/hongsua/UniOP/archive/{COMMIT}.tar.gz"
CODELOAD_PATH = f"/hongsua/UniOP/tar.gz/{COMMIT}"
CODELOAD_URL = f"https://codeload.github.com{CODELOAD_PATH}"
DB_PATH = f"uniop/UniOP-{COMMIT}.tar.gz"
FAKE_IMAGE_ID = "sha256:" + "5e" * 32
SKIP_DBS = ["--bakta", "skip", "--no-genomad", "--no-promotech-models"]
ACCEPT = "--accept-uniop-terms"


# ---------------------------------------------------------------- fixture archive


def fixture_archive(payload: bytes = b"not UniOP code\n") -> bytes:
    """A deterministic ``.tar.gz`` with one small file under ``UniOP-<commit>/``."""
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode="w") as tar:
        for name, data in (
            (f"UniOP-{COMMIT}/README.fixture", payload),
            (f"UniOP-{COMMIT}/src/UniOP", b"#!/usr/bin/env python\nprint('fixture')\n"),
        ):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 0
            info.mode = 0o644
            tar.addfile(info, io.BytesIO(data))
    return gzip.compress(raw.getvalue(), mtime=0)


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------- the pinned row


def uniop_source() -> Any:
    from gippa.setupdb import items

    try:
        return items.SOURCES["uniop_source"]
    except KeyError:
        pytest.fail("db_sources.yaml has no uniop_source row (14u1)")


def pin_uniop_source(monkeypatch: pytest.MonkeyPatch, sha_hex: str) -> Any:
    """Replace the pin of ``uniop_source`` by ``sha256:<sha_hex>`` wherever ``SOURCES`` is held.

    Without the row (before 14u1) nothing is patched and ``None`` is returned, so the test fails
    on the missing feature itself rather than in a fixture."""
    from gippa.setupdb import items

    original = items.SOURCES
    row = original.get("uniop_source")
    if row is None:
        return None
    patched = MappingProxyType(
        {**original, "uniop_source": dataclasses.replace(row, checksum=f"sha256:{sha_hex}")}
    )
    for name, mod in list(sys.modules.items()):
        if (name == "gippa" or name.startswith("gippa.")) and getattr(mod, "SOURCES", None) is (
            original
        ):
            monkeypatch.setattr(mod, "SOURCES", patched)
    return patched["uniop_source"]


# ---------------------------------------------------------------- fake http.client


@dataclass
class Route:
    status: int = 200
    headers: dict[str, str] = field(default_factory=dict)
    body: bytes = b""
    content_length: bool = True
    error: str = ""  # raise OSError(error) instead of answering (offline)


class _Sock:
    def __init__(self, raw: bytes) -> None:
        self.raw = raw

    def makefile(self, *_a: Any, **_k: Any) -> io.BytesIO:
        return io.BytesIO(self.raw)


class FakeNet:
    """Fake ``http.client.HTTPSConnection``; ``contacted`` lists ``(host, path)`` in order."""

    def __init__(self) -> None:
        self.routes: dict[tuple[str, str], Route] = {}
        self.contacted: list[tuple[str, str]] = []

    def conn(self) -> type:
        net = self

        class Conn(http.client.HTTPConnection):
            def __init__(self, host: str, port: int | None = None, **_kw: Any) -> None:
                super().__init__(host, port)
                self.target = ""

            def connect(self) -> None:
                raise AssertionError("FakeNet never opens a socket")

            def close(self) -> None: ...

            def request(self, method: str, url: str, *_a: Any, **_k: Any) -> None:  # type: ignore[override]
                self.target = url

            def getresponse(self) -> http.client.HTTPResponse:
                key = (str(self.host), self.target.split("?")[0])
                net.contacted.append(key)
                route = net.routes.get(key, Route(404, {}, b"not found"))
                if route.error:
                    raise OSError(route.error)
                head = "".join(f"{k}: {v}\r\n" for k, v in route.headers.items())
                if route.content_length:
                    head += f"Content-Length: {len(route.body)}\r\n"
                raw = f"HTTP/1.1 {route.status} X\r\n{head}Connection: close\r\n\r\n".encode()
                resp = http.client.HTTPResponse(_Sock(raw + route.body), method="GET")  # type: ignore[arg-type]
                resp.begin()
                return resp

        return Conn

    def serve_github(self, body: bytes) -> None:
        """GitHub's ``302 -> codeload.github.com``; codeload streams without Content-Length."""
        self.routes[("github.com", ARCHIVE_PATH)] = Route(302, {"Location": CODELOAD_URL})
        self.routes[("codeload.github.com", CODELOAD_PATH)] = Route(
            200, {"Content-Type": "application/x-gzip"}, body, content_length=False
        )


def install_fake_net(monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    for var in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY", "all_proxy",
                "ALL_PROXY", "no_proxy", "NO_PROXY"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    net = FakeNet()
    monkeypatch.setattr(http.client, "HTTPSConnection", net.conn())
    monkeypatch.setattr(http.client, "HTTPConnection", net.conn())
    return net


# ---------------------------------------------------------------- fake docker

FAKE_DOCKER = r"""#!/bin/sh
printf '%s\n' "docker $*" >> "$FAKE_DOCKER_LOG"
state="$FAKE_DOCKER_STATE"
mkdir -p "$state"
last=""; for a in "$@"; do last=$a; done
case "$1" in
build)
    f=""; t=""; prev=""
    for a in "$@"; do
        case "$prev" in -f) f=$a;; -t) t=$a;; esac
        prev=$a
    done
    if [ -n "$FAKE_BUILD_FAIL" ]; then
        echo "#7 [3/5] RUN micromamba install" >&2
        echo "$FAKE_BUILD_FAIL" >&2
        exit 1
    fi
    [ -f "$f" ] || { echo "fake docker: no Dockerfile $f" >&2; exit 9; }
    [ -f "$last/uniop-src.tar.gz" ] || { echo "fake docker: no archive in $last" >&2; exit 9; }
    rm -rf "$state/context"; cp -R "$last" "$state/context"
    echo "$t" >> "$state/tags"
    echo "Successfully tagged $t"
    exit 0;;
image|inspect)
    if [ "$1" = image ] && [ "$2" != inspect ]; then
        echo "fake docker: unexpected call: $*" >&2; exit 3
    fi
    if ! grep -qxF "$last" "$state/tags" 2>/dev/null; then
        echo "Error: No such image: $last" >&2; exit 1
    fi
    case "$*" in
    *"{{"*) echo "$FAKE_IMAGE_ID";;
    *) printf '[{"Id": "%s", "RepoTags": ["%s"]}]\n' "$FAKE_IMAGE_ID" "$last";;
    esac
    exit 0;;
pull)
    echo "${last%@*}: Pulling from fake"
    echo "Digest: ${last#*@}"
    echo "Status: Downloaded newer image for $last"
    exit 0;;
push) exit 0;;
version|info) echo "fake docker 0.0"; exit 0;;
esac
echo "fake docker: unexpected call: $*" >&2
exit 3
"""


def install_fake_docker(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """``docker`` first on ``PATH``; returns the state directory."""
    bindir = tmp_path / "fakebin"
    bindir.mkdir(parents=True, exist_ok=True)
    exe = bindir / "docker"
    exe.write_text(FAKE_DOCKER)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_DOCKER_LOG", str(tmp_path / "docker.log"))
    monkeypatch.setenv("FAKE_DOCKER_STATE", str(tmp_path / "docker_state"))
    monkeypatch.setenv("FAKE_IMAGE_ID", FAKE_IMAGE_ID)
    monkeypatch.delenv("FAKE_BUILD_FAIL", raising=False)
    return tmp_path / "docker_state"


def docker_calls(tmp_path: Path, sub: str = "") -> list[list[str]]:
    log = tmp_path / "docker.log"
    lines = log.read_text().splitlines() if log.exists() else []
    return [ln.split()[1:] for ln in lines if not sub or ln.split()[1:2] == [sub]]


def isolate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Never the real ``$GIPPA_DB_DIR``/``$GIPPA_ROOT``/caches/envs; a fixed umask of 022."""
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "GIPPA_REPO", "NXF_APPTAINER_CACHEDIR",
                "GIPPA_CONTAINER_RUNTIME", "UNIOP_HOME"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(tmp_path / "no_mamba"))


def setup_cli(args: list[str]) -> tuple[int, str]:
    """``gippa setup <args>`` (no TTY); returns the exit code and stdout + stderr."""
    from click.testing import CliRunner

    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["setup", *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    text = res.output
    with contextlib.suppress(ValueError):  # stderr not captured separately
        text += res.stderr
    return res.exit_code, text


def manifest_items(db: Path) -> dict[str, dict[str, Any]]:
    import json

    path = db / "db_manifest.json"
    if not path.is_file():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    return {item["name"]: item for item in data["items"]}


def uniop_key() -> str:
    from gippa.io.image_key import key

    return key("uniop")


def uniop_tag() -> str:
    return f"localhost/{PKG}-uniop:git-{COMMIT[:8]}-{uniop_key()}"


def mode(path: Path) -> int:
    return stat.S_IMODE(path.stat().st_mode)
