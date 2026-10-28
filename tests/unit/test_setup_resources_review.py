"""M2-09 lead decisions after code review (PR #67): items 1, 4, 5 and 6.

docs/specs/M2-09-gippa-setup-resources.md, "Lead decisions after code review". Every test runs
under ``setup_guard``: no network, no tool, no real $GIPPA_DB_DIR.

Ambiguities resolved here:
- (1) ``--check`` with antiSMASH: the item passes when ``db_manifest.json`` holds an ``ok`` record
  with ``url: "image"`` and ``antismash --check-prereqs`` (through the command runner) exits 0.
  Without such a record (and without an installed directory) ``--check`` still fails: it verifies
  what setup recorded. ``--check`` never downloads.
- (4) The cap on a body without Content-Length (close-delimited or chunked) is
  ``item.size_gib * gippa.setupdb.run.GIB`` bytes, or ``1 * GIB`` when ``size_gib`` is 0; a direct
  ``HttpsFetcher().fetch(url, dest)`` for a URL that is no item uses 1 GiB. Tests shrink
  ``gippa.setupdb.run.GIB`` so the cap is a few kB. Exceeding it is a ``SourceError`` naming the
  host, and the partial file is removed (install leaves the item ``absent``).
- (5) An unsafe ``path`` in ``db_sources.yaml`` makes ``_load_sources`` raise ``ValueError`` that
  names the offending path.
- (6) ``db_manifest.json`` is mode 0644 whatever the umask.
"""

from __future__ import annotations

import dataclasses
import http.client
import io
import os
import shutil
import stat
import sys
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any, NamedTuple

import pytest
from setup_guard import BLOCKED, real_https_fetch, setup_guard  # noqa: F401
from test_setupdb import NoFetcher, NoRunner

REPO = Path(__file__).resolve().parents[2]
REAL_GIB = 1024**3
SMALL_GIB = 4096
KO_URL = "https://rest.kegg.jp/link/pathway/ko"
KO_PATH = "/link/pathway/ko"


class _Usage(NamedTuple):
    total: int
    used: int
    free: int


@pytest.fixture(autouse=True)
def _env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Plenty of disk, inside Slurm, every GIPPA root (and $GIPPA_DB_DIR) in ``tmp_path``."""
    free = 10_000 * REAL_GIB
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(free, 0, free))
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    roots = {"GIPPA_ROOT": "root", "MAMBA_ROOT_PREFIX": "mamba",
             "NXF_APPTAINER_CACHEDIR": "cache", "GIPPA_DB_DIR": "not-the-real-db"}  # fmt: skip
    for var, sub in roots.items():
        monkeypatch.setenv(var, str(tmp_path / "_env" / sub))
    monkeypatch.setenv("GIPPA_REPO", str(REPO))


def sources() -> Any:
    from gippa.setupdb import SOURCES

    return SOURCES


def repin(monkeypatch: pytest.MonkeyPatch, name: str, **fields: Any) -> Any:
    """Replace ``SOURCES[name]`` in every ``gippa.setupdb`` module that holds it."""
    item = dataclasses.replace(sources()[name], **fields)
    patched = {**sources(), name: item}
    for mod_name, mod in list(sys.modules.items()):
        if mod_name.startswith("gippa.setupdb") and hasattr(mod, "SOURCES"):
            monkeypatch.setattr(mod, "SOURCES", patched)
    return item


def opts(db: Path, **kw: Any) -> Any:
    from gippa.setupdb import SetupOptions

    base: dict[str, Any] = {"db_dir": db, "bakta": "skip", "genomad": False, "antismash": False,
                            "promotech_models": False, "check": False}  # fmt: skip
    return SetupOptions(**(base | kw))


def setup(db: Path, fetcher: Any = None, runner: Any = None, **kw: Any) -> int:
    from gippa.setupdb import run_setup

    return int(run_setup(opts(db, **kw), runner or NoRunner(), fetcher or NoFetcher()))


def files_under(d: Path) -> list[Path]:
    return sorted(p for p in d.rglob("*") if p.is_file()) if d.exists() else []


# ------------------------------------------------- 1. antiSMASH --check from the image


class PrereqRunner:
    """``antismash --check-prereqs`` exits ``prereqs``; any download is a test failure."""

    def __init__(self, prereqs: int = 0, unavailable: bool = False) -> None:
        self.prereqs = prereqs
        self.unavailable = unavailable
        self.calls: list[list[str]] = []

    def run(self, tool: str, argv: Sequence[str]) -> int:
        from gippa.setupdb import ToolUnavailable

        self.calls.append(list(argv))
        assert tool == "antismash", tool
        if self.unavailable:
            raise ToolUnavailable("no gippa-antismash environment and no pinned antismash image")
        if "--check-prereqs" in argv:
            return self.prereqs
        raise AssertionError(f"no download expected: {list(argv)}")

    def prereq_calls(self) -> list[list[str]]:
        return [c for c in self.calls if "--check-prereqs" in c]


def install_from_image(db: Path) -> None:
    assert setup(db, runner=PrereqRunner(0), antismash=True) == 0
    from gippa.setupdb import read_records

    rec = read_records(db)["antismash"]
    assert (rec.status, rec.url) == ("ok", "image"), rec


def test_check_accepts_antismash_recorded_as_image(tmp_path: Path) -> None:
    """Reviewer's repro: install records ``url: image``; ``--check`` then exited 2."""
    install_from_image(tmp_path)
    runner = PrereqRunner(0)
    assert setup(tmp_path, runner=runner, antismash=True, check=True) == 0
    assert runner.prereq_calls(), "--check must ask antismash --check-prereqs in the image"
    assert BLOCKED == []


def test_check_image_record_fails_when_prereqs_fail(tmp_path: Path) -> None:
    install_from_image(tmp_path)
    runner = PrereqRunner(1)
    assert setup(tmp_path, runner=runner, antismash=True, check=True) == 2
    assert runner.prereq_calls()


def test_check_image_record_fails_when_image_is_unavailable(tmp_path: Path) -> None:
    install_from_image(tmp_path)
    assert setup(tmp_path, runner=PrereqRunner(unavailable=True), antismash=True, check=True) == 2


def test_check_without_image_record_still_fails(tmp_path: Path) -> None:
    """Prereqs alone do not pass ``--check``: nothing was recorded by setup."""
    assert setup(tmp_path, runner=PrereqRunner(0), antismash=True, check=True) == 2


def test_check_image_does_not_rewrite_or_download(tmp_path: Path) -> None:
    install_from_image(tmp_path)
    before = (tmp_path / "db_manifest.json").read_bytes()
    assert setup(tmp_path, runner=PrereqRunner(0), antismash=True, check=True) == 0
    assert (tmp_path / "db_manifest.json").read_bytes() == before
    assert not (tmp_path / sources()["antismash"].path).exists()


def test_rerun_after_image_install_stays_image(tmp_path: Path) -> None:
    install_from_image(tmp_path)
    assert setup(tmp_path, runner=PrereqRunner(0), antismash=True) == 0
    from gippa.setupdb import read_records

    assert read_records(tmp_path)["antismash"].url == "image"


# ---------------------------------------------------------------- 4. body without Content-Length


class _Sock:
    def __init__(self, raw: bytes) -> None:
        self.raw = raw

    def makefile(self, *_a: Any, **_k: Any) -> io.BytesIO:
        return io.BytesIO(self.raw)


def _chunked(body: bytes, size: int = 1000) -> bytes:
    parts = [body[i : i + size] for i in range(0, len(body), size)]
    return b"".join(b"%x\r\n%s\r\n" % (len(p), p) for p in parts) + b"0\r\n\r\n"


class NoLengthNet:
    """Fake ``http.client`` that serves 200 bodies with no Content-Length header: either
    close-delimited or ``Transfer-Encoding: chunked``."""

    def __init__(self, mode: str) -> None:
        self.mode = mode
        self.routes: dict[tuple[str, str], bytes] = {}

    def conn(self) -> type:
        net = self

        class Conn(http.client.HTTPConnection):
            def __init__(self, host: str, port: int | None = None, **_kw: Any) -> None:
                super().__init__(host, port)
                self.path = ""

            def connect(self) -> None:
                raise AssertionError("never opens a socket")

            def close(self) -> None: ...

            def request(self, method: str, url: str, *_a: Any, **_k: Any) -> None:  # type: ignore[override]
                self.path = url

            def getresponse(self) -> http.client.HTTPResponse:
                body = net.routes.get((self.host, self.path.split("?")[0]))
                if body is None:
                    raw = b"HTTP/1.1 404 X\r\nContent-Length: 0\r\n\r\n"
                elif net.mode == "chunked":
                    raw = (b"HTTP/1.1 200 OK\r\nTransfer-Encoding: chunked\r\n\r\n"
                           + _chunked(body))  # fmt: skip
                else:
                    raw = b"HTTP/1.1 200 OK\r\nConnection: close\r\n\r\n" + body
                resp = http.client.HTTPResponse(_Sock(raw), method="GET")  # type: ignore[arg-type]
                resp.begin()
                assert resp.getheader("Content-Length") is None
                return resp

        return Conn


@pytest.fixture(params=["close", "chunked"])
def nolength(
    request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch
) -> Iterator[NoLengthNet]:
    import gippa.setupdb.run as run_mod

    request.getfixturevalue("real_https_fetch")
    fake = NoLengthNet(request.param)
    monkeypatch.setattr(http.client, "HTTPSConnection", fake.conn())
    monkeypatch.setattr(http.client, "HTTPConnection", fake.conn())
    monkeypatch.setattr(run_mod, "GIB", SMALL_GIB)
    yield fake


def ko_body(lines: int) -> bytes:
    return b"".join(b"ko:K%05d\tpath:map00010\n" % (i % 100000) for i in range(lines))


def test_unbounded_body_over_one_gib_default_is_a_source_error(
    nolength: NoLengthNet, tmp_path: Path
) -> None:
    from gippa.setupdb import HttpsFetcher, SourceError

    nolength.routes[("example.org", "/big")] = b"x" * (3 * SMALL_GIB)
    dest = tmp_path / "f"
    with pytest.raises(SourceError) as exc:
        HttpsFetcher().fetch("https://example.org/big", dest)
    assert exc.value.host == "example.org"
    assert not dest.exists()


def test_body_without_length_under_the_cap_is_accepted(
    nolength: NoLengthNet, tmp_path: Path
) -> None:
    from gippa.setupdb import HttpsFetcher

    body = ko_body(50)  # 1200 bytes < 1 "GiB" (4096)
    nolength.routes[("rest.kegg.jp", KO_PATH)] = body
    HttpsFetcher().fetch(KO_URL, tmp_path / "f")
    assert (tmp_path / "f").read_bytes() == body


def test_install_stops_at_default_cap_when_size_gib_unset(
    nolength: NoLengthNet, tmp_path: Path
) -> None:
    from gippa.setupdb import HttpsFetcher, install

    item = sources()["kegg_ko_pathway"]
    assert item.size_gib == 0
    nolength.routes[("rest.kegg.jp", KO_PATH)] = ko_body(300)  # 7200 bytes > 4096
    rec = install(item, tmp_path, NoRunner(), HttpsFetcher())
    assert rec.status == "absent" and "rest.kegg.jp" in rec.message, rec
    assert [p for p in files_under(tmp_path) if not p.name.endswith(".lock")] == []


def test_install_cap_follows_item_size_gib(
    nolength: NoLengthNet, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.setupdb import HttpsFetcher, install

    body = ko_body(300)  # 7200 bytes
    nolength.routes[("rest.kegg.jp", KO_PATH)] = body
    big = repin(monkeypatch, "kegg_ko_pathway", size_gib=2)  # cap 8192
    rec = install(big, tmp_path, NoRunner(), HttpsFetcher())
    assert rec.status == "ok", rec
    assert (tmp_path / big.path).read_bytes() == body
    over = ko_body(400)  # 9600 bytes > 8192
    nolength.routes[("rest.kegg.jp", KO_PATH)] = over
    (tmp_path / big.path).unlink()
    (tmp_path / "db_manifest.json").unlink(missing_ok=True)
    rec = install(big, tmp_path, NoRunner(), HttpsFetcher())
    assert rec.status == "absent" and "rest.kegg.jp" in rec.message, rec
    assert not (tmp_path / big.path).exists()


def test_run_setup_with_unbounded_kegg_body_exits_2(nolength: NoLengthNet, tmp_path: Path) -> None:
    from gippa.setupdb import HttpsFetcher

    nolength.routes[("rest.kegg.jp", KO_PATH)] = ko_body(300)
    nolength.routes[("rest.kegg.jp", "/list/pathway")] = b"map00010\tGlycolysis\n"
    assert setup(tmp_path, HttpsFetcher(), kegg=True) == 2
    assert not (tmp_path / sources()["kegg_ko_pathway"].path).exists()
    assert [p for p in tmp_path.rglob("*") if ".partial" in p.name] == []


# ------------------------------------------------- 5. _load_sources rejects unsafe paths


class _FakeRes:
    def __init__(self, text: str) -> None:
        self.text = text

    def joinpath(self, _name: str) -> _FakeRes:
        return self

    def read_text(self, *_a: Any, **_k: Any) -> str:
        return self.text


def _yaml(path: str) -> str:
    return (
        "bad_item:\n  version: '1'\n  url: https://example.org/x\n  checksum: ''\n"
        f"  path: {path!r}\n  tool: null\n  check: sha256\n"
    )


@pytest.mark.parametrize(
    "path", ["../outside", "/etc/passwd", "go/../../escape", "a/../../b", "..", ""]
)
def test_load_sources_rejects_unsafe_path(path: str, monkeypatch: pytest.MonkeyPatch) -> None:
    import gippa.setupdb.items as items_mod

    monkeypatch.setattr(items_mod, "files", lambda _pkg: _FakeRes(_yaml(path)))
    with pytest.raises(ValueError) as exc:
        items_mod._load_sources()
    if path:
        assert path in str(exc.value)


@pytest.mark.parametrize("path", ["go/go-basic.obo", "bakta/db-light", "a/b/c.tsv"])
def test_load_sources_accepts_safe_path(path: str, monkeypatch: pytest.MonkeyPatch) -> None:
    import gippa.setupdb.items as items_mod

    monkeypatch.setattr(items_mod, "files", lambda _pkg: _FakeRes(_yaml(path)))
    assert items_mod._load_sources()["bad_item"].path == path


def test_shipped_sources_paths_are_all_safe() -> None:
    from pathlib import PurePosixPath

    for item in sources().values():
        p = PurePosixPath(item.path)
        assert item.path and not p.is_absolute() and ".." not in p.parts, item


# ---------------------------------------------------------------- 6. manifest mode 0644


@pytest.fixture(params=[0o022, 0o077, 0o002], ids=["022", "077", "002"])
def umask(request: pytest.FixtureRequest) -> Iterator[int]:
    old = os.umask(request.param)
    try:
        yield request.param
    finally:
        os.umask(old)


def _mode(p: Path) -> int:
    return stat.S_IMODE(p.stat().st_mode)


def test_write_manifest_is_0644(tmp_path: Path, umask: int) -> None:
    from gippa.setupdb import DbRecord, write_manifest

    rec = DbRecord("x", "1", "https://example.org/x", "", 0, "", "x", "absent", "")
    write_manifest(tmp_path, [rec])
    assert _mode(tmp_path / "db_manifest.json") == 0o644
    write_manifest(tmp_path, [dataclasses.replace(rec, version="2")])  # rewrite keeps the mode
    assert _mode(tmp_path / "db_manifest.json") == 0o644


def test_run_setup_manifest_is_0644(tmp_path: Path, umask: int) -> None:
    install_from_image(tmp_path)
    assert _mode(tmp_path / "db_manifest.json") == 0o644
