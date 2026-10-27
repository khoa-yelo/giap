"""`gippa setup` after code review (PR #57; docs/specs/M2-08 "Lead decisions after code review" 1-6).

Every test runs under ``setup_guard`` (decision 3): no network, no tool, no real $GIPPA_DB_DIR.

Ambiguities resolved here:
- Redirect handling is exercised through fake ``http.client`` connections (``HttpsFetcher`` may use
  ``urllib`` or ``http.client`` directly; either reaches these classes). A redirect that is not
  allowed must be refused *before* the new location is contacted.
- An unsafe archive (zip or tar) is never ``ok``: the item stays absent (``status`` ``mismatch`` or
  ``absent``) and nothing is written outside ``<path>.partial``.
- A previous manifest record is reused (idempotent reruns keep their date) only when it carries the
  pinned URL and checksum; a record with another checksum or URL is not carried forward (R6).
- The geNOMAD tarball is deleted after it is unpacked (§3.9 "delete the archive after its checksum
  passes").
"""

from __future__ import annotations

import dataclasses
import hashlib
import http.client
import io
import json
import re
import shutil
import socket
import subprocess
import tempfile
import zipfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any, NamedTuple

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from setup_guard import BLOCKED, GuardError, real_https_fetch, setup_guard  # noqa: F401
from test_setupdb import (
    GENOMAD_MD5,
    GENOMAD_TAR,
    FakeRunner,
    NoFetcher,
    NoRunner,
    ZipFetcher,
    _bakta_version,
    _tar,
    make_bakta,
    make_genomad,
    md5,
    setup,
)

GIB = 1024**3
REPO = Path(__file__).resolve().parents[2]
ZENODO_GENOMAD = "https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz"


class _Usage(NamedTuple):
    total: int
    used: int
    free: int


@pytest.fixture(autouse=True)
def _env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Plenty of disk, inside Slurm, and tool lookups pointed at empty temporary roots."""
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    roots = {"GIPPA_ROOT": "root", "MAMBA_ROOT_PREFIX": "mamba",
             "NXF_APPTAINER_CACHEDIR": "cache", "GIPPA_DB_DIR": "not-the-real-db"}  # fmt: skip
    for var, sub in roots.items():
        monkeypatch.setenv(var, str(tmp_path / "_env" / sub))
    monkeypatch.setenv("GIPPA_REPO", str(REPO))


def items(db: Path) -> dict[str, dict[str, Any]]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    return {i["name"]: i for i in data["items"]}


def leftovers(db: Path) -> list[Path]:
    return [p for p in db.rglob("*") if ".partial" in p.name or p.name.endswith((".tar.gz", ".gz"))]


def sources() -> Any:
    from gippa.setupdb import SOURCES

    return SOURCES


def genomad_item(payload: bytes) -> Any:
    return dataclasses.replace(sources()["genomad"], checksum=md5(payload))


# ---------------------------------------------------------------- decision 3: the guard


def test_guard_blocks_tools_subprocesses_sockets_and_https() -> None:
    from gippa.setupdb import HttpsFetcher, ToolRunner, ToolUnavailable

    with pytest.raises((GuardError, ToolUnavailable)):
        ToolRunner().run("bakta", ["bakta_db", "--help"])
    with pytest.raises(GuardError):
        HttpsFetcher().fetch(ZENODO_GENOMAD, Path("never-written"))
    with pytest.raises(GuardError):
        subprocess.run(["true"], check=False)
    with pytest.raises(GuardError):
        subprocess.Popen(["true"])
    with pytest.raises(GuardError):
        socket.create_connection(("zenodo.org", 443))
    assert not Path("never-written").exists()


def test_cli_genomad_download_goes_through_the_https_fetcher(tmp_path: Path) -> None:
    """The CLI wires geNOMAD to ``HttpsFetcher`` (decision 1); the guard stops it before the net."""
    from click.testing import CliRunner

    from gippa.cli import cli

    argv = ["setup", "--db-dir", str(tmp_path / "db"), "--bakta", "skip", "--genomad",
            "--no-antismash", "--no-promotech-models"]  # fmt: skip
    r = CliRunner().invoke(cli, argv)
    assert r.exit_code != 0
    assert any(b.startswith("HttpsFetcher.fetch") for b in BLOCKED), BLOCKED
    assert not any(b.startswith("ToolRunner.run") for b in BLOCKED), BLOCKED
    assert not (tmp_path / "db" / "genomad_db").exists()


# ---------------------------------------------------------------- decision 1: geNOMAD over HTTPS


def test_genomad_source_is_the_pinned_zenodo_tarball_without_a_tool() -> None:
    item = sources()["genomad"]
    assert item.tool is None
    assert item.url.startswith("https://zenodo.org/records/14886553/")
    assert item.host == "zenodo.org" and item.checksum == GENOMAD_MD5
    assert "nersc" not in item.url


def test_genomad_fetched_md5_checked_unpacked_archive_deleted(tmp_path: Path) -> None:
    from gippa.setupdb import install

    item = genomad_item(GENOMAD_TAR)
    fetcher = ZipFetcher(GENOMAD_TAR)
    rec = install(item, tmp_path, NoRunner(), fetcher)
    assert rec.status == "ok" and rec.version == "1.9", rec
    assert fetcher.urls == [item.url]
    assert rec.checksum == item.checksum and rec.url == item.url
    assert (tmp_path / "genomad_db" / "version.txt").read_text().strip() == "1.9"
    assert leftovers(tmp_path) == []


def test_genomad_md5_mismatch_is_mismatch_and_nothing_unpacked(tmp_path: Path) -> None:
    from gippa.setupdb import install

    rec = install(sources()["genomad"], tmp_path, NoRunner(), ZipFetcher(GENOMAD_TAR))
    assert rec.status == "mismatch", rec
    assert "zenodo.org" in rec.message
    assert not (tmp_path / "genomad_db").exists()
    assert leftovers(tmp_path) == []


def test_genomad_md5_mismatch_exits_2(tmp_path: Path) -> None:
    assert setup(tmp_path, NoRunner(), ZipFetcher(GENOMAD_TAR), bakta="skip", genomad=True) == 2
    assert items(tmp_path)["genomad"]["status"] == "mismatch"


def test_genomad_download_is_verified_after_unpacking(tmp_path: Path) -> None:
    """Checksum right, content wrong (v1.8): never installed."""
    from gippa.setupdb import install

    payload = _tar([("genomad_db/version.txt", b"1.8\n")])
    rec = install(genomad_item(payload), tmp_path, NoRunner(), ZipFetcher(payload))
    assert rec.status != "ok", rec
    assert not (tmp_path / "genomad_db").exists()
    assert leftovers(tmp_path) == []


def _unsafe_tars(base: Path) -> list[tuple[str, bytes]]:
    good = ("genomad_db/version.txt", b"1.9\n")
    return [
        ("dotdot", _tar([good, ("../../escape.txt", b"x")])),
        ("dotdot-inside", _tar([good, ("genomad_db/../../../escape.txt", b"x")])),
        ("absolute", _tar([good, (str(base / "abs_escape.txt"), b"x")])),
        ("symlink-out", _tar([good], links=[("genomad_db/link", "../../..")])),
        ("symlink-abs", _tar([good], links=[("genomad_db/link", str(base))])),
    ]


@pytest.mark.parametrize(
    "case", ["dotdot", "dotdot-inside", "absolute", "symlink-out", "symlink-abs"]
)
def test_genomad_unsafe_tarball_is_rejected(tmp_path: Path, case: str) -> None:
    from gippa.setupdb import install

    db = tmp_path / "a" / "db"
    db.mkdir(parents=True)
    payload = dict(_unsafe_tars(tmp_path))[case]
    rec = install(genomad_item(payload), db, NoRunner(), ZipFetcher(payload))
    assert rec.status != "ok", rec
    assert not (db / "genomad_db").exists()
    escaped = [p for p in tmp_path.rglob("*escape*")]
    assert escaped == [], escaped
    assert not [p for p in db.rglob("*") if p.is_symlink()]
    assert leftovers(db) == []


@settings(
    max_examples=30, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    st.lists(st.sampled_from(["genomad_db", "a", "..", "b"]), min_size=1, max_size=5).filter(
        lambda parts: ".." in parts
    )
)
def test_any_member_with_dotdot_is_rejected(parts: list[str]) -> None:
    """Property: a member path with a ``..`` component never yields an ``ok`` install."""
    from gippa.setupdb import install

    payload = _tar([("genomad_db/version.txt", b"1.9\n"), ("/".join([*parts, "f.txt"]), b"x")])
    with tempfile.TemporaryDirectory() as tmp:
        db = Path(tmp) / "x" / "y" / "db"
        db.mkdir(parents=True)
        rec = install(genomad_item(payload), db, NoRunner(), ZipFetcher(payload))
        assert rec.status != "ok", (parts, rec)
        assert not (db / "genomad_db").exists()
        assert [p for p in Path(tmp).rglob("f.txt")] == []


# ---------------------------------------------------------------- decision 5: zip-slip (Promotech)


def _zip(members: Sequence[tuple[str, bytes]]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, data in members:
            z.writestr(name, data)
    return buf.getvalue()


@pytest.mark.parametrize("bad", ["../evil.model", "RF-HOT/../../evil.model", "/ABS/evil.model"])
def test_promotech_zip_slip_is_rejected(tmp_path: Path, bad: str) -> None:
    from gippa.setupdb import install

    db = tmp_path / "a" / "db"
    db.mkdir(parents=True)
    name = str(tmp_path) + bad[4:] if bad.startswith("/ABS") else bad
    payload = _zip([("RF-HOT.model", b"model"), (name, b"evil")])
    item = sources()["promotech_models"]
    item = dataclasses.replace(item, checksum="sha256:" + hashlib.sha256(payload).hexdigest())
    rec = install(item, db, NoRunner(), ZipFetcher(payload))
    assert rec.status != "ok", rec
    assert [p for p in tmp_path.rglob("evil.model")] == []
    assert not (db / item.path / "RF-HOT.model").exists()


# ---------------------------------------------------------------- decision 2: prefetch adoption


def _prefetch(db: Path, checksum: str, source: str, date: str = "2025-03-02T01:08:26Z") -> None:
    (db / "prefetch_manifest.tsv").write_text(
        "item\tsource\tfile\tchecksum\tverified\tdate\n"
        f"genomad_db\t{source}\tgenomad_db_v1.9.tar.gz\t{checksum}\tzenodo md5\t{date}\n"
    )


BAD_MD5 = "md5:" + "0" * 32


@pytest.mark.parametrize(
    ("checksum", "source"),
    [
        (BAD_MD5, ZENODO_GENOMAD + "?download=1"),
        (GENOMAD_MD5, "http://evil.example/genomad_db_v1.9.tar.gz"),
        (BAD_MD5, "http://evil.example/genomad_db_v1.9.tar.gz"),
        (GENOMAD_MD5, "https://evil.example/records/14886553/files/genomad_db_v1.9.tar.gz"),
    ],
)
def test_prefetch_row_differing_from_the_pin_is_a_mismatch(
    tmp_path: Path, checksum: str, source: str
) -> None:
    _prefetch(tmp_path, checksum, source)
    make_genomad(tmp_path)
    assert setup(tmp_path, NoRunner(), NoFetcher(), bakta="skip", genomad=True) == 2
    rec = items(tmp_path)["genomad"]
    assert rec["status"] == "mismatch", rec
    assert rec["checksum"] == GENOMAD_MD5, "the pin is recorded, not the prefetch log's checksum"
    assert "evil.example" not in rec["url"]
    assert rec["retrieved"] != "2025-03-02"


def test_prefetch_row_matching_the_pin_is_adopted_with_its_old_date(tmp_path: Path) -> None:
    _prefetch(tmp_path, GENOMAD_MD5, ZENODO_GENOMAD + "?download=1")
    make_genomad(tmp_path)
    assert setup(tmp_path, NoRunner(), NoFetcher(), bakta="skip", genomad=True) == 0
    rec = items(tmp_path)["genomad"]
    assert (rec["status"], rec["checksum"], rec["retrieved"]) == ("ok", GENOMAD_MD5, "2025-03-02")


# ---------------------------------------------------------------- decision 5: previous provenance


def _old_manifest(db: Path, **kw: Any) -> None:
    from gippa.setupdb import DbRecord, write_manifest

    base = {"name": "genomad", "version": "1.9", "url": sources()["genomad"].url,
            "checksum": GENOMAD_MD5, "size_bytes": 1, "retrieved": "2025-01-02",
            "path": "genomad_db", "status": "ok", "message": ""}  # fmt: skip
    write_manifest(db, [DbRecord(**(base | kw))])


def test_rerun_keeps_previous_ok_provenance(tmp_path: Path) -> None:
    make_genomad(tmp_path)
    _old_manifest(tmp_path)
    assert setup(tmp_path, NoRunner(), NoFetcher(), bakta="skip", genomad=True) == 0
    assert items(tmp_path)["genomad"]["retrieved"] == "2025-01-02"


@pytest.mark.parametrize(
    "kw", [{"checksum": BAD_MD5}, {"url": "http://evil.example/genomad_db_v1.9.tar.gz"}]
)
def test_previous_record_off_the_pin_is_not_carried_forward(
    tmp_path: Path, kw: dict[str, str]
) -> None:
    make_genomad(tmp_path)
    _old_manifest(tmp_path, **kw)
    setup(tmp_path, NoRunner(), NoFetcher(), bakta="skip", genomad=True)
    rec = items(tmp_path)["genomad"]
    assert rec["checksum"] == GENOMAD_MD5 and rec["url"] == sources()["genomad"].url, rec


# ---------------------------------------------------------------- decision 5: remove-before-rename


def _broken_light(db: Path) -> Path:
    d = make_bakta(db, "light", major=5)
    (d / "old.txt").write_text("stale v5 file")
    return d


def test_verified_download_replaces_a_stale_copy(tmp_path: Path) -> None:
    d = _broken_light(tmp_path)
    assert setup(tmp_path, FakeRunner(), bakta="light") == 0
    assert json.loads((d / "version.json").read_text())["major"] == 6
    assert not (d / "old.txt").exists(), "the stale copy is removed, not merged"
    assert leftovers(tmp_path) == []


def test_failed_download_leaves_the_stale_copy_alone(tmp_path: Path) -> None:
    d = _broken_light(tmp_path)
    assert setup(tmp_path, FakeRunner(exit_code=1), bakta="light") == 2
    assert (d / "old.txt").exists(), "nothing is removed before a download verifies"
    assert json.loads((d / "version.json").read_text())["major"] == 5
    assert leftovers(tmp_path) == []


def test_download_that_does_not_verify_leaves_the_stale_copy_alone(tmp_path: Path) -> None:
    class WrongVersion(FakeRunner):
        def run(self, tool: str, argv: Sequence[str]) -> int:
            code = super().run(tool, argv)
            out = Path(argv[list(argv).index("--output") + 1])
            for p in out.rglob("version.json"):
                p.write_text(_bakta_version("light", major=7))
            return code

    d = _broken_light(tmp_path)
    assert setup(tmp_path, WrongVersion(), bakta="light") == 2
    assert (d / "old.txt").exists()
    assert leftovers(tmp_path) == []


# ---------------------------------------------------------------- decision 5: ToolRunner.command


ARGV_TAIL = ["--type", "light"]


def _image(cache: Path, tool: str) -> Path:
    conf = (REPO / "workflow/conf/containers.config").read_text(encoding="utf-8")
    m = re.search(rf"withName: '{tool.upper()}' \{{ container = '(?:docker://)?([^']+)'", conf)
    assert m is not None
    img = cache / (re.sub(r"[/:]", "-", m.group(1)) + ".img")
    img.parent.mkdir(parents=True, exist_ok=True)
    img.write_bytes(b"")
    return img


def _env_dir(mamba: Path, tool: str) -> None:
    (mamba / "envs" / f"gippa-{tool}").mkdir(parents=True)


def test_command_prefers_the_conda_environment(tmp_path: Path) -> None:
    from gippa.setupdb import ToolRunner

    _env_dir(tmp_path / "_env" / "mamba", "bakta")
    _image(tmp_path / "_env" / "cache", "bakta")
    argv = ["bakta_db", "download", "--output", str(tmp_path), *ARGV_TAIL]
    assert ToolRunner().command("bakta", argv) == ["micromamba", "run", "-n", "gippa-bakta", *argv]


def test_command_falls_back_to_the_pinned_image_with_binds(tmp_path: Path) -> None:
    from gippa.setupdb import ToolRunner

    img = _image(tmp_path / "_env" / "cache", "bakta")
    out = tmp_path / "out"
    out.mkdir()
    argv = ["bakta_db", "download", "--output", str(out), *ARGV_TAIL]
    cmd = ToolRunner().command("bakta", argv)
    assert cmd[:2] == ["apptainer", "exec"]
    assert cmd[-len(argv) :] == argv
    i = cmd.index(str(img))
    assert cmd[i + 1 :] == argv, "the image comes right before the tool's argv"
    assert any(cmd[k] == "--bind" and cmd[k + 1].split(":")[0] == str(out) for k in range(i)), cmd
    assert "light" not in cmd[:i], "relative arguments are not bound"


def _no_docker(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """14f2: ToolRunner falls back to ``docker run``; these cases mean a host without docker."""
    nodocker = tmp_path / "nodocker"
    nodocker.mkdir()
    monkeypatch.setenv("PATH", str(nodocker))


def test_command_without_environment_or_image_is_unavailable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _no_docker(tmp_path, monkeypatch)
    from gippa.setupdb import ToolRunner, ToolUnavailable

    with pytest.raises(ToolUnavailable, match="build_envs"):
        ToolRunner().command("bakta", ["bakta_db", "download"])


def test_install_with_no_environment_is_absent_and_starts_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _no_docker(tmp_path, monkeypatch)
    from gippa.setupdb import ToolRunner, install

    rec = install(sources()["bakta_light"], tmp_path, ToolRunner(), NoFetcher())
    assert rec.status == "absent" and "build_envs" in rec.message, rec
    assert BLOCKED == []
    assert leftovers(tmp_path) == []


# ---------------------------------------------------------------- decision 4: https-only, redirects


class _Sock:
    def __init__(self, raw: bytes) -> None:
        self.raw = raw

    def makefile(self, *_a: Any, **_k: Any) -> io.BytesIO:
        return io.BytesIO(self.raw)


class FakeNet:
    """Fake ``http.client``: ``routes[(scheme, host, path)] = (status, headers, body)``."""

    def __init__(self) -> None:
        self.routes: dict[tuple[str, str, str], tuple[int, dict[str, str], bytes]] = {}
        self.contacted: list[tuple[str, str, str]] = []

    def conn(self, scheme: str) -> type:
        net = self

        class Conn(http.client.HTTPConnection):  # only the network part is replaced
            def __init__(self, host: str, port: int | None = None, **_kw: Any) -> None:
                super().__init__(host, port)
                self.path = ""

            def connect(self) -> None:
                raise AssertionError("FakeNet never opens a socket")

            def close(self) -> None: ...

            def request(self, method: str, url: str, *_a: Any, **_k: Any) -> None:  # type: ignore[override]
                self.path = url

            def getresponse(self) -> http.client.HTTPResponse:
                key = (scheme, self.host, self.path.split("?")[0])
                net.contacted.append(key)
                status, headers, body = net.routes.get(key, (404, {}, b"not found"))
                head = "".join(f"{k}: {v}\r\n" for k, v in headers.items())
                raw = (f"HTTP/1.1 {status} X\r\n{head}Content-Length: {len(body)}\r\n"
                       "Connection: close\r\n\r\n").encode() + body  # fmt: skip
                resp = http.client.HTTPResponse(_Sock(raw), method="GET")  # type: ignore[arg-type]
                resp.begin()
                return resp

        return Conn


@pytest.fixture
def net(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    request.getfixturevalue("real_https_fetch")
    fake = FakeNet()
    monkeypatch.setattr(http.client, "HTTPSConnection", fake.conn("https"))
    monkeypatch.setattr(http.client, "HTTPConnection", fake.conn("http"))
    return fake


FILE = "/records/14886553/files/genomad_db_v1.9.tar.gz"


def test_fake_net_serves_a_plain_https_download(net: FakeNet, tmp_path: Path) -> None:
    from gippa.setupdb import HttpsFetcher

    net.routes[("https", "zenodo.org", FILE)] = (200, {}, b"DATA")
    HttpsFetcher().fetch("https://zenodo.org" + FILE, tmp_path / "f")
    assert (tmp_path / "f").read_bytes() == b"DATA"


def test_redirect_on_the_official_host_over_https_is_followed(net: FakeNet, tmp_path: Path) -> None:
    from gippa.setupdb import HttpsFetcher

    net.routes[("https", "zenodo.org", FILE)] = (
        302,
        {"Location": "https://zenodo.org/api/records/14886553/files/x/content"},
        b"",
    )
    net.routes[("https", "zenodo.org", "/api/records/14886553/files/x/content")] = (200, {}, b"OK")
    HttpsFetcher().fetch("https://zenodo.org" + FILE, tmp_path / "f")
    assert (tmp_path / "f").read_bytes() == b"OK"


@pytest.mark.parametrize(
    "location",
    ["http://zenodo.org/api/x", "https://evil.example/x", "http://evil.example/x",
     "https://zenodo.org.evil.example/x", "ftp://zenodo.org/x"],
)  # fmt: skip
def test_redirect_off_https_or_off_host_is_an_error(
    net: FakeNet, tmp_path: Path, location: str
) -> None:
    from gippa.setupdb import HttpsFetcher, SourceError

    net.routes[("https", "zenodo.org", FILE)] = (302, {"Location": location}, b"")
    evil = [("http", "zenodo.org", "/api/x"), ("https", "evil.example", "/x"),
            ("http", "evil.example", "/x"), ("https", "zenodo.org.evil.example", "/x")]  # fmt: skip
    for k in evil:
        net.routes[k] = (200, {}, b"EVIL")
    with pytest.raises(SourceError):
        HttpsFetcher().fetch("https://zenodo.org" + FILE, tmp_path / "f")
    assert net.contacted == [("https", "zenodo.org", FILE)], "the new location is never contacted"
    assert not (tmp_path / "f").exists() or b"EVIL" not in (tmp_path / "f").read_bytes()


@pytest.mark.parametrize("url", ["http://zenodo.org" + FILE, "ftp://zenodo.org" + FILE])
def test_fetch_refuses_non_https_urls_without_contacting(
    net: FakeNet, tmp_path: Path, url: str
) -> None:
    from gippa.setupdb import HttpsFetcher, SourceError

    with pytest.raises(SourceError):
        HttpsFetcher().fetch(url, tmp_path / "f")
    assert net.contacted == []


def test_install_redirect_off_host_leaves_item_absent(net: FakeNet, tmp_path: Path) -> None:
    from gippa.setupdb import HttpsFetcher, install

    item = genomad_item(GENOMAD_TAR)
    path = "/" + item.url.split("/", 3)[3]
    net.routes[("https", "zenodo.org", path.split("?")[0])] = (
        302,
        {"Location": "https://evil.example/genomad.tar.gz"},
        b"",
    )
    net.routes[("https", "evil.example", "/genomad.tar.gz")] = (200, {}, GENOMAD_TAR)
    rec = install(item, tmp_path, NoRunner(), HttpsFetcher())
    assert rec.status != "ok" and "zenodo.org" in rec.message, rec
    assert ("https", "evil.example", "/genomad.tar.gz") not in net.contacted
    assert not (tmp_path / "genomad_db").exists()


# ---------------------------------------------------------------- decision 6: amrfinderplus row


def _amr(db: Path, bakta_path: str, version: str = "2026-08-07.1") -> None:
    d = db / bakta_path / "amrfinderplus-db" / "latest"
    d.mkdir(parents=True, exist_ok=True)
    (d / "version.txt").write_text(version + "\n")


def test_amrfinderplus_row_from_bakta_light_db(tmp_path: Path) -> None:
    make_bakta(tmp_path, "light")
    _amr(tmp_path, "bakta/db-light")
    assert setup(tmp_path, NoRunner(), bakta="light") == 0
    amr = items(tmp_path)["amrfinderplus"]
    assert amr["status"] == "ok" and amr["version"] == "2026-08-07.1"
    assert amr["path"].startswith("bakta/db-light/amrfinderplus-db")
    assert amr["url"].startswith("https://ftp.ncbi.nlm.nih.gov/")


def test_amrfinderplus_row_absent_for_light_without_version_file(tmp_path: Path) -> None:
    make_bakta(tmp_path, "light")
    setup(tmp_path, NoRunner(), bakta="light")
    got = items(tmp_path)
    assert got["bakta_light"]["status"] == "ok"
    assert got.get("amrfinderplus", {"status": "absent"})["status"] == "absent"
