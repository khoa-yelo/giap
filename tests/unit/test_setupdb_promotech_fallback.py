"""Promotech RF-HOT falls back to the authors' GitHub LFS copy (issue #508, ADR-0069).

Interface pinned here (for the implementer):
- ``gippa/data/db_sources.yaml`` gains, for ``promotech_models`` only::

      alt_urls:
        - https://media.githubusercontent.com/media/BioinformaticsLabAtMUN/PromoTech/bdf2f3ce06ffd0e342607eed4bfa32e160816143/models/RF-HOT.zip

  and ``DbItem`` gains ``alt_urls: tuple[str, ...] = ()``. No other item has ``alt_urls``; the
  primary ``url``/``host`` stay ``www.cs.mun.ca``.
- ``install`` tries ``url`` first and each ``alt_urls`` entry only after a ``SourceError``
  (timeout, connection refused, HTTP 5xx/404, refused redirect...). A checksum mismatch (or an
  unsafe archive) is final: no further source is contacted. A full disk is ours, not the
  source's: it is not a reason to try another source either.
- Every source is checked against the same ADR-0014 pin.
- The record (``DbRecord.url`` in ``db_manifest.json``) and the ``downloaded from`` message name
  the URL/host that served the file.
- When every source fails, the message names every host tried and the manual route: put
  ``RF-HOT.zip`` and its unzipped ``RF-HOT.model`` in ``<db-dir>/promotech/`` (the real path is
  printed), then rerun ``gippa setup``. Only the hosts, the path, the two file names and
  "rerun" are pinned, not the wording.
- A previous ``ok`` record whose url is an official alt URL is carried forward like one with the
  primary URL (idempotent reruns); ``_redirect_hosts`` and ``_body_cap`` treat an alt URL as the
  item's own URL. The LFS URL answers directly, so it needs no redirect allowlist entry.

Ambiguities resolved here (test author):
- A prefetch row of ``promotech_rf_hot`` whose source host is the alt host and whose checksum is
  the pin is adopted (``ok``, the row's date), as for the primary host.
- The message after a fallback must name the alt host; it may also say why cs.mun.ca failed.
"""

from __future__ import annotations

import dataclasses
import errno
from collections.abc import Mapping
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pytest
import yaml
from setup_guard import real_https_fetch, setup_guard  # noqa: F401  (autouse guard)
from test_setupdb import NoFetcher, NoRunner, _zip, items, pin, setup, text
from test_setupdb_review import FakeNet, net  # noqa: F401  (fixture)

from gippa.exitcodes import ExitCode

GIB = 1024**3
PRIMARY = "https://www.cs.mun.ca/~lourdes/public/PromoTech_models/RF-HOT.zip"
PRIMARY_HOST = "www.cs.mun.ca"
ALT = (
    "https://media.githubusercontent.com/media/BioinformaticsLabAtMUN/PromoTech/"
    "bdf2f3ce06ffd0e342607eed4bfa32e160816143/models/RF-HOT.zip"
)
ALT_HOST = "media.githubusercontent.com"
PAYLOAD = _zip(b"the published RF-HOT model (fake)")
BAD = _zip(b"not the published RF-HOT model")


@pytest.fixture(autouse=True)
def _env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Plenty of disk, inside Slurm (as in test_setupdb)."""
    import shutil
    from typing import NamedTuple

    class _Usage(NamedTuple):
        total: int
        used: int
        free: int

    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    monkeypatch.setenv("SLURM_JOB_ID", "123")


@pytest.fixture
def item(monkeypatch: pytest.MonkeyPatch) -> Any:
    """``promotech_models`` re-pinned to the sha256 of the fake ``PAYLOAD`` (alt_urls kept)."""
    return pin(monkeypatch, "promotech_models", PAYLOAD)


class RouteFetcher:
    """``routes[url]`` is the body served or the exception raised; any other URL is a test
    failure (it must never be contacted)."""

    def __init__(self, routes: Mapping[str, bytes | BaseException]) -> None:
        self.routes = dict(routes)
        self.urls: list[str] = []

    def fetch(self, url: str, dest: Path) -> None:
        self.urls.append(url)
        if url not in self.routes:
            raise AssertionError(f"must not be contacted: {url}")
        got = self.routes[url]
        if isinstance(got, BaseException):
            raise got
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(got)


def down(url: str, reason: str) -> BaseException:
    from gippa.setupdb import SourceError

    return SourceError(urlparse(url).hostname or url, reason)


REASONS = ["timed out", "[Errno 111] Connection refused", "HTTP 503 for " + PRIMARY,
           "HTTP 404 for " + PRIMARY]  # fmt: skip


def run_promotech(db: Path, fetcher: Any) -> int:
    return setup(db, NoRunner(), fetcher, bakta="skip", promotech_models=True)


def promotech_dir(db: Path) -> Path:
    from gippa.setupdb import SOURCES

    return db / SOURCES["promotech_models"].path


def leftovers(db: Path) -> list[Path]:
    return [p for p in db.rglob("*") if ".partial" in p.name]


# ---------------------------------------------------------------- the source table


def test_promotech_has_the_authors_lfs_copy_as_its_only_alt_url() -> None:
    from gippa.setupdb import SOURCES

    row = SOURCES["promotech_models"]
    assert row.alt_urls == (ALT,)
    assert (row.url, row.host) == (PRIMARY, PRIMARY_HOST), "cs.mun.ca stays first"
    assert urlparse(ALT).scheme == "https" and urlparse(ALT).hostname == ALT_HOST


def test_no_other_item_has_alt_urls() -> None:
    from gippa.setupdb import SOURCES

    others = {n: i.alt_urls for n, i in SOURCES.items() if n != "promotech_models"}
    assert others and all(v == () for v in others.values()), others
    from importlib.resources import files

    raw = yaml.safe_load(files("gippa.data").joinpath("db_sources.yaml").read_text("utf-8"))
    assert [n for n, row in raw.items() if "alt_urls" in row] == ["promotech_models"]


def test_dbitem_alt_urls_defaults_to_empty() -> None:
    from gippa.setupdb import SOURCES, DbItem

    base = {f.name: getattr(SOURCES["genomad"], f.name) for f in dataclasses.fields(DbItem)
            if f.name != "alt_urls"}  # fmt: skip
    assert DbItem(**base).alt_urls == ()


# ---------------------------------------------------------------- behaviours 1-2: the fallback


@pytest.mark.parametrize("reason", REASONS)
def test_unreachable_cs_mun_ca_falls_back_to_the_lfs_copy(
    tmp_path: Path, item: Any, reason: str
) -> None:
    from gippa.setupdb import install

    fetcher = RouteFetcher({PRIMARY: down(PRIMARY, reason), ALT: PAYLOAD})
    rec = install(item, tmp_path, NoRunner(), fetcher)
    assert fetcher.urls == [PRIMARY, ALT], "cs.mun.ca first, then the LFS copy"
    assert rec.status == "ok", rec
    assert rec.url == ALT, "the record names the URL that served the file"
    assert "downloaded from" in rec.message and ALT_HOST in rec.message, rec.message
    d = promotech_dir(tmp_path)
    assert (d / "RF-HOT.zip").read_bytes() == PAYLOAD
    assert (d / "RF-HOT.model").is_file()
    assert leftovers(tmp_path) == []


def test_fallback_is_recorded_in_the_manifest_and_reported(
    tmp_path: Path,
    item: Any,
    capsys: pytest.CaptureFixture[str],
    caplog: pytest.LogCaptureFixture,
) -> None:
    fetcher = RouteFetcher({PRIMARY: down(PRIMARY, "timed out"), ALT: PAYLOAD})
    assert run_promotech(tmp_path, fetcher) == ExitCode.SUCCESS
    rec = items(tmp_path)["promotech_models"]
    assert (rec["status"], rec["url"], rec["checksum"]) == ("ok", ALT, item.checksum), rec
    out = text(capsys, caplog)
    assert "promotech_models" in out and "downloaded from" in out and ALT_HOST in out, out


@pytest.mark.parametrize("status", [404, 500, 503])
def test_http_error_from_cs_mun_ca_falls_back_over_real_https(
    net: FakeNet,  # noqa: F811
    tmp_path: Path,
    item: Any,
    status: int,
) -> None:
    from gippa.setupdb import HttpsFetcher, install

    net.routes[("https", PRIMARY_HOST, urlparse(PRIMARY).path)] = (status, {}, b"down")
    net.routes[("https", ALT_HOST, urlparse(ALT).path)] = (200, {}, PAYLOAD)
    rec = install(item, tmp_path, NoRunner(), HttpsFetcher())
    assert rec.status == "ok" and rec.url == ALT, rec
    assert net.contacted == [("https", PRIMARY_HOST, urlparse(PRIMARY).path),
                             ("https", ALT_HOST, urlparse(ALT).path)]  # fmt: skip


# ---------------------------------------------------------------- behaviour 3: primary first


def test_good_cs_mun_ca_never_contacts_the_lfs_copy(
    tmp_path: Path,
    item: Any,
    capsys: pytest.CaptureFixture[str],
    caplog: pytest.LogCaptureFixture,
) -> None:
    fetcher = RouteFetcher({PRIMARY: PAYLOAD})  # the alt URL would fail the test
    assert run_promotech(tmp_path, fetcher) == ExitCode.SUCCESS
    assert fetcher.urls == [PRIMARY]
    rec = items(tmp_path)["promotech_models"]
    assert (rec["status"], rec["url"]) == ("ok", PRIMARY), rec
    assert f"downloaded from {PRIMARY_HOST}" in rec["message"], rec
    assert ALT_HOST not in text(capsys, caplog)


# ---------------------------------------------------------------- behaviours 4-5: mismatch is final


def test_checksum_mismatch_from_cs_mun_ca_is_final(tmp_path: Path, item: Any) -> None:
    from gippa.setupdb import install

    fetcher = RouteFetcher({PRIMARY: BAD})  # contacting ALT fails the test
    rec = install(item, tmp_path, NoRunner(), fetcher)
    assert fetcher.urls == [PRIMARY], "a mismatch never falls through to another source"
    assert rec.status == "mismatch", rec
    assert "checksum mismatch" in rec.message and PRIMARY_HOST in rec.message, rec.message
    assert not list(tmp_path.rglob("RF-HOT.model")) and leftovers(tmp_path) == []


def test_checksum_mismatch_from_cs_mun_ca_exits_2_over_real_https(
    net: FakeNet,  # noqa: F811
    tmp_path: Path,
    item: Any,
) -> None:
    from gippa.setupdb import HttpsFetcher

    net.routes[("https", PRIMARY_HOST, urlparse(PRIMARY).path)] = (200, {}, BAD)
    net.routes[("https", ALT_HOST, urlparse(ALT).path)] = (200, {}, PAYLOAD)
    assert run_promotech(tmp_path, HttpsFetcher()) == ExitCode.ENVIRONMENT
    assert items(tmp_path)["promotech_models"]["status"] == "mismatch"
    assert [h for _s, h, _p in net.contacted] == [PRIMARY_HOST]


def test_checksum_mismatch_from_the_lfs_copy_is_a_mismatch(tmp_path: Path, item: Any) -> None:
    from gippa.setupdb import install

    fetcher = RouteFetcher({PRIMARY: down(PRIMARY, "timed out"), ALT: BAD})
    rec = install(item, tmp_path, NoRunner(), fetcher)
    assert fetcher.urls == [PRIMARY, ALT]
    assert rec.status == "mismatch", rec
    assert "checksum mismatch" in rec.message and ALT_HOST in rec.message, rec.message
    assert not list(tmp_path.rglob("RF-HOT.model")) and leftovers(tmp_path) == []
    assert not list(tmp_path.rglob("RF-HOT.zip"))


def test_full_disk_on_cs_mun_ca_does_not_try_the_lfs_copy(tmp_path: Path, item: Any) -> None:
    """Our disk, not the source (M11-15c1): the item fails and no other source is contacted."""
    from gippa.setupdb import install

    full = OSError(errno.ENOSPC, "No space left on device", str(tmp_path / "x"))
    fetcher = RouteFetcher({PRIMARY: full})
    rec = install(item, tmp_path, NoRunner(), fetcher)
    assert fetcher.urls == [PRIMARY]
    assert rec.status == "failed", rec


# ---------------------------------------------------------------- behaviour 6: both down


@pytest.mark.parametrize("reason", REASONS)
def test_both_sources_down_names_both_hosts_and_the_manual_route(
    tmp_path: Path,
    item: Any,
    reason: str,
    capsys: pytest.CaptureFixture[str],
    caplog: pytest.LogCaptureFixture,
) -> None:
    fetcher = RouteFetcher({PRIMARY: down(PRIMARY, reason), ALT: down(ALT, "timed out")})
    assert run_promotech(tmp_path, fetcher) == ExitCode.ENVIRONMENT
    assert fetcher.urls == [PRIMARY, ALT]
    rec = items(tmp_path)["promotech_models"]
    assert rec["status"] == "absent", rec
    msg = rec["message"]
    assert PRIMARY_HOST in msg and ALT_HOST in msg, msg
    assert "RF-HOT.zip" in msg and "RF-HOT.model" in msg, msg
    assert str(promotech_dir(tmp_path)) in msg, "the manual route names the real directory"
    assert "rerun" in msg.lower() and "gippa setup" in msg, msg
    out = text(capsys, caplog)
    assert PRIMARY_HOST in out and ALT_HOST in out and "rerun" in out.lower(), out
    assert leftovers(tmp_path) == [] and not promotech_dir(tmp_path).exists()


def test_manual_route_of_the_message_is_adopted_on_rerun(tmp_path: Path, item: Any) -> None:
    """The advice works: the two files dropped in ``<db-dir>/promotech/`` verify offline."""
    d = promotech_dir(tmp_path)
    d.mkdir(parents=True)
    (d / "RF-HOT.zip").write_bytes(PAYLOAD)
    (d / "RF-HOT.model").write_bytes(b"the published RF-HOT model (fake)")
    assert run_promotech(tmp_path, NoFetcher()) == ExitCode.SUCCESS
    assert items(tmp_path)["promotech_models"]["status"] == "ok"


# ---------------------------------------------------------------- behaviours 7-8: provenance


def _installed(db: Path) -> None:
    d = promotech_dir(db)
    d.mkdir(parents=True)
    (d / "RF-HOT.zip").write_bytes(PAYLOAD)
    (d / "RF-HOT.model").write_bytes(b"the published RF-HOT model (fake)")


def _old_record(db: Path, item: Any, url: str) -> None:
    from gippa.setupdb import DbRecord, write_manifest

    write_manifest(db, [DbRecord("promotech_models", item.version, url, item.checksum, 1,
                                 "2026-01-02", item.path, "ok")])  # fmt: skip


def test_previous_lfs_record_is_kept_on_rerun_and_passes_check(tmp_path: Path, item: Any) -> None:
    _installed(tmp_path)
    _old_record(tmp_path, item, ALT)
    assert run_promotech(tmp_path, NoFetcher()) == ExitCode.SUCCESS
    rec = items(tmp_path)["promotech_models"]
    assert (rec["status"], rec["url"], rec["retrieved"]) == ("ok", ALT, "2026-01-02"), rec
    code = setup(tmp_path, bakta="skip", promotech_models=True, check=True)
    assert code == ExitCode.SUCCESS
    assert items(tmp_path)["promotech_models"]["url"] == ALT, "--check writes nothing"


def test_previous_record_from_an_unlisted_url_is_not_carried_forward(
    tmp_path: Path, item: Any
) -> None:
    _installed(tmp_path)
    _old_record(tmp_path, item, "https://evil.example/RF-HOT.zip")
    assert run_promotech(tmp_path, NoFetcher()) == ExitCode.SUCCESS
    rec = items(tmp_path)["promotech_models"]
    assert rec["url"] == PRIMARY and rec["retrieved"] != "2026-01-02", rec


def test_prefetch_row_from_the_lfs_host_with_the_pin_is_adopted(tmp_path: Path, item: Any) -> None:
    _installed(tmp_path)
    (tmp_path / "prefetch_manifest.tsv").write_text(
        "item\tsource\tfile\tchecksum\tverified\tdate\n"
        f"promotech_rf_hot\t{ALT}\tRF-HOT.zip\t{item.checksum}\tlfs sha256\t2026-10-30T01:02:03Z\n"
    )
    assert run_promotech(tmp_path, NoFetcher()) == ExitCode.SUCCESS
    rec = items(tmp_path)["promotech_models"]
    assert (rec["status"], rec["retrieved"]) == ("ok", "2026-10-30"), rec


def test_lfs_url_gets_the_item_body_cap_and_no_redirect_allowlist() -> None:
    from gippa.setupdb import SOURCES
    from gippa.setupdb.run import _body_cap, _redirect_hosts

    assert _body_cap(ALT) == SOURCES["promotech_models"].size_gib * GIB == _body_cap(PRIMARY)
    assert _redirect_hosts(ALT, ALT_HOST) == frozenset(), "the LFS URL serves directly"


def test_same_host_redirect_on_the_lfs_host_is_followed(
    net: FakeNet,  # noqa: F811
    tmp_path: Path,
    item: Any,
) -> None:
    from gippa.setupdb import HttpsFetcher, install

    net.routes[("https", ALT_HOST, urlparse(ALT).path)] = (
        302, {"Location": f"https://{ALT_HOST}/served/RF-HOT.zip"}, b"")  # fmt: skip
    net.routes[("https", ALT_HOST, "/served/RF-HOT.zip")] = (200, {}, PAYLOAD)
    rec = install(item, tmp_path, NoRunner(), HttpsFetcher())  # cs.mun.ca: FakeNet 404
    assert rec.status == "ok" and rec.url == ALT, rec


def test_off_host_redirect_from_the_lfs_host_is_refused(
    net: FakeNet,  # noqa: F811
    tmp_path: Path,
    item: Any,
) -> None:
    from gippa.setupdb import HttpsFetcher, install

    net.routes[("https", ALT_HOST, urlparse(ALT).path)] = (
        302, {"Location": "https://evil.example/RF-HOT.zip"}, b"")  # fmt: skip
    net.routes[("https", "evil.example", "/RF-HOT.zip")] = (200, {}, PAYLOAD)
    rec = install(item, tmp_path, NoRunner(), HttpsFetcher())
    assert rec.status != "ok", rec
    assert ("https", "evil.example", "/RF-HOT.zip") not in net.contacted
    assert ("https", ALT_HOST, urlparse(ALT).path) in net.contacted, "the alt URL was tried"
    assert not promotech_dir(tmp_path).exists()
