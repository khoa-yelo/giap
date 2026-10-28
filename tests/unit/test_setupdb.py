"""`gippa setup` core (issue M2-08; GIPPA_SPEC.md §5.1, §7.9.1, §3.2, §3.9, R6).

Examples 1-8 and the edge cases of docs/specs/M2-08-gippa-setup-core.md, run with a fake command
runner and a fake fetcher: no test touches the network, a conda environment or $GIPPA_ROOT.

Ambiguities resolved here (the implementation may choose anything else the spec leaves open):
- The public names are importable from ``gippa.setupdb`` (re-exported from items/manifest/run).
- ``SetupOptions`` takes keyword fields named after the CLI options: ``db_dir``, ``bakta``,
  ``genomad``, ``antismash``, ``promotech_models``, ``check`` (``go``/``kegg``/``cog``/
  ``contaminants`` default to off).
- The fake runner does not know whether ``install`` points ``--output`` at ``<path>.partial`` or
  at its parent, so it writes the version file in both shapes (``OUT/version.json`` and
  ``OUT/db-light/version.json``); ``verify`` then decides.
- Where a database lives is ``SOURCES[name].path``; tests build layouts there, never at a
  hard-coded path (the prefetch put geNOMAD under ``genomad/genomad_db``, see the integration test).
- A checksum mismatch is recorded with ``status: "mismatch"`` (example 5) and never ``ok``; the
  item's files are left absent.
- Free disk is measured with ``shutil.disk_usage`` (patched here).
- Messages go to stderr; the tests read stderr and the log records together.
"""

from __future__ import annotations

import dataclasses
import hashlib
import io
import json
import re
import shutil
import tarfile
import zipfile
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, NamedTuple

import pytest
from setup_guard import setup_guard  # noqa: F401  (autouse: no network, no tools; decision 3)

from gippa.exitcodes import ExitCode

GIB = 1024**3
BAKTA_LIGHT_MD5 = "md5:4a6e059ded39e9c5537ef4137d2f5648"  # §3.9
GENOMAD_MD5 = "md5:67244b528bb8bed464d1ca147136d33e"  # Zenodo 14886553 record (prefetch log)
BAKTA_FULL_MD5 = "md5:4c1115e40abfa2b464ae5dd988bdd88e"  # ADR-0014 (Zenodo 14916843 metadata)
PROMOTECH_ZIP_SHA256 = (  # ADR-0014 first-download pin of RF-HOT.zip
    "sha256:fcd96ffe5a43b579af965839ee2bb896de94de3601978d88bcabd32c4fb5093f"
)
ISO_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")


# ---------------------------------------------------------------- fakes


class _Usage(NamedTuple):
    total: int
    used: int
    free: int


@pytest.fixture(autouse=True)
def big_disk(monkeypatch: pytest.MonkeyPatch) -> None:
    """Plenty of free disk unless a test says otherwise; inside a Slurm job by default."""
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    monkeypatch.setenv("SLURM_JOB_ID", "123")


def _out_dir(argv: Sequence[str]) -> Path:
    for flag in ("--output", "--database-dir"):
        if flag in argv:
            return Path(argv[list(argv).index(flag) + 1])
    return Path(argv[-1])  # genomad download-database DIR


def _bakta_version(kind: str, major: int = 6) -> str:
    return json.dumps({"major": major, "minor": 0, "type": kind})


class FakeRunner:
    """Records calls; fakes ``bakta_db download`` and ``genomad download-database``."""

    def __init__(self, exit_code: int = 0, hook: Callable[[], None] | None = None) -> None:
        self.calls: list[tuple[str, list[str]]] = []
        self.exit_code = exit_code
        self.hook = hook

    def run(self, tool: str, argv: Sequence[str]) -> int:
        self.calls.append((tool, list(argv)))
        if self.hook is not None:
            self.hook()
        out = _out_dir(argv)
        out.mkdir(parents=True, exist_ok=True)
        if self.exit_code != 0:
            (out / "half-written.tmp").write_text("partial")
            return self.exit_code
        if "bakta_db" in argv or any("bakta_db" in a for a in argv):
            kind = argv[list(argv).index("--type") + 1]
            sub = out / ("db-light" if kind == "light" else "db")
            sub.mkdir(parents=True, exist_ok=True)
            for d in (out, sub):
                (d / "version.json").write_text(_bakta_version(kind))
        elif "download-database" in argv:
            sub = out / "genomad_db"
            sub.mkdir(parents=True, exist_ok=True)
            for d in (out, sub):
                (d / "version.txt").write_text("1.9\n")
        return 0


class NoRunner:
    def run(self, tool: str, argv: Sequence[str]) -> int:
        raise AssertionError(f"runner must not be called: {tool} {argv}")


class NoFetcher:
    def fetch(self, url: str, dest: Path) -> None:
        raise AssertionError(f"fetcher must not be called: {url}")


def _zip(model_bytes: bytes) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("RF-HOT.model", model_bytes)
    return buf.getvalue()


class ZipFetcher:
    def __init__(self, payload: bytes) -> None:
        self.payload = payload
        self.urls: list[str] = []

    def fetch(self, url: str, dest: Path) -> None:
        self.urls.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(self.payload)


class DownFetcher:
    def fetch(self, url: str, dest: Path) -> None:
        from urllib.parse import urlparse

        from gippa.setupdb import SourceError

        raise SourceError(urlparse(url).hostname or url, "connection refused")


def _tar(members: Sequence[tuple[str, bytes]], links: Sequence[tuple[str, str]] = ()) -> bytes:
    """A gzipped tarball with the given regular files and symlinks (names kept verbatim)."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, data in members:
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
        for name, target in links:
            info = tarfile.TarInfo(name)
            info.type = tarfile.SYMTYPE
            info.linkname = target
            tar.addfile(info)
    return buf.getvalue()


GENOMAD_TAR = _tar([("genomad_db/version.txt", b"1.9\n"), ("genomad_db/genomad_db", b"mmseqs")])


def md5(payload: bytes) -> str:
    return "md5:" + hashlib.md5(payload).hexdigest()


def pin(monkeypatch: pytest.MonkeyPatch, name: str, payload: bytes) -> Any:
    """Re-pin ``SOURCES[name]`` to the checksum of a fake archive (same algorithm as the pin)."""
    import gippa.setupdb
    import gippa.setupdb.items
    import gippa.setupdb.run
    from gippa.setupdb import SOURCES

    algo = SOURCES[name].checksum.split(":", 1)[0]
    item = dataclasses.replace(
        SOURCES[name], checksum=f"{algo}:{hashlib.new(algo, payload).hexdigest()}"
    )
    patched = {**SOURCES, name: item}
    for mod in (gippa.setupdb, gippa.setupdb.items, gippa.setupdb.run):
        if hasattr(mod, "SOURCES"):
            monkeypatch.setattr(mod, "SOURCES", patched)
    return item


@pytest.fixture
def genomad_fetcher(monkeypatch: pytest.MonkeyPatch) -> ZipFetcher:
    """geNOMAD is an HTTPS fetch of the pinned Zenodo tarball (decision 1); fake one re-pinned."""
    pin(monkeypatch, "genomad", GENOMAD_TAR)
    return ZipFetcher(GENOMAD_TAR)


# ---------------------------------------------------------------- helpers


def opts(db: Path, **kw: Any) -> Any:
    from gippa.setupdb import SetupOptions

    base: dict[str, Any] = {
        "db_dir": db,
        "bakta": "light",
        "genomad": False,
        "antismash": False,
        "promotech_models": False,
        "check": False,
    }
    base.update(kw)
    return SetupOptions(**base)


def setup(db: Path, runner: Any = None, fetcher: Any = None, **kw: Any) -> int:
    from gippa.setupdb import run_setup

    return int(run_setup(opts(db, **kw), runner or NoRunner(), fetcher or NoFetcher()))


def manifest(db: Path) -> dict[str, Any]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def items(db: Path) -> dict[str, dict[str, Any]]:
    return {i["name"]: i for i in manifest(db)["items"]}


def src_path(name: str) -> str:
    from gippa.setupdb import SOURCES

    return SOURCES[name].path


def make_bakta(db: Path, kind: str, major: int = 6, type_: str | None = None) -> Path:
    d = db / src_path("bakta_light" if kind == "light" else "bakta_full")
    d.mkdir(parents=True, exist_ok=True)
    (d / "version.json").write_text(_bakta_version(type_ or kind, major))
    return d


def make_genomad(db: Path, version: str = "1.9") -> Path:
    d = db / src_path("genomad")
    d.mkdir(parents=True, exist_ok=True)
    (d / "version.txt").write_text(version + "\n")
    return d


def stray(db: Path) -> list[Path]:
    """Partial downloads and archives left behind."""
    return [p for p in db.rglob("*") if ".partial" in p.name or p.suffix in {".zip", ".gz", ".tmp"}]


def text(capsys: pytest.CaptureFixture[str], caplog: pytest.LogCaptureFixture) -> str:
    return capsys.readouterr().err + "\n" + caplog.text


# ---------------------------------------------------------------- sources and plan


def test_sources_are_the_official_pinned_ones() -> None:
    from urllib.parse import urlparse

    from gippa.annotate.tools import TOOLS
    from gippa.setupdb import SOURCES, DbItem

    assert {"bakta_light", "bakta_full", "genomad", "antismash", "promotech_models"} <= set(SOURCES)
    for name, item in SOURCES.items():
        assert isinstance(item, DbItem) and item.name == name
        assert re.fullmatch(r"(md5:[0-9a-f]{32}|sha256:[0-9a-f]{64})", item.checksum) or (
            name == "antismash" or name.startswith("kegg_")  # M2-09: KEGG is recorded, not pinned
        ), item
        assert item.url.startswith("https://"), item
        assert urlparse(item.url).hostname == item.host, item
        assert item.tool is None or item.tool in TOOLS, item
        assert not Path(item.path).is_absolute()
    assert SOURCES["bakta_light"].checksum == BAKTA_LIGHT_MD5
    assert SOURCES["genomad"].checksum == GENOMAD_MD5
    assert SOURCES["promotech_models"].checksum == PROMOTECH_ZIP_SHA256
    assert SOURCES["bakta_full"].checksum == BAKTA_FULL_MD5
    assert SOURCES["genomad"].path == "genomad_db"
    for name in ("bakta_light", "bakta_full"):
        assert SOURCES[name].host == "zenodo.org" and "14916843" in SOURCES[name].url
        assert SOURCES[name].tool == "bakta" and SOURCES[name].version == "6.0"
    assert SOURCES["bakta_light"].path == "bakta/db-light"
    assert SOURCES["bakta_full"].path == "bakta/db"
    assert SOURCES["genomad"].host == "zenodo.org" and "14886553" in SOURCES["genomad"].url
    assert SOURCES["genomad"].tool is None and SOURCES["genomad"].version == "1.9"  # decision 1
    assert SOURCES["promotech_models"].host == "www.cs.mun.ca"
    assert "/~lourdes/public/PromoTech_models/" in SOURCES["promotech_models"].url
    assert SOURCES["promotech_models"].tool is None
    assert SOURCES["antismash"].host == "dl.secondarymetabolites.org"
    assert SOURCES["antismash"].tool == "antismash"


@pytest.mark.parametrize(
    ("kw", "names"),
    [
        ({"bakta": "light", "genomad": True, "promotech_models": True},
         {"bakta_light", "genomad", "promotech_models"}),
        ({"bakta": "full"}, {"bakta_full"}),
        ({"bakta": "skip"}, set()),
        ({"bakta": "skip", "antismash": True}, {"antismash"}),
    ],
)  # fmt: skip
def test_plan_follows_the_options(tmp_path: Path, kw: dict[str, Any], names: set[str]) -> None:
    from gippa.setupdb import plan

    assert {i.name for i in plan(opts(tmp_path, **kw))} == names


# ---------------------------------------------------------------- verify


def test_verify_nothing_installed_is_absent_without_network(tmp_path: Path) -> None:
    from gippa.setupdb import SOURCES, verify

    for name in ("bakta_light", "bakta_full", "genomad", "promotech_models"):
        rec = verify(SOURCES[name], tmp_path)
        assert rec.name == name and rec.status == "absent", rec


def test_verify_bakta_reads_version_json(tmp_path: Path) -> None:
    from gippa.setupdb import SOURCES, verify

    make_bakta(tmp_path, "light")
    rec = verify(SOURCES["bakta_light"], tmp_path)
    assert rec.status == "ok" and rec.version == "6.0" and rec.path == "bakta/db-light"
    assert verify(SOURCES["bakta_full"], tmp_path).status == "absent"


@pytest.mark.parametrize(("major", "type_"), [(5, "light"), (6, "full")])
def test_verify_bakta_rejects_wrong_major_or_type(tmp_path: Path, major: int, type_: str) -> None:
    from gippa.setupdb import SOURCES, verify

    make_bakta(tmp_path, "light", major=major, type_=type_)
    assert verify(SOURCES["bakta_light"], tmp_path).status != "ok"


def test_verify_bakta_rejects_unreadable_version_json(tmp_path: Path) -> None:
    from gippa.setupdb import SOURCES, verify

    make_bakta(tmp_path, "light").joinpath("version.json").write_text("{not json")
    assert verify(SOURCES["bakta_light"], tmp_path).status != "ok"


@pytest.mark.parametrize(("version", "ok"), [("1.9", True), ("1.9.0", True), ("1.8", False)])
def test_verify_genomad_reads_version_txt(tmp_path: Path, version: str, ok: bool) -> None:
    from gippa.setupdb import SOURCES, verify

    make_genomad(tmp_path, version)
    assert (verify(SOURCES["genomad"], tmp_path).status == "ok") is ok


def test_verify_promotech_rejects_model_with_wrong_sha256(tmp_path: Path) -> None:
    from gippa.setupdb import SOURCES, verify

    d = tmp_path / src_path("promotech_models")
    d.mkdir(parents=True)
    (d / "RF-HOT.model").write_bytes(b"not the pinned model")
    assert verify(SOURCES["promotech_models"], tmp_path).status != "ok"


# ---------------------------------------------------------------- examples


def test_example1_bakta_light_installed(tmp_path: Path) -> None:
    runner = FakeRunner()
    code = setup(tmp_path, runner, bakta="light")
    assert code == ExitCode.SUCCESS == 0
    rec = items(tmp_path)["bakta_light"]
    assert rec["status"] == "ok" and rec["version"] == "6.0"
    assert rec["checksum"] == BAKTA_LIGHT_MD5 and rec["path"] == "bakta/db-light"
    assert "zenodo.org" in rec["url"] and ISO_DATE.fullmatch(rec["retrieved"])
    assert isinstance(rec["size_bytes"], int) and rec["size_bytes"] >= 0
    assert json.loads((tmp_path / "bakta/db-light/version.json").read_text())["type"] == "light"
    assert stray(tmp_path) == []
    [(tool, argv)] = runner.calls
    assert tool == "bakta"
    assert "download" in argv and "--output" in argv
    assert argv[argv.index("--type") + 1] == "light"
    assert any("bakta_db" in a for a in argv)


def test_example2_unreachable_source_exit2_and_nothing_left(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], caplog: pytest.LogCaptureFixture
) -> None:
    code = setup(tmp_path, FakeRunner(exit_code=1), bakta="light")
    assert code == ExitCode.ENVIRONMENT == 2
    assert "bakta_light: source zenodo.org unreachable" in text(capsys, caplog)
    assert not (tmp_path / "bakta/db-light").exists()
    assert stray(tmp_path) == []
    if (tmp_path / "db_manifest.json").exists():
        assert items(tmp_path).get("bakta_light", {}).get("status") != "ok"


def test_example3_check_with_nothing_installed(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], caplog: pytest.LogCaptureFixture
) -> None:
    code = setup(tmp_path, bakta="light", genomad=True, promotech_models=True, check=True)
    assert code == 2
    out = text(capsys, caplog)
    for name in ("bakta_light", "genomad", "promotech_models"):
        assert name in out
    if (tmp_path / "db_manifest.json").exists():
        assert all(i["status"] == "absent" for i in manifest(tmp_path)["items"])


def test_example4_check_after_install(tmp_path: Path) -> None:
    assert setup(tmp_path, FakeRunner(), bakta="light") == 0
    assert setup(tmp_path, bakta="light", check=True) == 0


def test_example5_promotech_checksum_mismatch(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], caplog: pytest.LogCaptureFixture
) -> None:
    fetcher = ZipFetcher(_zip(b"not the published RF-HOT model"))
    code = setup(tmp_path, fetcher=fetcher, bakta="skip", promotech_models=True)
    assert code == 2
    assert fetcher.urls and fetcher.urls[0] == _source_url("promotech_models")
    rec = items(tmp_path)["promotech_models"]
    assert rec["status"] == "mismatch"
    assert not list(tmp_path.rglob("*.zip")), "archive must be deleted"
    assert not list(tmp_path.rglob("RF-HOT.model"))
    assert stray(tmp_path) == []
    assert "promotech_models" in text(capsys, caplog)


def _source_url(name: str) -> str:
    from gippa.setupdb import SOURCES

    return SOURCES[name].url


def test_example6_prefetched_full_db_is_adopted(tmp_path: Path) -> None:
    make_bakta(tmp_path, "full")
    assert setup(tmp_path, NoRunner(), bakta="full") == 0
    rec = items(tmp_path)["bakta_full"]
    assert rec["status"] == "ok" and rec["version"] == "6.0" and rec["path"] == "bakta/db"


def test_example7_nothing_planned_leaves_manifest_unchanged(tmp_path: Path) -> None:
    assert setup(tmp_path, FakeRunner(), bakta="light") == 0
    before = (tmp_path / "db_manifest.json").read_bytes()
    assert setup(tmp_path, bakta="skip") == 0
    assert (tmp_path / "db_manifest.json").read_bytes() == before


def test_example8_two_setups_merge_sorted(tmp_path: Path, genomad_fetcher: ZipFetcher) -> None:
    assert setup(tmp_path, FakeRunner(), bakta="light") == 0
    assert setup(tmp_path, NoRunner(), genomad_fetcher, bakta="skip", genomad=True) == 0
    names = [i["name"] for i in manifest(tmp_path)["items"]]
    assert names == sorted(names) and {"bakta_light", "genomad"} <= set(names)
    assert all(i["status"] == "ok" for i in manifest(tmp_path)["items"])


def test_genomad_is_fetched_from_zenodo_not_by_its_downloader(
    tmp_path: Path, genomad_fetcher: ZipFetcher
) -> None:
    """Decision 1: no ``genomad download-database`` (portal.nersc.gov, unchecked)."""
    assert setup(tmp_path, NoRunner(), genomad_fetcher, bakta="skip", genomad=True) == 0
    assert genomad_fetcher.urls == [_source_url("genomad")]
    assert "zenodo.org/records/14886553/" in genomad_fetcher.urls[0]
    assert items(tmp_path)["genomad"]["version"] == "1.9"
    assert stray(tmp_path) == []


# ---------------------------------------------------------------- idempotence, R6


def test_second_identical_run_is_a_noop(tmp_path: Path, genomad_fetcher: ZipFetcher) -> None:
    assert setup(tmp_path, FakeRunner(), genomad_fetcher, bakta="light", genomad=True) == 0
    before = (tmp_path / "db_manifest.json").read_bytes()
    assert setup(tmp_path, NoRunner(), NoFetcher(), bakta="light", genomad=True) == 0
    assert (tmp_path / "db_manifest.json").read_bytes() == before


def test_same_setup_in_two_dirs_gives_identical_manifests(
    tmp_path: Path, genomad_fetcher: ZipFetcher
) -> None:
    a, b = tmp_path / "a", tmp_path / "b"
    for d in (a, b):
        assert setup(d, FakeRunner(), genomad_fetcher, bakta="light", genomad=True) == 0
    assert (a / "db_manifest.json").read_bytes() == (b / "db_manifest.json").read_bytes()


def test_unreachable_item_blocks_only_itself(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], caplog: pytest.LogCaptureFixture
) -> None:
    """R6: a failing source blocks that item; the others are still installed and recorded."""
    make_bakta(tmp_path, "light")
    code = setup(tmp_path, NoRunner(), DownFetcher(), bakta="light", promotech_models=True)
    assert code == 2
    assert "promotech_models: source www.cs.mun.ca unreachable" in text(capsys, caplog)
    got = items(tmp_path)
    assert got["bakta_light"]["status"] == "ok"
    assert got.get("promotech_models", {}).get("status") != "ok"
    assert stray(tmp_path) == []


def test_failing_item_does_not_stop_later_items(
    tmp_path: Path, genomad_fetcher: ZipFetcher
) -> None:
    runner = FakeRunner(exit_code=1)
    assert setup(tmp_path, runner, genomad_fetcher, bakta="light", genomad=True) == 2
    assert [tool for tool, _ in runner.calls] == ["bakta"]
    assert genomad_fetcher.urls == [_source_url("genomad")]
    assert items(tmp_path)["genomad"]["status"] == "ok"


def test_install_skips_verified_item(tmp_path: Path) -> None:
    from gippa.setupdb import SOURCES, install

    make_genomad(tmp_path)
    rec = install(SOURCES["genomad"], tmp_path, NoRunner(), NoFetcher())
    assert rec.status == "ok" and rec.version == "1.9"


def test_install_returns_absent_with_host_on_failure(tmp_path: Path) -> None:
    from gippa.setupdb import SOURCES, install

    rec = install(SOURCES["genomad"], tmp_path, NoRunner(), DownFetcher())
    assert rec.status == "absent" and "zenodo.org" in rec.message
    assert not (tmp_path / src_path("genomad")).exists() and stray(tmp_path) == []


# ---------------------------------------------------------------- prefetch adoption


PREFETCH = (
    "item\tsource\tfile\tchecksum\tverified\tdate\n"
    "genomad_db\thttps://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz?download=1\t"
    f"genomad_db_v1.9.tar.gz\t{GENOMAD_MD5}\tzenodo md5\t2025-03-02T01:08:26Z\n"
    "genome_k12_mg1655\thttps://ftp.ncbi.nlm.nih.gov/x.fna.gz\tx.fna.gz\tmd5:"
    "7e69874199f23fd21b060dc0b2b72321\tncbi md5checksums.txt\t2025-03-02T01:09:52Z\n"
)


def test_prefetched_item_adopted_with_its_logged_provenance(tmp_path: Path) -> None:
    (tmp_path / "prefetch_manifest.tsv").write_text(PREFETCH)
    make_genomad(tmp_path)
    assert setup(tmp_path, NoRunner(), NoFetcher(), bakta="skip", genomad=True) == 0
    rec = items(tmp_path)["genomad"]
    assert rec["status"] == "ok" and rec["version"] == "1.9"
    assert rec["checksum"] == GENOMAD_MD5
    assert rec["retrieved"] == "2025-03-02", "retrieval date comes from the prefetch log"
    assert (tmp_path / "prefetch_manifest.tsv").read_text() == PREFETCH, "prefetch log untouched"


# ---------------------------------------------------------------- edge cases


def test_low_disk_exits_2_before_downloading(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    caplog: pytest.LogCaptureFixture,
) -> None:
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(100 * GIB, 90 * GIB, 10 * GIB))
    assert setup(tmp_path, NoRunner(), bakta="full") == 2
    out = text(capsys, caplog)
    assert "free" in out.lower() and "10" in out


def test_low_disk_does_not_block_adoption(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(100 * GIB, 99 * GIB, 1 * GIB))
    make_bakta(tmp_path, "full")
    assert setup(tmp_path, NoRunner(), bakta="full") == 0


def test_full_outside_slurm_warns_but_proceeds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    caplog: pytest.LogCaptureFixture,
) -> None:
    for var in ("SLURM_JOB_ID", "SLURM_JOBID"):
        monkeypatch.delenv(var, raising=False)
    runner = FakeRunner()
    assert setup(tmp_path, runner, bakta="full") == 0
    assert runner.calls and "slurm" in text(capsys, caplog).lower()


def test_concurrent_setup_on_same_dir_exits_2(tmp_path: Path) -> None:
    make_genomad(tmp_path)
    inner: list[int] = []
    runner = FakeRunner(
        hook=lambda: inner.append(setup(tmp_path, NoRunner(), bakta="skip", genomad=True))
    )
    assert setup(tmp_path, runner, bakta="light") == 0
    assert inner == [2]
    # the lock is released afterwards
    assert setup(tmp_path, NoRunner(), bakta="skip", genomad=True) == 0
    assert items(tmp_path)["genomad"]["status"] == "ok"


def test_check_never_writes_item_files(tmp_path: Path) -> None:
    setup(tmp_path, bakta="full", genomad=True, promotech_models=True, check=True)
    assert {p.name for p in tmp_path.iterdir()} <= {"db_manifest.json", ".giap_setup.lock"}


def test_adr0014_pins_live_in_db_sources_yaml() -> None:
    import yaml

    raw = (Path(__file__).resolve().parents[2] / "src/gippa/data/db_sources.yaml").read_text()
    for pin in (BAKTA_FULL_MD5, GENOMAD_MD5, PROMOTECH_ZIP_SHA256):
        assert pin.split(":", 1)[1] in raw, pin
    assert isinstance(yaml.safe_load(raw), dict)


# ---------------------------------------------------------------- amrfinderplus row (Q2)


AMR_VERSION = "2026-08-07.1"


def make_amrfinder(db: Path, version: str = AMR_VERSION) -> None:
    d = db / src_path("bakta_full") / "amrfinderplus-db" / "latest"
    d.mkdir(parents=True, exist_ok=True)
    (d / "version.txt").write_text(version + "\n")


def test_amrfinderplus_row_from_bakta_full_db(tmp_path: Path) -> None:
    make_bakta(tmp_path, "full")
    make_amrfinder(tmp_path)
    assert setup(tmp_path, NoRunner(), bakta="full") == 0
    got = items(tmp_path)
    assert got["bakta_full"]["status"] == "ok"
    amr = got["amrfinderplus"]
    assert amr["status"] == "ok" and amr["version"] == AMR_VERSION
    assert amr["path"].startswith("bakta/db/amrfinderplus-db")
    names = [i["name"] for i in manifest(tmp_path)["items"]]
    assert names == sorted(names)


def test_amrfinderplus_row_absent_without_version_file(tmp_path: Path) -> None:
    make_bakta(tmp_path, "full")
    setup(tmp_path, NoRunner(), bakta="full")
    got = items(tmp_path)
    assert got["bakta_full"]["status"] == "ok"
    assert got.get("amrfinderplus", {"status": "absent"})["status"] == "absent"


def test_amrfinderplus_version_follows_the_file(tmp_path: Path) -> None:
    make_bakta(tmp_path, "full")
    make_amrfinder(tmp_path, "2025-01-01.2")
    assert setup(tmp_path, NoRunner(), bakta="full") == 0
    assert items(tmp_path)["amrfinderplus"]["version"] == "2025-01-01.2"


def test_amrfinderplus_row_after_full_download(tmp_path: Path) -> None:
    """The fake ``bakta_db download`` (which runs amrfinder_update) leaves the AMR version file."""

    class AmrRunner(FakeRunner):
        def run(self, tool: str, argv: Sequence[str]) -> int:
            code = super().run(tool, argv)
            out = _out_dir(argv)
            for d in (out, out / "db"):
                latest = d / "amrfinderplus-db" / "latest"
                latest.mkdir(parents=True, exist_ok=True)
                (latest / "version.txt").write_text(AMR_VERSION + "\n")
            return code

    assert setup(tmp_path, AmrRunner(), bakta="full") == 0
    assert items(tmp_path)["amrfinderplus"]["version"] == AMR_VERSION
