"""``gippa setup`` staging when TMPDIR and the DB dir sit on different filesystems (backlog
2026-10-27, "NFS TMPDIR partial-dir tests missing").

On FarmShare ``/tmp`` (TMPDIR) is node-local and ``$GIPPA_ROOT`` (the DB dir) is NFS: an
``os.rename``/``os.replace`` between them raises ``OSError(EXDEV)``, while ``shutil.move`` copies.
A simulator patches ``os.rename``/``os.replace`` so that any rename whose source and destination
lie on different "devices" (the ``local`` tree is TMPDIR, everything else is NFS) raises EXDEV.

What the code does (read from ``setupdb/run.py`` and ``setupdb/uniop_build.py``): every staging
path is a sibling of its target (``<path>.partial`` files and dirs, the ``local_tools.config`` temp
file made with ``mkstemp(dir=db_dir)``), so nothing is staged in TMPDIR and the final rename never
crosses a device. These tests pin that: setup ends with the item installed (claim (a)), with no
``.partial`` or temp leftovers in the DB dir or in TMPDIR (claim (b)), and an interrupted copy
leaves the previous version untouched (claim (c)).

When the final rename itself fails (forced EXDEV: the target is a mount point of another
filesystem; on NFS the same path is taken by ``ESTALE``/``EIO``), M2-08 says "on failure delete
the partial and return ``absent``/``mismatch``";
``test_run_setup_rename_failure_fails_the_item_cleanly`` asks for that.

No test touches the network, the real ``$HOME`` or ``$GIPPA_DB_DIR``.
"""

from __future__ import annotations

import contextlib
import errno
import hashlib
import io
import json
import os
import shutil
import tarfile
import tempfile
from collections.abc import Callable, Iterator, Sequence
from pathlib import Path
from typing import Any

import pytest

# ---------------------------------------------------------------- the two-device simulator


class TwoDevices:
    """``local`` (TMPDIR) and NFS (everything else). A rename across them raises EXDEV."""

    def __init__(self, local: Path, nfs: Path) -> None:
        self.local = local.resolve()
        self.nfs = nfs.resolve()
        self.renames: list[tuple[str, str]] = []
        self.forced: Callable[[Path, Path], bool] = lambda _s, _d: False

    def device(self, p: str | os.PathLike[str]) -> int:
        path = Path(os.path.abspath(os.fspath(p)))
        with contextlib.suppress(OSError):
            path = path.parent.resolve() / path.name
        return 1 if path.is_relative_to(self.local) else 2

    def check(self, src: Any, dst: Any) -> None:
        s, d = Path(os.fspath(src)), Path(os.fspath(dst))
        self.renames.append((str(s), str(d)))
        if self.device(s) != self.device(d) or self.forced(s, d):
            raise OSError(errno.EXDEV, os.strerror(errno.EXDEV), str(s), str(d))


@pytest.fixture
def devices(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[TwoDevices]:
    local, nfs = tmp_path / "local_tmp", tmp_path / "nfs"
    local.mkdir()
    nfs.mkdir()
    sim = TwoDevices(local, nfs)
    real_rename, real_replace = os.rename, os.replace

    def rename(src: Any, dst: Any, *a: Any, **kw: Any) -> None:
        sim.check(src, dst)
        real_rename(src, dst, *a, **kw)

    def replace(src: Any, dst: Any, *a: Any, **kw: Any) -> None:
        sim.check(src, dst)
        real_replace(src, dst, *a, **kw)

    monkeypatch.setattr(os, "rename", rename)
    monkeypatch.setattr(os, "replace", replace)
    monkeypatch.setenv("TMPDIR", str(local))
    monkeypatch.setattr(tempfile, "tempdir", None)  # re-read TMPDIR
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "MAMBA_ROOT_PREFIX"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    yield sim


def leftovers(*roots: Path) -> list[Path]:
    """Staging files or dirs left behind: ``.partial``, ``.tmp`` and mkstemp names."""
    return [
        p
        for r in roots
        for p in r.rglob("*")
        if ".partial" in p.name or p.name.endswith(".tmp") or p.name.startswith("tmp")
    ]


def nothing_in(d: Path) -> bool:
    return not any(d.iterdir())


# ---------------------------------------------------------------- fakes


def sha256(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


class BytesFetcher:
    """Writes ``payload`` to ``dest``; ``crash_after`` bytes then KeyboardInterrupt (a kill)."""

    def __init__(self, payload: bytes, crash_after: int | None = None) -> None:
        self.payload = payload
        self.crash_after = crash_after
        self.dests: list[Path] = []

    def fetch(self, url: str, dest: Path) -> None:
        self.dests.append(dest)
        with dest.open("wb") as fh:
            if self.crash_after is None:
                fh.write(self.payload)
                return
            fh.write(self.payload[: self.crash_after])
        raise KeyboardInterrupt("killed mid-copy")


class NoFetcher:
    def fetch(self, url: str, dest: Path) -> None:
        raise AssertionError(f"fetcher must not be called: {url}")


class NoRunner:
    def run(self, tool: str, argv: Sequence[str]) -> int:
        raise AssertionError(f"runner must not be called: {tool} {argv}")


class BaktaRunner:
    """Fakes ``bakta_db download --output DIR --type light`` (writes ``DIR/db-light``)."""

    def __init__(self, crash: bool = False) -> None:
        self.outs: list[Path] = []
        self.crash = crash

    def run(self, tool: str, argv: Sequence[str]) -> int:
        out = Path(argv[list(argv).index("--output") + 1])
        self.outs.append(out)
        sub = out / "db-light"
        sub.mkdir(parents=True, exist_ok=True)
        (sub / "half.bin").write_bytes(b"x" * 10)
        if self.crash:
            raise KeyboardInterrupt("killed mid-download")
        (sub / "version.json").write_text(json.dumps({"major": 6, "minor": 0, "type": "light"}))
        return 0


def _tar(members: Sequence[tuple[str, bytes]]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, data in members:
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


GENOMAD_TAR = _tar([("genomad_db/version.txt", b"1.9\n"), ("genomad_db/genomad_db", b"mmseqs")])
FILE_PAYLOAD = b">a\nACGT\n>b\nGGCC\n"


def file_item(**kw: Any) -> Any:
    from gippa.setupdb.items import DbItem

    base: dict[str, Any] = {
        "name": "contaminants",
        "version": "1",
        "url": "https://ftp.ncbi.nlm.nih.gov/x/contaminants.fa",
        "host": "ftp.ncbi.nlm.nih.gov",
        "checksum": sha256(FILE_PAYLOAD),
        "path": "contaminants/contaminants.fa",
        "tool": None,
        "check": "fasta:2",
    }
    return DbItem(**(base | kw))


def genomad_item() -> Any:
    from gippa.setupdb.items import DbItem

    return DbItem("genomad", "1.9", "https://zenodo.org/records/1/files/genomad_db_v1.9.tar.gz",
                  "zenodo.org", sha256(GENOMAD_TAR), "genomad_db", None, "genomad")  # fmt: skip


def bakta_item() -> Any:
    from gippa.setupdb.items import DbItem

    return DbItem("bakta_light", "6.0", "https://zenodo.org/records/2/files/db-light.tar.xz",
                  "zenodo.org", "", "bakta/db-light", "bakta", "bakta:light")  # fmt: skip


def install(item: Any, db: Path, runner: Any = None, fetcher: Any = None) -> Any:
    from gippa.setupdb.run import install as _install

    return _install(item, db, runner or NoRunner(), fetcher or NoFetcher())


# ---------------------------------------------------------------- the simulator itself


def test_simulator_raises_exdev_across_devices_and_shutil_move_copies(
    devices: TwoDevices,
) -> None:
    """Guards the tests below: a cross-device rename fails, a same-device one does not."""
    src = devices.local / "a.partial"
    src.write_bytes(b"x")
    with pytest.raises(OSError) as exc:
        os.replace(src, devices.nfs / "a")
    assert exc.value.errno == errno.EXDEV
    with pytest.raises(OSError):
        Path(src).rename(devices.nfs / "a")
    shutil.move(str(src), str(devices.nfs / "a"))  # falls back to copy + unlink
    assert (devices.nfs / "a").read_bytes() == b"x" and not src.exists()
    (devices.nfs / "a").replace(devices.nfs / "b")
    assert tempfile.gettempdir() == str(devices.local)


# ------------------------------------------------ (a) + (b): installs across devices


def test_file_item_installs_with_tmpdir_on_another_device(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    fetcher = BytesFetcher(FILE_PAYLOAD)
    rec = install(file_item(), db, fetcher=fetcher)
    target = db / "contaminants" / "contaminants.fa"
    assert rec.status == "ok", rec.message
    assert target.read_bytes() == FILE_PAYLOAD
    assert [d.parent for d in fetcher.dests] == [target.parent]  # staged next to the target
    assert leftovers(db, devices.local) == []
    assert nothing_in(devices.local)


def test_archive_dir_item_installs_with_tmpdir_on_another_device(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    fetcher = BytesFetcher(GENOMAD_TAR)
    rec = install(genomad_item(), db, fetcher=fetcher)
    assert rec.status == "ok", rec.message
    assert (db / "genomad_db" / "version.txt").read_text() == "1.9\n"
    assert all(devices.device(d) == devices.device(db) for d in fetcher.dests)
    assert leftovers(db, devices.local) == []
    assert nothing_in(devices.local)


def test_downloader_dir_item_installs_with_tmpdir_on_another_device(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    runner = BaktaRunner()
    rec = install(bakta_item(), db, runner=runner)
    assert rec.status == "ok", rec.message
    assert json.loads((db / "bakta" / "db-light" / "version.json").read_text())["major"] == 6
    assert [o.parent for o in runner.outs] == [db / "bakta"]  # the downloader writes on NFS
    assert leftovers(db, devices.local) == []
    assert nothing_in(devices.local)


def test_no_rename_in_setup_crosses_a_device(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    install(file_item(), db, fetcher=BytesFetcher(FILE_PAYLOAD))
    install(genomad_item(), db, fetcher=BytesFetcher(GENOMAD_TAR))
    install(bakta_item(), db, runner=BaktaRunner())
    assert devices.renames, "the installs rename their staging copies into place"
    assert all(devices.device(s) == devices.device(d) for s, d in devices.renames)


def test_db_dir_symlinked_onto_nfs_installs(devices: TwoDevices) -> None:
    """``--db-dir`` given as a link from the local disk onto NFS (``~/db -> $GIPPA_ROOT/db``)."""
    real = devices.nfs / "db"
    real.mkdir()
    link = devices.local.parent / "home_db"
    link.symlink_to(real)
    rec = install(file_item(), link, fetcher=BytesFetcher(FILE_PAYLOAD))
    assert rec.status == "ok", rec.message
    assert (real / "contaminants" / "contaminants.fa").read_bytes() == FILE_PAYLOAD
    assert leftovers(real, devices.local) == []


def test_run_setup_installs_with_tmpdir_on_another_device(
    devices: TwoDevices, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.exitcodes import ExitCode
    from gippa.setupdb import run as run_mod

    db = devices.nfs / "db"
    item = file_item()
    monkeypatch.setattr(run_mod, "plan", lambda _o: [item])
    opts = _opts(db)
    code = run_mod.run_setup(opts, NoRunner(), BytesFetcher(FILE_PAYLOAD))
    assert code == ExitCode.SUCCESS
    recs = {i["name"]: i for i in json.loads((db / "db_manifest.json").read_text())["items"]}
    assert recs["contaminants"]["status"] == "ok"
    assert leftovers(db, devices.local) == []
    assert nothing_in(devices.local)


def _opts(db: Path) -> Any:
    from gippa.setupdb import SetupOptions

    return SetupOptions(db_dir=db, bakta="skip", genomad=False, antismash=False,
                        promotech_models=False, check=False, go=False, kegg=False, cog=False,
                        contaminants=False)  # fmt: skip


# ---------------------------------------------------------------- (b): leftovers of a killed run


def test_stale_partial_file_from_a_killed_run_is_replaced(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    target = db / "contaminants" / "contaminants.fa"
    target.parent.mkdir(parents=True)
    target.with_name(target.name + ".partial").write_bytes(b"garbage from a SIGKILL" * 100)
    rec = install(file_item(), db, fetcher=BytesFetcher(FILE_PAYLOAD))
    assert rec.status == "ok", rec.message
    assert target.read_bytes() == FILE_PAYLOAD
    assert leftovers(db, devices.local) == []


def test_stale_partial_dir_from_a_killed_run_is_cleared(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    stale = db / "genomad_db.partial"
    (stale / "genomad_db").mkdir(parents=True)
    (stale / "genomad_db" / "junk").write_bytes(b"half")
    (stale / "genomad_db" / "version.txt").write_text("0.1\n")
    rec = install(genomad_item(), db, fetcher=BytesFetcher(GENOMAD_TAR))
    assert rec.status == "ok", rec.message
    assert not (db / "genomad_db" / "junk").exists()
    assert leftovers(db, devices.local) == []


# ---------------------------------------------------------------- (c): a crash mid-copy


def test_crash_mid_copy_of_a_file_leaves_the_old_version(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    target = db / "contaminants" / "contaminants.fa"
    target.parent.mkdir(parents=True)
    target.write_bytes(b">old\nAAAA\n")  # an older release: does not verify, so it is re-fetched
    with pytest.raises(KeyboardInterrupt):
        install(file_item(), db, fetcher=BytesFetcher(FILE_PAYLOAD, crash_after=5))
    assert target.read_bytes() == b">old\nAAAA\n"
    assert leftovers(db, devices.local) == []
    assert nothing_in(devices.local)


def test_crash_mid_copy_of_an_archive_leaves_the_old_dir(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    old = db / "genomad_db"
    old.mkdir(parents=True)
    (old / "version.txt").write_text("1.8\n")  # older: not adopted, so re-fetched
    (old / "genomad_db").write_bytes(b"old mmseqs")
    with pytest.raises(KeyboardInterrupt):
        install(genomad_item(), db, fetcher=BytesFetcher(GENOMAD_TAR, crash_after=20))
    assert (old / "version.txt").read_text() == "1.8\n"
    assert (old / "genomad_db").read_bytes() == b"old mmseqs"
    assert leftovers(db, devices.local) == []
    assert nothing_in(devices.local)


def test_crash_mid_download_by_a_tool_leaves_the_old_dir(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    old = db / "bakta" / "db-light"
    old.mkdir(parents=True)
    (old / "version.json").write_text(json.dumps({"major": 5, "minor": 1, "type": "light"}))
    with pytest.raises(KeyboardInterrupt):
        install(bakta_item(), db, runner=BaktaRunner(crash=True))
    assert json.loads((old / "version.json").read_text())["major"] == 5
    assert sorted(p.name for p in old.iterdir()) == ["version.json"]
    assert leftovers(db, devices.local) == []


# ---------------------------------------------------------------- the final rename fails


def test_failed_final_rename_of_a_file_keeps_the_old_version_and_no_partial(
    devices: TwoDevices,
) -> None:
    db = devices.nfs / "db"
    target = db / "contaminants" / "contaminants.fa"
    target.parent.mkdir(parents=True)
    target.write_bytes(b">old\nAAAA\n")
    devices.forced = lambda s, d: d == target
    with pytest.raises(OSError):  # install itself may raise; run_setup must not (below)
        install(file_item(), db, fetcher=BytesFetcher(FILE_PAYLOAD))
    assert target.read_bytes() == b">old\nAAAA\n"
    assert leftovers(db, devices.local) == []


def test_failed_final_rename_of_a_dir_leaves_no_partial(devices: TwoDevices) -> None:
    db = devices.nfs / "db"
    devices.forced = lambda s, d: d == db / "genomad_db"
    with pytest.raises(OSError):
        install(genomad_item(), db, fetcher=BytesFetcher(GENOMAD_TAR))
    assert leftovers(db, devices.local) == []
    assert nothing_in(devices.local)


def test_run_setup_rename_failure_fails_the_item_cleanly(
    devices: TwoDevices, monkeypatch: pytest.MonkeyPatch
) -> None:
    """M2-08: "on failure delete the partial and return absent/mismatch with the host"; §7.9.1:
    a failed item is exit 2 and the later items still install. An ``OSError`` at the rename
    (EXDEV here; ESTALE/EIO on NFS) must not escape ``run_setup`` as exit 5."""
    from gippa.exitcodes import ExitCode
    from gippa.setupdb import run as run_mod

    db = devices.nfs / "db"
    first, second = file_item(), genomad_item()
    target = db / first.path
    monkeypatch.setattr(run_mod, "plan", lambda _o: [first, second])
    devices.forced = lambda s, d: d == target

    class Both:
        def fetch(self, url: str, dest: Path) -> None:
            dest.write_bytes(GENOMAD_TAR if url == second.url else FILE_PAYLOAD)

    code = run_mod.run_setup(_opts(db), NoRunner(), Both())
    assert code == ExitCode.ENVIRONMENT
    recs = {i["name"]: i for i in json.loads((db / "db_manifest.json").read_text())["items"]}
    assert recs["contaminants"]["status"] in ("absent", "failed")
    message = recs["contaminants"]["message"]
    assert first.path in message or str(target) in message
    assert recs["genomad"]["status"] == "ok"
    assert leftovers(db, devices.local) == []


# ---------------------------------------------------------------- uniop_build


def uniop_item() -> Any:
    payload = _tar([("UniOP-3ed14dfc/UniOP.py", b"print(1)\n")])
    from gippa.setupdb.items import DbItem

    item = DbItem("uniop_source", "git:3ed14dfc", "https://github.com/hongsua/UniOP/x.tar.gz",
                  "github.com", sha256(payload), "uniop/UniOP-3ed14dfc.tar.gz", None,
                  "sha256")  # fmt: skip
    return item, payload


def test_uniop_archive_fetch_across_devices(devices: TwoDevices) -> None:
    from gippa.setupdb.uniop_build import fetch_archive

    db = devices.nfs / "db"
    item, payload = uniop_item()
    fetcher = BytesFetcher(payload)
    rec = fetch_archive(item, db, fetcher)
    assert rec.status == "ok", rec.message
    assert (db / item.path).read_bytes() == payload
    assert [d.parent for d in fetcher.dests] == [(db / item.path).parent]
    assert leftovers(db, devices.local) == []
    assert nothing_in(devices.local)


def test_uniop_archive_crash_mid_copy_leaves_no_partial(devices: TwoDevices) -> None:
    from gippa.setupdb.uniop_build import fetch_archive

    db = devices.nfs / "db"
    item, payload = uniop_item()
    with pytest.raises(KeyboardInterrupt):
        fetch_archive(item, db, BytesFetcher(payload, crash_after=7))
    assert not (db / item.path).exists()
    assert leftovers(db, devices.local) == []


def test_uniop_verified_archive_is_kept_without_a_fetch(devices: TwoDevices) -> None:
    from gippa.setupdb.uniop_build import fetch_archive

    db = devices.nfs / "db"
    item, payload = uniop_item()
    (db / item.path).parent.mkdir(parents=True)
    (db / item.path).write_bytes(payload)
    assert fetch_archive(item, db, NoFetcher()).status == "ok"
    assert (db / item.path).read_bytes() == payload


def test_uniop_failed_archive_rename_leaves_no_partial(devices: TwoDevices) -> None:
    from gippa.setupdb.uniop_build import fetch_archive

    db = devices.nfs / "db"
    item, payload = uniop_item()
    devices.forced = lambda s, d: d == db / item.path
    with pytest.raises(OSError):  # setup_uniop turns an OSError into an absent record
        fetch_archive(item, db, BytesFetcher(payload))
    assert leftovers(db, devices.local) == []


def test_local_tools_config_is_written_across_devices(devices: TwoDevices) -> None:
    from gippa.setupdb.uniop_build import _write_local_tools

    db = devices.nfs / "db"
    db.mkdir()
    _write_local_tools(db, "docker", "gippa-local/uniop:abc")
    assert "gippa-local/uniop:abc" in (db / "local_tools.config").read_text()
    assert leftovers(db, devices.local) == []
    assert nothing_in(devices.local)
    assert all(devices.device(s) == devices.device(d) for s, d in devices.renames)


@pytest.mark.parametrize("how", ["crash", "rename"])
def test_local_tools_config_interrupted_keeps_the_old_file(
    devices: TwoDevices, monkeypatch: pytest.MonkeyPatch, how: str
) -> None:
    from gippa.setupdb import uniop_build
    from gippa.setupdb.uniop_build import _write_local_tools

    db = devices.nfs / "db"
    db.mkdir()
    _write_local_tools(db, "docker", "gippa-local/uniop:old")
    before = (db / "local_tools.config").read_bytes()
    if how == "crash":
        real = uniop_build.local_tools_text

        def boom(*a: Any, **kw: Any) -> str:
            real(*a, **kw)
            raise KeyboardInterrupt("killed mid-write")

        monkeypatch.setattr(uniop_build, "local_tools_text", boom)
        expected: type[BaseException] = KeyboardInterrupt
    else:
        devices.forced = lambda s, d: d == db / "local_tools.config"
        expected = OSError
    with pytest.raises(expected):
        _write_local_tools(db, "docker", "gippa-local/uniop:new")
    assert (db / "local_tools.config").read_bytes() == before
    assert [p for p in db.iterdir() if p.name != "local_tools.config"] == []
    assert nothing_in(devices.local)
