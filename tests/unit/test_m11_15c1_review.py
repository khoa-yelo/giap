"""M11-15c1 review round 1 (PR #440): the test side of B1, B3-B7, the prompt order and three nits.

docs/specs/M11-15-portable-clusters.md item 6 ("Before any download (15c1)", "Quotas (15c1)",
"``--db-dir`` becomes optional (15c1)"), examples 6.4, 6.8, 6.10; the review of 7f3e896e.

Readings and lead decisions pinned here:

- B1: ``tests/conftest.py`` isolates every test (``_isolated_user_dirs``): ``XDG_CONFIG_HOME``
  and ``XDG_DATA_HOME`` in a per-test temporary directory, ``space.MOUNTINFO`` at
  ``tests/data/mountinfo/root_ext4.txt``, ``space.user_quota`` answering ``None`` unless the
  test faked ``subprocess.run``. The guard tests below fail if a setup run writes under the home
  config or runs a real ``quota``.
- B3: a downloader run through the real ``ToolRunner`` (here ``download-antismash-databases`` in
  ``gippa-antismash``) that exits nonzero after "No space left on device" / "Disk quota exceeded"
  on its stderr is ``failed`` with the ``--db-dir`` advice, not ``absent`` "source unreachable".
  The fake ``subprocess.run`` honours ``capture_output``/``stderr=PIPE``/``stderr=<file>``, so
  the implementation may capture or tee the stream; inherited stderr gets it on fd 2.
- B4 (lead decision): setup keeps reading the pointer, but the home-quota warning (example 6.10)
  appears whenever the resolved DB directory is under ``$HOME`` on a shared filesystem and
  neither ``--db-dir`` nor ``GIPPA_DB_DIR`` was given, whatever the source (pointer included).
- B6: the quota applies per place: the image cache on another quota'd filesystem is checked
  against that filesystem's row. ``quota -w`` runs once per setup and is filtered per place.
- B7: an item without a recorded size prints "size not recorded" in its plan line, not 0.0 GiB.
- Prompt order (lead decision): ``Proceed? [Y/n]`` first, then the UniOP terms prompt, both on a
  TTY; "n" to Proceed? stops before any download, pull, fetch or build, UniOP path included.
- Nits: an unreadable mount table applies no quota; two roles sharing one path count once; a
  read-only config directory makes ``write_pointer`` warn, not crash.

The B5 tests kill mutants that survived the round-1 tests; they pass on 7f3e896e by design.
No test touches the network (``no_network``; fake ``http.client``), a real tool, the real
``$HOME`` or a real DB directory.
"""

from __future__ import annotations

import contextlib
import errno
import http.client
import io
import json
import os
import pwd
import re
import shutil
import stat
import subprocess
import sys
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any, NamedTuple

import pytest
from click.testing import CliRunner
from helpers.prereqs import run_fake_prereq
from helpers.setup_opt_outs import NO_RESOURCE_FLAGS, opt_outs

from gippa.cli import cli
from gippa.exitcodes import ExitCode
from gippa.setupdb import (
    SOURCES,
    HttpsFetcher,
    SetupOptions,
    SourceError,
    ToolRunner,
    run_setup,
    space,
)
from gippa.setupdb.dbdir import pointer_file

pytestmark = pytest.mark.usefixtures("no_network")

# captured at import (collection), before any fixture patches them
REAL_RUN = subprocess.run
REAL_FETCH = HttpsFetcher.__dict__["fetch"]
REAL_TOOLRUNNER_RUN = ToolRunner.__dict__["run"]

GIB = 1024**3
KIB = 1024
# M11-15c2 turns the resource files and the images on by default; these tests are about genomad
# + promotech_models (and images only where they pass --images), so the rest is opted out
# (opt_outs keeps only the flags the current setup accepts). Without --no-images a real
# subprocess.run would start real image pulls.
_OTHERS = tuple(opt_outs(*NO_RESOURCE_FLAGS))
HTTPS_DBS = ("--bakta", "skip", "--genomad", "--promotech-models", *_OTHERS)
HTTPS_ITEMS = (*HTTPS_DBS, *opt_outs("images"))
NO_DBS = ("--bakta", "skip", "--no-genomad", "--no-promotech-models", *_OTHERS)
PLANNED = ("genomad", "promotech_models")
HOME_WARNING = (
    "the default DB directory is in your home directory, which usually has a small quota on "
    "clusters; pass --db-dir on a project or scratch filesystem"
)
ADVICE = "disk quota exceeded or disk full at"
MARK = "<<FIRST-FETCH-OR-PULL>>"
ENOSPC_TEXT = "OSError: [Errno 28] No space left on device: '/db/antismash.partial/pfam.hmm'"
EDQUOT_TEXT = "error: write failed: Disk quota exceeded"


class _Usage(NamedTuple):
    total: int
    used: int
    free: int


def _done(argv: Sequence[str], rc: int, out: str = "", err: str = "") -> Any:
    return subprocess.CompletedProcess(list(argv), rc, out, err)


class Calls:
    """Fake ``subprocess.run`` for ``quota``/``lfs``/``apptainer pull``; anything else fails."""

    def __init__(self) -> None:
        self.argv: list[list[str]] = []
        self.quota_out = ""
        self.lfs_out = ""

    def __call__(self, argv: Sequence[str], *a: Any, **k: Any) -> Any:
        argv = [str(x) for x in argv]
        if (done := run_fake_prereq(argv, *a, **k)) is not None:  # 15d: java -version
            return done
        self.argv.append(argv)
        if argv[0] == "quota":
            return _done(argv, 0, self.quota_out)
        if argv[0] == "lfs":
            return _done(argv, 0, self.lfs_out)
        if argv[:2] == ["apptainer", "pull"]:
            print(MARK, file=sys.stderr, flush=True)
            return _done(argv, 1, "", "FATAL: fake pull refused")
        raise AssertionError(f"unexpected subprocess.run {argv}")

    def named(self, prog: str) -> list[list[str]]:
        return [a for a in self.argv if a[0] == prog]

    @property
    def pulls(self) -> list[list[str]]:
        return [a for a in self.argv if a[:2] == ["apptainer", "pull"]]


@pytest.fixture(autouse=True)
def _base(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USER", "alice")
    monkeypatch.setenv("LOGNAME", "alice")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(tmp_path / "no_mamba"))
    for var in (
        "GIPPA_DB_DIR",
        "GIPPA_ROOT",
        "GIPPA_REPO",
        "NXF_APPTAINER_CACHEDIR",
        "APPTAINER_CACHEDIR",
        "APPTAINER_TMPDIR",
        "SINGULARITY_CACHEDIR",
        "GIPPA_CONTAINER_RUNTIME",
        "UNIOP_HOME",
    ):
        monkeypatch.delenv(var, raising=False)
    return home


@pytest.fixture
def env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, _base: Path) -> dict[str, Any]:
    """CLI fakes: 10 TB free everywhere, ``subprocess.run`` = ``Calls``, a fetcher that records
    the URL and refuses (so "setup proceeded" is "the fetcher was called")."""
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    calls = Calls()
    monkeypatch.setattr(subprocess, "run", calls)
    fetched: list[str] = []

    def fetch(_self: Any, url: str, dest: Path) -> None:
        print(MARK, file=sys.stderr, flush=True)
        fetched.append(url)
        raise SourceError("fake.example", "refused by the test")

    monkeypatch.setattr(HttpsFetcher, "fetch", fetch)
    return {"home": _base, "calls": calls, "fetched": fetched}


def _setup(*args: str, stdin: str | None = None) -> Any:
    return CliRunner().invoke(cli, ["setup", *args], input=stdin, catch_exceptions=False)


def _plan(output: str) -> str:
    return output.split(MARK, 1)[0]


def _tty(monkeypatch: pytest.MonkeyPatch) -> None:
    from gippa.setupdb import uniop_build

    monkeypatch.setattr(uniop_build, "stdin_is_tty", lambda: True, raising=False)
    monkeypatch.setattr(space, "stdin_is_tty", lambda: True, raising=False)


def _mountinfo(rows: Sequence[tuple[str, str, str]], root: str = "/") -> str:
    """``/proc/self/mountinfo`` lines for (mount point as written, fstype, source)."""
    return "".join(
        f"{i} 1 0:{i} {root} {mnt} rw,relatime shared:{i} - {fstype} {src} rw\n"
        for i, (mnt, fstype, src) in enumerate(rows, start=20)
    )


def _set_mounts(tmp: Path, mp: pytest.MonkeyPatch, rows: Sequence[tuple[str, str, str]]) -> None:
    f = tmp / "mountinfo"
    f.write_text(_mountinfo([("/", "ext4", "/dev/sda1"), *rows]))
    mp.setattr(space, "MOUNTINFO", f)


def _quota_rows(*rows: tuple[str, int, int]) -> str:
    """``quota -w`` output; rows of (filesystem, used KiB, limit KiB)."""
    head = (
        "Disk quotas for user alice (uid 1001): \n"
        "     Filesystem  blocks   quota   limit   grace   files   quota   limit   grace\n"
    )
    return head + "".join(
        f"{fs} {used} {lim} {lim}            10       0       0\n" for fs, used, lim in rows
    )


def _dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path.resolve()


def _items(db: Path) -> dict[str, dict[str, Any]]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    return {i["name"]: i for i in data["items"]}


# ================================================================ B1: isolation guards


def test_conftest_isolates_xdg_homes(tmp_path_factory: pytest.TempPathFactory) -> None:
    base = tmp_path_factory.getbasetemp().resolve()
    real_home = Path(pwd.getpwuid(os.getuid()).pw_dir).resolve()
    for var in ("XDG_CONFIG_HOME", "XDG_DATA_HOME"):
        value = os.environ.get(var)
        assert value, f"tests/conftest.py must set {var} for every test (review B1)"
        where = Path(value).resolve()
        assert where.is_relative_to(base), f"{var}={value} is not a pytest temporary directory"
        assert where != real_home / ".config" and where != real_home / ".local" / "share"


def test_conftest_points_mountinfo_at_a_fixture() -> None:
    assert Path(space.MOUNTINFO) != Path("/proc/self/mountinfo"), "the host's mount table"
    assert Path(space.MOUNTINFO).is_file()
    assert Path(space.MOUNTINFO).resolve().is_relative_to(Path(__file__).resolve().parents[1])


def _isolated_config(factory: pytest.TempPathFactory) -> str:
    """``$XDG_CONFIG_HOME``, asserted to be a pytest temporary directory before setup runs."""
    xdg = os.environ.get("XDG_CONFIG_HOME")
    assert xdg and Path(xdg).resolve().is_relative_to(factory.getbasetemp().resolve()), (
        "XDG_CONFIG_HOME is not isolated: running setup now would write the real home config"
    )
    return xdg


def test_setup_run_never_writes_under_the_home_config(
    tmp_path: Path, tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Guard (B1): a plain ``gippa setup`` writes its pointer under the isolated config home,
    nothing under ``$HOME/.config`` (here a stand-in home; the real one in the next test)."""
    xdg = _isolated_config(tmp_path_factory)
    fake_home = tmp_path / "standin_home"
    fake_home.mkdir()
    monkeypatch.setenv("HOME", str(fake_home))
    monkeypatch.setattr(HttpsFetcher, "fetch", _refuse)
    monkeypatch.setattr(subprocess, "run", Calls())
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    assert not (fake_home / ".config").exists(), r.output
    assert pointer_file(Path(xdg)).is_file(), "setup did not record the DB directory at all"


def test_setup_run_leaves_the_real_home_config_alone(
    tmp_path: Path, tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    _isolated_config(tmp_path_factory)
    real = Path(pwd.getpwuid(os.getuid()).pw_dir) / ".config" / "gippa" / "db_dir"
    before = real.stat().st_mtime_ns if real.exists() else None
    monkeypatch.setenv("HOME", pwd.getpwuid(os.getuid()).pw_dir)
    monkeypatch.setattr(HttpsFetcher, "fetch", _refuse)
    monkeypatch.setattr(subprocess, "run", Calls())
    _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    after = real.stat().st_mtime_ns if real.exists() else None
    assert after == before, f"setup wrote {real}"


def _refuse(_self: Any, url: str, dest: Path) -> None:
    raise SourceError("fake.example", "refused by the test")


QUOTA_SPY = """#!/bin/sh
echo "$0 $*" >> "$QUOTA_SPY_LOG"
exit 0
"""


def test_setup_run_never_starts_a_real_quota(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Guard (B1): with the real ``subprocess.run`` no ``quota``/``lfs`` process starts (a spy
    ``quota`` first on ``PATH`` logs any call; the real one is never reached)."""
    bindir = tmp_path / "spybin"
    bindir.mkdir()
    for prog in ("quota", "lfs"):
        exe = bindir / prog
        exe.write_text(QUOTA_SPY)
        exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    log = tmp_path / "quota_spy.log"
    monkeypatch.setenv("QUOTA_SPY_LOG", str(log))
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setattr(subprocess, "run", REAL_RUN)
    monkeypatch.setattr(HttpsFetcher, "fetch", _refuse)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    assert "GiB" in r.output, r.output  # the plan ran
    assert not log.exists(), f"a quota command ran: {log.read_text()}"


# ================================================================ B3: full disk in a downloader


def _emit(argv: list[str], kw: dict[str, Any], rc: int, err: str) -> Any:
    """Return or write ``err`` the way the caller asked for stderr."""
    text = bool(
        kw.get("text") or kw.get("universal_newlines") or kw.get("encoding") or kw.get("errors")
    )
    payload: str | bytes = err if text else err.encode()
    empty: str | bytes = "" if text else b""
    target = kw.get("stderr")
    if kw.get("capture_output") or target == subprocess.PIPE:
        return subprocess.CompletedProcess(argv, rc, empty, payload)
    if target == subprocess.STDOUT and kw.get("stdout") == subprocess.PIPE:
        return subprocess.CompletedProcess(argv, rc, payload, None)
    if isinstance(target, int) and target >= 0:
        os.write(target, err.encode())
    elif target is not None and hasattr(target, "write"):
        try:
            target.write(payload)
        except TypeError:
            target.write(err if isinstance(payload, bytes) else err.encode())
        if hasattr(target, "flush"):
            target.flush()
    else:
        os.write(2, err.encode())
    return subprocess.CompletedProcess(argv, rc, None, None)


class Downloader:
    """Fake ``subprocess.run`` behind the real ``ToolRunner``: ``antismash --check-prereqs``
    fails (databases not in the image), ``download-antismash-databases`` writes half a file,
    prints ``err`` on stderr and exits 1."""

    def __init__(self, err: str) -> None:
        self.err = err
        self.argv: list[list[str]] = []

    def __call__(self, argv: Sequence[str], *a: Any, **kw: Any) -> Any:
        argv = [str(x) for x in argv]
        if (done := run_fake_prereq(argv, *a, **kw)) is not None:  # 15d: java -version
            return done
        self.argv.append(argv)
        if argv[0] in ("quota", "lfs"):
            return _done(argv, 0)
        if "--check-prereqs" in argv:
            return _emit(argv, kw, 1, "antiSMASH: databases missing\n")
        if "download-antismash-databases" in argv:
            out = Path(argv[argv.index("--database-dir") + 1])
            out.mkdir(parents=True, exist_ok=True)
            (out / "pfam.hmm").write_bytes(b"half")
            return _emit(argv, kw, 1, f"Downloading Pfam-A.hmm.gz\n{self.err}\n")
        raise AssertionError(f"unexpected subprocess.run {argv}")


@pytest.mark.parametrize("err", [ENOSPC_TEXT, EDQUOT_TEXT], ids=["ENOSPC", "EDQUOT"])
def test_b3_full_disk_inside_the_downloader_is_failed_with_the_advice(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, err: str
) -> None:
    (tmp_path / "mamba" / "envs" / "gippa-antismash").mkdir(parents=True)
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(tmp_path / "mamba"))
    monkeypatch.setattr(ToolRunner, "run", REAL_TOOLRUNNER_RUN)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    fake = Downloader(err)
    monkeypatch.setattr(subprocess, "run", fake)
    monkeypatch.setattr(HttpsFetcher, "fetch", _refuse)
    db = tmp_path / "db"
    opts = SetupOptions(db, bakta="skip", genomad=False, antismash=True, promotech_models=False)
    got = run_setup(opts, ToolRunner(db_dir=db), HttpsFetcher())
    assert any("download-antismash-databases" in a for a in fake.argv), fake.argv
    assert got == ExitCode.ENVIRONMENT
    item = _items(db)["antismash"]
    assert item["status"] == "failed", item
    assert ADVICE in item["message"] and "pass --db-dir" in item["message"], item
    assert str(db) in item["message"], item
    assert "unreachable" not in item["message"], "a full disk is not the source's fault"
    assert not (db / SOURCES["antismash"].path).exists()
    assert [p for p in db.rglob("*") if ".part" in p.name] == []


def test_b3_other_downloader_failures_stay_absent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "mamba" / "envs" / "gippa-antismash").mkdir(parents=True)
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(tmp_path / "mamba"))
    monkeypatch.setattr(ToolRunner, "run", REAL_TOOLRUNNER_RUN)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    monkeypatch.setattr(subprocess, "run", Downloader("ConnectionError: Max retries exceeded"))
    db = tmp_path / "db"
    opts = SetupOptions(db, bakta="skip", genomad=False, antismash=True, promotech_models=False)
    run_setup(opts, ToolRunner(db_dir=db), HttpsFetcher())
    item = _items(db)["antismash"]
    assert item["status"] == "absent", item
    assert ADVICE not in item["message"]


# ================================================================ B4: home warning, any source


def _home_on_nfs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, home: Path) -> None:
    _set_mounts(tmp_path, monkeypatch, [(str(home.resolve()), "nfs4", "srv:/export/home")])


def _write_pointer(target: Path) -> None:
    ptr = pointer_file(Path(os.environ["XDG_CONFIG_HOME"]))
    ptr.parent.mkdir(parents=True, exist_ok=True)
    ptr.write_text(f"{target}\n", encoding="utf-8")


@pytest.mark.parametrize("where", ["default", "other"])
def test_b4_pointer_into_shared_home_still_warns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any], where: str
) -> None:
    """A second plain ``gippa setup`` (pointer recorded) in a shared home warns as the first."""
    home = env["home"]
    monkeypatch.delenv("XDG_DATA_HOME")  # the data dir is ~/.local/share/gippa
    _home_on_nfs(tmp_path, monkeypatch, home)
    target = home / ".local" / "share" / "gippa" / "db" if where == "default" else home / "mydb"
    target.mkdir(parents=True)
    _write_pointer(target)
    r = _setup(*HTTPS_ITEMS)
    plan = _plan(r.output)
    assert HOME_WARNING in " ".join(plan.split()), r.output
    assert env["fetched"], "without a TTY setup proceeds"


def test_b4_pointer_outside_home_does_not_warn(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    _home_on_nfs(tmp_path, monkeypatch, env["home"])
    target = _dir(tmp_path / "scratch" / "db")
    _write_pointer(target)
    r = _setup(*HTTPS_ITEMS)
    assert "small quota" not in r.output
    assert env["fetched"]


# ================================================================ B5: surviving mutants


def test_b5_check_writes_no_pointer(tmp_path: Path, env: dict[str, Any]) -> None:
    ptr = pointer_file(Path(os.environ["XDG_CONFIG_HOME"]))
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS, "--check")
    assert not ptr.exists(), f"--check recorded the DB directory:\n{r.output}"
    _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    assert ptr.is_file(), "control: a plain setup records it"


# ---------------------------------------------------------------- real HttpsFetcher re-raises


class _Resp:
    status = 200

    def __init__(self, body: bytes) -> None:
        self._body = io.BytesIO(body)
        self._len = len(body)

    def getheader(self, name: str, default: str | None = None) -> str | None:
        return str(self._len) if name.lower() == "content-length" else default

    def read(self, n: int = -1) -> bytes:
        return self._body.read(n)


class _Conn:
    """Fake ``http.client.HTTPSConnection``: 200 with 64 KiB; never opens a socket."""

    def __init__(self, host: str, port: int | None = None, **_kw: Any) -> None:
        self.host = host

    def request(self, *_a: Any, **_k: Any) -> None: ...

    def getresponse(self) -> _Resp:
        return _Resp(b"x" * 65536)

    def close(self) -> None: ...


class _FullFile:
    def __init__(self, real: Any, code: int) -> None:
        self.real, self.code = real, code

    def write(self, _b: bytes) -> int:
        raise OSError(self.code, os.strerror(self.code))

    def __enter__(self) -> _FullFile:
        return self

    def __exit__(self, *_a: Any) -> None:
        self.real.close()

    def close(self) -> None:
        self.real.close()


def _full_disk_under(monkeypatch: pytest.MonkeyPatch, root: Path, code: int, at: str) -> None:
    """``Path.open(..., "wb")`` under ``root`` hits ``code``: on open, or on the first write
    (the file then exists, half-written, as on a real full disk)."""
    real_open = Path.open

    def fake_open(self: Path, mode: str = "r", *a: Any, **k: Any) -> Any:
        if "w" in mode and Path(self).absolute().is_relative_to(root.absolute()):
            if at == "open":
                raise OSError(code, os.strerror(code), str(self))
            return _FullFile(real_open(self, mode, *a, **k), code)
        return real_open(self, mode, *a, **k)

    monkeypatch.setattr(Path, "open", fake_open)


CODES = [pytest.param(errno.ENOSPC, id="ENOSPC"), pytest.param(errno.EDQUOT, id="EDQUOT")]


@pytest.mark.parametrize("at", ["open", "write"])
@pytest.mark.parametrize("code", CODES)
def test_b5_real_https_fetcher_reraises_disk_full(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, code: int, at: str
) -> None:
    monkeypatch.setattr(HttpsFetcher, "fetch", REAL_FETCH)
    for var in (
        "http_proxy",
        "https_proxy",
        "HTTP_PROXY",
        "HTTPS_PROXY",
        "all_proxy",
        "ALL_PROXY",
        "no_proxy",
        "NO_PROXY",
    ):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setattr(http.client, "HTTPSConnection", _Conn)
    out = tmp_path / "out"
    out.mkdir()
    _full_disk_under(monkeypatch, out, code, at)
    dest = out / "genomad_db.tar.gz"
    with pytest.raises(OSError) as info:
        HttpsFetcher().fetch(SOURCES["genomad"].url, dest)
    assert not isinstance(info.value, SourceError), "our disk, not the source"
    assert info.value.errno == code
    assert not dest.exists(), "no partial file is left"


@pytest.mark.parametrize("code", CODES)
def test_b5_real_fetcher_disk_full_makes_the_file_item_failed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, code: int
) -> None:
    monkeypatch.setattr(HttpsFetcher, "fetch", REAL_FETCH)
    monkeypatch.setattr(http.client, "HTTPSConnection", _Conn)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    db = tmp_path / "db"
    _full_disk_under(monkeypatch, db / "go", code, "write")
    opts = SetupOptions(db, bakta="skip", genomad=False, promotech_models=False, go=True)
    run_setup(opts, ToolRunner(db_dir=db), HttpsFetcher())
    item = _items(db)["go_basic"]
    assert item["status"] == "failed", item
    assert ADVICE in item["message"], item


# ---------------------------------------------------------------- quota device filtering

PROJ_DEV = "nfs.example.org:/export/proj"
HOME_DEV = "nfs.example.org:/export/home"


def test_b5_quota_remaining_filters_by_device() -> None:
    text = _quota_rows((PROJ_DEV, 0, 1000 * GIB // KIB), (HOME_DEV, 0, 1024))
    assert space.quota_remaining(text, [PROJ_DEV]) == 1000 * GIB
    assert space.quota_remaining(text, [HOME_DEV]) == 1024 * KIB
    assert space.quota_remaining(text, ["other:/x"]) is None


def test_b5_another_filesystems_quota_does_not_refuse_setup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    proj = _dir(tmp_path / "proj")
    _set_mounts(tmp_path, monkeypatch, [(str(proj), "nfs4", PROJ_DEV)])
    env["calls"].quota_out = _quota_rows((PROJ_DEV, 0, 1000 * GIB // KIB), (HOME_DEV, 0, 1))
    r = _setup("--db-dir", str(proj / "db"), *HTTPS_ITEMS)
    assert env["calls"].named("quota"), "setup did not ask quota"
    assert "on a larger disk" not in r.output, r.output
    assert env["fetched"]


# ---------------------------------------------------------------- mount_of


def test_b5_mount_of_overmount_takes_the_later_mount(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = _dir(tmp_path / "data")
    _dir(data / "x")
    _set_mounts(
        tmp_path, monkeypatch, [(str(data), "nfs4", "srv:/a"), (str(data), "ext4", "/dev/sdb1")]
    )
    assert space.mount_of(data / "x") == (str(data), "ext4", "/dev/sdb1")


def test_b5_mount_of_decodes_octal_escapes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    disk = _dir(tmp_path / "my disk")
    _dir(disk / "db")
    written = str(disk).replace(" ", "\\040")
    _set_mounts(tmp_path, monkeypatch, [(written, "nfs4", "srv:/my")])
    assert space.mount_of(disk / "db") == (str(disk), "nfs4", "srv:/my")


def test_b5_mount_of_bind_mount_uses_the_mount_point_not_the_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bind = _dir(tmp_path / "bind")
    f = tmp_path / "mountinfo"
    f.write_text(
        _mountinfo([("/", "ext4", "/dev/sda1")])
        + f"40 20 8:17 /export/sub {bind} rw,relatime shared:40 - xfs /dev/sdb1 rw\n"
    )
    monkeypatch.setattr(space, "MOUNTINFO", f)
    assert space.mount_of(bind) == (str(bind), "xfs", "/dev/sdb1")


# ---------------------------------------------------------------- per-filesystem grouping


def _review(lines: Sequence[Any], places: dict[str, Path], db: Path) -> int | None:
    return space.review_plan(
        lines, places, db, home_warning=False, yes=True, environ=dict(os.environ)
    )


@pytest.mark.parametrize(("free_gib", "refused"), [(25, True), (34, False)])
def test_b5_needs_on_one_filesystem_add_up(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, free_gib: int, refused: bool
) -> None:
    """DB dir, image cache and temp dir on one device: 10 + 10 + 10 GiB against one figure."""
    db, cache, tmp = _dir(tmp_path / "db"), _dir(tmp_path / "cache"), _dir(tmp_path / "tmp")
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(100 * GIB, 0, free_gib * GIB))
    lines = [
        space.PlanLine("genomad", 10.0, db),
        space.PlanLine("image x", 10.0, cache),
        space.PlanLine("temporary space", 10.0, tmp, temporary=True),
    ]
    places = {"DB directory": db, "image cache": cache, "temporary space": tmp}
    got = _review(lines, places, db)
    assert got == (int(ExitCode.ENVIRONMENT) if refused else None)


def test_b5_lfs_wrapped_name_under_a_device_filter() -> None:
    name = "/very/long/lustre/mount/point/for/alice"
    text = (
        "Disk quotas for usr alice (uid 1001):\n"
        "     Filesystem  kbytes   quota   limit   grace   files   quota   limit   grace\n"
        f"{name}\n"
        "                1048576       0 104857600       -     120       0       0       -\n"
    )
    assert space.quota_remaining(text, [name]) == (104857600 - 1048576) * KIB
    assert space.quota_remaining(text, ["/scratch"]) is None


def test_b5_in_shared_home_needs_a_path_boundary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = _dir(tmp_path / "homes")
    al, alice = _dir(root / "al"), _dir(root / "alice")
    _dir(alice / "db")
    _dir(al / "db")
    _set_mounts(tmp_path, monkeypatch, [(str(root), "nfs4", "srv:/export/home")])
    environ = {"HOME": str(al)}
    assert space.in_shared_home(alice / "db", environ) is False, "/home/alice is not in /home/al"
    assert space.in_shared_home(al / "db", environ) is True


# ---------------------------------------------------------------- plan lines: images

R1 = "quay.io/biocontainers/one@sha256:" + "a" * 64
R2 = "quay.io/biocontainers/two@sha256:" + "b" * 64


@pytest.fixture
def three_tools(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """``a`` and ``c`` share image R1 (3 GiB), ``b`` is R2 (5 GiB)."""
    from gippa.setupdb import images

    monkeypatch.setattr(images, "pinned_images", lambda *_a, **_k: {"a": R1, "b": R2, "c": R1})
    monkeypatch.setattr(images, "image_sizes", lambda *_a, **_k: {"a": 3.0, "b": 5.0, "c": 3.0})
    return ["a", "b", "c"]


def _gib(line: str) -> float:
    m = re.search(r"(\d+(?:\.\d+)?) GiB", line)
    assert m is not None, line
    return float(m.group(1))


def test_b5_temp_space_is_the_largest_image_and_shared_images_count_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    three_tools: list[str],
) -> None:
    db, cache, tmp = _dir(tmp_path / "db"), tmp_path / "cache", _dir(tmp_path / "tmp")
    monkeypatch.setenv("APPTAINER_TMPDIR", str(tmp))
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(100 * GIB, 0, 90 * GIB))
    lines, places = space.setup_plan(db, [], three_tools, "apptainer", cache, dict(os.environ))
    _review(lines, places, db)
    out = capsys.readouterr().err
    temp = [ln for ln in out.splitlines() if "plan:" in ln and str(tmp) in ln]
    assert temp, out
    assert _gib(temp[0]) == pytest.approx(5.0), "temp space = the largest image, not the sum"
    total = [ln for ln in out.splitlines() if re.search(r"plan: total\b", ln)]
    assert total, out
    assert _gib(total[0]) == pytest.approx(8.0), "a shared image is counted once: 3 + 5"


# ================================================================ B6: quota per place


def test_b6_image_cache_on_another_quotad_filesystem_is_checked(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    proj, home = _dir(tmp_path / "proj"), env["home"].resolve()
    _set_mounts(
        tmp_path, monkeypatch, [(str(proj), "nfs4", PROJ_DEV), (str(home), "nfs4", HOME_DEV)]
    )
    monkeypatch.setenv("APPTAINER_TMPDIR", str(_dir(proj / "tmp")))
    left_kib = int(0.3 * GIB) // KIB  # the images need more than 0.3 GiB
    limit_kib = 20 * GIB // KIB
    env["calls"].quota_out = _quota_rows(
        (PROJ_DEV, 0, 1000 * GIB // KIB), (HOME_DEV, limit_kib - left_kib, limit_kib)
    )
    cache = home / "cache"
    r = _setup(
        "--db-dir",
        str(proj / "db"),
        "--cache-dir",
        str(cache),
        *NO_DBS,
        "--images",
        "--profile",
        "apptainer",
        "--no-uniop",
    )
    assert r.exit_code == 2, r.output
    assert env["calls"].pulls == [], "nothing may be pulled when a filesystem is too small"
    assert "on a larger disk" in r.output
    cache_lines = [ln for ln in r.output.splitlines() if "free" in ln and str(cache) in ln]
    assert any("0.3 GiB" in ln for ln in cache_lines), r.output
    assert len(env["calls"].named("quota")) == 1, "one quota -w call, filtered per place"


def test_b6_one_quota_call_for_three_places(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    proj, home = _dir(tmp_path / "proj"), env["home"].resolve()
    _set_mounts(
        tmp_path, monkeypatch, [(str(proj), "nfs4", PROJ_DEV), (str(home), "nfs4", HOME_DEV)]
    )
    monkeypatch.setenv("APPTAINER_TMPDIR", str(_dir(home / "tmp")))
    env["calls"].quota_out = _quota_rows(
        (PROJ_DEV, 0, 1000 * GIB // KIB), (HOME_DEV, 0, 1000 * GIB // KIB)
    )
    r = _setup(
        "--db-dir",
        str(proj / "db"),
        "--cache-dir",
        str(home / "cache"),
        *HTTPS_DBS,
        "--images",
        "--profile",
        "apptainer",
        "--no-uniop",
    )
    assert env["calls"].pulls, r.output
    assert len(env["calls"].named("quota")) == 1, env["calls"].named("quota")


# ================================================================ B7: size not recorded


def test_b7_unrecorded_size_is_said_not_printed_as_zero(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    db = _dir(tmp_path / "db")
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(100 * GIB, 0, 90 * GIB))
    items = [("go_basic", 0.0, "go/go-basic.obo"), ("genomad", 4.0, "genomad_db")]
    lines, places = space.setup_plan(db, items, [], None, tmp_path / "c", dict(os.environ))
    _review(lines, places, db)
    out = capsys.readouterr().err
    (go,) = [ln for ln in out.splitlines() if "plan:" in ln and "go_basic" in ln]
    assert "size not recorded" in go, go
    assert "0.0 GiB" not in go, go
    (gm,) = [ln for ln in out.splitlines() if "plan:" in ln and "genomad" in ln]
    assert "4.0 GiB" in gm and "not recorded" not in gm


def test_b7_cli_plan_for_items_without_a_size(tmp_path: Path, env: dict[str, Any]) -> None:
    r = _setup("--db-dir", str(tmp_path / "db"), *NO_DBS, "--go", "--cog")
    plan = _plan(r.output)
    zero = [
        n
        for n in ("go_basic", "goslim_generic", "cog2020_def", "cog2020_fun")
        if not SOURCES[n].size_gib
    ]
    for name in zero:
        rows = [ln for ln in plan.splitlines() if re.search(rf"\b{name}\b", ln)]
        assert rows, f"no plan line for {name}:\n{plan}"
        assert all("size not recorded" in ln and "0.0 GiB" not in ln for ln in rows), rows


# ================================================================ prompt order (lead decision)


@pytest.fixture
def uniop_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Any]:
    """The 14u2 world: fake docker/apptainer/micromamba on PATH, a fake GitHub, the real
    ``subprocess.run`` and ``HttpsFetcher.fetch``."""
    from helpers.uniop_14u1 import (
        fixture_archive,
        install_fake_net,
        isolate,
        pin_uniop_source,
        sha256_hex,
    )
    from helpers.uniop_14u2 import install_fake_runtimes

    isolate(tmp_path, monkeypatch)
    install_fake_runtimes(tmp_path, monkeypatch)
    monkeypatch.setattr(subprocess, "run", REAL_RUN)
    monkeypatch.setattr(HttpsFetcher, "fetch", REAL_FETCH)
    net = install_fake_net(monkeypatch)
    data = fixture_archive()
    pin_uniop_source(monkeypatch, sha256_hex(data))
    net.serve_github(data)
    old = os.umask(0o022)
    yield net
    os.umask(old)


def _invoke(args: list[str], answer: str) -> tuple[int, str]:
    res = CliRunner().invoke(cli, ["setup", *args], input=answer)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    text = res.output
    with contextlib.suppress(ValueError):  # stderr not captured separately
        text += res.stderr
    return res.exit_code, text


def _uniop_case(db: Path, case: str) -> list[str]:
    # M11-15c2 turns antiSMASH and the resource files on by default; these cases are about the
    # images and UniOP, so the rest is opted out (opt_outs: only flags setup accepts).
    skip = ["--bakta", "skip", "--no-genomad", "--no-promotech-models",
            *opt_outs(*NO_RESOURCE_FLAGS)]  # fmt: skip
    if case == "images":
        return ["--db-dir", str(db), *skip, "--images", "--profile", "docker"]
    if case == "uniop+images":
        return ["--db-dir", str(db), *skip, "--uniop", "--images", "--profile", "docker"]
    return [
        "--db-dir",
        str(db),
        "--bakta",
        "skip",
        "--genomad",
        "--no-promotech-models",
        *opt_outs(*NO_RESOURCE_FLAGS),
        "--uniop",
        "--profile",
        "docker",
    ]


def test_prompt_order_proceed_then_terms(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, uniop_env: Any
) -> None:
    from helpers.uniop_14u2 import PROMPT

    _tty(monkeypatch)
    code, out = _invoke(_uniop_case(tmp_path / "db", "images"), "y\nn\n")
    assert "Proceed?" in out, out
    assert PROMPT in out, out
    assert out.index("Proceed?") < out.index(PROMPT), "Proceed? first, then the UniOP terms"
    assert code == 0, out


# `--uniop` installs UniOP only (Q-K27): with a database flag too, nothing else is planned, so
# Proceed? is never asked there; test_example_18c_tty_refusal_under_uniop_exits_2 covers that path.
@pytest.mark.parametrize("case", ["images", "uniop+images"])
def test_no_to_proceed_stops_before_any_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, uniop_env: Any, case: str
) -> None:
    from helpers.uniop_14u1 import docker_calls
    from helpers.uniop_14u2 import PROMPT

    _tty(monkeypatch)
    db = tmp_path / "db"
    code, out = _invoke(_uniop_case(db, case), "n\n")
    assert "Proceed?" in out, out
    assert code == 0, out
    assert PROMPT not in out, "the terms prompt never comes after a no"
    assert uniop_env.contacted == [], "nothing fetched (UniOP archive included)"
    assert docker_calls(tmp_path, "pull") == [] and docker_calls(tmp_path, "build") == []
    assert not (db / "db_manifest.json").exists() or all(
        i["status"] != "ok" for i in _items(db).values()
    )


# ================================================================ nits


def test_nit_unreadable_mountinfo_applies_no_quota(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    """Mount unknown: an unrelated filesystem's quota row must not refuse setup."""
    monkeypatch.setattr(space, "MOUNTINFO", tmp_path / "no-such-mountinfo")
    env["calls"].quota_out = _quota_rows((HOME_DEV, 0, 1))
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    assert "on a larger disk" not in r.output, r.output
    assert env["fetched"]


@pytest.mark.parametrize("shared", ["tmp=cache", "cache=db"])
def test_nit_two_roles_on_one_path_count_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, shared: str
) -> None:
    db, cache = _dir(tmp_path / "db"), _dir(tmp_path / "cache")
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(100 * GIB, 0, 25 * GIB))
    if shared == "tmp=cache":  # 10 + 10 = 20 GiB; 22 < 25
        lines = [
            space.PlanLine("image x", 10.0, cache),
            space.PlanLine("temporary space", 10.0, cache, temporary=True),
        ]
        places = {"DB directory": db, "image cache": cache, "temporary space": cache}
    else:  # 10 + 10 = 20 GiB; 22 < 25
        lines = [space.PlanLine("genomad", 10.0, db), space.PlanLine("image x", 10.0, db)]
        places = {"DB directory": db, "image cache": db}
    assert _review(lines, places, db) is None, "a path shared by two roles is counted once"


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores directory permissions")
def test_nit_read_only_config_dir_warns_and_does_not_crash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    cfg = tmp_path / "ro_config"
    cfg.mkdir()
    cfg.chmod(0o500)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(cfg))
    try:
        res = CliRunner().invoke(cli, ["setup", "--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS])
    finally:
        cfg.chmod(0o700)
    assert res.exception is None or isinstance(res.exception, SystemExit), repr(res.exception)
    assert "warning" in res.output.lower() and str(cfg) in res.output, res.output
    assert env["fetched"], "setup still downloads"
