"""M11-15c1: the ``gippa setup`` plan printed before any download, free space and quotas.

docs/specs/M11-15-portable-clusters.md item 6 ("Before any download (15c1)", "Quotas (15c1)"),
examples 6.4 and 6.10, the edge case "a DB directory on a filesystem whose free space cannot be
read (warning, proceed)", and the 15c1 acceptance criteria (plan/size/temp/total/free-space tests
with an injected ``disk_usage`` and quota runner, ``quota``/``lfs`` output fixtures, timeout).
The two licence lines and ``--minimal`` belong to 15c2 and are not tested here.

Interface pinned here (readings chosen where the spec leaves room):

- New module ``gippa.setupdb.space`` (not ``gippa.setupdb.plan``: that name is the re-exported
  ``plan()`` function of ``gippa.setupdb``). It reads the mount table from the module attribute
  ``MOUNTINFO`` (a ``Path``, ``/proc/self/mountinfo`` in production; looked up at call time, the
  tests point it at a fixture) to find the filesystem type and mount point of a path (longest
  mount-point prefix; item 9's ``fs_kind`` classes: ``nfs4`` is ``shared``).
- ``quota_remaining(text) -> int | None``: bytes left under the user's quota parsed from
  ``quota -w`` or ``lfs quota -u <user> <mount>`` output (KiB columns ``used soft hard``; a
  trailing ``*`` on ``used`` means over the soft limit; the limit is the soft one when nonzero,
  else the hard one; never negative; ``lfs`` may wrap a long filesystem name onto its own line).
  ``None`` when no row has a limit (``none``, both limits 0, unparseable).
- Quota commands run through ``subprocess.run(argv, ..., timeout=10)``: ``argv[0] == "quota"``
  with ``-w`` in ``argv``; on Lustre also ``["lfs", "quota", "-u", <user>, ..., <mount>]`` with
  the mount point last. ``<user>`` is ``$USER``/``$LOGNAME`` (both set here). A timeout, a missing
  command or a nonzero exit is ignored (best effort).
- Free space: ``shutil.disk_usage`` of the DB directory's, the image cache's and the temp
  directory's filesystems; each figure is printed as ``<x.y> GiB`` (one decimal, as the spec's
  ``12.0 GiB free``). A parsed quota lowers the DB directory's figure. Free < total x 1.1: exit 2
  before anything is fetched or pulled, with ``needs <t> GiB, <f> GiB free on <where>: pass
  --db-dir (or set GIPPA_DB_DIR) on a larger disk`` (the ``or use --minimal`` tail is 15c2's).
- The plan: one line per planned item with its name and its size in GiB (``size_gib`` of
  ``db_sources.yaml``; images: their sizes), a temporary-space line naming ``$APPTAINER_TMPDIR``
  (default ``/tmp``) when images are pulled with apptainer, a total line. All printed before the
  first fetch or pull.
- Prompt: ``Proceed? [Y/n]`` only on a terminal (detected through
  ``gippa.setupdb.uniop_build.stdin_is_tty``, or a ``stdin_is_tty`` in ``gippa.setupdb.space``;
  the tests patch both), default yes; ``--yes`` (new option) or no terminal: proceed, no prompt.
- Home warning: default DB directory (no ``--db-dir``, no ``GIPPA_DB_DIR``) under ``$HOME`` on a
  ``shared`` filesystem prints the spec's sentence; without a terminal setup proceeds.
- ``--check`` downloads nothing, so it is never refused for lack of space.

No test touches the network, a real tool, the real ``$HOME`` or the real ``$GIPPA_DB_DIR``:
``setup_guard`` blocks them; the fetcher is a fake that records the URL and raises
``SourceError`` (so "setup proceeded" is "the fetcher was called").
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, NamedTuple

import pytest
from click.testing import CliRunner
from helpers.prereqs import run_fake_prereq
from hypothesis import given, settings
from hypothesis import strategies as st
from setup_guard import setup_guard  # noqa: F401  (autouse: no network, no tools)

from gippa.cli import cli
from gippa.setupdb import SOURCES, HttpsFetcher, SourceError

GIB = 1024**3
KIB = 1024
HTTPS_ITEMS = ("--bakta", "skip", "--genomad", "--promotech-models", "--no-antismash")
NO_DBS = ("--bakta", "skip", "--no-genomad", "--no-promotech-models", "--no-antismash")
PLANNED = ("genomad", "promotech_models")
HOME_WARNING = (
    "the default DB directory is in your home directory, which usually has a small quota on "
    "clusters; pass --db-dir on a project or scratch filesystem"
)
MARK = "<<FIRST-FETCH-OR-PULL>>"


class _Usage(NamedTuple):
    total: int
    used: int
    free: int


def _total_gib() -> float:
    return float(sum(SOURCES[n].size_gib for n in PLANNED))


def _space() -> Any:
    try:
        from gippa.setupdb import space
    except ImportError as exc:
        pytest.fail(f"gippa.setupdb.space is not implemented ({exc})")
    return space


def _mountinfo(rows: Sequence[tuple[str, str, str]]) -> str:
    """``/proc/self/mountinfo`` lines for (mount point, fstype, source)."""
    out = []
    for i, (mnt, fstype, src) in enumerate(rows, start=20):
        out.append(f"{i} 1 0:{i} / {mnt} rw,relatime shared:{i} - {fstype} {src} rw\n")
    return "".join(out)


def _set_mounts(tmp: Path, mp: pytest.MonkeyPatch, rows: Sequence[tuple[str, str, str]]) -> None:
    f = tmp / "mountinfo"
    f.write_text(_mountinfo([("/", "ext4", "/dev/sda1"), *rows]))
    mp.setattr(_space(), "MOUNTINFO", f)


class Calls:
    """Fake ``subprocess.run`` for ``quota``/``lfs``/``apptainer``; everything else is refused."""

    def __init__(self) -> None:
        self.argv: list[list[str]] = []
        self.kwargs: list[dict[str, Any]] = []
        self.quota: Callable[[list[str]], Any] = lambda a: _done(a, 0, "")
        self.lfs: Callable[[list[str]], Any] = lambda a: _done(a, 0, "")

    def __call__(self, argv: Sequence[str], *a: Any, **k: Any) -> Any:
        argv = [str(x) for x in argv]
        if (done := run_fake_prereq(argv, *a, **k)) is not None:  # 15d: java -version
            return done
        self.argv.append(argv)
        self.kwargs.append(k)
        if argv[0] == "quota":
            return self.quota(argv)
        if argv[0] == "lfs":
            return self.lfs(argv)
        if argv[:2] == ["apptainer", "pull"]:
            print(MARK, file=sys.stderr, flush=True)
            return _done(argv, 1, "", "FATAL: fake pull refused")
        raise AssertionError(f"unexpected subprocess.run {argv}")

    def named(self, prog: str) -> list[tuple[list[str], dict[str, Any]]]:
        return [(a, k) for a, k in zip(self.argv, self.kwargs, strict=True) if a[0] == prog]


def _done(argv: Sequence[str], rc: int, out: str, err: str = "") -> Any:
    return subprocess.CompletedProcess(list(argv), rc, out, err)


@pytest.fixture(autouse=True)
def env(
    setup_guard: list[str],  # noqa: F811  (patched after the guard, so the fakes win)
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> dict[str, Any]:
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USER", "alice")
    monkeypatch.setenv("LOGNAME", "alice")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    for var in ("XDG_DATA_HOME", "GIPPA_DB_DIR", "GIPPA_ROOT", "NXF_APPTAINER_CACHEDIR",
                "APPTAINER_CACHEDIR", "APPTAINER_TMPDIR", "SINGULARITY_CACHEDIR"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    calls = Calls()
    monkeypatch.setattr(subprocess, "run", calls)
    fetched: list[str] = []

    def fetch(_self: Any, url: str, dest: Path) -> None:
        print(MARK, file=sys.stderr, flush=True)
        fetched.append(url)
        raise SourceError("fake.example", "refused by the test")

    monkeypatch.setattr(HttpsFetcher, "fetch", fetch)
    try:
        from gippa.setupdb import space
    except ImportError:
        pass
    else:  # deterministic mounts unless a test sets its own
        f = tmp_path / "mountinfo.default"
        f.write_text(_mountinfo([("/", "ext4", "/dev/sda1")]))
        monkeypatch.setattr(space, "MOUNTINFO", f, raising=False)
    return {"home": home, "calls": calls, "fetched": fetched}


def _tty(monkeypatch: pytest.MonkeyPatch) -> None:
    from gippa.setupdb import uniop_build

    monkeypatch.setattr(uniop_build, "stdin_is_tty", lambda: True)
    try:
        from gippa.setupdb import space
    except ImportError:
        return
    monkeypatch.setattr(space, "stdin_is_tty", lambda: True, raising=False)


def _setup(*args: str, stdin: str | None = None) -> Any:
    return CliRunner().invoke(cli, ["setup", *args], input=stdin, catch_exceptions=False)


def _before_first_fetch(output: str) -> str:
    return output.split(MARK, 1)[0]


def _gib_line(text: str, name: str) -> bool:
    return re.search(rf"(?m)^.*\b{re.escape(name)}\b.*\d+(\.\d+)? ?GiB", text) is not None


# ---------------------------------------------------------------- quota output parsing (pure)

QUOTA_W = """Disk quotas for user alice (uid 1001):
     Filesystem  blocks   quota   limit   grace   files   quota   limit   grace
nfs.example.org:/export/home  15728640  20971520 20971520            1234       0       0
"""

QUOTA_W_OVER = """Disk quotas for user alice (uid 1001):
     Filesystem  blocks   quota   limit   grace   files   quota   limit   grace
nfs.example.org:/export/home  20971600*  20971520 20971520   6days    1234       0       0
"""

QUOTA_NONE = "Disk quotas for user alice (uid 1001): none\n"

LFS = """Disk quotas for usr alice (uid 1001):
     Filesystem  kbytes   quota   limit   grace   files   quota   limit   grace
       /scratch 1048576       0 104857600       -     120       0       0       -
"""

LFS_WRAPPED = """Disk quotas for usr alice (uid 1001):
     Filesystem  kbytes   quota   limit   grace   files   quota   limit   grace
/very/long/lustre/mount/point/for/alice
                1048576       0 104857600       -     120       0       0       -
"""

LFS_NO_LIMIT = """Disk quotas for usr alice (uid 1001):
     Filesystem  kbytes   quota   limit   grace   files   quota   limit   grace
       /scratch 1048576       0       0       -     120       0       0       -
"""


@pytest.mark.parametrize(
    ("text", "want"),
    [
        (QUOTA_W, (20971520 - 15728640) * KIB),
        (QUOTA_W_OVER, 0),
        (QUOTA_NONE, None),
        (LFS, (104857600 - 1048576) * KIB),
        (LFS_WRAPPED, (104857600 - 1048576) * KIB),
        (LFS_NO_LIMIT, None),
        ("", None),
        ("quota: Cannot resolve mountpoint path /x\n", None),
    ],
    ids=[
        "quota-w",
        "quota-w-over",
        "quota-none",
        "lfs",
        "lfs-wrapped",
        "lfs-no-limit",
        "empty",
        "garbage",
    ],
)
def test_quota_remaining_fixtures(text: str, want: int | None) -> None:
    assert _space().quota_remaining(text) == want


@settings(max_examples=200, deadline=None)
@given(
    used=st.integers(0, 1 << 40),
    soft=st.integers(0, 1 << 40),
    hard=st.integers(0, 1 << 40),
    lfs=st.booleans(),
)
def test_quota_remaining_property(used: int, soft: int, hard: int, lfs: bool) -> None:
    """max(0, limit - used) KiB, limit = soft if soft else hard; None without a limit."""
    limit = soft or hard
    star = "*" if soft and used > soft else ""
    grace = "6days" if star else ("-" if lfs else "")
    head = "Disk quotas for usr alice (uid 1001):" if lfs else "Disk quotas for user alice (uid 1):"
    fs = "/scratch" if lfs else "server:/export/home"
    cols = "     Filesystem  blocks   quota   limit   grace   files   quota   limit   grace"
    tail = "-" if lfs else ""
    text = f"{head}\n{cols}\n{fs} {used}{star} {soft} {hard} {grace} 7 0 0 {tail}\n"
    got = _space().quota_remaining(text)
    if limit == 0:
        assert got is None
    else:
        assert got == max(0, limit - used) * KIB
        assert got >= 0


# ---------------------------------------------------------------- the plan


def test_plan_lists_items_sizes_and_total_before_fetching(
    tmp_path: Path, env: dict[str, Any]
) -> None:
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    assert env["fetched"], f"setup did not proceed to the downloads:\n{r.output}"
    plan = _before_first_fetch(r.output)
    for name in PLANNED:
        assert _gib_line(plan, name), f"no plan line with {name} and its size:\n{plan}"
    total = _total_gib()
    m = re.search(r"(?mi)^.*\btotal\b.*?(\d+(?:\.\d+)?) ?GiB", plan)
    assert m is not None, f"no total line before the first download:\n{plan}"
    assert float(m.group(1)) == pytest.approx(total, abs=0.05)


def test_plan_states_free_space_of_db_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                          env: dict[str, Any]) -> None:  # fmt: skip
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(900 * GIB, 0, 777 * GIB))
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    assert "777.0 GiB" in _before_first_fetch(r.output), r.output
    assert env["fetched"]


@pytest.mark.parametrize(("factor", "refused"), [(1.05, True), (0.5, True), (1.15, False)])
def test_free_space_below_total_times_1_1_exits_2_before_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any], factor: float,
    refused: bool,
) -> None:  # fmt: skip
    """Example 6.4: not enough space: exit 2 before any download, with the advice."""
    free = int(_total_gib() * factor * GIB)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(free * 2, free, free))
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    if refused:
        assert r.exit_code == 2, r.output
        assert env["fetched"] == [], "nothing may be fetched when the space check fails"
        assert f"{free / GIB:.1f} GiB free on " in r.output
        assert "pass --db-dir (or set GIPPA_DB_DIR) on a larger disk" in r.output
        m = re.search(r"needs (\d+(?:\.\d+)?) GiB", r.output)
        assert m is not None, r.output
        t = _total_gib()
        assert float(m.group(1)) in {round(t, 1), round(t * 1.1, 1)}
    else:
        assert env["fetched"], r.output
        assert "on a larger disk" not in r.output


def test_check_is_never_refused_for_space(tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                          env: dict[str, Any]) -> None:  # fmt: skip
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(GIB, GIB, 0))
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS, "--check")
    assert "on a larger disk" not in r.output
    assert env["fetched"] == []


def test_unreadable_free_space_warns_and_proceeds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    def broken(_p: Any) -> Any:
        raise OSError(5, "Input/output error")

    monkeypatch.setattr(shutil, "disk_usage", broken)
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    assert "warning" in _before_first_fetch(r.output).lower(), r.output
    assert env["fetched"] == [SOURCES[n].url for n in PLANNED]
    assert r.exit_code == 2  # the fake sources refuse; no crash


# ---------------------------------------------------------------- quotas


def _proj_mount(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fstype: str, src: str) -> Path:
    proj = (tmp_path / "proj").resolve()
    proj.mkdir()
    _set_mounts(tmp_path, monkeypatch, [(str(proj), fstype, src)])
    return proj


def _quota_text(device: str, used_kib: int, limit_kib: int) -> str:
    return (
        "Disk quotas for user alice (uid 1001): \n"
        "     Filesystem  blocks   quota   limit   grace   files   quota   limit   grace\n"
        f"{device} {used_kib} {limit_kib} {limit_kib}            10       0       0\n"
    )


def test_quota_lowers_free_figure_and_refuses(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    dev = "nfs.example.org:/export/proj"
    proj = _proj_mount(tmp_path, monkeypatch, "nfs4", dev)
    left = int(_total_gib() * 0.5 * GIB)  # filesystem: 10 TB free; quota: half the total
    limit_kib = 100 * GIB // KIB
    used_kib = limit_kib - left // KIB
    env["calls"].quota = lambda a: _done(a, 0, _quota_text(dev, used_kib, limit_kib))
    r = _setup("--db-dir", str(proj / "db"), *HTTPS_ITEMS)
    quota = env["calls"].named("quota")
    assert quota, "setup did not ask quota"
    assert "-w" in quota[0][0]
    assert quota[0][1].get("timeout") == 10
    assert r.exit_code == 2, r.output
    assert env["fetched"] == []
    assert f"{(limit_kib - used_kib) * KIB / GIB:.1f} GiB free on " in r.output


def test_quota_with_room_proceeds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    dev = "nfs.example.org:/export/proj"
    proj = _proj_mount(tmp_path, monkeypatch, "nfs4", dev)
    env["calls"].quota = lambda a: _done(a, 0, _quota_text(dev, 0, 1000 * GIB // KIB))
    _setup("--db-dir", str(proj / "db"), *HTTPS_ITEMS)
    assert env["calls"].named("quota")
    assert env["fetched"]


@pytest.mark.parametrize("failure", ["timeout", "missing", "exit1"])
def test_quota_failure_is_best_effort(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any], failure: str
) -> None:
    def quota(argv: list[str]) -> Any:
        if failure == "timeout":
            raise subprocess.TimeoutExpired(argv, 10)
        if failure == "missing":
            raise FileNotFoundError(2, "No such file or directory", "quota")
        return _done(argv, 1, "", "quota: error while getting quota from server")

    env["calls"].quota = quota
    _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    q = env["calls"].named("quota")
    assert q, "setup did not try quota"
    assert all(k.get("timeout") == 10 for _a, k in q)
    assert env["fetched"] == [SOURCES[n].url for n in PLANNED]


def test_lustre_asks_lfs_quota_for_the_mount(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    proj = _proj_mount(tmp_path, monkeypatch, "lustre", "10.0.0.1@o2ib:/scratch")
    left = int(_total_gib() * 0.5 * GIB)
    limit = 100 * GIB // KIB
    used = limit - left // KIB
    lfs = (
        "Disk quotas for usr alice (uid 1001):\n"
        "     Filesystem  kbytes   quota   limit   grace   files   quota   limit   grace\n"
        f"{proj} {used} 0 {limit} - 5 0 0 -\n"
    )
    env["calls"].quota = lambda a: _done(a, 0, QUOTA_NONE)
    env["calls"].lfs = lambda a: _done(a, 0, lfs)
    r = _setup("--db-dir", str(proj / "db"), *HTTPS_ITEMS)
    calls = env["calls"].named("lfs")
    assert calls, "setup did not ask lfs quota on Lustre"
    argv, kw = calls[0]
    assert argv[:4] == ["lfs", "quota", "-u", "alice"]
    assert argv[-1] == str(proj)
    assert kw.get("timeout") == 10
    assert r.exit_code == 2, r.output
    assert env["fetched"] == []


def test_no_lfs_off_lustre(tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                           env: dict[str, Any]) -> None:  # fmt: skip
    proj = _proj_mount(tmp_path, monkeypatch, "nfs4", "nfs.example.org:/export/proj")
    _setup("--db-dir", str(proj / "db"), *HTTPS_ITEMS)
    assert env["calls"].named("lfs") == []
    assert env["fetched"]


# ---------------------------------------------------------------- home warning (example 6.10)


def _home_on(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, home: Path, fstype: str) -> None:
    _set_mounts(tmp_path, monkeypatch, [(str(home.resolve()), fstype, "srv:/export/home")])


def test_default_db_dir_in_home_on_nfs_warns_and_proceeds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    """Example 6.10: default DB dir under ``$HOME`` on ``nfs4``, no TTY."""
    _home_on(tmp_path, monkeypatch, env["home"], "nfs4")
    r = _setup(*HTTPS_ITEMS)
    plan = _before_first_fetch(r.output)
    assert HOME_WARNING in " ".join(plan.split()), r.output
    assert env["fetched"] == [SOURCES[n].url for n in PLANNED], "without a TTY setup proceeds"
    assert "Proceed?" not in r.output


def test_no_home_warning_on_local_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                       env: dict[str, Any]) -> None:  # fmt: skip
    _home_on(tmp_path, monkeypatch, env["home"], "ext4")
    r = _setup(*HTTPS_ITEMS)
    assert "small quota" not in r.output
    assert env["fetched"]


@pytest.mark.parametrize("how", ["flag", "env"])
def test_no_home_warning_for_an_explicit_dir_in_home(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any], how: str
) -> None:
    _home_on(tmp_path, monkeypatch, env["home"], "nfs4")
    db = env["home"] / "mydb"
    if how == "env":
        monkeypatch.setenv("GIPPA_DB_DIR", str(db))
        r = _setup(*HTTPS_ITEMS)
    else:
        r = _setup("--db-dir", str(db), *HTTPS_ITEMS)
    assert "small quota" not in r.output
    assert env["fetched"]


# ---------------------------------------------------------------- the prompt


def test_terminal_asks_and_no_stops_before_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    _tty(monkeypatch)
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS, stdin="n\n")
    assert "Proceed? [Y/n]" in r.output
    assert env["fetched"] == []


def test_terminal_default_answer_is_yes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    _tty(monkeypatch)
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS, stdin="\n")
    assert "Proceed? [Y/n]" in r.output
    assert env["fetched"] == [SOURCES[n].url for n in PLANNED]


def test_yes_skips_the_prompt(tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                              env: dict[str, Any]) -> None:  # fmt: skip
    _tty(monkeypatch)
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS, "--yes")
    assert "no such option" not in r.output.lower(), r.output
    assert "Proceed?" not in r.output
    assert env["fetched"] == [SOURCES[n].url for n in PLANNED]


def test_no_terminal_no_prompt(tmp_path: Path, env: dict[str, Any]) -> None:
    r = _setup("--db-dir", str(tmp_path / "db"), *HTTPS_ITEMS)
    assert "Proceed?" not in r.output
    assert env["fetched"]
    assert _gib_line(_before_first_fetch(r.output), "genomad")


# ---------------------------------------------------------------- images: temp space, 3 filesystems


def _image_tools() -> list[str]:
    from gippa.config import resolve_config
    from gippa.setupdb.images import image_set

    return image_set(resolve_config(None, []), False)


def test_apptainer_plan_lists_images_temp_space_and_three_filesystems(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: dict[str, Any]
) -> None:
    db, cache, tmp = tmp_path / "A" / "db", tmp_path / "B" / "cache", tmp_path / "C" / "tmp"
    for d in (db, cache, tmp):
        d.mkdir(parents=True)
    monkeypatch.setenv("APPTAINER_TMPDIR", str(tmp))
    free = {str(db): 111, str(cache): 222, str(tmp): 333}

    def usage(p: Any) -> _Usage:
        s = str(Path(p))
        for root, gib in free.items():
            if s == root or s.startswith(root + "/"):
                return _Usage(1000 * GIB, 0, gib * GIB)
        return _Usage(10_000 * GIB, 0, 10_000 * GIB)

    monkeypatch.setattr(shutil, "disk_usage", usage)
    r = _setup("--db-dir", str(db), "--cache-dir", str(cache), *NO_DBS, "--images",
               "--profile", "apptainer", "--no-uniop")  # fmt: skip
    pulls = [a for a in env["calls"].argv if a[:2] == ["apptainer", "pull"]]
    assert pulls, f"setup did not proceed to the pulls:\n{r.output}"
    plan = _before_first_fetch(r.output)
    for tool in _image_tools():
        assert _gib_line(plan, tool), f"no plan line for image {tool}:\n{plan}"
    temp = [ln for ln in plan.splitlines() if str(tmp) in ln and "GiB" in ln]
    assert temp, f"no temporary-space line naming APPTAINER_TMPDIR:\n{plan}"
    for gib in (111, 222, 333):
        assert f"{gib}.0 GiB" in plan, f"free space {gib}.0 GiB not in the plan:\n{plan}"


def test_apptainer_temp_space_defaults_to_tmp(tmp_path: Path, env: dict[str, Any]) -> None:
    r = _setup("--db-dir", str(tmp_path / "db"), *NO_DBS, "--images", "--profile", "apptainer",
               "--no-uniop")  # fmt: skip
    plan = _before_first_fetch(r.output)
    lines = [ln for ln in plan.splitlines()
             if re.search(r"(?<![\w/.-])/tmp(?![\w/.-])", ln) and "GiB" in ln]  # fmt: skip
    assert lines, f"no temporary-space line for /tmp:\n{plan}"
