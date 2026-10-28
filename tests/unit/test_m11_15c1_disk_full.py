"""M11-15c1: ``EDQUOT``/``ENOSPC`` during a download or pull, no partial files; ``--disable-cache``.

docs/specs/M11-15-portable-clusters.md item 6 ("Quotas (15c1)" third bullet, "Apptainer's own
cache (15c1)"), examples 6.8 and 6.9, the ``.part`` rule of the 15c1 acceptance criteria.

Readings chosen where the spec leaves room:

- Example 6.8 says "fake fetcher raises ``OSError(EDQUOT)``" while downloading antiSMASH, but
  antiSMASH's databases come from its downloader through the command runner, not the fetcher.
  Both paths are pinned: the runner raising ``OSError(ENOSPC|EDQUOT)`` (antiSMASH), and the
  fetcher raising it (geNOMAD archive, GO file).
- The failed item is recorded in ``db_manifest.json`` with ``status: "failed"`` (a new
  ``Status`` value) and a message containing ``disk quota exceeded or disk full at <path>`` (a
  path inside the DB directory) and ``pass --db-dir``. The next items still run; exit 2.
- "The ``.part`` file": the code already downloads to ``<name>.partial`` and existing tests pin
  that name (tests/unit/test_m11_14f2_setup_images.py, test_setupdb.py); the spec's ``.part`` is
  read as that rule. Here: nothing whose name contains ``.part`` (so neither ``.part`` nor
  ``.partial``) is left, the final path is never written before completion.
- A pull that fails because the disk is full is detected from the runtime's message
  (``no space left on device`` / ``disk quota exceeded``, case-insensitive): item ``failed`` with
  the same advice and the partial ``.img`` removed.
- Example 6.9: ``apptainer pull --disable-cache <cache>/<name>.img docker://<ref>`` exactly when
  ``APPTAINER_CACHEDIR`` is unset or empty; with it set the argv has no ``--disable-cache``.
  ``singularity pull`` (``SINGULARITY_CACHEDIR``) needs the 15f runtime and is not tested here.
"""

from __future__ import annotations

import errno
import json
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Any, NamedTuple

import pytest
from setup_guard import setup_guard  # noqa: F401  (autouse: no network, no tools)

from gippa.exitcodes import ExitCode
from gippa.setupdb import SOURCES, SetupOptions, SourceError, run_setup

GIB = 1024**3
ADVICE = "disk quota exceeded or disk full at"
ERRNOS = [pytest.param(errno.EDQUOT, id="EDQUOT"), pytest.param(errno.ENOSPC, id="ENOSPC")]


class _Usage(NamedTuple):
    total: int
    used: int
    free: int


@pytest.fixture(autouse=True)
def _env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    for var in ("XDG_DATA_HOME", "GIPPA_DB_DIR", "GIPPA_ROOT", "NXF_APPTAINER_CACHEDIR",
                "APPTAINER_CACHEDIR", "SINGULARITY_CACHEDIR"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))


def _full(code: int, where: Path) -> OSError:
    return OSError(code, "Disk quota exceeded" if code == errno.EDQUOT else
                   "No space left on device", str(where))  # fmt: skip


class FullFetcher:
    """Writes half a file to ``dest`` for the ``full`` URLs, then raises the disk error; other
    URLs are refused by the source (recorded, so "the next item still ran" is visible)."""

    def __init__(self, code: int, full: Sequence[str]) -> None:
        self.code = code
        self.full = set(full)
        self.calls: list[tuple[str, Path]] = []

    def fetch(self, url: str, dest: Path) -> None:
        self.calls.append((url, dest))
        if url in self.full:
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(b"half a download")
            raise _full(self.code, dest)
        raise SourceError("fake.example", "refused by the test")


class FullRunner:
    """antiSMASH not in the image; its downloader writes a file and then hits the disk error."""

    def __init__(self, code: int) -> None:
        self.code = code
        self.calls: list[tuple[str, list[str]]] = []

    def run(self, tool: str, argv: Sequence[str]) -> int:
        self.calls.append((tool, list(argv)))
        if "--check-prereqs" in argv:
            return 1
        out = Path(argv[list(argv).index("--database-dir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / "pfam.hmm").write_bytes(b"half")
        raise _full(self.code, out / "pfam.hmm")


def _items(db: Path) -> dict[str, dict[str, Any]]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    return {i["name"]: i for i in data["items"]}


def _partials(db: Path) -> list[Path]:
    return [p for p in db.rglob("*") if ".part" in p.name]


def _assert_failed(item: dict[str, Any], db: Path) -> None:
    assert item["status"] == "failed", item
    assert ADVICE in item["message"], item
    assert "pass --db-dir" in item["message"], item
    assert str(db) in item["message"], item


@pytest.mark.parametrize("code", ERRNOS)
def test_example_6_8_quota_hit_in_antismash_and_genomad(tmp_path: Path, code: int) -> None:
    """Example 6.8: the item fails with the advice, no partial left, the next items still run."""
    db = tmp_path / "db"
    fetcher = FullFetcher(code, [SOURCES["genomad"].url])
    runner = FullRunner(code)
    opts = SetupOptions(db, bakta="skip", genomad=True, antismash=True, promotech_models=True)
    got = run_setup(opts, runner, fetcher)
    assert got == ExitCode.ENVIRONMENT
    items = _items(db)
    _assert_failed(items["genomad"], db)
    _assert_failed(items["antismash"], db)
    assert [u for u, _d in fetcher.calls] == [SOURCES["genomad"].url,
                                               SOURCES["promotech_models"].url]  # fmt: skip
    assert items["promotech_models"]["status"] == "absent"  # ran, refused by the fake source
    assert _partials(db) == []
    assert not (db / SOURCES["genomad"].path).exists()
    assert not (db / SOURCES["antismash"].path).exists()


@pytest.mark.parametrize("code", ERRNOS)
def test_file_item_disk_full_leaves_nothing_and_continues(tmp_path: Path, code: int) -> None:
    db = tmp_path / "db"
    fetcher = FullFetcher(code, [SOURCES["go_basic"].url])
    opts = SetupOptions(db, bakta="skip", genomad=False, promotech_models=False, go=True)
    got = run_setup(opts, FullRunner(code), fetcher)
    assert got == ExitCode.ENVIRONMENT
    items = _items(db)
    _assert_failed(items["go_basic"], db)
    assert [u for u, _d in fetcher.calls] == [SOURCES["go_basic"].url,
                                               SOURCES["goslim_generic"].url]  # fmt: skip
    target = db / SOURCES["go_basic"].path
    dest = fetcher.calls[0][1]
    assert dest != target, "downloads go to a partial name, renamed only when complete"
    assert ".part" in dest.name
    assert not target.exists()
    assert _partials(db) == []


def test_other_items_unaffected_by_a_disk_error_record(tmp_path: Path) -> None:
    """A disk error is not a source error: the message never blames the official host."""
    db = tmp_path / "db"
    fetcher = FullFetcher(errno.ENOSPC, [SOURCES["genomad"].url])
    run_setup(SetupOptions(db, bakta="skip", genomad=True, promotech_models=False),
              FullRunner(errno.ENOSPC), fetcher)  # fmt: skip
    msg = _items(db)["genomad"]["message"]
    assert "unreachable" not in msg
    assert SOURCES["genomad"].host not in msg


# ---------------------------------------------------------------- pulls (example 6.9)

REF = "quay.io/biocontainers/fake@sha256:" + "ab" * 32


class PullCalls:
    def __init__(self, rc: int = 0, err: str = "") -> None:
        self.argv: list[list[str]] = []
        self.rc = rc
        self.err = err

    def __call__(self, argv: Sequence[str], *a: Any, **k: Any) -> Any:
        argv = [str(x) for x in argv]
        self.argv.append(argv)
        if argv[:2] != ["apptainer", "pull"] and argv[:2] != ["docker", "pull"]:
            raise AssertionError(f"unexpected subprocess.run {argv}")
        if argv[0] == "apptainer":
            Path(argv[-2]).write_text("sif" if self.rc == 0 else "half a sif")
            return subprocess.CompletedProcess(argv, self.rc, "", self.err)
        return subprocess.CompletedProcess(argv, 0, f"Digest: {REF.split('@')[1]}\n", "")


@pytest.fixture
def pull(monkeypatch: pytest.MonkeyPatch) -> Any:
    from gippa.setupdb import images

    monkeypatch.setattr(images, "pinned_images", lambda *_a, **_k: {"faketool": REF})
    monkeypatch.setattr(images, "image_versions", lambda *_a, **_k: {"faketool": "1.0"})

    def go(tmp_path: Path, runtime: str, calls: PullCalls) -> tuple[int, Path, Path]:
        monkeypatch.setattr(subprocess, "run", calls)
        cache, db = tmp_path / "cache", tmp_path / "db"
        code = images.install_images(["faketool"], runtime, cache, db)
        return code, cache, db

    return go


def _simple(ref: str) -> str:
    return ref.replace("/", "-").replace(":", "-") + ".img"


@pytest.mark.parametrize("cachedir", [None, ""])
def test_example_6_9_apptainer_pull_disables_its_cache(
    tmp_path: Path, pull: Any, monkeypatch: pytest.MonkeyPatch, cachedir: str | None
) -> None:
    if cachedir is not None:
        monkeypatch.setenv("APPTAINER_CACHEDIR", cachedir)
    calls = PullCalls()
    code, cache, _db = pull(tmp_path, "apptainer", calls)
    assert code == ExitCode.SUCCESS
    assert calls.argv == [["apptainer", "pull", "--disable-cache", str(cache / _simple(REF)),
                           f"docker://{REF}"]]  # fmt: skip


def test_example_6_9_user_cachedir_left_alone(
    tmp_path: Path, pull: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("APPTAINER_CACHEDIR", str(tmp_path / "x"))
    calls = PullCalls()
    code, cache, _db = pull(tmp_path, "apptainer", calls)
    assert code == ExitCode.SUCCESS
    (argv,) = calls.argv
    assert "--disable-cache" not in argv
    assert argv[:2] == ["apptainer", "pull"]
    assert argv[-2:] == [str(cache / _simple(REF)), f"docker://{REF}"]


def test_docker_pull_unchanged(tmp_path: Path, pull: Any) -> None:
    calls = PullCalls()
    code, _cache, _db = pull(tmp_path, "docker", calls)
    assert code == ExitCode.SUCCESS
    assert calls.argv == [["docker", "pull", REF]]


@pytest.mark.parametrize(
    "err",
    [
        "FATAL: While making image from oci registry: write /x/blob: no space left on device",
        "FATAL: unable to write /x/sif: Disk quota exceeded",
    ],
    ids=["ENOSPC", "EDQUOT"],
)
def test_pull_on_full_disk_fails_with_advice_and_no_partial_image(
    tmp_path: Path, pull: Any, err: str
) -> None:
    calls = PullCalls(rc=255, err=err)
    code, cache, db = pull(tmp_path, "apptainer", calls)
    assert code == ExitCode.ENVIRONMENT
    item = _items(db)["image:faketool"]
    assert item["status"] == "failed", item
    assert ADVICE in item["message"] and "pass --db-dir" in item["message"], item
    assert not (cache / _simple(REF)).exists(), "the partial image must be removed"
