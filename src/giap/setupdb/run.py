"""``giap setup``: adopt or download each planned item, then record it (GIAP_SPEC.md §7.9.1, R6).

Downloads come only from the official source of each item (``db_sources.yaml``), into
``<path>.partial``; the copy is verified before it is renamed into place, and a failure leaves
the item absent with the host named. ``--check`` only verifies. One setup per ``db_dir`` at a
time (``.giap_setup.lock``).
"""

from __future__ import annotations

import dataclasses
import fcntl
import http.client
import os
import posixpath
import re
import shutil
import ssl
import subprocess
import sys
import tarfile
import zipfile
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
from typing import Protocol
from urllib.parse import urljoin, urlparse

from giap.exitcodes import ExitCode
from giap.setupdb.items import (
    DbItem,
    DbRecord,
    SetupOptions,
    check_dir,
    file_digest,
    plan,
    today,
    verify,
)
from giap.setupdb.manifest import read_prefetch, read_records, write_manifest

__all__ = [
    "LOCK",
    "CommandRunner",
    "Fetcher",
    "HttpsFetcher",
    "SourceError",
    "ToolRunner",
    "ToolUnavailable",
    "amrfinder_record",
    "install",
    "run_setup",
]  # fmt: skip

LOCK = ".giap_setup.lock"
GIB = 1024**3
AMR_URL = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/"


class SourceError(Exception):
    """An official source could not be reached; no other source is ever tried (R6)."""

    def __init__(self, host: str, reason: str) -> None:
        super().__init__(f"{host}: {reason}")
        self.host = host
        self.reason = reason


class CommandRunner(Protocol):
    def run(self, tool: str, argv: Sequence[str]) -> int: ...


class Fetcher(Protocol):
    def fetch(self, url: str, dest: Path) -> None: ...


class ToolUnavailable(Exception):
    """Neither a conda environment nor the pinned container image of a tool is installed."""


def _repo() -> Path:
    return Path(os.environ.get("GIAP_REPO") or Path(__file__).resolve().parents[3])


class ToolRunner:
    """Runs a downloader in ``giap-<tool>`` (micromamba) or else the image pinned in
    ``workflow/conf/containers.config`` via Apptainer (M2-08 lead decision 5)."""

    def command(self, tool: str, argv: Sequence[str]) -> list[str]:
        root = Path(os.environ.get("GIAP_ROOT") or _repo() / "local")
        mamba = Path(os.environ.get("MAMBA_ROOT_PREFIX") or root / "micromamba")
        if (mamba / "envs" / f"giap-{tool}").is_dir():
            return ["micromamba", "run", "-n", f"giap-{tool}", *argv]
        conf = (_repo() / "workflow/conf/containers.config").read_text(encoding="utf-8")
        m = re.search(rf"withName: '{tool.upper()}' \{{ container = 'docker://([^']+)'", conf)
        cache = Path(os.environ.get("NXF_APPTAINER_CACHEDIR") or root / "containers")
        img = cache / (re.sub(r"[/:]", "-", m.group(1)) + ".img") if m else None
        if img is None or not img.is_file():
            raise ToolUnavailable(
                f"no giap-{tool} environment and no pinned {tool} image ({img}); "
                "run sbatch scripts/slurm/build_envs.sbatch " + tool
            )
        binds = [b for a in argv if os.path.isabs(a) and Path(a).exists() for b in ("--bind", a)]
        return ["apptainer", "exec", *binds, str(img), *argv]

    def run(self, tool: str, argv: Sequence[str]) -> int:
        return subprocess.run(self.command(tool, argv), check=False).returncode


_REDIRECTS = frozenset({301, 302, 303, 307, 308})


class HttpsFetcher:
    """HTTPS download of one official URL. Redirects are followed only to https URLs on the same
    host, checked before the new location is contacted; no other scheme, host or mirror (R6)."""

    max_redirects = 5

    def fetch(self, url: str, dest: Path) -> None:
        host = urlparse(url).hostname or url
        current = url
        for _ in range(self.max_redirects + 1):
            parts = urlparse(current)
            if parts.scheme != "https" or parts.hostname != host:
                raise SourceError(host, f"refused {current}: only https on {host} is allowed")
            conn = http.client.HTTPSConnection(
                host, parts.port, timeout=60, context=ssl.create_default_context()
            )
            try:
                target = (parts.path or "/") + (f"?{parts.query}" if parts.query else "")
                conn.request("GET", target, headers={"User-Agent": "giap-setup"})
                resp = conn.getresponse()
                if resp.status in _REDIRECTS:
                    current = urljoin(current, resp.getheader("Location") or "")
                    continue
                if resp.status != 200:
                    raise SourceError(host, f"HTTP {resp.status} for {current}")
                with dest.open("wb") as out:
                    shutil.copyfileobj(resp, out, 8 << 20)
                return
            except (OSError, http.client.HTTPException) as exc:
                raise SourceError(host, str(exc)) from exc
            finally:
                conn.close()
        raise SourceError(host, f"more than {self.max_redirects} redirects")


def _say(message: str) -> None:
    print(f"giap setup: {message}", file=sys.stderr, flush=True)


def _argv(item: DbItem, out: Path) -> list[str]:
    if item.check.startswith("bakta:"):
        return ["bakta_db", "download", "--output", str(out), "--type", item.check.split(":")[1]]
    return ["download-antismash-databases", "--database-dir", str(out)]


def _absent(item: DbItem, message: str, status: str = "absent") -> DbRecord:
    rec = DbRecord(
        item.name, item.version, item.url, item.checksum, 0, "", item.path, "absent", message
    )
    return dataclasses.replace(rec, status="mismatch") if status == "mismatch" else rec


def _safe_name(name: str) -> bool:
    p = PurePosixPath(name)
    return bool(name) and not p.is_absolute() and ".." not in p.parts


def _safe_member(m: tarfile.TarInfo) -> bool:
    """Regular files and directories with relative names; links only within the archive."""
    if not _safe_name(m.name):
        return False
    if m.issym():
        dest = posixpath.normpath(posixpath.join(posixpath.dirname(m.name), m.linkname))
        return not posixpath.isabs(m.linkname) and dest != ".." and not dest.startswith("../")
    return m.isfile() or m.isdir() or (m.islnk() and _safe_name(m.linkname))


def _unpack(archive: Path, out: Path) -> str | None:
    """Extract ``archive`` into ``out``, or name the first unsafe member (nothing extracted)."""
    if archive.suffix == ".zip":
        with zipfile.ZipFile(archive) as z:
            bad = next((n for n in z.namelist() if not _safe_name(n)), None)
            if bad is None:
                z.extractall(out)
    else:
        with tarfile.open(archive) as tar:
            bad = next((m.name for m in tar.getmembers() if not _safe_member(m)), None)
            if bad is None:
                tar.extractall(out, filter="data")
    return None if bad is None else f"unsafe member {bad!r} in {archive.name}"


def _fetch_archive(item: DbItem, out: Path, fetcher: Fetcher) -> DbRecord | None:
    """HTTPS fetch, checksum, safe unpack. The Promotech zip is kept so ``verify`` can re-hash
    it (lead decision 4); other archives are deleted once unpacked (§3.9)."""
    archive = out / PurePosixPath(urlparse(item.url).path).name
    fetcher.fetch(item.url, archive)
    got = file_digest(archive, item.checksum.split(":", 1)[0])
    if got != item.checksum:
        return _absent(item, f"checksum mismatch from {item.host}: {got}, pinned {item.checksum}",
                       "mismatch")  # fmt: skip
    unsafe = _unpack(archive, out)
    if unsafe is not None:
        return _absent(item, f"{unsafe} from {item.host}", "mismatch")
    if item.check != "promotech":
        archive.unlink()
    return None


def install(item: DbItem, db_dir: Path, runner: CommandRunner, fetcher: Fetcher) -> DbRecord:
    """Adopt ``item`` if it verifies, else download it from its official source and verify it."""
    rec = verify(item, db_dir)
    if rec.status == "ok":
        return rec
    free = shutil.disk_usage(db_dir).free
    if free < item.size_gib * GIB:
        return _absent(
            item, f"only {free / GIB:.1f} GiB free in {db_dir}, need {item.size_gib} GiB"
        )
    target = db_dir / item.path
    partial = target.with_name(target.name + ".partial")
    shutil.rmtree(partial, ignore_errors=True)
    partial.mkdir(parents=True)
    try:
        if item.tool is None:
            failed = _fetch_archive(item, partial, fetcher)
            if failed is not None:
                return failed
        else:
            code = runner.run(item.tool, _argv(item, partial))
            if code != 0:
                return _absent(item, f"source {item.host} unreachable (downloader exit {code}; it "
                                     "also fails on a checksum mismatch)")  # fmt: skip
        built = partial / target.name if (partial / target.name).is_dir() else partial
        got = check_dir(item, built)
        if got.status != "ok":
            return _absent(item, f"download from {item.host} does not verify: {got.message}",
                           got.status)  # fmt: skip
        if target.exists():
            shutil.rmtree(target)
        built.rename(target)
    except ToolUnavailable as exc:
        return _absent(item, str(exc))
    except SourceError as exc:
        return _absent(item, f"source {exc.host} unreachable ({exc.reason})")
    finally:
        shutil.rmtree(partial, ignore_errors=True)
    return dataclasses.replace(verify(item, db_dir), message=f"downloaded from {item.host}")


def amrfinder_record(bakta: DbItem, db_dir: Path) -> DbRecord:
    """The AMRFinderPlus DB that ``bakta_db download`` pulls (lead decision Q2)."""
    path = f"{bakta.path}/amrfinderplus-db"
    try:
        version = (db_dir / path / "latest" / "version.txt").read_text(encoding="utf-8").strip()
    except OSError:
        return DbRecord("amrfinderplus", "", AMR_URL, "", 0, "", path, "absent",
                        f"no {path}/latest/version.txt")  # fmt: skip
    return DbRecord("amrfinderplus", version, AMR_URL, "", 0, today(), path, "ok")


def _provenance(
    rec: DbRecord,
    old: dict[str, DbRecord],
    prefetch: dict[str, tuple[str, str, str]],
    prefetch_key: str | None,
) -> DbRecord:
    """Provenance of an item that is already installed: a prefetch row is adopted only when its
    checksum and host equal the pin (else ``mismatch``); a previous ``ok`` record is kept
    (idempotent reruns) only when it carries the pinned URL and checksum (R6)."""
    row = prefetch.get(prefetch_key) if prefetch_key is not None else None
    host = urlparse(rec.url).hostname or ""
    if row is not None and (row[0], row[2]) != (rec.checksum, host):
        return dataclasses.replace(rec, size_bytes=0, retrieved="", status="mismatch", message=(
            f"prefetched copy is {row[0]} from {row[2] or 'an unknown host'}, "
            f"pinned {rec.checksum} from {host}"))  # fmt: skip
    prev = old.get(rec.name)
    pinned = (rec.version, rec.path, rec.url, rec.checksum)
    if prev is not None and prev.status == "ok" and (
        prev.version, prev.path, prev.url, prev.checksum) == pinned:  # fmt: skip
        return prev
    if row is not None:
        return dataclasses.replace(rec, retrieved=row[1])
    return rec


@contextmanager
def _locked(db_dir: Path) -> Iterator[bool]:
    with (db_dir / LOCK).open("a") as fh:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            yield False
            return
        try:
            yield True
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def _report(rec: DbRecord) -> None:
    detail = f" ({rec.message})" if rec.message else ""
    _say(f"{rec.name}: {rec.status} v{rec.version}{detail}" if rec.status == "ok" else
         f"{rec.name}: {rec.message or rec.status}")  # fmt: skip


def run_setup(opts: SetupOptions, runner: CommandRunner, fetcher: Fetcher) -> ExitCode:
    """Verify (``--check``) or adopt/download every planned item and update the manifest."""
    items = plan(opts)
    extra = [f for f in ("go", "kegg", "cog", "contaminants") if getattr(opts, f)]
    if extra:
        _say(f"--{', --'.join(extra)}: not handled by this version (M2-09)")
    if not items:
        return ExitCode.SUCCESS
    db = opts.db_dir
    if opts.check:
        checked = [verify(i, db) for i in items]
        for r in checked:
            _report(r)
        return ExitCode.SUCCESS if all(r.status == "ok" for r in checked) else ExitCode.ENVIRONMENT
    if any(i.name == "bakta_full" for i in items) and not (
        os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOBID")
    ):
        _say("warning: --bakta full downloads about 84 GB; run it as a Slurm batch job (R14)")
    db.mkdir(parents=True, exist_ok=True)
    with _locked(db) as acquired:
        if not acquired:
            _say(f"another giap setup is running on {db} ({LOCK}); try again later")
            return ExitCode.ENVIRONMENT
        old, prefetch = read_records(db), read_prefetch(db)
        recs: list[DbRecord] = []
        for item in items:
            rec = install(item, db, runner, fetcher)
            if rec.status == "ok" and not rec.message:  # adopted, not downloaded now
                rec = _provenance(rec, old, prefetch, item.prefetch)
            _report(rec)
            recs.append(rec)
        status_ok = all(r.status == "ok" for r in recs)
        for item in items:
            if item.check.startswith("bakta:"):
                amr = amrfinder_record(item, db)
                if amr.status == "ok":
                    recs.append(_provenance(amr, old, {}, None))
                elif item.name == "bakta_full":  # the light DB may ship without it
                    recs.append(amr)
        write_manifest(db, recs)
    return ExitCode.SUCCESS if status_ok else ExitCode.ENVIRONMENT
