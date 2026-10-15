"""``giap setup``: adopt or download each planned item, then record it (GIAP_SPEC.md §7.9.1, R6).

Downloads come only from the official source of each item (``db_sources.yaml``), into
``<path>.partial``; the copy is verified before it is renamed into place, and a failure leaves
the item absent with the host named. ``--check`` only verifies. One setup per ``db_dir`` at a
time (``.giap_setup.lock``).
"""

from __future__ import annotations

import dataclasses
import fcntl
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
import zipfile
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Protocol
from urllib.parse import urlparse

from giap.exitcodes import ExitCode
from giap.setupdb.items import (
    PROMOTECH_ARCHIVE,
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


class HttpsFetcher:
    """Plain HTTPS download of one official URL; no redirects to other schemes, no mirrors."""

    def fetch(self, url: str, dest: Path) -> None:
        host = urlparse(url).hostname or url
        if not url.startswith("https://"):
            raise SourceError(host, "not an https URL")
        try:
            with urllib.request.urlopen(url, timeout=60) as resp, dest.open("wb") as out:
                shutil.copyfileobj(resp, out, 8 << 20)
        except (urllib.error.URLError, OSError) as exc:
            raise SourceError(host, str(exc)) from exc


def _say(message: str) -> None:
    print(f"giap setup: {message}", file=sys.stderr, flush=True)


def _argv(item: DbItem, out: Path) -> list[str]:
    if item.check.startswith("bakta:"):
        return ["bakta_db", "download", "--output", str(out), "--type", item.check.split(":")[1]]
    if item.check == "genomad":
        return ["genomad", "download-database", str(out)]
    return ["download-antismash-databases", "--database-dir", str(out)]


def _absent(item: DbItem, message: str, status: str = "absent") -> DbRecord:
    rec = DbRecord(
        item.name, item.version, item.url, item.checksum, 0, "", item.path, "absent", message
    )
    return dataclasses.replace(rec, status="mismatch") if status == "mismatch" else rec


def _fetch_archive(item: DbItem, out: Path, fetcher: Fetcher) -> DbRecord | None:
    """HTTPS fetch, checksum, unzip; the archive is kept so ``verify`` can re-hash it."""
    archive = out / PROMOTECH_ARCHIVE
    fetcher.fetch(item.url, archive)
    got = file_digest(archive, item.checksum.split(":", 1)[0])
    if got != item.checksum:
        return _absent(item, f"checksum mismatch from {item.host}: {got}, pinned {item.checksum}",
                       "mismatch")  # fmt: skip
    with zipfile.ZipFile(archive) as z:
        for member in z.namelist():
            if Path(member).is_absolute() or ".." in Path(member).parts:
                return _absent(item, f"unsafe path {member!r} in archive", "mismatch")
        z.extractall(out)
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


def _provenance(rec: DbRecord, old: dict[str, DbRecord], prefetch: dict[str, tuple[str, str]],
                prefetch_key: str | None) -> DbRecord:  # fmt: skip
    """Keep the recorded provenance of an item that is already installed (idempotent reruns)."""
    prev = old.get(rec.name)
    if (
        prev is not None
        and prev.status == "ok"
        and (prev.version, prev.path) == (rec.version, rec.path)
    ):
        return prev
    if prefetch_key is not None and prefetch_key in prefetch:
        checksum, date = prefetch[prefetch_key]
        return dataclasses.replace(rec, checksum=checksum, retrieved=date)
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
            if item.name == "bakta_full":
                amr = amrfinder_record(item, db)
                recs.append(_provenance(amr, old, {}, None) if amr.status == "ok" else amr)
        write_manifest(db, recs)
    return ExitCode.SUCCESS if status_ok else ExitCode.ENVIRONMENT
