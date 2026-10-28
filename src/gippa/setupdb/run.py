"""``gippa setup``: adopt or download each planned item, then record it (GIPPA_SPEC.md §7.9.1, R6).

Downloads come only from the official source of each item (``db_sources.yaml``), into
``<path>.partial``; the copy is verified before it is renamed into place, and a failure leaves
the item absent with the host named. ``--check`` only verifies. One setup per ``db_dir`` at a
time (``.giap_setup.lock``).
"""

from __future__ import annotations

import codecs
import dataclasses
import fcntl
import http.client
import os
import posixpath
import shutil
import ssl
import subprocess
import sys
import tarfile
import threading
import zipfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Protocol
from urllib.parse import urljoin, urlparse

from gippa.exitcodes import ExitCode
from gippa.io.image_key import image_cache
from gippa.setupdb.images import pinned_images, simple_name
from gippa.setupdb.items import (
    FILE_CHECKS,
    SOURCES,
    DbItem,
    DbRecord,
    SetupOptions,
    Status,
    _safe_name,
    check_dir,
    check_file,
    file_digest,
    plan,
    today,
    verify,
)
from gippa.setupdb.manifest import read_prefetch, read_records, write_manifest
from gippa.setupdb.space import DISK_FULL, disk_full, is_disk_full

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
    "setup_nextflow",
]  # fmt: skip

LOCK = ".giap_setup.lock"
GIB = 1024**3
_TAIL = 65536  # bytes of a downloader's stderr kept to tell a full disk
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


def _checkout() -> Path | None:
    """``$GIPPA_REPO``, else this source checkout; ``None`` for an installed wheel."""
    if os.environ.get("GIPPA_REPO"):
        return Path(os.environ["GIPPA_REPO"])
    here = Path(__file__).resolve().parents[3]
    return here if (here / "workflow" / "main.nf").is_file() else None


def _containers_config() -> Path:
    """The pinned images: the checkout's, else the workflow bundled in the package (M11-01)."""
    repo = _checkout()
    if repo is not None:
        return repo / "workflow" / "conf" / "containers.config"
    from gippa.runner.nextflow import workflow_main

    return workflow_main().parent / "conf" / "containers.config"


# Images whose ENTRYPOINT is a wrapper: docker runs the downloader as the entrypoint instead.
_DOCKER_ENTRYPOINT = frozenset({"antismash"})


def _mamba_root() -> Path | None:
    """``$MAMBA_ROOT_PREFIX``, else ``$GIPPA_ROOT/micromamba``, else the checkout's
    ``local/micromamba``; ``None`` for an installed package without either variable."""
    if os.environ.get("MAMBA_ROOT_PREFIX"):
        return Path(os.environ["MAMBA_ROOT_PREFIX"])
    if os.environ.get("GIPPA_ROOT"):
        return Path(os.environ["GIPPA_ROOT"]) / "micromamba"
    repo = _checkout()
    return repo / "local" / "micromamba" if repo is not None else None


class ToolRunner:
    """Runs a downloader in ``gippa-<tool>`` (micromamba), else in the image pinned in the bundled
    ``containers.config``: the apptainer cache file (``image_cache``), else ``docker run`` (M2-08
    lead decision 5; M11-14 14f2). Neither ``$GIPPA_ROOT`` nor a checkout is required. With
    ``cache_dir`` the image is looked up there and nowhere else (``gippa setup`` passes the
    directory it pulls into)."""

    def __init__(self, *, db_dir: Path | None = None, cache_dir: Path | None = None) -> None:
        self.db_dir = db_dir
        self.cache_dir = cache_dir

    def _cache(self) -> Path | None:
        if self.cache_dir is not None:
            return self.cache_dir
        if self.db_dir is not None:
            return image_cache(os.environ, self.db_dir)
        if os.environ.get("NXF_APPTAINER_CACHEDIR"):
            return Path(os.environ["NXF_APPTAINER_CACHEDIR"])
        if os.environ.get("GIPPA_ROOT"):
            return Path(os.environ["GIPPA_ROOT"]) / "containers"
        return None

    def command(self, tool: str, argv: Sequence[str]) -> list[str]:
        mamba = _mamba_root()
        if mamba is not None and (mamba / "envs" / f"gippa-{tool}").is_dir():
            return ["micromamba", "run", "-n", f"gippa-{tool}", *argv]
        ref = pinned_images(_containers_config().parent).get(tool)
        cache = self._cache()
        img = cache / simple_name(ref) if ref is not None and cache is not None else None
        paths = list(dict.fromkeys(a for a in argv if os.path.isabs(a) and Path(a).exists()))
        if img is not None and img.is_file():
            binds = [b for a in paths for b in ("--bind", a)]
            return ["apptainer", "exec", *binds, str(img), *argv]
        if ref is not None and shutil.which("docker") is not None:
            vols = [v for a in paths for v in ("-v", f"{a}:{a}")]
            head = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", *vols]
            if tool in _DOCKER_ENTRYPOINT and argv:
                return [*head, "--entrypoint", argv[0], ref, *argv[1:]]
            return [*head, ref, *argv]
        where = f"no pinned {tool} image" if ref is None else f"no pinned {tool} image in {img}"
        raise ToolUnavailable(
            f"no gippa-{tool} environment, {where} and no docker on PATH; run gippa setup --images "
            f"--profile apptainer|docker (in a checkout: sbatch scripts/slurm/build_envs.sbatch "
            f"{tool})"
        )

    last_stderr: str = ""  # the tail of the last run's stderr (a full disk is told from it)

    def run(self, tool: str, argv: Sequence[str]) -> int:
        """Run ``argv``; its stderr still reaches the user live and its tail is kept."""
        env = tool_env(os.environ)
        read_fd, write_fd = os.pipe()
        tail = bytearray()
        sink = sys.stderr

        def tee() -> None:
            decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
            with os.fdopen(read_fd, "rb", buffering=0) as src:
                while chunk := src.read(65536):
                    tail.extend(chunk)
                    del tail[:-_TAIL]
                    sink.write(decoder.decode(chunk))
                    sink.flush()

        reader = threading.Thread(target=tee, daemon=True)
        reader.start()
        try:
            cmd = self.command(tool, argv)
            code = subprocess.run(cmd, check=False, env=env, stderr=write_fd).returncode
        finally:
            os.close(write_fd)
            reader.join(timeout=10)
        self.last_stderr = tail.decode("utf-8", errors="replace")
        return code


def tool_env(environ: Mapping[str, str]) -> dict[str, str]:
    """A copy of ``environ`` with ``PYTHONNOUSERSITE=1``: a tool's Python must not import
    packages from the user site (``~/.local``) over its pinned environment."""
    return {**environ, "PYTHONNOUSERSITE": "1"}


_REDIRECTS = frozenset({301, 302, 303, 307, 308})
# Official hosts of the same organisation a source may redirect to (M2-09 lead decision 1; R6).
REDIRECT_HOSTS = {
    # the UniOP archive: GitHub answers 302 to codeload, which streams it (M11-14 item 1b)
    "github.com": frozenset({"codeload.github.com"}),
    "purl.obolibrary.org": frozenset({"release.geneontology.org", "current.geneontology.org"}),
}


class HttpsFetcher:
    """HTTPS download of one official URL. Redirects are followed only to https URLs on the same
    host or on that source's allowlisted official hosts, checked before the new location is
    contacted; no other scheme, host or mirror (R6). A short body is an error; a body without
    Content-Length stops at the item's ``size_gib`` (1 GiB when unset or for a URL that is no
    item) and is an error too; nothing is left."""

    max_redirects = 5

    def fetch(self, url: str, dest: Path) -> None:
        try:
            self._fetch(url, dest)
        except BaseException:
            dest.unlink(missing_ok=True)
            raise

    def _fetch(self, url: str, dest: Path) -> None:
        host = urlparse(url).hostname or url
        allowed = {host} | _redirect_hosts(url, host)
        current = url
        for _ in range(self.max_redirects + 1):
            parts = urlparse(current)
            if parts.scheme != "https" or parts.hostname not in allowed:
                raise SourceError(host, f"refused {current}: only https on {sorted(allowed)}")
            conn = http.client.HTTPSConnection(
                parts.hostname, parts.port, timeout=60, context=ssl.create_default_context()
            )
            try:
                target = (parts.path or "/") + (f"?{parts.query}" if parts.query else "")
                conn.request("GET", target, headers={"User-Agent": "gippa-setup"})
                resp = conn.getresponse()
                if resp.status in _REDIRECTS:
                    current = urljoin(current, resp.getheader("Location") or "")
                    continue
                if resp.status != 200:
                    raise SourceError(host, f"HTTP {resp.status} for {current}")
                want = resp.getheader("Content-Length")
                cap = None if want is not None and want.isdigit() else _body_cap(url)
                with dest.open("wb") as out:
                    got = _copy(resp, out, cap)
                if cap is not None and got > cap:
                    raise SourceError(host, f"body without Content-Length exceeds {cap} bytes "
                                            f"from {current}")  # fmt: skip
                if want is not None and want.isdigit() and got != int(want):
                    raise SourceError(host, f"truncated body: {got} of {want} bytes from {current}")
                return
            except (OSError, http.client.HTTPException) as exc:
                if isinstance(exc, OSError) and exc.errno in DISK_FULL:
                    raise  # our disk, not the source: the item is failed (M11-15 15c1)
                raise SourceError(host, str(exc)) from exc
            finally:
                conn.close()
        raise SourceError(host, f"more than {self.max_redirects} redirects")


def _redirect_hosts(url: str, host: str) -> frozenset[str]:
    """The redirect allowlist of the item(s) at ``url``: their ``redirect_hosts`` when a row
    names some (M11-15 item 5), else the ``REDIRECT_HOSTS`` default of ``host``."""
    rows = [h for i in SOURCES.values() if i.url == url for h in i.redirect_hosts]
    return frozenset(rows) if rows else REDIRECT_HOSTS.get(host, frozenset())


def _body_cap(url: str) -> int:
    """Byte cap for a body of unknown length: ``size_gib`` of the item(s) at ``url``, else 1 GiB
    (M2-09 review decision 4)."""
    gib = max((i.size_gib for i in SOURCES.values() if i.url == url), default=0)
    return (gib or 1) * GIB


def _copy(src: BinaryIO, out: BinaryIO, cap: int | None) -> int:
    """Copy ``src`` to ``out`` in 8 MiB blocks; with a ``cap``, stop after ``cap + 1`` bytes."""
    got = 0
    while cap is None or got <= cap:
        n = 8 << 20 if cap is None else min(8 << 20, cap + 1 - got)
        block = src.read(n)
        if not block:
            break
        out.write(block)
        got += len(block)
    return got


def _say(message: str) -> None:
    print(f"gippa setup: {message}", file=sys.stderr, flush=True)


def _argv(item: DbItem, out: Path) -> list[str]:
    if item.check.startswith("bakta:"):
        return ["bakta_db", "download", "--output", str(out), "--type", item.check.split(":")[1]]
    return ["download-antismash-databases", "--database-dir", str(out)]


def _absent(item: DbItem, message: str, status: Status = "absent") -> DbRecord:
    return DbRecord(
        item.name, item.version, item.url, item.checksum, 0, "", item.path, status, message
    )


def _full(item: DbItem, db_dir: Path) -> DbRecord:
    """``EDQUOT``/``ENOSPC`` while writing ``item``: failed, with the advice (M11-15 15c1); the
    message names the item's final path, not its ``.partial``."""
    return _absent(item, disk_full(db_dir / item.path), "failed")


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


def _install_file(item: DbItem, db_dir: Path, fetcher: Fetcher) -> DbRecord:
    """Fetch one file to ``<path>.partial``, verify it, rename it into place (M2-09). Unpinned
    files (KEGG) are recorded with the sha256 of the bytes and the retrieval date as version."""
    target = db_dir / item.path
    partial = target.with_name(target.name + ".partial")
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        fetcher.fetch(item.url, partial)
        algo = item.checksum.partition(":")[0]
        if item.checksum and (got := file_digest(partial, algo)) != item.checksum:
            return _absent(item, f"checksum mismatch from {item.host}: pulled {got}, pinned "
                                 f"{item.checksum}", "mismatch")  # fmt: skip
        rec = check_file(item, partial, item.checksum, item.version or today())
        if rec.status != "ok":
            return _absent(item, f"download from {item.host} does not verify: {rec.message}",
                           rec.status)  # fmt: skip
        if item.executable:
            partial.chmod(0o755)
        partial.replace(target)
    except SourceError as exc:
        return _absent(item, f"source {exc.host} unreachable ({exc.reason})")
    except OSError as exc:
        if exc.errno not in DISK_FULL:
            raise
        return _full(item, db_dir)
    finally:
        partial.unlink(missing_ok=True)
    return dataclasses.replace(rec, message=f"downloaded from {item.host}")


def _antismash_image(item: DbItem, runner: CommandRunner) -> DbRecord | None:
    """``ok`` with url "image" when ``antismash --check-prereqs`` passes on the image's own
    databases (M2-09 lead decision 3), else None."""
    try:
        if runner.run("antismash", ["antismash", "--check-prereqs"]) != 0:
            return None
    except ToolUnavailable:
        return None
    return DbRecord(item.name, item.version, "image", "", 0, today(), item.path, "ok",
                    "databases in the antiSMASH image")  # fmt: skip


def install(item: DbItem, db_dir: Path, runner: CommandRunner, fetcher: Fetcher) -> DbRecord:
    """Adopt ``item`` if it verifies, else download it from its official source and verify it."""
    rec = verify(item, db_dir)
    if rec.status == "ok":
        return rec
    if item.check == "antismash" and (image := _antismash_image(item, runner)) is not None:
        return image
    try:
        free: int | None = shutil.disk_usage(db_dir).free
    except OSError:
        free = None  # unreadable: the plan has warned; the download decides
    if free is not None and free < item.size_gib * GIB:
        return _absent(
            item, f"only {free / GIB:.1f} GiB free in {db_dir}, need {item.size_gib} GiB"
        )
    if item.check.partition(":")[0] in FILE_CHECKS:
        return _install_file(item, db_dir, fetcher)
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
            if code != 0 and is_disk_full(getattr(runner, "last_stderr", "")):
                return _full(item, db_dir)
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
    except OSError as exc:
        if exc.errno not in DISK_FULL:
            raise
        return _full(item, db_dir)
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


def _pin_changed(item: DbItem, rec: DbRecord) -> bool:
    """An installed copy of another release than this one's pin (Bakta: another major; geNOMAD:
    another version), as read from its version file by ``verify``."""
    if rec.status == "ok" or not rec.version or rec.version == item.version:
        return False
    kind = item.check.partition(":")[0]
    if kind == "bakta":
        return rec.version.split(".")[0] != item.version.split(".")[0]
    return kind == "genomad"


def _check(item: DbItem, db_dir: Path, old: dict[str, DbRecord], runner: CommandRunner) -> DbRecord:
    """``--check`` of one item. antiSMASH recorded as installed from the image (url "image")
    passes when ``antismash --check-prereqs`` passes there (M2-09 review decision 1)."""
    rec = verify(item, db_dir)
    if _pin_changed(item, rec):
        return dataclasses.replace(rec, status="mismatch", message=(
            f"installed {rec.version}, this release pins {item.version}"))  # fmt: skip
    prev = old.get(item.name)
    if rec.status == "ok" or item.check != "antismash" or prev is None:
        return rec
    if (prev.status, prev.url) != ("ok", "image"):
        return rec
    image = _antismash_image(item, runner)
    if image is None:
        return dataclasses.replace(rec, message="antismash --check-prereqs fails in the image")
    return dataclasses.replace(prev, message=image.message)


def run_setup(opts: SetupOptions, runner: CommandRunner, fetcher: Fetcher) -> ExitCode:
    """Verify (``--check``) or adopt/download every planned item and update the manifest."""
    items = plan(opts)
    if not items:
        return ExitCode.SUCCESS
    db = opts.db_dir
    if opts.check:
        old = read_records(db)
        checked = [_check(i, db, old, runner) for i in items]
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
            _say(f"another gippa setup is running on {db} ({LOCK}); try again later")
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
        write_manifest(db, recs, setup={"set": "minimal" if opts.minimal else "default"})
    return ExitCode.SUCCESS if status_ok else ExitCode.ENVIRONMENT


def setup_nextflow(db_dir: Path, fetcher: Fetcher, environ: Mapping[str, str]) -> ExitCode:
    """The Nextflow launcher (M11-15 item 5): a ``nextflow`` on ``PATH`` is reported (version and
    path) and nothing is downloaded; else the pinned ``-dist`` launcher is adopted or fetched into
    ``<db_dir>/bin/nextflow`` (mode 0755) and recorded in ``db_manifest.json``."""
    from gippa.runner.prereqs import nextflow_version

    found = shutil.which("nextflow", path=environ.get("PATH", ""))
    if found is not None:
        version = nextflow_version(Path(found), environ) or "of unknown version"
        _say(f"nextflow: found {version} at {found}; the launcher is not downloaded")
        return ExitCode.SUCCESS
    item = SOURCES["nextflow"]
    db_dir.mkdir(parents=True, exist_ok=True)
    with _locked(db_dir) as acquired:
        if not acquired:
            _say(f"another gippa setup is running on {db_dir} ({LOCK}); try again later")
            return ExitCode.ENVIRONMENT
        old = read_records(db_dir)
        rec = install(item, db_dir, ToolRunner(db_dir=db_dir), fetcher)
        if rec.status == "ok":
            (db_dir / item.path).chmod(0o755)
            if not rec.message:  # adopted, not downloaded now
                rec = _provenance(rec, old, {}, None)
        _report(rec)
        write_manifest(db_dir, [rec])
    return ExitCode.SUCCESS if rec.status == "ok" else ExitCode.ENVIRONMENT
