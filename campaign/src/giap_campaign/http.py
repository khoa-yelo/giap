"""https-only JSON GET with cache, checked downloads and safe archive extraction (campaign SPEC §2).

The only module (with ``driver``) that touches the network; every call goes through an injectable
``opener`` / ``run`` so tests stay offline. ``NCBI_API_KEY`` is sent as header ``api-key`` only to
``api.ncbi.nlm.nih.gov`` and is never logged or written (SPEC §1).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import tarfile
import urllib.error
import urllib.request
import zipfile
from collections.abc import Callable
from pathlib import Path, PurePosixPath
from typing import IO, Any
from urllib.parse import urlsplit

__all__ = ["DownloadError", "HttpError", "NotFound", "download", "get_json", "safe_extract"]

NCBI_API_HOST = "api.ncbi.nlm.nih.gov"
_HEX = re.compile(r"[0-9a-fA-F]+")
_DIGEST_BY_LENGTH = {32: "md5", 40: "sha1", 64: "sha256"}

Opener = Callable[..., Any]
Runner = Callable[..., "subprocess.CompletedProcess[Any]"]


class NotFound(LookupError):
    """HTTP 404: the record does not exist (never retried)."""


class HttpError(OSError):
    """A failed HTTP request (4xx other than 404, or still failing after the retries)."""

    def __init__(self, message: str, status: int | None = None) -> None:
        super().__init__(message)
        self.status = status


class DownloadError(OSError):
    """A download that failed, did not match its size or checksum, or an unsafe archive."""


def _require_https(url: str) -> None:
    parts = urlsplit(url)
    if parts.scheme.lower() != "https" or not parts.hostname:
        raise ValueError(f"only https:// URLs are allowed: {url!r}")


def _request(url: str) -> urllib.request.Request:
    headers = {"Accept": "application/json"}
    key = os.environ.get("NCBI_API_KEY")
    if key and urlsplit(url).hostname == NCBI_API_HOST:
        headers["api-key"] = key
    return urllib.request.Request(url, headers=headers)


def get_json(
    url: str,
    *,
    cache_dir: str | Path | None = None,
    timeout: float = 60,
    retries: int = 1,
    opener: Opener | None = None,
) -> Any:
    """Parsed JSON of ``url``; with ``cache_dir`` stored as ``<sha256(url)>.json`` and reused."""
    _require_https(url)
    cached: Path | None = None
    if cache_dir is not None:
        cached = Path(cache_dir) / f"{hashlib.sha256(url.encode()).hexdigest()}.json"
        if cached.is_file():
            return json.loads(cached.read_bytes())
    open_ = opener if opener is not None else urllib.request.urlopen
    last: BaseException | None = None
    for _ in range(max(0, retries) + 1):
        try:
            with open_(_request(url), timeout=timeout) as response:
                body: bytes = response.read()
        except urllib.error.HTTPError as exc:  # before URLError: HTTPError subclasses it
            if exc.code == 404:
                raise NotFound(f"not found: {url}") from exc
            if exc.code < 500:
                raise HttpError(f"HTTP {exc.code} for {url}", exc.code) from exc
            last = exc
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
        else:
            data = json.loads(body)  # parse first: errors are never cached
            if cached is not None:
                cached.parent.mkdir(parents=True, exist_ok=True)
                tmp = cached.with_name(cached.name + ".tmp")
                tmp.write_bytes(body)
                tmp.replace(cached)
            return data
    status = last.code if isinstance(last, urllib.error.HTTPError) else None
    raise HttpError(f"request failed after {retries} retries: {url}: {last}", status) from last


def _digest_name(checksum: str | None) -> str | None:
    """hashlib name for a 32/40/64-hex checksum; None for anything else (then ignored)."""
    if not checksum or not _HEX.fullmatch(checksum):
        return None
    return _DIGEST_BY_LENGTH.get(len(checksum))


def _check(path: Path, expected_bytes: int | None, checksum: str | None) -> tuple[bool, str]:
    """(passes size and checksum checks, SHA-256 of the file)."""
    if expected_bytes is not None and path.stat().st_size != expected_bytes:
        return False, ""
    algo = _digest_name(checksum)
    sha256 = hashlib.sha256()
    other = hashlib.new(algo) if algo is not None and algo != "sha256" else None
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            sha256.update(chunk)
            if other is not None:
                other.update(chunk)
    digest = sha256.hexdigest()
    if algo is not None and checksum is not None:
        got = digest if other is None else other.hexdigest()
        return got == checksum.lower(), digest
    return True, digest


def download(
    url: str,
    dest: str | Path,
    *,
    expected_bytes: int | None = None,
    checksum: str | None = None,
    run: Runner = subprocess.run,
) -> str:
    """Fetch ``url`` to ``dest`` with curl (https only), check, rename; returns its SHA-256."""
    _require_https(url)
    dest = Path(dest)
    if dest.is_file():
        ok, digest = _check(dest, expected_bytes, checksum)
        if ok:
            return digest
    part = Path(str(dest) + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["curl", "-fsSL", "--proto", "=https", "--proto-redir", "=https", "--retry", "1"]
    cmd += ["-o", str(part), url]
    proc = run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        part.unlink(missing_ok=True)
        detail = (proc.stderr or "").strip()[-500:]
        raise DownloadError(f"curl exited {proc.returncode} for {url}: {detail}")
    if not part.is_file():
        raise DownloadError(f"curl wrote no file for {url}")
    ok, digest = _check(part, expected_bytes, checksum)
    if not ok:
        part.unlink(missing_ok=True)
        raise DownloadError(f"size or checksum mismatch for {url}")
    part.replace(dest)
    return digest


def _safe_target(out_dir: Path, name: str) -> Path:
    """``out_dir / name``, refusing absolute paths and ``..`` (also after resolving symlinks)."""
    member = PurePosixPath(name)
    if member.is_absolute() or ".." in member.parts or not member.parts:
        raise DownloadError(f"unsafe archive member: {name!r}")
    target = out_dir / Path(*member.parts)
    if not target.resolve().is_relative_to(out_dir.resolve()):
        raise DownloadError(f"archive member escapes the output directory: {name!r}")
    return target


def _write(target: Path, source: IO[bytes]) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as fh:
        shutil.copyfileobj(source, fh)


def _extract_zip(archive: Path, out_dir: Path) -> list[Path]:
    with zipfile.ZipFile(archive) as zf:
        members = zf.infolist()
        plan: list[tuple[zipfile.ZipInfo, Path]] = []
        for info in members:  # validate everything before writing anything
            if stat.S_ISLNK(info.external_attr >> 16):
                raise DownloadError(f"symlink in archive: {info.filename!r}")
            plan.append((info, _safe_target(out_dir, info.filename)))
        files: list[Path] = []
        for info, target in plan:
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            with zf.open(info) as src:
                _write(target, src)
            files.append(target)
    return files


def _extract_tar(archive: Path, out_dir: Path) -> list[Path]:
    with tarfile.open(archive, "r:*") as tf:
        plan: list[tuple[tarfile.TarInfo, Path]] = []
        for info in tf.getmembers():
            if info.issym() or info.islnk():
                raise DownloadError(f"link in archive: {info.name!r}")
            if not (info.isfile() or info.isdir()):
                raise DownloadError(f"special file in archive: {info.name!r}")
            plan.append((info, _safe_target(out_dir, info.name)))
        files: list[Path] = []
        for info, target in plan:
            if info.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            src = tf.extractfile(info)
            if src is None:
                raise DownloadError(f"unreadable archive member: {info.name!r}")
            with src:
                _write(target, src)
            files.append(target)
    return files


def safe_extract(archive: str | Path, out_dir: str | Path) -> list[Path]:
    """Extract a ``.zip``/``.tar``/``.tar.gz``/``.tgz``, refusing unsafe members; files written."""
    archive, out_dir = Path(archive), Path(out_dir)
    name = archive.name.lower()
    out_dir.mkdir(parents=True, exist_ok=True)
    if name.endswith(".zip"):
        return _extract_zip(archive, out_dir)
    if name.endswith((".tar", ".tar.gz", ".tgz")):
        return _extract_tar(archive, out_dir)
    raise ValueError(f"unsupported archive type: {archive.name}")
