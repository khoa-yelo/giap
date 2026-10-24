"""Download, verify and unpack the PXD078459 runs a manifest lists (M10-01, R6, R11, §3.9).

Stdlib only; ``scripts/slurm/fetch_pxd078459.sbatch`` runs it as ``python3 -I``. Downloads go
through ``curl`` (https or ftps only, TLS always verified); nothing is written outside ``--dest``.

    pxd_fetch.py https DRAFT --out MANIFEST     # ftp:// rows -> same-host https when HEAD agrees
    pxd_fetch.py download --dest DIR --manifest MANIFEST [--files "A B C" --base-url URL]

Exit codes: 0 done; 1 a download, verification or unpack failed (after one retry); 2 a
precondition (no file list, scheme, disk, format, an https route not confirmed).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import zipfile
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from functools import partial
from pathlib import Path, PurePosixPath
from typing import IO

COLUMNS = (
    "run",
    "file_name",
    "url",
    "url_record",
    "bytes",
    "checksum_record",
    "sha256",
    "sample",
    "selected_because",
)
NA = "NA"
ASSUMED_RUN_BYTES = 5 * 10**9  # spec: 5 GB per run when the record gives no size
TLS_SCHEMES = ("https://", "ftps://")
ARCHIVES = (".zip", ".tar", ".tar.gz", ".tgz")
CURL_TLS = ["--proto", "=https,ftps", "--proto-redir", "=https,ftps"]
FILLED_NAME = "pxd078459_files.tsv"

Row = dict[str, str]


class Unsafe(Exception):
    """An archive member that could write outside its run folder."""


class Precondition(Exception):
    """Exit 2: nothing is downloaded."""


# ------------------------------------------------------------------ manifest


def read_manifest(path: Path) -> list[Row]:
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader, None)
        if header is None or tuple(header) != COLUMNS:
            raise Precondition(f"{path}: the header is not {' '.join(COLUMNS)}")
        rows = [dict(zip(COLUMNS, r, strict=True)) for r in reader if r]
    return rows


def write_manifest(path: Path, rows: Sequence[Row]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n",
                           quoting=csv.QUOTE_NONE)  # fmt: skip
        w.writeheader()
        w.writerows(rows)


def run_label(name: str) -> str:
    parts = PurePosixPath(name).parts
    if len(parts) > 1 and parts[0].endswith(".d"):  # a file inside a run folder: <run>.d/<file>
        return parts[0].removesuffix(".d")
    stem = parts[-1]
    for ext in (*ARCHIVES, ".d"):
        stem = stem.removesuffix(ext)
    return stem.removesuffix(".d")


def override_rows(files: str, base_url: str) -> list[Row]:
    """``PXD078459_FILES`` (space separated, never globbed) under ``PXD078459_BASE_URL``."""
    if not base_url:
        raise Precondition("PXD078459_FILES is set without PXD078459_BASE_URL")
    rows = []
    for name in files.split():
        url = f"{base_url.rstrip('/')}/{name}"
        rows.append({**dict.fromkeys(COLUMNS, NA), "run": run_label(name), "file_name": name,
                     "url": url, "url_record": url,
                     "selected_because": "override:PXD078459_FILES"})  # fmt: skip
    return rows


def kind_of(name: str) -> str:
    """``archive`` or ``member`` (a file inside a ``<run>.d/`` folder); else a precondition."""
    if name.lower().endswith(ARCHIVES):
        return "archive"
    if re.match(r"^[^/]+\.d/.+", name, re.I):
        return "member"
    raise Precondition(
        f"unsupported archive format: {name} (zip, tar, tar.gz only; question for the lead)"
    )


def check_rows(rows: Sequence[Row]) -> None:
    if not rows:
        raise Precondition("the file list is empty")
    for r in rows:
        name, url = r["file_name"], r["url"]
        if not url.startswith(TLS_SCHEMES):
            raise Precondition(
                f"{name}: {url}: R6 / M3-24 decision 11: plain FTP or HTTP is refused"
            )
        parts = PurePosixPath(name).parts
        if name.startswith("/") or ".." in parts or "\\" in name:
            raise Precondition(f"{name}: unsafe file name in the manifest")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", r["run"]):
            raise Precondition(f"{name}: unsafe run label {r['run']!r}")
        kind_of(name)


# ------------------------------------------------------------------ files and checksums


def digest(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _int(value: str) -> int | None:
    return int(value) if value.isdigit() else None


def verified(path: Path, row: Row) -> bool:
    """Size and record checksum match whenever the record states them."""
    if not path.is_file():
        return False
    size = _int(row["bytes"])
    if size is not None and path.stat().st_size != size:
        return False
    if row["checksum_record"] != NA:
        algorithm, _, want = row["checksum_record"].partition(":")
        return digest(path, algorithm) == want.lower()
    return True


def target_of(dest: Path, row: Row) -> Path:
    if kind_of(row["file_name"]) == "archive":
        return dest / row["file_name"]
    return dest / f"{row['run']}.d" / row["file_name"].split("/", 1)[1]


def unpacked(dest: Path, row: Row) -> bool:
    run_d = dest / f"{row['run']}.d"
    return (
        kind_of(row["file_name"]) == "archive"
        and row["sha256"] != NA
        and (run_d / "analysis.tdf").is_file()
        and (run_d / "analysis.tdf_bin").is_file()
    )


def free_bytes(path: Path) -> int:
    out = subprocess.run(["df", "-Pk", str(path)], capture_output=True, text=True, check=True)
    return int(out.stdout.strip().splitlines()[-1].split()[3]) * 1024


def curl_get(url: str, target: Path) -> int:
    target.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["curl", "--fail", "--location", "--silent", "--show-error", *CURL_TLS,
           "--continue-at", "-", "--max-time", "86400", "--output", str(target), url]  # fmt: skip
    return subprocess.run(cmd, check=False).returncode


def fetch(row: Row, target: Path, get: Callable[[str, Path], int]) -> str | None:
    """Download ``row`` to ``target`` (resuming a partial file); an error message or None."""
    name = row["file_name"]
    known = row["bytes"] != NA or row["checksum_record"] != NA
    if known and verified(target, row):
        print(f"already there, skipping: {name}")
        return None
    error = ""
    for attempt in (1, 2):
        size = _int(row["bytes"])
        if target.is_file() and size is not None and target.stat().st_size > size:
            target.unlink()
        print(f"downloading {row['url']}" + (" (retry)" if attempt == 2 else ""))
        rc = get(row["url"], target)
        if rc != 0:
            error = f"download failed after one retry: {name} (curl exit {rc})"
        elif verified(target, row):
            return None
        else:
            target.unlink()
            error = f"checksum mismatch after one retry: {name}"
            print(f"checksum mismatch: {name}", file=sys.stderr)
    return error


# ------------------------------------------------------------------ archives


@dataclass(frozen=True)
class Member:
    name: str
    size: int
    is_dir: bool
    open: Callable[[], io.BufferedIOBase | IO[bytes] | None]


def _check_name(name: str) -> None:
    n = name.replace("\\", "/")
    if n.startswith("/") or re.match(r"^[A-Za-z]:", n) or ".." in n.split("/"):
        raise Unsafe(f"unsafe archive member: {name}")


def _zip_members(zf: zipfile.ZipFile) -> list[Member]:
    infos = zf.infolist()
    for i in infos:
        _check_name(i.filename)
        if stat.S_ISLNK(i.external_attr >> 16):
            raise Unsafe(f"unsafe archive member: {i.filename} (symbolic link)")
    bad = zf.testzip()  # the CRC test of every member, as ``unzip -t`` does
    if bad is not None:
        raise ValueError(f"corrupt zip member: {bad}")
    return [Member(i.filename, i.file_size, i.is_dir(), partial(zf.open, i)) for i in infos]


def _tar_members(tf: tarfile.TarFile) -> list[Member]:
    infos = tf.getmembers()
    for t in infos:
        _check_name(t.name)
        if not (t.isfile() or t.isdir()):
            raise Unsafe(f"unsafe archive member: {t.name} (link or special file)")
    return [Member(t.name, t.size, t.isdir(), partial(tf.extractfile, t)) for t in infos]


def unpack(archive: Path, run: str, dest: Path, max_bytes: int) -> str:
    """Unpack the one ``.d`` folder of ``archive`` to ``dest/<run>.d``; its archive folder name.

    Every member is checked (absolute, ``..``, links, special files) before the first byte is
    written; the unpacked size must fit in ``max_bytes``.
    """
    if archive.name.lower().endswith(".zip"):
        if not zipfile.is_zipfile(archive):
            raise ValueError(f"not a readable zip: {archive.name}")
        with zipfile.ZipFile(archive) as zf:
            return _extract(archive, _zip_members(zf), run, dest, max_bytes)
    with tarfile.open(archive, "r:*") as tf:
        return _extract(archive, _tar_members(tf), run, dest, max_bytes)


def _extract(archive: Path, found: list[Member], run: str, dest: Path, max_bytes: int) -> str:
    files = [m for m in found if not m.is_dir]
    tdfs = [m for m in files if PurePosixPath(m.name).name == "analysis.tdf"]
    if len(tdfs) != 1:
        raise ValueError(f"{archive.name}: {len(tdfs)} analysis.tdf files, expected exactly one")
    root = PurePosixPath(tdfs[0].name).parent
    top = root == PurePosixPath(".")
    keep = [m for m in files if top or root in PurePosixPath(m.name).parents]
    total = sum(m.size for m in keep)
    if total > max_bytes:
        raise ValueError(f"{archive.name} unpacks to {total} bytes, more than {max_bytes} free")
    tmp = dest / f".{run}.d.partial"
    shutil.rmtree(tmp, ignore_errors=True)
    for m in keep:
        out = (tmp / PurePosixPath(m.name).relative_to(root)).resolve()
        if tmp.resolve() not in out.parents:
            raise Unsafe(f"unsafe archive member: {m.name}")
        out.parent.mkdir(parents=True, exist_ok=True)
        src = m.open()
        if src is None:
            raise ValueError(f"unreadable archive member: {m.name}")
        with src, out.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    final = dest / f"{run}.d"
    shutil.rmtree(final, ignore_errors=True)
    tmp.rename(final)
    source = root.name or archive.name
    (dest / f"{run}.d.source_folder").write_text(source + "\n")
    return source


# ------------------------------------------------------------------ commands


def download(
    rows: list[Row],
    dest: Path,
    *,
    get: Callable[[str, Path], int] = curl_get,
    free: Callable[[Path], int] = free_bytes,
) -> int:
    check_rows(rows)
    dest.mkdir(parents=True, exist_ok=True)
    todo = [r for r in rows if not unpacked(dest, r)]
    fetched = {id(r) for r in todo}
    need = 3 * sum(_int(r["bytes"]) or ASSUMED_RUN_BYTES for r in todo)
    avail = free(dest)
    if avail < need:
        raise Precondition(f"{avail} bytes free under {dest}, need {need} (3 x the download)")
    errors, failed = [], set()
    for r in rows:
        if id(r) in fetched:
            target = target_of(dest, r)
            err = fetch(r, target, get)
            if err:
                errors.append(err)
                failed.add(id(r))
                print(err, file=sys.stderr)
                continue
            r["sha256"] = digest(target, "sha256")
        else:
            print(f"already unpacked, skipping: {r['run']}.d ({r['file_name']})")
    done = [r for r in rows if r["sha256"] != NA and id(r) not in failed]
    sums = "".join(f"{r['sha256']}  {r['file_name']}\n" for r in done)
    (dest / "SHA256SUMS").write_text(sums)
    write_manifest(dest / FILLED_NAME, rows)
    for r in done:
        archive = dest / r["file_name"]
        if id(r) not in fetched or kind_of(r["file_name"]) != "archive":
            continue
        try:
            source = unpack(archive, r["run"], dest, free(dest))
        except (Unsafe, ValueError, tarfile.TarError, zipfile.BadZipFile, OSError) as err:
            errors.append(str(err))
            print(err, file=sys.stderr)
            shutil.rmtree(dest / f".{r['run']}.d.partial", ignore_errors=True)
            continue
        archive.unlink()
        print(f"unpacked {r['file_name']} ({source}) -> {r['run']}.d")
    return 1 if errors else 0


def curl_head(url: str) -> tuple[str, str]:
    cmd = ["curl", "--head", "--silent", "--location", *CURL_TLS, "--output", "/dev/null",
           "--write-out", "%{http_code} %header{content-length}", url]  # fmt: skip
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    code, _, length = res.stdout.strip().partition(" ")
    return code, length.strip()


def https_route(rows: list[Row], head: Callable[[str], tuple[str, str]] = curl_head) -> list[str]:
    """Q1: an ftp:// row downloads over https on the same host and path when HEAD gives 200
    and the record's size; the rows without a confirmed route are returned."""
    refused = []
    for r in rows:
        if r["url"].startswith(TLS_SCHEMES):
            continue
        if r["url"].startswith("ftp://") and r["bytes"] != NA:
            https = "https://" + r["url"].removeprefix("ftp://")
            code, length = head(https)
            if code == "200" and length == r["bytes"]:
                r["url"] = https
                continue
        refused.append(r["file_name"])
    return refused


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("https")
    h.add_argument("draft", type=Path)
    h.add_argument("--out", type=Path, required=True)
    d = sub.add_parser("download")
    d.add_argument("--dest", type=Path, required=True)
    d.add_argument("--manifest", type=Path, required=True)
    d.add_argument("--files", default="")
    d.add_argument("--base-url", default="")
    args = ap.parse_args(argv)
    try:
        if args.cmd == "https":
            rows = read_manifest(args.draft)
            refused = https_route(rows)
            write_manifest(args.out, rows)
            if refused:
                print(f"no confirmed https route (R6 / M3-24 decision 11; Q1): "
                      f"{', '.join(refused)}", file=sys.stderr)  # fmt: skip
                return 2
            return 0
        if args.files.strip():
            rows = override_rows(args.files, args.base_url)
        elif args.manifest.is_file():
            rows = read_manifest(args.manifest)
            if args.base_url:
                for r in rows:
                    r["url"] = f"{args.base_url.rstrip('/')}/{r['file_name']}"
        else:
            raise Precondition("no file list: run --resolve first")
        return download(rows, args.dest)
    except Precondition as err:
        print(err, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
