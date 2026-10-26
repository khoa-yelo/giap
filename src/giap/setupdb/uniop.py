"""``giap setup --uniop``: fetch the authors' UniOP archive, verify it, build it locally (M11-14
item 1b, 14u1).

UniOP has no licence file, so GIAP never distributes it (R5). Setup fetches the pinned archive
(``db_sources.yaml`` ``uniop_source``) from GitHub only, following its one allowlisted redirect to
codeload (R6), checks its sha256 before anything is built, stages a build context
``<cache>/uniop/context/`` from the bundled files and the archive, and builds the image under the
runtime the profile list resolves to (``R(P)``, :func:`giap.preflight.environment.uniop_runtime`).
The built image is named in ``<db_dir>/local_tools.config`` (one ``profiles`` block per runtime)
and recorded as the ``image:uniop`` item of ``db_manifest.json``. Nothing is ever pushed.

14u1 builds with docker; the apptainer and micromamba builders are 14u2.
"""

from __future__ import annotations

import contextlib
import dataclasses
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Literal

from giap.exitcodes import ExitCode
from giap.io.image_key import UNIOP_ARCHIVE, key, key_inputs, tag
from giap.preflight.environment import (
    LOCAL_TOOLS,
    UNIOP_SETUP_ITEM,
    _blocks,
    _strip_comments,
    uniop_runtime,
    uniop_runtime_problem,
)
from giap.setupdb import items
from giap.setupdb.items import DbItem, DbRecord, Status, check_file, file_digest, today
from giap.setupdb.manifest import write_manifest
from giap.setupdb.run import LOCK, Fetcher, SourceError, _absent, _locked, _report, _say

__all__ = [
    "TERMS_NOT_ACCEPTED",
    "docker_build_argv",
    "fetch_archive",
    "local_tools_text",
    "setup_uniop",
    "stage_context",
    "uniop_ref",
]

Terms = Literal["flag", "prompt", ""]
SOURCE_ITEM = "uniop_source"
TERMS_NOT_ACCEPTED = "terms not accepted; rerun with --accept-uniop-terms"
DOCKERFILE = "workflow/containers/uniop/Dockerfile"
UNIOP_DEF = "workflow/containers/uniop/uniop.def"
_PKG = __name__.split(".", 1)[0]
_IMAGE_ID = re.compile(r"sha256:[0-9a-f]{64}")


def _source() -> DbItem:
    """The pinned archive row, read at call time."""
    return items.SOURCES[SOURCE_ITEM]


def _workflow_root() -> Path:
    """The directory holding the bundled ``workflow/`` (checkout or installed package)."""
    from giap.runner.nextflow import workflow_main

    return workflow_main().parent.parent


def _private_dir(path: Path) -> None:
    """``path`` (and its parents) exists and is ``0700``: never a lab-wide copy of UniOP."""
    path.mkdir(parents=True, exist_ok=True)
    path.chmod(0o700)


# ---------------------------------------------------------------- pure parts


def uniop_ref(version: str, image_key: str) -> str:
    """``localhost/<pkg>-uniop:<version tag>-<key>``: a local name, never a registry one."""
    return f"localhost/{_PKG}-uniop:{tag(version, image_key)}"


def docker_build_argv(context: Path, ref: str) -> list[str]:
    """The docker build of spec item 1b, exactly."""
    return ["docker", "build", "-f", f"{context}/Dockerfile", "-t", ref, str(context)]


def local_tools_text(existing: str, runtime: str, container: str) -> str:
    """``local_tools.config`` with the ``runtime`` block naming ``container`` (``conda = ''``),
    every other runtime's block of ``existing`` kept. Never sets ``conda.enabled``."""
    kept: list[tuple[str, str]] = []
    for top, body in _blocks(_strip_comments(existing))[0]:
        if top == "profiles":
            kept += [(n, b.strip()) for n, b in _blocks(body)[0] if n and n != runtime]
    line = f"process {{ withName: 'UNIOP' {{ container = '{container}'; conda = '' }} }}"
    blocks = sorted([*kept, (runtime, line)])
    head = "// Written by giap setup --uniop (M11-14 item 1b): one block per runtime built.\n"
    return head + "profiles {\n" + "".join(f"    {n} {{ {b} }}\n" for n, b in blocks) + "}\n"


def _record(item: DbItem, status: Status, message: str, checksum: str = "") -> DbRecord:
    """The ``image:uniop`` item (``path`` ``""``: docker keeps the image in its own store)."""
    ok = status == "ok"
    return DbRecord(UNIOP_SETUP_ITEM, item.version, item.url, checksum or item.checksum, 0,
                    today() if ok else "", "", status, message)  # fmt: skip


def _last_line(text: str) -> str:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    return lines[-1] if lines else "no output"


# ---------------------------------------------------------------- fetch and stage


def fetch_archive(item: DbItem, db_dir: Path, fetcher: Fetcher) -> DbRecord:
    """The verified archive at ``<db_dir>/<item.path>`` (``0600`` in a ``0700`` directory): a
    copy already in place that verifies is adopted, else one download from the pinned URL.
    A sha256 mismatch deletes the bytes (``uniop: archive sha256 <got> != pinned <want>``)."""
    target = db_dir / item.path
    _private_dir(target.parent)
    rec = check_file(item, target, item.checksum, item.version)
    if rec.status == "ok":
        target.chmod(0o600)
        return rec
    target.unlink(missing_ok=True)
    partial = target.with_name(target.name + ".partial")
    try:
        fetcher.fetch(item.url, partial)
        got = file_digest(partial, "sha256")
        if got != item.checksum:
            return _absent(item, f"uniop: archive sha256 {got} != pinned {item.checksum}",
                           "mismatch")  # fmt: skip
        partial.chmod(0o600)
        partial.replace(target)
    except SourceError as exc:
        return _absent(item, f"source {exc.host} unreachable ({exc.reason})")
    finally:
        partial.unlink(missing_ok=True)
    return dataclasses.replace(check_file(item, target, item.checksum, item.version),
                               message=f"downloaded from {item.host}")  # fmt: skip


def stage_context(context: Path, archive: Path, root: Path) -> None:
    """``context`` holds the key inputs at their repo paths, the ``Dockerfile`` and ``uniop.def``
    at its root, and the verified archive as ``uniop-src.tar.gz`` (``0600``)."""
    if context.exists():
        shutil.rmtree(context)
    context.mkdir(parents=True)
    for rel in key_inputs("uniop", root):
        dest = context / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / rel, dest)
    shutil.copyfile(root / DOCKERFILE, context / "Dockerfile")
    shutil.copyfile(root / UNIOP_DEF, context / "uniop.def")
    shutil.copyfile(archive, context / UNIOP_ARCHIVE)
    (context / UNIOP_ARCHIVE).chmod(0o600)


# ---------------------------------------------------------------- docker


def _docker_build(context: Path, ref: str) -> tuple[str, str]:
    """``(image ID, "")`` after the build and ``docker image inspect``, else ``("", why)``."""
    docker = shutil.which("docker")
    if docker is None:
        return "", "docker not found on PATH"
    argv = [docker, *docker_build_argv(context, ref)[1:]]
    try:
        proc = subprocess.run(argv, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            return "", f"docker build failed: {_last_line(proc.stderr or proc.stdout)}"
        inspect = subprocess.run([docker, "image", "inspect", "--format", "{{.Id}}", ref],
                                 capture_output=True, text=True, check=False)  # fmt: skip
    except OSError as exc:
        return "", f"docker build failed: {exc}"
    image_id = inspect.stdout.strip()
    if inspect.returncode != 0 or not _IMAGE_ID.fullmatch(image_id):
        why = _last_line(inspect.stderr or inspect.stdout)
        return "", f"docker image inspect {ref} failed: {why}"
    return image_id, ""


def _write_local_tools(db_dir: Path, runtime: str, container: str) -> None:
    path = db_dir / LOCAL_TOOLS
    try:
        existing = path.read_text(encoding="utf-8")
    except OSError:
        existing = ""
    fd, tmp = tempfile.mkstemp(prefix=f".{LOCAL_TOOLS}.", suffix=".tmp", dir=db_dir)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(local_tools_text(existing, runtime, container))
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


# ---------------------------------------------------------------- the step


def _build(item: DbItem, db_dir: Path, cache: Path, fetcher: Fetcher,
           terms: Terms) -> tuple[list[DbRecord], dict[str, str] | None]:  # fmt: skip
    """Fetch, verify, stage and build under docker; the records to write and the
    ``local_build`` object of a successful build."""
    src = fetch_archive(item, db_dir, fetcher)
    _report(src)
    if src.status != "ok":
        return [src, _record(item, src.status, src.message)], None
    image_key = key("uniop")
    ref = uniop_ref(item.version, image_key)
    context = cache / "uniop" / "context"
    _private_dir(cache / "uniop")
    stage_context(context, db_dir / item.path, _workflow_root())
    image_id, why = _docker_build(context, ref)
    if why:
        return [src, _record(item, "absent", why, src.checksum)], None
    _write_local_tools(db_dir, "docker", ref)
    build = {"runtime": "docker", "ref": ref, "digest": image_id, "key": image_key,
             "terms": terms}  # fmt: skip
    return [src, _record(item, "ok", "", src.checksum)], build


def setup_uniop(profile: str, db_dir: Path, cache: Path, fetcher: Fetcher, *, accepted: bool,
                explicit: bool) -> ExitCode:  # fmt: skip
    """The UniOP step of ``setup --uniop`` (``explicit``) or ``setup --images``. Without
    accepted terms nothing is fetched: exit 2 under ``--uniop``, a warning and exit 0 under
    ``--images``. Once accepted, a failed fetch, checksum or build is exit 2."""
    item = _source()
    runtime = uniop_runtime(profile)
    db_dir.mkdir(parents=True, exist_ok=True)
    with _locked(db_dir) as acquired:
        if not acquired:
            _say(f"another giap setup is running on {db_dir} ({LOCK}); try again later")
            return ExitCode.ENVIRONMENT
        code = ExitCode.ENVIRONMENT
        build: dict[str, str] | None = None
        if runtime is None:
            recs = [_record(item, "absent", f"uniop: {uniop_runtime_problem(profile)}")]
        elif not accepted:
            recs = [_record(item, "absent", TERMS_NOT_ACCEPTED)]
            code = ExitCode.ENVIRONMENT if explicit else ExitCode.SUCCESS
        elif runtime != "docker":
            recs = [_record(item, "absent", f"uniop: the local {runtime} build is not available "
                                            "in this version; use --profile docker")]  # fmt: skip
        else:
            recs, build = _build(item, db_dir, cache, fetcher, "flag")
            code = ExitCode.SUCCESS if recs[-1].status == "ok" else ExitCode.ENVIRONMENT
        image = recs[-1]
        if image.status == "ok":
            _report(image)
        elif code == ExitCode.SUCCESS:
            _say(f"warning: {image.name}: {image.message}")
        else:
            _say(f"{image.name}: {image.message}")
        write_manifest(db_dir, recs, {UNIOP_SETUP_ITEM: build} if build is not None else None)
    return code
