"""``gippa setup --uniop``: fetch the authors' UniOP archive, verify it, build it locally (M11-14
item 1b, 14u1).

UniOP has no licence file, so GIPPA never distributes it (R5). Setup fetches the pinned archive
(``db_sources.yaml`` ``uniop_source``) from GitHub only, following its one allowlisted redirect to
codeload (R6), checks its sha256 before anything is built, stages a build context
``<cache>/uniop/context/`` from the bundled files and the archive, and builds the image under the
runtime the profile list resolves to (``R(P)``, :func:`gippa.preflight.environment.uniop_runtime`).
The built image or env is named in ``<db_dir>/local_tools.config`` (one ``profiles`` block per
runtime) and recorded as the ``image:uniop`` item of ``db_manifest.json`` with its ``local_build``
object (one per item: the last build wins). Nothing is ever pushed. Builders: docker (14u1),
apptainer with the fakeroot fallback, micromamba under plain ``conda``, and the check of the
prebuilt env under ``farmshare*`` (14u2). UniOP has no licence, so the terms notice is always
printed and nothing is fetched until the user accepts (``--accept-uniop-terms`` or the prompt).
"""

from __future__ import annotations

import contextlib
import dataclasses
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Literal

from gippa.exitcodes import ExitCode
from gippa.io.image_key import UNIOP_ARCHIVE, key, key_inputs, tag
from gippa.preflight.environment import (
    _BUILD_BLOCK,
    _FARMSHARE_PROFILES,
    LOCAL_TOOLS,
    UNIOP_SETUP_ITEM,
    _blocks,
    _profiles,
    _strip_comments,
    _uniop_build,
    uniop_local_problem,
    uniop_runtime,
    uniop_runtime_problem,
)
from gippa.setupdb import items
from gippa.setupdb.items import DbItem, DbRecord, Status, check_file, file_digest, today
from gippa.setupdb.manifest import write_manifest
from gippa.setupdb.run import (
    LOCK,
    Fetcher,
    SourceError,
    _absent,
    _locked,
    _report,
    _say,
    _unpack,
)
from gippa.setupdb.space import DISK_FULL, disk_full

__all__ = [
    "NOTICE",
    "TERMS_NOT_ACCEPTED",
    "apptainer_build_argvs",
    "docker_build_argv",
    "fetch_archive",
    "local_tools_text",
    "setup_uniop",
    "stage_context",
    "stdin_is_tty",
    "uniop_ref",
]

Terms = Literal["flag", "prompt", ""]
SOURCE_ITEM = "uniop_source"
TERMS_NOT_ACCEPTED = "terms not accepted; rerun with --accept-uniop-terms"
NOTICE = (
    "UniOP (github.com/hongsua/UniOP @3ed14dfc) is fetched from its authors and built on this "
    f"machine; it has no licence file, so {__name__.split('.', 1)[0]} does not distribute it and "
    "whether your use is permitted is yours to decide."
)
PROMPT = "Build UniOP for your own use? [y/N]"
APPTAINER_FAILED = (
    "apptainer build failed (fakeroot unavailable or build error): {line}; use --profile docker "
    "or conda, or set annotation.uniop.enabled=auto"
)
PREBUILT_HINT = "sbatch scripts/slurm/build_envs.sbatch uniop"
LOCK_FILE = "workflow/envs/locks/uniop.explicit.txt"
DOCKERFILE = "workflow/containers/uniop/Dockerfile"
UNIOP_DEF = "workflow/containers/uniop/uniop.def"
_PKG = __name__.split(".", 1)[0]
_IMAGE_ID = re.compile(r"sha256:[0-9a-f]{64}")


def _source() -> DbItem:
    """The pinned archive row, read at call time."""
    return items.SOURCES[SOURCE_ITEM]


def _workflow_root() -> Path:
    """The directory holding the bundled ``workflow/`` (checkout or installed package)."""
    from gippa.runner.nextflow import workflow_main

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


def local_tools_text(
    existing: str, runtime: str, container: str, drop: frozenset[str] = frozenset()
) -> str:
    """``local_tools.config`` with the ``runtime`` block (``docker``/``apptainer``: naming
    ``container`` with ``conda = ''``; ``conda``: ``conda = '<prefix>'``); the other blocks of
    ``existing`` are kept, except those in ``drop``. Never sets ``conda.enabled``."""
    kept: list[tuple[str, str]] = []
    for top, body in _blocks(_strip_comments(existing))[0]:
        if top == "profiles":
            kept += [(n, b.strip()) for n, b in _blocks(body)[0]
                     if n and n != runtime and n not in drop]  # fmt: skip
    inner = (f"conda = '{container}'" if runtime == "conda"
             else f"container = '{container}'; conda = ''")  # fmt: skip
    line = f"process {{ withName: 'UNIOP' {{ {inner} }} }}"
    blocks = sorted([*kept, (runtime, line)])
    head = "// Written by gippa setup --uniop (M11-14 item 1b): one block per runtime built.\n"
    return head + "profiles {\n" + "".join(f"    {n} {{ {b} }}\n" for n, b in blocks) + "}\n"


def apptainer_build_argvs(sif: Path, definition: Path, root: bool) -> list[list[str]]:
    """The apptainer build of spec item 1b: as root one plain build; otherwise ``--fakeroot``,
    then once ``--fakeroot --ignore-fakeroot-command`` (FarmShare's SUID apptainer)."""
    if root:
        return [["apptainer", "build", str(sif), str(definition)]]
    first = ["apptainer", "build", "--fakeroot", str(sif), str(definition)]
    return [first, [*first[:3], "--ignore-fakeroot-command", *first[3:]]]


def _record(
    item: DbItem, status: Status, message: str, checksum: str = "", path: str = ""
) -> DbRecord:
    """The ``image:uniop`` item; ``path``: the ``.sif`` or env prefix (``""`` for docker, which
    keeps the image in its own store)."""
    ok = status == "ok"
    return DbRecord(UNIOP_SETUP_ITEM, item.version, item.url, checksum or item.checksum, 0,
                    today() if ok else "", path, status, message)  # fmt: skip


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


def _apptainer_build(context: Path, sif: Path) -> tuple[str, str]:
    """``(sha256 digest of the .sif, "")`` after the build (fallback included), else
    ``("", why)``. Runs in ``context``: ``%files`` sources are relative to it."""
    exe = shutil.which("apptainer") or shutil.which("singularity")
    if exe is None:
        return "", APPTAINER_FAILED.format(line="apptainer not found on PATH")
    sif.unlink(missing_ok=True)
    why = "no output"
    for argv in apptainer_build_argvs(sif, context / "uniop.def", os.geteuid() == 0):
        try:
            proc = subprocess.run([exe, *argv[1:]], cwd=context, stdin=subprocess.DEVNULL,
                                  capture_output=True, text=True, check=False)  # fmt: skip
        except OSError as exc:
            why = str(exc)
            continue
        if proc.returncode == 0 and sif.is_file():
            sif.chmod(0o600)
            return file_digest(sif, "sha256"), ""
        why = _last_line(proc.stderr or proc.stdout)
        sif.unlink(missing_ok=True)
    return "", APPTAINER_FAILED.format(line=why)


def _micromamba_build(context: Path, prefix: Path, version: str) -> str:
    """``micromamba create`` of ``prefix`` from the lock, then the verified archive unpacked
    (first path component stripped) into ``<prefix>/share/uniop`` with ``GIPPA_COMMIT``; ``""``
    or why it failed."""
    exe = os.environ.get("MAMBA_EXE") or shutil.which("micromamba")
    if exe is None:
        return "micromamba not found on PATH"
    if prefix.exists():
        shutil.rmtree(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    argv = [exe, "create", "-y", "-p", str(prefix), "-f", str(context / LOCK_FILE)]
    try:
        proc = subprocess.run(argv, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                              check=False)  # fmt: skip
    except OSError as exc:
        return f"micromamba create failed: {exc}"
    if proc.returncode != 0:
        return f"micromamba create failed: {_last_line(proc.stderr or proc.stdout)}"
    prefix.chmod(0o700)
    staging = prefix / "share" / ".uniop-src"
    staging.mkdir(parents=True)
    unsafe = _unpack(context / UNIOP_ARCHIVE, staging)
    tops = [] if unsafe else list(staging.iterdir())
    if unsafe or len(tops) != 1 or not tops[0].is_dir():
        return f"uniop: {unsafe or 'the archive does not hold one top-level directory'}"
    home = tops[0].replace(prefix / "share" / "uniop")
    staging.rmdir()
    (home / "GIPPA_COMMIT").write_text(version.removeprefix("git:") + "\n", encoding="utf-8")
    return ""


def _write_local_tools(
    db_dir: Path, runtime: str, container: str, drop: frozenset[str] = frozenset()
) -> None:
    path = db_dir / LOCAL_TOOLS
    try:
        existing = path.read_text(encoding="utf-8")
    except OSError:
        existing = ""
    fd, tmp = tempfile.mkstemp(prefix=f".{LOCAL_TOOLS}.", suffix=".tmp", dir=db_dir)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(local_tools_text(existing, runtime, container, drop))
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


# ---------------------------------------------------------------- the step


def _build(item: DbItem, db_dir: Path, cache: Path, fetcher: Fetcher, runtime: str,
           terms: Terms) -> tuple[list[DbRecord], dict[str, str] | None]:  # fmt: skip
    """Fetch, verify, stage and build under ``runtime``; the records to write and the
    ``local_build`` object of a successful build."""
    cache = cache.absolute()  # the apptainer build runs in the context: no relative paths
    src = fetch_archive(item, db_dir, fetcher)
    _report(src)
    if src.status != "ok":
        return [src, _record(item, src.status, src.message)], None
    image_key = key("uniop")
    ref = uniop_ref(item.version, image_key)
    context = cache / "uniop" / "context"
    _private_dir(cache / "uniop")
    stage_context(context, db_dir / item.path, _workflow_root())
    path, digest, why = "", "", ""
    if runtime == "docker":
        digest, why = _docker_build(context, ref)
    elif runtime == "apptainer":
        ref = path = str(cache.absolute() / (ref.replace("/", "-").replace(":", "-") + ".sif"))
        digest, why = _apptainer_build(context, Path(path))
    else:
        ref = path = str(cache.absolute() / "envs" / "gippa-uniop")
        why = _micromamba_build(context, Path(path), item.version)
    if why:
        return [src, _record(item, "absent", why, src.checksum)], None
    old = _uniop_build(db_dir) or {}
    drop = frozenset({_BUILD_BLOCK.get(old.get("runtime", ""), "")} - {_BUILD_BLOCK[runtime], ""})
    _write_local_tools(db_dir, _BUILD_BLOCK[runtime], ref, drop)
    build = {"runtime": runtime, "ref": ref, "digest": digest, "key": image_key,
             "terms": terms}  # fmt: skip
    return [src, _record(item, "ok", "", src.checksum, path)], build


def _prebuilt(item: DbItem) -> tuple[list[DbRecord], dict[str, str] | None]:
    """``farmshare*``: the prebuilt ``$GIPPA_ROOT/micromamba/envs/gippa-uniop`` is checked and
    recorded (``prebuilt_envs.config`` owns UNIOP there); nothing is fetched or built."""
    root = os.environ.get("GIPPA_ROOT")
    prefix = Path(root) / "micromamba" / "envs" / "gippa-uniop" if root else None
    if prefix is None or not (prefix / "share" / "uniop").is_dir():
        where = f"{prefix} has no share/uniop" if prefix else "$GIPPA_ROOT is not set"
        return [_record(item, "absent", f"uniop: prebuilt env {where}; run: {PREBUILT_HINT}")], None
    build = {"runtime": "prebuilt", "ref": str(prefix), "digest": "", "key": key("uniop"),
             "terms": ""}  # fmt: skip
    return [_record(item, "ok", "", "", str(prefix))], build


def stdin_is_tty() -> bool:
    """Whether setup can ask (a terminal on stdin)."""
    try:
        return sys.stdin is not None and sys.stdin.isatty()
    except (AttributeError, ValueError):
        return False


def _terms(accepted: bool) -> Terms:
    """``flag`` with ``--accept-uniop-terms``; else, on a TTY, ``prompt`` when the user answers
    yes; ``""`` otherwise (the default is no; end of input is a no)."""
    if accepted:
        return "flag"
    if not stdin_is_tty():
        return ""
    print(f"{PROMPT} ", end="", file=sys.stderr, flush=True)
    answer = sys.stdin.readline()
    return "prompt" if answer.strip().lower() in ("y", "yes") else ""


def setup_uniop(profile: str, db_dir: Path, cache: Path, fetcher: Fetcher, *, accepted: bool,
                explicit: bool) -> ExitCode:  # fmt: skip
    """The UniOP step of ``setup --uniop`` (``explicit``) or ``setup --images``. A current
    build is kept as it is. Without accepted terms nothing is fetched: exit 2 under ``--uniop``,
    a warning and exit 0 under ``--images``; a good ``image:uniop`` record is then kept. Once
    accepted, a failed fetch, checksum or build is exit 2."""
    item = _source()
    runtime = uniop_runtime(profile)
    prebuilt = runtime == "micromamba" and bool(_profiles(profile) & _FARMSHARE_PROFILES)
    if runtime is not None and not prebuilt:
        print(f"gippa setup: {NOTICE}", flush=True)  # stdout: the one line a user must see
    db_dir.mkdir(parents=True, exist_ok=True)
    with _locked(db_dir) as acquired:
        if not acquired:
            _say(f"another gippa setup is running on {db_dir} ({LOCK}); try again later")
            return ExitCode.ENVIRONMENT
        code = ExitCode.ENVIRONMENT
        build: dict[str, str] | None = None
        recs: list[DbRecord] = []
        if runtime is None:
            recs = [_record(item, "absent", f"uniop: {uniop_runtime_problem(profile)}")]
        elif prebuilt:
            recs, build = _prebuilt(item)
        elif _uniop_build(db_dir) and uniop_local_problem(profile, db_dir / LOCAL_TOOLS) == "":
            _say(f"{UNIOP_SETUP_ITEM}: ok (the {runtime} build is current)")
            return ExitCode.SUCCESS
        elif not (terms := _terms(accepted)):
            code = ExitCode.ENVIRONMENT if explicit else ExitCode.SUCCESS
            _say(("warning: " if code == ExitCode.SUCCESS else "")
                 + f"{UNIOP_SETUP_ITEM}: {TERMS_NOT_ACCEPTED}")  # fmt: skip
            if _uniop_build(db_dir) is None:  # a good record is kept (#396 review S3)
                write_manifest(db_dir, [_record(item, "absent", TERMS_NOT_ACCEPTED)])
            return code
        else:
            try:
                recs, build = _build(item, db_dir, cache, fetcher, runtime, terms)
            except OSError as exc:  # e.g. a shared cache dir owned by another user (review S2)
                full = exc.errno in DISK_FULL  # our disk: failed, with the advice (15c1)
                at = str(exc.filename or cache).removesuffix(".partial")
                why = disk_full(at, " or --cache-dir") if full else f"uniop: {exc}"
                recs, build = [_record(item, "failed" if full else "absent", why)], None
        if recs[-1].status == "ok":
            code = ExitCode.SUCCESS
            _report(recs[-1])
        else:
            _say(f"{recs[-1].name}: {recs[-1].message}")
        write_manifest(db_dir, recs, {UNIOP_SETUP_ITEM: build} if build is not None else None)
    return code
