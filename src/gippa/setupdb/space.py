"""Disk space for ``gippa setup``: the plan, free space, user quotas, a full disk (M11-15 item 6).

Before any download setup prints one line per planned item with its size, the temporary space
of an apptainer build, the total and the free space of each filesystem it writes to; free space
below the need x 1.1 is exit 2. ``shutil.disk_usage`` sees the filesystem, not the user's quota,
so ``quota -w`` (and ``lfs quota`` on Lustre) are asked too, best effort with a 10 s timeout.
"""

from __future__ import annotations

import errno
import os
import re
import shutil
import subprocess
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import NamedTuple

import click

from gippa.exitcodes import ExitCode

__all__ = ["MOUNTINFO", "PlanLine", "disk_full", "fs_kind", "is_disk_full", "mount_of",
           "quota_remaining", "review_plan", "user_quota"]  # fmt: skip

MOUNTINFO = Path("/proc/self/mountinfo")
GIB = 1024**3
MARGIN = 1.1
QUOTA_TIMEOUT = 10
DISK_FULL = frozenset({errno.ENOSPC, errno.EDQUOT})
_FULL_TEXT = re.compile(r"no space left on device|disk quota exceeded", re.IGNORECASE)
SHARED = frozenset({"nfs", "nfs4", "lustre", "gpfs", "beegfs", "cephfs", "panfs", "wekafs",
                    "fuse.glusterfs"})  # fmt: skip
LOCAL = frozenset({"ext2", "ext3", "ext4", "xfs", "btrfs", "tmpfs", "overlay", "zfs"})
HOME_WARNING = (
    "warning: the default DB directory is in your home directory, which usually has a small "
    "quota on clusters; pass --db-dir on a project or scratch filesystem"
)
ADVICE = "pass --db-dir (or set GIPPA_DB_DIR) on a larger disk"


def disk_full(where: Path | str, also: str = "") -> str:
    """The message of an item that hit ``EDQUOT``/``ENOSPC`` at ``where``."""
    return (f"disk quota exceeded or disk full at {where}: pass --db-dir{also} (or set "
            "GIPPA_DB_DIR) on a larger disk")  # fmt: skip


def is_disk_full(text: str) -> bool:
    """Whether a runtime's error text says the disk or the quota is full."""
    return _FULL_TEXT.search(text) is not None


def fs_kind(fstype: str) -> str:
    """``shared``, ``local`` or ``unknown`` (item 9)."""
    return "shared" if fstype in SHARED else "local" if fstype in LOCAL else "unknown"


def _existing(path: Path) -> Path:
    p = Path(os.path.realpath(path))
    while not p.exists() and p != p.parent:
        p = p.parent
    return p


def mount_of(path: Path) -> tuple[str, str, str] | None:
    """(mount point, fstype, source) of ``path`` from ``MOUNTINFO``: the longest prefix."""
    try:
        text = MOUNTINFO.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    where, best = str(_existing(path)), None
    for line in text.splitlines():
        head, sep, tail = line.partition(" - ")
        f, t = head.split(), tail.split()
        if not sep or len(f) < 5 or len(t) < 2:
            continue
        mnt = re.sub(r"\\([0-7]{3})", lambda m: chr(int(m.group(1), 8)), f[4])
        inside = mnt == "/" or where == mnt or where.startswith(mnt.rstrip("/") + "/")
        if inside and (best is None or len(mnt) >= len(best[0])):
            best = (mnt, t[0], t[1])
    return best


def quota_remaining(text: str, devices: Iterable[str] | None = None) -> int | None:
    """Bytes left under the tightest limit in ``quota -w`` / ``lfs quota`` output (KiB columns
    ``used soft hard``; limit = soft, else hard); only rows of ``devices`` when given."""
    want = set(devices) if devices is not None else None
    left: list[int] = []
    pending = ""
    for line in text.splitlines():
        tok = line.split()
        if not tok or line.startswith("Disk quotas") or tok[0] == "Filesystem":
            continue
        if len(tok) == 1:
            pending = tok[0]  # lfs wraps a long filesystem name onto its own line
            continue
        numeric = tok[0].rstrip("*").isdigit()
        fs, nums = (pending, tok[:3]) if numeric else (tok[0], tok[1:4])
        pending = ""
        try:
            used, soft, hard = (int(n.rstrip("*")) for n in nums)
        except ValueError:
            continue
        limit = soft or hard
        if limit and (want is None or fs in want):
            left.append(max(0, limit - used) * 1024)
    return min(left) if left else None


def _ask(argv: Sequence[str]) -> str:
    try:
        proc = subprocess.run(argv, capture_output=True, text=True, timeout=QUOTA_TIMEOUT,
                              check=False)  # fmt: skip
    except (OSError, subprocess.SubprocessError, ValueError):  # best effort: never stops setup
        return ""
    return proc.stdout if isinstance(proc.stdout, str) else ""  # over quota exits 1: still read


def user_quota(paths: Iterable[Path], environ: Mapping[str, str]) -> dict[Path, int]:
    """Bytes the user may still write on each path's filesystem, where a quota is known (best
    effort): one ``quota -w`` filtered per place, one ``lfs quota`` per Lustre mount. A path
    whose mount is unknown gets no quota (another filesystem's row must not count)."""
    mounts = {p: m for p in paths if (m := mount_of(p)) is not None}
    if not mounts:
        return {}
    text = _ask(["quota", "-w"])
    user = environ.get("USER") or environ.get("LOGNAME")
    lfs = {m[0]: _ask(["lfs", "quota", "-u", user, m[0]])
           for m in mounts.values() if m[1] == "lustre" and user}  # fmt: skip
    out: dict[Path, int] = {}
    for path, (mnt, _fstype, source) in mounts.items():
        found = [quota_remaining(t, (mnt, source)) for t in (text, lfs.get(mnt, ""))]
        known = [q for q in found if q is not None]
        if known:
            out[path] = min(known)
    return out


class PlanLine(NamedTuple):
    label: str
    gib: float  # 0: the size is not recorded
    where: Path | None  # the directory it is written to; None: nothing to write (present)
    note: str = ""
    temporary: bool = False  # build space, freed afterwards: not in the total


def _say(message: str) -> None:
    print(f"gippa setup: {message}", file=sys.stderr, flush=True)


def _free(where: Path) -> int | None:
    try:
        return shutil.disk_usage(_existing(where)).free
    except OSError as exc:
        _say(f"warning: cannot read the free space of {where} ({exc}); proceeding")
        return None


def _device(where: Path) -> tuple[int | str, str]:
    """The filesystem of ``where``: its device and mount point (a quota is per mount)."""
    mount = mount_of(where)
    try:
        dev: int | str = _existing(where).stat().st_dev
    except OSError:
        dev = str(where)
    return dev, "" if mount is None else mount[0]


def _size(gib: float) -> str:
    return f"{gib:.1f} GiB" if gib else "size not recorded"


def review_plan(
    lines: Sequence[PlanLine], places: Mapping[str, Path], db_dir: Path, *,
    home_warning: bool, yes: bool, environ: Mapping[str, str],
) -> int | None:  # fmt: skip
    """Print the plan; exit code 2 when a filesystem is too small, 0 when the user says no,
    else ``None`` (proceed). ``places`` names the directories written to (role -> path); a path
    shared by two roles counts once, and each path's quota is its own filesystem's."""
    for ln in lines:
        _say(f"plan: {ln.label} {_size(ln.gib)}{f' ({ln.note})' if ln.note else ''}")
    total = sum(ln.gib for ln in lines if ln.where is not None and not ln.temporary)
    _say(f"plan: total {total:.1f} GiB")
    roles: dict[Path, list[str]] = {}
    for role, where in places.items():
        roles.setdefault(where, []).append(role)
    quotas = user_quota(list(roles), environ) or {}  # one quota -w for every place
    groups: dict[tuple[int | str, str], tuple[Path, float, int | None]] = {}
    for where, names in roles.items():
        free = _free(where)
        if free is not None and where in quotas:
            free = min(free, quotas[where])
        if free is not None:
            # a directory not created yet is measured on its nearest existing ancestor (#503)
            on = where if where.exists() else (
                f"{_existing(where)} (nearest existing directory of {where})")  # fmt: skip
            _say(f"free: {free / GIB:.1f} GiB on {on} ({', '.join(names)})")
        need = sum(ln.gib for ln in lines if ln.where == where)
        key = _device(where)
        first, had, low = groups.get(key, (where, 0.0, None))
        known = [f for f in (low, free) if f is not None]
        groups[key] = (first, had + need, min(known) if known else None)
    if home_warning:
        _say(HOME_WARNING)
    for where, need, free in groups.values():
        if free is not None and need and free < need * MARGIN * GIB:
            _say(f"needs {need:.1f} GiB, {free / GIB:.1f} GiB free on {where}: {ADVICE}")
            return int(ExitCode.ENVIRONMENT)
    from gippa.setupdb import uniop_build

    asked = not yes and uniop_build.stdin_is_tty()
    if asked and not click.confirm("gippa setup: Proceed?", default=True, err=True):
        _say("nothing downloaded")
        return int(ExitCode.SUCCESS)
    return None


def setup_plan(
    db_dir: Path, items: Sequence[tuple[str, float, str]], tools: Sequence[str],
    runtime: str | None, cache: Path, environ: Mapping[str, str],
) -> tuple[list[PlanLine], dict[str, Path]]:  # fmt: skip
    """Plan lines and the directories written to (role -> path) for the DB ``items`` (name, GiB,
    path under ``db_dir``) and the images of ``tools``; what already exists is ``present``."""
    from gippa.setupdb.images import image_sizes, pinned_images, simple_name

    lines = [PlanLine(n, g, None, "present") if (db_dir / p).exists() else PlanLine(n, g, db_dir)
             for n, g, p in items]  # fmt: skip
    places = {"DB directory": db_dir}
    pins, sizes, seen = pinned_images(), image_sizes(), set[str]()
    apptainer = runtime == "apptainer"
    for tool in tools:
        ref = pins[tool]
        same = [sizes[t] for t, r in pins.items() if r == ref and t in sizes]
        gib = sizes.get(tool, max(same, default=0.0))
        if apptainer and (cache / simple_name(ref)).is_file():
            line = PlanLine(f"image {tool}", gib, None, "present")
        elif ref in seen:
            line = PlanLine(f"image {tool}", gib, None, "same image as another tool")
        else:
            line = PlanLine(f"image {tool}", gib, cache if apptainer else None)
        seen.add(ref)
        lines.append(line)
    pulled = [ln.gib for ln in lines if apptainer and ln.where == cache]
    if pulled:
        tmp = Path(environ.get("APPTAINER_TMPDIR") or "/tmp")
        lines.append(PlanLine(f"temporary space in {tmp} (APPTAINER_TMPDIR) for the largest "
                              "image build", max(pulled), tmp, temporary=True))  # fmt: skip
        places |= {"image cache": cache, "temporary space": tmp}
    return lines, places


def in_shared_home(path: Path, environ: Mapping[str, str]) -> bool:
    """Whether ``path`` is under ``$HOME`` on a ``shared`` filesystem (the home-quota warning)."""
    home = environ.get("HOME")
    if not home:
        return False
    where, root = str(_existing(path)), os.path.realpath(home)
    mount = mount_of(path)
    inside = where == root or where.startswith(root.rstrip("/") + "/")
    return inside and mount is not None and fs_kind(mount[1]) == "shared"
