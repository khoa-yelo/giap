"""The DB directory when ``--db-dir`` is not given (docs/specs/M11-15 item 6, slice 15c1).

Order: ``--db-dir``, ``$GIPPA_DB_DIR``, the pointer file ``gippa setup`` writes
(``<config home>/gippa/db_dir``, one line), ``<data dir>/db`` with ``<data dir>`` =
``$XDG_DATA_HOME/gippa``, else ``$HOME/.local/share/gippa``. Empty variables count as unset.
"""

from __future__ import annotations

import contextlib
import errno
import os
from collections.abc import Mapping
from pathlib import Path
from typing import Literal, NamedTuple

__all__ = [
    "DbDir",
    "config_home",
    "creatable",
    "pointer_file",
    "resolve_db_dir",
    "writable",
    "write_pointer",
]

Source = Literal["flag", "env", "pointer", "default"]


class DbDir(NamedTuple):
    path: Path
    source: Source
    warning: str = ""  # an empty or unreadable pointer file, for the caller to print


def config_home(environ: Mapping[str, str]) -> Path:
    """``$XDG_CONFIG_HOME``, else ``$HOME/.config``."""
    xdg = environ.get("XDG_CONFIG_HOME")
    return Path(xdg) if xdg else Path(environ.get("HOME") or Path.home()) / ".config"


def pointer_file(config: Path) -> Path:
    return config / "gippa" / "db_dir"


def _data_dir(environ: Mapping[str, str]) -> Path:
    xdg = environ.get("XDG_DATA_HOME")
    if xdg:
        return Path(xdg) / "gippa"
    return Path(environ.get("HOME") or Path.home()) / ".local" / "share" / "gippa"


def resolve_db_dir(flag: Path | None, environ: Mapping[str, str], config: Path) -> DbDir:
    """The DB directory and where it came from (the pointer file is the only thing read)."""
    if flag is not None:
        return DbDir(Path(flag), "flag")
    if env := environ.get("GIPPA_DB_DIR"):
        return DbDir(Path(env), "env")
    ptr = pointer_file(config)
    warning = ""
    if ptr.exists() or ptr.is_symlink():
        try:
            line = ptr.read_text(encoding="utf-8").strip()
        except (OSError, UnicodeDecodeError) as exc:
            line, warning = "", f"warning: ignoring unreadable {ptr} ({exc})"
        else:
            warning = "" if line else f"warning: ignoring empty {ptr}"
        if line:
            return DbDir(Path(line), "pointer")
    return DbDir(_data_dir(environ) / "db", "default", warning)


def write_pointer(config: Path, db_dir: Path) -> str | None:
    """Record ``db_dir`` (absolute) for the readers, atomically (a temporary file renamed over
    the pointer); ``None``, or why it could not be written (the caller warns, setup goes on)."""
    ptr = pointer_file(config)
    tmp = ptr.with_name(f".{ptr.name}.{os.getpid()}.tmp")
    try:
        ptr.parent.mkdir(parents=True, exist_ok=True)
        tmp.write_text(f"{db_dir.absolute()}\n", encoding="utf-8")
        os.replace(tmp, ptr)
    except OSError as exc:
        with contextlib.suppress(OSError):
            tmp.unlink()
        return f"cannot record the DB directory in {ptr} ({exc}); pass --db-dir to the readers"
    return None


def writable(directory: Path) -> OSError | None:
    """Why the existing ``directory`` cannot be written into (``EROFS`` on a read-only filesystem,
    else ``EACCES``); ``None`` when it can."""
    if os.access(directory, os.W_OK | os.X_OK):
        return None
    try:
        read_only = bool(os.statvfs(directory).f_flag & os.ST_RDONLY)
    except OSError:
        read_only = False
    code = errno.EROFS if read_only else errno.EACCES
    return OSError(code, os.strerror(code), str(directory))


def creatable(path: Path) -> OSError | None:
    """Why ``path`` cannot be (or is not) a directory ``gippa setup`` can create, without creating
    anything (#503); ``None`` when it is a directory or its nearest existing ancestor is a
    directory the user can write into. Never raises: ``os.path`` answers ``False`` for what an
    unsearchable parent hides, so that parent is the ancestor reported."""
    if os.path.isdir(path):
        return None
    if os.path.lexists(path):
        return FileExistsError(errno.EEXIST, os.strerror(errno.EEXIST), str(path))
    ancestor = path.absolute().parent
    while not os.path.lexists(ancestor) and ancestor != ancestor.parent:
        ancestor = ancestor.parent
    if not os.path.isdir(ancestor):
        return NotADirectoryError(errno.ENOTDIR, os.strerror(errno.ENOTDIR), str(ancestor))
    return writable(ancestor)
