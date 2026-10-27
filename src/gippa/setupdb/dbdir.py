"""The DB directory when ``--db-dir`` is not given (docs/specs/M11-15 item 6, slice 15c1).

Order: ``--db-dir``, ``$GIPPA_DB_DIR``, the pointer file ``gippa setup`` writes
(``<config home>/gippa/db_dir``, one line), ``<data dir>/db`` with ``<data dir>`` =
``$XDG_DATA_HOME/gippa``, else ``$HOME/.local/share/gippa``. Empty variables count as unset.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Literal, NamedTuple

__all__ = ["DbDir", "config_home", "pointer_file", "resolve_db_dir", "write_pointer"]

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


def write_pointer(config: Path, db_dir: Path) -> Path:
    """Record ``db_dir`` (absolute) for the readers; returns the pointer file."""
    ptr = pointer_file(config)
    ptr.parent.mkdir(parents=True, exist_ok=True)
    ptr.write_text(f"{db_dir.absolute()}\n", encoding="utf-8")
    return ptr
