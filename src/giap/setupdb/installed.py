"""Installed ``giap setup`` items that pipeline steps read (GIAP_SPEC.md §7.9.1)."""

from __future__ import annotations

from pathlib import Path

__all__ = ["cog_table"]


def cog_table(db_dir: Path | None) -> Path | None:
    """The COG 2020 ``cog-20.def.tab`` (path from ``db_sources.yaml`` ``cog2020_def``) when
    ``db_dir``'s setup manifest lists ``cog2020_def`` (M4-15 decisions 3, 17, 22), else ``None``."""
    from giap.preflight.environment import read_manifest
    from giap.setupdb.items import SOURCES

    if db_dir is None:
        return None
    path = db_dir / SOURCES["cog2020_def"].path
    return path if "cog2020_def" in (read_manifest(str(db_dir)) or {}) and path.is_file() else None
