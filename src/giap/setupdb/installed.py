"""Installed ``giap setup`` items that pipeline steps read (GIAP_SPEC.md §7.9.1)."""

from __future__ import annotations

from pathlib import Path

__all__ = ["FUNCTION_RESOURCES", "cog_table", "resource"]

#: the setup items the ``function`` step reads (M8-06 item 7), in ``read_resources`` key form
FUNCTION_RESOURCES: tuple[str, ...] = ("go_basic", "goslim_generic", "kegg_ko_pathway",
                                       "kegg_pathway_names", "cog2020_fun")  # fmt: skip


def resource(db_dir: Path | None, name: str) -> Path | None:
    """``db_dir / SOURCES[name].path`` when ``db_dir``'s setup manifest lists ``name`` and the
    file exists (M4-15 decision 5, generalised in M8-06 item 7), else ``None``."""
    from giap.preflight.environment import read_manifest
    from giap.setupdb.items import SOURCES

    if db_dir is None:
        return None
    path = db_dir / SOURCES[name].path
    return path if name in (read_manifest(str(db_dir)) or {}) and path.is_file() else None


def cog_table(db_dir: Path | None) -> Path | None:
    """The COG 2020 ``cog-20.def.tab`` (``db_sources.yaml`` ``cog2020_def``) when installed
    (M4-15 decisions 3, 17, 22), else ``None``."""
    return resource(db_dir, "cog2020_def")
