"""Database and model setup (`giap setup`; GIAP_SPEC.md §7.9.1, docs/specs/M2-08)."""

from giap.setupdb.items import SOURCES, DbItem, DbRecord, SetupOptions, plan, verify
from giap.setupdb.manifest import MANIFEST, read_records, write_manifest
from giap.setupdb.run import (
    CommandRunner,
    Fetcher,
    HttpsFetcher,
    SourceError,
    ToolRunner,
    ToolUnavailable,
    install,
    run_setup,
)

__all__ = [
    "MANIFEST",
    "SOURCES",
    "CommandRunner",
    "DbItem",
    "DbRecord",
    "Fetcher",
    "HttpsFetcher",
    "SetupOptions",
    "SourceError",
    "ToolRunner",
    "ToolUnavailable",
    "install",
    "plan",
    "read_records",
    "run_setup",
    "verify",
    "write_manifest",
]
