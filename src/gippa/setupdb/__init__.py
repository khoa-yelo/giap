"""Database and model setup (`gippa setup`; GIPPA_SPEC.md §7.9.1, docs/specs/M2-08)."""

from gippa.setupdb.items import (
    SOURCES,
    DbItem,
    DbRecord,
    SetupOptions,
    count_fasta_records,
    licence_lines,
    plan,
    verify,
    verify_go,
    verify_kegg_links,
)
from gippa.setupdb.manifest import MANIFEST, read_records, read_setup, write_manifest
from gippa.setupdb.run import (
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
    "count_fasta_records",
    "install",
    "licence_lines",
    "plan",
    "read_records",
    "read_setup",
    "run_setup",
    "verify",
    "verify_go",
    "verify_kegg_links",
    "write_manifest",
]
