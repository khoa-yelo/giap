"""OSTIR command line (GIPPA_SPEC.md §7.1.1; M2-16).

OSTIR is GPL-3.0: it is never imported, only run as a separate process from its own environment.
Kept in sync with ``workflow/modules/local/ostir/main.nf``.
"""

from __future__ import annotations

from pathlib import Path


def ostir_command(inp: Path, out: Path, threads: int) -> list[str]:
    """``ostir -i <inp> -o <out> -j <threads>``."""
    return ["ostir", "-i", str(inp), "-o", str(out), "-j", str(threads)]
