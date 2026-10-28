"""geNOMAD command line (GIPPA_SPEC.md §7.1.1; M2-18).

geNOMAD (DB non-commercial) is never imported, only run as a separate process from its own
environment. Kept in sync with ``workflow/modules/local/genomad/main.nf``.
"""

from __future__ import annotations

from pathlib import Path


def genomad_command(genome: Path, outdir: Path, db: Path, threads: int) -> list[str]:
    """``genomad end-to-end --cleanup --splits 8 --threads <n> <genome> <outdir> <db>``."""
    return [
        "genomad",
        "end-to-end",
        "--cleanup",
        "--splits",
        "8",
        "--threads",
        str(threads),
        str(genome),
        str(outdir),
        str(db),
    ]
