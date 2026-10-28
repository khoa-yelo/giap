"""antiSMASH command line (GIPPA_SPEC.md §7.1.1; M2-18).

antiSMASH is AGPL-3.0: it is never imported, only run as a separate process (R5). Kept in sync
with ``workflow/modules/local/antismash/main.nf``, which always passes
``--databases <db_dir>/antismash`` (M2-18 lead decision 2).
"""

from __future__ import annotations

from pathlib import Path


def antismash_command(gbff: Path, outdir: Path, threads: int, databases: Path) -> list[str]:
    """``antismash --genefinding-tool none --minimal --databases <db> --cpus <n> --output-dir
    <outdir> <gbff>``.

    ``databases`` is required (M2-18 lead decision 4): real runs always need the databases
    ``gippa setup`` installed, and ``main.nf`` always passes them.
    """
    return [
        "antismash",
        "--genefinding-tool",
        "none",
        "--minimal",
        "--databases",
        str(databases),
        "--cpus",
        str(threads),
        "--output-dir",
        str(outdir),
        str(gbff),
    ]
