"""antiSMASH command line (GIPPA_SPEC.md §7.1.1; M2-18).

antiSMASH is AGPL-3.0: it is never imported, only run as a separate process (R5). Kept in sync
with ``workflow/modules/local/antismash/main.nf``, whose task shell passes
``--databases <db_dir>/antismash`` only when that path is a directory with content and otherwise
leaves antiSMASH on the image's bundled databases. This supersedes M2-18 lead decisions 2 and 4
(``--databases`` always passed and required) per #515.
"""

from __future__ import annotations

from pathlib import Path


def antismash_command(gbff: Path, outdir: Path, threads: int, databases: Path | None) -> list[str]:
    """``antismash --genefinding-tool none --minimal [--databases <db>] --cpus <n> --output-dir
    <outdir> <gbff>``.

    ``databases`` stays a required argument; ``None`` omits the flag (the databases of the image,
    which is what ``gippa setup`` records under a container profile). Supersedes M2-18 lead
    decisions 2 and 4 per #515.
    """
    flag = [] if databases is None else ["--databases", str(databases)]
    return [
        "antismash",
        "--genefinding-tool",
        "none",
        "--minimal",
        *flag,
        "--cpus",
        str(threads),
        "--output-dir",
        str(outdir),
        str(gbff),
    ]
