"""UniOP command line (GIPPA_SPEC.md §7.1.1; M2-12).

UniOP has no licence (ADR-0008): it is never imported, only run as a separate process from its
own environment. Kept in sync with ``workflow/modules/local/uniop/main.nf``.
"""

from __future__ import annotations

from pathlib import Path


def uniop_command(gff: Path, outdir: Path, uniop_home: Path) -> list[str]:
    """``python <uniop_home>/src/UniOP -f <gff> --gff_source prodigal_gff -t <outdir>``."""
    return [
        "python",
        str(uniop_home / "src" / "UniOP"),
        "-f",
        str(gff),
        "--gff_source",
        "prodigal_gff",
        "-t",
        str(outdir),
    ]
