"""TransTermHP command line (GIAP_SPEC.md §7.1.1; M2-15).

``expterm.dat`` ships with the bioconda ``transtermhp=2.09`` package (``$CONDA_PREFIX/data``,
M2-15 Q1). Kept in sync with ``workflow/modules/local/transtermhp/main.nf``.
"""

from __future__ import annotations

from pathlib import Path


def transterm_command(genome_fa: Path, coords: Path, expterm: Path) -> list[str]:
    """``transterm -p <expterm.dat> <genome_fa> <genes.coords>`` (output on stdout)."""
    return ["transterm", "-p", str(expterm), str(genome_fa), str(coords)]
