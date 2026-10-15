"""TerminatorNet command line (GIAP_SPEC.md §7.1.1, G12; M2-15).

TerminatorNet's argument parser matches flags by substring: pass only the documented flags, always
``-s`` (both strands, forward coordinates) and staged file names without whitespace. Kept in sync
with ``workflow/modules/local/terminatornet/main.nf``.
"""

from __future__ import annotations

import re
from pathlib import Path

_SPACE = re.compile(r"\s")


def terminatornet_command(genome_fa: Path, out_tsv: Path, threads: int, home: Path) -> list[str]:
    """``python <home>/TerminatorNet.py -f <genome_fa> -s -n <threads> -o <out_tsv>``."""
    for p in (genome_fa, out_tsv):
        if _SPACE.search(str(p)):
            raise ValueError(f"TerminatorNet file names may not contain whitespace: {str(p)!r}")
    if threads < 1:
        raise ValueError(f"threads must be >= 1, got {threads}")
    return [
        "python", str(home / "TerminatorNet.py"), "-f", str(genome_fa), "-s",
        "-n", str(threads), "-o", str(out_tsv),
    ]  # fmt: skip
