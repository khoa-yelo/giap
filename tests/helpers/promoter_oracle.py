"""Independent oracle and API lookup for the M2-13 promoter-window tests.

The spec places the functions in ``giap.annotate.promoters`` and ``giap.annotate.promotech_cmd``
without saying which file holds the Promotech parser; ``api(name)`` looks in both (promoters first).

Window geometry (docs/specs/M2-13, S = first base of the start codon 0-based, E = interval end):
``+`` window = genome ``[S + k - 40, S + k)``, ``end3p = S + k``; ``-`` window = revcomp of genome
``[E - k, E - k + 40)``, ``end3p = E - k + 1``. Circular replicons wrap (end3p reported in
``[1, L]``); on linear contigs windows leaving ``[0, L)`` are dropped.
"""

from __future__ import annotations

import importlib
from typing import Any

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
MODULES = ("giap.annotate.promoters", "giap.annotate.promotech_cmd")


def api(name: str) -> Any:
    for mod in MODULES:
        m = importlib.import_module(mod)
        if hasattr(m, name):
            return getattr(m, name)
    raise AttributeError(f"{name} not found in {MODULES}")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def oracle_windows(
    start: int,
    end: int,
    strand: str,
    seq: str,
    circular: bool,
    k_range: tuple[int, int] = (-200, 20),
) -> list[tuple[int, int, str]]:
    """``(k, end3p, window_seq)`` in k-ascending order, by brute force on the genome string."""
    n = len(seq)
    out = []
    for k in range(k_range[0], k_range[1] + 1):
        if strand == "+":
            lo, hi, p0 = start + k - 40, start + k, start + k - 1  # p0: 3'-most base, 0-based
        else:
            lo, hi, p0 = end - k, end - k + 40, end - k
        if not circular and (lo < 0 or hi > n):
            continue
        w = "".join(seq[i % n] for i in range(lo, hi))
        out.append((k, p0 % n + 1, w if strand == "+" else revcomp(w)))
    return out
