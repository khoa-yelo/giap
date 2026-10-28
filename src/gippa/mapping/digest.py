"""In-silico tryptic digestion (GIPPA_SPEC.md §11.6 TE1/TE8; ``abundance.enzyme``; M3-22).

``trypsin_p`` cuts after every K/R (also before P, the DIA-NN default); ``trypsin`` does not cut
when the next residue is P. The protein end closes the last product. A product spans at most
``missed`` internal cut sites and has ``min_len..max_len`` residues (inclusive). The protein is
digested as written, with its initiator M (lead decision 1): no Met-excised products. Products
are ``(start, end, peptide)`` in 0-based half-open protein coordinates, sorted by (start, end).
"""

from __future__ import annotations

from gippa.config.model import Enzyme

__all__ = ["cut_sites", "digest"]

_ENZYMES: frozenset[str] = frozenset({"trypsin_p", "trypsin"})


def cut_sites(protein: str, enzyme: Enzyme = "trypsin_p") -> list[int]:
    """Product end positions: after each cleaved K/R, plus ``len(protein)`` (none if empty).

    Raises:
        ValueError: for an unknown ``enzyme``.
    """
    if enzyme not in _ENZYMES:
        raise ValueError(f"unknown enzyme {enzyme!r}; expected one of {sorted(_ENZYMES)}")
    n = len(protein)
    blocked = "P" if enzyme == "trypsin" else ""
    ends = [
        i + 1
        for i in range(n - 1)
        if protein[i] in "KR" and (not blocked or protein[i + 1] != blocked)
    ]
    return [*ends, n] if n else []


def digest(
    protein: str,
    enzyme: Enzyme = "trypsin_p",
    missed: int = 1,
    min_len: int = 7,
    max_len: int = 30,
) -> list[tuple[int, int, str]]:
    """Tryptic products of ``protein`` with at most ``missed`` missed cleavages.

    Raises:
        ValueError: for an unknown ``enzyme`` or a negative ``missed``.
    """
    if missed < 0:
        raise ValueError(f"missed must be >= 0, got {missed}")
    ends = cut_sites(protein, enzyme)
    starts = [0, *ends[:-1]]
    out: list[tuple[int, int, str]] = []
    for i, s in enumerate(starts):
        for e in ends[i : i + missed + 1]:
            n = e - s
            if n > max_len:
                break
            if n >= min_len:
                out.append((s, e, protein[s:e]))
    return out
