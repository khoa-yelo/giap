"""Modification stripping (GIPPA_SPEC.md §5.4).

Removes ``(...)``, ``[...]``, ``{...}`` and ``<...>`` groups repeatedly until none remain (so
nested names such as ``(Acetyl (Protein N-term))`` go entirely), then flanking residues
(``K.PEPTIDE.R`` keeps ``PEPTIDE``), bare mass deltas (``+15.995``), and finally every character
that is not ``A``-``Z`` (lower-case residues, ``n``/``c`` terminus markers, ``.``, ``_``).
"""

from __future__ import annotations

import re

__all__ = ["strip_modifications"]

# An innermost group: one bracket pair with no bracket of any kind inside.
_GROUP = re.compile(r"\([^()\[\]{}<>]*\)|\[[^()\[\]{}<>]*\]|\{[^()\[\]{}<>]*\}|<[^()\[\]{}<>]*>")
_FLANK = re.compile(r"^[A-Z-]\.(.+)\.[A-Z-]$")
_DELTA = re.compile(r"[+-]\d+(?:\.\d+)?")
_NON_UPPER = re.compile(r"[^A-Z]")


def strip_modifications(seq: str) -> str:
    """Return the plain uppercase residue sequence of a modified peptide string."""
    n = 1
    while n:
        seq, n = _GROUP.subn("", seq)
    flank = _FLANK.match(seq)
    if flank is not None:
        seq = flank.group(1)
    return _NON_UPPER.sub("", _DELTA.sub("", seq))
