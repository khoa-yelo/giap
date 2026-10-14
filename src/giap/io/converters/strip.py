"""Modification stripping (GIAP_SPEC.md §5.4).

Removes bracketed or parenthesised modifications (``(UniMod:35)``, ``[+15.995]``, ``[UNIMOD:35]``,
``n[42]``), bare mass deltas (``+15.995``), and then every character that is not ``A``-``Z``
(lower-case residues, ``n``/``c`` terminus markers, ``.``, ``_``).
"""

from __future__ import annotations

import re

__all__ = ["strip_modifications"]

_GROUP = re.compile(r"\([^()]*\)|\[[^\[\]]*\]")
_DELTA = re.compile(r"[+-]\d+(?:\.\d+)?")
_NON_UPPER = re.compile(r"[^A-Z]")


def strip_modifications(seq: str) -> str:
    """Return the plain uppercase residue sequence of a modified peptide string."""
    return _NON_UPPER.sub("", _DELTA.sub("", _GROUP.sub("", seq)))
