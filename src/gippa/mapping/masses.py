"""Residue and modification masses and the generated mass-equivalence list.

GIPPA_SPEC.md Appendix B.5 (masses, the checks the list must pass), §7.2.4 (block exchanges) and
§7.5.3 (the ``mass_equivalent`` flag). The masses live in the bundled TSVs
``gippa/data/residue_masses.tsv`` and ``gippa/data/modifications.tsv``; this module reads them once
at import (package resources, not user input) and enumerates equivalences purely.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from importlib.resources import files
from itertools import combinations, combinations_with_replacement
from types import MappingProxyType
from typing import Literal

__all__ = [
    "BLOCK_DELTA_DA",
    "MOD_DELTA",
    "RESIDUE_MASS",
    "STANDARD_RESIDUES",
    "EquivalenceKind",
    "MassEquivalence",
    "composition_mass",
    "mass_equivalents",
]

EquivalenceKind = Literal["block", "near_block", "substitution", "modification"]

STANDARD_RESIDUES: str = "ACDEFGHKLMNPQRSTVWY"  # the 20 table-11 residues, I merged into L
BLOCK_DELTA_DA: float = 0.001  # composition exchanges below this are "block", else "near_block"
DEFAULT_TOL_DA: float = 0.05


def _parse_tsv(text: str) -> dict[str, float]:
    """Parse a two-column ``key<TAB>value`` table; ``#`` comments and the header row skipped."""
    rows = [ln.split("\t") for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    return {key: float(value) for key, value, *_ in rows[1:]}


def _load(name: str) -> Mapping[str, float]:
    text = files("gippa.data").joinpath(name).read_text(encoding="utf-8")
    return MappingProxyType(_parse_tsv(text))


RESIDUE_MASS: Mapping[str, float] = _load("residue_masses.tsv")
MOD_DELTA: Mapping[str, float] = _load("modifications.tsv")


@dataclass(frozen=True, slots=True)
class MassEquivalence:
    """Two residue strings (or ``"A+Oxidation"`` modified forms) of near-equal mass."""

    left: str
    right: str
    delta_da: float  # |m(left) - m(right)|, rounded to 6 decimals
    kind: EquivalenceKind


DELTA_DECIMALS = 6


def composition_mass(side: str) -> float:
    """Mass of a residue string, or of a ``"<residue>+<modification>"`` form."""
    if "+" in side:
        residue, mod = side.split("+", 1)
        return RESIDUE_MASS[residue] + MOD_DELTA[mod]
    return sum(RESIDUE_MASS[c] for c in side)


def _entry(left: str, right: str, kind: EquivalenceKind, tol_da: float) -> MassEquivalence | None:
    raw = abs(composition_mass(left) - composition_mass(right))
    delta = round(raw, DELTA_DECIMALS)
    # Inclusive tolerance on either the raw or the rounded delta (lead decisions 1 and 4).
    return MassEquivalence(left, right, delta, kind) if min(raw, delta) <= tol_da else None


def _composition_kind(left: str, right: str) -> EquivalenceKind:
    if len(left) == len(right) == 1:
        return "substitution"
    delta = abs(composition_mass(left) - composition_mass(right))
    return "block" if delta < BLOCK_DELTA_DA else "near_block"


def mass_equivalents(tol_da: float = DEFAULT_TOL_DA) -> list[MassEquivalence]:
    """Enumerate near-equal-mass readings with ``|delta| <= tol_da`` (Appendix B.5).

    Over the 20 standard residues (I and L merged; U and O excluded): 1<->1 substitutions,
    1<->2 and 2<->2 composition exchanges (2<->2 pairs sharing a residue are reducible to a
    1<->1 substitution and are skipped), and residue + modification <-> residue. Compositions are
    written with sorted letters, the shorter side on the left. ``delta_da`` is rounded to
    ``DELTA_DECIMALS`` (6) decimals before sorting and output, so the order carries no
    floating-point noise. The result is sorted by (kind, delta_da, left, right).

    The ``kind`` here (``block``/``near_block`` by delta against ``BLOCK_DELTA_DA``) is a label
    of this list only; it is independent of the DP's ``BLOCKS`` set (M3-06).
    """
    comps = list(STANDARD_RESIDUES) + [
        "".join(p) for p in combinations_with_replacement(STANDARD_RESIDUES, 2)
    ]
    found: list[MassEquivalence] = []
    for a, b in combinations(comps, 2):
        if len(a) == len(b) == 2 and set(a) & set(b):
            continue
        left, right = sorted((a, b), key=lambda s: (len(s), s))
        entry = _entry(left, right, _composition_kind(left, right), tol_da)
        if entry is not None:
            found.append(entry)
    for residue in STANDARD_RESIDUES:
        for mod in MOD_DELTA:
            for target in STANDARD_RESIDUES:
                entry = _entry(f"{residue}+{mod}", target, "modification", tol_da)
                if entry is not None:
                    found.append(entry)
    return sorted(found, key=lambda e: (e.kind, e.delta_da, e.left, e.right))
