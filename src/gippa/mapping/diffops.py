"""Cost-1 difference operations of the DP verifier and their numba tables (GIPPA_SPEC.md §7.2.4).

``BLOCKS`` and ``NEAR_BLOCKS`` are fixed by §7.2.4 (not derived from the generated list); the
tests check them against :func:`gippa.mapping.masses.mass_equivalents`. :func:`block_tables`
expands every composition into all residue orders over the M3-02 letter codes (I -> L).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from itertools import permutations
from typing import Literal

import numpy as np

from gippa.mapping.sixframe import ALPHABET, encode_target

__all__ = [
    "BLOCKS",
    "BLOCK_KIND",
    "DIFF_TYPES",
    "NEAR_BLOCKS",
    "NEAR_BLOCK_KIND",
    "BlockTables",
    "DiffType",
    "block_tables",
]

DiffType = Literal["substitution", "swap", "block", "near_block", "frameshift"]
DIFF_TYPES: tuple[DiffType, ...] = ("substitution", "swap", "block", "near_block", "frameshift")

BLOCKS: tuple[tuple[str, str], ...] = (
    ("N", "GG"),
    ("Q", "AG"),
    ("GT", "AS"),
    ("AV", "GL"),
    ("AN", "GQ"),
)
NEAR_BLOCKS: tuple[tuple[str, str], ...] = (("W", "AD"), ("W", "EG"), ("R", "GV"))

BLOCK_KIND: int = 0
NEAR_BLOCK_KIND: int = 1
_KIND_NAMES: tuple[Literal["block", "near_block"], ...] = ("block", "near_block")
_PAD: int = 255  # unused slot of a length-1 side


@dataclass(frozen=True, slots=True)
class BlockTables:
    """One row per ordered exchange ``left <-> right`` (each residue order listed once).

    ``left``/``right`` are ``(n, 2)`` M3-02 codes padded with 255; ``left_len``/``right_len`` hold
    1 or 2; ``kind`` is ``BLOCK_KIND`` or ``NEAR_BLOCK_KIND``. All arrays are C-contiguous uint8.

    Rows are directional (``left`` -> ``right``): each exchange is stored in one orientation only,
    and the M3-06 DP tries both orientations when matching a block.
    """

    left: np.ndarray
    left_len: np.ndarray
    right: np.ndarray
    right_len: np.ndarray
    kind: np.ndarray

    def pairs(self) -> Iterator[tuple[str, str, str]]:
        """Decoded ``(left, right, kind)`` rows in M3-02 letters."""
        for i in range(self.kind.shape[0]):
            left = "".join(ALPHABET[c] for c in self.left[i, : self.left_len[i]])
            right = "".join(ALPHABET[c] for c in self.right[i, : self.right_len[i]])
            yield left, right, _KIND_NAMES[int(self.kind[i])]


def _orders(comp: str) -> list[str]:
    return sorted({"".join(p) for p in permutations(comp)})


def _row(side: str) -> list[int]:
    codes = [int(c) for c in encode_target(side)]
    return codes + [_PAD] * (2 - len(codes))


def block_tables() -> BlockTables:
    """Expand ``BLOCKS`` and ``NEAR_BLOCKS`` into every residue order as uint8 arrays."""
    rows: list[tuple[str, str, int]] = []
    for comps, kind in ((BLOCKS, BLOCK_KIND), (NEAR_BLOCKS, NEAR_BLOCK_KIND)):
        for a, b in comps:
            rows.extend((x, y, kind) for x in _orders(a) for y in _orders(b))

    def arr(values: list[list[int]] | list[int]) -> np.ndarray:
        return np.ascontiguousarray(np.array(values, dtype=np.uint8))

    return BlockTables(
        left=arr([_row(x) for x, _, _ in rows]),
        left_len=arr([len(x) for x, _, _ in rows]),
        right=arr([_row(y) for _, y, _ in rows]),
        right_len=arr([len(y) for _, y, _ in rows]),
        kind=arr([k for _, _, k in rows]),
    )
