"""Decoy translations of the six-frame index (GIAP_SPEC.md §7.2.1, §8.2, §8.5 attack 5, §7.6.4).

- ``reverse`` (database input): each frame's amino-acid string reversed.
- ``shuffle`` (de novo input): ``*`` and ``X`` stay in place; within each maximal stop-free
  segment the non-``X`` residues are cut into consecutive windows of ``window`` residues from the
  segment start (the last may be shorter) and each window is permuted uniformly. Frames are
  processed in index order from one generator, so the stream is deterministic (§6.6); the step
  seeds it with ``giap.rng.for_step("sixframe_index", config.seed)``.

Both keep the residue multiset, frame length and frame identity, so every decoy position maps back
to a mirrored target position (:func:`mirror_aa`). Decoy hits are only ever counted, never placed.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numba import njit

from giap.mapping.sixframe import AMBIG, STOP, Frame, SixFrameIndex

DecoyType = Literal["reverse", "shuffle"]


def reverse_frame(aa: np.ndarray) -> np.ndarray:
    """The frame's amino-acid codes reversed (a new uint8 array)."""
    return np.ascontiguousarray(aa[::-1], dtype=np.uint8)


@njit(cache=True, nogil=True)
def _window_slots(
    aa: np.ndarray, window: int
) -> tuple[np.ndarray, np.ndarray]:  # pragma: no cover - compiled by numba; via shuffle_frame
    """Positions of the shuffled residues and each one's index inside its window."""
    n = aa.shape[0]
    pos = np.empty(n, dtype=np.int64)
    local = np.empty(n, dtype=np.int64)
    m = 0
    k = 0
    for i in range(n):
        c = aa[i]
        if c == STOP:
            k = 0
        elif c != AMBIG:
            pos[m] = i
            local[m] = k
            m += 1
            k += 1
            if k == window:
                k = 0
    return pos[:m], local[:m]


@njit(cache=True, nogil=True)
def _apply_swaps(
    out: np.ndarray, pos: np.ndarray, local: np.ndarray, draws: np.ndarray
) -> None:  # pragma: no cover - compiled by numba; via shuffle_frame
    """Inside-out Fisher-Yates per window: slot ``l`` swaps with slot ``draws[s]`` in ``[0, l]``."""
    for s in range(pos.shape[0]):
        j = s - local[s] + draws[s]
        a = pos[s]
        b = pos[j]
        tmp = out[a]
        out[a] = out[b]
        out[b] = tmp


def shuffle_frame(aa: np.ndarray, window: int, rng: np.random.Generator) -> np.ndarray:
    """Window-shuffle one frame (new array); ``*`` and ``X`` keep their positions."""
    if window < 1:
        raise ValueError(f"shuffle window must be >= 1, got {window}")
    out = np.array(aa, dtype=np.uint8, copy=True)
    pos, local = _window_slots(out, window)
    draws = rng.integers(0, local + 1, dtype=np.int64)
    _apply_swaps(out, pos, local, draws)
    return out


def build_decoys(
    index: SixFrameIndex, kind: DecoyType, window: int, rng: np.random.Generator
) -> tuple[Frame, ...]:
    """Decoy frames in ``index.frames`` order; ``reverse`` draws nothing from ``rng``."""
    if kind == "reverse":
        return tuple(_like(f, reverse_frame(f.aa)) for f in index.frames)
    if kind == "shuffle":
        return tuple(_like(f, shuffle_frame(f.aa, window, rng)) for f in index.frames)
    raise ValueError(f"unknown decoy type {kind!r}")


def _like(frame: Frame, aa: np.ndarray) -> Frame:
    return Frame(frame.contig, frame.strand, frame.offset, aa, frame.ext_nt)


def mirror_aa(kind: DecoyType, n_aa: int, aa_start: int, aa_end: int) -> tuple[int, int]:
    """Target window mirrored by decoy window ``[aa_start, aa_end)`` of a frame of ``n_aa``."""
    if kind == "reverse":
        return n_aa - aa_end, n_aa - aa_start
    if kind == "shuffle":
        return aa_start, aa_end
    raise ValueError(f"unknown decoy type {kind!r}")
