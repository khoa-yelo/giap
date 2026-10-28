"""Per-step random streams (GIPPA_SPEC.md §6.6).

All randomness flows from ``config.seed``; each step gets its own stream keyed by its name, so
adding a step never changes another step's draws.
"""

from __future__ import annotations

import zlib

import numpy as np


def for_step(step_name: str, seed: int) -> np.random.Generator:
    """Generator for ``step_name``: ``default_rng(SeedSequence([seed, crc32(step_name)]))``."""
    return np.random.default_rng(np.random.SeedSequence([seed, zlib.crc32(step_name.encode())]))
