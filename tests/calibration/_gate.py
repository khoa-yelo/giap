"""One calibration criterion at several levels through ``_criteria.calibrate`` (M4-05 decision 7).

Suites check a criterion at every alpha / q in ``levels`` on the *same* replicates (one simulation
per replicate, one value per level). ``calibrate`` takes one scalar per replicate and one level, so
``calibrate_levels`` calls it once per level with a memoised replicate function: replicate ``i`` is
simulated once (keyed by its seed ``suite_seed(stream, i)``) and every level reads its column. The
two-block confirmation (ADR-0015), sharding and the jsonl log are those of ``calibrate``; nothing
here evaluates a criterion.

``Gate.values`` holds the per-replicate outputs of every replicate this process ran, in run order
(block 1, then block 2 if a level needed it; only the shard's range under ``--shard``; empty for
checks skipped under ``GIAP_CALIBRATION_CONFIRM``). Suites use it for non-criterion side checks
(e.g. "something is accepted"), never for the gate itself.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
from numpy.typing import ArrayLike

#: statuses that are not a calibration failure: verdicts pass / unconfirmed (ADR-0015), and the
#: no-verdict records of a sharded run (``shard``) or of a check not listed for block 2
#: (``skipped``)
NOT_FAILED = frozenset({"pass", "unconfirmed", C.SHARD_STATUS, C.SKIPPED_STATUS})


@dataclass
class Gate:
    """Records of one ``calibrate_levels`` call and the replicate outputs it computed."""

    records: list[dict[str, Any]]
    values: list[np.ndarray] = field(default_factory=list)

    @property
    def failed(self) -> list[dict[str, Any]]:
        return [r for r in self.records if r["status"] not in NOT_FAILED]

    def detail(self) -> list[str]:
        return [
            f"{r['stream']} {r['criterion']} level={r['level']:g} block={r['block']}"
            f" status={r['status']} observed={r['observed']} bound={r['bound']} R={r['R']}"
            for r in self.records
        ]

    def assert_ok(self) -> None:
        """Fail the test iff some level is a confirmed ``fail`` (both blocks, ADR-0015)."""
        assert not self.failed, self.detail()

    def mean(self) -> np.ndarray | None:
        """Mean replicate output per level over the replicates run here (None if none ran)."""
        return np.vstack(self.values).mean(axis=0) if self.values else None


def _column(v: np.ndarray, j: int, n_levels: int) -> float:
    if v.size == 1:
        return float(v[0])
    if v.size != n_levels:
        raise ValueError(f"replicate returned {v.size} values for {n_levels} levels")
    return float(v[j])


def calibrate_levels(
    suite: str,
    criterion: str,
    levels: Sequence[float],
    fn: Callable[[np.random.Generator], ArrayLike],
    R: int,
    *,
    stream: str,
    log_dir: Path | None = None,
) -> Gate:
    """``calibrate(suite, criterion, level, ., R, stream=stream)`` for every level in ``levels``.

    ``fn(rng)`` simulates one replicate and returns one value per level (or one scalar used at
    every level, e.g. a p-value). It runs once per replicate however many levels read it.
    """
    cache: dict[int, np.ndarray] = {}
    order: list[int] = []

    def value(rng: np.random.Generator) -> np.ndarray:
        key = int(rng.bit_generator.seed_seq.entropy)  # = suite_seed(stream, i)
        if key not in cache:
            cache[key] = np.atleast_1d(np.asarray(fn(rng), dtype=float)).ravel()
            order.append(key)
        return cache[key]

    n = len(levels)
    records = [
        C.calibrate(
            suite,
            criterion,
            float(a),
            lambda g, j=j: _column(value(g), j, n),
            R,
            stream=stream,
            log_dir=log_dir,
        )
        for j, a in enumerate(levels)
    ]
    return Gate(records, [cache[k] for k in order])
