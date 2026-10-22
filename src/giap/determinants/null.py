"""The TU-size-preserving permutation null of the determinants models (GIAP_SPEC.md §7.8.2
step 4, §7.8.3, §8.1, §8.2; docs/specs/M8-09-determinants-attribution-null.md rules 5-7).

The response moves among TUs with the same number of model rows, whole TU vectors at a time
(gene order inside a TU kept); a TU whose size is unique stays in place (conservative). Each
permutation refits one CV repeat on the repeat-1 folds (decision 16); ``p = (1 + b)/(1 + B)``.
Pure functions, no file I/O.
"""

from __future__ import annotations

import zlib
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.config.model import DeterminantsConfig
from giap.determinants.models import MODELS, ModelName, refit_r2
from giap.stats.permutation import perm_p

__all__ = ["NULL_STREAM", "NullResult", "null_p", "permute_response", "tu_size_blocks"]

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]

#: the stream of ``giap.rng.for_step``; permutation b draws from its b-th spawned child
NULL_STREAM = "determinants.null"
DEGENERATE = "null_degenerate"


@dataclass(frozen=True, slots=True)
class NullResult:
    """The permutation null of one model (rule 6)."""

    model: ModelName
    observed: float  # observed mean CV R² over the repeats
    null_r2: F64  # the B null R², in b order
    b_count: int  # #{null R² >= observed}
    n_perm: int  # B
    p: float  # (1 + b_count)/(1 + B)
    flags: str = ""  # ";"-joined; ``null_degenerate`` when no TU can move


def tu_size_blocks(tu: Sequence[str]) -> list[I64]:
    """Rule 5: one ``(k, s)`` array of row indices per distinct TU size ``s`` (ascending); its
    rows are the ``k`` TUs of that size sorted by TU ID, each TU's rows in input order. Groups of
    one TU are included (the fixed blocks); every row index appears exactly once."""
    ids = np.asarray(tu, dtype=str)
    if ids.size == 0:
        return []
    _, codes, counts = np.unique(ids, return_inverse=True, return_counts=True)
    rows = np.argsort(codes, kind="stable").astype(np.int64)  # by TU ID, input order within
    starts = np.cumsum(counts) - counts
    out: list[I64] = []
    for s in np.unique(counts).tolist():
        first = starts[counts == s]
        out.append(np.asarray(rows[first[:, None] + np.arange(s)], dtype=np.int64))
    return out


def _movable(blocks: Sequence[I64]) -> bool:
    return any(b.shape[0] > 1 for b in blocks)


def _permute(y: F64, blocks: Sequence[I64], rng: np.random.Generator) -> F64:
    out = y.copy()
    for b in blocks:
        if b.shape[0] > 1:
            out[b] = y[b[rng.permutation(b.shape[0])]]
    return out


def permute_response(y: F64, tu: Sequence[str], rng: np.random.Generator) -> F64:
    """Rule 5: ``y`` with whole TU response vectors permuted uniformly among the TUs of each size
    group; TUs of a unique size stay put. ``y`` itself is not modified."""
    return _permute(np.asarray(y, dtype=np.float64), tu_size_blocks(tu), rng)


def null_p(
    X: pd.DataFrame, y: F64, tu: Sequence[str], folds_rep1: I64,
    observed_mean_r2: Mapping[str, float], cfg: DeterminantsConfig, B: int, seed: int,
    threads: int,
) -> dict[str, NullResult]:  # fmt: skip
    """Rule 6: for every model named in ``observed_mean_r2``, B TU-size-preserving permutations
    of ``y``, each refitted on one CV repeat with the repeat-1 folds (``cross_validate`` with
    ``folds_rep1[None, :]`` and ``tu=tu``); ``p = perm_p(observed, null, "greater")``.

    Permutation b draws from child b of ``SeedSequence([seed, crc32(NULL_STREAM)]).spawn(B)``
    (independent of B), so results are reproducible and identical for any ``threads``: the
    permutations run on ``threads`` workers, each fit single-threaded, collected in b order.
    Every TU size unique: every permutation is the identity, so no refit is made, each null R²
    is the observed value, b = B, p = 1 and the result is flagged ``null_degenerate``."""
    if B < 1:
        raise ValueError(f"B must be >= 1, got {B}")
    unknown = set(observed_mean_r2) - set(MODELS)
    if unknown:
        raise ValueError(f"unknown models {sorted(unknown)}")
    models = [m for m in MODELS if m in observed_mean_r2]
    yv = np.asarray(y, dtype=np.float64)
    tus = [str(t) for t in tu]
    blocks = tu_size_blocks(tus)
    if not _movable(blocks):
        return {m: _degenerate(m, float(observed_mean_r2[m]), B) for m in models}
    folds = np.asarray(folds_rep1, dtype=np.int64)[None, :]
    children = np.random.SeedSequence([seed, zlib.crc32(NULL_STREAM.encode())]).spawn(B)

    def one(child: np.random.SeedSequence) -> list[float]:
        yb = _permute(yv, blocks, np.random.default_rng(child))
        return [float(refit_r2(m, X, yb, folds, cfg, seed, 1, tu=tus)[0]) for m in models]

    if threads > 1:
        with ThreadPoolExecutor(max_workers=threads) as pool:
            null = np.array(list(pool.map(one, children)), dtype=np.float64)
    else:
        null = np.array([one(c) for c in children], dtype=np.float64)
    out: dict[str, NullResult] = {}
    for j, m in enumerate(models):
        obs = float(observed_mean_r2[m])
        p = perm_p(obs, null[:, j], alternative="greater")
        n = int(np.count_nonzero(~np.isnan(null[:, j])))
        b = round(p * (1 + n)) - 1 if np.isfinite(p) else 0
        out[m] = NullResult(m, obs, null[:, j].copy(), b, B, p)
    return out


def _degenerate(model: ModelName, observed: float, B: int) -> NullResult:
    return NullResult(model, observed, np.full(B, observed), B, B, 1.0, DEGENERATE)
