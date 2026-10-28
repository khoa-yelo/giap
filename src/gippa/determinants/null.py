"""The TU-size-preserving permutation null of the determinants models (GIPPA_SPEC.md §7.8.2
step 4, §7.8.3, §8.1, §8.2; docs/specs/M8-09-determinants-attribution-null.md rules 5-7).

The response moves among TUs with the same number of model rows, whole TU vectors at a time
(gene order inside a TU kept); a TU whose size is unique stays in place (conservative). Each
permutation refits one CV repeat on the repeat-1 folds (decision 16); ``p = (1 + b)/(1 + B)``.
Pure functions, no file I/O.
"""

from __future__ import annotations

import multiprocessing
import zlib
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from threadpoolctl import threadpool_limits  # type: ignore[import-untyped]

from gippa.config.model import DeterminantsConfig
from gippa.determinants.models import MODELS, ModelName, refit_r2
from gippa.stats.permutation import perm_p

__all__ = ["NULL_STREAM", "NullResult", "null_p", "permute_response", "tu_size_blocks"]

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]

#: the stream of ``gippa.rng.for_step``; permutation b draws from its b-th spawned child
NULL_STREAM = "determinants.null"
DEGENERATE = "null_degenerate"
#: pass 1 (ADR-0034): the permutation p is valid only for exchangeable TUs; on every result
UNCALIBRATED = "null_uncalibrated"
UNESTIMATED = "unestimated"
#: with threads > 1 the refits run on spawned worker processes, never threads: the elastic net's
#: path holds the GIL, so threading it -- or overlapping it with a LightGBM thread pool (GIL convoy
#: on a loaded node) -- is slower than one thread (2.5x at 4 threads; ADR-0053, ADR-0063)
#: below this B the spawn of the worker processes (about 2 s each, in parallel) outweighs the
#: refits: LightGBM's refits run on a thread pool and the elastic net's serially (ADR-0053)
_PROCESS_MIN_B = 200
#: permutations per task sent to a worker process (amortises the pickling of the response)
_CHUNK = 8

#: the fixed refit arguments of one worker process (set once by ``_init_worker``)
_WORKER: dict[str, Any] = {}


def _init_worker(X: pd.DataFrame, folds: I64, cfg: DeterminantsConfig, seed: int,
                 tus: list[str]) -> None:  # fmt: skip
    _WORKER.update(X=X, folds=folds, cfg=cfg, seed=seed, tus=tus)
    # one fit per worker process: BLAS/OpenMP pools of one thread (threads workers in all)
    threadpool_limits(1)


def _worker_r2(args: tuple[ModelName, F64]) -> float:
    m, yb = args
    w = _WORKER
    return float(refit_r2(m, w["X"], yb, w["folds"], w["cfg"], w["seed"], 1, tu=w["tus"])[0])


def _process_columns(
    models: Sequence[ModelName], ys: Sequence[F64], X: pd.DataFrame, folds: I64,
    cfg: DeterminantsConfig, seed: int, tus: list[str], threads: int,
) -> list[list[float]]:  # fmt: skip
    """``[[one(m, yb) for yb in ys] for m in models]`` on one pool of ``threads`` spawned worker
    processes, model by model, each column in b order."""
    ctx = multiprocessing.get_context("spawn")  # never fork a process with live BLAS threads
    pool = ProcessPoolExecutor(
        max_workers=min(threads, len(ys)), mp_context=ctx, initializer=_init_worker,
        initargs=(X, folds, cfg, seed, tus),
    )  # fmt: skip
    try:
        return [
            list(pool.map(_worker_r2, [(m, yb) for yb in ys], chunksize=_CHUNK)) for m in models
        ]
    except BaseException:
        pool.shutdown(wait=True, cancel_futures=True)  # a failed refit: drop the queued chunks
        raise
    finally:
        pool.shutdown(wait=True)


@dataclass(frozen=True, slots=True)
class NullResult:
    """The permutation null of one model (rule 6)."""

    model: ModelName
    observed: float  # observed mean CV R² over the repeats
    null_r2: F64  # the B null R², in b order
    b_count: int  # #{null R² >= observed}
    n_perm: int  # B
    p: float  # (1 + b_count)/(1 + B)
    flags: str = ""  # ";"-joined: ``unestimated`` (NaN observed), ``null_degenerate``, then
    # ``null_uncalibrated`` (always, ADR-0034)


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
    (independent of B), so results are reproducible and identical for any ``threads``: each fit
    is single-threaded and independent; with ``threads`` > 1 and B >= ``_PROCESS_MIN_B`` they
    run model by model on one pool of ``threads`` spawned worker processes (one BLAS/OpenMP
    thread each); with a smaller B LightGBM's run on ``threads`` threads and the elastic net's
    serially; with ``threads`` = 1 all serially; collected in b order.
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

    ys = [_permute(yv, blocks, np.random.default_rng(c)) for c in children]

    def one(m: ModelName, yb: F64) -> float:
        return float(refit_r2(m, X, yb, folds, cfg, seed, 1, tu=tus)[0])

    def column(m: ModelName) -> list[float]:
        return [one(m, yb) for yb in ys]

    null = np.empty((B, len(models)), dtype=np.float64)
    # Bit identity across the three paths rests on each refit being deterministic and
    # independent of the BLAS/OpenMP thread count: the parent's serial and threaded refits use
    # the process's default BLAS threads, the workers one; the tests and ADR-0063 pin it.
    if threads > 1 and B >= _PROCESS_MIN_B:
        for j, col in enumerate(_process_columns(models, ys, X, folds, cfg, seed, tus, threads)):
            null[:, j] = col
    else:
        for j, m in enumerate(models):
            if threads > 1 and m != "elastic_net":
                with ThreadPoolExecutor(max_workers=threads) as pool:
                    null[:, j] = list(pool.map(lambda yb, m=m: one(m, yb), ys))
            else:
                null[:, j] = column(m)
    out: dict[str, NullResult] = {}
    for j, m in enumerate(models):
        obs = float(observed_mean_r2[m])
        p = perm_p(obs, null[:, j], alternative="greater")
        n = int(np.count_nonzero(~np.isnan(null[:, j])))
        b = round(p * (1 + n)) - 1 if np.isfinite(p) else 0
        out[m] = NullResult(m, obs, null[:, j].copy(), b, B, p, _flags(obs))
    return out


def _degenerate(model: ModelName, observed: float, B: int) -> NullResult:
    return NullResult(
        model, observed, np.full(B, observed), B, B, 1.0, _flags(observed, DEGENERATE)
    )


def _flags(observed: float, *extra: str) -> str:
    nan = [UNESTIMATED] if np.isnan(observed) else []
    return ";".join([*nan, *extra, UNCALIBRATED])
