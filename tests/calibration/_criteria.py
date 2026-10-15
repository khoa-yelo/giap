"""Shared calibration harness (issue M4-05; GIAP_SPEC.md §8.4, §11.2, §12).

Every constant comes from the protected ``tests/calibration/tolerances.yaml``; this module only
reads it. Criteria, with k = ``se_multiplier`` and R the number of independent null replicates:

- level:   rate of p <= alpha  <=  alpha + k * sqrt(alpha (1 - alpha) / R)
- uniform: KS p-value against U(0, 1) > ``ks_min_p`` (continuous p only; > 5 % ties refused)
- FDR:     mean FDP <= q + k * sd(FDP, ddof=1) / sqrt(R)
- FWER:    share of replicates with >= 1 false discovery <= alpha + k * sqrt(alpha (1 - alpha) / R)
- pooled level (``m=``, PR #69 lead decision 1): each replicate pools m independent p-values and
  reports its rejection share s_i = #{p <= alpha} / m; sum_i s_i / R <= level_bound(alpha, R m)

Replicate ``rep`` of stream ``suite`` draws from ``default_rng(suite_seed(suite, rep))`` with
``suite_seed = SeedSequence([calibration_base, crc32(suite), rep])`` as one uint64, so any shard
reproduces a replicate alone. Results are appended as JSON lines to
``local/logs/calibration/<suite>.jsonl`` (``GIAP_CALIBRATION_LOG_DIR`` overrides the directory).

Fail closed (PR #54 lead decision 1): ``check_level`` / ``check_uniform`` / ``check_fdr`` refuse
non-finite values or values outside [0, 1], ``check_fwer`` refuses anything but 0/1 indicators and
``check_fdr`` refuses R < 2 (no SE).

``calibrate`` is the gate (ADR-0015, two-block confirmation): replicates 0..R-1 are block 1; if the
criterion fails there, replicates R..2R-1 (same seed formula) are block 2 and the gate fails only if
both blocks fail. Records carry ``block`` (1 | 2), ``status`` (pass | unconfirmed | fail) and the
per-replicate ``values``. Under ``GIAP_CALIBRATION_REPLICATES=start:stop`` (set by ``--shard K/N``)
``calibrate`` runs only that range of block ``GIAP_CALIBRATION_BLOCK`` (default 1, ``--block``) and
logs a ``shard`` record; ``merge_shards`` evaluates the criterion on the full R. A merged check
whose block 1 fails and whose block 2 was not logged is ``pending`` (no verdict; never ``fail``
from block 1 alone). ``--merge`` (PR #54 lead decision 5) then runs block 2 itself: it reruns the
suite's pytest command (``pytest_args``) with ``GIAP_CALIBRATION_CONFIRM`` listing the pending
checks, so ``calibrate`` runs replicates R..2R-1 of exactly those checks (same seeds as unsharded
``calibrate``) and skips every other check, then merges again. Exit codes of ``--merge``: 0 every
check pass or unconfirmed, 1 a confirmed fail (both blocks), 2 merge error (missing or malformed
records, block 2 not obtained).

Run as a script, this module is the nightly job-array entry point (one suite per task, task i = the
i-th key of ``replicates_full``): ``python tests/calibration/_criteria.py --task-id I [--fast]
[--shard K/N] [--block B] [--merge] [--dry-run]``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shlex
import subprocess
import sys
import time
import zlib
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, TypeVar

import numpy as np
import yaml
from numpy.typing import ArrayLike
from scipy import stats as sps

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TOLERANCES_PATH = HERE / "tolerances.yaml"
LOG_DIR_ENV = "GIAP_CALIBRATION_LOG_DIR"
SHARD_ENV = "GIAP_CALIBRATION_REPLICATES"
BLOCK_ENV = "GIAP_CALIBRATION_BLOCK"
CONFIRM_ENV = "GIAP_CALIBRATION_CONFIRM"  # JSON [[stream, criterion, level, R], ...]: block 2 only
CRITERIA = ("level", "uniform", "fdr", "fwer")
STATUSES = ("pass", "unconfirmed", "fail")
SHARD_STATUS = "shard"  # partial record of one shard; only ``merge_shards`` gives a verdict
PENDING_STATUS = "pending"  # merged check whose block 1 fails and whose block 2 is not logged yet
SKIPPED_STATUS = "skipped"  # check not listed in GIAP_CALIBRATION_CONFIRM (not run, not logged)
MERGE_FAIL, MERGE_ERROR = 1, 2  # --merge exit codes: confirmed calibration failure / merge error
MAX_TIE_SHARE = 0.05  # §8.4 discrete-p rule: more repeated values than this means discrete p

T = TypeVar("T")


# --------------------------------------------------------------------------- tolerances


@dataclass(frozen=True)
class Tolerances:
    """Parsed ``tolerances.yaml`` (read-only views)."""

    levels: tuple[float, ...]
    se_multiplier: float
    ks_min_p: float
    calibration_base: int
    replicates_fast: Mapping[str, int]
    replicates_full: Mapping[str, int]


@dataclass(frozen=True)
class CheckResult:
    """Outcome of one criterion: ``ok`` iff ``observed`` lies on the passing side of ``bound``."""

    ok: bool
    observed: float
    bound: float
    R: int
    detail: str


def _counts(table: Mapping[str, Any]) -> Mapping[str, int]:
    return MappingProxyType({str(k): int(v) for k, v in table.items()})


def _parse(path: Path) -> Tolerances:
    with path.open() as fh:
        raw: dict[str, Any] = yaml.safe_load(fh)
    return Tolerances(
        levels=tuple(float(x) for x in raw["levels"]),
        se_multiplier=float(raw["se_multiplier"]),
        ks_min_p=float(raw["ks_min_p"]),
        calibration_base=int(raw["seeds"]["calibration_base"]),
        replicates_fast=_counts(raw["replicates_fast"]),
        replicates_full=_counts(raw["replicates_full"]),
    )


@cache
def _default() -> Tolerances:
    return _parse(TOLERANCES_PATH)


def load_tolerances(path: Path | None = None) -> Tolerances:
    """Tolerances from ``path`` (default: the protected ``tolerances.yaml``)."""
    return _default() if path is None else _parse(Path(path))


def replicates(suite: str, full: bool) -> int:
    """``replicates_full[suite]`` or ``replicates_fast[suite]``; KeyError if unknown."""
    tol = load_tolerances()
    table = tol.replicates_full if full else tol.replicates_fast
    if suite not in table:
        raise KeyError(f"unknown calibration suite {suite!r} (not in tolerances.yaml)")
    return table[suite]


# --------------------------------------------------------------------------- seeds


def suite_seed(suite: str, rep: int) -> int:
    """``SeedSequence([calibration_base, crc32(suite), rep])`` drawn as one uint64."""
    base = load_tolerances().calibration_base
    ss = np.random.SeedSequence([base, zlib.crc32(suite.encode()), rep])
    return int(ss.generate_state(1, np.uint64)[0])


def replicate_rng(suite: str, rep: int) -> np.random.Generator:
    """Independent generator of replicate ``rep`` of stream ``suite``."""
    return np.random.default_rng(suite_seed(suite, rep))


def run_replicates(
    suite: str, fn: Callable[[np.random.Generator], T], R: int, start: int = 0
) -> list[T]:
    """``[fn(replicate_rng(suite, i)) for i in range(start, start + R)]``."""
    return [fn(replicate_rng(suite, i)) for i in range(start, start + R)]


# --------------------------------------------------------------------------- criteria


def _vector(x: ArrayLike, what: str) -> np.ndarray:
    a = np.asarray(x, dtype=float).ravel()
    if a.size == 0:
        raise ValueError(f"{what}: no replicates")
    return a


def _unit_vector(x: ArrayLike, what: str) -> np.ndarray:
    """``_vector`` that fails closed on non-finite values or values outside [0, 1]."""
    a = _vector(x, what)
    bad = ~np.isfinite(a) | (a < 0.0) | (a > 1.0)
    if bad.any():
        i = int(np.flatnonzero(bad)[0])
        raise ValueError(
            f"{what}: {int(bad.sum())} value(s) non-finite or outside [0, 1] (first: index {i}"
            f" = {a[i]!r}); a method bug must not pass the gate"
        )
    return a


def level_bound(alpha: float, R: int) -> float:
    """``alpha + se_multiplier * sqrt(alpha (1 - alpha) / R)``."""
    if R < 1:
        raise ValueError(f"R must be >= 1, got {R}")
    k = load_tolerances().se_multiplier
    return alpha + k * math.sqrt(alpha * (1 - alpha) / R)


def check_level(p: ArrayLike, alpha: float) -> CheckResult:
    """Null rejection rate (share of p <= alpha) against ``level_bound(alpha, R)``."""
    a = _unit_vector(p, "check_level")
    r = int(a.size)
    observed = float(np.count_nonzero(a <= alpha)) / r
    bound = level_bound(alpha, r)
    ok = observed <= bound
    return CheckResult(
        ok,
        observed,
        bound,
        r,
        f"level alpha={alpha}: rate {observed:.6g} vs bound {bound:.6g} (R={r})",
    )


def check_uniform(p: ArrayLike) -> CheckResult:
    """KS test of continuous p against U(0, 1); passes when its p-value exceeds ``ks_min_p``."""
    a = _unit_vector(p, "check_uniform")
    r = int(a.size)
    tie_share = 1.0 - np.unique(a).size / r
    if tie_share > MAX_TIE_SHARE:
        raise ValueError(
            f"check_uniform: {tie_share:.1%} repeated values (> {MAX_TIE_SHARE:.0%} ties); discrete"
            " p-values use the level check only (§8.4)"
        )
    observed = float(sps.kstest(a, "uniform").pvalue)
    bound = load_tolerances().ks_min_p
    ok = observed > bound
    detail = f"KS vs U(0,1): p {observed:.6g} vs > {bound} (R={r})"
    return CheckResult(ok, observed, bound, r, detail)


def check_fdr(fdp: ArrayLike, q: float) -> CheckResult:
    """Mean FDP against ``q + se_multiplier * sd(FDP, ddof=1) / sqrt(R)``; needs R >= 2."""
    a = _unit_vector(fdp, "check_fdr")
    r = int(a.size)
    if r < 2:
        raise ValueError(f"check_fdr: R = {r} < 2 replicates, so SE(FDP) is undefined")
    observed = float(a.mean())
    se = float(a.std(ddof=1)) / math.sqrt(r)
    bound = q + load_tolerances().se_multiplier * se
    ok = observed <= bound
    return CheckResult(
        ok, observed, bound, r, f"FDR q={q}: mean FDP {observed:.6g} vs bound {bound:.6g} (R={r})"
    )


def _valid_m(m: int) -> int:
    if isinstance(m, bool) or int(m) != m or m < 1:
        raise ValueError(f"m must be an integer >= 1, got {m!r}")
    return int(m)


def check_level_pooled(shares: ArrayLike, alpha: float, m: int) -> CheckResult:
    """Pooled level check: replicate i rejects ``shares[i] * m`` of its ``m`` independent p-values.

    The rejection rate over all R m p-values is gated by ``level_bound(alpha, R m)``; ``R`` of the
    result is the number of replicates. Shares outside [0, 1] or not multiples of 1 / m are refused.
    """
    m = _valid_m(m)
    a = _unit_vector(shares, "check_level_pooled")
    counts = a * m
    rounded = np.round(counts)
    if not np.allclose(counts, rounded, rtol=0.0, atol=1e-9 * m):
        i = int(np.flatnonzero(~np.isclose(counts, rounded, rtol=0.0, atol=1e-9 * m))[0])
        raise ValueError(
            f"check_level_pooled: share {a[i]!r} (index {i}) is not a multiple of 1/m = 1/{m}"
        )
    r = int(a.size)
    n = r * m
    observed = float(rounded.sum()) / n
    bound = level_bound(alpha, n)
    ok = observed <= bound
    return CheckResult(
        ok,
        observed,
        bound,
        r,
        f"level alpha={alpha} pooled m={m}: rate {observed:.6g} vs bound {bound:.6g}"
        f" (R={r}, R*m={n})",
    )


def check_fwer(any_false: ArrayLike, alpha: float = 0.05) -> CheckResult:
    """Share of replicates with >= 1 false discovery against ``level_bound(alpha, R)``."""
    a = _unit_vector(any_false, "check_fwer")
    if not np.isin(a, (0.0, 1.0)).all():
        raise ValueError("check_fwer: indicators must be 0/1 (or bool)")
    r = int(a.size)
    observed = float(np.count_nonzero(a)) / r
    bound = level_bound(alpha, r)
    ok = observed <= bound
    return CheckResult(
        ok,
        observed,
        bound,
        r,
        f"FWER alpha={alpha}: rate {observed:.6g} vs bound {bound:.6g} (R={r})",
    )


# --------------------------------------------------------------------------- results log


def _default_log_dir() -> Path:
    env = os.environ.get(LOG_DIR_ENV)
    return Path(env) if env else REPO / "local" / "logs" / "calibration"


def log_dir() -> Path:
    """``$GIAP_CALIBRATION_LOG_DIR`` or ``<repo>/local/logs/calibration``."""
    return _default_log_dir()


def log_result(
    suite: str,
    criterion: str,
    level: float,
    result: CheckResult,
    seed: int,
    duration: float,
    log_dir: Path | None = None,
) -> Path:
    """Append one JSON line for ``result`` to ``<log_dir>/<suite>.jsonl``; return the path."""
    directory = Path(log_dir) if log_dir is not None else _default_log_dir()
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{suite}.jsonl"
    record = {
        "suite": suite,
        "criterion": criterion,
        "level": float(level),
        "observed": float(result.observed),
        "bound": float(result.bound),
        "R": int(result.R),
        "seed": int(seed),
        "duration": float(duration),
        "ok": bool(result.ok),
    }
    with path.open("a") as fh:  # one write per line: safe for concurrent xdist workers
        fh.write(json.dumps(record) + "\n")
    return path


def read_results(path: Path) -> list[dict[str, Any]]:
    """Records of one ``<suite>.jsonl`` file, in order."""
    with Path(path).open() as fh:
        return [json.loads(line) for line in fh if line.strip()]


def _fmt(x: Any, spec: str) -> str:
    return "-" if x is None else format(x, spec)


def summary_markdown(records: Sequence[Mapping[str, Any]]) -> str:
    """Markdown table of results for the PR summary.

    ``ok`` is the block's own result (``pass`` / ``FAIL``); ``status`` is the gate verdict of the
    check (pass | unconfirmed | fail, ADR-0015; ``-`` for records without one).
    """
    lines = [
        "| suite | criterion | level | block | observed | bound | R | seed | duration (s) | ok |"
        " status |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in records:
        ok = r.get("ok")
        lines.append(
            f"| {r['suite']} | {r['criterion']} | {r['level']:g} | {r.get('block', '-')} | "
            f"{_fmt(r.get('observed'), '.4g')} | {_fmt(r.get('bound'), '.4g')} | {r['R']} | "
            f"{r['seed']} | {_fmt(r.get('duration'), '.2f')} | "
            f"{'-' if ok is None else 'pass' if ok else 'FAIL'} | {r.get('status', '-')} |"
        )
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- gate (ADR-0015)


def replicate_range(R: int) -> range:
    """``range(start, stop)`` from ``GIAP_CALIBRATION_REPLICATES=start:stop``, else ``range(R)``.

    Requires ``0 <= start < stop <= R``; anything else raises ValueError (fail closed).
    """
    if R < 1:
        raise ValueError(f"R must be >= 1, got {R}")
    raw = os.environ.get(SHARD_ENV)
    if raw is None:
        return range(R)
    parts = raw.split(":")
    try:
        if len(parts) != 2:
            raise ValueError("expected start:stop")
        start, stop = int(parts[0]), int(parts[1])
    except ValueError as exc:
        raise ValueError(f"{SHARD_ENV}={raw!r}: {exc}") from None
    if not 0 <= start < stop <= R:
        raise ValueError(f"{SHARD_ENV}={raw!r}: need 0 <= start < stop <= R = {R}")
    return range(start, stop)


def replicate_block() -> int:
    """Block of a sharded run: ``GIAP_CALIBRATION_BLOCK`` (1 or 2), default 1."""
    raw = os.environ.get(BLOCK_ENV, "1")
    if raw not in ("1", "2"):
        raise ValueError(f"{BLOCK_ENV}={raw!r}: must be 1 or 2")
    return int(raw)


def check(criterion: str, values: ArrayLike, level: float, m: int | None = None) -> CheckResult:
    """Apply ``criterion`` (level | uniform | fdr | fwer) at ``level`` to the replicate values.

    With ``m`` (level only) the values are per-replicate rejection shares of ``m`` pooled p-values
    (``check_level_pooled``).
    """
    if m is not None:
        if criterion != "level":
            raise ValueError(f"m= applies to the level criterion only, not {criterion!r}")
        return check_level_pooled(values, level, m)
    if criterion == "level":
        return check_level(values, level)
    if criterion == "uniform":
        return check_uniform(values)
    if criterion == "fdr":
        return check_fdr(values, level)
    if criterion == "fwer":
        return check_fwer(values, level)
    raise ValueError(f"unknown criterion {criterion!r}; choose from {', '.join(CRITERIA)}")


def _scalar(x: Any, rep: int) -> float:
    a = np.asarray(x, dtype=float)
    if a.shape != ():
        raise ValueError(f"replicate {rep}: fn must return one scalar, got shape {a.shape}")
    return float(a)


def _record(
    suite: str,
    stream: str,
    criterion: str,
    level: float,
    R: int,
    block: int,
    reps: range,
    values: Sequence[float],
    duration: float,
    result: CheckResult | None,
    status: str,
    m: int | None = None,
) -> dict[str, Any]:
    return {
        "suite": suite,
        "stream": stream,
        "criterion": criterion,
        "level": float(level),
        "observed": None if result is None else float(result.observed),
        "bound": None if result is None else float(result.bound),
        "R": int(R),
        "seed": suite_seed(stream, (block - 1) * R),  # seed of the block's first replicate
        "duration": float(duration),
        "ok": None if result is None else bool(result.ok),
        "block": int(block),
        "status": status,
        "start": int(reps.start),
        "stop": int(reps.stop),
        "values": [float(v) for v in values],
        "m": m,
    }


def _append(records: Sequence[Mapping[str, Any]], suite: str, log_dir: Path | None) -> None:
    directory = Path(log_dir) if log_dir is not None else _default_log_dir()
    directory.mkdir(parents=True, exist_ok=True)
    text = "".join(json.dumps(dict(r)) + "\n" for r in records)
    with (directory / f"{suite}.jsonl").open("a") as fh:  # one write: safe for xdist workers
        fh.write(text)


def _verdict(r1: CheckResult, r2: CheckResult | None) -> str:
    """ADR-0015 status: fail only when both blocks fail; ``pending`` while block 2 is missing."""
    if r1.ok:
        return "pass"
    if r2 is None:
        return PENDING_STATUS  # confirmation block not run: no verdict (never fail from block 1)
    return "fail" if not r2.ok else "unconfirmed"


CheckKey = tuple[str, str, float, int]  # (stream, criterion, level, R)


def confirm_checks() -> frozenset[CheckKey] | None:
    """Checks listed in ``GIAP_CALIBRATION_CONFIRM`` (block-2 run of ``--merge``), else None."""
    raw = os.environ.get(CONFIRM_ENV)
    if raw is None:
        return None
    try:
        items = json.loads(raw)
        return frozenset((str(a), str(b), float(c), int(d)) for a, b, c, d in items)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{CONFIRM_ENV}={raw!r}: {exc}") from None


def calibrate(
    suite: str,
    criterion: str,
    level: float,
    fn: Callable[[np.random.Generator], Any],
    R: int,
    *,
    stream: str | None = None,
    log_dir: Path | None = None,
    m: int | None = None,
) -> dict[str, Any]:
    """Run one calibration check with two-block confirmation (ADR-0015) and log it.

    Replicate ``i`` is ``fn(replicate_rng(stream or suite, i))`` (one scalar). Unsharded: block 1 is
    replicates 0..R-1; if ``criterion`` fails there, block 2 (R..2R-1) is run and the check fails
    only if both blocks fail. Both block records are appended to ``<log_dir>/<suite>.jsonl`` with
    the final ``status``; the final record (last block run) is returned. Sharded
    (``GIAP_CALIBRATION_REPLICATES`` set): only that range of block ``GIAP_CALIBRATION_BLOCK`` is
    run and logged as a ``shard`` record (no verdict; see ``merge_shards``). Under
    ``GIAP_CALIBRATION_CONFIRM`` (set by ``--merge``) a listed check runs all of block 2 as a
    ``shard`` record and any other check is ``skipped`` (nothing run or logged). With ``m``
    (level only) each replicate value is the rejection share of ``m`` pooled p-values and the
    criterion is ``check_level_pooled`` (bound ``level_bound(level, R m)``); records carry ``m``.
    """
    if suite not in array_suites():
        raise KeyError(f"unknown calibration suite {suite!r} (not in tolerances.yaml)")
    if criterion not in CRITERIA:
        raise ValueError(f"unknown criterion {criterion!r}; choose from {', '.join(CRITERIA)}")
    if m is not None:
        m = _valid_m(m)
        if criterion != "level":
            raise ValueError(f"m= applies to the level criterion only, not {criterion!r}")
    name = stream if stream is not None else suite
    confirm = confirm_checks()
    reps = range(R) if confirm is not None else replicate_range(R)

    def run(block: int, local: range) -> tuple[range, list[float], float]:
        offset = (block - 1) * R
        absolute = range(offset + local.start, offset + local.stop)
        t0 = time.perf_counter()
        values = [_scalar(fn(replicate_rng(name, i)), i) for i in absolute]
        return absolute, values, time.perf_counter() - t0

    if confirm is not None:
        if (name, criterion, float(level), int(R)) not in confirm:
            return _record(
                suite, name, criterion, level, R, 1, reps, [], 0.0, None, SKIPPED_STATUS, m
            )
        absolute, values, dt = run(2, reps)
        rec = _record(suite, name, criterion, level, R, 2, absolute, values, dt, None, "shard", m)
        _append([rec], suite, log_dir)
        return rec

    if reps != range(R):
        block = replicate_block()
        absolute, values, dt = run(block, reps)
        rec = _record(
            suite, name, criterion, level, R, block, absolute, values, dt, None, "shard", m
        )
        _append([rec], suite, log_dir)
        return rec

    reps1, v1, dt1 = run(1, reps)
    r1 = check(criterion, v1, level, m)
    if r1.ok:
        rec1 = _record(suite, name, criterion, level, R, 1, reps1, v1, dt1, r1, "pass", m)
        _append([rec1], suite, log_dir)
        return rec1
    reps2, v2, dt2 = run(2, reps)
    r2 = check(criterion, v2, level, m)
    status = _verdict(r1, r2)
    rec1 = _record(suite, name, criterion, level, R, 1, reps1, v1, dt1, r1, status, m)
    rec2 = _record(suite, name, criterion, level, R, 2, reps2, v2, dt2, r2, status, m)
    _append([rec1, rec2], suite, log_dir)
    return rec2


def _collect(
    group: Sequence[Mapping[str, Any]], block: int, R: int, key: tuple[Any, ...]
) -> tuple[list[float], float] | None:
    """Values of replicates ``(block-1)R .. block R - 1`` from the records of one check."""
    offset = (block - 1) * R
    got: dict[int, float] = {}
    duration = 0.0
    seen: set[tuple[int, int]] = set()
    for r in group:
        if int(r.get("block", 1)) != block or "values" not in r:
            continue
        start, stop, values = int(r["start"]), int(r["stop"]), r["values"]
        if len(values) != stop - start or not (offset <= start < stop <= offset + R):
            raise ValueError(f"{key}: malformed record for replicates {start}:{stop}")
        for i, v in zip(range(start, stop), values, strict=True):
            if i in got and got[i] != float(v):
                raise ValueError(f"{key}: replicate {i} logged twice with different values")
            got[i] = float(v)
        if (start, stop) not in seen:
            seen.add((start, stop))
            duration += float(r.get("duration", 0.0))
    if not got:
        return None
    missing = [i for i in range(offset, offset + R) if i not in got]
    if missing:
        raise ValueError(
            f"{key}: block {block} misses {len(missing)} of {R} replicates (first {missing[0]});"
            " a verdict needs the full R"
        )
    return [got[i] for i in range(offset, offset + R)], duration


def merge_shards(suite: str, log_dir: Path | None = None) -> list[dict[str, Any]]:
    """One verdict record per check from the (shard) records in ``<log_dir>/<suite>.jsonl``.

    A check is (stream, criterion, level, R). Block 1 must cover replicates 0..R-1 (else
    ValueError); if it fails, block 2 (R..2R-1) decides per ADR-0015; if block 2 was not logged
    the check is ``pending`` (block-1 record, no verdict: never ``fail`` from block 1 alone).
    Records written by ``log_result`` (no ``values``) are skipped. The records of one check must
    agree on ``m`` (pooled level, else ValueError).
    """
    directory = Path(log_dir) if log_dir is not None else _default_log_dir()
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for r in read_results(directory / f"{suite}.jsonl"):
        if "values" not in r:
            continue
        key = (r.get("stream", suite), r["criterion"], float(r["level"]), int(r["R"]))
        groups.setdefault(key, []).append(r)
    out: list[dict[str, Any]] = []
    for key, group in groups.items():
        stream, criterion, level, R = key
        b1 = _collect(group, 1, R, key)
        if b1 is None:
            raise ValueError(f"{key}: no block-1 replicates logged")
        ms = {r.get("m") for r in group}
        if len(ms) != 1:
            raise ValueError(f"{key}: records disagree on m ({sorted(map(str, ms))})")
        (m,) = ms
        r1 = check(criterion, b1[0], level, m)
        b2 = None if r1.ok else _collect(group, 2, R, key)
        r2 = None if b2 is None else check(criterion, b2[0], level, m)
        status = _verdict(r1, r2)
        block, values, dt, res = (1, *b1, r1) if b2 is None or r2 is None else (2, *b2, r2)
        reps = range((block - 1) * R, block * R)
        rec = _record(suite, stream, criterion, level, R, block, reps, values, dt, res, status, m)
        out.append(rec)
    return out


# --------------------------------------------------------------------------- job array


def array_suites() -> tuple[str, ...]:
    """Suites of the nightly job array: keys of ``replicates_full`` in file order."""
    return tuple(load_tolerances().replicates_full)


def task_suite(task_id: int) -> str:
    """Suite of array task ``task_id``."""
    suites = array_suites()
    if not 0 <= task_id < len(suites):
        raise ValueError(f"task id {task_id} outside 0..{len(suites) - 1}")
    return suites[task_id]


def shard_replicates(R: int, n: int, k: int) -> range:
    """Block ``k`` of ``n`` contiguous, near-equal blocks of ``range(R)``."""
    if n < 1 or not 0 <= k < n:
        raise ValueError(f"shard {k}/{n} out of range")
    q, extra = divmod(R, n)
    start = k * q + min(k, extra)
    return range(start, start + q + (1 if k < extra else 0))


def pytest_args(suite: str, full: bool) -> list[str]:
    """pytest arguments that run one suite (``calsuite_<suite>`` keyword set by conftest.py)."""
    marker = "calibration_full" if full else "calibration and not calibration_full"
    return [str(REPO / "tests" / "calibration"), "-m", marker, "-k", f"calsuite_{suite}"]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse the job-array options; adds ``suite``, ``full`` and ``replicates`` (a range)."""
    ap = argparse.ArgumentParser(description="Run one calibration suite (nightly job array task).")
    ap.add_argument(
        "--task-id", type=int, default=None, help="array index (else SLURM_ARRAY_TASK_ID)"
    )
    ap.add_argument("--suite", default=None, help="suite name instead of a task id")
    ap.add_argument("--fast", action="store_true", help="use replicates_fast")
    ap.add_argument("--shard", default=None, help="K/N: run block K of N replicate blocks")
    ap.add_argument(
        "--block", type=int, choices=(1, 2), default=1, help="replicate block of a sharded run"
    )
    ap.add_argument(
        "--merge",
        action="store_true",
        help="merge shard records into verdicts; runs block 2 of failing checks only",
    )
    ap.add_argument("--dry-run", action="store_true", help="print the pytest command only")
    ns = ap.parse_args(list(argv) if argv is not None else None)
    if ns.suite is None:
        task = ns.task_id
        if task is None and os.environ.get("SLURM_ARRAY_TASK_ID"):
            task = int(os.environ["SLURM_ARRAY_TASK_ID"])
        if task is None:
            ap.error("give --task-id, --suite or SLURM_ARRAY_TASK_ID")
        try:
            ns.suite = task_suite(task)
        except ValueError as exc:
            ap.error(str(exc))
    elif ns.suite not in array_suites():
        ap.error(f"unknown suite {ns.suite!r}; choose from {', '.join(array_suites())}")
    ns.full = not ns.fast
    r = replicates(ns.suite, full=ns.full)
    if ns.shard is None:
        ns.replicates = range(r)
    else:
        try:
            k_s, n_s = ns.shard.split("/")
            ns.replicates = shard_replicates(r, int(n_s), int(k_s))
        except ValueError as exc:
            ap.error(f"--shard {ns.shard!r}: {exc}")
    return ns


def _check_key(r: Mapping[str, Any]) -> list[Any]:
    return [r["stream"], r["criterion"], float(r["level"]), int(r["R"])]


def _run_block2(suite: str, full: bool, pending: Sequence[Mapping[str, Any]]) -> int:
    """Run block 2 of the ``pending`` checks through the suite's pytest command; its exit code."""
    cmd = [sys.executable, "-m", "pytest", *pytest_args(suite, full)]
    checks = json.dumps([_check_key(r) for r in pending])
    print(f"{CONFIRM_ENV}={shlex.quote(checks)} {shlex.join(cmd)}", file=sys.stderr)
    base = {k: v for k, v in os.environ.items() if k not in (SHARD_ENV, BLOCK_ENV, CONFIRM_ENV)}
    return subprocess.run(cmd, env={**base, CONFIRM_ENV: checks}, check=False).returncode


def _merge_main(suite: str, full: bool) -> int:
    """Merge the shard records of ``suite``, running block 2 of failing checks (decision 5).

    0 if every check is pass or unconfirmed, ``MERGE_FAIL`` (1) if one fails in both blocks,
    ``MERGE_ERROR`` (2) if the records cannot be merged or block 2 could not be obtained.
    """
    try:
        merged = merge_shards(suite)
        pending = [r for r in merged if r["status"] == PENDING_STATUS]
        if pending:
            code = _run_block2(suite, full, pending)
            merged = merge_shards(suite)
            still = [r for r in merged if r["status"] == PENDING_STATUS]
            if still:
                raise ValueError(
                    f"block 2 not logged for {len(still)} failing check(s) (first"
                    f" {_check_key(still[0])}; pytest exit {code}); no verdict"
                )
    except (OSError, ValueError) as exc:
        print(f"merge {suite}: {exc}", file=sys.stderr)
        return MERGE_ERROR
    if not merged:
        print(f"merge {suite}: no calibrate records in {log_dir()}", file=sys.stderr)
        return MERGE_ERROR
    path = log_dir() / f"{suite}.merged.jsonl"
    path.write_text("".join(json.dumps(r) + "\n" for r in merged))
    print(summary_markdown(merged), end="")
    return MERGE_FAIL if any(r["status"] == "fail" for r in merged) else 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run (or with ``--dry-run`` print) the pytest command of one suite; ``--merge`` merges."""
    ns = parse_args(argv)
    if ns.merge:
        return _merge_main(ns.suite, ns.full)
    cmd = [sys.executable, "-m", "pytest", *pytest_args(ns.suite, ns.full)]
    env: dict[str, str] = {}
    if ns.shard is not None:  # unsharded runs leave the variables unset (behaviour unchanged)
        env = {SHARD_ENV: f"{ns.replicates.start}:{ns.replicates.stop}", BLOCK_ENV: str(ns.block)}
    print(" ".join([*(f"{k}={v}" for k, v in env.items()), shlex.join(cmd)]))
    if ns.dry_run:
        return 0
    base = {k: v for k, v in os.environ.items() if k not in (SHARD_ENV, BLOCK_ENV)}
    return subprocess.run(cmd, env={**base, **env}, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
