"""Shared calibration harness (issue M4-05; GIAP_SPEC.md §8.4, §11.2, §12).

Every constant comes from the protected ``tests/calibration/tolerances.yaml``; this module only
reads it. Criteria, with k = ``se_multiplier`` and R the number of independent null replicates:

- level:   rate of p <= alpha  <=  alpha + k * sqrt(alpha (1 - alpha) / R)
- uniform: KS p-value against U(0, 1) > ``ks_min_p`` (continuous p only; > 5 % ties refused)
- FDR:     mean FDP <= q + k * sd(FDP, ddof=1) / sqrt(R)
- FWER:    share of replicates with >= 1 false discovery <= alpha + k * sqrt(alpha (1 - alpha) / R)
- pooled level (``m=``, PR #69 lead decision 1): each replicate pools m independent p-values and
  reports its rejection share s_i = #{p <= alpha} / m; sum_i s_i / R <= level_bound(alpha, R m)
- pooled F1 (``calibrate_pooled``, M7-13 rule 2): replicate i reports counts (tp_p, n_pred,
  tp_r, n_truth); F1 of P = sum tp_p / sum n_pred and Rc = sum tp_r / sum n_truth >= bound (point
  estimate, optional count guards and vacuity floors; jackknife SE recorded only)
- FDR from counts (``calibrate_fdr_counts``, M7-13 rule 6): replicate i reports (n_called,
  n_false, ...); FDP_i = n_false / max(n_called, 1), then the FDR criterion; optional vacuity
  floors on block-1 count means

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
per-replicate ``values``. Under ``GIAP_CALIBRATION_REPLICATES=start:stop@K/N`` (set by ``--shard
K/N``) ``calibrate`` runs only ``shard_replicates(R, N, K)`` of the check's own R (any R; an empty
shard logs nothing) of block ``GIAP_CALIBRATION_BLOCK`` (default 1, ``--block``) and logs a
``shard`` record (a bare ``start:stop`` is that range of the given R); ``merge_shards``
evaluates the criterion on the full R. A merged check
whose block 1 fails and whose block 2 was not logged is ``pending`` (no verdict; never ``fail``
from block 1 alone). ``--merge`` (PR #54 lead decision 5) then runs block 2 itself: it reruns the
suite's pytest command (``pytest_args``) with ``GIAP_CALIBRATION_CONFIRM`` listing the pending
checks, so ``calibrate`` runs replicates R..2R-1 of exactly those checks (same seeds as unsharded
``calibrate``) and skips every other check, then merges again. Exit codes of ``--merge``: 0 every
check pass or unconfirmed, 1 a confirmed fail (both blocks), 2 merge error (missing or malformed
records, block 2 not obtained).

Tests in a shard (``SHARD_SKIP``, ``skip_without_verdict``, ``owns_direct_checks``): a test logs its
checks, then skips (no verdict); its direct parts (assertions or records over every replicate) run
only where ``owns_direct_checks()`` (unsharded or shard 0 of block 1). ADR-0057: records carry the
owning ``test`` and its ``xfail`` mark (``set_current_test``, called by conftest.py) and
``--merge`` judges each test as pytest would: a confirmed fail of a test with an xfail mark is an
expected failure (``outcome`` xfail, exit unchanged), a strict-xfail test whose checks all pass is
a strict XPASS (exit 1), a non-strict mark is reported only.

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
POOLED_F1 = "pooled_f1"  # criterion of ``calibrate_pooled`` (count records, not scalars)
F1_COUNTS = ("tp_p", "n_pred", "tp_r", "n_truth")
FDR_COUNTS = "fdr_counts"  # criterion of ``calibrate_fdr_counts`` (count records, not scalars)
COUNT_CRITERIA = (POOLED_F1, FDR_COUNTS)  # count-record checks: merged per (stream, key, level, R)
VACUOUS_STATUS = "vacuous"  # pooled check whose block-1 mean of a count is below its floor
STATUSES = ("pass", "unconfirmed", "fail")
SHARD_STATUS = "shard"  # partial record of one shard; only ``merge_shards`` gives a verdict
PENDING_STATUS = "pending"  # merged check whose block 1 fails and whose block 2 is not logged yet
SKIPPED_STATUS = "skipped"  # check not listed in GIAP_CALIBRATION_CONFIRM (not run, not logged)
#: records without a verdict: a shard's range, a merged check awaiting block 2, an unlisted check
NO_VERDICT_STATUSES = frozenset({SHARD_STATUS, PENDING_STATUS, SKIPPED_STATUS})
SHARD_SKIP = "shard or confirm run: the verdict comes from --merge"
XFAILED = "xfailed"  # outcome record of a test that failed as expected in a shard (ADR-0057)
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
    """Pooled level check: replicate i rejects ``shares[i] * m`` of its ``m`` p-values
    (for the ORA arms, categories that are negatively dependent given N, not independent).

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
    check (pass | unconfirmed | fail, ADR-0015; vacuous, M7-13 rule 7; ``-`` for records without
    one).
    """
    lines = [
        "| suite | criterion | level | block | observed | bound | R | seed | duration (s) | ok |"
        " status |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in records:
        ok = r.get("ok")
        status = r.get("status", "-")
        if r.get("outcome"):  # ADR-0057: the owning test's xfail outcome
            strict = (r.get("xfail") or {}).get("strict")
            status = f"{status} ({r['outcome']}{', strict' if strict else ''})"
        lines.append(
            f"| {r['suite']} | {r['criterion']} | {r['level']:g} | {r.get('block', '-')} | "
            f"{_fmt(r.get('observed'), '.4g')} | {_fmt(r.get('bound'), '.4g')} | {r['R']} | "
            f"{r['seed']} | {_fmt(r.get('duration'), '.2f')} | "
            f"{'-' if ok is None else 'pass' if ok else 'FAIL'} | {status} |"
        )
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- gate (ADR-0015)


def _split_shard(raw: str) -> tuple[str, tuple[int, int] | None]:
    """``start:stop`` and the ``(K, N)`` of an optional ``@K/N`` suffix (ValueError if bad)."""
    body, sep, spec = raw.partition("@")
    if not sep:
        return body, None
    try:
        k_s, n_s = spec.split("/")
        k, n = int(k_s), int(n_s)
    except ValueError:
        raise ValueError(f"{SHARD_ENV}={raw!r}: suffix must be @K/N") from None
    if n < 1 or not 0 <= k < n:
        raise ValueError(f"{SHARD_ENV}={raw!r}: need 0 <= K < N")
    return body, (k, n)


def shard_spec() -> tuple[int, int] | None:
    """``(K, N)`` of ``--shard K/N`` (the ``@K/N`` suffix of ``GIAP_CALIBRATION_REPLICATES``)."""
    raw = os.environ.get(SHARD_ENV)
    return None if raw is None else _split_shard(raw)[1]


def replicate_range(R: int) -> range:
    """Replicates of a check with ``R`` replicates that this process runs.

    ``GIAP_CALIBRATION_REPLICATES=start:stop@K/N`` (set by ``--shard K/N``): ``shard_replicates(R,
    N, K)`` for every R, so a check whose R is not the suite's R is sharded too (empty when
    R < N for some K). ``start:stop`` without the suffix: ``range(start, stop)``, requiring
    ``0 <= start < stop <= R``. Unset: ``range(R)``. Anything else raises ValueError (fail closed).
    """
    if R < 1:
        raise ValueError(f"R must be >= 1, got {R}")
    raw = os.environ.get(SHARD_ENV)
    if raw is None:
        return range(R)
    body, spec = _split_shard(raw)
    if spec is not None:
        if len(body.split(":")) != 2:
            raise ValueError(f"{SHARD_ENV}={raw!r}: expected start:stop@K/N")
        return shard_replicates(R, spec[1], spec[0])
    parts = body.split(":")
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


def owns_direct_checks() -> bool:
    """True where a test's direct (non-``calibrate``) parts run: unsharded, or shard 0 of block 1.

    Shard 0 is ``@0/N`` (or, for a bare ``start:stop``, the range starting at 0). False in every
    other shard, in any block-2 shard and in the ``--merge`` confirm run, so exactly one process
    evaluates a direct assertion or writes an all-replicate record, as the unsharded run does.
    """
    if os.environ.get(CONFIRM_ENV) is not None:
        return False
    raw = os.environ.get(SHARD_ENV)
    if raw is None:
        return True
    if replicate_block() != 1:
        return False
    body, spec = _split_shard(raw)
    if spec is not None:
        return spec[0] == 0
    return body.split(":")[0].strip() == "0"


def skip_without_verdict(records: Mapping[str, Any] | Sequence[Mapping[str, Any]]) -> None:
    """``pytest.skip(SHARD_SKIP)`` when any record has no verdict (``NO_VERDICT_STATUSES``).

    Call it after the test's last harness call (so every shard record reaches ``--merge``) and
    before asserting on a status."""
    recs = [records] if isinstance(records, Mapping) else list(records)
    if any(r["status"] in NO_VERDICT_STATUSES for r in recs):
        import pytest

        pytest.skip(SHARD_SKIP)


#: test that is running (node ID) and its evaluated xfail mark; set by conftest.py around each
#: test call and written into every record (ADR-0057)
_OWNER: dict[str, Any] = {"test": None, "xfail": None}


def set_current_test(test: str | None, xfail: Mapping[str, Any] | None) -> None:
    """Name the owner of the records written from now on (conftest.py ``pytest_runtest_call``).

    ``xfail`` is ``{"strict": bool, "reason": str}`` when the test carries an xfail mark that
    applies to an assertion failure, else None."""
    _OWNER["test"] = test
    _OWNER["xfail"] = None if xfail is None else dict(xfail)


def log_test_outcome(
    suite: str,
    test: str,
    outcome: str,
    xfail: Mapping[str, Any] | None,
    log_dir: Path | None = None,
) -> Path:
    """Append the outcome of ``test`` in a shard to ``<log_dir>/<suite>.outcomes.jsonl``.

    conftest.py writes ``xfailed`` when a test with an xfail mark failed as expected in a shard
    (a direct assertion of shard 0), so ``--merge`` does not call the test an XPASS."""
    directory = Path(log_dir) if log_dir is not None else _default_log_dir()
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{suite}.outcomes.jsonl"
    rec = {"suite": suite, "test": test, "outcome": outcome,
           "xfail": None if xfail is None else dict(xfail)}  # fmt: skip
    with path.open("a") as fh:
        fh.write(json.dumps(rec) + "\n")
    return path


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
        "test": _OWNER["test"],  # owning test and its xfail mark (ADR-0057)
        "xfail": _OWNER["xfail"],
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


#: (stream, criterion, level, R), plus the count-record key for pooled checks
CheckKey = tuple[str, str, float, int] | tuple[str, str, float, int, str]


def _parse_check_key(item: Sequence[Any]) -> CheckKey:
    if len(item) == 4:
        return (str(item[0]), str(item[1]), float(item[2]), int(item[3]))
    if len(item) == 5:
        return (str(item[0]), str(item[1]), float(item[2]), int(item[3]), str(item[4]))
    raise ValueError(f"check key {item!r}: expected 4 or 5 items")


def confirm_checks() -> frozenset[CheckKey] | None:
    """Checks listed in ``GIAP_CALIBRATION_CONFIRM`` (block-2 run of ``--merge``), else None."""
    raw = os.environ.get(CONFIRM_ENV)
    if raw is None:
        return None
    try:
        items = json.loads(raw)
        return frozenset(_parse_check_key(item) for item in items)
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
        if values:  # an empty shard (R < N) logs nothing
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
    group: Sequence[Mapping[str, Any]],
    block: int,
    R: int,
    key: tuple[Any, ...],
    conv: Callable[[Any], Any] = float,
) -> tuple[list[Any], float] | None:
    """Values (through ``conv``) of replicates ``(block-1)R .. block R - 1`` of one check."""
    offset = (block - 1) * R
    got: dict[int, Any] = {}
    duration = 0.0
    seen: set[tuple[int, int]] = set()
    for r in group:
        if int(r.get("block", 1)) != block or "values" not in r:
            continue
        start, stop, values = int(r["start"]), int(r["stop"]), r["values"]
        if len(values) != stop - start or not (offset <= start < stop <= offset + R):
            raise ValueError(f"{key}: malformed record for replicates {start}:{stop}")
        for i, v in zip(range(start, stop), values, strict=True):
            if i in got and got[i] != conv(v):
                raise ValueError(f"{key}: replicate {i} logged twice with different values")
            got[i] = conv(v)
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
    agree on ``m`` (pooled level, else ValueError). ``pooled_f1`` and ``fdr_counts`` records
    (``calibrate_pooled``, ``calibrate_fdr_counts``) are merged per (stream, key, level, R) by
    ``_merge_pooled``; their verdicts come last.
    """
    directory = Path(log_dir) if log_dir is not None else _default_log_dir()
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    pooled: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for r in read_results(directory / f"{suite}.jsonl"):
        if "values" not in r:
            continue
        if r["criterion"] in COUNT_CRITERIA:
            pkey = (r["criterion"], r["stream"], str(r["key"]), float(r["level"]), int(r["R"]))
            pooled.setdefault(pkey, []).append(r)
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
        out.append(_owned(rec, group, key))
    out.extend(_owned(_merge_pooled(suite, pkey, group), group, pkey)
               for pkey, group in pooled.items())  # fmt: skip
    return out


def _owned(
    rec: dict[str, Any], group: Sequence[Mapping[str, Any]], key: tuple[Any, ...]
) -> dict[str, Any]:
    """``rec`` with the owning ``test`` and ``xfail`` of its records (ADR-0057). Records of one
    check that disagree on the xfail state, or an xfail check written by two tests, are a merge
    error (fail closed: the verdict would depend on which record is read)."""
    xfails = {json.dumps(r.get("xfail"), sort_keys=True) for r in group}
    tests = {r.get("test") for r in group}
    if len(xfails) != 1:
        raise ValueError(f"{key}: records disagree on the xfail mark ({sorted(xfails)})")
    xfail = group[0].get("xfail")
    if xfail is not None and len(tests) != 1:
        raise ValueError(f"{key}: xfail check written by several tests ({sorted(map(str, tests))})")
    rec["test"] = next(iter(tests)) if len(tests) == 1 else None
    rec["xfail"] = xfail
    return rec


def xfailed_tests(suite: str, log_dir: Path | None = None) -> frozenset[str]:
    """Tests with an ``xfailed`` outcome record in ``<log_dir>/<suite>.outcomes.jsonl``."""
    directory = Path(log_dir) if log_dir is not None else _default_log_dir()
    path = directory / f"{suite}.outcomes.jsonl"
    if not path.exists():
        return frozenset()
    return frozenset(str(r["test"]) for r in read_results(path) if r.get("outcome") == XFAILED)


def apply_xfail(merged: Sequence[dict[str, Any]], xfailed: frozenset[str]) -> list[str]:
    """Set ``outcome`` on the merged records as pytest would judge their owning test (ADR-0057);
    return the strict XPASS tests.

    Per test with an xfail mark: it failed if one of its checks is a confirmed ``fail`` or it has
    an ``xfailed`` outcome record (a direct assertion in shard 0); then every check of it is
    ``xfail``. If every check is ``pass`` / ``unconfirmed`` and it did not fail, every check is
    ``xpass``, and a strict mark makes the test a strict XPASS. Records of unmarked tests get
    ``outcome`` None; a ``vacuous`` check leaves its test's outcome None."""
    by_test: dict[str, list[dict[str, Any]]] = {}
    for r in merged:
        r["outcome"] = None
        if r.get("xfail") is not None:
            by_test.setdefault(str(r["test"]), []).append(r)
    xpass: list[str] = []
    for test, recs in by_test.items():
        if test in xfailed or any(r["status"] == "fail" for r in recs):
            outcome = "xfail"
        elif all(r["status"] in ("pass", "unconfirmed") for r in recs):
            outcome = "xpass"
            if recs[0]["xfail"].get("strict"):
                xpass.append(test)
        else:
            continue
        for r in recs:
            r["outcome"] = outcome
    return xpass


# --------------------------------------------------------------------------- pooled F1 (M7-13)


def _f1_totals(counts: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    """Sums of ``F1_COUNTS`` over replicates; fails closed on bad counts."""
    if len(counts) == 0:
        raise ValueError("pooled_f1: no replicates")
    tot = dict.fromkeys(F1_COUNTS, 0)
    for i, c in enumerate(counts):
        for k in F1_COUNTS:
            v = c[k]
            if isinstance(v, bool) or not isinstance(v, int | np.integer) or v < 0:
                raise ValueError(f"pooled_f1: replicate {i} {k} = {v!r} is not a count")
            tot[k] += int(v)
        if c["tp_p"] > c["n_pred"] or c["tp_r"] > c["n_truth"]:
            raise ValueError(f"pooled_f1: replicate {i} has tp > n ({dict(c)})")
    return tot


def pooled_f1(counts: Sequence[Mapping[str, Any]]) -> tuple[float, float, float]:
    """(precision, recall, F1) from the summed counts (M7-13 rule 2); F1 = 0 when either
    denominator sum is 0."""
    t = _f1_totals(counts)
    if t["n_pred"] == 0 or t["n_truth"] == 0:
        p = t["tp_p"] / t["n_pred"] if t["n_pred"] else 0.0
        r = t["tp_r"] / t["n_truth"] if t["n_truth"] else 0.0
        return p, r, 0.0
    p, r = t["tp_p"] / t["n_pred"], t["tp_r"] / t["n_truth"]
    return p, r, (2 * p * r / (p + r) if p + r > 0 else 0.0)


def jackknife_se(counts: Sequence[Mapping[str, Any]]) -> float:
    """Leave-one-replicate-out jackknife SE of the pooled F1 (recorded only); NaN when R < 2."""
    n = len(counts)
    _f1_totals(counts)
    if n < 2:
        return math.nan
    loo = np.array([pooled_f1([c for j, c in enumerate(counts) if j != i])[2] for i in range(n)])
    return math.sqrt((n - 1) / n * float(((loo - loo.mean()) ** 2).sum()))


def check_pooled_min(
    counts: Sequence[Mapping[str, Any]], bound: float, *, guard: Mapping[str, int] | None = None
) -> CheckResult:
    """Pooled F1 >= ``bound`` and every ``guard`` count sum >= its minimum (rule 7)."""
    p, r, f1 = pooled_f1(counts)
    sums = {k: sum(int(c[k]) for c in counts) for k in (guard or {})}
    miss = {k: sums[k] for k, m in (guard or {}).items() if sums[k] < m}
    detail = f"pooled F1 {f1:.6g} (P {p:.4g}, Rc {r:.4g}) vs >= {bound} (R={len(counts)})"
    if miss:
        detail += f"; guard miss {miss} vs {dict(guard or {})}"
    return CheckResult(f1 >= bound and not miss, f1, float(bound), len(counts), detail)


def _is_vacuous(values: Sequence[Mapping[str, Any]], vacuous: Mapping[str, float] | None) -> bool:
    return bool(vacuous) and any(
        float(np.mean([float(v[k]) for v in values])) < float(m) for k, m in (vacuous or {}).items()
    )


def _counts_of(v: Mapping[str, Any]) -> dict[str, int]:
    """A count record as plain ints; refuses non-integer or non-finite values (fail closed: a
    ``ValueError``, so ``--merge`` gives ``MERGE_ERROR``, never an ``OverflowError``)."""
    out: dict[str, int] = {}
    for k, x in v.items():
        if (
            isinstance(x, bool)
            or not isinstance(x, int | np.integer | float)
            or not math.isfinite(x)
            or int(x) != x
        ):
            raise ValueError(f"count {k} = {x!r} is not an integer")
        out[str(k)] = int(x)
    return out


def _validate(criterion: str, values: Sequence[Mapping[str, Any]]) -> None:
    """Fail closed on bad count records before anything is logged."""
    if criterion == POOLED_F1:
        _f1_totals(values)
    else:
        for v in values:
            fdp_of(v)


def _count_check(
    criterion: str,
    values: Sequence[Mapping[str, Any]],
    level: float,
    guard: Mapping[str, int] | None,
) -> CheckResult:
    if criterion == POOLED_F1:
        return check_pooled_min(values, level, guard=guard)
    return check_fdr_counts(values, level)


def _pooled_record(
    criterion: str,
    suite: str,
    stream: str,
    key: str,
    bound: float,
    R: int,
    block: int,
    reps: range,
    values: Sequence[Mapping[str, Any]],
    duration: float,
    result: CheckResult | None,
    status: str,
    guard: Mapping[str, int] | None,
    vacuous: Mapping[str, float] | None,
) -> dict[str, Any]:
    rec = _record(suite, stream, criterion, bound, R, block, reps, [], duration, result, status)
    rec["values"] = [_counts_of(v) for v in values]
    rec["key"] = key
    rec["guard"] = None if guard is None else dict(guard)
    rec["vacuous"] = None if vacuous is None else dict(vacuous)
    if criterion == FDR_COUNTS:
        if values and result is None and status == VACUOUS_STATUS:
            rec["observed"] = float(np.mean([fdp_of(v) for v in values]))
        return rec
    rec["precision"] = rec["recall"] = rec["se"] = None
    if values:
        p, r, f1 = pooled_f1(values)
        se = jackknife_se(values)
        rec.update(precision=p, recall=r, se=se if math.isfinite(se) else None)
        if result is None and status == VACUOUS_STATUS:
            rec["observed"] = f1
    return rec


def _pooled_verdict(
    criterion: str,
    v1: Sequence[Mapping[str, Any]],
    v2: Sequence[Mapping[str, Any]] | None,
    bound: float,
    guard: Mapping[str, int] | None,
    vacuous: Mapping[str, float] | None,
) -> tuple[str, CheckResult | None, CheckResult | None]:
    """(status, block-1 result, block-2 result); vacuity is judged on block 1 (rule 6)."""
    if _is_vacuous(v1, vacuous):
        return VACUOUS_STATUS, None, None
    r1 = _count_check(criterion, v1, bound, guard)
    r2 = None if r1.ok or v2 is None else _count_check(criterion, v2, bound, guard)
    return _verdict(r1, r2), r1, r2


def _calibrate_counts(
    criterion: str,
    suite: str,
    stream: str,
    key: str,
    bound: float,
    records: Callable[[int], Mapping[str, Any]],
    R: int,
    guard: Mapping[str, int] | None,
    vacuous: Mapping[str, float] | None,
    log_dir: Path | None,
) -> dict[str, Any]:
    """Shared gate of the count-record criteria (``calibrate_pooled``, ``calibrate_fdr_counts``)."""
    if suite not in array_suites():
        raise KeyError(f"unknown calibration suite {suite!r} (not in tolerances.yaml)")
    confirm = confirm_checks()
    reps = range(R) if confirm is not None else replicate_range(R)

    def run(block: int, local: range) -> tuple[range, list[dict[str, int]], float]:
        offset = (block - 1) * R
        absolute = range(offset + local.start, offset + local.stop)
        t0 = time.perf_counter()
        values = [_counts_of(records(i)[key]) for i in absolute]
        if values:  # an empty shard (R < N) runs and logs nothing
            _validate(criterion, values)  # fail closed before anything is logged
        return absolute, values, time.perf_counter() - t0

    def rec(block: int, rr: range, v: Sequence[Mapping[str, Any]], dt: float,
            res: CheckResult | None, status: str) -> dict[str, Any]:  # fmt: skip
        return _pooled_record(
            criterion, suite, stream, key, bound, R, block, rr, v, dt, res, status, guard, vacuous
        )

    if confirm is not None:
        if (stream, criterion, float(bound), int(R), key) not in confirm:
            return rec(1, reps, [], 0.0, None, SKIPPED_STATUS)
        out = rec(2, *run(2, reps), None, SHARD_STATUS)
        _append([out], suite, log_dir)
        return out
    if reps != range(R):
        block = replicate_block()
        out = rec(block, *run(block, reps), None, SHARD_STATUS)
        if out["values"]:
            _append([out], suite, log_dir)
        return out

    reps1, v1, dt1 = run(1, reps)
    status, r1, _ = _pooled_verdict(criterion, v1, None, bound, guard, vacuous)
    if r1 is None or r1.ok:
        out = rec(1, reps1, v1, dt1, r1, status)
        _append([out], suite, log_dir)
        return out
    reps2, v2, dt2 = run(2, reps)
    status, r1, r2 = _pooled_verdict(criterion, v1, v2, bound, guard, vacuous)
    recs = [rec(1, reps1, v1, dt1, r1, status), rec(2, reps2, v2, dt2, r2, status)]
    _append(recs, suite, log_dir)
    return recs[1]


def calibrate_pooled(
    suite: str,
    stream: str,
    key: str,
    bound: float,
    records: Callable[[int], Mapping[str, Any]],
    R: int,
    *,
    guard: Mapping[str, int] | None = None,
    vacuous: Mapping[str, float] | None = None,
    log_dir: Path | None = None,
) -> dict[str, Any]:
    """Pooled-F1 gate (M7-13 rule 2, ADR-0015) over the count records ``records(i)[key]``.

    ``records`` takes the absolute replicate index (block 2 is R..2R-1). Pass when the pooled F1
    of block 1 is >= ``bound`` and the ``guard`` sums are met; else block 2 decides (``fail``
    only if both blocks fail). ``vacuous`` (``{count: floor}``) makes block 1 ``vacuous`` (no
    block 2, not a pass) when a block mean is below its floor. Sharding and
    ``GIAP_CALIBRATION_CONFIRM`` behave as in ``calibrate``; the confirm key is
    ``_check_key(record)`` (it includes ``key``).
    """
    return _calibrate_counts(
        POOLED_F1, suite, stream, key, bound, records, R, guard, vacuous, log_dir
    )


def _merge_pooled(
    suite: str, pkey: tuple[Any, ...], group: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """Verdict record of one count-record check from its (shard) records (``merge_shards``)."""
    criterion, stream, key, bound, R = pkey
    settings = {(json.dumps(r.get("guard"), sort_keys=True),
                 json.dumps(r.get("vacuous"), sort_keys=True)) for r in group}  # fmt: skip
    if len(settings) != 1:
        raise ValueError(f"{pkey}: records disagree on guard / vacuous ({sorted(settings)})")
    guard, vacuous = group[0].get("guard"), group[0].get("vacuous")
    b1 = _collect(group, 1, R, pkey, conv=_counts_of)
    if b1 is None:
        raise ValueError(f"{pkey}: no block-1 replicates logged")
    status, r1, _ = _pooled_verdict(criterion, b1[0], None, bound, guard, vacuous)
    b2 = None if r1 is None or r1.ok else _collect(group, 2, R, pkey, conv=_counts_of)
    if b2 is not None:
        status, r1, r2 = _pooled_verdict(criterion, b1[0], b2[0], bound, guard, vacuous)
        return _pooled_record(criterion, suite, stream, key, bound, R, 2, range(R, 2 * R), b2[0],
                              b2[1], r2, status, guard, vacuous)  # fmt: skip
    return _pooled_record(criterion, suite, stream, key, bound, R, 1, range(R), b1[0], b1[1], r1,
                          status, guard, vacuous)  # fmt: skip


# --------------------------------------------------------------------------- FDR from counts


def fdp_of(counts: Mapping[str, Any]) -> float:
    """``n_false / max(n_called, 1)`` (M7-13 rule 6: 0 when nothing is called); fails closed on a
    missing key (KeyError), a negative or non-integer count or ``n_false > n_called``."""
    got: dict[str, int] = {}
    for k in ("n_called", "n_false"):
        v = counts[k]
        if isinstance(v, bool) or not isinstance(v, int | np.integer) or v < 0:
            raise ValueError(f"fdr_counts: {k} = {v!r} is not a count")
        got[k] = int(v)
    if got["n_false"] > got["n_called"]:
        raise ValueError(f"fdr_counts: n_false > n_called ({dict(counts)})")
    return got["n_false"] / max(got["n_called"], 1)


def check_fdr_counts(values: Sequence[Mapping[str, Any]], q: float) -> CheckResult:
    """``check_fdr`` on the per-replicate FDPs of the count records."""
    return check_fdr([fdp_of(v) for v in values], q)


def calibrate_fdr_counts(
    suite: str,
    stream: str,
    key: str,
    q: float,
    records: Callable[[int], Mapping[str, Any]],
    R: int,
    *,
    vacuous: Mapping[str, float] | None = None,
    log_dir: Path | None = None,
) -> dict[str, Any]:
    """FDR gate (§8.4, M7-13 rule 6, ADR-0015) over the count records ``records(i)[key]``.

    As ``calibrate_pooled`` with criterion ``fdr_counts``: mean FDP <= q + k SE(FDP) at
    ``level = q``; ``vacuous`` (``{count: floor}``, e.g. ``{"n_targets": 1 / q}``) makes block 1
    ``vacuous`` when a block mean is below its floor (no block 2, not a pass).
    """
    return _calibrate_counts(FDR_COUNTS, suite, stream, key, q, records, R, None, vacuous, log_dir)


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
        ns.shard_kn = None
    else:
        try:
            k_s, n_s = ns.shard.split("/")
            ns.replicates = shard_replicates(r, int(n_s), int(k_s))
            ns.shard_kn = (int(k_s), int(n_s))
        except ValueError as exc:
            ap.error(f"--shard {ns.shard!r}: {exc}")
    return ns


def _check_key(r: Mapping[str, Any]) -> list[Any]:
    """Confirm-list entry of a check; count-record checks (one stream, many keys) add ``key``."""
    out = [r["stream"], r["criterion"], float(r["level"]), int(r["R"])]
    return [*out, str(r["key"])] if r["criterion"] in COUNT_CRITERIA else out


def _run_block2(suite: str, full: bool, pending: Sequence[Mapping[str, Any]]) -> int:
    """Run block 2 of the ``pending`` checks through the suite's pytest command; its exit code."""
    cmd = [sys.executable, "-m", "pytest", *pytest_args(suite, full)]
    checks = json.dumps([_check_key(r) for r in pending])
    print(f"{CONFIRM_ENV}={shlex.quote(checks)} {shlex.join(cmd)}", file=sys.stderr)
    base = {k: v for k, v in os.environ.items() if k not in (SHARD_ENV, BLOCK_ENV, CONFIRM_ENV)}
    return subprocess.run(cmd, env={**base, CONFIRM_ENV: checks}, check=False).returncode


def _merge_main(suite: str, full: bool) -> int:
    """Merge the shard records of ``suite``, running block 2 of failing checks (decision 5).

    0 if every check is pass, unconfirmed or ``vacuous`` (rule 7: vacuous is neither a pass nor
    a fail; the summary names it), ``MERGE_FAIL`` (1) if one fails in both blocks, ``MERGE_ERROR``
    (2) if the records cannot be merged, block 2 could not be obtained, or the confirm run's
    pytest exited non-zero with no confirmed fail (a test outside the harness failed; M7-13 lead
    decisions 7-8). A confirmed fail keeps ``MERGE_FAIL`` whatever the confirm run's exit.
    """
    code = 0
    if not (log_dir() / f"{suite}.jsonl").exists() and _shards_complete(suite):
        # no harness check in the suite (e.g. denovo_precision): its verdicts are the shard exits
        print(f"merge {suite}: no calibrate records; every shard finished (verdicts: shard exits)")
        return 0
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
    try:
        xpass = apply_xfail(merged, xfailed_tests(suite))
    except (OSError, ValueError, KeyError) as exc:
        print(f"merge {suite}: outcome records: {exc}", file=sys.stderr)
        return MERGE_ERROR
    path = log_dir() / f"{suite}.merged.jsonl"
    path.write_text("".join(json.dumps(r) + "\n" for r in merged))
    print(summary_markdown(merged), end="")
    for test in xpass:
        print(f"merge {suite}: XPASS(strict) {test}: every check passed under a strict xfail")
    if any(r["status"] == "fail" and r.get("outcome") != "xfail" for r in merged):
        return MERGE_FAIL
    if xpass:
        return MERGE_FAIL
    if code != 0:
        print(f"merge {suite}: the block-2 confirm run's pytest exited {code}", file=sys.stderr)
        return MERGE_ERROR
    return 0


WORKTREE_ENV = "GIAP_WORKTREE"


def worktree_giap_error() -> str | None:
    """Under ``GIAP_WORKTREE``: None when ``giap`` imports from ``$GIAP_WORKTREE/src``, else the
    reason. Without it (the main checkout) None. Prints ``giap.__file__`` to stderr (Slurm 1885838
    ran a worktree's tests against the main checkout's src)."""
    wt = os.environ.get(WORKTREE_ENV)
    if not wt:
        return None
    import importlib.util

    spec = importlib.util.find_spec("giap")
    origin = None if spec is None or spec.origin is None else Path(spec.origin).resolve()
    print(f"giap.__file__={origin}", file=sys.stderr)
    src = (Path(wt) / "src").resolve()
    if origin is None or not origin.is_relative_to(src):
        return f"giap.__file__={origin} is not under {WORKTREE_ENV}/src = {src}"
    return None


def main(argv: Sequence[str] | None = None) -> int:
    """Run (or with ``--dry-run`` print) the pytest command of one suite; ``--merge`` merges."""
    ns = parse_args(argv)
    bad = worktree_giap_error()
    if bad is not None:
        print(f"main: {bad}; set PYTHONPATH to the worktree (calibration_full.sbatch does)",
              file=sys.stderr)  # fmt: skip
        return MERGE_ERROR
    if ns.merge:
        return _merge_main(ns.suite, ns.full)
    cmd = [sys.executable, "-m", "pytest", *pytest_args(ns.suite, ns.full)]
    env: dict[str, str] = {}
    if ns.shard is not None:  # unsharded runs leave the variables unset (behaviour unchanged)
        k, n = ns.shard_kn  # @K/N: every check's R is sharded K of N (replicate_range)
        env = {SHARD_ENV: f"{ns.replicates.start}:{ns.replicates.stop}@{k}/{n}",
               BLOCK_ENV: str(ns.block)}  # fmt: skip
    print(" ".join([*(f"{k}={v}" for k, v in env.items()), shlex.join(cmd)]))
    if ns.dry_run:
        return 0
    # a leaked GIAP_CALIBRATION_CONFIRM would make every gate test skip (statistician note 4)
    base = {k: v for k, v in os.environ.items() if k not in (SHARD_ENV, BLOCK_ENV, CONFIRM_ENV)}
    if CONFIRM_ENV in os.environ:
        print(f"main: ignoring {CONFIRM_ENV} (set only by --merge's block-2 run)", file=sys.stderr)
    rc = subprocess.run(cmd, env={**base, **env}, check=False).returncode
    if ns.shard is not None:
        _log_shard_done(ns.suite, ns.shard_kn, ns.block, rc)
    return rc


def _log_shard_done(suite: str, kn: tuple[int, int], block: int, rc: int) -> None:
    """Append ``{suite, shard, n, block, pytest_exit}`` to ``<log_dir>/<suite>.shards.jsonl``."""
    d = log_dir()
    d.mkdir(parents=True, exist_ok=True)
    rec = {"suite": suite, "shard": kn[0], "n": kn[1], "block": block, "pytest_exit": rc}
    with (d / f"{suite}.shards.jsonl").open("a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _shards_complete(suite: str) -> bool:
    """True when ``<suite>.shards.jsonl`` holds every block-1 shard 0..N-1 of one N."""
    path = log_dir() / f"{suite}.shards.jsonl"
    if not path.exists():
        return False
    done = [r for r in read_results(path) if int(r["block"]) == 1]
    ns = {int(r["n"]) for r in done}
    return len(ns) == 1 and {int(r["shard"]) for r in done} == set(range(ns.pop()))


if __name__ == "__main__":
    raise SystemExit(main())
