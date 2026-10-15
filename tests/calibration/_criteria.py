"""Shared calibration harness (issue M4-05; GIAP_SPEC.md §8.4, §11.2, §12).

Every constant comes from the protected ``tests/calibration/tolerances.yaml``; this module only
reads it. Criteria, with k = ``se_multiplier`` and R the number of independent null replicates:

- level:   rate of p <= alpha  <=  alpha + k * sqrt(alpha (1 - alpha) / R)
- uniform: KS p-value against U(0, 1) > ``ks_min_p`` (continuous p only; > 5 % ties refused)
- FDR:     mean FDP <= q + k * sd(FDP, ddof=1) / sqrt(R)
- FWER:    share of replicates with >= 1 false discovery <= alpha + k * sqrt(alpha (1 - alpha) / R)

Replicate ``rep`` of stream ``suite`` draws from ``default_rng(suite_seed(suite, rep))`` with
``suite_seed = SeedSequence([calibration_base, crc32(suite), rep])`` as one uint64, so any shard
reproduces a replicate alone. Results are appended as JSON lines to
``local/logs/calibration/<suite>.jsonl`` (``GIAP_CALIBRATION_LOG_DIR`` overrides the directory).

Run as a script, this module is the nightly job-array entry point (one suite per task, task i = the
i-th key of ``replicates_full``): ``python tests/calibration/_criteria.py --task-id I [--fast]
[--shard K/N] [--dry-run]``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shlex
import subprocess
import sys
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


def level_bound(alpha: float, R: int) -> float:
    """``alpha + se_multiplier * sqrt(alpha (1 - alpha) / R)``."""
    if R < 1:
        raise ValueError(f"R must be >= 1, got {R}")
    k = load_tolerances().se_multiplier
    return alpha + k * math.sqrt(alpha * (1 - alpha) / R)


def check_level(p: ArrayLike, alpha: float) -> CheckResult:
    """Null rejection rate (share of p <= alpha) against ``level_bound(alpha, R)``."""
    a = _vector(p, "check_level")
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
    a = _vector(p, "check_uniform")
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
    """Mean FDP against ``q + se_multiplier * sd(FDP, ddof=1) / sqrt(R)``."""
    a = _vector(fdp, "check_fdr")
    r = int(a.size)
    observed = float(a.mean())
    se = float(a.std(ddof=1)) / math.sqrt(r) if r > 1 else 0.0
    bound = q + load_tolerances().se_multiplier * se
    ok = observed <= bound
    return CheckResult(
        ok, observed, bound, r, f"FDR q={q}: mean FDP {observed:.6g} vs bound {bound:.6g} (R={r})"
    )


def check_fwer(any_false: ArrayLike, alpha: float = 0.05) -> CheckResult:
    """Share of replicates with >= 1 false discovery against ``level_bound(alpha, R)``."""
    a = _vector(any_false, "check_fwer")
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


def summary_markdown(records: Sequence[Mapping[str, Any]]) -> str:
    """Markdown table of results for the PR summary (``pass`` / ``FAIL``)."""
    lines = [
        "| suite | criterion | level | observed | bound | R | seed | duration (s) | ok |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in records:
        lines.append(
            f"| {r['suite']} | {r['criterion']} | {r['level']:g} | {r['observed']:.4g} | "
            f"{r['bound']:.4g} | {r['R']} | {r['seed']} | {r['duration']:.2f} | "
            f"{'pass' if r['ok'] else 'FAIL'} |"
        )
    return "\n".join(lines) + "\n"


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


def main(argv: Sequence[str] | None = None) -> int:
    """Run (or with ``--dry-run`` print) the pytest command of one suite."""
    ns = parse_args(argv)
    cmd = [sys.executable, "-m", "pytest", *pytest_args(ns.suite, ns.full)]
    env = {SHARD_ENV: f"{ns.replicates.start}:{ns.replicates.stop}"}
    print(" ".join(f"{k}={v}" for k, v in env.items()), shlex.join(cmd))
    if ns.dry_run:
        return 0
    return subprocess.run(cmd, env={**os.environ, **env}, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
