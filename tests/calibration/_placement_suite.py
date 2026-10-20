"""Genome-level placement FDR suite runner (issue M3-21; lead decisions 1-5).

``run_fdr_suite`` simulates replicate i once (genome seed from ``replicate_rng(base_stream, i)``)
and gates, through ``_gate.calibrate_levels`` (ADR-0015 two-block confirmation, sharding):

* every list (decoy type x m in {0, 1, 2}) -- an unexercised list has FDP 0 and passes, so "every
  list the suite exercises" is covered whatever the truth says (decision 4);
* the pooled FDP over the union, each peptide once at its best-cost list (decision 5).

Each check has its own record stream ``<base>.<list>``; the replicate index is recovered from the
check stream's seed, so all checks of replicate i read one simulation. After the gate, from the
replicates actually run: every claimed list must be exercised (mean >= 1 wrong or non-origin
peptide with a target hit at that cost per replicate; else a suite-configuration error), the
unique guard holds in every replicate, and with planted near-paralogs the mutant "adds all loci of
cost <= m_T + 1" is caught by the guard in every replicate (decision 3). FDP / accepted /
exercise tables are appended to ``<log_dir>/placement_m3_21_diagnostics.jsonl`` for the
statistician's sign-off.
"""

from __future__ import annotations

import json
import time
from collections.abc import Callable, Sequence
from pathlib import Path

import _criteria as C
import numpy as np
from _gate import calibrate_levels
from helpers.placement_truth import LISTS, POOLED, Replicate, list_name

SUITE = "placement"
CHECKS: tuple[str, ...] = (*(list_name(x) for x in LISTS), POOLED)


def write_diagnostics(record: dict[str, object]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "placement_m3_21_diagnostics.jsonl").open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


def run_fdr_suite(
    base: str,
    replicate: Callable[[int, Path], Replicate],
    tmp: Path,
    claimed: Sequence[tuple[str, int]],
    *,
    full: bool = False,
    mutant: bool = False,
) -> dict[int, Replicate]:
    """Gate one suite; returns the replicates run (index -> result). ``replicate(seed, dir)``."""
    levels = C.load_tolerances().levels
    R = C.replicates(SUITE, full)
    index = {C.suite_seed(f"{base}.{c}", i): i for c in CHECKS for i in range(2 * R)}
    memo: dict[int, Replicate] = {}

    def rep(i: int) -> Replicate:
        if i not in memo:
            seed = int(C.replicate_rng(base, i).integers(0, 2**31 - 1))
            memo[i] = replicate(seed, tmp / f"r{i}")
        return memo[i]

    t0 = time.perf_counter()
    gates = {}
    for c in CHECKS:
        gates[c] = calibrate_levels(
            SUITE,
            "fdr",
            levels,
            lambda g, c=c: rep(index[int(g.bit_generator.seed_seq.entropy)]).fdp[c],
            R,
            stream=f"{base}.{c}",
        )
    secs = time.perf_counter() - t0
    failed = [line for g in gates.values() if g.failed for line in g.detail()]
    reps = [memo[i] for i in sorted(memo)]
    if reps:
        ex = {c: float(np.mean([r.exercise[c] for r in reps])) for c in CHECKS if c != POOLED}
        write_diagnostics(
            {
                "suite": base,
                "R_run": len(reps),
                "secs": secs,
                "levels": list(levels),
                "mean_fdp": {c: np.mean([r.fdp[c] for r in reps], axis=0).tolist() for c in CHECKS},
                "mean_accepted": {
                    c: np.mean([r.accepted[c] for r in reps], axis=0).tolist() for c in CHECKS
                },
                "mean_exercise": ex,
                "loci_false_share": float(np.mean([r.loci_false_share for r in reps])),
                "mutant_pooled_fdp": np.mean([r.fdp_mutant for r in reps], axis=0).tolist(),
                "mutant_guard_violations": int(sum(len(r.guard_mutant) for r in reps)),
                "gate": [rec for g in gates.values() for rec in g.detail()],
            }
        )
    assert not failed, failed
    if reps:
        bad = [list_name(x) for x in claimed if ex[list_name(x)] < 1.0]
        assert not bad, (
            f"suite-configuration error: claimed lists {bad} not exercised (mean wrong or"
            f" non-origin peptides with a target hit per replicate: {ex})"
        )
        guard = [(i, v) for i, r in memo.items() for v in r.guard]
        assert not guard, f"unique guard (decision 3) violated: {guard[:20]}"
        if mutant:
            missed = [i for i, r in memo.items() if not r.guard_mutant]
            assert not missed, (
                f"the cost <= m_T + 1 mutant survives the guard in replicates {missed}"
            )
    return memo
