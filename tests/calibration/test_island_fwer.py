"""M5-07: island FWER calibration on the giapsim AR(1) backgrounds, the M5 gate (GIAP_SPEC.md §13
row M5, §7.3.6, §8.4 "Islands | AR(1) backgrounds without islands | FWER | 50 | 500"; spec
docs/specs/M5-07-island-calibration-power.md; ADR-0021 with Amendment 1; ADR-0015).

Replicate i (spec "Public interface"): ``scenarios.ar1_background(n_genes, rho)`` drawn with
``replicate_rng(stream, i)``; the gene table is the simulator truth (``detected_true``,
``p_detect_true``, ``abundance_log2_true``, TUs from the layout; ``islands_data.scenario_contigs``);
the real M5-06 step runs with ``stats.bootstrap = 1000``, ``stats.permutations = 1000`` and the
default ``context.islands``. No island is planted, so the indicator is "at least one island
reported" and the criterion is ``fwer`` through ``_gate.calibrate_levels`` (one ``calibrate`` per
tolerances level, two-block confirmation, the suite log ``islands.jsonl``).

- Fast (per PR): 1,000-gene chromosome, rho 0.5, R = ``replicates_fast.islands`` (Example 1:
  bound 0.05 + 2 sqrt(0.05 * 0.95 / 50) = 0.1116).
- Full (the M5 gate): 4,500-gene chromosome plus one 60-gene plasmid, rho in {0, 0.5, 0.9},
  R = ``replicates_full.islands`` per rho, streams ``islands.rho0`` / ``islands.rho05`` /
  ``islands.rho09`` (Example 2: bound 0.0695).
- Recorded, not gated (lead decision 28, ADR-0021 Amendment 1): rho 0.9 at 1,000 and 1,800 genes,
  logged to ``islands_recorded.jsonl``.
- The per-rho distribution of l-hat and l' (methods note) is logged to ``islands_blocklen.jsonl``.
- Teeth: the iid bootstrap (``context.islands.block_len = 1``) on the fast cell must be a
  confirmed ``fail`` of the same criterion.
- Example 6: the same seeds twice give identical indicator vectors.
"""

from __future__ import annotations

import dataclasses
import json
import sys
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import islands_data as idat

SUITE = "islands"
B, PERM = 1000, 1000  # spec: stats.bootstrap = 1000, stats.permutations = 1000
FAST_GENES, FAST_RHO = 1000, 0.5
FULL_GENES, PLASMID_GENES = 4500, 60
FULL_RHOS = (0.0, 0.5, 0.9)
STREAM = {0.0: "islands.rho0", 0.5: "islands.rho05", 0.9: "islands.rho09"}
FAST_STREAM = "islands.fast.rho05"

pytestmark = pytest.mark.calsuite(SUITE)


def gate_pair(n_genes: int, rho: float, plasmid: int | None = None) -> Any:
    """``ar1_background(n_genes, rho)``, with one circular plasmid of ``plasmid`` genes."""
    from giapsim import scenarios as sc

    lcfg, ecfg = sc.ar1_background(n_genes, rho)
    if plasmid is not None:
        lcfg = dataclasses.replace(lcfg, plasmids=(plasmid,))
    return lcfg, ecfg


def replicate(
    rng: np.random.Generator,
    pair: Any,
    levels: Sequence[float],
    *,
    block_len: int | str = "auto",
    sink: list[dict[str, Any]] | None = None,
) -> np.ndarray:
    """One null replicate: the indicator "at least one island reported" at every level.

    The step runs once at ``alpha = max(levels)``; an island is reported at level a iff its
    window p <= a, so the indicator at a is ``min p <= a`` (checked against the table at the top
    level). ``sink`` collects l-hat, l' and the scanned flags per replicon.
    """
    contigs = idat.scenario_contigs(rng, pair)
    cfg = idat.config(bootstrap=B, permutations=PERM, alpha=max(levels), block_len=block_len)
    genes, tus, replicons = idat.frames(contigs)
    out, scan = idat.step(genes, idat.build_series(genes, tus, replicons, cfg), cfg, rng)
    p = idat.min_p(scan)
    assert (len(out) > 0) == (p <= max(levels)), (len(out), p)
    assert all(ell <= 20 for ell in scan.block_len), scan.block_len  # l' <= 20 (ADR-0021)
    if sink is not None:
        sink.append(
            {
                "contigs": [s.contig for s in scan.series],
                "n_tus": [len(s.tu_ids) for s in scan.series],
                "scanned": list(scan.scanned),
                "l_hat": list(scan.l_hat),
                "block_len": list(scan.block_len),
                "min_p": p,
            }
        )
    return np.array([float(p <= a) for a in levels])


def _log_block_lengths(stream: str, rho: float, n_genes: int, sink: list[dict[str, Any]]) -> None:
    """Methods note: the observed l-hat / l' distribution per rho (chromosome = replicon 0)."""
    if not sink:
        return
    l_hat = np.array([s["l_hat"][0] for s in sink], dtype=float)
    ell = np.array([s["block_len"][0] for s in sink], dtype=float)
    vals, n = np.unique(ell, return_counts=True)
    counts = {str(int(k)): int(v) for k, v in zip(vals, n, strict=True)}
    path = C.log_dir() / "islands_blocklen.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    rec = {
        "suite": SUITE,
        "stream": stream,
        "rho": rho,
        "n_genes": n_genes,
        "R": len(sink),
        "chr_tus_median": float(np.median([s["n_tus"][0] for s in sink])),
        "l_hat_quantiles": np.quantile(l_hat, [0.0, 0.1, 0.5, 0.9, 1.0]).tolist(),
        "block_len_quantiles": np.quantile(ell, [0.0, 0.1, 0.5, 0.9, 1.0]).tolist(),
        "block_len_counts": counts,
        "scanned_share": np.mean([s["scanned"] for s in sink], axis=0).tolist(),
    }  # fmt: skip
    with path.open("a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"block lengths {stream}: {rec}")


# --------------------------------------------------------------------------- the gate


@pytest.mark.calibration
def test_fast_fwer_rho05_1000_genes() -> None:
    """Example 1: R = 50, rho 0.5, 1,000 genes; at most 5 of 50 replicates with an island."""
    levels = C.load_tolerances().levels
    pair = gate_pair(FAST_GENES, FAST_RHO)
    sink: list[dict[str, Any]] = []
    gate = calibrate_levels(
        SUITE, "fwer", levels, lambda rng: replicate(rng, pair, levels, sink=sink),
        C.replicates(SUITE, full=False), stream=FAST_STREAM,
    )  # fmt: skip
    _log_block_lengths(FAST_STREAM, FAST_RHO, FAST_GENES, sink)
    print(f"fast gate: {gate.detail()} (block 1 {gate.block_seconds} s)")
    gate.assert_ok()


@pytest.mark.calibration_full
@pytest.mark.parametrize("rho", FULL_RHOS)
def test_full_fwer_gate(rho: float) -> None:
    """Example 2: 4,500-gene chromosome + 60-gene plasmid, R = 500 per rho, rate <= 0.0695."""
    levels = C.load_tolerances().levels
    pair = gate_pair(FULL_GENES, rho, plasmid=PLASMID_GENES)
    sink: list[dict[str, Any]] = []
    gate = calibrate_levels(
        SUITE, "fwer", levels, lambda rng: replicate(rng, pair, levels, sink=sink),
        C.replicates(SUITE, full=True), stream=STREAM[rho],
    )  # fmt: skip
    _log_block_lengths(STREAM[rho], rho, FULL_GENES, sink)
    print(f"full gate rho {rho}: {gate.detail()}")
    if sink:  # the layout is the spec's: one chromosome and one plasmid in one family
        assert all(s["contigs"] == ["chr1", "p1"] for s in sink)
    gate.assert_ok()


# --------------------------------------------------------- recorded cells (not gated)


@pytest.mark.calibration_full
@pytest.mark.parametrize("n_genes", [1000, 1800])
def test_recorded_rho09_smaller_chromosomes(n_genes: int) -> None:
    """ADR-0021 Amendment 1 / lead decision 28: rho 0.9 at 1,000 and 1,800 genes is recorded,
    not gated (outside the claimed dependence regime at these sizes). Logged, never asserted."""
    stream = f"islands.recorded.rho09.n{n_genes}"
    R = C.replicates(SUITE, full=True)
    pair = gate_pair(n_genes, 0.9)
    sink: list[dict[str, Any]] = []
    reps = C.replicate_range(R)
    t0 = time.perf_counter()
    ind = [replicate(C.replicate_rng(stream, i), pair, [0.05], sink=sink)[0] for i in reps]
    res = C.check_fwer(ind, 0.05)
    C.log_result("islands_recorded", "fwer", 0.05, res, C.suite_seed(stream, reps.start),
                 time.perf_counter() - t0)  # fmt: skip
    _log_block_lengths(stream, 0.9, n_genes, sink)
    print(f"recorded {stream}: {res.detail}")
    assert 0.0 <= res.observed <= 1.0


# --------------------------------------------------------------------------- teeth


@pytest.mark.calibration
def test_teeth_iid_bootstrap_fails_the_gate(teeth_log: Path) -> None:
    """A known-bad null: the iid bootstrap (l' = 1) ignores the AR(1) dependence, so the null
    maximum is too small and the same fast cell must be a confirmed ``fail`` at alpha 0.05."""
    pair = gate_pair(FAST_GENES, FAST_RHO)
    rec = C.calibrate(
        SUITE, "fwer", 0.05, lambda rng: replicate(rng, pair, [0.05], block_len=1)[0],
        C.replicates(SUITE, full=False), stream="islands.teeth.iid_bootstrap",
        log_dir=teeth_log,
    )  # fmt: skip
    print(f"teeth iid bootstrap: {rec}")
    assert rec["status"] == "fail", f"the iid bootstrap is not a confirmed fail: {rec}"


# --------------------------------------------------------------------------- Example 6


@pytest.mark.calibration
def test_same_seeds_same_indicators() -> None:
    """Example 6: replicates 0..2 of the fast stream twice give identical indicator vectors
    (and identical block lengths and minimum p)."""
    levels = C.load_tolerances().levels
    pair = gate_pair(FAST_GENES, FAST_RHO)
    runs = []
    for _ in range(2):
        sink: list[dict[str, Any]] = []
        ind = [replicate(C.replicate_rng(FAST_STREAM, i), pair, levels, sink=sink)
               for i in range(3)]  # fmt: skip
        runs.append((np.stack(ind), sink))
    np.testing.assert_array_equal(runs[0][0], runs[1][0])
    assert runs[0][1] == runs[1][1]
