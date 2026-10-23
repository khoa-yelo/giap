"""Coverage of the ``giap compare`` TU-block bootstrap CI for Spearman rho (M9-22 acceptance
"Calibration"; suite ``compare``; GIAP_SPEC.md §7.9.3, §8.2, §8.3, §8.4; ADR-0042).

ADR-0042 (accepted by Khoa 2026-10-22, merged in #201) adds ``replicates_fast/full.compare`` =
100 / 1000 to ``tolerances.yaml`` and states that this file **gates**; M9-22 lead decision 2
("recorded, not gated, until a Khoa-approved ADR adds the counts") is therefore past its
condition. C1 gates (anti-conservative side, two-block confirmation through
``_criteria.calibrate``) and every claim is also recorded to
``local/logs/calibration/compare_coverage.jsonl``.

World (one replicate, ``replicate_rng("compare.rho_ci.<arm>", i)``): ``N_TU`` = 120 TUs of 1-5
genes (uniform sizes, as giapsim draws operon sizes; about 360 genes). A known common abundance
``mu = u_TU + v_gene`` (SD 1 and 0.5) is measured twice with independent TU-correlated noise:
``a = mu + e_TU^A + eps^A``, ``b = mu + e_TU^B + eps^B`` (TU noise SD 0.7, gene noise SD 0.3). Every
gene is detected and measured, so all are shared. Marginally ``(a, b)`` is bivariate normal with
``r = 1.25 / 1.83``; the target is the population Spearman ``rho_s = (6 / pi) asin(r / 2)``
(Kruskal 1958), which the sample rho estimates with O(1/n) bias.

Claims:

C1 (**gated**, suite ``compare``): the 95 % CI of ``abundance_pairs`` (``tus`` = the TUs of the
   first run, B = 1000 resamples from ``giap.rng.for_step("compare", seed)``) misses ``rho_s`` in
   at most ``level_bound(0.05, R)`` of R replicates (R = ``replicates_fast/full.compare``).
   Over-coverage beyond the binomial band is recorded and warned, not failed.
C2 (**recorded**): coverage, mean CI width and the mean bias of rho-hat.
C3 (**teeth**): the same world with each gene its own resampling unit (``tus`` = locus tags)
   must under-cover: positive within-TU correlation makes gene-level resampling too narrow.
   Without this the CI would not depend on resampling whole TUs.
"""

from __future__ import annotations

import json
import math
import time
import warnings
from collections.abc import Callable

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from helpers.compare_data import gene, genes_frame

from giap.rng import for_step

SUITE = "compare"
pytestmark = [pytest.mark.calibration, pytest.mark.calsuite(SUITE)]
# ~0.7 s per replicate (B = 1000 on ~360 genes): the fast mode passes 60 s, so it is nightly.
MODES = [pytest.param(False, marks=pytest.mark.calibration_nightly, id="fast"),
         pytest.param(True, marks=pytest.mark.calibration_full, id="full")]  # fmt: skip
LEVEL = 0.95
B = 1000
N_TU = 120
SD_U, SD_V, SD_TU_NOISE, SD_GENE_NOISE = 1.0, 0.5, 0.7, 0.3
VAR_MU = SD_U**2 + SD_V**2
R_PEARSON = VAR_MU / (VAR_MU + SD_TU_NOISE**2 + SD_GENE_NOISE**2)
RHO_S = 6.0 / math.pi * math.asin(R_PEARSON / 2.0)


def _world(rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    sizes = rng.integers(1, 6, N_TU)
    tu = np.repeat(np.arange(N_TU), sizes)
    n = tu.size
    mu = rng.normal(0.0, SD_U, N_TU)[tu] + rng.normal(0.0, SD_V, n)
    a = mu + rng.normal(0.0, SD_TU_NOISE, N_TU)[tu] + rng.normal(0.0, SD_GENE_NOISE, n)
    b = mu + rng.normal(0.0, SD_TU_NOISE, N_TU)[tu] + rng.normal(0.0, SD_GENE_NOISE, n)

    def frame(x: np.ndarray) -> pd.DataFrame:
        pairs = enumerate(zip(x, tu, strict=True))
        rows = [gene(f"g{i:04d}", start=300 * i, log2=20.0 + float(v), tu=f"TU_{t}")
                for i, (v, t) in pairs]  # fmt: skip
        return genes_frame(rows)

    return frame(a), frame(b)


def _replicate(rng: np.random.Generator, gene_level: bool) -> tuple[float, float, float]:
    """(covered 0/1 or NaN, CI width, rho-hat - rho_s) of one world."""
    from giap.compare.abundance import abundance_pairs

    fa, fb = _world(rng)
    units = fa["locus_tag"] if gene_level else fa["tu_id"]
    tus = pd.Series(units.to_numpy(), index=fa["locus_tag"].to_numpy())
    seed = int(rng.integers(0, 2**31 - 1))
    row = abundance_pairs({"A": fa, "B": fb}, tus, B, for_step("compare", seed)).iloc[0]
    lo, hi, rho = float(row["ci_low"]), float(row["ci_high"]), float(row["spearman_rho"])
    assert int(row["n_shared"]) == len(fa)
    if not (math.isfinite(lo) and math.isfinite(hi)):
        return math.nan, math.nan, rho - RHO_S
    return float(lo <= RHO_S <= hi), hi - lo, rho - RHO_S


def _record(arm: str, R: int, out: list[tuple[float, float, float]], t0: float) -> dict:
    cov = np.array([c for c, _, _ in out], dtype=float)
    ok = np.isfinite(cov)
    record = {
        "suite": SUITE,
        "stream": f"{SUITE}.rho_ci.{arm}",
        "claim": "ci_coverage",
        "statistic": "spearman_rho_tu_block_bootstrap",
        "R": R,
        "B": B,
        "n_tu": N_TU,
        "rho_s": RHO_S,
        "estimable": float(ok.mean()),
        "coverage": float(cov[ok].mean()) if ok.any() else math.nan,
        "mean_ci_width": float(np.nanmean([w for _, w, _ in out])) if ok.any() else math.nan,
        "mean_bias": float(np.mean([d for _, _, d in out])),
        "seconds": time.perf_counter() - t0,
    }
    path = C.log_dir() / "compare_coverage.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(record) + "\n")
    print(f"recorded compare rho CI {arm}: {record}")
    return record


def _arm(arm: str, R: int, gene_level: bool) -> dict:
    t0 = time.perf_counter()
    out = C.run_replicates(f"{SUITE}.rho_ci.{arm}", lambda r: _replicate(r, gene_level), R)
    return _record(arm, R, out, t0)


def test_target_is_the_population_spearman() -> None:
    """Sanity of the oracle: r = 0.683, rho_s = 0.666 (not the Pearson r)."""
    assert pytest.approx(1.25 / 1.83) == R_PEARSON
    assert pytest.approx(0.6660, abs=5e-4) == RHO_S


@pytest.mark.parametrize("full", MODES)
def test_c1_tu_block_ci_covers_rho_at_level(full: bool) -> None:
    """C1 gated + C2 recorded (from the same replicates)."""
    R = C.replicates(SUITE, full=full)
    seen: list[tuple[float, float, float]] = []
    t0 = time.perf_counter()

    def miss(rng: np.random.Generator) -> float:  # "p" = 0 when the CI misses rho_s or is NA
        out = _replicate(rng, gene_level=False)
        seen.append(out)
        return 1.0 if out[0] == 1.0 else 0.0

    fn: Callable[[np.random.Generator], float] = miss
    res = C.calibrate(SUITE, "level", 1.0 - LEVEL, fn, R, stream=f"{SUITE}.rho_ci.tu_block")
    rec = _record("tu_block", len(seen), seen, t0)
    assert rec["estimable"] == 1.0, rec  # 120 TUs: every replicate has a CI
    assert abs(rec["mean_bias"]) < 0.02, rec
    assert res["status"] != "fail", res
    half = float(C.load_tolerances().se_multiplier) * math.sqrt(LEVEL * (1 - LEVEL) / R)
    if rec["coverage"] > LEVEL + half:
        warnings.warn(f"compare rho CI over-covers: {rec['coverage']:.3f} > {LEVEL + half:.3f}"
                      f" over R={R} (recorded)", stacklevel=1)  # fmt: skip


@pytest.mark.parametrize("full", MODES)
def test_c3_teeth_gene_level_resampling_under_covers(full: bool) -> None:
    R = C.replicates(SUITE, full=full)
    rec = _arm("gene_level", R, gene_level=True)
    half = float(C.load_tolerances().se_multiplier) * math.sqrt(LEVEL * (1 - LEVEL) / R)
    assert rec["coverage"] < LEVEL - half, rec
