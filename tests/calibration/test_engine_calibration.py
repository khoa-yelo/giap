"""Engine calibration through the harness (issue M4-05; GIAP_SPEC.md §8.4).

Each test runs fast (``calibration``, ``replicates_fast``) and full (``calibration_full``,
``replicates_full``). Engine tests reuse the R of the closest consumer suite (M4-index lead
decision 6): ``gene_fdr`` for BH, ``regions`` for permutation, circular shifts and Besag-Clifford.
Replicate i of stream ``<suite>.<test>`` uses ``_criteria.replicate_rng`` (seed
``suite_seed(stream, i)``). Every level is gated by ``_criteria.calibrate`` (ADR-0015 two-block
confirmation, sharding) through ``_gate.calibrate_levels``, which logs to
``local/logs/calibration/<suite>.jsonl``.

- ``test_bh_fdr``: 1,000 p per replicate, 900 U(0, 1) nulls + 100 Beta(0.1, 1) alternatives;
  mean FDP of BH at q <= q + k SE for q in ``levels``.
- ``test_perm_p_level``: two groups of 30 iid N(0, 1), stat = difference of means, B = 199 label
  permutations; level check (discrete p, no KS).
- ``test_circular_shift_level_iid``: 1,000 iid N(0, 1) on a circular index, region = 20 consecutive
  indices, stat = sum, exhaustive shifts; level check.
- ``test_besag_clifford_level``: as ``test_perm_p_level`` with |difference|, h = 10, max_B = 1,000.

Autocorrelation engine (issue M5-02, suite ``autocorrelation``; docs/specs/M5-02 acceptance,
§7.3.6 "autocorrelation tests at level under the TU-block null", ADR-0020 decision 3):

- ``test_moran_lag1_level_label_null``: iid N(0, 1), n = 500, circular; ``lag_test`` of
  ``morans_i`` lag 1 under B = 199 label permutations; level at every alpha.
- ``test_moran_lag1_level_tu_null``: TU-structured noise without inter-TU structure (iid TU sizes
  1 + Geometric(0.55), TU effect N(0, 1) + gene noise N(0, 0.5^2), ~500 genes, circular) under
  B = 199 TU-block permutations; level.
- ``test_residual_join_count_level_label_null``: p_i ~ U(0.05, 0.95), d_i ~ Bernoulli(p_i), n = 500;
  residual BB lag 1 with p permuted together with d (decision 7); level.
- ``test_example11_iid_block_len_one`` / ``test_example12_ar09_block_len_twenty``: spec examples
  11-12 on 100 seeds (``replicate_rng``), n = 2,000, B = 199, circular (example 12 as amended by
  lead decision 9: l-hat >= 15 in >= 95 %).
- ``test_decision8_giapsim_tu_scores_block_len`` (suite ``islands``, lead decision 8): giapsim
  chromosome at AR(1) rho 0.9 -> Blom scores (left-censored) -> TU means -> 10 % LOESS detrend
  (circular padding, statsmodels as M5-01 specifies) -> ``choose_block_len`` with B = 500 gives
  mean l-hat ≈ 8 at 4,500 genes and ≈ 1.3 at 1,000 genes (statistician, Slurm 1860281 / 1860318:
  8.2-8.3 over 200 replicates, 1.34 over 300), pinning the dilution ADR-0021 decision 2 rests on.
"""

from __future__ import annotations

import importlib
import math
from collections.abc import Callable
from types import ModuleType

import _criteria as C
import numpy as np
import pytest
from _gate import Gate, calibrate_levels

from giap.stats.multitest import bh
from giap.stats.permutation import (
    apply_shift,
    besag_clifford,
    choose_shifts,
    circular_shifts,
    perm_p,
)

FAST_SECONDS = 60.0
N_GROUP = 30

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def _diff(x: np.ndarray, labels: np.ndarray) -> float:
    return float(x[labels].mean() - x[~labels].mean())


def _level(
    suite: str, stream: str, full: bool, replicate: Callable[[np.random.Generator], float]
) -> None:
    r = C.replicates(suite, full=full)

    def checked(rng: np.random.Generator) -> float:
        p = replicate(rng)
        assert 0 < p <= 1, p
        return p

    gate = calibrate_levels(suite, "level", C.load_tolerances().levels, checked, r, stream=stream)
    _finish(gate, r, full)


def _finish(gate: Gate, r: int, full: bool) -> None:
    assert all(rec["R"] == r for rec in gate.records)
    gate.assert_ok()
    if not full:  # one simulation per replicate; block-1 time kept if a level needs block 2
        duration = gate.block_seconds.get(1, 0.0)
        assert duration < FAST_SECONDS, duration


@pytest.mark.parametrize("full", MODES)
def test_bh_fdr(full: bool) -> None:
    suite, stream = "gene_fdr", "gene_fdr.engine_bh"
    levels = C.load_tolerances().levels
    r = C.replicates(suite, full=full)

    def rep(rng: np.random.Generator) -> np.ndarray:
        p = np.concatenate([rng.uniform(size=900), rng.beta(0.1, 1.0, size=100)])
        null = np.arange(1000) < 900
        q = bh(p)
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = q <= a
            n = int(called.sum())
            out[j] = (called & null).sum() / n if n else 0.0
        return out

    _finish(calibrate_levels(suite, "fdr", levels, rep, r, stream=stream), r, full)


@pytest.mark.parametrize("full", MODES)
def test_perm_p_level(full: bool) -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP
        null = np.array([_diff(x, rng.permutation(labels)) for _ in range(199)])
        return perm_p(_diff(x, labels), null)

    _level("regions", "regions.engine_perm_p", full, rep)


@pytest.mark.parametrize("full", MODES)
def test_circular_shift_level_iid(full: bool) -> None:
    n, width = 1000, 20

    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=n)
        idx = (int(rng.integers(0, n)) + np.arange(width)) % n
        shifts, exhaustive = choose_shifts(
            circular_shifts(n, idx, circular=True), 20_000, n, 1_000, rng
        )
        assert bool(exhaustive) and len(shifts) == n - 1
        null = np.array([x[apply_shift(idx, int(k), n)].sum() for k in shifts])
        return perm_p(float(x[idx].sum()), null)

    _level("regions", "regions.engine_circular_shift", full, rep)


@pytest.mark.parametrize("full", MODES)
def test_besag_clifford_level(full: bool) -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP

        def stat(g: np.random.Generator) -> float:
            return abs(_diff(x, g.permutation(labels)))

        p, b_used = besag_clifford(stat, abs(_diff(x, labels)), 10, 1_000, rng)
        assert 10 <= b_used <= 1_000
        return p

    _level("regions", "regions.engine_besag_clifford", full, rep)


# ------------------------------------------------------------------ M5-02 autocorrelation

AUTOCORR = pytest.mark.calsuite("autocorrelation")


def _ac() -> ModuleType:
    """giap.stats.autocorr, imported per test so the M4 engine tests above still collect."""
    return importlib.import_module("giap.stats.autocorr")


N_AC, B_AC = 500, 199


@AUTOCORR
@pytest.mark.parametrize("full", MODES)
def test_moran_lag1_level_label_null(full: bool) -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=N_AC)
        orders = _ac().label_permutations(N_AC, B_AC, rng)
        null = _ac().morans_i_matrix(x[orders], [1], True)[:, 0]
        return _ac().lag_test(_ac().morans_i(x, 1, True), null).p_value

    _level("autocorrelation", "autocorrelation.engine_moran_label", full, rep)


def _tu_noise(rng: np.random.Generator, n_genes: int) -> tuple[np.ndarray, list[str]]:
    sizes: list[int] = []
    while sum(sizes) < n_genes:
        sizes.append(int(rng.geometric(0.55)))
    x = np.concatenate([rng.normal() + 0.5 * rng.normal(size=s) for s in sizes])
    labels = [f"TU{t}" for t, s in enumerate(sizes) for _ in range(s)]
    return x, labels


@AUTOCORR
@pytest.mark.parametrize("full", MODES)
def test_moran_lag1_level_tu_null(full: bool) -> None:
    def rep(rng: np.random.Generator) -> float:
        x, labels = _tu_noise(rng, N_AC)
        orders = _ac().tu_block_permutations(_ac().tu_blocks(labels, True), B_AC, rng)
        null = _ac().morans_i_matrix(x[orders], [1], True)[:, 0]
        return _ac().lag_test(_ac().morans_i(x, 1, True), null).p_value

    _level("autocorrelation", "autocorrelation.engine_moran_tu", full, rep)


@AUTOCORR
@pytest.mark.parametrize("full", MODES)
def test_residual_join_count_level_label_null(full: bool) -> None:
    def rep(rng: np.random.Generator) -> float:
        p = rng.uniform(0.05, 0.95, size=N_AC)
        d = rng.uniform(size=N_AC) < p
        orders = _ac().label_permutations(N_AC, B_AC, rng)
        null = _ac().join_count_matrix(d[orders], [1], True, P=p[orders])[:, 0]
        return _ac().lag_test(_ac().join_count(d, 1, True, p=p), null).p_value

    _level("autocorrelation", "autocorrelation.engine_join_count_residual", full, rep)


def _block_lens(stream: str, series: Callable[[np.random.Generator], np.ndarray]) -> np.ndarray:
    out = []
    for i in range(100):
        rng = C.replicate_rng(stream, i)
        out.append(_ac().choose_block_len(series(rng), True, 199, rng))
    return np.array(out)


@AUTOCORR
@pytest.mark.calibration
def test_example11_iid_block_len_one() -> None:
    """Spec example 11: iid N(0, 1), n = 2,000 -> l = 1 in >= 90 % of 100 seeds."""
    ells = _block_lens("autocorrelation.engine_block_len_iid", lambda r: r.normal(size=2_000))
    assert np.mean(ells == 1) >= 0.90, np.bincount(ells)


def _ar1(rng: np.random.Generator, n: int, rho: float) -> np.ndarray:
    from scipy.signal import lfilter

    e = rng.normal(size=n)
    e[0] /= math.sqrt(1.0 - rho * rho)
    return np.asarray(lfilter([1.0], [1.0, -rho], e))


@AUTOCORR
@pytest.mark.calibration
def test_example12_ar09_block_len_twenty() -> None:
    """Spec example 12 as amended by lead decision 9: AR(1) rho 0.9, n = 2,000 -> l-hat >= 15 in
    >= 95 % of 100 seeds (l-hat = 20 only in 0.85-0.92: the sampling SD of I_17 ... I_20 is about
    rho^20)."""
    ells = _block_lens("autocorrelation.engine_block_len_ar09", lambda r: _ar1(r, 2_000, 0.9))
    assert np.mean(ells >= 15) >= 0.95, np.bincount(ells, minlength=21)[1:]


def _tu_score_series(rng: np.random.Generator, n_genes: int, rho: float) -> np.ndarray:
    """Detrended TU scores of the giapsim chromosome (M5-06 steps 1-3, built here from scipy and
    statsmodels so the test does not depend on M5-01)."""
    from scipy.special import ndtri
    from scipy.stats import rankdata
    from statsmodels.nonparametric.smoothers_lowess import lowess

    import giapsim.expression as ex
    import giapsim.layout as la

    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=n_genes), rng)
    expr = ex.simulate_expression(lay, ex.ExpressionConfig(ar1_rho=rho), rng)
    g = lay.genes[["locus_tag", "contig", "start", "partial"]].merge(
        expr.genes[["locus_tag", "detected_true", "abundance_log2_true"]], on="locus_tag"
    )
    tu_of = {
        gene: t for t, gs in zip(lay.tus["tu_id"], lay.tus["genes"], strict=True) for gene in gs
    }
    g["tu"] = g["locus_tag"].map(tu_of)
    chrom = str(lay.replicons["contig"].iloc[0])
    g = g[(g["contig"] == chrom) & ~g["partial"].astype(bool)]
    g = g.sort_values(["start", "locus_tag"], kind="mergesort")
    cens = ~g["detected_true"].to_numpy(dtype=bool)
    n, n0 = len(g), int(cens.sum())
    r = np.empty(n)
    r[cens] = (n0 + 1) / 2.0
    r[~cens] = n0 + rankdata(g["abundance_log2_true"].to_numpy()[~cens])
    g["score"] = ndtri((r - 0.375) / (n + 0.25))
    tus = g.groupby("tu", sort=False).agg(z=("score", "mean"), first=("start", "min"))
    z = tus.sort_values("first", kind="mergesort")["z"].to_numpy()
    T = z.size
    k = min(T, max(3, math.ceil(0.10 * T - 1e-9)))
    pos = np.arange(T, dtype=float)
    pos_pad = np.concatenate([pos[-k:] - T, pos, pos[:k] + T])
    z_pad = np.concatenate([z[-k:], z, z[:k]])
    fit = lowess(z_pad, pos_pad, frac=k / pos_pad.size, it=0, delta=0.0, return_sorted=False)
    return np.asarray(z - fit[k : k + T])


@pytest.mark.calsuite("islands")
@pytest.mark.calibration
def test_decision8_giapsim_tu_scores_block_len() -> None:
    """Lead decision 8: l-hat ≈ 8 at 4,500 genes and ≈ 1.3 at 1,000 genes (rho 0.9). Bands are
    ±3.5 SE of a 30-replicate mean (SD ≈ 2.0 and ≈ 0.65 in the statistician's runs)."""
    means = {}
    for n_genes in (1_000, 4_500):
        ells = []
        for i in range(30):
            rng = C.replicate_rng(f"islands.engine_block_len_giapsim_{n_genes}", i)
            det = _tu_score_series(rng, n_genes, 0.9)
            ells.append(_ac().choose_block_len(det, True, 500, rng))
        means[n_genes] = float(np.mean(ells))
    assert 6.7 <= means[4_500] <= 9.5, means
    assert 1.0 <= means[1_000] <= 1.8, means
    assert means[4_500] - means[1_000] >= 4.0, means
