"""Calibration of giapsim evidence (issue M4-09; GIAP_SPEC.md §11.3, §8.4, §8.7).

Suite ``gene_fdr`` (SUITE_REGISTRY file entry): this generator is the gene-FDR null of M4-12.
Every statistical claim of docs/specs/M4-09-sim-evidence.md becomes a per-replicate value gated
through the M4-05 harness (``_criteria.calibrate``, ADR-0015 two-block confirmation). Replicate
``i`` draws a 2,000-gene layout, its expression and its evidence from
``replicate_rng(stream, i)``. pi_wrong in {0.01, 0.05} (lead Q2 proposed default: the fast suite
covers both; nightly adds 0.10).

* q-values (example 9): FDP among target rows with q <= alpha (truth ``is_false``), criterion
  ``fdr`` at every ``levels`` alpha;
* PEP: the sum of PEP over targets has, given the labels, mean n_false E0 + n_true E1 and variance
  n_false V0 + n_true V1 (exact moments of PEP under each component; oracle); z-test p-values
  are uniform (``uniform``) and reject at rate <= alpha (``level``);
* scores: true ~ N(mu1, 1), false targets and decoys ~ N(mu0, 1): one-sample KS per component;
* false targets land on gene g with probability proportional to n_observable_g (§8.4 "false
  peptides on random genes, size-proportional"; §8.7 pile-up on undetected genes): binomial test
  of the count on ``detected_true = False`` genes (discrete p, so ``level`` only); same for the
  decoys' paired genes;
* intensities: ln I - abundance_log2_true ln 2 ~ N(0, 0.5^2) for true peptides (one-sample KS).
"""

from __future__ import annotations

import math
import sys
from collections.abc import Callable
from functools import cache
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import evidence_oracle as eo

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
PIS = [0.01, 0.05]
SCORED = frozenset({"q_value", "pep", "score", "is_decoy", "intensity", "protein_ids"})
SUITE = "gene_fdr"
_MEMO: dict[tuple[str, int], Any] = {}


def _two_sided_z(z: float) -> float:
    return float(2 * stats.norm.sf(abs(z)))


def _binom_p(k: int, n: int, p: float) -> float:
    return float(stats.binomtest(k, n, p).pvalue) if n else 1.0


def _replicate(rng: np.random.Generator, pi: float) -> tuple[Any, Any, Any]:
    import giapsim.evidence as ev
    import giapsim.expression as ex
    import giapsim.layout as la

    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=2000), rng)
    expr = ex.simulate_expression(lay, ex.ExpressionConfig(), rng)
    e = ev.simulate_evidence(lay, expr, ev.EvidenceConfig(pi_wrong=pi, emit=SCORED), rng)
    return lay, expr, e


def _memo(stream: str, fn: Callable[[np.random.Generator], Any]) -> Callable[[Any], Any]:
    def run(rng: np.random.Generator) -> Any:
        key = (stream, int(rng.bit_generator.seed_seq.entropy))  # type: ignore[attr-defined]
        if key not in _MEMO:
            _MEMO[key] = fn(rng)
        return _MEMO[key]

    return run


def _labels(e: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(score, pep, q, is_false) of the target rows; decoy scores separately via ``_decoy``."""
    t = eo.targets(e.peptides)
    bad = eo.false_set(e.truth)
    return (t["score"].to_numpy(float), t["pep"].to_numpy(float), t["q_value"].to_numpy(float),
            t["peptide"].isin(bad).to_numpy())  # fmt: skip


def _gate_p(stream: str, full: bool, fn: Callable[[np.random.Generator], float],
            uniform: bool = True) -> None:  # fmt: skip
    R = C.replicates(SUITE, full)
    statuses = []
    if uniform:
        statuses.append(C.calibrate(SUITE, "uniform", 0.0, fn, R, stream=stream)["status"])
    for a in C.load_tolerances().levels:
        statuses.append(C.calibrate(SUITE, "level", a, fn, R, stream=stream)["status"])
    assert "fail" not in statuses, (stream, statuses)


@cache
def _moments(pi: float) -> tuple[float, float, float, float]:
    return eo.pep_moments(pi, 0.0, 3.0)


# ----------------------------------------------------------------------------- q-values


@pytest.mark.parametrize("pi", PIS)
@pytest.mark.parametrize("full", MODES)
def test_q_values_control_fdp(full: bool, pi: float) -> None:
    """Example 9 at every level: mean FDP of target rows with q <= alpha within alpha + k SE."""
    stream = f"sim_evidence.rep.pi{pi}"
    rep = _memo(stream, lambda rng: _labels(_replicate(rng, pi)[2]))
    R = C.replicates(SUITE, full)
    statuses = []
    for a in C.load_tolerances().levels:

        def fdp(rng: np.random.Generator, a: float = a) -> float:
            _, _, q, bad = rep(rng)
            acc = q <= a
            return float(bad[acc].mean()) if acc.any() else 0.0

        statuses.append(C.calibrate(SUITE, "fdr", a, fdp, R, stream=stream)["status"])
    assert "fail" not in statuses, statuses


@pytest.mark.parametrize("pi", PIS)
@pytest.mark.parametrize("full", MODES)
def test_pep_sum_matches_false_count(full: bool, pi: float) -> None:
    stream = f"sim_evidence.rep.pi{pi}"
    rep = _memo(stream, lambda rng: _labels(_replicate(rng, pi)[2]))
    e0, v0, e1, v1 = _moments(pi)

    def fn(rng: np.random.Generator) -> float:
        _, pep, _, bad = rep(rng)
        nf, nt = int(bad.sum()), int((~bad).sum())
        z = (float(pep.sum()) - nf * e0 - nt * e1) / math.sqrt(nf * v0 + nt * v1)
        return _two_sided_z(z)

    _gate_p(stream, full, fn)


# ----------------------------------------------------------------------------- scores


@pytest.mark.parametrize("pi", PIS)
@pytest.mark.parametrize("full", MODES)
def test_score_components(full: bool, pi: float) -> None:
    """True ~ N(3, 1); false targets and decoys ~ N(0, 1): one-sample KS p-value per component."""

    def true_p(rng: np.random.Generator) -> float:
        s, _, _, bad = _labels(_replicate(rng, pi)[2])
        return float(stats.kstest(s[~bad], "norm", args=(3.0, 1.0)).pvalue)

    def null_p(rng: np.random.Generator) -> float:
        _, _, e = _replicate(rng, pi)
        s, _, _, bad = _labels(e)
        x = np.concatenate([s[bad], eo.decoys(e.peptides)["score"].to_numpy(float)])
        return float(stats.kstest(x, "norm", args=(0.0, 1.0)).pvalue)

    _gate_p(f"sim_evidence.scores_true.pi{pi}", full, true_p)
    _gate_p(f"sim_evidence.scores_null.pi{pi}", full, null_p)


# ----------------------------------------------------------------------------- size-proportional


@pytest.mark.parametrize("full", MODES)
def test_false_targets_size_proportional(full: bool) -> None:
    """§8.4 / §8.7: share of false targets on undetected genes = their n_observable share."""
    stream = "sim_evidence.false_genes"

    def parts(rng: np.random.Generator) -> tuple[float, float]:
        lay, expr, e = _replicate(rng, 0.05)
        undet = sorted(str(t) for t, d in zip(expr.genes["locus_tag"], expr.genes["detected_true"],
                                              strict=True) if not d)  # fmt: skip
        share = eo.n_obs_share(lay, undet)
        bad = eo.false_set(e.truth)
        h = e.peptide_hits
        genes = [str(g) for p, g in zip(h["peptide"], h["gene"], strict=True) if p in bad]
        u = set(undet)
        p_false = _binom_p(sum(g in u for g in genes), len(genes), share)
        dec = e.truth.set_index("peptide").loc[eo.decoys(e.peptides)["peptide"], "decoy_gene"]
        dgenes = [str(g) for g in dec]
        p_dec = _binom_p(sum(g in u for g in dgenes), len(dgenes), share)
        return p_false, p_dec

    rep = _memo(stream, parts)
    _gate_p(f"{stream}.false", full, lambda r: rep(r)[0], uniform=False)
    _gate_p(f"{stream}.decoy", full, lambda r: rep(r)[1], uniform=False)


# ----------------------------------------------------------------------------- intensities


@pytest.mark.parametrize("full", MODES)
def test_intensity_lognormal(full: bool) -> None:
    """ln I - abundance_log2_true ln 2 ~ N(0, 0.5^2) for true peptides (one-sample KS)."""

    def fn(rng: np.random.Generator) -> float:
        _, expr, e = _replicate(rng, 0.01)
        ab = dict(zip(expr.genes["locus_tag"], expr.genes["abundance_log2_true"], strict=True))
        t = e.truth
        bad = eo.false_set(t)
        gene = {p: g for p, g in zip(t["peptide"], t["gene_true"], strict=True)
                if p not in bad and not eo.is_na(g)}  # fmt: skip
        p = eo.targets(e.peptides)
        p = p[p["peptide"].isin(set(gene))]
        x = np.log(p["intensity"].to_numpy(float))
        r = x - np.array([ab[gene[q]] * math.log(2) for q in p["peptide"]])
        return float(stats.kstest(r, "norm", args=(0.0, 0.5)).pvalue)

    _gate_p("sim_evidence.intensity", full, fn)
