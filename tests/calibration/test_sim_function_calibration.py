"""Calibration of the gippasim function and determinants generators (issue M8-14; GIPPA_SPEC.md
§8.4 rows ORA, GSEA, Determinants null; §8.5 attack 2; docs/specs/M8-14-sim-function-
determinants.md examples 1, 2, 3 and 6; interface as assumed in tests/unit/test_sim_function.py).

Each claim is a null that must hold *in the generator*; replicate ``i`` draws one world from
``replicate_rng(stream, i)`` and returns one p-value; gated through the M4-05 harness (§8.4 level
at every alpha of ``levels``, ADR-0015 two-block confirmation) in the suite whose nulls the
generator feeds (registered by ``calsuite`` markers): ``ora`` (R 100 / 1,000), ``gsea``
(50 / 300), ``determinants`` (20 / 100).

- Example 1 (``ora_null(4000)``): detection independent of category membership given p_true.
  For the largest category S, Z = Σ_S (D - p) / sqrt(Σ_S p (1 - p)) is N(0, 1) under the null
  (D ~ Bernoulli(p) independently, membership independent of D given p); p = 2 Φ(-|Z|) (the
  spec's χ²₁ = Z²).
- Example 2 (``ora_attack_confounded(4000, -1.5)``): detection still depends only on p_true — the
  same Z over the confounded members. Teeth: the unweighted 2x2 χ² (members' detection rate vs
  the rest, ignoring p_true) must reject (a confirmed ``fail``) — the attack bites.
- Example 3 (``gsea_null(4000)``): set membership independent of abundance. TU-level rank test
  (Mann-Whitney, two-sided) of the TU mean ``abundance_log2_true`` of the first set's TUs (by
  ``set_id``) against the other TUs: exact conditional on the abundances whenever the set's TUs
  are exchangeable given its size (TU-correlated abundance does not break it; gene-level ranking
  would).
- Example 6 (``determinants_null(3000)``): every feature independent of the response. TU-level
  Spearman (first gene of each TU, genome-order ``genes[0]``) of each feature against the TU mean
  response (features of distinct TUs exchangeable); the three p-values of one world are pooled
  into one level check per replicate (``check_level_pooled``, m = 3; lead decision 4, B3: nine
  separate full-tier checks failed by multiplicity, about 0.6 % each, at fixed seeds).

Cost (fast tier, 1 core): ~0.4 s per world at 4,000 genes; ~1 min for ``ora`` per arm, ~20 s for
``gsea``, ~10 s per feature for ``determinants``; the teeth arm runs both blocks. The full tier
(``calibration_full``) is for the nightly Slurm array only.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from scipy import stats

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
FEATURES = ("promoter", "rbs", "cai")
STREAM_DET_POOLED = "sim_function.determinants_null.pooled"


def _draw(pair: Any, rng: np.random.Generator) -> Any:
    import gippasim.function_sim as fs

    return fs.draw(pair, rng)


def _sc() -> Any:
    import gippasim.scenarios

    return gippasim.scenarios


def _members(world: Any) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for s, g in zip(world.sets["set_id"], world.sets["locus_tag"], strict=True):
        out.setdefault(str(s), set()).add(str(g))
    return out


def _z_given_p(genes: pd.DataFrame, members: set[str]) -> float:
    """Two-sided p of Σ_S (D - p) = 0 with Var Σ_S p (1 - p) (normal approximation)."""
    s = genes[genes["locus_tag"].isin(members)]
    p = s["p_detect_true"].to_numpy(float)
    d = s["detected_true"].to_numpy(bool).astype(float)
    z = float((d - p).sum() / np.sqrt((p * (1 - p)).sum()))
    return float(2 * stats.norm.sf(abs(z)))


def _tu_frame(world: Any) -> pd.DataFrame:
    tus = world.layout.tus
    rows = [
        (str(t), [str(g) for g in gs]) for t, gs in zip(tus["tu_id"], tus["genes"], strict=True)
    ]
    return pd.DataFrame(rows, columns=["tu_id", "genes"])


def _gate(suite: str, stream: str, full: bool, fn: Callable[[np.random.Generator], float],
          expect: str = "pass", log_dir: Path | None = None) -> None:  # fmt: skip
    R = C.replicates(suite, full)
    statuses = [C.calibrate(suite, "level", a, fn, R, stream=stream, log_dir=log_dir)["status"]
                for a in C.load_tolerances().levels]  # fmt: skip
    C.skip_without_verdict([{"status": s} for s in statuses])  # a shard: --merge decides
    if expect == "pass":
        assert "fail" not in statuses, (stream, statuses)
    else:
        assert "fail" in statuses, (stream, statuses)


# ----------------------------------------------------------------------------- ora


@pytest.mark.calsuite("ora")
@pytest.mark.parametrize("full", MODES)
def test_example1_ora_null_detection_independent_of_membership(full: bool) -> None:
    def fn(rng: np.random.Generator) -> float:
        w = _draw(_sc().ora_null(4000), rng)
        biggest = max(_members(w).values(), key=len)
        return _z_given_p(w.expression.genes, biggest)

    _gate("ora", "sim_function.ora_null.largest_category", full, fn)


@pytest.mark.calsuite("ora")
@pytest.mark.parametrize("full", MODES)
def test_example2_attack_detection_depends_only_on_p_true(full: bool) -> None:
    def fn(rng: np.random.Generator) -> float:
        w = _draw(_sc().ora_attack_confounded(4000, -1.5), rng)
        s = w.sets
        conf = set(s.loc[s["confounded"].astype(bool), "locus_tag"].astype(str))
        return _z_given_p(w.expression.genes, conf)

    _gate("ora", "sim_function.ora_attack.confounded_given_p", full, fn)


@pytest.mark.calsuite("ora")
@pytest.mark.parametrize("full", MODES)
def test_example2_teeth_unweighted_test_rejects(full: bool, teeth_log: Path) -> None:
    """Without p_true (ω = 1) the confounded members' detection rate differs: confirmed fail."""

    def fn(rng: np.random.Generator) -> float:
        w = _draw(_sc().ora_attack_confounded(4000, -1.5), rng)
        s = w.sets
        conf = set(s.loc[s["confounded"].astype(bool), "locus_tag"].astype(str))
        g = w.expression.genes
        inside = g["locus_tag"].isin(conf).to_numpy()
        det = g["detected_true"].to_numpy(bool)
        table = [[int((inside & det).sum()), int((inside & ~det).sum())],
                 [int((~inside & det).sum()), int((~inside & ~det).sum())]]  # fmt: skip
        return float(stats.chi2_contingency(table, correction=False)[1])

    _gate(
        "ora",
        "sim_function.ora_attack.teeth_unweighted",
        full,
        fn,
        expect="fail",
        log_dir=teeth_log,
    )  # teeth: unsharded in shard 0 (conftest), private log


# ----------------------------------------------------------------------------- gsea


@pytest.mark.calsuite("gsea")
@pytest.mark.parametrize("full", MODES)
def test_example3_gsea_null_membership_independent_of_abundance(full: bool) -> None:
    def fn(rng: np.random.Generator) -> float:
        w = _draw(_sc().gsea_null(4000), rng)
        ab = dict(zip(w.expression.genes["locus_tag"].astype(str),
                      w.expression.genes["abundance_log2_true"], strict=True))  # fmt: skip
        tus = _tu_frame(w)
        tu_ab = np.array([np.mean([ab[g] for g in gs]) for gs in tus["genes"]])
        first = _members(w)[min(_members(w))]
        inside = np.array([any(g in first for g in gs) for gs in tus["genes"]])
        res = stats.mannwhitneyu(tu_ab[inside], tu_ab[~inside], alternative="two-sided")
        return float(res.pvalue)

    _gate("gsea", "sim_function.gsea_null.first_set_tu_rank", full, fn)


# ----------------------------------------------------------------------------- determinants


def _determinants_null_ps(rng: np.random.Generator) -> tuple[float, float, float]:
    """The three feature p-values of one ``determinants_null(3000)`` world."""
    w = _draw(_sc().determinants_null(3000), rng)
    y = dict(zip(w.expression.genes["locus_tag"].astype(str),
                 w.expression.genes["abundance_log2_true"], strict=True))  # fmt: skip
    tus = _tu_frame(w)
    r = np.array([np.mean([y[g] for g in gs]) for gs in tus["genes"]])
    out = []
    for feature in FEATURES:
        f = dict(zip(w.features["locus_tag"].astype(str), w.features[feature], strict=True))
        x = np.array([f[gs[0]] for gs in tus["genes"]], dtype=float)
        out.append(float(stats.spearmanr(x, r).pvalue))
    return out[0], out[1], out[2]


@pytest.mark.calsuite("determinants")
@pytest.mark.parametrize("full", MODES)
def test_example6_determinants_null_features_independent_of_response(full: bool) -> None:
    """B3 (lead decision 4): one pooled check per replicate over the three features
    (``check_level_pooled``, m = 3; the planted features are independent of each other), in
    place of one check per feature. The world of replicate ``i`` is drawn once and reused by the
    three alphas."""
    cache: dict[Any, tuple[float, float, float]] = {}

    def ps(rng: np.random.Generator) -> tuple[float, float, float]:
        key = rng.bit_generator.seed_seq.entropy  # type: ignore[attr-defined]
        if key not in cache:
            cache[key] = _determinants_null_ps(rng)
        return cache[key]

    R = C.replicates("determinants", full)
    statuses = []
    for a in C.load_tolerances().levels:

        def share(rng: np.random.Generator, a: float = a) -> float:
            return sum(p <= a for p in ps(rng)) / len(FEATURES)

        rec = C.calibrate("determinants", "level", a, share, R, stream=STREAM_DET_POOLED,
                          m=len(FEATURES))  # fmt: skip
        assert rec["m"] == len(FEATURES)
        statuses.append(rec["status"])
    assert "fail" not in statuses, statuses
