"""Attacks on the M8-14 giapsim function and determinants generators (GIAP_SPEC.md §7.7.5, §8.4
rows ORA / GSEA / Determinants null, §8.5 attacks 2 and 3; statistician review of PR #251).
Statistician-owned. These test the *generators* M8-15 gates on, through the shipped engines where
that is the only way to say whether an attack is a fair null.

Arms (probe evidence: ``$GIAP_ROOT``-free scratch probes of the PR #251 review, seeds below):

* **Attack 2, fitted p_detect** (strict xfail, review B1). ``ora_attack_confounded(4000, -1.5)``
  moves the members' detection logit directly; their covariates are untouched (members are random,
  lead decision 3 gap 7). An M4-13 fit on the covariates therefore cannot see the shift, and the
  weighted ORA with the fitted p_detect rejects the confounded category in every replicate
  (probe: 30 of 30 at alpha 0.01). M8-15 gates this arm with the fitted p ("never the truth"), so
  the gate fails by construction: from the data the category *is* under-detected beyond its
  covariates, which no detectability-weighted method may explain away.
* **Attack 2, oracle p_true** (gated): the same world with ``p_detect = p_detect_true`` is at level
  (probe 0 / 0.067 / 0.10 at R = 30). This is the only arm the current generator supports.
* **ORA null, fitted p_detect** (gated, nightly): ``ora_null(4000)`` with the M4-13 fit; the
  largest in-range category, two-sided. M8-15's main arm in its real configuration.
* **Attack 3, testable after the TU collapse** (strict xfail, review B2): every
  ``gsea_attack_operon_packed`` set holds 1-3 TUs, below ``function.gsea_size`` = [10, 500], so
  the TU-collapsed engine writes no row and M8-15's "level after TU collapse" gate has m = 0.
* **Attack 3, packed TUs exchangeable at TU level** (strict xfail, review B2): the extra
  N(0, 3^2) effect makes the packed TUs' scores more dispersed than every other TU's, so a set
  built from them has a genuine TU-level effect; it is not a null for any TU-level test.
* **GSEA null with the default position decay** (gated): the scenario sets
  ``position_decay_log2 = 0``; with the M4-08 default 0.2 one random in-range set per replicate
  is still at level (probe, all sets pooled, R = 80: 0.0084 / 0.047 / 0.096). Engine:
  ``SIMPLE_ENGINE`` (structure only, no R-engine level claim; ADR-0038).
* **Determinants truth is the estimand** (gated): OLS of the response on the three planted
  features among *detected* genes recovers the recorded true R² at 0.2 and 0.5 (mean over
  20 worlds within 0.02). The unit test "recorded = analytic to 1e-6" restates the formula; this
  arm checks that the formula describes the data (features uncorrelated with each other, with the
  TU effect and with detection).
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest

pytestmark = [pytest.mark.calibration]

LEVELS = (0.01, 0.05, 0.10)
GSEA_SIZE = (10, 500)


def _fs() -> Any:
    import giapsim.function_sim as fs

    return fs


def _sc() -> Any:
    import giapsim.scenarios as sc

    return sc


def _coll(name: str, sets: dict[str, set[str]]) -> Any:
    from giap.function.collections import GeneSets

    s = {t: frozenset(g) for t, g in sorted(sets.items())}
    return GeneSets(name, s, {t: t for t in s}, {}, 0.5, 0)  # type: ignore[arg-type]


def _tu_maps(w: Any) -> tuple[dict[str, str], dict[str, int]]:
    tu: dict[str, str] = {}
    pos: dict[str, int] = {}
    for t, gs in zip(w.layout.tus["tu_id"], w.layout.tus["genes"], strict=True):
        for k, g in enumerate(gs, start=1):
            tu[str(g)], pos[str(g)] = str(t), k
    return tu, pos


def _members(w: Any) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for s, g in zip(w.sets["set_id"], w.sets["locus_tag"], strict=True):
        out.setdefault(str(s), set()).add(str(g))
    return out


def _gate(suite: str, stream: str, fn: Any, R: int, expect: str = "pass") -> None:
    statuses = [C.calibrate(suite, "level", a, fn, R, stream=stream)["status"] for a in LEVELS]
    if expect == "pass":
        assert "fail" not in statuses, (stream, statuses)
    else:
        assert "fail" in statuses, (stream, statuses)


# ----------------------------------------------------------------------------- ORA


def _ora_inputs(w: Any, fitted: bool, rng: np.random.Generator) -> tuple[pd.DataFrame, Any]:
    """The ORA gene frame (eligible genes: a protein sequence) and M4-13 folds."""
    from giap.config.model import DetectabilityConfig
    from giap.genetable.detectability import fit_detectability, tu_folds

    g = w.expression.genes
    tu, _ = _tu_maps(w)
    tags = g["locus_tag"].astype(str)
    length = w.layout.genes.set_index("locus_tag").loc[tags, "aa"].str.len().to_numpy()
    ok = length > 0
    det = pd.Series(g["detected_true"].to_numpy(bool), index=tags, dtype="boolean")
    det[~ok] = pd.NA
    tus = pd.Series([tu[t] for t in tags], index=tags)
    frame = pd.DataFrame({"locus_tag": tags[ok].to_numpy(), "tu_id": tus[ok].to_numpy(),
                          "partial": False,
                          "detected": det[ok].astype(bool).to_numpy()})  # fmt: skip
    if fitted:
        cov = pd.DataFrame({"locus_tag": tags, "length_aa": length,
                            "n_observable_peptides": g["n_obs_true"].to_numpy(),
                            "gravy": g["gravy_true"].to_numpy()})  # fmt: skip
        cfg = DetectabilityConfig(covariates=["n_observable_peptides", "length_aa", "gravy"])
        fit = fit_detectability(det, cov, tus, cfg, rng)
        assert not fit.constant
        return frame.assign(p_detect=fit.p_detect.reindex(frame["locus_tag"]).to_numpy()), (
            fit.fold_of_tu
        )
    p = g["p_detect_true"].to_numpy(float)[ok]
    return frame.assign(p_detect=p), tu_folds(list(frame["tu_id"]), 5, rng)


def _ora_two_sided(genes: pd.DataFrame, members: set[str], folds: Any) -> float:
    from giap.function import ora

    tab = ora.detection_enrichment(genes, [_coll("kegg", {"S": members})], (5, 10**6), 0.05,
                                   folds=folds)  # fmt: skip
    assert "tu_recal_skipped" not in tab.attrs["run_flags"]
    p = tab["p_value"].to_numpy(float)
    assert p.size == 2 and np.all(np.isfinite(p))
    return float(min(1.0, 2 * p.min()))


def _attack2_rep(fitted: bool) -> Any:
    def fn(rng: np.random.Generator) -> float:
        w = _fs().draw(_sc().ora_attack_confounded(4000, -1.5), rng)
        s = w.sets
        conf = set(s.loc[s["confounded"].astype(bool), "locus_tag"].astype(str))
        genes, folds = _ora_inputs(w, fitted, rng)
        return _ora_two_sided(genes, conf & set(genes["locus_tag"]), folds)

    return fn


@pytest.mark.calsuite("ora")
@pytest.mark.calibration_nightly
@pytest.mark.xfail(strict=True, reason="PR #251 B1: the member shift is invisible to covariates")
def test_attack2_fitted_p_detect_confounded_category_at_level() -> None:
    _gate("ora", "attack.sim_function.ora_attack.fitted", _attack2_rep(True), 20)


@pytest.mark.calsuite("ora")
@pytest.mark.calibration_nightly
def test_attack2_oracle_p_true_confounded_category_at_level() -> None:
    _gate("ora", "attack.sim_function.ora_attack.oracle", _attack2_rep(False),
          C.replicates("ora", False))  # fmt: skip


@pytest.mark.calsuite("ora")
@pytest.mark.calibration_nightly
@pytest.mark.slow
def test_ora_null_fitted_p_detect_largest_category_at_level() -> None:
    def fn(rng: np.random.Generator) -> float:
        w = _fs().draw(_sc().ora_null(4000), rng)
        genes, folds = _ora_inputs(w, True, rng)
        tu = dict(zip(genes["locus_tag"], genes["tu_id"], strict=True))
        sets = {k: v & set(tu) for k, v in _members(w).items()}
        ok = {k: v for k, v in sets.items() if 5 <= len({tu[g] for g in v}) <= 500}
        biggest = ok[max(sorted(ok), key=lambda k: len(ok[k]))]
        return _ora_two_sided(genes, biggest, folds)

    _gate("ora", "attack.sim_function.ora_null.fitted", fn, C.replicates("ora", False))


# ----------------------------------------------------------------------------- GSEA


def _gsea_genes(w: Any) -> pd.DataFrame:
    g = w.expression.genes
    tu, pos = _tu_maps(w)
    tags = g["locus_tag"].astype(str)
    det = g["detected_true"].to_numpy(bool)
    return pd.DataFrame({
        "locus_tag": tags, "tu_id": [tu[t] for t in tags], "tu_position": [pos[t] for t in tags],
        "detected": det, "abundance_log2": np.where(det, g["abundance_log2_true"], np.nan),
        "abundance_imputed": False, "abundance_scale": "molar", "partial": False,
        "p_detect": g["p_detect_true"].to_numpy(float),
    })  # fmt: skip


@pytest.mark.calsuite("gsea")
@pytest.mark.xfail(strict=True, reason="PR #251 B2: packed sets are 1-3 TUs, below gsea_size")
def test_attack3_packed_sets_are_tested_after_tu_collapse() -> None:
    from giap.config.model import GseaConfig
    from giap.function import gsea
    from giap.rng import for_step

    w = _fs().draw(_sc().gsea_attack_operon_packed(4000, 20), np.random.default_rng(7))
    tab = gsea.gsea_table(_gsea_genes(w), [_coll("kegg", _members(w))], GSEA_SIZE, GseaConfig(),
                          0.05, for_step("function.gsea", 1))  # fmt: skip
    assert int(np.isfinite(tab["p_value"].to_numpy(float)).sum()) >= 10, len(tab)


@pytest.mark.calsuite("gsea")
@pytest.mark.xfail(strict=True, reason="PR #251 B2: packed TUs carry an extra N(0, 3^2) effect")
def test_attack3_packed_tus_exchangeable_with_other_tus() -> None:
    """Fligner-Killeen test of TU-score dispersion, packed TUs vs the rest (one p per world)."""
    from scipy import stats

    from giap.function.gsea import tu_ranking

    def fn(rng: np.random.Generator) -> float:
        w = _fs().draw(_sc().gsea_attack_operon_packed(4000, 20), rng)
        rank = tu_ranking(_gsea_genes(w))
        tu, _ = _tu_maps(w)
        packed = {tu[g] for g in w.sets["locus_tag"].astype(str)}
        inside = np.array([t in packed for t in rank.tu_ids])
        return float(stats.fligner(rank.scores[inside], rank.scores[~inside]).pvalue)

    _gate("gsea", "attack.sim_function.packed.exchangeable", fn, C.replicates("gsea", False))


@pytest.mark.calsuite("gsea")
def test_gsea_null_with_default_position_decay_at_level() -> None:
    from giap.config.model import GseaConfig
    from giap.function import gsea
    from giap.rng import for_step

    def fn(rng: np.random.Generator) -> float:
        lcfg, ecfg = _sc().gsea_null(4000)
        w = _fs().draw((lcfg, replace(ecfg, position_decay_log2=0.2)), rng)
        genes = _gsea_genes(w)
        rank = gsea.tu_ranking(genes)
        tu, _ = _tu_maps(w)
        scoring = {g for gs in rank.tu_genes.values() for g in gs}  # membership: scoring genes
        sets = {
            k: v
            for k, v in _members(w).items()
            if GSEA_SIZE[0] <= len({tu[g] for g in v & scoring}) <= GSEA_SIZE[1]
        }
        name = sorted(sets)[int(rng.integers(0, len(sets)))]
        step = for_step("function.gsea", int(rng.integers(0, 2**31 - 1)))
        tab = gsea.gsea_table(genes, [_coll("kegg", {name: sets[name]})], GSEA_SIZE,
                              GseaConfig(), 0.05, step)  # fmt: skip
        p = tab["p_value"].to_numpy(float)
        assert p.size == 1 and np.isfinite(p[0])
        return float(p[0])

    _gate("gsea", "attack.sim_function.gsea_null.decay0.2", fn, C.replicates("gsea", False))


# ----------------------------------------------------------------------------- determinants


@pytest.mark.calsuite("determinants")
@pytest.mark.parametrize(
    ("coefs", "r2"),
    [({"promoter": 0.5, "rbs": 0.4, "cai": 0.3}, 0.2),
     ({"promoter": 1.0, "rbs": 0.8, "cai": 0.6}, 0.5)],
)  # fmt: skip
def test_determinants_true_r2_is_recovered_among_detected_genes(
    coefs: dict[str, float], r2: float
) -> None:
    fs = _fs()
    est = []
    for i in range(20):
        w = fs.draw(_sc().determinants_planted(3000, coefs, 1.0),
                    C.replicate_rng("attack.sim_function.determinants.r2", i))  # fmt: skip
        true = float(w.truth["determinants"]["true_r2"].iloc[0])
        assert abs(true - r2) < 1e-12
        g = w.expression.genes.set_index("locus_tag")
        f = w.features.set_index("locus_tag")
        det = g.loc[f.index, "detected_true"].to_numpy(bool)
        y = g.loc[f.index, "abundance_log2_true"].to_numpy(float)[det]
        x = np.column_stack([np.ones(det.sum())]
                            + [f[c].to_numpy(float)[det] for c in fs.FEATURES])  # fmt: skip
        b, *_ = np.linalg.lstsq(x, y, rcond=None)
        est.append(1.0 - float(np.var(y - x @ b)) / float(np.var(y)))
    assert abs(float(np.mean(est)) - r2) < 0.02, (r2, np.mean(est), np.std(est))
