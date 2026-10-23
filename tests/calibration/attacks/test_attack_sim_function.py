"""Attacks on the M8-14 giapsim function and determinants generators (GIAP_SPEC.md §7.7.5, §8.4
rows ORA / GSEA / Determinants null, §8.5 attacks 2 and 3; statistician review of PR #251).
Statistician-owned. These test the *generators* M8-15 gates on, through the shipped engines where
that is the only way to say whether an attack is a fair null.

Arms (probe evidence: ``$GIAP_ROOT``-free scratch probes of the PR #251 review, seeds below):

* **Attack 2, fitted p_detect** (gated, nightly; review B1, lead decision 4). Before the fix
  ``ora_attack_confounded(4000, -1.5)`` moved the members' detection logit directly, leaving their
  covariates untouched, so an M4-13 fit could not see the shift and the weighted ORA with the
  fitted p_detect rejected the confounded category in every replicate (probe: 30 of 30 at alpha
  0.01). Decision 4 carries the shift through an M4-13 covariate (default GRAVY, recorded as
  ``function_sim.CONFOUND_COVARIATE``; members still random, gap 7), so the weighted ORA with the
  fitted p_detect must now stay at level. Teeth: the unweighted 2x2 test of the calibration file
  (``test_example2_teeth_unweighted_test_rejects``).
* **Attack 2, oracle p_true** (gated): the same world with ``p_detect = p_detect_true`` is at level
  (probe 0 / 0.067 / 0.10 at R = 30). This is the only arm the current generator supports.
* **ORA null, fitted p_detect** (gated, nightly): ``ora_null(4000)`` with the M4-13 fit; the
  largest in-range category, two-sided. M8-15's main arm in its real configuration.
* **Attack 3, testable after the TU collapse** (gated; review B2, lead decision 4): before the
  fix every ``gsea_attack_operon_packed`` set held 1-3 TUs, below ``function.gsea_size`` =
  [10, 500], so the TU-collapsed engine wrote no row. Now each set is 1-3 packed operons of
  10-30 genes plus single-gene TUs, >= 10 TUs: every set gets a finite p.
* **Attack 3, packed TUs exchangeable at TU level** (gated; B2): the extra N(0, 3^2) effect is
  dropped. Fligner-Killeen test of the true TU level (``tu_log2``), packed operons against the
  other chr1 TUs (the TU *score* would mix in the gene-noise average, which depends on TU size).
* **Attack 3, TU-level GSEA at level** (gated; B2): one random packed set per world through
  ``gsea_table`` (``SIMPLE_ENGINE``). Teeth: the same set ranked at gene level (genes permuted
  as if independent, as in ``test_gsea_tu_calibration.py``) is a confirmed fail.
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
from pathlib import Path
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
@pytest.mark.slow
def test_attack2_fitted_p_detect_confounded_category_at_level() -> None:
    """B1 gate (lead decision 4): with the shift carried by a covariate, the weighted ORA with
    the fitted M4-13 p_detect is at level on the confounded category."""
    _gate("ora", "attack.sim_function.ora_attack.fitted", _attack2_rep(True),
          C.replicates("ora", False))  # fmt: skip


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
def test_attack3_packed_sets_are_tested_after_tu_collapse() -> None:
    """B2: every packed set holds >= 10 scoring TUs, so the TU-collapsed engine tests all 20."""
    from giap.config.model import GseaConfig
    from giap.function import gsea
    from giap.rng import for_step

    w = _fs().draw(_sc().gsea_attack_operon_packed(4000, 20), np.random.default_rng(7))
    tab = gsea.gsea_table(_gsea_genes(w), [_coll("kegg", _members(w))], GSEA_SIZE, GseaConfig(),
                          0.05, for_step("function.gsea", 1))  # fmt: skip
    assert int(np.isfinite(tab["p_value"].to_numpy(float)).sum()) == 20, len(tab)


@pytest.mark.calsuite("gsea")
def test_attack3_packed_tus_exchangeable_with_other_tus() -> None:
    """Fligner-Killeen test of the true TU level ``tu_log2``: packed operons (TUs of >= 10 genes
    in a packed set) against the other chr1 TUs, one p per world."""
    from scipy import stats

    def fn(rng: np.random.Generator) -> float:
        w = _fs().draw(_sc().gsea_attack_operon_packed(4000, 20), rng)
        tu, _ = _tu_maps(w)
        tus = w.layout.tus
        size = {str(t): len(gs) for t, gs in zip(tus["tu_id"], tus["genes"], strict=True)}
        chr1 = set(tus.loc[tus["contig"] == "chr1", "tu_id"].astype(str))
        packed = {tu[g] for g in w.sets["locus_tag"].astype(str) if size[tu[g]] >= 10}
        g = w.expression.genes
        lvl = g.assign(tu=g["locus_tag"].astype(str).map(tu)).groupby("tu")["tu_log2"].first()
        inside = lvl.index.isin(packed)
        rest = lvl.index.isin(chr1) & ~inside
        assert inside.sum() >= 20
        return float(stats.fligner(lvl[inside].to_numpy(float),
                                   lvl[rest].to_numpy(float)).pvalue)  # fmt: skip

    _gate("gsea", "attack.sim_function.packed.exchangeable", fn, C.replicates("gsea", False))


def _packed_pick(rng: np.random.Generator) -> tuple[pd.DataFrame, set[str]]:
    w = _fs().draw(_sc().gsea_attack_operon_packed(4000, 20), rng)
    sets = _members(w)
    name = sorted(sets)[int(rng.integers(0, len(sets)))]
    return _gsea_genes(w), sets[name]


def _packed_tu_rep(rng: np.random.Generator) -> float:
    from giap.config.model import GseaConfig
    from giap.function import gsea
    from giap.rng import for_step

    genes, members = _packed_pick(rng)
    step = for_step("function.gsea", int(rng.integers(0, 2**31 - 1)))
    tab = gsea.gsea_table(genes, [_coll("kegg", {"S": members})], GSEA_SIZE, GseaConfig(), 0.05,
                          step)  # fmt: skip
    p = tab["p_value"].to_numpy(float)
    assert p.size == 1 and np.isfinite(p[0]), tab
    return float(p[0])


def _packed_gene_level_rep(rng: np.random.Generator) -> float:
    """The arm GIAP rejects: detected genes ranked individually, permuted as if independent."""
    from giap.config.model import GseaConfig
    from giap.rng import for_step
    from giap.stats.fgsea import fgsea_simple, prepare

    genes, members = _packed_pick(rng)
    d = genes[genes["detected"]]
    x = d["abundance_log2"].to_numpy(float)
    cfg = GseaConfig()
    prep = prepare(x - np.median(x), list(d["locus_tag"]), {"S": sorted(members)}, *GSEA_SIZE)
    step = for_step("function.gsea", int(rng.integers(0, 2**31 - 1)))
    res = fgsea_simple(prep, cfg.n_perm_simple, cfg.sample_size, step)
    return float(res.p[0])


@pytest.mark.calsuite("gsea")
def test_attack3_tu_level_gsea_at_level_on_packed_sets() -> None:
    """B2 gate: after the TU collapse a packed set (1-3 operons + single-gene TUs) is at level."""
    _gate("gsea", "attack.sim_function.packed.tu_level", _packed_tu_rep,
          C.replicates("gsea", False))  # fmt: skip


@pytest.mark.calsuite("gsea")
def test_teeth_attack3_gene_level_gsea_fails_on_packed_sets(teeth_log: Path) -> None:
    """B2 teeth: gene-label permutation sees 10-30 correlated genes per packed operon and is
    anti-conservative on the same sets: a confirmed fail at 0.05."""
    rec = C.calibrate("gsea", "level", 0.05, _packed_gene_level_rep, C.replicates("gsea", False),
                      stream="attack.sim_function.packed.teeth_gene_level",
                      log_dir=teeth_log)  # fmt: skip
    assert rec["status"] == "fail", f"the gene-level arm is not a confirmed fail: {rec}"


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
