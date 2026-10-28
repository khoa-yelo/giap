"""M6-02 detection-limit model: examples, degenerate cases and interface
(docs/specs/M6-02-detection-limit-model.md; GIPPA_SPEC.md §7.4.3 "One detected, one not", §7.3.3).

Interface pinned here (the spec's sketch plus what its acceptance criteria need to be testable):

- ``gippa.operons.detection_limit.PMissModel`` (frozen, slots): the six spec fields ``kind``,
  ``coef``, ``k_min``, ``n_fit``, ``auc``, ``usable`` plus, for cross-fitting (spec 2, "a partner's
  peptides never train the model that scores it"), ``fold_of_tu: Mapping[str, int]`` and
  ``fold_coef: tuple[Mapping[str, float], ...]`` (``fold_coef[j]`` trained on the TUs whose fold is
  not ``j``), both defaulting to empty so the toy model is built from the six spec fields.
- Peptide-model ``coef`` keys: ``intercept``, ``abundance_log2``, ``length``, ``gravy`` (linear
  logit, as in the spec's toy model ``logit p_j = -12 + 0.6 a - 0.02 L_j``).
- ``model_for_tu(model, tu_id) -> PMissModel``: the out-of-fold model of ``tu_id`` (coef =
  ``fold_coef[fold_of_tu[tu_id]]``); ``model`` itself for a TU that trained nothing.
- ``fit_detection_limit(genes, peptide_hits, peptide_runs, cfg, rng) -> DetectionLimitFit`` with
  ``model`` (the one in force), ``peptide``, ``replicate`` (each ``PMissModel | None``) and
  ``summary`` (detection_summary rows, section ``detection_limit``, one ``auc`` row per fitted
  model with ``group`` = its kind: lead decision 3, "record both AUCs").
- ``fit_replicate_model(genes, peptide_runs, cfg, rng, *, peptide_hits)``: C16 says which peptide
  was seen in which run, C11 says which gene each peptide is assigned to.
- ``p_detect`` / ``p_miss`` take a keyword-only ``cov`` (the gene's covariates, by name) for the
  replicate model; the peptide model ignores it.
- ``genes`` frame columns read: ``locus_tag``, ``tu_id``, ``aa``, ``detected``, ``abundance_log2``,
  ``detection_tier`` (+ ``cfg.detectability.covariates`` for the replicate model).
"""

from __future__ import annotations

import dataclasses
import math
import re
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import detection_limit_data as D

from gippa.config import GippaConfig
from gippa.rng import for_step

SEED = 20261021


def _dl() -> Any:
    from gippa.operons import detection_limit

    return detection_limit


def _toy(k_min: int = 1, **kw: Any) -> Any:
    return _dl().PMissModel(
        kind="peptide", coef=dict(D.TOY_COEF), k_min=k_min, n_fit=1000, auc=0.8, usable=True, **kw
    )


def _cfg(**kw: Any) -> GippaConfig:
    return GippaConfig(seed=SEED, **kw)


def _rng(seed: int = SEED) -> np.random.Generator:
    return for_step(D.STEP, seed)


BETA = {"intercept": -10.0, "abundance_log2": 0.5, "length": -0.03, "gravy": 0.8}


def _sim(n: int = 300, seed: int = 1, **kw: Any) -> D.PeptideSim:
    return D.simulate_peptides(np.random.default_rng(seed), n, BETA, **kw)


# ---------------------------------------------------------------------------- interface


def test_pmissmodel_is_frozen_with_spec_fields() -> None:
    m = _toy()
    names = {f.name for f in dataclasses.fields(m)}
    assert {"kind", "coef", "k_min", "n_fit", "auc", "usable"} <= names
    assert {"fold_of_tu", "fold_coef"} <= names
    with pytest.raises(dataclasses.FrozenInstanceError):
        m.k_min = 2  # type: ignore[misc]
    assert not hasattr(m, "__dict__")  # slots
    assert dict(m.fold_of_tu) == {} and tuple(m.fold_coef) == ()


def test_no_q_value_in_the_public_interface() -> None:
    """Lead decision 4: p_miss is a probability, never a p-value; nothing here makes a q."""
    pub = getattr(_dl(), "__all__", [n for n in dir(_dl()) if not n.startswith("_")])
    assert "p_miss" in pub and "p_detect" in pub
    bad = [n for n in pub if re.search(r"(^|_)(q|qvalue|q_value|qvalues|bh)(_|$)", n.lower())]
    assert bad == []


# ---------------------------------------------------------------------------- examples


def test_example1_per_peptide_probabilities() -> None:
    """Example 1: per-peptide p at a = 16 are 0.0691, 0.0573, 0.0474 (one-peptide p_detect)."""
    m = _toy()
    got = [_dl().p_detect(m, 16.0, [pep]) for pep in D.TOY_PEPS]
    assert got == pytest.approx([0.0691, 0.0573, 0.0474], abs=5e-5)
    oracle = [D.peptide_p(D.TOY_COEF, 16.0, n, g) for n, g in D.TOY_PEPS]
    assert got == pytest.approx(oracle, rel=1e-12)


def test_example2_p_detect_k1() -> None:
    p = _dl().p_detect(_toy(1), 17.0, list(D.TOY_PEPS))
    assert p == pytest.approx(0.2730, abs=5e-5)
    assert p == pytest.approx(D.p_detect_oracle(D.TOY_COEF, 17.0, D.TOY_PEPS, 1), rel=1e-12)


def test_example3_p_detect_k2() -> None:
    p = _dl().p_detect(_toy(2), 17.0, list(D.TOY_PEPS))
    assert p == pytest.approx(0.0281, abs=5e-5)
    assert p == pytest.approx(D.p_detect_oracle(D.TOY_COEF, 17.0, D.TOY_PEPS, 2), rel=1e-12)


def test_example4_p_miss() -> None:
    """Example 4: per-r 1 - P_det = 0.8359, 0.7524, 0.6853, 0.4820 -> mean 0.6889."""
    dl = _dl()
    ratios = [-1.0, -0.2, 0.3, 1.5]
    per_r = [1.0 - dl.p_detect(_toy(), 17.0 + r, list(D.TOY_PEPS)) for r in ratios]
    assert per_r == pytest.approx([0.8359, 0.7524, 0.6853, 0.4820], abs=5e-5)
    pm = dl.p_miss(_toy(), 17.0, ratios, list(D.TOY_PEPS))
    assert pm == pytest.approx(0.6889, abs=5e-5)
    assert pm == pytest.approx(D.p_miss_oracle(D.TOY_COEF, 17.0, ratios, D.TOY_PEPS, 1), rel=1e-12)


def test_example4_ratios_as_numpy_array_not_mutated() -> None:
    ratios = np.array([-1.0, -0.2, 0.3, 1.5])
    before = ratios.copy()
    pm = _dl().p_miss(_toy(), 17.0, ratios, list(D.TOY_PEPS))
    assert pm == pytest.approx(0.6889, abs=5e-5)
    np.testing.assert_array_equal(ratios, before)


def test_example5_partner_expected_far_higher() -> None:
    """Example 5: ratios = [3.0] -> small p_miss (< 0.25); hand value 1 - P_det(20) = 0.2125."""
    pm = _dl().p_miss(_toy(), 17.0, [3.0], list(D.TOY_PEPS))
    assert pm < 0.25
    assert pm == pytest.approx(0.21254, abs=5e-5)


def test_example6_engine_check() -> None:
    from gippa.stats.poisson_binomial import poisson_binomial_sf

    assert poisson_binomial_sf([0.3, 0.4, 0.6], 1) == pytest.approx(0.8320, abs=1e-12)
    assert poisson_binomial_sf([0.3, 0.4, 0.6], 2) == pytest.approx(0.3960, abs=1e-12)
    # and p_detect is that engine applied to the per-peptide probabilities
    m = _dl().PMissModel(
        kind="peptide",
        coef={"intercept": 0.0, "abundance_log2": 0.0, "length": 0.0, "gravy": 1.0},
        k_min=2,
        n_fit=100,
        auc=0.7,
        usable=True,
    )
    g = [math.log(q / (1 - q)) for q in (0.3, 0.4, 0.6)]  # gravy coefficient 1: logit = gravy
    assert _dl().p_detect(m, 5.0, [(10, x) for x in g]) == pytest.approx(0.3960, abs=1e-12)


def test_example8_no_observable_peptide_is_nan() -> None:
    dl = _dl()
    assert math.isnan(dl.p_detect(_toy(), 17.0, []))
    assert math.isnan(dl.p_miss(_toy(), 17.0, [-1.0, 0.5], []))


def test_example8_gene_without_observable_peptides_has_none() -> None:
    """A protein whose tryptic fragments are all < 7 or > 30 aa has no observable peptide."""
    cfg = _cfg()
    assert _dl().observable_peptides("MKAGKRSSK", cfg) == []
    assert _dl().observable_peptides("A" * 40 + "K", cfg) == []


def test_example9_same_seed_identical_folds_and_coefficients() -> None:
    s = _sim(200, seed=9)
    dl = _dl()
    a = dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng(7))
    b = dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng(7))
    assert a.usable and b.usable
    assert dict(a.fold_of_tu) == dict(b.fold_of_tu)
    assert dict(a.coef) == dict(b.coef)
    assert [dict(c) for c in a.fold_coef] == [dict(c) for c in b.fold_coef]
    assert a.auc == b.auc and a.n_fit == b.n_fit
    c = dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng(8))
    assert dict(c.fold_of_tu) != dict(a.fold_of_tu)  # 100 TUs: the folds follow the seed


@pytest.mark.slow
def test_example7_known_beta_recovered_within_10pct() -> None:
    """Example 7: ~5,000 training peptides from a known logistic truth; every coefficient within
    10 % and AUC >= 0.75. (The spec names M6-03's evidence simulator, not merged yet; the truth
    here is the same logistic in (a, L, GRAVY) on hand-built tryptic proteins. Design chosen so the
    unpenalised MLE is within 5 % on this seed; nearly every gene is detected, so selection on
    detection is negligible.)"""
    beta = {"intercept": -8.0, "abundance_log2": 0.8, "length": -0.15, "gravy": 1.0}
    s = D.simulate_peptides(
        for_step("test.ex7", 0), 0, beta, a_mu=14.0, a_sd=3.0, n_peps=(20, 30), min_total=5000
    )
    m = _dl().fit_peptide_model(s.genes, s.hits, _cfg(), _rng())
    assert m.usable and m.kind == "peptide"
    n_rows = sum(len(s.peptides[t]) for t in s.genes.loc[s.genes.detected, "locus_tag"])
    assert m.n_fit == n_rows >= 5000
    for k, v in beta.items():
        assert m.coef[k] == pytest.approx(v, rel=0.10), k
    assert m.auc >= 0.75


# ---------------------------------------------------------------------------- observable peptides


def test_observable_peptides_lengths_and_gravy() -> None:
    peps = ["ACDEFGHK", "MNQSTVWYK", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGK"]  # 8, 9, 30 aa
    aa = "MK" + "".join(peps) + "AK" + "A" * 30 + "K"  # "MK", "AK" too short; 31 aa too long
    got = sorted(_dl().observable_peptides(aa, _cfg()))
    assert [n for n, _ in got] == [8, 9, 30]
    want = sorted((len(p), D.gravy(p)) for p in peps)
    assert [g for _, g in got] == pytest.approx([g for _, g in want], abs=1e-9)
    assert all(isinstance(n, int) for n, _ in got)


def test_observable_peptides_dedup_by_sequence_not_by_tuple() -> None:
    """Deduplicated within the gene by peptide sequence: a repeated peptide counts once; two
    different peptides with the same (length, GRAVY) both count."""
    p1, p2 = "ADEGHSTK", "DAEGHSTK"  # anagrams: same length and GRAVY
    got = _dl().observable_peptides(p1 + p2 + p1 + p1, _cfg())
    assert len(got) == 2
    assert got[0] == pytest.approx(got[1])


def test_observable_peptides_uses_configured_range() -> None:
    cfg = _cfg(abundance={"observable_len": [9, 30]})
    got = _dl().observable_peptides("ACDEFGHK" + "MNQSTVWYK", cfg)
    assert [n for n, _ in got] == [9]


# ---------------------------------------------------------------------------- the peptide model


def _small_world() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    """Hand-built: 3 detected genes (5, 4, 3 observable peptides), 1 undetected (4), 1 detected
    with NaN abundance (2). Spacers 'GK' are never observable."""
    rng = np.random.default_rng(5)
    peps = {t: D.make_peptides(rng, n) for t, n in zip("ABCDE", (5, 4, 3, 4, 2), strict=True)}
    genes = D.genes_frame(
        list("ABCDE"),
        ["T1", "T1", "T2", "T3", "T4"],
        [D.protein(peps[t], spacers=True) for t in "ABCDE"],
        [True, True, True, False, True],
        [16.0, 15.0, 14.0, 12.0, math.nan],
    )
    hits = D.hits_frame(
        [D.hit_row(peps["A"][0], "A"), D.hit_row(peps["A"][1], "A"), D.hit_row(peps["B"][0], "B"),
         D.hit_row(peps["C"][2], "C"), D.hit_row(peps["E"][0], "E")]
    )  # fmt: skip
    return genes, hits, peps


def test_training_rows_are_observable_peptides_of_detected_genes_with_abundance() -> None:
    """Spec 2: 5 + 4 + 3 = 12 rows (D undetected, E NaN abundance); 12 < 50 -> not usable."""
    genes, hits, _ = _small_world()
    m = _dl().fit_peptide_model(genes, hits, _cfg(), _rng())
    assert m.n_fit == 12
    assert m.usable is False
    assert m.kind == "peptide"


def test_fewer_than_50_rows_not_usable_boundary() -> None:
    """49 training rows (7 detected genes x 7 peptides, both classes present) -> usable False."""
    rng = np.random.default_rng(11)
    tags = [f"g{i}" for i in range(7)]
    peps = {t: D.make_peptides(rng, 7) for t in tags}
    genes = D.genes_frame(
        tags, [f"T{i}" for i in range(7)], [D.protein(peps[t]) for t in tags], [True] * 7,
        [12.0 + i for i in range(7)],
    )  # fmt: skip
    hits = D.hits_frame([D.hit_row(p, t) for t in tags for p in peps[t][: 1 + int(t[1:]) % 6]])
    m = _dl().fit_peptide_model(genes, hits, _cfg(), _rng())
    assert m.n_fit == 49 and m.usable is False


def test_usable_on_a_realistic_world() -> None:
    s = _sim(300, seed=2)
    m = _dl().fit_peptide_model(s.genes, s.hits, _cfg(), _rng())
    assert m.usable is True and m.kind == "peptide" and m.k_min == 1
    assert m.n_fit == sum(len(s.peptides[t]) for t in s.genes.loc[s.genes.detected, "locus_tag"])
    assert set(m.coef) == {"intercept", "abundance_log2", "length", "gravy"}
    assert all(math.isfinite(v) for v in m.coef.values())
    assert m.coef["abundance_log2"] > 0  # truth 0.5
    assert 0.5 < m.auc <= 1.0
    k = _cfg().detectability.folds
    assert len(m.fold_coef) == k
    assert set(m.fold_of_tu.values()) <= set(range(k))
    train_tus = set(s.genes.loc[s.genes.detected, "tu_id"])
    assert train_tus <= set(m.fold_of_tu)


def test_one_class_response_is_constant_not_usable() -> None:
    """No observed peptide at all (empty C11): the logistic has one class -> usable False."""
    s = _sim(200, seed=3)
    m = _dl().fit_peptide_model(s.genes, D.hits_frame([]), _cfg(), _rng())
    assert m.usable is False


def test_every_peptide_observed_is_constant_not_usable() -> None:
    s = _sim(150, seed=4)
    rows = [D.hit_row(p, t) for t in s.genes.locus_tag for p in s.peptides[t]]
    genes = s.genes.assign(detected=True)
    m = _dl().fit_peptide_model(genes, D.hits_frame(rows), _cfg(), _rng())
    assert m.usable is False


def test_separated_logistic_never_returns_nonfinite() -> None:
    """Observed iff a > 15 (perfect separation on abundance): M4-13 handling; either a finite,
    usable model with probabilities in [0, 1], or usable False."""
    s = _sim(200, seed=12)
    rows = [
        D.hit_row(p, t)
        for t, a in zip(s.genes.locus_tag, s.genes.abundance_log2, strict=True)
        if a > 15
        for p in s.peptides[t]
    ]
    genes = s.genes.assign(detected=True)
    m = _dl().fit_peptide_model(genes, D.hits_frame(rows), _cfg(), _rng())
    if m.usable:
        assert all(math.isfinite(v) for v in m.coef.values())
        for a in (10.0, 15.0, 20.0):
            p = _dl().p_detect(m, a, s.peps("G00000"))
            assert 0.0 <= p <= 1.0


def test_observed_means_any_assignment_uniqueness_ignored() -> None:
    """Lead decision 2 (Q2): 'observed' = any C11 assignment of the peptide to the gene;
    ``shared`` / ``multi_identical`` count exactly like ``unique``."""
    s = _sim(250, seed=13)
    dl = _dl()
    base = dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng())
    for u in ("shared", "multi_identical"):
        h = s.hits.assign(uniqueness=u)
        m = dl.fit_peptide_model(s.genes, h, _cfg(), _rng())
        assert dict(m.coef) == pytest.approx(dict(base.coef), rel=1e-9, abs=1e-12)
        assert m.n_fit == base.n_fit


def test_paralogues_both_counted() -> None:
    """Degenerate row 'every peptide of a gene is also in another gene (shared)': still counted.
    A duplicated gene G'' with the same sequence and the same peptides assigned (as ``shared``)
    adds its observable peptides to the training rows."""
    s = _sim(250, seed=14)
    dl = _dl()
    base = dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng())
    src = s.genes[s.genes.detected].iloc[0]
    dup = src.copy()
    dup["locus_tag"], dup["tu_id"] = "PARALOG", "TU_PARALOG"
    genes = pd.concat([s.genes, dup.to_frame().T], ignore_index=True)
    genes = genes.astype({"detected": bool, "abundance_log2": "float64"})
    own = s.hits[s.hits.gene == src.locus_tag]
    hits = pd.concat(
        [s.hits.assign(uniqueness=np.where(s.hits.gene == src.locus_tag, "shared", "unique")),
         own.assign(gene="PARALOG", uniqueness="shared")],
        ignore_index=True,
    )  # fmt: skip
    m = dl.fit_peptide_model(genes, hits, _cfg(), _rng())
    assert m.n_fit == base.n_fit + len(s.peptides[src.locus_tag])


def test_unassigned_hits_and_unobservable_peptides_do_not_count() -> None:
    """A C11 row with ``gene`` NA does not make the peptide observed for any gene; a hit on a
    peptide outside the gene's observable set adds no training row."""
    s = _sim(250, seed=15)
    dl = _dl()
    base = dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng())
    # every unobserved peptide gets an unassigned (gene NA) hit; plus a hit on a 6-mer
    extra = [
        D.hit_row(p, None)
        for t in s.genes.locus_tag
        for p, o in zip(s.peptides[t], s.observed[t], strict=True)
        if not o
    ]
    extra.append(D.hit_row("ADEGHK", str(s.genes.locus_tag.iloc[0])))
    hits = pd.concat([s.hits, D.hits_frame(extra)], ignore_index=True)
    m = dl.fit_peptide_model(s.genes, hits, _cfg(), _rng())
    assert m.n_fit == base.n_fit
    assert dict(m.coef) == pytest.approx(dict(base.coef), rel=1e-9, abs=1e-12)


def test_two_peptide_tier_sets_k_min_2() -> None:
    """Spec 3: k_min = 2 under the two_peptide tier (detection.min_peptides_sequences_only)."""
    s = _sim(300, seed=16, tier="two_peptide", k_min=2)
    m = _dl().fit_peptide_model(s.genes, s.hits, _cfg(), _rng())
    assert m.k_min == _cfg().detection.min_peptides_sequences_only == 2
    for tier in ("mayu", "decoy_picked"):
        g = s.genes.assign(detection_tier=tier)
        assert _dl().fit_peptide_model(g, s.hits, _cfg(), _rng()).k_min == 1


def test_k_min_2_with_one_observable_peptide() -> None:
    """Degenerate row: k_min = 2 with one observable peptide -> p_detect 0, p_miss 1."""
    dl = _dl()
    assert dl.p_detect(_toy(2), 25.0, [(10, 0.0)]) == 0.0
    assert dl.p_miss(_toy(2), 25.0, [-0.5, 0.0, 2.0], [(10, 0.0)]) == 1.0


def test_empty_ratios_is_nan() -> None:
    assert math.isnan(_dl().p_miss(_toy(), 17.0, [], list(D.TOY_PEPS)))
    assert math.isnan(_dl().p_miss(_toy(), 17.0, np.array([]), list(D.TOY_PEPS)))


def test_unusable_model_gives_nan() -> None:
    """No silent fallback: a non-usable model scores nothing (M6-05 then emits
    ``no_detection_limit_model``)."""
    m = dataclasses.replace(_toy(), usable=False)
    assert math.isnan(_dl().p_detect(m, 17.0, list(D.TOY_PEPS)))
    assert math.isnan(_dl().p_miss(m, 17.0, [0.0], list(D.TOY_PEPS)))


def test_extrapolation_far_outside_training_range_is_a_probability() -> None:
    """Degenerate row: a_v + r far outside the training range -> the logistic extrapolates."""
    dl = _dl()
    lo = dl.p_detect(_toy(), -50.0, list(D.TOY_PEPS))
    hi = dl.p_detect(_toy(), 80.0, list(D.TOY_PEPS))
    assert 0.0 <= lo < 1e-6 and 1.0 - 1e-6 < hi <= 1.0


# ---------------------------------------------------------------------------- cross-fitting


def test_model_for_tu_uses_the_out_of_fold_coefficients() -> None:
    s = _sim(250, seed=17)
    dl = _dl()
    m = dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng())
    for tu, j in list(m.fold_of_tu.items())[:10]:
        mt = dl.model_for_tu(m, tu)
        assert dict(mt.coef) == dict(m.fold_coef[j])
        assert mt.k_min == m.k_min and mt.usable == m.usable and mt.kind == m.kind
    assert dl.model_for_tu(m, "TU_NOT_TRAINED") == m


@pytest.mark.parametrize("j", [0, 2, 4])
def test_cross_fit_fold_model_ignores_its_own_tus(j: int) -> None:
    """Acceptance: no gene is scored by a model trained on its own TU. Rewriting the C11 evidence
    of every gene whose TU is in fold j leaves ``fold_coef[j]`` bit-identical (same seed, same
    folds) and changes the other folds' models."""
    s = _sim(250, seed=18)
    dl = _dl()
    a = dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng(3))
    in_j = {tu for tu, f in a.fold_of_tu.items() if f == j}
    tags_j = set(s.genes.loc[s.genes.tu_id.isin(in_j), "locus_tag"])
    assert tags_j
    keep = s.hits[~s.hits.gene.isin(tags_j)]
    flip = [D.hit_row(p, t) for t in sorted(tags_j) for p in s.peptides[t]]  # all observed
    b = dl.fit_peptide_model(
        s.genes, pd.concat([keep, D.hits_frame(flip)], ignore_index=True), _cfg(), _rng(3)
    )
    assert dict(b.fold_of_tu) == dict(a.fold_of_tu)
    assert dict(b.fold_coef[j]) == dict(a.fold_coef[j])
    others = [i for i in range(len(a.fold_coef)) if i != j]
    assert any(dict(b.fold_coef[i]) != dict(a.fold_coef[i]) for i in others)


def test_folds_from_the_step_stream() -> None:
    """Spec 7: the caller's ``for_step("operons.detection_limit", seed)`` generator decides the
    folds; the TU-fold map is a pure function of the seed (not of data order)."""
    s = _sim(200, seed=19)
    dl = _dl()
    a = dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng(5))
    perm = s.genes.sample(frac=1.0, random_state=0).reset_index(drop=True)
    b = dl.fit_peptide_model(perm, s.hits.iloc[::-1].reset_index(drop=True), _cfg(), _rng(5))
    assert dict(a.fold_of_tu) == dict(b.fold_of_tu)
    assert dict(b.coef) == pytest.approx(dict(a.coef), rel=1e-3, abs=1e-6)


# ---------------------------------------------------------------------------- inputs untouched


def test_inputs_not_mutated() -> None:
    s = _sim(200, seed=20)
    genes, hits = s.genes.copy(deep=True), s.hits.copy(deep=True)
    runs = _runs_for(s, 3, seed=20)
    runs0 = runs.copy(deep=True)
    dl = _dl()
    dl.fit_peptide_model(s.genes, s.hits, _cfg(), _rng())
    dl.fit_detection_limit(s.genes, s.hits, runs, _cfg(), _rng())
    dl.fit_replicate_model(s.genes, runs, _cfg(), _rng(), peptide_hits=s.hits)
    pd.testing.assert_frame_equal(s.genes, genes)
    pd.testing.assert_frame_equal(s.hits, hits)
    pd.testing.assert_frame_equal(runs, runs0)
    peps = list(D.TOY_PEPS)
    dl.p_miss(_toy(), 17.0, [0.1], peps)
    assert peps == list(D.TOY_PEPS)


# ---------------------------------------------------------------------------- dispatcher, C16


def _runs_for(s: D.PeptideSim, n_runs: int, seed: int) -> pd.DataFrame:
    """C16: each observed peptide is seen in each run with probability 0.6 (at least one run)."""
    rng = np.random.default_rng(seed)
    rows: list[tuple[str, str]] = []
    for t in s.genes.locus_tag:
        for p, o in zip(s.peptides[t], s.observed[t], strict=True):
            if not o:
                continue
            seen = rng.random(n_runs) < 0.6
            seen[int(rng.integers(n_runs))] = True
            rows.extend((p, f"run{r}") for r in range(n_runs) if seen[r])
    return D.runs_frame(rows)


def test_no_peptides_and_one_run_no_fit() -> None:
    """Degenerate row 1: no peptide input and one run -> usable False, no fit attempted."""
    s = _sim(200, seed=21)
    fit = _dl().fit_detection_limit(s.genes, None, None, _cfg(), _rng())
    assert fit.model.usable is False
    one = D.runs_frame([(p, "only") for p in s.hits.peptide])
    fit1 = _dl().fit_detection_limit(s.genes, D.hits_frame([]), one, _cfg(), _rng())
    assert fit1.model.usable is False
    r = _dl().fit_replicate_model(s.genes, one, _cfg(), _rng(), peptide_hits=s.hits)
    assert r.usable is False and r.n_fit == 0


def test_c16_absent_peptide_model_still_fits() -> None:
    """C16 is optional (lead decision 1): ``peptide_runs = None`` -> the peptide model alone."""
    s = _sim(300, seed=22)
    fit = _dl().fit_detection_limit(s.genes, s.hits, None, _cfg(), _rng())
    assert fit.model.kind == "peptide" and fit.model.usable
    assert fit.replicate is None
    direct = _dl().fit_peptide_model(s.genes, s.hits, _cfg(), _rng())
    assert dict(fit.model.coef) == dict(direct.coef)


def test_peptides_and_runs_peptide_model_wins_both_aucs_recorded() -> None:
    """Degenerate row 2 + lead decision 3: peptides and >= 2 runs -> the peptide model is in
    force; both AUCs are recorded in the detection summary."""
    s = _sim(400, seed=23)
    runs = _runs_for(s, 4, seed=23)
    fit = _dl().fit_detection_limit(s.genes, s.hits, runs, _cfg(), _rng())
    assert fit.model.kind == "peptide" and fit.model.usable
    assert fit.peptide is not None and fit.replicate is not None
    assert fit.model == fit.peptide
    summ = fit.summary
    assert set(summ["section"]) == {"detection_limit"}
    aucs = summ[summ["metric"] == "auc"].set_index("group")["value"]
    assert set(aucs.index) == {"peptide", "replicate"}
    assert aucs["peptide"] == pytest.approx(fit.peptide.auc, nan_ok=True)
    assert aucs["replicate"] == pytest.approx(fit.replicate.auc, nan_ok=True)


def test_replicate_model_training_rows_and_direction() -> None:
    """Spec 4: rows = (gene, run) for genes seen in some runs but not all; regressed on the gene's
    abundance_log2 (covariates off here). Truth: per-run detection logit = -8 + 0.5 a."""
    rng = np.random.default_rng(24)
    n, n_runs = 600, 4
    tags = [f"R{i:04d}" for i in range(n)]
    a = rng.normal(16.0, 2.0, n)
    peps = {t: D.make_peptides(rng, 3) for t in tags}
    p = 1.0 / (1.0 + np.exp(-(-8.0 + 0.5 * a)))
    seen = rng.random((n, n_runs)) < p[:, None]
    det = seen.any(axis=1)
    genes = D.genes_frame(
        tags, [f"T{i}" for i in range(n)], [D.protein(peps[t]) for t in tags], det.tolist(), a
    )
    hits = D.hits_frame([D.hit_row(peps[t][0], t) for t, d in zip(tags, det, strict=True) if d])
    runs = D.runs_frame(
        [(peps[t][0], f"run{r}") for i, t in enumerate(tags) for r in range(n_runs) if seen[i, r]]
    )
    cfg = _cfg(detectability={"covariates": []})
    m = _dl().fit_replicate_model(genes, runs, cfg, _rng(), peptide_hits=hits)
    mixed = int(((seen.sum(1) > 0) & (seen.sum(1) < n_runs)).sum())
    assert m.kind == "replicate" and m.usable
    assert m.n_fit == mixed * n_runs
    assert set(m.coef) == {"intercept", "abundance_log2"}
    assert m.coef["abundance_log2"] > 0
    # replicate p_detect is the model's probability directly (no Poisson-binomial, peps ignored)
    want = 1.0 / (1.0 + math.exp(-(m.coef["intercept"] + m.coef["abundance_log2"] * 16.0)))
    dl = _dl()
    assert dl.p_detect(m, 16.0, [(10, 0.0)]) == pytest.approx(want, rel=1e-12)
    assert dl.p_detect(m, 16.0, list(D.TOY_PEPS)) == pytest.approx(want, rel=1e-12)
    assert dl.p_detect(m, 18.0, []) >= dl.p_detect(m, 16.0, [])


def test_replicate_model_with_default_covariates_takes_cov() -> None:
    rng = np.random.default_rng(25)
    n, n_runs = 400, 4
    tags = [f"R{i:04d}" for i in range(n)]
    a = rng.normal(16.0, 2.0, n)
    peps = {t: D.make_peptides(rng, 3) for t in tags}
    seen = rng.random((n, n_runs)) < (1.0 / (1.0 + np.exp(-(-8.0 + 0.5 * a))))[:, None]
    det = seen.any(axis=1)
    cov = {
        "n_observable_peptides": [3] * n,
        "length_aa": [len(D.protein(peps[t])) for t in tags],
        "gravy": rng.normal(-0.3, 0.3, n).tolist(),
        "pi": rng.normal(7.0, 1.5, n).tolist(),
        "tm_segments": rng.integers(0, 3, n).tolist(),
    }
    genes = D.genes_frame(
        tags, [f"T{i}" for i in range(n)], [D.protein(peps[t]) for t in tags], det.tolist(), a,
        extra=cov,
    )  # fmt: skip
    hits = D.hits_frame([D.hit_row(peps[t][0], t) for t, d in zip(tags, det, strict=True) if d])
    runs = D.runs_frame(
        [(peps[t][0], f"run{r}") for i, t in enumerate(tags) for r in range(n_runs) if seen[i, r]]
    )
    m = _dl().fit_replicate_model(genes, runs, _cfg(), _rng(), peptide_hits=hits)
    assert m.usable and {"intercept", "abundance_log2"} <= set(m.coef)
    g0 = {k: v[0] for k, v in cov.items()}
    p = _dl().p_detect(m, 16.0, [], cov=g0)
    assert 0.0 <= p <= 1.0
    assert _dl().p_miss(m, 16.0, [0.0, 1.0], [], cov=g0) == pytest.approx(
        np.mean([1 - _dl().p_detect(m, 16.0 + r, [], cov=g0) for r in (0.0, 1.0)]), rel=1e-12
    )
