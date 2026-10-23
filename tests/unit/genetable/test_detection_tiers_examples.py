"""M4-12 detection tiers: examples 1-10 and edge cases (docs/specs/M4-12-detection-tiers.md).

Spec: GIAP_SPEC.md §2.1 decision 4, §7.3 (eligible = non-partial), §7.3.2 (tier table, picked steps
1-4, MAYU), §9.3 (``detection_status`` enum, ``detection_summary.tsv`` section ``tier``); ADR-0009
(picked decoy wins rescaled by 1/f); M4-index lead decisions 4 (decoys pass q <= ``q_value_max``;
the pairing cut-off reuses ``mapping.proteins.min_identity`` / ``min_coverage``, 0.80); M4-02
statistician items forwarded to M4-12 (input validation; ``fp`` and ``n_genes`` next to each q
threshold; unspecific-search warning; ``fdr_2pep`` never an FDR of the >= 2-peptide set); M4-11 lead
decision 3 (M4-12 decides whether contaminant evidence counts).

Inputs come from the merged producers (tests/helpers/detection_tiers_fixtures.py). Target and
decoy keys in examples 1, 2 and 10 are encoded with equal q and PEP and the score as tie-breaker,
so the key order is the order of the M4-01 example scores.

Resolved ambiguities (pinned here):
- ADR-0009 rescaling: FDR(t) = (D(t)/f + 1) / T(t), f = paired / (paired + unpaired) decoys among
  decoys with ``protein_ids`` passing q <= ``q_value_max`` ("each decoy win counts as 1/f"; the +1
  is not a decoy win). Summary values report f and 1/f.
- Contaminant evidence counts for every tier: accepted is accepted (§7.3.1), and the
  ``two_peptide`` tier must agree with M4-11 ``detected_2pep``, which counts contaminant peptides;
  contaminant prefixes only keep protein IDs out of the protein -> gene map.
- Decoys failing q <= ``q_value_max`` are neither paired nor unpaired (not in any count).
- Protein-group IDs ``A;B`` inside one ``protein_ids`` element are split; each member is stripped
  and looked up; a decoy pairs only when every member maps to the same one gene.
- Unspecific search (statistician, PR #34): in the ``mayu`` tier, when most accepted target
  peptides are not tryptic (C-terminal residue not K/R), ``warnings`` carries one entry containing
  "unspecific"; never for tryptic peptides.
- Unmapped proteins may be absent from the map or present with ``gene`` NA. The map covers every C2
  ``locus_tag`` (method ``locus_tag``), so a decoy ``rev_<tag>`` pairs even with a gene that has no
  target peptide (decoy-only genes enter the picked list as decoy wins, M4-01).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_fixtures as fx
from detection_tiers_fixtures import NAN, Builder, cfg, tag


def _det() -> Any:
    from giap.genetable import detection

    return detection


def _alpha(a: float, **kw: Any) -> Any:
    return cfg(stats={"alpha": {"identification": a}}, **kw)


def assert_table_shape(r: fx.Run) -> None:
    t = r.det.table
    assert set(fx.TABLE_COLUMNS) <= set(t.columns), list(t.columns)
    assert t["locus_tag"].astype(str).tolist() == r.genes["locus_tag"].astype(str).tolist()
    assert set(t["detection_status"].astype(str)) <= set(fx.STATUSES)
    assert (t["detection_tier"].astype(str) == r.det.tier).all()
    assert r.det.tier in fx.TIERS
    det = t["detected"].astype(bool).to_numpy()
    assert (det == (t["detection_status"].astype(str) == "detected").to_numpy()).all()


def _picked(scores_t: list[float], scores_d: list[float]) -> Builder:
    """Genes G1..Gn with one target (score 10 s) and, when not NaN, one paired decoy (10 d)."""
    b = Builder(len(scores_t))
    assert b.tags is not None
    for g, s, d in zip(b.tags, scores_t, scores_d, strict=True):
        if not math.isnan(s):
            b.target(g, 0.001, 0.01, 10.0 * s)
        if not math.isnan(d):
            b.decoy([f"rev_{g}"], 0.001, 0.01, 10.0 * d)
    return b


# ------------------------------------------------------------------------- example 1 (M4-01 ex 8)

EX1_T = [5, 4, 3, 3, 2, 1]
EX1_D = [1, NAN, 3.6, NAN, 2, NAN]
EX1_Q = [0.5, 0.5, NAN, 2 / 3, NAN, 0.75]


def test_example1_picked_gene_q_values() -> None:
    r = _picked(EX1_T, EX1_D).run()
    assert_table_shape(r)
    assert r.det.tier == "decoy_picked"
    for i, exp in enumerate(EX1_Q):
        g = tag(i)
        if math.isnan(exp):
            assert math.isnan(r.q(g)), (g, r.q(g))
        else:
            assert r.q(g) == pytest.approx(exp, rel=1e-12), g
        assert r.status(g) == "not_detected", g  # every q > 0.01
        assert not r.detected(g)


@pytest.mark.parametrize(
    ("alpha", "detected"),
    [(0.6, {0, 1}), (0.7, {0, 1, 3}), (0.75, {0, 1, 3, 5})],
)
def test_example1_detected_at_identification_alpha(alpha: float, detected: set[int]) -> None:
    r = _picked(EX1_T, EX1_D).run(_alpha(alpha))
    for i in range(6):
        assert r.detected(tag(i)) == (i in detected), (alpha, i)
        assert r.status(tag(i)) == ("detected" if i in detected else "not_detected")
    # decoy winners keep q NA whatever alpha (G3, G5)
    assert math.isnan(r.q(tag(2))) and math.isnan(r.q(tag(4)))


# ------------------------------------------------------------------------- example 2 (B.6 shape)


def test_example2_150_genes_above_one_decoy() -> None:
    b = Builder(150)
    assert b.tags is not None
    for i, g in enumerate(b.tags):
        b.target(g, 0.001, 0.01, 100.0 + i)
    b.decoy([f"rev_{b.tags[0]}"], 0.001, 0.01, 1.0)  # below every target
    r = b.run()
    assert r.det.tier == "decoy_picked"
    q = np.array([r.q(g) for g in b.tags])
    np.testing.assert_allclose(q, 1 / 150, rtol=1e-12)
    assert all(r.status(g) == "detected" for g in b.tags)
    assert fx.tier_metric(r.det, "n_detected") == 150
    assert fx.curve_metric(r.det, 0.01, "n_detected") == 150
    assert fx.curve_metric(r.det, 0.005, "n_detected") == 0  # 1/150 > 0.005


# ------------------------------------------------------------------------- examples 3-5 (pairing)


def _with_decoys(n_paired: int, n_unpaired: int, unpaired_q: float = 0.001) -> Builder:
    b = Builder(12)
    assert b.tags is not None
    for i, g in enumerate(b.tags):
        b.target(g, 0.001, 0.01, 100.0 + i)
    for i in range(n_paired):
        b.decoy([f"rev_{b.tags[i]}"], 0.001, 0.01, 1.0 + i)
    for i in range(n_unpaired):
        b.decoy([f"rev_UNKNOWN{i}"], unpaired_q, 0.01, 1.0 + i)
    return b


def test_example3_nine_of_ten_paired_is_decoy_picked() -> None:
    b = _with_decoys(9, 1)
    det = _det()
    pmap = det.protein_gene_map(b.genes_frame(), b.peptides_frame(), fx.eo.hits_frame(b.hits),
                                fx.ao.protein_hits_frame([]))  # fmt: skip
    pairing = det.pair_decoys(b.peptides_frame(), pmap, cfg().input.decoy_prefixes)
    assert {"peptide", "q_value", "pep", "score", "gene"} <= set(pairing.columns)
    assert len(pairing) == 10, "one row per decoy peptide"
    assert int(pairing["gene"].isna().sum()) == 1
    tier, warns = det.choose_tier(b.peptides_frame(), fx.eo.hits_frame(b.hits),
                                  fx.ao.protein_hits_frame([]), pairing)  # fmt: skip
    assert tier == "decoy_picked"
    assert "decoy_pairing_low" not in warns
    r = b.run()
    assert r.det.tier == "decoy_picked"
    assert fx.tier_metric(r.det, "n_decoys_paired") == 9
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 1


def test_example4_seven_of_ten_paired_falls_back_to_mayu() -> None:
    b = _with_decoys(7, 3)
    r = b.run()
    assert r.det.tier == "mayu"
    assert "decoy_pairing_low" in r.det.warnings
    assert set(r.det.table["detection_tier"].astype(str)) == {"mayu"}


def test_exactly_eighty_percent_paired_is_decoy_picked() -> None:
    r = _with_decoys(8, 2).run()
    assert r.det.tier == "decoy_picked"
    assert "decoy_pairing_low" not in r.det.warnings


def test_pairing_cutoff_reuses_min_identity_and_coverage() -> None:
    """Lead decision 4: the cut-off is mapping.proteins.min_identity / min_coverage (0.80)."""
    c = cfg(mapping={"proteins": {"min_identity": 0.95, "min_coverage": 0.95}})
    r = _with_decoys(9, 1).run(c)
    assert r.det.tier == "mayu"
    assert "decoy_pairing_low" in r.det.warnings


def test_decoys_above_q_value_max_are_not_counted() -> None:
    """Lead decision 4 (M4-12 Q1): decoys pass the same q <= q_value_max as targets."""
    r = _with_decoys(7, 3, unpaired_q=0.5).run()
    assert r.det.tier == "decoy_picked"
    assert fx.tier_metric(r.det, "n_decoys_paired") == 7
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0


def test_example5_decoy_on_two_genes_is_unpaired() -> None:
    b = Builder(2)
    g, h = b.tags or []
    b.protein_hit("P1", g)
    b.protein_hit("P2", h)
    b.protein_hit("P3", g)
    b.decoy(["rev_P1", "rev_P2"])
    b.decoy(["rev_P1;rev_P3"])  # one group element, both members on G -> paired with G
    det = _det()
    peps = b.peptides_frame()
    pmap = det.protein_gene_map(b.genes_frame(), peps, fx.eo.hits_frame([]),
                                fx.ao.protein_hits_frame(b.protein_hit_rows))  # fmt: skip
    assert fx.mapped(pmap, "P1") == (g, "protein_hit")
    assert fx.mapped(pmap, "P2") == (h, "protein_hit")
    pairing = det.pair_decoys(peps, pmap, cfg().input.decoy_prefixes).set_index("peptide")
    p_two, p_group = peps["peptide"].tolist()
    assert fx.is_na(pairing.loc[p_two, "gene"])
    assert pairing.loc[p_group, "gene"] == g


def test_first_matching_decoy_prefix_is_stripped() -> None:
    b = Builder(3)
    t = b.tags or []
    b.decoy([f"DECOY_{t[0]}"])
    b.decoy([f"##{t[1]}"])
    b.decoy([f"XXX_{t[2]}"])
    det = _det()
    peps = b.peptides_frame()
    pmap = det.protein_gene_map(b.genes_frame(), peps, fx.eo.hits_frame([]),
                                fx.ao.protein_hits_frame([]))  # fmt: skip
    pairing = det.pair_decoys(peps, pmap, cfg().input.decoy_prefixes)
    assert pairing.set_index("peptide").loc[peps["peptide"], "gene"].tolist() == t


# ------------------------------------------------------------------------- ADR-0009 rescaling


def test_adr0009_decoy_wins_rescaled_by_paired_share() -> None:
    """f = 8/10: T0..T4 (scores 100..96) above a decoy-only gene X (60), T5..T9 (50..46) below;
    seven more paired decoy-only genes below every target; two unpaired decoys.
    FDR at T9 = (1/0.8 + 1)/10 = 0.225 (not 0.2 unscaled, not 2/(0.8*10) = 0.25)."""
    b = Builder(18)
    t = b.tags or []
    for i in range(5):
        b.target(t[i], 0.001, 0.01, 100.0 - i)
    for i in range(5, 10):
        b.target(t[i], 0.001, 0.01, 55.0 - i)
    b.decoy([f"rev_{t[10]}"], 0.001, 0.01, 60.0)
    for i in range(11, 18):
        b.decoy([f"rev_{t[i]}"], 0.001, 0.01, float(i - 10))
    b.decoy(["rev_UNKNOWN1"], 0.001, 0.01, 1.0)
    b.decoy(["rev_UNKNOWN2"], 0.001, 0.01, 1.0)
    r = b.run()
    assert r.det.tier == "decoy_picked"
    np.testing.assert_allclose([r.q(g) for g in t[:5]], 0.2, rtol=1e-12)
    np.testing.assert_allclose([r.q(g) for g in t[5:10]], 0.225, rtol=1e-12)
    assert all(math.isnan(r.q(g)) for g in t[10:])
    s = fx.summary(r.det)
    vals = s.loc[s["group_type"] == "tier", "value"].astype(float).to_numpy()
    assert np.isclose(vals, 0.8).any(), "paired share f reported (ADR-0009)"
    assert np.isclose(vals, 1.25).any(), "rescaling factor 1/f reported (ADR-0009)"


# ------------------------------------------------------------------------- example 6 (B.3 MAYU)


@pytest.fixture(scope="module")
def b3_run() -> fx.Run:
    """5,000 eligible genes (uniform sizes), 3,000 hit by 5,000 peptides at q = 0.01, no PEP:
    F(0.01) = 0.01 * 5,000 = 50 (Appendix B.3)."""
    b = Builder(5000)
    t = b.tags or []
    for i in range(3000):
        b.target(t[i], 0.01, None, 50.0)
        if i < 2000:
            b.target(t[i], 0.01, None, 50.0)
    return b.run(sizes=10)


def test_example6_mayu_b3(b3_run: fx.Run) -> None:
    r = b3_run
    assert r.det.tier == "mayu"
    assert fx.curve_metric(r.det, 0.01, "n_detected") == 3000
    assert fx.curve_metric(r.det, 0.01, "fdr") == pytest.approx(0.0067008, rel=1e-4)
    q = r.det.table["gene_q_value"].astype(float).to_numpy()
    np.testing.assert_allclose(q[:3000], 0.0067008, rtol=1e-4)
    assert np.isnan(q[3000:]).all()
    st = r.det.table["detection_status"].astype(str).to_numpy()
    assert (st[:3000] == "detected").all() and (st[3000:] == "not_detected").all()


def test_example6_reports_fp_n_genes_and_fdr_2pep(b3_run: fx.Run) -> None:
    """Statistician items (PR #34): fp and n_genes next to the threshold; fdr_2pep (M4-02 ex 2)
    is reported under its own name and is not the gene FDR."""
    det = b3_run.det
    assert fx.curve_metric(det, 0.01, "fp") == pytest.approx(20.1024, rel=1e-4)
    assert fx.curve_metric(det, 0.01, "n_genes") == 3000
    assert fx.curve_metric(det, 0.01, "fdr_2pep") == pytest.approx(3.3114e-5, rel=1e-3)
    thresholds = sorted(set(fx.curve_rows(det)["group"].astype(float)))
    np.testing.assert_allclose(thresholds, fx.CURVE_THRESHOLDS)


def test_fdr_2pep_only_in_mayu_tier() -> None:
    r = _picked(EX1_T, EX1_D).run()
    assert "fdr_2pep" not in set(fx.curve_rows(r.det)["metric"])


# ------------------------------------------------------------------------- examples 7-10


def _sequences_only() -> Builder:
    b = Builder(3)
    g, h, _k = b.tags or []
    b.target(g, None, None, None)
    b.target(g, None, None, None)
    b.target(h, None, None, None)
    return b


def test_example7_two_peptide_tier() -> None:
    b = _sequences_only()
    g, h, k = b.tags or []
    r = b.run()
    assert_table_shape(r)
    assert r.det.tier == "two_peptide"
    assert r.status(g) == "detected" and r.detected(g)
    assert r.status(h) == "single_peptide_unestimated" and not r.detected(h)
    assert r.status(k) == "not_detected" and not r.detected(k)
    assert all(math.isnan(r.q(x)) for x in (g, h, k))
    assert fx.tier_metric(r.det, "n_detected") == 1
    assert fx.tier_metric(r.det, "n_detected_2pep") == 1


def test_example8_protein_list_tier() -> None:
    b = Builder(2)
    g, h = b.tags or []
    b.protein_hit("UP_G", g)
    r = b.run(no_peptides=True)
    assert_table_shape(r)
    assert r.det.tier == "protein_list"
    assert r.status(g) == "detected" and r.detected(g)
    assert r.status(h) == "not_detected" and not r.detected(h)
    assert math.isnan(r.q(g)) and math.isnan(r.q(h))


def test_protein_list_needs_rank_one_hit() -> None:
    b = Builder(2)
    g, h = b.tags or []
    b.protein_hit("UP_G", g, rank=1)
    b.protein_hit("UP_H", h, rank=2)
    r = b.run(no_peptides=True)
    assert r.detected(g) and not r.detected(h)


@pytest.mark.parametrize("tier", ["decoy_picked", "two_peptide"])
def test_example9_partial_gene_excluded(tier: str) -> None:
    b = Builder(3, partial={tag(0)})
    t = b.tags or []
    seq_only = tier == "two_peptide"
    for _ in range(5):
        b.target(t[0], *((None, None, None) if seq_only else (0.001, 0.01, 90.0)))
    for _ in range(2):
        b.target(t[1], *((None, None, None) if seq_only else (0.001, 0.01, 80.0)))
    if not seq_only:
        b.decoy([f"rev_{t[2]}"], 0.001, 0.01, 1.0)
    r = b.run(_alpha(1.0))
    assert r.det.tier == tier
    assert r.status(t[0]) == "partial_excluded"
    assert not r.detected(t[0])
    assert math.isnan(r.q(t[0]))
    assert r.detected(t[1])


def test_example10_tie_goes_to_the_decoy() -> None:
    b = Builder(2)
    g, h = b.tags or []
    b.target(g, 0.001, 0.01, 50.0)
    b.decoy([f"rev_{g}"], 0.001, 0.01, 50.0)
    b.target(h, 0.001, 0.01, 80.0)
    r = b.run(_alpha(1.0))
    assert r.det.tier == "decoy_picked"
    assert r.status(g) == "not_detected" and math.isnan(r.q(g))
    assert r.detected(h)


@pytest.mark.parametrize(
    ("target", "decoy", "target_wins"),
    [
        ((0.001, 0.01, 50.0), (0.002, 0.0, 99.0), True),  # q first
        ((0.001, 0.01, 50.0), (0.001, 0.02, 99.0), True),  # then lower PEP
        ((0.001, 0.01, 50.0), (0.001, 0.01, 49.0), True),  # then higher score
        ((0.001, 0.01, 50.0), (0.001, 0.005, 1.0), False),
        # PEP only on the decoy side: not a key column, so the score decides (M4-24 rule 5)
        ((0.001, None, 50.0), (0.001, 0.5, 1.0), True),
    ],
)
def test_lexicographic_score_key(
    target: tuple[float, float | None, float], decoy: tuple[float, float, float], target_wins: bool
) -> None:
    b = Builder(1)
    g = (b.tags or [])[0]
    b.target(g, *target)
    b.decoy([f"rev_{g}"], *decoy)
    r = b.run(_alpha(1.0))
    assert r.detected(g) == target_wins
    assert math.isnan(r.q(g)) != target_wins


def test_best_target_of_gene_is_its_best_peptide() -> None:
    """s_g is the BEST accepted target peptide: a weak second peptide does not lose to the decoy."""
    b = Builder(1)
    g = (b.tags or [])[0]
    b.target(g, 0.0001, 0.01, 50.0)
    b.target(g, 0.009, 0.01, 50.0)
    b.decoy([f"rev_{g}"], 0.001, 0.01, 50.0)
    assert b.run(_alpha(1.0)).detected(g)


def test_gene_without_accepted_peptide_not_detected_q_na() -> None:
    b = Builder(3)
    t = b.tags or []
    b.target(t[0], 0.001, 0.01, 50.0)
    b.target(t[1], 0.5, 0.9, 1.0)  # rejected: q > q_value_max
    b.decoy([f"rev_{t[2]}"], 0.001, 0.01, 1.0)
    r = b.run(_alpha(1.0))
    for g in t[1:]:
        assert r.status(g) == "not_detected" and math.isnan(r.q(g)), g


# ------------------------------------------------------------------------- tier choice


def test_no_decoys_with_q_is_mayu() -> None:
    b = Builder(2)
    b.target(tag(0), 0.001, None, None)
    assert b.run().det.tier == "mayu"


def test_no_decoys_with_pep_only_is_mayu() -> None:
    b = Builder(2)
    b.target(tag(0), None, 0.01, None)
    r = b.run()
    assert r.det.tier == "mayu"
    assert not fx.is_na(r.row(tag(0))["gene_q_value"]), "lead decision 5: q derived from PEP"


def test_decoys_without_protein_ids_do_not_give_picked() -> None:
    b = Builder(2)
    b.target(tag(0), 0.001, 0.01, 50.0)
    b.decoy([], 0.001, 0.01, 1.0)
    assert b.run().det.tier == "mayu"


def test_score_only_is_two_peptide() -> None:
    b = Builder(2)
    b.target(tag(0), None, None, 50.0)
    b.target(tag(0), None, None, 40.0)
    assert b.run().det.tier == "two_peptide"


def test_c11_only_q_values_give_mayu() -> None:
    """Without C9 the accepted q comes from C11 ``best_q_input`` (M4-11 decision 4)."""
    b = Builder(2)
    b.target(tag(0), 0.001, 0.01, 50.0)
    b.target(tag(1), 0.002, 0.01, 50.0)
    r = b.run(no_peptides=True)
    assert r.det.tier == "mayu"


def test_nothing_at_all_is_protein_list() -> None:
    r = Builder(2).run(no_peptides=True)
    assert r.det.tier == "protein_list"
    assert not r.det.table["detected"].astype(bool).any()


# ------------------------------------------------------------------------- edge cases


def test_zero_eligible_genes() -> None:
    b = Builder(2, partial={tag(0), tag(1)})
    b.target(tag(0), 0.001, 0.01, 50.0)
    b.decoy([f"rev_{tag(1)}"], 0.001, 0.01, 1.0)
    r = b.run()
    assert r.det.tier == "decoy_picked", "tier still chosen"
    assert "no_eligible_genes" in r.det.warnings
    assert fx.curve_rows(r.det).empty
    assert set(r.det.table["detection_status"]) == {"partial_excluded"}


def test_q_zero_handled_by_ranks() -> None:
    b = Builder(6)
    t = b.tags or []
    for i in range(5):
        b.target(t[i], 0.0, 0.001 * (i + 1), 50.0)
    b.decoy([f"rev_{t[5]}"], 0.0, 0.5, 50.0)  # worse PEP than every target
    r = b.run()
    q = np.array([r.q(g) for g in t[:5]])
    assert np.isfinite(q).all()
    np.testing.assert_allclose(q, 0.2, rtol=1e-12)


def test_protein_map_priority_and_methods() -> None:
    b = Builder(3)
    g1, g2, g3 = b.tags or []
    b.protein_hit(g1, g2)  # (1) locus tag wins over the C12 hit
    b.protein_hit("UP_1", g2)  # (2) rank-1 C12 gene wins over the vote
    b.protein_hit("UP_2", g1, rank=2)  # rank 2 ignored
    for _ in range(2):
        b.target(g3, protein_ids=["UP_1", "UP_3"])
    b.target(g1, protein_ids=["UP_3", "UP_4", g1])
    b.target(g2, protein_ids=["UP_4"])  # UP_4: one peptide on g1, one on g2 -> tie
    det = _det()
    pmap = det.protein_gene_map(b.genes_frame(), b.peptides_frame(), fx.eo.hits_frame(b.hits),
                                fx.ao.protein_hits_frame(b.protein_hit_rows))  # fmt: skip
    assert fx.mapped(pmap, g1) == (g1, "locus_tag")
    assert fx.mapped(pmap, "UP_1") == (g2, "protein_hit")
    assert fx.mapped(pmap, "UP_2") == (None, None)
    assert fx.mapped(pmap, "UP_3") == (g3, "peptide_vote")  # 2 peptides on g3, 1 on g1
    assert fx.mapped(pmap, "UP_4") == (None, None)
    methods = set(pmap.loc[pmap["gene"].notna(), "method"].astype(str))
    assert methods <= {"locus_tag", "protein_hit", "peptide_vote"}


def test_contaminant_ids_never_mapped() -> None:
    b = Builder(2)
    g = (b.tags or [])[0]
    b.protein_hit("CON__P02768", g)
    for _ in range(3):
        b.target(g, protein_ids=["Cont_TRYP", "CON__P02768"], contaminant=True)
    det = _det()
    pmap = det.protein_gene_map(b.genes_frame(), b.peptides_frame(), fx.eo.hits_frame(b.hits),
                                fx.ao.protein_hits_frame(b.protein_hit_rows))  # fmt: skip
    assert fx.mapped(pmap, "CON__P02768") == (None, None)
    assert fx.mapped(pmap, "Cont_TRYP") == (None, None)


def test_contaminant_evidence_counts_for_detection() -> None:
    b = _sequences_only()
    _g, h, _k = b.tags or []
    b.target(h, None, None, None, contaminant=True)  # H now has 2 accepted, one contaminant
    r = b.run()
    assert bool(r.evidence.set_index("locus_tag").loc[h, "contaminant_evidence"])
    assert r.status(h) == "detected"


def test_pep_partial_uses_t_times_p_and_warns() -> None:
    from giap.stats.mayu import mayu_gene_fdr

    b = Builder(40)
    t = b.tags or []
    for i in range(20):
        b.target(t[i], 0.001 * (1 + i % 5), 0.01 if i % 2 else None, 50.0)
    r = b.run(sizes=5)
    assert r.det.tier == "mayu"
    assert "pep_partial" in r.det.warnings
    best = r.evidence["best_peptide_q"].astype(float).to_numpy()
    acc_q = r.accepted.drop_duplicates("peptide")["q_value"].astype(float).to_numpy()
    ref = mayu_gene_fdr(best, acc_q, None, np.full(40, 5.0), cfg().detection.mayu_size_bins)
    q = r.det.table["gene_q_value"].astype(float).to_numpy()
    np.testing.assert_allclose(q, ref.gene_q, rtol=1e-12, equal_nan=True)


def test_mayu_uses_pep_when_every_peptide_has_one() -> None:
    from giap.stats.mayu import mayu_gene_fdr

    b = Builder(40)
    t = b.tags or []
    for i in range(20):
        b.target(t[i], 0.001 * (1 + i % 5), 0.002 * (1 + i % 3), 50.0)
    r = b.run(sizes={g: 1 + i % 7 for i, g in enumerate(t)})
    assert "pep_partial" not in r.det.warnings
    acc = r.accepted.drop_duplicates("peptide")
    sizes = np.array([1 + i % 7 for i in range(40)], dtype=float)
    ref = mayu_gene_fdr(r.evidence["best_peptide_q"].astype(float).to_numpy(),
                        acc["q_value"].astype(float).to_numpy(),
                        acc["pep"].astype(float).to_numpy(), sizes,
                        cfg().detection.mayu_size_bins)  # fmt: skip
    q = r.det.table["gene_q_value"].astype(float).to_numpy()
    np.testing.assert_allclose(q, ref.gene_q, rtol=1e-12, equal_nan=True)


# ------------------------------------------------------------------------- input validation


@pytest.mark.parametrize("bad_pep", [-0.1, 1.5])
def test_pep_outside_unit_interval_rejected(bad_pep: float) -> None:
    b = Builder(3)
    b.target(tag(0), 0.001, bad_pep, 50.0)
    b.target(tag(1), 0.001, 0.01, 50.0)
    with pytest.raises(ValueError, match=r"(?i)pep"):
        b.run()


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_non_finite_sizes_rejected(bad: float) -> None:
    b = Builder(3)
    b.target(tag(0), 0.001, 0.01, 50.0)
    with pytest.raises(ValueError):
        b.run(sizes={tag(0): 4, tag(1): bad, tag(2): 4})


def test_bad_mayu_bins_rejected() -> None:
    b = Builder(3)
    b.target(tag(0), 0.001, 0.01, 50.0)
    with pytest.raises(ValueError):
        b.run(cfg(detection={"mayu_size_bins": 0}))


def test_unspecific_search_warning() -> None:
    b = Builder(30)
    t = b.tags or []
    for i, g in enumerate(t):
        b.target(g, 0.001, 0.01, 50.0, peptide=fx.peptide_name(1000 + i, last="L"))
    warns = b.run().det.warnings
    assert any("unspecific" in w for w in warns), warns
    b2 = Builder(30)
    for g in b2.tags or []:
        b2.target(g, 0.001, 0.01, 50.0)
    assert not any("unspecific" in w for w in b2.run().det.warnings)


def test_detection_is_frozen_dataclass() -> None:
    r = _picked(EX1_T, EX1_D).run()
    with pytest.raises((AttributeError, TypeError)):
        r.det.tier = "mayu"
    assert isinstance(r.det.warnings, tuple)
    assert isinstance(r.det.table, pd.DataFrame)
