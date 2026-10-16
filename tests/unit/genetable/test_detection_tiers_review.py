"""M4-12 lead decisions 5-9 after the statistician review (PR #92, M4-12-detection-tiers spec).

5. PEP-only input (spec §2.1 decision 4: "q-values or PEPs only -> MAYU"): gene thresholds use the
   peptide q derived from PEP, the running mean of the sorted PEPs of the distinct accepted
   peptides (q_i = mean of the PEPs <= PEP_i, the quantity the MAYU-with-PEP path sums); never
   0 genes silently.
6. Contaminant decoys are neither paired nor unpaired. Decision 12 (correcting decision 11): a decoy
   is a contaminant decoy only when ALL its stripped IDs carry an ``input.contaminant_prefixes``
   entry; a mixed decoy (``rev_P;rev_CON__P``) pairs on its non-contaminant IDs (paired when P
   maps, unpaired otherwise).
7. Decoy-prefix stripping removes the longest matching prefix first (``REV__X`` -> ``X``) whatever
   the order of ``input.decoy_prefixes``.
8. ``Detection.protein_gene_map`` is the protein -> gene map ``detect`` used (columns
   ``protein_id``, ``gene``, ``method``).
9. A decoy pairs only when every listed protein maps (to one gene); ``unspecific_search`` fires
   at more than 50 % non-tryptic peptides, not at exactly 50 %.

Inputs come from the merged producers (tests/helpers/detection_tiers_fixtures.py: M4-11
``accepted_peptides`` / ``gene_evidence``, M4-10 ``covariate_table``); the giapsim (M4-09) versions
of 5, 6 and 8 are in test_detection_tiers_contract.py.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_fixtures as fx
from detection_tiers_fixtures import Builder, cfg, tag


def _det() -> Any:
    from giap.genetable import detection

    return detection


def _pmap_rows(pmap: pd.DataFrame) -> set[tuple[str, str, str]]:
    m = pmap[pmap["gene"].notna()]
    return set(zip(m["protein_id"].astype(str), m["gene"].astype(str), m["method"].astype(str),
                   strict=True))  # fmt: skip


# ------------------------------------------------------------------ decision 5: PEP-only input


def test_pep_only_reviewer_reproduction_detects_genes() -> None:
    """Review finding 5: 300 genes x 2 peptides at PEP 1e-4, no q -> every gene detected."""
    b = Builder(300)
    for g in b.tags or []:
        b.target(g, None, 1e-4, None)
        b.target(g, None, 1e-4, None)
    r = b.run(sizes=20)
    assert r.det.tier == "mayu"
    q = r.det.table["gene_q_value"].astype(float)
    assert q.notna().all(), "PEP-only gene q must come from the PEP-derived peptide q"
    assert r.det.table["detected"].astype(bool).all()
    assert fx.tier_metric(r.det, "n_detected") == 300


def test_pep_only_gene_q_is_mayu_on_running_mean_of_pep() -> None:
    from giap.stats.mayu import mayu_gene_fdr

    b = Builder(40)
    t = b.tags or []
    peps = [0.0005 * (i + 1) for i in range(20)]  # distinct, so no tie convention is needed
    for i in range(20):
        b.target(t[i], None, peps[i], None)
    sizes = {g: 1 + i % 7 for i, g in enumerate(t)}
    r = b.run(sizes=sizes)
    assert r.det.tier == "mayu"
    derived = np.cumsum(peps) / np.arange(1, 21)
    best = np.full(40, np.nan)
    best[:20] = derived
    ref = mayu_gene_fdr(best, derived, np.asarray(peps), np.array([sizes[g] for g in t], float),
                        cfg().detection.mayu_size_bins)  # fmt: skip
    q = r.det.table["gene_q_value"].astype(float).to_numpy()
    np.testing.assert_allclose(q, ref.gene_q, rtol=1e-12, equal_nan=True)
    assert np.isfinite(q[:20]).all()


def test_pep_only_q_threshold_curve_is_populated() -> None:
    b = Builder(300)
    for g in b.tags or []:
        b.target(g, None, 1e-4, None)
        b.target(g, None, 1e-4, None)
    det = b.run(sizes=20).det
    assert fx.curve_metric(det, 0.01, "n_detected") == 300


# ------------------------------------------------------------------ decisions 6 and 7: decoy IDs


def _eight_gene_decoys(contaminant_ids: list[str]) -> Builder:
    """10 genes with a strong target each, 8 decoys ``rev_<tag>`` below every target, plus one
    decoy per entry of ``contaminant_ids``."""
    b = Builder(10)
    t = b.tags or []
    for i, g in enumerate(t):
        b.target(g, 0.001, 0.01, 100.0 + i)
    for i in range(8):
        b.decoy([f"rev_{t[i]}"], 0.001, 0.01, 1.0 + i)
    for i, cid in enumerate(contaminant_ids):
        b.decoy([cid], 0.001, 0.01, 20.0 + i)
    return b


def test_contaminant_decoys_not_counted_8_plus_3() -> None:
    """Decision 6: 8 gene + 3 pure ``rev_CON__`` decoys (every stripped ID a contaminant,
    decision 12) stay in the picked tier (8/8, not 8/11)."""
    ids = ["rev_CON__P02768", "rev_CON__P00761", "rev_CON__P04264"]
    r = _eight_gene_decoys(ids).run()
    assert r.det.tier == "decoy_picked"
    assert "decoy_pairing_low" not in r.det.warnings
    assert fx.tier_metric(r.det, "n_decoys_paired") == 8
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0
    assert fx.tier_metric(r.det, "decoy_paired_share") == pytest.approx(1.0)


def test_contaminant_decoys_do_not_change_gene_q() -> None:
    ids = ["rev_CON__P02768", "rev_CON__P00761", "rev_CON__P04264"]
    base = _eight_gene_decoys([]).run()
    with_con = _eight_gene_decoys(ids).run()
    pd.testing.assert_frame_equal(with_con.det.table, base.det.table)


@pytest.mark.parametrize("prefix", ["CON__", "Cont_", "CONT_", "cRAP"])
def test_every_contaminant_prefix_excluded(prefix: str) -> None:
    r = _eight_gene_decoys([f"rev_{prefix}X{i}" for i in range(3)]).run()
    assert r.det.tier == "decoy_picked"
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0


def test_contaminant_decoys_follow_configured_prefixes() -> None:
    c = cfg(input={"contaminant_prefixes": ["MYCON_"]})
    r = _eight_gene_decoys([f"rev_MYCON_{i}" for i in range(3)]).run(c)
    assert r.det.tier == "decoy_picked"
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0
    r2 = _eight_gene_decoys([f"rev_CON__X{i}" for i in range(3)]).run(c)
    assert fx.tier_metric(r2.det, "n_decoys_unpaired") == 3  # CON__ is not a contaminant here
    assert r2.det.tier == "mayu"


def test_unknown_decoys_still_count_as_unpaired() -> None:
    """Only all-contaminant decoys leave the counts: 8 paired + 3 unknown -> 8/11 < 0.8 -> mayu."""
    r = _eight_gene_decoys([f"rev_UNKNOWN{i}" for i in range(3)]).run()
    assert r.det.tier == "mayu"
    assert "decoy_pairing_low" in r.det.warnings


@pytest.mark.parametrize("order", [["REV_", "REV__"], ["REV__", "REV_"]])
def test_longest_decoy_prefix_stripped_first(order: list[str]) -> None:
    """Decision 7: ``REV__X`` -> ``X`` (MaxQuant decoys), never ``_X``."""
    b = Builder(3)
    t = b.tags or []
    for g in t:
        b.decoy([f"REV__{g}"])
    det = _det()
    peps = b.peptides_frame()
    pmap = det.protein_gene_map(b.genes_frame(), peps, fx.eo.hits_frame([]),
                                fx.ao.protein_hits_frame([]))  # fmt: skip
    pairing = det.pair_decoys(peps, pmap, order)
    assert pairing.set_index("peptide").loc[peps["peptide"], "gene"].tolist() == t


def test_maxquant_contaminant_decoy_excluded_via_longest_prefix() -> None:
    """``REV__CON__P02768`` strips to ``CON__P02768`` (decision 7), a contaminant (decision 6)."""
    prefixes = [*cfg().input.decoy_prefixes, "REV__"]
    c = cfg(input={"decoy_prefixes": prefixes})
    b = Builder(10)
    t = b.tags or []
    for i, g in enumerate(t):
        b.target(g, 0.001, 0.01, 100.0 + i)
    for i in range(8):
        b.decoy([f"REV__{t[i]}"], 0.001, 0.01, 1.0 + i)
    for i in range(3):
        b.decoy([f"REV__CON__P0276{i}"], 0.001, 0.01, 20.0 + i)
    r = b.run(c)
    assert r.det.tier == "decoy_picked"
    assert fx.tier_metric(r.det, "n_decoys_paired") == 8
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0


def test_all_contaminant_multi_id_decoys_not_counted() -> None:
    """Decision 12: ``rev_CON__A;rev_Cont_B`` has only contaminant IDs -> out of the counts."""
    ids = [f"rev_CON__P0{i};rev_Cont_T{i}" for i in range(3)]
    r = _eight_gene_decoys(ids).run()
    assert r.det.tier == "decoy_picked"
    assert fx.tier_metric(r.det, "n_decoys_paired") == 8
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0


# ------------------------------------------------------------------ decision 12: mixed decoys


def _pair(b: Builder) -> pd.DataFrame:
    det = _det()
    peps = b.peptides_frame()
    pmap = det.protein_gene_map(b.genes_frame(), peps, fx.eo.hits_frame([]),
                                fx.ao.protein_hits_frame(b.protein_hit_rows))  # fmt: skip
    return det.pair_decoys(peps, pmap, cfg().input.decoy_prefixes).set_index("peptide")


@pytest.mark.parametrize("ids", [["rev_P1;rev_CON__P1"], ["rev_P1", "rev_CON__P02768"]])
def test_mixed_decoy_pairs_when_its_target_maps(ids: list[str]) -> None:
    """``rev_P1;rev_CON__P1`` pairs on P1 (it maps to g): a paired row, not dropped."""
    b = Builder(2)
    g = (b.tags or [])[0]
    b.protein_hit("P1", g)
    b.decoy(ids)
    pairing = _pair(b)
    p = b.peptides_frame()["peptide"].iloc[0]
    assert p in pairing.index, "a mixed decoy is not a contaminant decoy (decision 12)"
    assert pairing.loc[p, "gene"] == g


@pytest.mark.parametrize("ids", [["rev_UNMAPPED;rev_CON__P1"], ["rev_UNMAPPED", "rev_CON__P1"]])
def test_mixed_decoy_unpaired_when_its_target_does_not_map(ids: list[str]) -> None:
    b = Builder(2)
    b.decoy(ids)
    pairing = _pair(b)
    p = b.peptides_frame()["peptide"].iloc[0]
    assert p in pairing.index, "a mixed decoy is not a contaminant decoy (decision 12)"
    assert fx.is_na(pairing.loc[p, "gene"])


def test_mixed_decoys_paired_in_detect() -> None:
    """8 gene decoys + 3 mixed ``rev_<tag>;rev_CON__..`` whose tag maps -> 11 paired, 0 unpaired."""
    b = _eight_gene_decoys([])
    t = b.tags or []
    for i in range(3):
        b.decoy([f"rev_{t[7 + i]};rev_CON__P0276{i}"], 0.001, 0.01, 20.0 + i)
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 11
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0
    assert r.det.tier == "decoy_picked"


def test_mixed_decoys_unpaired_in_detect() -> None:
    """8 gene decoys + 3 mixed decoys whose non-contaminant ID does not map: 8/11 < 0.8 -> mayu
    (under "any" they would vanish and keep ``decoy_picked``: the attack of decision 12)."""
    ids = [f"rev_UNMAPPED{i};rev_CON__P0276{i}" for i in range(3)]
    r = _eight_gene_decoys(ids).run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 8
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 3
    assert fx.tier_metric(r.det, "decoy_paired_share") == pytest.approx(8 / 11)
    assert r.det.tier == "mayu"
    assert "decoy_pairing_low" in r.det.warnings


def test_pairing_share_on_non_contaminant_ids_with_many_mixed_decoys() -> None:
    """10 decoys, 4 mixed (40 % >= 20 %):
    5 plain ``rev_<tag>`` -> paired; 2 ``rev_<tag>;rev_CON__..`` -> paired (tag maps);
    2 ``rev_UNMAPPED;rev_Cont_..`` -> unpaired; 1 ``rev_CON__..;rev_cRAP..`` -> excluded.
    Share = (5 + 2) / (5 + 2 + 2) = 7/9 = 0.778 < 0.8 -> mayu. ("Any" gives 5/5 = 1 -> picked.)"""
    b = Builder(10)
    t = b.tags or []
    for i, g in enumerate(t):
        b.target(g, 0.001, 0.01, 100.0 + i)
    for i in range(5):
        b.decoy([f"rev_{t[i]}"], 0.001, 0.01, 1.0 + i)
    for i in range(2):
        b.decoy([f"rev_{t[5 + i]};rev_CON__P0{i}"], 0.001, 0.01, 10.0 + i)
    for i in range(2):
        b.decoy([f"rev_UNMAPPED{i};rev_Cont_T{i}"], 0.001, 0.01, 20.0 + i)
    b.decoy(["rev_CON__P09;rev_cRAP009"], 0.001, 0.01, 30.0)
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 7
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 2
    assert fx.tier_metric(r.det, "decoy_paired_share") == pytest.approx(7 / 9)
    assert r.det.tier == "mayu"
    assert "decoy_pairing_low" in r.det.warnings


# ------------------------------------------------------------------ decision 8: exposed map


def test_detection_exposes_the_map_detect_used() -> None:
    b = Builder(3)
    g1, g2, g3 = b.tags or []
    b.protein_hit("UP_1", g2)
    for _ in range(2):
        b.target(g3, protein_ids=["UP_3"])
    b.target(g1, protein_ids=[g1])
    r = b.run()
    pm = r.det.protein_gene_map
    assert isinstance(pm, pd.DataFrame)
    assert {"protein_id", "gene", "method"} <= set(pm.columns)
    c = cfg()
    ref = _det().protein_gene_map(r.genes, r.peptides, r.peptide_hits, r.protein_hits,
                                  accepted=r.accepted,
                                  contaminant_prefixes=c.input.contaminant_prefixes)  # fmt: skip
    assert _pmap_rows(pm) == _pmap_rows(ref)
    assert fx.mapped(pm, "UP_1") == (g2, "protein_hit")
    assert fx.mapped(pm, "UP_3") == (g3, "peptide_vote")
    assert fx.mapped(pm, g1) == (g1, "locus_tag")


def test_exposed_map_uses_configured_contaminant_prefixes() -> None:
    b = Builder(2)
    g = (b.tags or [])[0]
    b.protein_hit("MYCON_1", g)
    b.protein_hit("CON__P02768", g)
    b.target(g)
    r = b.run(cfg(input={"contaminant_prefixes": ["MYCON_"]}))
    assert fx.mapped(r.det.protein_gene_map, "MYCON_1") == (None, None)
    assert fx.mapped(r.det.protein_gene_map, "CON__P02768") == (g, "protein_hit")


def test_exposed_map_votes_only_over_accepted_peptides() -> None:
    """Peptides failing q <= q_value_max do not vote in the map ``detect`` used."""
    b = Builder(2)
    g, h = b.tags or []
    b.target(g, 0.001, protein_ids=["UP_9"])
    for _ in range(3):
        b.target(h, 0.5, protein_ids=["UP_9"])  # rejected: q > q_value_max
    r = b.run()
    assert fx.mapped(r.det.protein_gene_map, "UP_9") == (g, "peptide_vote")


def test_exposed_map_present_in_every_tier() -> None:
    b = Builder(2)
    b.protein_hit("UP_1", tag(0))
    r = b.run(no_peptides=True)
    assert r.det.tier == "protein_list"
    assert fx.mapped(r.det.protein_gene_map, "UP_1") == (tag(0), "protein_hit")


# ------------------------------------------------------------------ decision 9: surviving mutants


@pytest.mark.parametrize("ids", [["rev_P1;rev_UNKNOWN"], ["rev_P1", "rev_UNKNOWN"]])
def test_decoy_with_one_unmapped_protein_is_unpaired(ids: list[str]) -> None:
    b = Builder(2)
    g = (b.tags or [])[0]
    b.protein_hit("P1", g)
    b.decoy(ids)
    b.decoy(["rev_P1"])
    det = _det()
    peps = b.peptides_frame()
    pmap = det.protein_gene_map(b.genes_frame(), peps, fx.eo.hits_frame([]),
                                fx.ao.protein_hits_frame(b.protein_hit_rows))  # fmt: skip
    pairing = det.pair_decoys(peps, pmap, cfg().input.decoy_prefixes).set_index("peptide")
    p_mixed, p_one = peps["peptide"].tolist()
    assert fx.is_na(pairing.loc[p_mixed, "gene"])
    assert pairing.loc[p_one, "gene"] == g


def test_decoy_with_one_unmapped_protein_counts_unpaired_in_detect() -> None:
    b = Builder(10)
    t = b.tags or []
    for i, g in enumerate(t):
        b.target(g, 0.001, 0.01, 100.0 + i)
    for i in range(7):
        b.decoy([f"rev_{t[i]}"], 0.001, 0.01, 1.0 + i)
    for i in range(3):
        b.decoy([f"rev_{t[7 + i]};rev_UNKNOWN{i}"], 0.001, 0.01, 10.0 + i)
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 7
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 3
    assert r.det.tier == "mayu"


def _unspecific_run(n_nontryptic: int, n: int = 30) -> Any:
    b = Builder(n)
    for i, g in enumerate(b.tags or []):
        last = "L" if i < n_nontryptic else "K"
        b.target(g, 0.001, 0.01, 50.0, peptide=fx.peptide_name(2000 + i, last=last))
    return b.run()


def test_unspecific_search_not_flagged_at_exactly_half() -> None:
    r = _unspecific_run(15)
    assert r.det.tier == "mayu"
    assert not any("unspecific" in w for w in r.det.warnings), r.det.warnings


def test_unspecific_search_flagged_just_above_half() -> None:
    r = _unspecific_run(16)
    assert any("unspecific" in w for w in r.det.warnings), r.det.warnings
