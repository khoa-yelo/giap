"""M4-18 ``giap.context.regions_table``: plasmid rows, prophage subset rows, the table assembler.

docs/specs/M4-18-regions-stats-table.md examples 1-7 and 9, edge cases, and the lead decisions of
2026-10-20 (binding): raw log2 HL for plasmids without imputed genes (2), the reference-replicon
rule (3), the depth offset (4), NA + ``small_region`` for 0 eligible genes and ``observed`` vs
``hl_shift_log2`` (5), ``null_model = none`` with NA p / q and ``n_perm = 0`` and no BH recompute
(6), the plasmid's own RNG stream (7) and descriptive subset rows (8).

Inputs are hand-built §9.3 gene tables (tests/helpers/regions_stats_data.py on top of the M4-17
helper); expected values come from the independent oracles there.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import region_tests_data as rd
import regions_stats_data as sd
from region_tests_data import Contig

B = 200


_T, _rng, _plasmids, _tests, _table = sd.module, sd.rng, sd.plasmid_rows, sd.tests, sd.table
_a, _pairs, chr_contig, plasmid_contig = sd.abund, sd.pairs, sd.chr_contig, sd.plasmid_contig
example_1 = sd.example_1


def _assert_no_test(row: pd.Series) -> None:
    assert sd.na(row["p_value"]) and sd.na(row["q_value"]), row.to_dict()
    assert row["null_model"] == "none"
    assert int(row["n_perm"]) == 0


# --------------------------------------------------------------------------- example 1


def test_example_1_two_plasmid_rows() -> None:
    g = example_1()
    res = _plasmids(g)
    assert len(res) == 2
    assert sorted(res["stat"]) == ["hl_shift", "o_minus_e"]
    for _, row in res.iterrows():
        assert row["region_id"] == "plasmid_1" and row["region_type"] == "plasmid"
        assert row["contig"] == "p1"
        _assert_no_test(row)
        assert row["method"] == sd.PLASMID_METHOD
        assert row["label"] == sd.PLASMID_LABEL
        assert "dosage_plus_regulation" in sd.flags(row["flags"])
        if "subset" in res:
            assert row["subset"] == "all"


def test_example_1_hl_is_raw_shift_against_chromosome() -> None:
    g = example_1()
    hl = sd.unadjusted(_plasmids(g), "plasmid_1", "hl_shift")
    want = sd.plasmid_hl(g.genes, "p1", ["chr1"])
    assert abs(want - 4.0) < 0.3  # the fixture: plasmid median 4 units above
    assert hl["effect"] == pytest.approx(want, abs=1e-12)
    assert hl["hl_shift_log2"] == pytest.approx(want, abs=1e-12)
    assert hl["ci_low"] <= hl["effect"] <= hl["ci_high"]
    assert hl["ci_high"] - hl["ci_low"] > 0


def test_example_1_o_minus_e_is_the_plasmids_own() -> None:
    """Lead decision 1: the plasmid's own O - E (E = sum of p_detect), no chromosome-rate row."""
    g = example_1()
    res = _plasmids(g)
    oe = sd.unadjusted(res, "plasmid_1", "o_minus_e")
    want = sd.detection_counts(g.genes, ["p1"])
    for c in ("observed", "expected", "o_minus_e", "o_over_e", "n_genes", "n_eligible"):
        assert float(oe[c]) == pytest.approx(want[c], abs=1e-9), c
    assert oe["effect"] == pytest.approx(want["o_minus_e"], abs=1e-9)
    assert oe["ci_low"] <= oe["effect"] <= oe["ci_high"]


def test_detection_columns_filled_on_both_plasmid_rows() -> None:
    """Decision 5: on the hl_shift row ``observed`` is the detected count, the HL value is in
    ``hl_shift_log2``; every row carries all detection columns."""
    g = example_1()
    res = _plasmids(g)
    oe = sd.unadjusted(res, "plasmid_1", "o_minus_e")
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    want = sd.detection_counts(g.genes, ["p1"])
    for c in ("observed", "expected", "o_minus_e", "o_over_e"):
        assert float(hl[c]) == pytest.approx(want[c], abs=1e-9), c
    assert oe["hl_shift_log2"] == pytest.approx(hl["effect"], abs=1e-12)
    assert hl["observed"] != pytest.approx(hl["effect"])


# --------------------------------------------------------------------------- decision 2 (mutants)


def test_hl_not_residualised_per_replicon() -> None:
    """Per-replicon residuals centre each replicon and give HL ~ 0 whatever the dosage; the raw
    contrast must keep the 4 log2 shift (statistician's mutant)."""
    g = example_1()
    hl = sd.unadjusted(_plasmids(g), "plasmid_1", "hl_shift")["effect"]
    assert hl > 3.5


def test_hl_not_residualised_on_covariates() -> None:
    """Abundance tracks length_aa and plasmid genes are longer: a residualised HL (pooled or per
    replicon) differs from the raw oracle."""
    n, m = 120, 60
    la_c = [60 + 3 * (i % 20) for i in range(n)]
    la_p = [150 + 3 * (i % 20) for i in range(m)]
    a_c = [8.0 + 0.03 * x + 0.01 * i for i, x in enumerate(la_c)]
    a_p = [8.0 + 0.03 * x + 0.01 * i for i, x in enumerate(la_p)]
    c = Contig("chr1", [1] * n, [0.6] * n, a=a_c, tu=_pairs("chr1", n), length_aa=la_c)
    p = Contig("p1", [1] * m, [0.7] * m, a=a_p, tu=_pairs("p1", m), kind="plasmid", length_aa=la_p)
    g = sd.genome([c, p], plasmids=["p1"])
    hl = sd.unadjusted(_plasmids(g), "plasmid_1", "hl_shift")
    want = sd.plasmid_hl(g.genes, "p1", ["chr1"])
    assert want > 2.0
    assert hl["effect"] == pytest.approx(want, abs=1e-12)


def test_imputed_genes_excluded_from_hl() -> None:
    """Decision 2: imputed genes (M4-14 ``abundance_imputed``; the genome-wide fill) drop out of
    both sides of the plasmid contrast; they still count as detected."""
    n, m = 120, 60
    a_p = _a(m, 14.3)
    imputed_idx = [i for i in range(m) if i % 4 != 0][:20]  # detected genes
    for i in imputed_idx:
        a_p[i] = 10.5  # the chromosome-dominated fill value
    a_c = _a(n, 10.0)
    c_imp = [i for i in range(n) if i % 3 != 0][:10]
    for i in c_imp:
        a_c[i] = 10.5
    c = Contig("chr1", [0 if i % 3 == 0 else 1 for i in range(n)], [0.6] * n, a=a_c,
               tu=_pairs("chr1", n))  # fmt: skip
    p = plasmid_contig(n=m)
    p.a = a_p
    g0 = sd.genome([c, p], plasmids=["p1"])
    tags = [f"p1_g{i:05d}" for i in imputed_idx] + [f"chr1_g{i:05d}" for i in c_imp]
    g = sd.genome([c, p], plasmids=["p1"], imputed=tags)
    res = _plasmids(g)
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    want = sd.plasmid_hl(g.genes, "p1", ["chr1"])
    with_imputed = sd.plasmid_hl(g0.genes, "p1", ["chr1"])
    assert abs(want - with_imputed) > 0.3  # the fixture separates the two readings
    assert hl["effect"] == pytest.approx(want, abs=1e-12)
    oe = sd.unadjusted(res, "plasmid_1", "o_minus_e")
    assert float(oe["observed"]) == sd.detection_counts(g.genes, ["p1"])["observed"]


def test_partial_and_unquantified_genes_left_out() -> None:
    """Partial genes are not eligible; detected genes without abundance drop out of HL."""
    m = 60
    part = [i in (0, m - 1) for i in range(m)]
    p = plasmid_contig(n=m, partial=part)
    a = list(p.a or [])
    a[5] = math.nan
    p.a = a
    g = sd.genome([chr_contig(), p], plasmids=["p1"])
    res = _plasmids(g)
    oe = sd.unadjusted(res, "plasmid_1", "o_minus_e")
    want = sd.detection_counts(g.genes, ["p1"])
    assert int(oe["n_genes"]) == m and int(oe["n_eligible"]) == m - 2
    assert float(oe["expected"]) == pytest.approx(want["expected"], abs=1e-9)
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    assert hl["effect"] == pytest.approx(sd.plasmid_hl(g.genes, "p1", ["chr1"]), abs=1e-12)


# --------------------------------------------------------------------------- CIs


def test_hl_ci_resamples_plasmid_and_chromosome_as_strata() -> None:
    """One TU per replicon: with strata (plasmid / chromosome) every resample is the data, so the
    CI collapses onto HL; pooled resampling would draw two chromosome TUs (empty plasmid side)."""
    n, m = 40, 20
    c = chr_contig(n=n, tu=["TU_c"] * n)
    p = plasmid_contig(n=m, tu=["TU_p"] * m)
    g = sd.genome([c, p], plasmids=["p1"])
    res = _plasmids(g)
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    assert hl["ci_low"] == pytest.approx(hl["effect"], abs=1e-12)
    assert hl["ci_high"] == pytest.approx(hl["effect"], abs=1e-12)
    oe = sd.unadjusted(res, "plasmid_1", "o_minus_e")
    assert oe["ci_low"] == pytest.approx(oe["effect"], abs=1e-12)
    assert oe["ci_high"] == pytest.approx(oe["effect"], abs=1e-12)


def test_o_minus_e_ci_resamples_the_plasmids_tus() -> None:
    """Three plasmid TUs with sums of d - p of 1.0, -1.0 and 0.0: every replicate is a sum of 3
    draws of the TU sums, so the CI lies in [-3, 3] (and is not degenerate)."""
    d = [1, 1, 0, 0, 1, 0]
    p = [0.5] * 6
    tu = ["TU_a", "TU_a", "TU_b", "TU_b", "TU_c", "TU_c"]  # sums 1.0, -1.0, 0.0
    pl = Contig("p1", d, p, a=[14.0, 15.0, 0, 0, 16.0, 0], tu=tu, kind="plasmid")
    g = sd.genome([chr_contig(), pl], plasmids=["p1"])
    oe = sd.unadjusted(_plasmids(g), "plasmid_1", "o_minus_e")
    assert oe["effect"] == pytest.approx(0.0, abs=1e-12)
    assert -3.0 - 1e-9 <= oe["ci_low"] <= oe["ci_high"] <= 3.0 + 1e-9
    assert oe["ci_high"] > oe["ci_low"]


# ----------------------------------------------------------------------- examples 2 and 7: depth


def _adjusted(res: pd.DataFrame) -> pd.DataFrame:
    return res[res["label"].astype(str) == sd.ADJUSTED_LABEL]


def test_example_2_depth_adjusted_row() -> None:
    g = example_1()
    res = _plasmids(g, sd.depth_frame({"chr1": 50.0, "p1": 400.0}))
    adj = _adjusted(res)
    assert len(adj) == 1
    row = adj.iloc[0]
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    assert row["stat"] == "hl_shift" and row["region_id"] == "plasmid_1"
    assert row["effect"] == pytest.approx(hl["effect"] - 3.0, abs=1e-12)
    assert row["ci_low"] == pytest.approx(hl["ci_low"] - 3.0, abs=1e-12)
    assert row["ci_high"] == pytest.approx(hl["ci_high"] - 3.0, abs=1e-12)
    _assert_no_test(row)
    assert sd.DEPTH_FLAG in sd.flags(row["flags"])  # lead decision 13
    assert len(res) == 3
    for stat in ("o_minus_e", "hl_shift"):
        assert sd.DEPTH_FLAG not in sd.flags(sd.unadjusted(res, "plasmid_1", stat)["flags"])


def test_depth_does_not_change_the_unadjusted_rows() -> None:
    g = example_1()
    plain = _plasmids(g)
    res = _plasmids(g, sd.depth_frame({"chr1": 50.0, "p1": 400.0}))
    for stat in ("o_minus_e", "hl_shift"):
        a, b = sd.unadjusted(plain, "plasmid_1", stat), sd.unadjusted(res, "plasmid_1", stat)
        for c in ("effect", "ci_low", "ci_high", "observed", "expected"):
            assert float(a[c]) == float(b[c]), (stat, c)


def test_chromosome_depth_is_length_weighted() -> None:
    c1 = chr_contig("chr1", n=100)
    c2 = chr_contig("chr2", n=20, kind="chromosome")
    g = sd.genome([c1, c2, plasmid_contig()], plasmids=["p1"])
    depth = {"chr1": 50.0, "chr2": 150.0, "p1": 400.0}
    res = _plasmids(g, sd.depth_frame(depth))
    dc = sd.weighted_depth(g.replicons, depth, ["chr1", "chr2"])
    off = math.log2(400.0 / dc)
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    assert hl["effect"] == pytest.approx(sd.plasmid_hl(g.genes, "p1", ["chr1", "chr2"]), abs=1e-12)
    adj = _adjusted(res).iloc[0]
    assert adj["effect"] == pytest.approx(hl["effect"] - off, abs=1e-12)


@pytest.mark.parametrize("missing", ["p1", "chr1"])
def test_example_7_depth_missing_contig_omits_adjusted_row(missing: str) -> None:
    g = example_1()
    depth = {"chr1": 50.0, "p1": 400.0}
    del depth[missing]
    res = _plasmids(g, sd.depth_frame(depth))
    assert _adjusted(res).empty
    assert len(res) == 2


# --------------------------------------------------------------------- examples 5, 6: reference


def test_example_5_no_plasmid_no_rows() -> None:
    g = sd.genome([chr_contig()])
    res = _plasmids(g)
    assert res.empty


def test_plasmid_type_switched_off_gives_no_rows() -> None:
    from giap.config.model import GiapConfig

    cfg = GiapConfig.model_validate({"stats": {"bootstrap": B},
                                     "context": {"region_types": ["prophage", "bgc"]}})  # fmt: skip
    assert _plasmids(example_1(), cfg=cfg).empty


def test_reference_is_every_chromosome_replicon() -> None:
    c1 = chr_contig("chr1", n=90)
    c2 = chr_contig("chr2", n=40, base=12.0, kind="chromosome")
    ctg = chr_contig("ctg1", n=30, base=30.0, kind="unknown")  # not a chromosome: excluded
    g = sd.genome([c1, c2, ctg, plasmid_contig()], plasmids=["p1"])
    hl = sd.unadjusted(_plasmids(g), "plasmid_1", "hl_shift")
    want = sd.plasmid_hl(g.genes, "p1", ["chr1", "chr2"])
    assert want != pytest.approx(sd.plasmid_hl(g.genes, "p1", ["chr1", "chr2", "ctg1"]), abs=1e-3)
    assert hl["effect"] == pytest.approx(want, abs=1e-12)
    assert "no_reference_replicon" not in sd.flags(hl["flags"])


def test_example_6_no_chromosome_falls_back_with_flag() -> None:
    """Decision 3 (supersedes example 6's 'hl_shift NA'): without a ``kind = chromosome`` replicon
    the reference is every replicon that is neither ``kind = plasmid`` nor a plasmid region's
    contig; the row is flagged ``no_reference_replicon``."""
    ctg = chr_contig("ctg1", kind="unknown")
    other = plasmid_contig("p2", n=30, base=40.0)  # kind plasmid, no plasmid region: excluded
    pl = plasmid_contig("p1", kind="unknown")  # plasmid by its region (geNOMAD), kind unknown
    g = sd.genome([ctg, other, pl], plasmids=["p1"])
    ref, fallback = sd.reference_contigs(g.replicons, ["p1"])
    assert ref == ["ctg1"] and fallback
    res = _plasmids(g)
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    assert hl["effect"] == pytest.approx(sd.plasmid_hl(g.genes, "p1", ["ctg1"]), abs=1e-12)
    assert "no_reference_replicon" in sd.flags(hl["flags"])
    _assert_no_test(hl)


def test_no_reference_gene_gives_na_hl_with_flag() -> None:
    """Only the plasmid replicon: HL NA, flag ``no_reference_replicon``; the plasmid's own O - E is
    still reported."""
    g = sd.genome([plasmid_contig("p1", kind="unknown")], plasmids=["p1"])
    res = _plasmids(g)
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    assert sd.na(hl["effect"]) and sd.na(hl["hl_shift_log2"])
    assert "no_reference_replicon" in sd.flags(hl["flags"])
    oe = sd.unadjusted(res, "plasmid_1", "o_minus_e")
    assert float(oe["effect"]) == pytest.approx(sd.detection_counts(g.genes, ["p1"])["o_minus_e"])


def test_fallback_reference_without_detected_gene_gives_na_hl() -> None:
    n = 50
    ctg = Contig("ctg1", [0] * n, [0.3] * n, kind="unknown", tu=_pairs("ctg1", n))
    g = sd.genome([ctg, plasmid_contig("p1", kind="unknown")], plasmids=["p1"])
    hl = sd.unadjusted(_plasmids(g), "plasmid_1", "hl_shift")
    assert sd.na(hl["effect"])
    assert "no_reference_replicon" in sd.flags(hl["flags"])


def test_fallback_depth_uses_the_same_reference() -> None:
    ctg = chr_contig("ctg1", kind="unknown")
    g = sd.genome([ctg, plasmid_contig("p1", kind="unknown")], plasmids=["p1"])
    res = _plasmids(g, sd.depth_frame({"ctg1": 25.0, "p1": 100.0}))
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    adj = _adjusted(res)
    assert len(adj) == 1
    assert adj.iloc[0]["effect"] == pytest.approx(hl["effect"] - 2.0, abs=1e-12)


# --------------------------------------------------------------------------- decision 5 edges


def test_plasmid_without_eligible_gene_all_na_small_region() -> None:
    m = 6
    pl = Contig("p1", [0] * m, [0.5] * m, kind="plasmid", partial=[True] * m, circular=True)
    g = sd.genome([chr_contig(), pl], plasmids=["p1"])
    res = _plasmids(g)
    assert set(res["stat"]) == {"o_minus_e", "hl_shift"}
    for _, row in res.iterrows():
        assert int(row["n_eligible"]) == 0 and int(row["n_genes"]) == m
        for c in ("observed", "expected", "o_minus_e", "o_over_e", "hl_shift_log2", "effect",
                  "ci_low", "ci_high"):  # fmt: skip
            assert sd.na(row[c]), (row["stat"], c, row[c])
        assert "small_region" in sd.flags(row["flags"])
        _assert_no_test(row)


def test_o_over_e_na_when_expected_zero() -> None:
    m = 10
    pl = Contig("p1", [0] * m, [0.0] * m, kind="plasmid", tu=_pairs("p1", m))
    g = sd.genome([chr_contig(), pl], plasmids=["p1"])
    oe = sd.unadjusted(_plasmids(g), "plasmid_1", "o_minus_e")
    assert float(oe["expected"]) == 0.0 and float(oe["observed"]) == 0.0
    assert sd.na(oe["o_over_e"])


def test_plasmid_without_detected_gene_has_na_hl_and_real_o_minus_e() -> None:
    m = 10
    pl = Contig("p1", [0] * m, [0.4] * m, kind="plasmid", tu=_pairs("p1", m))
    g = sd.genome([chr_contig(), pl], plasmids=["p1"])
    res = _plasmids(g)
    assert sd.na(sd.unadjusted(res, "plasmid_1", "hl_shift")["effect"])
    assert float(sd.unadjusted(res, "plasmid_1", "o_minus_e")["effect"]) == pytest.approx(-4.0)


def test_no_adjusted_row_when_hl_is_na() -> None:
    """Lead decision 13: HL NA (no detected plasmid gene) -> no depth-adjusted row."""
    m = 10
    pl = Contig("p1", [0] * m, [0.4] * m, kind="plasmid", tu=_pairs("p1", m))
    g = sd.genome([chr_contig(), pl], plasmids=["p1"])
    res = _plasmids(g, sd.depth_frame({"chr1": 50.0, "p1": 400.0}))
    assert _adjusted(res).empty
    assert len(res) == 2


def _small() -> int:
    from giap.context import regions as R

    return int(R.SMALL_REGION)


def test_small_plasmid_flagged_with_na_statistics() -> None:
    """Lead decision 13: fewer eligible genes than M4-17's small-region threshold -> flag
    ``small_region`` and NA statistics (the plasmid rows are not computed)."""
    m = _small() - 1
    assert m >= 1
    pl = Contig("p1", [1] * m, [0.5] * m, a=[20.0] * m, kind="plasmid", tu=_pairs("p1", m))
    g = sd.genome([chr_contig(), pl], plasmids=["p1"])
    res = _plasmids(g, sd.depth_frame({"chr1": 50.0, "p1": 400.0}))
    assert _adjusted(res).empty
    assert set(res["stat"]) == {"o_minus_e", "hl_shift"} and len(res) == 2
    for _, row in res.iterrows():
        assert int(row["n_eligible"]) == m
        assert "small_region" in sd.flags(row["flags"])
        for c in ("effect", "ci_low", "ci_high", "o_minus_e", "o_over_e", "hl_shift_log2"):
            assert sd.na(row[c]), (row["stat"], c, row[c])
        _assert_no_test(row)


def test_plasmid_at_threshold_is_not_small() -> None:
    m = _small()
    pl = Contig("p1", [1] * m, [0.5] * m, a=[20.0 + i for i in range(m)], kind="plasmid",
                tu=_pairs("p1", m))  # fmt: skip
    g = sd.genome([chr_contig(), pl], plasmids=["p1"])
    res = _plasmids(g)
    for _, row in res.iterrows():
        assert "small_region" not in sd.flags(row["flags"])
    assert not sd.na(sd.unadjusted(res, "plasmid_1", "hl_shift")["effect"])


# --------------------------------------------------------------------------- decision 7: RNG


def test_plasmid_rows_use_their_own_stream_per_region() -> None:
    """Adding plasmid_1 leaves plasmid_2's CIs unchanged (child stream keyed by region_id)."""
    contigs = [chr_contig(), plasmid_contig("p1"), plasmid_contig("p2", n=40, base=16.0)]
    both = sd.genome(contigs, plasmids=["p1", "p2"])
    only2 = sd.genome(contigs, {"plasmid_2": ("plasmid", both.genes.loc[
        both.genes["contig"] == "p2", "locus_tag"].tolist())})  # fmt: skip
    a = _plasmids(both)
    b = _plasmids(only2)
    assert set(b["region_id"]) == {"plasmid_2"}
    for stat in ("o_minus_e", "hl_shift"):
        ra, rb = sd.unadjusted(a, "plasmid_2", stat), sd.unadjusted(b, "plasmid_2", stat)
        for c in ("effect", "ci_low", "ci_high"):
            assert float(ra[c]) == float(rb[c]), (stat, c)


def test_m4_17_rows_unchanged_when_a_plasmid_is_added() -> None:
    """Decision 7: the M4-17 rows (and their CIs) in the assembled table do not depend on whether
    the genome has a plasmid region."""
    contigs = [chr_contig(), plasmid_contig()]
    gtmp = sd.genome(contigs)
    tags = gtmp.genes.loc[gtmp.genes["contig"] == "chr1", "locus_tag"].tolist()
    mem = {"bgc_1": ("bgc", tags[10:25]), "prophage_1": ("prophage", tags[40:70])}
    with_p = sd.genome(contigs, mem, plasmids=["p1"])
    without = sd.genome(contigs, mem)
    a, b = _table(with_p), _table(without)
    keep = ["region_id", "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "n_perm"]
    ka = a[a["region_type"] != "plasmid"][keep].reset_index(drop=True)
    kb = b[b["region_type"] != "plasmid"][keep].reset_index(drop=True)
    pd.testing.assert_frame_equal(ka, kb, check_exact=True)
    assert (a["region_type"] == "plasmid").any() and not (b["region_type"] == "plasmid").any()


# --------------------------------------------------------------------------- example 3: subsets


def test_example_3_prophage_subset_rows() -> None:
    g = sd.prophage_genome()
    res = _T().subset_rows(g.genes, g.regions, g.region_genes)
    assert set(res["region_id"]) == {"prophage_1"}
    assert set(res["subset"]) == {"structural", "lysis"}
    p = g.genes.set_index("locus_tag")["p_detect"].astype(float)
    rg = g.region_genes.set_index("locus_tag")
    for name, n_genes, obs in (("structural", 16, 8), ("lysis", 4, 0)):
        row = sd.one(res, "prophage_1", "o_minus_e", subset=name)
        tags = rg.index[rg["subset"] == name]
        e = float(p[tags].sum())
        assert int(row["n_genes"]) == n_genes and int(row["n_eligible"]) == n_genes
        assert float(row["observed"]) == obs
        assert float(row["expected"]) == pytest.approx(e, abs=1e-12)
        assert float(row["o_minus_e"]) == pytest.approx(obs - e, abs=1e-12)
        assert float(row["effect"]) == pytest.approx(obs - e, abs=1e-12)
        assert float(row["o_over_e"]) == pytest.approx(obs / e, abs=1e-12)
        assert row["region_type"] == "prophage" and row["contig"] == "chr1"
        assert "descriptive" in sd.flags(row["flags"])
        assert row["method"] == sd.DESCRIPTIVE and row["label"] == sd.DESCRIPTIVE  # decision 13
        assert sd.na(row["ci_low"]) and sd.na(row["ci_high"])
        _assert_no_test(row)


def test_subset_rows_skip_absent_subsets_and_non_prophages() -> None:
    g = sd.prophage_genome()
    res = _T().subset_rows(g.genes, g.regions, g.region_genes)
    assert len(res) == 2
    assert not set(res["subset"]) & {"all", "integrase", "repressor", "cargo"}


def test_subset_row_with_partial_member() -> None:
    n = 40
    part = [i == 0 for i in range(n)]
    c = Contig("chr1", [1] * n, [0.5] * n, circular=False, partial=part, tu=_pairs("chr1", n),
               kind="chromosome")  # fmt: skip
    gtmp = sd.genome([c])
    tags = gtmp.genes["locus_tag"].tolist()
    g = sd.genome([c], {"prophage_1": ("prophage", tags[0:6])},
                  subsets={t: "integrase" for t in tags[0:6]})  # fmt: skip
    row = sd.one(_T().subset_rows(g.genes, g.regions, g.region_genes), "prophage_1", "o_minus_e",
                 subset="integrase")  # fmt: skip
    assert int(row["n_genes"]) == 6 and int(row["n_eligible"]) == 5
    assert float(row["observed"]) == 5.0 and float(row["expected"]) == pytest.approx(2.5)


# --------------------------------------------------------------------------- example 4 and assembly


def test_regions_stats_columns() -> None:
    T = _T()
    assert tuple(T.REGIONS_STATS_COLUMNS) == sd.REGIONS_STATS_COLUMNS
    res = _table(sd.mixed_genome())
    assert tuple(res.columns) == sd.REGIONS_STATS_COLUMNS


def test_example_4_m4_17_rows_carried_unchanged() -> None:
    g = sd.mixed_genome()
    tests = _tests(g)
    res = _table(g)
    for stat in ("o_minus_e", "hl_shift"):
        src = rd.pick(tests, "bgc_1", stat)
        row = sd.one(res, "bgc_1", stat)
        assert row["subset"] == "all"
        assert row["null_model"] == "circular_shift"
        assert int(row["n_perm"]) == int(src["n_perm"])
        assert row["method"] == src["method"]
        for c in ("effect", "ci_low", "ci_high", "p_value", "q_value"):
            assert float(row[c]) == float(src[c]) or (sd.na(row[c]) and sd.na(src[c])), c
        assert sd.flags(row["flags"]) == rd.flag_set(src["flags"])
    assert len(sd.rows(res, "bgc_1")) == 2


def test_test_rows_self_contained_detection_columns() -> None:
    g = sd.mixed_genome()
    tests = _tests(g)
    res = _table(g)
    for rid in ("prophage_1", "bgc_1", "crispr_1", "amr_1"):
        oe_src = rd.pick(tests, rid, "o_minus_e")
        hl_src = rd.pick(tests, rid, "hl_shift")
        for stat in ("o_minus_e", "hl_shift"):
            row = sd.one(res, rid, stat, subset="all")
            assert float(row["observed"]) == float(oe_src["observed"]), (rid, stat)
            assert float(row["expected"]) == pytest.approx(float(oe_src["expected"]), abs=1e-12)
            assert float(row["o_minus_e"]) == pytest.approx(float(oe_src["effect"]), abs=1e-12)
            e = float(oe_src["expected"])
            assert float(row["o_over_e"]) == pytest.approx(float(oe_src["observed"]) / e)
            h = float(hl_src["effect"])
            assert (sd.na(row["hl_shift_log2"]) and math.isnan(h)) or float(
                row["hl_shift_log2"]
            ) == pytest.approx(h, abs=1e-12)


def test_no_bh_recompute() -> None:
    """Decision 6: ``regions_stats`` copies q-values from the M4-17 rows, never recomputes BH."""
    T = _T()
    g = sd.mixed_genome()
    tests = _tests(g)
    tests = tests.copy()
    tests["q_value"] = np.where(tests["p_value"].notna(), 0.123, np.nan)
    res = T.regions_stats(tests, _plasmids(g), T.subset_rows(g.genes, g.regions, g.region_genes),
                          g.regions)  # fmt: skip
    t = res[res["null_model"] == "circular_shift"]
    q = t["q_value"].to_numpy(float)
    assert np.all((q == 0.123) | np.isnan(q))
    assert (q == 0.123).any()


def test_coordinates_from_c7_and_amr_na() -> None:
    g = sd.mixed_genome()
    res = _table(g)
    c7 = g.regions.set_index("region_id")
    for rid in ("prophage_1", "prophage_2", "bgc_1", "crispr_1", "plasmid_1"):
        for _, row in sd.rows(res, rid).iterrows():
            assert int(row["start"]) == int(c7.loc[rid, "start"])
            assert int(row["end"]) == int(c7.loc[rid, "end"])
            assert row["contig"] == c7.loc[rid, "contig"]
    for _, row in sd.rows(res, "amr_1").iterrows():
        assert sd.na(row["start"]) and sd.na(row["end"])


def test_sort_order() -> None:
    """region_type order of ``context.region_types``, then contig, start, region_id, subset
    (``all`` first)."""
    res = _table(sd.mixed_genome())
    types = list(dict.fromkeys(res["region_type"]))
    assert types == [t for t in sd.TYPE_ORDER if t in types]
    assert list(dict.fromkeys(res["region_id"])) == [
        "prophage_2", "prophage_1", "plasmid_1", "bgc_1", "crispr_1", "amr_1"]  # fmt: skip
    p1 = res[res["region_id"] == "prophage_1"]
    subs = list(p1["subset"])
    assert subs[0] == "all" and subs[1] == "all"
    assert set(subs[2:]) == {"structural", "cargo"}
    for rid, grp in res.groupby("region_id", sort=False):  # rows of a region are contiguous
        idx = grp.index.to_numpy()
        assert (np.diff(idx) == 1).all(), rid
    sd.check_row_order(res.reset_index(drop=True))  # lead decision 11


def test_row_order_with_depth_row() -> None:
    """Decision 11 with the adjusted plasmid row: stat enum order within the plasmid region, and
    the plasmid region still between the prophages and bgc_1."""
    res = _table(sd.mixed_genome(), sd.depth_frame({"chr1": 40.0, "p1": 160.0}))
    res = res.reset_index(drop=True)
    sd.check_row_order(res)
    pl = res[res["region_id"] == "plasmid_1"]
    order = sd.stat_order()
    assert list(pl["stat"]) == sorted(pl["stat"], key=order.index)
    assert len(pl) == 3


def test_prophage_subsets_in_table_have_prophage_coordinates() -> None:
    g = sd.mixed_genome()
    res = _table(g)
    c7 = g.regions.set_index("region_id")
    sub = res[(res["region_id"] == "prophage_1") & (res["subset"] != "all")]
    assert len(sub) == 2
    assert (sub["start"].astype(int) == int(c7.loc["prophage_1", "start"])).all()
    assert (sub["stat"] == "o_minus_e").all()


def test_every_row_null_model_none_iff_untested() -> None:
    """Decision 6 / acceptance: no row has ``p_value`` with ``null_model = none``; no row has a
    ``q_value`` without ``p_value``; ``null_model = none`` rows have NA p, q and ``n_perm = 0``."""
    g = sd.mixed_genome()
    res = _table(g, sd.depth_frame({"chr1": 40.0, "p1": 160.0}))
    for _, row in res.iterrows():
        if row["null_model"] == "none":
            _assert_no_test(row)
        if not sd.na(row["p_value"]):
            assert row["null_model"] != "none"
        if not sd.na(row["q_value"]):
            assert not sd.na(row["p_value"])
    pl = res[res["region_type"] == "plasmid"]
    assert len(pl) == 3 and (pl["null_model"] == "none").all()
    assert (res["subset"] != "all").sum() == 2


def test_empty_genome_gives_empty_table() -> None:
    T = _T()
    g = sd.genome([chr_contig()])
    res = T.regions_stats(_tests(g), _plasmids(g), T.subset_rows(g.genes, g.regions,
                          g.region_genes), g.regions)  # fmt: skip
    assert res.empty
    assert tuple(res.columns) == sd.REGIONS_STATS_COLUMNS


# --------------------------------------------------------------------------- example 9


def test_example_9_same_seed_same_table(tmp_path: Path) -> None:
    from giap.io.regions_tsv import write_regions_stats

    g = sd.mixed_genome()
    depth = sd.depth_frame({"chr1": 40.0, "p1": 160.0})
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    write_regions_stats(a, _table(g, depth, seed=5))
    write_regions_stats(b, _table(g, depth, seed=5))
    assert a.read_bytes() == b.read_bytes()
    c = tmp_path / "c.tsv"
    write_regions_stats(c, _table(g, depth, seed=6))
    assert c.read_bytes() != a.read_bytes()  # the CIs do use the seed
