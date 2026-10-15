"""M4-16 lead decisions 4-9 after the statistician's review (PR #87, changes requested).

Spec: docs/specs/M4-16-detection-summary.md ("Lead decisions after the statistician's review"),
GIAP_SPEC.md §7.5.1 (``id_mixups``, genome-inconsistent peptides), §7.3.5(a), §8.2;
docs/methods/detection-summary.md (statistician sign-off items 1-3).

Resolved ambiguities (pinned here):
- Decision 4: ``ci_few_units`` on ``o_minus_e`` iff ``n_tu < 10`` or min(E, n_genes - E) < 10,
  with E = sum of ``p_detect`` (expected, not observed, detections); strict ``<`` (E = 10 is not
  flagged). The flag never appears on any other basic row.
- Decision 5: ``genome_inconsistent_rate`` = share of distinct target (non-decoy),
  non-contaminant (C9 ``is_contaminant``) ``search_db = other`` peptides without an exact C11
  placement; ``n`` = that denominator. NA when there are none (e.g. all peptides ``same``).
- Decision 6: ``id_mixups`` row: value NA, flag ``requires_variants``. New QC row
  ``engine_protein_mismatch_rate``: distinct C11 peptides flagged ``engine_protein_mismatch`` over
  distinct placed peptides (pinned only where every placed peptide is a same-db target, so the
  denominator is unambiguous); NA without C11 rows.
- Decision 7: ValueError naming the genes for an eligible gene with NaN ``p_detect`` or without
  ``tu_id``; partial genes may lack both.
- Decision 8: COG ``-`` and ``""`` count in ``none``; a ``detected`` column of strings raises
  ValueError (``"False"`` would otherwise read as true).
- Decision 9 (surviving mutants): groups are processed in sorted (group_type, group) order, each
  drawing ``rng.integers(0, n_tu, size=(B, n_tu))`` (M4-04 ``tu_block_bootstrap``); the oracle
  below uses TU ids whose sorted order equals their order of appearance, so it does not depend on
  how the implementation orders TUs inside a group.
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
import detection_summary_fixture as fx
from detection_summary_fixture import gene, metric, qc_hit, qc_peptide

SEED = 20261015


def _summary() -> Any:
    from giap.context import summary

    return summary


def _basic(genes: pd.DataFrame, b: int = 200, seed: int = SEED) -> pd.DataFrame:
    from giap.rng import for_step

    return _summary().basic_rows(genes, b, for_step("context_stats", seed))  # type: ignore[no-any-return]


def _flags(out: pd.DataFrame, gt: str = "replicon", g: str = "chr1") -> set[str]:
    return fx.flag_set(metric(out, gt, g, "o_minus_e")["flags"])


def _singletons(n: int, p: float, det: Any = None) -> pd.DataFrame:
    d = det if det is not None else [i % 2 == 0 for i in range(n)]
    return fx.genes_frame(gene(f"G{i:03d}", f"TU{i:03d}", bool(d[i]), p) for i in range(n))


# --------------------------------------------------------------------------- decision 4


@pytest.mark.parametrize(
    ("n", "p", "flagged"),
    [
        (20, 0.4, True),  # E = 8 < 10 although n_tu = 20
        (20, 0.5, False),  # E = n - E = 10: strict <
        (20, 0.6, True),  # n - E = 8 < 10
        (30, 0.2, True),  # E = 6
        (30, 0.5, False),  # E = 15, n - E = 15
        (100, 0.05, True),  # sORF-like: E = 5, n_tu = 100
        (100, 0.95, True),  # near saturation: n - E = 5
        (100, 0.15, False),  # E = 15
    ],
)
def test_flag_rule_min_expected(n: int, p: float, flagged: bool) -> None:
    out = _basic(_singletons(n, p))
    assert metric(out, "replicon", "chr1", "n_tu")["value"] == n >= 10
    assert (fx.FEW_UNITS_FLAG in _flags(out)) is flagged, (n, p)


def test_flag_rule_uses_expected_not_observed() -> None:
    # 20 TUs, p = 0.5 (E = 10), nothing detected: O = 0 but E and n - E are 10 -> no flag
    out = _basic(_singletons(20, 0.5, det=[False] * 20))
    assert metric(out, "replicon", "chr1", "n_detected")["value"] == 0
    assert fx.FEW_UNITS_FLAG not in _flags(out)
    # 20 TUs, p = 0.3 (E = 6), every gene detected: O = 20 but E < 10 -> flag
    out = _basic(_singletons(20, 0.3, det=[True] * 20))
    assert fx.FEW_UNITS_FLAG in _flags(out)


def test_flag_rule_per_group() -> None:
    """A low-E sORF group is flagged while the cds group of the same frame is not."""
    rows = [gene(f"C{i:03d}", f"TC{i:03d}", i % 2 == 0, 0.5) for i in range(40)]
    rows += [gene(f"S{i:03d}", f"TS{i:03d}", False, 0.03, feature_type="sorf") for i in range(60)]
    out = _basic(fx.genes_frame(rows))
    assert metric(out, "feature_class", "sorf", "expected")["value"] == pytest.approx(1.8)
    assert fx.FEW_UNITS_FLAG in _flags(out, "feature_class", "sorf")
    assert fx.FEW_UNITS_FLAG not in _flags(out, "feature_class", "cds")
    # chr1: E = 21.8, n - E = 78.2, 100 TUs
    assert fx.FEW_UNITS_FLAG not in _flags(out)


@pytest.mark.parametrize("n", [5, 20])
def test_flag_only_on_o_minus_e(n: int) -> None:
    """Mutant ``flag_on_all``: no basic row other than ``o_minus_e`` carries ``ci_few_units``."""
    rows = [gene(f"G{i:03d}", f"TU{i:03d}", i % 3 == 0, 0.2, cog_category="K",
                 hypothetical=True) for i in range(n)]  # fmt: skip
    out = _basic(fx.genes_frame(rows))
    for r in out.itertuples(index=False):
        has = fx.FEW_UNITS_FLAG in fx.flag_set(r.flags)
        assert has is (r.metric == "o_minus_e"), (r.group_type, r.group, r.metric, r.flags)
        if r.metric != "o_minus_e":
            assert fx.flag_set(r.flags) == set(), (r.group_type, r.group, r.metric, r.flags)


# --------------------------------------------------------------------------- decision 9: order


def _oracle_cis(
    genes: pd.DataFrame, b: int, seed: int
) -> dict[tuple[str, str], tuple[float, float]]:
    from giap.rng import for_step

    rng = for_step("context_stats", seed)
    out = {}
    for key in sorted(fx.group_members(genes)):
        mask = fx.group_members(genes)[key]
        sums = fx.tu_sums(genes, mask)
        v = np.array([sums[t] for t in sorted(sums)])
        k = v.size
        draws = rng.integers(0, k, size=(b, k))
        boot = v[draws].sum(axis=1)
        lo, hi = np.quantile(boot, [0.025, 0.975], method="linear")
        out[key] = (float(lo), float(hi))
    return out


def _order_frame() -> pd.DataFrame:
    rows = [gene(f"G{i:03d}", f"TU{i // 2:03d}", i % 3 != 0, 0.05 * (i % 17),
                 contig="chr1" if i < 30 else "p1",
                 feature_type=["cds", "sorf", "cds", "pseudogene"][i % 4],
                 hypothetical=i % 5 == 0,
                 cog_category=["K", "T", None, "KT", "J"][i % 5]) for i in range(44)]  # fmt: skip
    return fx.genes_frame(rows)


def test_groups_draw_in_sorted_order() -> None:
    """Mutant ``unsorted_groups``: every group's CI equals the oracle that draws groups in sorted
    (group_type, group) order from one ``for_step("context_stats", seed)`` stream."""
    g = _order_frame()
    b = 300
    out = _basic(g, b=b, seed=7)
    want = _oracle_cis(g, b, 7)
    assert len(want) >= 8
    for key, (lo, hi) in want.items():
        r = metric(out, *key, "o_minus_e")
        assert float(r["ci_low"]) == pytest.approx(lo, abs=1e-9), key
        assert float(r["ci_high"]) == pytest.approx(hi, abs=1e-9), key


def test_reversed_group_order_would_change_cis() -> None:
    """Sanity of the oracle: drawing the groups in reverse order gives different CIs here."""
    from giap.rng import for_step

    g = _order_frame()
    groups = fx.group_members(g)
    fwd = _oracle_cis(g, 300, 7)
    rng = for_step("context_stats", 7)
    rev = {}
    for key in sorted(groups, reverse=True):
        sums = fx.tu_sums(g, groups[key])
        v = np.array([sums[t] for t in sorted(sums)])
        boot = v[rng.integers(0, v.size, size=(300, v.size))].sum(axis=1)
        rev[key] = tuple(np.quantile(boot, [0.025, 0.975], method="linear"))
    assert any(not np.allclose(fwd[k], rev[k]) for k in fwd)


# --------------------------------------------------------------------------- decision 7


def test_nan_p_detect_on_eligible_gene_raises_naming_it() -> None:
    rows = [gene(f"G{i}", f"TU{i}", True, 0.5) for i in range(4)]
    rows.append(gene("BADP1", "TU9", False, None))
    rows.append(gene("BADP2", "TU8", True, None))
    with pytest.raises(ValueError, match="BADP1") as exc:
        _basic(fx.genes_frame(rows))
    assert "BADP2" in str(exc.value)


def test_missing_tu_on_eligible_gene_raises_naming_it() -> None:
    for missing in (None, math.nan):
        rows = [gene(f"G{i}", f"TU{i}", True, 0.5) for i in range(4)]
        bad = gene("NOTU1", "x", False, 0.2)
        bad["tu_id"] = missing
        rows.append(bad)
        with pytest.raises(ValueError, match="NOTU1"):
            _basic(fx.genes_frame(rows))


def test_partial_gene_without_p_or_tu_is_fine() -> None:
    rows = [gene(f"G{i}", f"TU{i}", True, 0.5) for i in range(4)]
    part = gene("P1", "x", False, None, status="partial_excluded")
    part["tu_id"] = None
    rows.append(part)
    out = _basic(fx.genes_frame(rows))
    assert metric(out, "replicon", "chr1", "n_genes")["value"] == 4


# --------------------------------------------------------------------------- decision 8


def test_cog_dash_and_blank_map_to_none() -> None:
    rows = [
        gene("G1", "A", True, 0.5, cog_category="-"),
        gene("G2", "B", False, 0.5, cog_category=""),
        gene("G3", "C", True, 0.5),
        gene("G4", "D", True, 0.5, cog_category="K"),
    ]
    out = _basic(fx.genes_frame(rows), b=50)
    cog = set(out.loc[out["group_type"] == "cog_category", "group"].astype(str))
    assert cog == {"K", "none"}
    assert metric(out, "cog_category", "none", "n_genes")["value"] == 3
    assert metric(out, "cog_category", "none", "n_detected")["value"] == 2


@pytest.mark.parametrize(
    "values", [["True", "False", "False", "True"], ["true", "false", "no", "yes"]]
)
def test_string_detected_raises(values: list[str]) -> None:
    g = fx.genes_frame(gene(f"G{i}", f"TU{i}", True, 0.5) for i in range(4))
    g["detected"] = pd.Series(values, dtype=object)
    with pytest.raises(ValueError):
        _basic(g, b=20)


def test_bool_detected_accepted() -> None:
    g = fx.genes_frame(gene(f"G{i}", f"TU{i}", i < 2, 0.5) for i in range(4))
    assert g["detected"].dtype == bool
    assert metric(_basic(g, b=20), "replicon", "chr1", "n_detected")["value"] == 2


# --------------------------------------------------------------------------- decisions 5-6 (QC)


def _qc(peps: pd.DataFrame | None, hits: pd.DataFrame) -> pd.DataFrame:
    return _summary().qc_rows(peps, hits, None)  # type: ignore[no-any-return]


def _row(df: pd.DataFrame, name: str) -> pd.Series:
    return fx.metric(df, "run", "all", name)


def test_genome_inconsistent_rate_over_other_db_excluding_contaminants() -> None:
    names = [f"OTH{i:03d}K" for i in range(100)]
    hits = [qc_hit(p) for p in names[:97]]
    hits += [qc_hit(names[97], route="tolerant"), qc_hit(names[98], route="tolerant")]
    peps = [qc_peptide(p, search_db="other") for p in names]  # names[99]: no C11 row
    # excluded: contaminants (unplaced; would raise the rate), decoys, same-db peptides (unplaced)
    peps += [dict(qc_peptide(f"CON{i}K", search_db="other"), is_contaminant=True)
             for i in range(10)]  # fmt: skip
    peps += [qc_peptide(f"DEC{i}K", search_db="other", is_decoy=True) for i in range(10)]
    peps += [qc_peptide(f"SAM{i}K") for i in range(10)]
    out = _qc(fx.peptides_frame(peps), fx.hits_frame(hits))
    r = _row(out, "genome_inconsistent_rate")
    assert float(r["value"]) == pytest.approx(0.03)
    assert int(r["n"]) == 100


def test_genome_inconsistent_rate_na_for_same_db_and_all_contaminant() -> None:
    same = [qc_peptide(f"SAM{i}K") for i in range(10)]  # none placed: old rule gave 1.0
    out = _qc(fx.peptides_frame(same), fx.hits_frame([qc_hit("ZZZZK")]))
    assert math.isnan(float(_row(out, "genome_inconsistent_rate")["value"]))
    cont = [dict(qc_peptide(f"CON{i}K", search_db="other"), is_contaminant=True) for i in range(3)]
    out = _qc(fx.peptides_frame(cont), fx.hits_frame([qc_hit("ZZZZK")]))
    assert math.isnan(float(_row(out, "genome_inconsistent_rate")["value"]))


def test_id_mixups_na_with_requires_variants() -> None:
    hits = fx.hits_frame([qc_hit("AAAAAAAK", flags=["engine_protein_mismatch"]),
                          qc_hit("CCCCCCCK")])  # fmt: skip
    peps = fx.peptides_frame([qc_peptide("AAAAAAAK"), qc_peptide("CCCCCCCK")])
    out = _qc(peps, hits)
    r = _row(out, "id_mixups")
    assert math.isnan(float(r["value"]))
    assert fx.flag_set(r["flags"]) == {fx.REQUIRES_VARIANTS_FLAG}
    others = out[out["metric"] != "id_mixups"]
    assert all(fx.REQUIRES_VARIANTS_FLAG not in fx.flag_set(f) for f in others["flags"])


def test_engine_protein_mismatch_rate() -> None:
    hits = fx.hits_frame([
        qc_hit("AAAAAAAK", flags=["engine_protein_mismatch"]),
        qc_hit("AAAAAAAK", flags=["engine_protein_mismatch"]),  # same peptide, second locus
        qc_hit("CCCCCCCK", route="tolerant", flags=["engine_protein_mismatch"]),
        qc_hit("DDDDDDDK", flags=["same_db"]),
        qc_hit("EEEEEEEK"),
    ])  # fmt: skip
    peps = fx.peptides_frame(qc_peptide(p) for p in sorted(set(hits["peptide"])))
    out = _qc(peps, hits)
    assert float(_row(out, "engine_protein_mismatch_rate")["value"]) == pytest.approx(0.5)
    none = _qc(peps, fx.hits_frame([]))
    assert math.isnan(float(_row(none, "engine_protein_mismatch_rate")["value"]))


def test_qc_metric_set_after_review() -> None:
    out = _qc(fx.peptides_frame([qc_peptide("AAAAAAAK")]), fx.hits_frame([qc_hit("AAAAAAAK")]))
    assert {"id_mixups", "engine_protein_mismatch_rate", "genome_inconsistent_rate"} <= set(
        out["metric"].astype(str)
    )
