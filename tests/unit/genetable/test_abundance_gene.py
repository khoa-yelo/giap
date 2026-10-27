"""M4-14 gene-level abundance (docs/specs/M4-14-abundance.md): examples 1-10, edge cases, routes.

Spec: GIPPA_SPEC.md §2.1 decisions 2-3, §7.3.4 (iBAQ, Top3, protein-level route, transform, missing
values), §7.2.6 (two proteins on one gene, cross-strain proteins), §9.3 abundance columns, §9.4
flags, §5.2 ``abundance`` config. M4-11 lead decisions: only ``abundance_eligible`` accepted rows
(unique, not an N/C-terminal extension, exact or tolerant under ``related``) are quantified.
Every scenario is cross-checked against tests/helpers/abundance_oracle.py where it applies.

Resolved ambiguities (pinned here; questions 1-3 of the spec use the proposed defaults):
- Top3 is the linear mean of the (up to) 3 largest intensities, then log2 (Q1).
- The fill statistic uses observed values of eligible (non-partial) detected genes only (Q2); a
  partial gene is excluded even if a stub marks it detected.
- Undetected genes keep an observed value (Q3).
- ``detection`` is an M4-12 ``Detection``-like object; only ``.table`` (locus_tag, detected)
  is read.
- PR #82 lead decision 1: C10 has NO ``scale`` column (M4-06 decision on PR #51); the measurement
  scale is ``abundance.scale`` when not ``auto``, else the preflight scale passed as
  ``gene_abundance(..., c10_scale=...)`` (``protein_abundance(..., scale)``). The tests here keep
  ``abundance.scale = auto`` and pass ``c10_scale``; precedence is in test_abundance_review.py.
- PR #82 lead decision 2: the protein -> gene map has only protein_id, gene, method;
  ``cross_strain_protein`` comes from the C12 protein-hit flags (``protein_hits=`` keyword of
  ``gene_abundance`` / ``protein_abundance``) of a protein mapped by ``protein_hit``.
- On the protein route ``top3_log2`` is NA; ``abundance_raw`` keeps the input value (log2 input
  stays log2). With source ``none`` ``abundance_imputed`` is NA or False (no abundance at all).
- ``n_quant_peptides`` counts distinct quantified peptides; genes without one are either absent
  from ``peptide_abundance`` or carry 0 and NA values.
- The ``no_observed_abundance`` warning is a Python warning or a log record carrying that code.
- The declared abundance type is ``abundance.type`` (backlog: "record the declared abundance type
  in the input metadata" is wiring for the step, not part of the pure interface tested here).
"""

from __future__ import annotations

import logging
import math
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import abundance_oracle as ao
import evidence_oracle as eo
from abundance_oracle import (
    acc,
    accepted_frame,
    c10,
    covariates_frame,
    detection,
    pmap_frame,
    protein_hits_frame,
)

G, H, K, L, M = "GSIM_00005", "GSIM_00010", "GSIM_00015", "GSIM_00020", "GSIM_00025"
TAGS = [G, H, K, L, M]
TOL = 1e-6


def _ab() -> Any:
    from gippa.genetable import abundance

    return abundance


def cfg(**abundance: Any) -> Any:
    from gippa.config import GippaConfig

    relation = abundance.pop("relation", "related")
    return GippaConfig(genome_relation=relation, abundance=abundance)


def genes_frame(tags: list[str] = TAGS, partial: tuple[str, ...] = ()) -> pd.DataFrame:
    return eo.genes_frame(
        [
            eo.gene(t, 1000 + 1000 * i, 1300 + 1000 * i, partial=t in partial)
            for i, t in enumerate(tags)
        ]
    )


def empty_pmap() -> pd.DataFrame:
    return pmap_frame([])


def run(
    accepted: pd.DataFrame | None = None,
    *,
    n_obs: dict[str, int] | None = None,
    detected: tuple[str, ...] = (),
    partial: tuple[str, ...] = (),
    abundance: pd.DataFrame | None = None,
    pmap: pd.DataFrame | None = None,
    policy: str | None = None,
    tags: list[str] = TAGS,
    c10_scale: str | None = None,
    protein_hits: pd.DataFrame | None = None,
    **conf: Any,
) -> Any:
    ab = _ab()
    genes = genes_frame(tags, partial)
    acc_df = accepted if accepted is not None else accepted_frame([])
    cov = covariates_frame(n_obs or {}, order=tags)
    det = detection(tags, detected, partial)
    return ab.gene_abundance(
        genes,
        acc_df,
        cov,
        det,
        abundance,
        pmap if pmap is not None else empty_pmap(),
        cfg(**conf),
        policy=policy,
        c10_scale=c10_scale,
        protein_hits=protein_hits if protein_hits is not None else protein_hits_frame([]),
    )


def row(result: Any, tag: str) -> dict[str, Any]:
    t = result.table
    sel = t[t["locus_tag"] == tag]
    assert len(sel) == 1, f"{tag}: {len(sel)} rows"
    return sel.iloc[0].to_dict()


def flags(r: dict[str, Any]) -> set[str]:
    cell = r["abundance_flags"]
    if ao.is_na(cell):
        return set()
    items = [str(f) for f in cell]
    assert len(items) == len(set(items)), f"duplicate flags {items}"
    return set(items)


def na(x: Any) -> bool:
    return ao.is_na(x)


def approx(x: Any, y: float) -> bool:
    return not na(x) and abs(float(x) - y) <= TOL


def ex1_rows(gene: str = G) -> list[dict[str, Any]]:
    return [acc(gene, f"PEPTIDE{i}K", v, start=30 * i) for i, v in enumerate((100, 300, 600, 1000))]


# ------------------------------------------------------------------------------ table shape


def test_table_has_one_row_per_c2_gene_in_order_with_9_3_columns() -> None:
    res = run(accepted_frame(ex1_rows()), n_obs={G: 4}, detected=(G,))
    t = res.table
    assert t["locus_tag"].tolist() == TAGS
    for col in ("locus_tag", *ao.ABUNDANCE_COLUMNS, "abundance_flags"):
        assert col in t.columns, col
    assert res.source == "peptide_intensity"


# ------------------------------------------------------------------------------ example 1


def test_example1_ibaq_and_top3_values() -> None:
    res = run(accepted_frame(ex1_rows()), n_obs={G: 4}, detected=(G,))
    r = row(res, G)
    assert approx(r["abundance_raw"], 500.0)
    assert approx(r["abundance_log2"], 8.965784)
    assert approx(r["top3_log2"], 9.306821)
    assert r["abundance_type"] == "ibaq"
    assert r["abundance_scale"] == "molar"
    assert not bool(r["abundance_imputed"])
    assert "top3_lt3" not in flags(r)


def test_example1_method_top3_selects_top3() -> None:
    res = run(accepted_frame(ex1_rows()), n_obs={G: 4}, detected=(G,), method="top3")
    r = row(res, G)
    assert approx(r["abundance_raw"], 1900.0 / 3.0)
    assert approx(r["abundance_log2"], 9.306821)
    assert approx(r["top3_log2"], 9.306821)
    assert r["abundance_type"] == "top3"
    assert r["abundance_scale"] == "molar"


def test_example1_peptide_abundance_direct() -> None:
    ab = _ab()
    out = ab.peptide_abundance(accepted_frame(ex1_rows()), pd.Series({G: 4}))
    sel = out[out["locus_tag"] == G].iloc[0]
    assert sel["ibaq"] == pytest.approx(500.0, abs=TOL)
    assert sel["top3"] == pytest.approx(1900.0 / 3.0, abs=TOL)
    assert int(sel["n_quant_peptides"]) == 4


# ------------------------------------------------------------------------------ example 2


def test_example2_n_observable_zero_and_top3_lt3() -> None:
    rows = [acc(G, "AAAAAAAK", 100), acc(G, "CCCCCCCK", 300, start=30)]
    res = run(accepted_frame(rows), n_obs={G: 0}, detected=(G,))
    r = row(res, G)
    assert approx(r["abundance_raw"], 400.0)
    assert approx(r["abundance_log2"], math.log2(400.0))
    assert approx(r["top3_log2"], math.log2(200.0))
    assert "top3_lt3" in flags(r)


def test_example2_peptide_abundance_direct() -> None:
    ab = _ab()
    rows = [acc(G, "AAAAAAAK", 100), acc(G, "CCCCCCCK", 300, start=30)]
    out = ab.peptide_abundance(accepted_frame(rows), pd.Series({G: 0}))
    sel = out[out["locus_tag"] == G].iloc[0]
    assert sel["ibaq"] == pytest.approx(400.0)
    assert sel["top3"] == pytest.approx(200.0)
    assert int(sel["n_quant_peptides"]) == 2


def test_top3_lt3_only_for_one_or_two_quantified_peptides() -> None:
    rows = [acc(G, "AAAAAAAK", 100)]  # 1 peptide
    rows += [acc(H, f"PEPH{i}AAAK", 10.0 * (i + 1), start=30 * i) for i in range(3)]  # 3
    res = run(accepted_frame(rows), n_obs={G: 3, H: 3}, detected=(G, H, K))
    assert "top3_lt3" in flags(row(res, G))
    assert approx(row(res, G)["top3_log2"], math.log2(100.0))
    assert "top3_lt3" not in flags(row(res, H))
    assert "top3_lt3" not in flags(row(res, K)), "no quantified peptide: no top3_lt3"


# ------------------------------------------------------------------------------ example 3


def test_example3_multi_identical_peptide_not_used() -> None:
    rows = [
        acc(G, "MULTIPEPK", 1e6, uniqueness="multi_identical"),
        acc(G, "UNIQUEPEPK", 100, start=60),
    ]
    res = run(accepted_frame(rows), n_obs={G: 1}, detected=(G,))
    r = row(res, G)
    assert approx(r["abundance_raw"], 100.0)
    assert approx(r["top3_log2"], math.log2(100.0))
    assert "top3_lt3" in flags(r)


@pytest.mark.parametrize("cls", ["nterm_extension", "cterm_extension"])
def test_extension_placements_not_abundance_eligible(cls: str) -> None:
    rows = [acc(G, "EXTPEPTIDEK", 1e6, locus_class=cls), acc(G, "UNIQUEPEPK", 100, start=60)]
    res = run(accepted_frame(rows), n_obs={G: 1}, detected=(G,))
    assert approx(row(res, G)["abundance_raw"], 100.0)


# ------------------------------------------------------------------------------ example 4


def _ex4(relation: str) -> Any:
    from gippa.genetable import evidence

    hits = eo.hits_frame(
        [
            eo.hit("TOLERANTPEPK", G, 1000, route="tolerant", flags=["strain_difference"]),
            eo.hit("EXACTPEPTIDEK", G, 1060),
        ]
    )
    peps = eo.peptides_frame(
        [eo.pep("TOLERANTPEPK", intensity=50.0), eo.pep("EXACTPEPTIDEK", intensity=150.0)]
    )
    c = cfg(relation=relation)
    accepted = evidence.accepted_peptides(hits, peps, c)
    return run(accepted, n_obs={G: 2}, detected=(G,), relation=relation)


def test_example4_tolerant_hit_included_under_related() -> None:
    r = row(_ex4("related"), G)
    assert approx(r["abundance_raw"], (50.0 + 150.0) / 2)
    assert approx(r["top3_log2"], math.log2(100.0))


def test_example4_tolerant_hit_not_accepted_under_same() -> None:
    r = row(_ex4("same"), G)
    assert approx(r["abundance_raw"], 150.0 / 2)
    assert approx(r["top3_log2"], math.log2(150.0))


# ------------------------------------------------------------------------------ example 5


def test_example5_two_proteins_one_gene_takes_larger_log2() -> None:
    ab = c10({"P1": 20.0, "P2": 22.0})
    pm = pmap_frame([("P1", G, "protein_hit"), ("P2", G, "peptide_vote")])
    res = run(None, detected=(G,), abundance=ab, pmap=pm, type="ibaq", c10_scale="log2")
    assert res.source == "protein_table"
    r = row(res, G)
    assert approx(r["abundance_log2"], 22.0)
    assert approx(r["abundance_raw"], 22.0)
    assert "two_proteins_one_gene" in flags(r)
    assert r["abundance_scale"] == "molar"
    assert r["abundance_type"] == "ibaq"
    assert na(r["top3_log2"])


def test_example5_protein_abundance_direct() -> None:
    abm = _ab()
    ab = c10({"P1": 20.0, "P2": 22.0, "P3": 18.0})
    pm = pmap_frame([("P1", G, "protein_hit"), ("P2", G, "protein_hit"), ("P3", H, "locus_tag")])
    out = abm.protein_abundance(ab, pm, "ibaq", "log2")
    g = out[out["locus_tag"] == G].iloc[0]
    h = out[out["locus_tag"] == H].iloc[0]
    assert g["value_log2"] == pytest.approx(22.0)
    assert int(g["n_proteins"]) == 2
    assert "two_proteins_one_gene" in [str(f) for f in g["flags"]]
    assert h["value_log2"] == pytest.approx(18.0)
    assert int(h["n_proteins"]) == 1
    assert "two_proteins_one_gene" not in [str(f) for f in h["flags"]]


def test_protein_route_linear_values_to_log2_and_raw_kept() -> None:
    ab = c10({"P1": 1024.0, "P2": 4096.0})
    pm = pmap_frame([("P1", G, "protein_hit"), ("P2", H, "protein_hit")])
    res = run(None, detected=(G, H), abundance=ab, pmap=pm, type="top3", c10_scale="linear")
    assert approx(row(res, G)["abundance_log2"], 10.0)
    assert approx(row(res, G)["abundance_raw"], 1024.0)
    assert approx(row(res, H)["abundance_log2"], 12.0)
    assert row(res, G)["abundance_scale"] == "molar"
    assert row(res, G)["abundance_type"] == "top3"


def test_protein_route_unmapped_protein_ignored() -> None:
    ab = c10({"P1": 20.0, "PX": 30.0})
    pm = pmap_frame([("P1", G, "protein_hit"), ("PX", None, "protein_hit")])
    res = run(None, detected=(G,), abundance=ab, pmap=pm, type="ibaq", c10_scale="log2")
    assert approx(row(res, G)["abundance_log2"], 20.0)
    for t in (H, K, L, M):
        assert na(row(res, t)["abundance_log2"])


def test_cross_strain_protein_carries_abundance_and_flag() -> None:
    ab = c10({"Q9XYZ1": 21.0, "P2": 19.0})
    pm = pmap_frame([("Q9XYZ1", H, "protein_hit"), ("P2", G, "protein_hit")])
    hits = protein_hits_frame([("Q9XYZ1", H, ["cross_strain_protein"]), ("P2", G, [])])
    res = run(None, detected=(G, H), abundance=ab, pmap=pm, type="ibaq", c10_scale="log2",
              protein_hits=hits)  # fmt: skip
    assert approx(row(res, H)["abundance_log2"], 21.0)
    assert "cross_strain_protein" in flags(row(res, H))
    assert "cross_strain_protein" not in flags(row(res, G))


# ------------------------------------------------------------------------------ example 6


@pytest.mark.parametrize(
    ("declared", "scale"),
    [
        ("maxlfq", "relative"),
        ("intensity", "relative"),
        ("unknown", "relative"),
        ("ibaq", "molar"),
        ("top3", "molar"),
    ],
)
def test_example6_protein_route_scale_from_declared_type(declared: str, scale: str) -> None:
    ab = c10({"P1": 20.0})
    pm = pmap_frame([("P1", G, "protein_hit")])
    res = run(None, detected=(G,), abundance=ab, pmap=pm, type=declared, c10_scale="log2")
    r = row(res, G)
    assert r["abundance_type"] == declared
    assert r["abundance_scale"] == scale
    for t in TAGS:
        assert row(res, t)["abundance_type"] == declared
        assert row(res, t)["abundance_scale"] == scale


# ------------------------------------------------------------------------------ example 7


EX7_FILL = {"median": 13.0, "mean": 14.0, "min": 10.0, "zero": 9.0, "exclude": None}


def _ex7_log2() -> pd.Series:
    # G..L detected; M detected and missing; "U" undetected with a large value (not in the stat)
    return pd.Series(
        [10.0, 12.0, 14.0, 20.0, math.nan, 100.0, math.nan],
        index=[G, H, K, L, M, "U1", "U2"],
    )


def _ex7_detected() -> pd.Series:
    return pd.Series([True, True, True, True, True, False, False], index=_ex7_log2().index)


@pytest.mark.parametrize("policy", list(EX7_FILL))
def test_example7_impute_policies(policy: str) -> None:
    ab = _ab()
    filled, mask, fill = ab.impute(_ex7_log2(), _ex7_detected(), policy)
    exp = EX7_FILL[policy]
    filled = pd.Series(filled, index=_ex7_log2().index)
    mask = pd.Series(mask, index=_ex7_log2().index).astype(bool)
    if exp is None:
        assert fill is None
        assert na(filled[M])
        assert not mask.any()
    else:
        assert fill == pytest.approx(exp, abs=TOL)
        assert filled[M] == pytest.approx(exp, abs=TOL)
        assert mask.tolist() == [False, False, False, False, True, False, False]
    assert na(filled["U2"]), "undetected genes are never imputed"
    assert filled["U1"] == 100.0
    for t, v in ((G, 10.0), (H, 12.0), (K, 14.0), (L, 20.0)):
        assert filled[t] == v
    o_filled, o_mask, o_fill = ao.impute(_ex7_log2().to_dict(), _ex7_detected().to_dict(), policy)
    assert (o_fill is None) == (fill is None)
    assert mask.to_dict() == o_mask
    for t, v in o_filled.items():
        assert (na(v) and na(filled[t])) or filled[t] == pytest.approx(v, abs=TOL)


def _ex7_gene_run(policy: str | None, **conf: Any) -> Any:
    # single-peptide genes with n_observable 1: iBAQ = intensity
    rows = [acc(t, f"PEP{t}K", 2.0**v) for t, v in ((G, 10), (H, 12), (K, 14), (L, 20))]
    return run(accepted_frame(rows), n_obs={t: 1 for t in TAGS}, detected=tuple(TAGS),
               policy=policy, **conf)  # fmt: skip


@pytest.mark.parametrize("policy", list(EX7_FILL))
def test_example7_gene_abundance_policy_override(policy: str) -> None:
    res = _ex7_gene_run(policy, missing="median")
    r = row(res, M)
    exp = EX7_FILL[policy]
    if exp is None:
        assert na(r["abundance_log2"])
        assert not bool(r["abundance_imputed"])
        assert "imputed_abundance" not in flags(r)
        assert res.fill_value is None
    else:
        assert approx(r["abundance_log2"], exp)
        assert bool(r["abundance_imputed"])
        assert "imputed_abundance" in flags(r)
        assert res.fill_value == pytest.approx(exp, abs=TOL)
    assert set(res.table["imputation_policy"]) == {policy}
    for t in (G, H, K, L):
        assert not bool(row(res, t)["abundance_imputed"])
        assert "imputed_abundance" not in flags(row(res, t))


@pytest.mark.parametrize("policy", list(EX7_FILL))
def test_example7_config_missing_used_without_override(policy: str) -> None:
    res = _ex7_gene_run(None, missing=policy)
    assert set(res.table["imputation_policy"]) == {policy}
    exp = EX7_FILL[policy]
    r = row(res, M)
    assert na(r["abundance_log2"]) if exp is None else approx(r["abundance_log2"], exp)


def test_default_policy_is_median() -> None:
    res = _ex7_gene_run(None)
    assert approx(row(res, M)["abundance_log2"], 13.0)
    assert set(res.table["imputation_policy"]) == {"median"}


# ------------------------------------------------------------------------------ example 8


def test_example8_undetected_missing_gene_stays_na() -> None:
    rows = [acc(t, f"PEP{t}K", 2.0**v) for t, v in ((G, 10), (H, 12))]
    res = run(accepted_frame(rows), n_obs={t: 1 for t in TAGS}, detected=(G, H), policy="median")
    for t in (K, L, M):
        r = row(res, t)
        assert na(r["abundance_log2"])
        assert not bool(r["abundance_imputed"])
        assert "imputed_abundance" not in flags(r)
        assert r["imputation_policy"] == "median"


def test_undetected_gene_keeps_observed_value() -> None:
    # Q3: a single-peptide (undetected) gene keeps its observed value
    rows = [acc(G, "PEPGK", 1024.0), acc(H, "PEPHK", 4096.0)]
    res = run(accepted_frame(rows), n_obs={G: 1, H: 1}, detected=(G, K))
    assert approx(row(res, H)["abundance_log2"], 12.0)
    assert not bool(row(res, H)["abundance_imputed"])
    # fill statistic uses detected G only
    assert approx(row(res, K)["abundance_log2"], 10.0)
    assert res.fill_value == pytest.approx(10.0)


# ------------------------------------------------------------------------------ example 9


def test_example9_linear_zero_is_missing_then_imputed_if_detected() -> None:
    ab = c10({"P1": 0.0, "P2": 1024.0, "P3": 4096.0, "P4": 0.0})
    pm = pmap_frame([("P1", G, "protein_hit"), ("P2", H, "protein_hit"),
                     ("P3", K, "protein_hit"), ("P4", L, "protein_hit")])  # fmt: skip
    res = run(None, detected=(G, H, K), abundance=ab, pmap=pm, type="ibaq", c10_scale="linear")
    g = row(res, G)
    assert approx(g["abundance_log2"], 11.0), "median of 10 and 12"
    assert bool(g["abundance_imputed"])
    assert "imputed_abundance" in flags(g)
    l_ = row(res, L)
    assert na(l_["abundance_log2"]), "undetected linear zero stays NA"
    assert not bool(l_["abundance_imputed"])


def test_example9_protein_abundance_zero_is_na() -> None:
    out = _ab().protein_abundance(
        c10({"P1": 0.0, "P2": -5.0}),
        pmap_frame([("P1", G, "protein_hit"), ("P2", H, "protein_hit")]),
        "ibaq",
        "linear",
    )
    for t in (G, H):
        sel = out[out["locus_tag"] == t]
        assert sel.empty or na(sel.iloc[0]["value_log2"])


# ------------------------------------------------------------------------------ example 10


@pytest.mark.parametrize("with_rows", [False, True])
def test_example10_no_intensities_no_c10(with_rows: bool) -> None:
    rows = [acc(G, "PEPGK", None), acc(G, "PEPGGK", None, start=30)] if with_rows else []
    res = run(accepted_frame(rows), n_obs={G: 2}, detected=(G, H))
    assert res.source == "none"
    assert res.fill_value is None
    t = res.table
    assert len(t) == len(TAGS)
    assert (t["abundance_type"] == "none").all()
    for col in ("abundance_raw", "abundance_log2", "top3_log2", "abundance_scale",
                "imputation_policy"):  # fmt: skip
        assert t[col].isna().all(), col
    assert all(na(v) or v is False or v == False for v in t["abundance_imputed"])  # noqa: E712
    for tag in TAGS:
        assert "imputed_abundance" not in flags(row(res, tag))


# ------------------------------------------------------------------------------ source choice


def test_peptide_intensities_win_over_c10() -> None:
    ab = c10({"P1": 30.0})
    pm = pmap_frame([("P1", G, "protein_hit")])
    res = run(accepted_frame(ex1_rows()), n_obs={G: 4}, detected=(G,), abundance=ab, pmap=pm,
              type="maxlfq", c10_scale="log2")  # fmt: skip
    assert res.source == "peptide_intensity"
    r = row(res, G)
    assert approx(r["abundance_log2"], 8.965784)
    assert r["abundance_type"] == "ibaq"
    assert r["abundance_scale"] == "molar"


def test_ineligible_intensities_only_fall_back_to_c10() -> None:
    rows = [acc(G, "MULTIPEPK", 1e6, uniqueness="multi_identical")]
    ab = c10({"P1": 20.0})
    pm = pmap_frame([("P1", G, "protein_hit")])
    res = run(accepted_frame(rows), n_obs={G: 1}, detected=(G,), abundance=ab, pmap=pm,
              type="maxlfq", c10_scale="log2")  # fmt: skip
    assert res.source == "protein_table"
    assert approx(row(res, G)["abundance_log2"], 20.0)
    assert row(res, G)["abundance_scale"] == "relative"


# ------------------------------------------------------------------------------ edge cases


def test_all_detected_missing_no_fill_and_warning(caplog: pytest.LogCaptureFixture) -> None:
    rows = [acc(H, "PEPHK", 4096.0)]  # only the undetected H is quantified
    with warnings.catch_warnings(record=True) as rec, caplog.at_level(logging.DEBUG):
        warnings.simplefilter("always")
        res = run(accepted_frame(rows), n_obs={H: 1}, detected=(G, K))
    assert res.source == "peptide_intensity"
    assert res.fill_value is None
    assert not res.table["abundance_imputed"].fillna(False).astype(bool).any()
    assert na(row(res, G)["abundance_log2"])
    msgs = [str(w.message) for w in rec] + [r.getMessage() for r in caplog.records]
    assert any("no_observed_abundance" in m for m in msgs), msgs


def test_na_intensities_skipped_not_zero() -> None:
    rows = [acc(G, "AAAAAAAK", 100), acc(G, "CCCCCCCK", None, start=30),
            acc(G, "DDDDDDDK", 300, start=60)]  # fmt: skip
    res = run(accepted_frame(rows), n_obs={G: 4}, detected=(G,))
    r = row(res, G)
    assert approx(r["abundance_raw"], 100.0)
    assert approx(r["top3_log2"], math.log2(200.0))
    assert "top3_lt3" in flags(r)


def test_zero_intensity_skipped() -> None:
    rows = [acc(G, "AAAAAAAK", 100), acc(G, "CCCCCCCK", 0.0, start=30),
            acc(G, "DDDDDDDK", 300, start=60)]  # fmt: skip
    res = run(accepted_frame(rows), n_obs={G: 4}, detected=(G,))
    r = row(res, G)
    assert approx(r["top3_log2"], math.log2(200.0))
    assert "top3_lt3" in flags(r)
    out = _ab().peptide_abundance(accepted_frame(rows), pd.Series({G: 4}))
    assert int(out[out["locus_tag"] == G].iloc[0]["n_quant_peptides"]) == 2


def test_partial_gene_reported_but_not_in_fill_statistic() -> None:
    rows = [acc(t, f"PEP{t}K", 2.0**v) for t, v in ((G, 10), (H, 12), (L, 30))]
    res = run(accepted_frame(rows), n_obs={t: 1 for t in TAGS}, detected=(G, H, K),
              partial=(L,))  # fmt: skip
    assert approx(row(res, L)["abundance_log2"], 30.0), "partial gene value reported"
    assert approx(row(res, K)["abundance_log2"], 11.0), "fill = median(10, 12), not with 30"


def test_partial_gene_marked_detected_still_not_in_fill_statistic() -> None:
    ab = _ab()
    genes = genes_frame(TAGS, partial=(L,))
    rows = [acc(t, f"PEP{t}K", 2.0**v) for t, v in ((G, 10), (H, 12), (L, 30))]
    det = detection(TAGS, (G, H, K, L))  # a stub that wrongly marks the partial L detected
    res = ab.gene_abundance(genes, accepted_frame(rows), covariates_frame({t: 1 for t in TAGS},
                            TAGS), det, None, empty_pmap(), cfg())  # fmt: skip
    assert approx(row(res, K)["abundance_log2"], 11.0)


def test_genes_without_evidence_in_accepted_are_na() -> None:
    res = run(accepted_frame(ex1_rows()), n_obs={G: 4}, detected=(G,))
    for t in (H, K, L, M):
        r = row(res, t)
        assert na(r["abundance_raw"]) and na(r["abundance_log2"]) and na(r["top3_log2"])
        assert r["abundance_type"] == "ibaq"
        assert r["imputation_policy"] == "median"


def test_imputed_gene_top3_log2_not_filled() -> None:
    rows = [acc(G, "PEPGK", 1024.0)]
    res = run(accepted_frame(rows), n_obs={G: 1}, detected=(G, H))
    assert bool(row(res, H)["abundance_imputed"])
    assert na(row(res, H)["top3_log2"])


# ------------------------------------------------------------------------------ determinism


def _rich_rows() -> list[dict[str, Any]]:
    rng = np.random.default_rng(7)
    rows = []
    for t in (G, H, K):
        for j in range(int(rng.integers(1, 6))):
            rows.append(acc(t, f"P{t}{j}K", float(rng.lognormal(10, 2)), start=30 * j))
    rows.append(acc(L, "PLMULTIK", 5.0, uniqueness="multi_identical"))
    return rows


def test_deterministic_and_row_order_invariant() -> None:
    rows = _rich_rows()
    n_obs = {G: 3, H: 5, K: 2, L: 1, M: 0}
    a = run(accepted_frame(rows), n_obs=n_obs, detected=(G, H, L, M))
    b = run(accepted_frame(rows), n_obs=n_obs, detected=(G, H, L, M))
    c = run(accepted_frame(rows[::-1]), n_obs=n_obs, detected=(G, H, L, M))
    pd.testing.assert_frame_equal(a.table, b.table)
    pd.testing.assert_frame_equal(a.table, c.table)
    assert a.fill_value == b.fill_value == c.fill_value


def test_inputs_not_mutated() -> None:
    accepted = accepted_frame(_rich_rows())
    before = accepted.copy(deep=True)
    run(accepted, n_obs={G: 3, H: 5, K: 2}, detected=(G,))
    pd.testing.assert_frame_equal(accepted, before)


def test_matches_oracle_on_rich_example() -> None:
    rows = _rich_rows()
    n_obs = {G: 3, H: 5, K: 2, L: 1, M: 0}
    res = run(accepted_frame(rows), n_obs=n_obs, detected=(G, H, L, M))
    vals = ao.peptide_values(accepted_frame(rows))
    log2 = {t: (math.log2(ao.ibaq(vals[t], n_obs[t])) if t in vals else math.nan) for t in TAGS}
    det = {t: t in (G, H, L, M) for t in TAGS}
    filled, mask, fill = ao.impute(log2, det, "median")
    assert res.fill_value == pytest.approx(fill)
    for t in TAGS:
        r = row(res, t)
        assert (na(r["abundance_log2"]) and na(filled[t])) or approx(r["abundance_log2"], filled[t])
        assert bool(r["abundance_imputed"]) == mask[t]
        if t in vals:
            assert approx(r["top3_log2"], math.log2(ao.top3(vals[t])))
