"""Pure helpers of validation/proxies/te_proxies.py (issue M3-22; PR #115 review blockers 1-3).

Hand-built tables whose exact expected values distinguish the documented rules from near
variants (the reviewer's surviving mutants):

* ``pooled_decoy_fdr``: sum_m (D_m(e*) + 1) / A, with e* the largest accepted E per (decoy type,
  stratum), D counting decoys of the same source with best cost == m and E <= e*, A = the number
  of accepted placements; 0 when nothing is accepted (lead decision 12: unit decoy weights).
* ``inframe_inside`` / ``te1_result``: lead decision 10 (exact hits in frame inside an own gene
  count even when the locus is ranked ``rna_overlap`` or ``pseudogene``; reported as
  ``exact_inframe_masked``). ``peptide_hits.tsv`` coordinates are 1-based closed.
* ``te8_result``: ``nissle_exact`` / ``nissle_tolerant`` count accepted placements by best target
  cost, ``nissle_rejected`` placed but not accepted, ``nissle_unmapped`` without a placement row;
  ``passed`` needs both the K-12 share and the FDR condition.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from gippa.mapping.placement_fdr import PER_PEPTIDE_COLUMNS
from gippa.model.coords import Interval

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
te: Any = importlib.import_module("validation.proxies.te_proxies")


# ------------------------------------------------------------------ tables


def _per(rows: list[tuple[str, int, float, bool]], source: str = "database") -> pd.DataFrame:
    """Placement rows (peptide, m_target, e_target, accepted)."""
    dtype = {"database": "reverse", "denovo": "shuffle"}[source]
    return pd.DataFrame(
        [
            {
                "peptide": p,
                "source": source,
                "decoy_type": dtype,
                "length": 8,
                "m_target": m,
                "e_target": e,
                "placement_q": 0.0,
                "accepted": acc,
                "flags": "",
            }
            for p, m, e, acc in rows
        ],
        columns=list(PER_PEPTIDE_COLUMNS),
    ).astype({"m_target": int, "e_target": float, "accepted": bool})


def _costs(rows: list[tuple[str, str, int, float]]) -> pd.DataFrame:
    """Best decoy costs (peptide, source, m_decoy, e_decoy)."""
    return pd.DataFrame(rows, columns=["peptide", "source", "m_decoy", "e_decoy"]).astype(
        {"m_decoy": int, "e_decoy": float}
    )


# ------------------------------------------------------------------ pooled FDR (blocker 1)

# 13 accepted database placements: 12 in stratum 0 (e* = 1e-4), 1 in stratum 1 (e* = 1e-3);
# one rejected stratum-0 target with a large E that must not raise e*.
PER = _per(
    [("T01", 0, 1e-6, True), ("T02", 0, 1e-4, True)]
    + [(f"T{i:02d}", 0, 1e-7, True) for i in range(3, 13)]
    + [("T13", 1, 1e-3, True), ("R01", 0, 1e-2, False)]
)
COSTS = _costs(
    [
        ("D1", "database", 0, 1e-4),  # tied at e*_0: counted (<=)
        ("D2", "database", 0, 1e-5),  # below e*_0: counted
        ("D3", "database", 0, 1e-3),  # above e*_0: not counted (counted if e* used the reject)
        ("D4", "database", 1, 1e-3),  # stratum 1, tied at e*_1: counted
        ("D5", "database", 2, 1e-9),  # stratum 2 has no accepted target: never counted
        ("D6", "denovo", 0, 1e-9),  # other source, no accepted de novo target: not counted
    ]
)


def test_pooled_fdr_exact_value() -> None:
    # stratum 0: D1, D2 -> 2 + 1; stratum 1: D4 -> 1 + 1; A = 13 accepted (not 14 rows)
    assert te.pooled_decoy_fdr(COSTS, PER) == pytest.approx(5 / 13)


@pytest.mark.parametrize(
    ("variant", "value"),
    [
        ("no +1 per stratum", 3 / 13),
        ("e* = smallest accepted E (nanmin)", 3 / 13),
        ("e_decoy < e* (strict)", 3 / 13),
        ("m_decoy >= m", 7 / 13),
        ("denominator = all placement rows", 5 / 14),
        ("e* from all rows incl. rejected", 6 / 13),
    ],
)
def test_pooled_fdr_differs_from_near_variants(variant: str, value: float) -> None:
    assert te.pooled_decoy_fdr(COSTS, PER) != pytest.approx(value), variant


def test_pooled_fdr_nan_e_target_is_ignored_for_e_star() -> None:
    per = pd.concat([PER, _per([("T14", 0, float("nan"), True)])], ignore_index=True)
    assert te.pooled_decoy_fdr(COSTS, per) == pytest.approx(5 / 14)


def test_pooled_fdr_without_decoys_is_one_per_stratum_over_accepted() -> None:
    assert te.pooled_decoy_fdr(_costs([]), PER) == pytest.approx(2 / 13)


def test_pooled_fdr_is_capped_at_one() -> None:
    per = _per([("T1", 0, 1e-3, True)])
    costs = _costs([(f"D{i}", "database", 0, 1e-4) for i in range(5)])
    assert te.pooled_decoy_fdr(costs, per) == 1.0


def test_pooled_fdr_nothing_accepted_is_zero() -> None:
    assert te.pooled_decoy_fdr(COSTS, _per([("R1", 0, 1e-9, False)])) == 0.0
    assert te.pooled_decoy_fdr(COSTS, _per([])) == 0.0


def test_pooled_fdr_pools_decoy_types_separately() -> None:
    per = pd.concat(
        [_per([("T1", 0, 1e-4, True)]), _per([("N1", 0, 1e-6, True)], source="denovo")],
        ignore_index=True,
    )
    costs = _costs([("D1", "denovo", 0, 1e-5), ("D2", "database", 0, 1e-5)])
    # database stratum 0 (e* 1e-4): D2 -> 2; denovo stratum 0 (e* 1e-6): D1 above -> 1
    assert te.pooled_decoy_fdr(costs, per) == pytest.approx(min(1.0, 3 / 2))
    many = pd.concat([per, _per([(f"X{i}", 0, 1e-9, True) for i in range(8)])], ignore_index=True)
    assert te.pooled_decoy_fdr(costs, many) == pytest.approx(3 / 10)


# ------------------------------------------------------------------ inframe_inside (blocker 2)

L = 3000
PLUS = Interval("chr", 100, 400, "+")
MINUS_PARTIAL = Interval("chr", 500, 801, "-")  # 301 bp: frames from 5' and 3' ends differ
ORIGIN = Interval("chr", 2900, 3100, "+")  # unwrapped across the origin of a 3000-bp replicon


@pytest.mark.parametrize(
    ("start", "end", "strand", "gene", "circ", "expected"),
    [
        (103, 124, "+", PLUS, L, True),  # in frame
        (100, 400, "+", PLUS, L, True),  # the whole gene
        (104, 125, "+", PLUS, L, False),  # off frame
        (103, 124, "-", PLUS, L, False),  # antisense
        (97, 118, "+", PLUS, L, False),  # starts before the gene
        (382, 403, "+", PLUS, L, False),  # in frame but ends 3 nt past the gene end
        (379, 400, "+", PLUS, L, True),  # ends exactly at the gene end
        (777, 798, "-", MINUS_PARTIAL, L, True),  # minus strand: frame counted from the end
        (776, 797, "-", MINUS_PARTIAL, L, False),  # in frame from the start, not from the end
        (2990, 3011, "+", ORIGIN, L, True),  # crossing hit, unwrapped (end > L)
        (50, 71, "+", ORIGIN, L, True),  # crossing gene, hit given wrapped: needs +L
        (5951, 5972, "+", ORIGIN, L, True),  # hit given one turn ahead: needs -L
        (50, 71, "+", ORIGIN, None, False),  # linear contig: no shift
        (51, 72, "+", ORIGIN, L, False),  # wrapped but off frame
    ],
)
def test_inframe_inside(
    start: int, end: int, strand: str, gene: Interval, circ: int | None, expected: bool
) -> None:
    assert te.inframe_inside("chr", start, end, strand, gene, circ) is expected


def test_inframe_inside_other_contig() -> None:
    assert te.inframe_inside("plasmid", 103, 124, "+", PLUS, L) is False


# ------------------------------------------------------------------ te1_result (blocker 2)

GENES = {
    "G1": PLUS,
    "G2": MINUS_PARTIAL,
    "G3": ORIGIN,
    "G9": Interval("chr", 1500, 1800, "+"),  # a foreign gene
}


def _hit(pep: str, route: str, cls: str, gene: str, s0: int, e0: int, strand: str) -> dict:
    """One ``peptide_hits.tsv`` row; ``s0`` 0-based, written 1-based closed like the file."""
    return {"peptide": pep, "route": route, "locus_class": cls, "gene": gene, "contig": "chr",
            "start": str(s0 + 1), "end": str(e0), "strand": strand}  # fmt: skip


SOURCES = {p: {"G3" if p in ("pH", "pI") else "G2" if p == "pE" else "G1"} for p in [
    "pA", "pB", "pC", "pD", "pE", "pF", "pG", "pH", "pI", "pJ", "pK", "pL", "pM", "pN",
]}  # fmt: skip
EXCLUDED = {"pL": "G1"}
HITS = pd.DataFrame(
    [
        _hit("pA", "exact", "cds_inframe", "G1", 103, 124, "+"),  # good
        _hit("pB", "exact", "cds_inframe", "G9", 1503, 1524, "+"),  # foreign gene: miss
        _hit("pC", "tolerant", "cds_inframe", "G1", 103, 124, "+"),  # tolerant: miss
        _hit("pD", "exact", "rna_overlap", "", 106, 127, "+"),  # masked in frame
        _hit("pE", "exact", "pseudogene", "G2", 777, 798, "-"),  # masked, minus strand
        _hit("pF", "exact", "rna_overlap", "", 104, 125, "+"),  # masked off frame: miss
        _hit("pG", "exact", "rna_overlap", "", 103, 124, "-"),  # masked antisense: miss
        _hit("pH", "exact", "rna_overlap", "", 2990, 3011, "+"),  # masked, crossing (end > L)
        _hit("pI", "exact", "rna_overlap", "", 50, 71, "+"),  # masked, crossing, wrapped
        _hit("pJ", "exact", "rna_overlap", "", 109, 130, "+"),  # masked row first ...
        _hit("pJ", "exact", "cds_inframe", "G1", 109, 130, "+"),  # ... and a good row
        _hit("pK", "exact", "intergenic", "", 112, 133, "+"),  # not a masking class: miss
        _hit("pN", "exact", "rna_overlap", "", 382, 403, "+"),  # 3 nt past the gene end: miss
        _hit("pM", "tolerant", "rna_overlap", "", 115, 136, "+"),  # tolerant masked: miss
    ]
)
CIRC = {"chr": L}


@pytest.fixture(scope="module")
def te1() -> Any:
    return te.te1_result(SOURCES, EXCLUDED, HITS, GENES, CIRC)


def test_te1_result_counts(te1: Any) -> None:
    assert te1.counts == {
        "peptides": 14,
        "non_atg_nterm": 1,
        "exact_inframe": 6,  # pA pD pE pH pI pJ
        "misses": 7,  # pB pC pF pG pK pM pN
        "exact_inframe_masked": 4,  # pD pE pH pI (pJ has a good row)
        "sec_site": 0,
    }
    assert te1.shares["exact_inframe"] == pytest.approx(6 / 13)
    assert te1.passed is False


def test_te1_result_miss_and_exclusion_rows(te1: Any) -> None:
    kinds = {(r["kind"], r["peptide"]) for r in te1.exceptions}
    misses = {"pB", "pC", "pF", "pG", "pK", "pM", "pN"}
    assert {p for k, p in kinds if k == "miss"} == misses
    assert {p for k, p in kinds if k == "non_atg_nterm"} == {"pL"}
    for r in te1.exceptions:
        assert r["locus_tag"] == ("G1" if r["peptide"] != "pE" else "G2")


@pytest.mark.parametrize("pep", ["pD", "pE", "pH", "pI"])
def test_te1_each_masked_case_alone(pep: str) -> None:
    r = te.te1_result({pep: SOURCES[pep]}, {}, HITS[HITS["peptide"] == pep], GENES, CIRC)
    assert r.counts["exact_inframe"] == 1 and r.counts["exact_inframe_masked"] == 1
    assert r.passed is True


@pytest.mark.parametrize("pep", ["pB", "pC", "pF", "pG", "pK", "pM", "pN"])
def test_te1_each_miss_case_alone(pep: str) -> None:
    r = te.te1_result({pep: SOURCES[pep]}, {}, HITS[HITS["peptide"] == pep], GENES, CIRC)
    assert r.counts["exact_inframe"] == 0 and r.counts["misses"] == 1
    assert r.passed is False


def test_te1_good_and_masked_rows_count_once_as_good_in_either_order() -> None:
    rows = HITS[HITS["peptide"] == "pJ"]
    for order in (rows, rows.iloc[::-1]):
        r = te.te1_result({"pJ": {"G1"}}, {}, order, GENES, CIRC)
        assert r.counts["exact_inframe"] == 1
        assert r.counts["exact_inframe_masked"] == 0


def test_te1_sorf_at_own_gene_is_good() -> None:
    hits = pd.DataFrame([_hit("pS", "exact", "sorf", "G1", 103, 124, "+")])
    r = te.te1_result({"pS": {"G1"}}, {}, hits, GENES, CIRC)
    assert r.counts["exact_inframe"] == 1 and r.counts["exact_inframe_masked"] == 0


def test_te1_masked_on_linear_contig_is_not_shifted() -> None:
    hits = HITS[HITS["peptide"] == "pI"]
    r = te.te1_result({"pI": {"G3"}}, {}, hits, GENES, {"chr": None})
    assert r.counts["misses"] == 1


def test_te1_unplaced_peptide_is_listed_as_unplaced() -> None:
    r = te.te1_result({"pZ": {"G1"}}, {}, HITS.iloc[0:0], GENES, CIRC)
    assert r.counts["misses"] == 1
    assert r.exceptions == [
        {"kind": "miss", "locus_tag": "G1", "peptide": "pZ", "hits": "unplaced"}
    ]


def test_te1_exceptions_capped_at_50() -> None:
    sources = {f"p{i:03d}": {"G1"} for i in range(80)}
    r = te.te1_result(sources, {}, HITS.iloc[0:0], GENES, CIRC)
    assert r.counts["misses"] == 80 and len(r.exceptions) == 50


# ------------------------------------------------------------------ te8_result (blocker 3)


def test_te8_result_counts_and_shares() -> None:
    sources = {f"p{i}": {f"P{i}"} for i in range(1, 7)}
    per = _per(
        [
            ("p1", 0, 1e-6, True),  # accepted exact
            ("p2", 1, 1e-5, True),  # accepted tolerant (stratum 1)
            ("p3", 0, 1e-1, False),  # placed, rejected
            ("p4", 2, 1e-5, True),  # accepted tolerant (stratum 2)
        ]
    )  # p5, p6: no placement row (unmapped)
    r = te.te8_result(sources, {"p6"}, {"p1", "p2", "p3", "p6"}, per, _costs([]))
    assert r.counts == {
        "peptides": 6,
        "non_atg_nterm": 1,
        "k12_exact": 3,  # p1 p2 p3; p6 is excluded although it hits
        "nissle_exact": 1,  # accepted only (p3 has m = 0 but is rejected)
        "nissle_tolerant": 2,
        "nissle_rejected": 1,
        "nissle_unmapped": 2,  # no row at all (not 6 - 3 accepted)
    }
    assert r.shares["k12_exact"] == pytest.approx(3 / 5)
    assert r.shares["nissle_exact"] == pytest.approx(1 / 6)
    assert r.shares["nissle_tolerant"] == pytest.approx(2 / 6)
    assert r.shares["nissle_unmapped"] == pytest.approx(2 / 6)
    assert r.shares["nissle_placement_fdr"] == pytest.approx(1.0)  # 3 strata, 3 accepted
    assert {e["peptide"] for e in r.exceptions} == {"p4", "p5"}
    assert all(e["kind"] == "k12_miss" for e in r.exceptions)
    assert r.passed is False


def _pass_case(n_decoys: int, k12_hits: int = 300) -> Any:
    sources = {f"q{i:03d}": {f"Q{i}"} for i in range(300)}
    per = _per([(f"q{i:03d}", 0, 1e-6, True) for i in range(200)])
    costs = _costs([(f"D{i}", "database", 0, 1e-7) for i in range(n_decoys)])
    hits = set(sorted(sources)[:k12_hits])
    return te.te8_result(sources, set(), hits, per, costs)


def test_te8_passes_when_both_conditions_hold() -> None:
    r = _pass_case(0)  # FDR = 1/200 = 0.005, K-12 share 1.0
    assert r.shares["nissle_placement_fdr"] == pytest.approx(0.005)
    assert r.passed is True


def test_te8_fails_on_fdr_alone() -> None:
    r = _pass_case(2)  # FDR = 3/200 = 0.015 > 0.01, K-12 share still 1.0
    assert r.shares["k12_exact"] == 1.0
    assert r.shares["nissle_placement_fdr"] == pytest.approx(0.015)
    assert r.passed is False


def test_te8_fdr_boundary_is_inclusive() -> None:
    r = _pass_case(1)  # FDR = 2/200 = 0.01 exactly
    assert r.shares["nissle_placement_fdr"] == pytest.approx(0.01)
    assert r.passed is True


def test_te8_fails_on_k12_share_alone() -> None:
    r = _pass_case(0, k12_hits=298)  # 298/300 = 0.99333 < 0.995
    assert r.shares["nissle_placement_fdr"] <= 0.01
    assert r.passed is False


def test_te8_tolerant_stratum_enters_the_fdr() -> None:
    sources = {f"q{i:03d}": {f"Q{i}"} for i in range(300)}
    per = _per(
        [(f"q{i:03d}", 0, 1e-6, True) for i in range(250)]
        + [(f"q{i:03d}", 1, 1e-4, True) for i in range(250, 300)]
    )
    costs = _costs([("D1", "database", 1, 1e-4), ("D2", "database", 1, 1e-3)])
    r = te.te8_result(sources, set(), set(sources), per, costs)
    # stratum 0: 0 + 1; stratum 1 (tolerant, e* 1e-4): D1 + 1 -> (1 + 2) / 300
    assert r.counts["nissle_tolerant"] == 50 and r.counts["nissle_exact"] == 250
    assert r.shares["nissle_placement_fdr"] == pytest.approx(3 / 300)
    assert r.passed is True


def test_te8_empty_inputs_do_not_pass() -> None:
    r = te.te8_result({}, set(), set(), _per([]), _costs([]))
    assert r.counts["peptides"] == 0
    assert r.shares["k12_exact"] == 0.0
    assert r.passed is False


def test_te8_all_excluded_does_not_pass() -> None:
    r = te.te8_result({"p1": {"P1"}}, {"p1"}, {"p1"}, _per([]), _costs([]))
    assert r.counts["non_atg_nterm"] == 1 and r.counts["k12_exact"] == 0
    assert r.passed is False


def test_proxy_result_to_json_round_trip() -> None:
    import json

    r = _pass_case(0)
    d = json.loads(json.dumps(r.to_json()))
    assert d["proxy"] == "TE8" and d["passed"] is True
    assert d["counts"] == r.counts
    assert np.isclose(d["shares"]["nissle_placement_fdr"], 0.005)


# ------------------------------------------------------------------ gene helpers (decisions 5, 11)


def _gene(tag: str, kind: str, aa: str, nt: str | None) -> Any:
    from gippa.model.genes import Gene

    return Gene(tag, kind, Interval("chr", 0, 3 * len(aa) + 3, "+"), aa, nt, None)  # type: ignore[arg-type]


def test_strict_non_atg_needs_both_conditions() -> None:
    assert te.is_strict_non_atg(_gene("a", "cds", "MAAK", "GTGGCGGCGAAATAA")) is True
    assert te.is_strict_non_atg(_gene("b", "cds", "VAAK", "GTGGCGGCGAAATAA")) is False  # Edge
    assert te.is_strict_non_atg(_gene("c", "cds", "MAAK", "ATGGCGGCGAAATAA")) is False
    assert te.is_strict_non_atg(_gene("d", "cds", "MAAK", "atgGCGGCGAAATAA")) is False
    assert te.is_strict_non_atg(_gene("e", "cds", "MAAK", None)) is False


def test_coding_genes_keep_cds_sorf_and_pseudogene() -> None:
    genes = [
        _gene("c", "cds", "MAAK", "ATG"),
        _gene("s", "sorf", "MAAK", "ATG"),
        _gene("p", "pseudogene", "MAAK", "ATG"),
        _gene("n", "cds", "", "ATG"),  # no protein
    ]
    assert [g.locus_tag for g in te.coding_genes(genes)] == ["c", "s", "p"]


def test_nterm_products_start_at_zero_only() -> None:
    # position 0 is a cut site: AAAAAAAR (start 1) is a product but not N-terminal
    assert te.nterm_products("KAAAAAAAR") == ["KAAAAAAAR"]


def test_non_atg_nterm_maps_products_to_smallest_locus() -> None:
    aa = "MAAAAAAAKWWWWWWWR"
    genes = [_gene("g2", "cds", aa, "TTG"), _gene("g1", "cds", aa, "GTG"),
             _gene("g3", "cds", aa, "ATG")]  # fmt: skip
    assert te.non_atg_nterm(genes) == {"MAAAAAAAK": "g1", "MAAAAAAAKWWWWWWWR": "g1"}


def test_te1_excluded_peptide_with_a_good_hit_is_not_counted_exact() -> None:
    hits = pd.DataFrame([_hit("pX", "exact", "cds_inframe", "G1", 103, 124, "+")])
    r = te.te1_result({"pX": {"G1"}, "pY": {"G1"}}, {"pX": "G1"}, hits, GENES, CIRC)
    assert r.counts["non_atg_nterm"] == 1
    assert r.counts["exact_inframe"] == 0 and r.counts["misses"] == 1
    assert r.counts["peptides"] == 2


def test_te1_all_excluded_does_not_pass() -> None:
    r = te.te1_result({"pX": {"G1"}}, {"pX": "G1"}, HITS.iloc[0:0], GENES, CIRC)
    assert r.shares["exact_inframe"] == 0.0 and r.passed is False
