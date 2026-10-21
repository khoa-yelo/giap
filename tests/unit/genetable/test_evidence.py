"""M4-11 accepted evidence and per-gene evidence counts (docs/specs/M4-11-accepted-evidence.md).

Spec: GIAP_SPEC.md §7.3.1 (five conditions, distinct peptides), §7.2.6 (uniqueness and locus
classes), §7.2.7 (genome relation), §9.3 (`genes.tsv` evidence columns), §9.4 (gene flags), §6.5
(origin-crossing intervals). Lead decisions: pseudogene placements assign to the pseudogene
(M4-index decision 1); `cterm_extension` is an assigning class (ADR-0013, M4-06 review decision 2).

Every scenario is also checked against the independent oracle in tests/helpers/evidence_oracle.py.
Resolved ambiguities (pinned here):
- `gene_evidence` returns one row per C2 gene in C2 order (a `locus_tag` column is used if present).
- `evidence_flags` is compared as a set (order not specified), without duplicates.
- NA `placement_q` without `same_db` is not an accepted placement (`NA <= alpha` is false).
- A placement stored after the origin (`start < end <= length`) inside an origin-crossing gene adds
  coverage (circular intersection); only accepted placements add coverage (Q2 default).
- `is_contaminant` is not a §7.3.1 condition: a placed, assigned contaminant peptide is accepted.
- C9 decoy rows never appear in C11; they produce nothing and raise nothing.
"""

from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import evidence_oracle as eo
from evidence_oracle import hit, pep

G, G2, H = "GSIM_00005", "GSIM_00010", "GSIM_00015"
GENES = [eo.gene(G, 1000, 1300), eo.gene(G2, 5000, 5300), eo.gene(H, 9000, 9300, strand="-")]
REPLICONS = (("chr1", 100_000, "circular"),)


def _ev() -> Any:
    from giap.genetable import evidence

    return evidence


def _cfg(q: float = 0.01, alpha: float = 0.01, relation: str = "related") -> Any:
    from giap.config import GiapConfig

    return GiapConfig(
        genome_relation=relation,
        input={"q_value_max": q},
        stats={"alpha": {"placement": alpha}},
    )


@dataclasses.dataclass
class Run:
    accepted: pd.DataFrame
    evidence: pd.DataFrame
    genes: pd.DataFrame

    def row(self, tag: str) -> dict[str, Any]:
        ev = self.evidence
        if "locus_tag" in ev.columns:
            sel = ev[ev["locus_tag"] == tag]
            assert len(sel) == 1, f"{tag}: {len(sel)} rows"
            return sel.iloc[0].to_dict()
        pos = self.genes["locus_tag"].tolist().index(tag)
        return ev.iloc[pos].to_dict()


def run(
    hit_rows: list[dict[str, Any]],
    pep_rows: list[dict[str, Any]] | None,
    *,
    genes: list[dict[str, Any]] | None = None,
    replicons: tuple[tuple[str, int, str], ...] = REPLICONS,
    q: float = 0.01,
    alpha: float = 0.01,
    relation: str = "related",
    min_peptides: int = 2,
    check_oracle: bool = True,
) -> Run:
    ev = _ev()
    gdf = eo.genes_frame(GENES if genes is None else genes)
    hits = eo.hits_frame(hit_rows)
    peps = eo.peptides_frame(pep_rows) if pep_rows is not None else None
    reps = eo.replicons_frame(replicons)
    acc = ev.accepted_peptides(hits, peps, _cfg(q, alpha, relation))
    out = ev.gene_evidence(gdf, hits, acc, reps, min_peptides=min_peptides)
    r = Run(acc, out, gdf)
    if check_oracle:
        rules = eo.Rules(q, alpha, relation)
        assert_matches_oracle(r, eo.gene_evidence(gdf, hits, peps, reps, rules, min_peptides))
        pairs = set(zip(acc["gene"].astype(str), acc["peptide"].astype(str), strict=True))
        assert pairs == eo.accepted_pairs(hits, peps, rules)
        assert len(acc) == len(pairs), "one row per accepted (gene, peptide) pair"
    return r


def _is_na(x: Any) -> bool:
    return x is None or (not isinstance(x, (list, tuple, np.ndarray)) and bool(pd.isna(x)))


def flags_of(row: dict[str, Any]) -> set[str]:
    cell = row["evidence_flags"]
    if _is_na(cell):
        return set()
    items = [str(f) for f in cell]
    assert len(items) == len(set(items)), f"duplicate flags {items}"
    return set(items)


def assert_matches_oracle(r: Run, expected: dict[str, dict[str, Any]]) -> None:
    assert len(r.evidence) == len(r.genes), "one row per C2 gene"
    if "locus_tag" in r.evidence.columns:
        assert r.evidence["locus_tag"].tolist() == r.genes["locus_tag"].tolist(), "C2 order"
    for tag, exp in expected.items():
        got = r.row(tag)
        for col in (
            "n_peptides_unique",
            "n_peptides_multi_identical",
            "n_peptides_shared",
            "n_peptides_tolerant",
            "n_accepted",
        ):
            assert int(got[col]) == exp[col], (tag, col, got[col], exp[col])
        assert bool(got["detected_2pep"]) == exp["detected_2pep"], (tag, "detected_2pep")
        assert float(got["coverage_frac"]) == pytest.approx(exp["coverage_frac"], abs=1e-12), tag
        if math.isnan(exp["best_peptide_q"]):
            assert _is_na(got["best_peptide_q"]), (tag, got["best_peptide_q"])
        else:
            assert float(got["best_peptide_q"]) == pytest.approx(exp["best_peptide_q"]), tag
        assert flags_of(got) == exp["evidence_flags"], (tag, got["evidence_flags"])


# ------------------------------------------------------------------------------- examples 1-10

P1, P2 = "PEPTIDEK", "ANOTHERPEPK"
EX1_HIT = hit(P1, G, 1030, locus_class="cds_inframe", uniqueness="unique", route="exact")


@dataclasses.dataclass(frozen=True)
class Example:
    hits: list[dict[str, Any]]
    peps: list[dict[str, Any]]
    expect: dict[str, dict[str, Any]]
    relation: str = "related"


EXAMPLES: dict[str, Example] = {
    "ex1_one_peptide": Example(
        [EX1_HIT],
        [pep(P1, 0.001)],
        {
            G: {
                "n_peptides_unique": 1,
                "n_accepted": 1,
                "detected_2pep": False,
                "best_peptide_q": 0.001,
                "flags": set(),
            }
        },
    ),
    "ex2_two_peptides": Example(
        [EX1_HIT, hit(P2, G, 1100)],
        [pep(P1, 0.001), pep(P2, 0.004)],
        {
            G: {
                "n_peptides_unique": 2,
                "n_accepted": 2,
                "detected_2pep": True,
                "best_peptide_q": 0.001,
            }
        },
    ),
    "ex3_c9_q_too_high": Example(
        [EX1_HIT],
        [pep(P1, 0.02)],
        {
            G: {
                "n_peptides_unique": 0,
                "n_accepted": 0,
                "best_peptide_q": None,
                "detected_2pep": False,
                "coverage_frac": 0.0,
            }
        },
    ),
    "ex4_same_db": Example(
        [hit(P1, G, 1030, placement_q=0.03, flags=["same_db"])],
        [pep(P1, 0.001)],
        {G: {"n_peptides_unique": 1, "n_accepted": 1, "best_peptide_q": 0.001}},
    ),
    "ex5_charges_runs_once": Example(
        [EX1_HIT],
        [pep(P1, 0.001, n_rows=4, n_runs=2, modified=[P1, "PEPTIDEK[+16]"])],
        {G: {"n_peptides_unique": 1, "n_accepted": 1, "detected_2pep": False}},
    ),
    "ex6_multi_identical": Example(
        [
            hit(P1, G, 1030, uniqueness="multi_identical"),
            hit(P1, G2, 5030, uniqueness="multi_identical"),
        ],
        [pep(P1, 0.001)],
        {
            t: {
                "n_peptides_multi_identical": 1,
                "n_peptides_unique": 0,
                "n_accepted": 1,
                "flags": {"multi_identical"},
            }
            for t in (G, G2)
        },
    ),
    "ex7_shared": Example(
        [hit(P1, G, 1030, uniqueness="shared"), hit(P1, H, 9030, strand="-", uniqueness="shared")],
        [pep(P1, 0.001)],
        {
            G: {
                "n_peptides_shared": 1,
                "n_accepted": 0,
                "n_peptides_unique": 0,
                "best_peptide_q": None,
                "flags": {"shared_evidence"},
            }
        },
    ),
    "ex8_tolerant_related": Example(
        [hit(P1, G, 1030, route="tolerant", n_diff=1, placement_q=0.004)],
        [pep(P1, 0.001)],
        {
            G: {
                "n_accepted": 1,
                "n_peptides_tolerant": 1,
                "n_peptides_unique": 1,
                "flags": {"tolerant_evidence", "strain_difference"},
            }
        },
    ),
    "ex8_tolerant_same": Example(
        [hit(P1, G, 1030, route="tolerant", n_diff=1, placement_q=0.004)],
        [pep(P1, 0.001)],
        {G: {"n_accepted": 0, "n_peptides_tolerant": 0, "flags": set()}},
        relation="same",
    ),
    "ex9_nterm_extension": Example(
        [hit(P1, G, 982, locus_class="nterm_extension")],
        [pep(P1, 0.001)],
        {G: {"n_accepted": 1, "flags": {"nterm_extension_evidence"}}},
    ),
    "ex10_coverage": Example(
        [hit("AAAAAAAAAK", G, 1000, 1030), hit("CCCCCCCCCCCCCCK", G, 1015, 1060)],
        [pep("AAAAAAAAAK", 0.001), pep("CCCCCCCCCCCCCCK", 0.002)],
        {G: {"coverage_frac": 60 / 300, "n_accepted": 2, "detected_2pep": True}},
    ),
}


@pytest.mark.parametrize("name", list(EXAMPLES))
def test_examples(name: str) -> None:
    ex = EXAMPLES[name]
    r = run(ex.hits, ex.peps, relation=ex.relation)
    for tag, exp in ex.expect.items():
        row = r.row(tag)
        for col, val in exp.items():
            if col == "flags":
                assert flags_of(row) == val, (name, tag, row["evidence_flags"])
            elif val is None:
                assert _is_na(row[col]), (name, tag, col, row[col])
            elif isinstance(val, bool):
                assert bool(row[col]) is val, (name, tag, col)
            elif isinstance(val, float):
                assert float(row[col]) == pytest.approx(val), (name, tag, col)
            else:
                assert int(row[col]) == val, (name, tag, col, row[col])


def test_example3_c9_wins_over_c11_best_q_input() -> None:
    r = run([hit(P1, G, 1030, best_q=0.001)], [pep(P1, 0.02)])
    assert r.accepted.empty


def test_example5_duplicate_c11_rows_counted_once() -> None:
    r = run([EX1_HIT, dict(EX1_HIT)], [pep(P1, 0.001, n_rows=4, n_runs=2)])
    assert len(r.accepted) == 1
    assert int(r.row(G)["n_accepted"]) == 1


def test_example6_multi_identical_not_abundance_eligible() -> None:
    r = run(EXAMPLES["ex6_multi_identical"].hits, EXAMPLES["ex6_multi_identical"].peps)
    assert sorted(r.accepted["gene"]) == [G, G2]
    assert not r.accepted["abundance_eligible"].astype(bool).any()


def test_example7_shared_not_accepted_anywhere() -> None:
    r = run(EXAMPLES["ex7_shared"].hits, EXAMPLES["ex7_shared"].peps)
    assert r.accepted.empty
    assert int(r.row(H)["n_peptides_shared"]) == 1


def test_example8_tolerant_related_is_abundance_eligible() -> None:
    r = run(EXAMPLES["ex8_tolerant_related"].hits, EXAMPLES["ex8_tolerant_related"].peps)
    assert r.accepted["abundance_eligible"].astype(bool).tolist() == [True]
    assert r.accepted["route"].tolist() == ["tolerant"]


def test_example9_nterm_extension_not_abundance_eligible() -> None:
    r = run(EXAMPLES["ex9_nterm_extension"].hits, EXAMPLES["ex9_nterm_extension"].peps)
    assert r.accepted["abundance_eligible"].astype(bool).tolist() == [False]
    assert r.accepted["locus_class"].tolist() == ["nterm_extension"]


# ------------------------------------------------------------------------------- accepted rows


def test_accepted_columns_and_values_from_c9() -> None:
    r = run(
        [hit(P1, G, 1030, best_q=0.002)],
        [pep(P1, 0.002, pep_value=0.03, score=77.0, intensity=5.0e5)],
    )
    acc = r.accepted
    assert set(eo.ACCEPTED_COLUMNS) <= set(acc.columns)
    row = acc.iloc[0]
    assert (row["gene"], row["peptide"]) == (G, P1)
    assert float(row["q_value"]) == pytest.approx(0.002)
    assert float(row["pep"]) == pytest.approx(0.03)
    assert float(row["score"]) == pytest.approx(77.0)
    assert float(row["intensity"]) == pytest.approx(5.0e5)
    assert (row["uniqueness"], row["route"], row["locus_class"]) == (
        "unique",
        "exact",
        "cds_inframe",
    )
    assert (row["contig"], int(row["start"]), int(row["end"])) == ("chr1", 1030, 1054)
    assert not bool(row["crosses_origin"])
    assert bool(row["abundance_eligible"]) is True


def test_na_q_passes_engine_threshold() -> None:
    r = run([EX1_HIT], [pep(P1, None)])
    assert int(r.row(G)["n_accepted"]) == 1
    assert _is_na(r.row(G)["best_peptide_q"])


@pytest.mark.parametrize(("q", "ok"), [(0.01, True), (0.0100001, False), (0.0, True)])
def test_engine_threshold_boundary(q: float, ok: bool) -> None:
    r = run([EX1_HIT], [pep(P1, q)])
    assert int(r.row(G)["n_accepted"]) == int(ok)


@pytest.mark.parametrize(
    ("pq", "flags", "ok"),
    [
        (0.01, [], True),
        (0.0100001, [], False),
        (0.5, ["same_db"], True),
        (None, [], False),
        (None, ["same_db"], True),
        (0.03, ["strain_difference"], False),
    ],
)
def test_placement_acceptance(pq: float | None, flags: list[str], ok: bool) -> None:
    r = run([hit(P1, G, 1030, placement_q=pq, flags=flags)], [pep(P1, 0.001)])
    assert int(r.row(G)["n_accepted"]) == int(ok)


def test_config_thresholds_are_read() -> None:
    hits = [hit(P1, G, 1030, placement_q=0.03)]
    assert int(run(hits, [pep(P1, 0.02)]).row(G)["n_accepted"]) == 0
    assert int(run(hits, [pep(P1, 0.02)], q=0.05, alpha=0.05).row(G)["n_accepted"]) == 1
    assert int(run(hits, [pep(P1, 0.02)], q=0.05, alpha=0.01).row(G)["n_accepted"]) == 0
    assert int(run(hits, [pep(P1, 0.02)], q=0.01, alpha=0.05).row(G)["n_accepted"]) == 0


def test_same_relation_keeps_exact_hits() -> None:
    r = run(
        [EX1_HIT, hit(P2, G, 1100, route="tolerant", n_diff=1)],
        [pep(P1, 0.001), pep(P2, 0.001)],
        relation="same",
    )
    row = r.row(G)
    assert int(row["n_accepted"]) == 1 and int(row["n_peptides_tolerant"]) == 0
    assert r.accepted["peptide"].tolist() == [P1]


@pytest.mark.parametrize(
    "cls", ["cds_inframe", "sorf", "nterm_extension", "cterm_extension", "pseudogene"]
)
def test_assigning_classes_accepted(cls: str) -> None:
    genes = [
        eo.gene(
            G, 1000, 1300, feature_type={"sorf": "sorf", "pseudogene": "pseudogene"}.get(cls, "cds")
        )
    ]
    r = run([hit(P1, G, 1030, locus_class=cls)], [pep(P1, 0.001)], genes=genes)
    assert int(r.row(G)["n_accepted"]) == 1
    assert {"cds_inframe", "sorf", "nterm_extension", "cterm_extension", "pseudogene"} <= set(
        _ev().ASSIGNING_CLASSES
    )


@pytest.mark.parametrize("cls", ["cds_outframe", "antisense", "intergenic", "rna_overlap"])
def test_non_assigning_class_is_second_guard(cls: str) -> None:
    # M3 sets gene NA on non-assigning rows; a gene set anyway must not count (class check).
    r = run([hit(P1, G, 1031, locus_class=cls)], [pep(P1, 0.001)])
    assert int(r.row(G)["n_accepted"]) == 0
    assert r.accepted.empty


def test_outframe_row_on_second_gene_only_assigning_counts() -> None:
    r = run([hit(P1, G, 1030), hit(P1, None, 5031, locus_class="cds_outframe")], [pep(P1, 0.001)])
    assert int(r.row(G)["n_accepted"]) == 1
    assert int(r.row(G2)["n_accepted"]) == 0
    assert r.accepted["gene"].tolist() == [G]


def test_gene_na_rows_never_accepted() -> None:
    r = run([hit(P1, None, 20_000, locus_class="intergenic")], [pep(P1, 0.001)])
    assert r.accepted.empty


def test_decoy_rows_in_c9_produce_nothing() -> None:
    r = run(
        [EX1_HIT], [pep(P1, 0.001), pep("KEDITPEP", 0.0, is_decoy=True, protein_ids=["rev_" + G])]
    )
    assert r.accepted["peptide"].tolist() == [P1]
    assert int(r.row(G)["n_accepted"]) == 1


def test_contaminant_flag_is_not_a_condition() -> None:
    r = run([EX1_HIT], [pep(P1, 0.001, is_contaminant=True, protein_ids=["CONT_P02769"])])
    assert int(r.row(G)["n_accepted"]) == 1


def test_missing_c9_row_raises_naming_peptide() -> None:
    with pytest.raises(ValueError, match="MISSINGPEPK"):
        run([EX1_HIT, hit("MISSINGPEPK", G, 1100)], [pep(P1, 0.001)], check_oracle=False)


# ------------------------------------------------------------------------------- per-gene counts


def test_multi_and_unique_on_one_gene_detect() -> None:
    r = run(
        [
            EX1_HIT,
            hit(P2, G, 1100, uniqueness="multi_identical"),
            hit(P2, G2, 5100, uniqueness="multi_identical"),
        ],
        [pep(P1, 0.001), pep(P2, 0.003)],
    )
    row = r.row(G)
    assert (int(row["n_peptides_unique"]), int(row["n_peptides_multi_identical"])) == (1, 1)
    assert int(row["n_accepted"]) == 2 and bool(row["detected_2pep"])
    assert flags_of(row) == {"multi_identical"}


def test_shared_counts_distinct_peptides_and_ignores_q() -> None:
    hits = [
        hit(P1, G, 1030, uniqueness="shared"),
        hit(P1, G, 1130, uniqueness="shared"),
        hit(P2, G, 1200, uniqueness="shared"),
        hit(P1, H, 9030, uniqueness="shared"),
    ]
    r = run(hits, [pep(P1, 0.001), pep(P2, 0.5)])
    assert int(r.row(G)["n_peptides_shared"]) == 2
    assert int(r.row(H)["n_peptides_shared"]) == 1
    assert int(r.row(G2)["n_peptides_shared"]) == 0


def test_partial_gene_counts_but_not_detected() -> None:
    genes = [eo.gene(G, 1000, 1300, partial=True), *GENES[1:]]
    r = run([EX1_HIT, hit(P2, G, 1100)], [pep(P1, 0.001), pep(P2, 0.002)], genes=genes)
    row = r.row(G)
    assert int(row["n_accepted"]) == 2 and int(row["n_peptides_unique"]) == 2
    assert bool(row["detected_2pep"]) is False


@pytest.mark.parametrize(
    ("k", "n", "detected"), [(1, 1, True), (2, 1, False), (3, 2, False), (3, 3, True)]
)
def test_min_peptides(k: int, n: int, detected: bool) -> None:
    seqs = ["AAAAAAAK", "CCCCCCCK", "DDDDDDDK"][:n]
    r = run(
        [hit(s, G, 1030 + 30 * i) for i, s in enumerate(seqs)],
        [pep(s, 0.001) for s in seqs],
        min_peptides=k,
    )
    assert bool(r.row(G)["detected_2pep"]) is detected


def test_best_q_over_accepted_only() -> None:
    r = run([EX1_HIT, hit(P2, G, 1100, placement_q=0.5)], [pep(P1, 0.001), pep(P2, 0.0001)])
    assert float(r.row(G)["best_peptide_q"]) == pytest.approx(0.001)


def test_coverage_uses_accepted_only() -> None:
    r = run(
        [hit("AAAAAAAAAK", G, 1000, 1030), hit("CCCCCCCCCCK", G, 1100, 1133)],
        [pep("AAAAAAAAAK", 0.001), pep("CCCCCCCCCCK", 0.5)],
    )
    assert float(r.row(G)["coverage_frac"]) == pytest.approx(30 / 300)


def test_coverage_clipped_to_gene_for_nterm_extension() -> None:
    # placement [970, 1030) overlaps the gene [1000, 1300) on 30 bp.
    r = run(
        [hit("MAAAAAAAAAAAAAAAAAAK", G, 970, 1030, locus_class="nterm_extension")],
        [pep("MAAAAAAAAAAAAAAAAAAK", 0.001)],
    )
    assert float(r.row(G)["coverage_frac"]) == pytest.approx(30 / 300)


def test_coverage_minus_strand_gene() -> None:
    r = run([hit("AAAAAAAAAK", H, 9270, 9300, strand="-")], [pep("AAAAAAAAAK", 0.001)])
    assert float(r.row(H)["coverage_frac"]) == pytest.approx(30 / 300)


# ------------------------------------------------------------------------------- origin crossing

L = 10_000
OGENE = eo.gene("GSIM_00100", L - 150, L + 150, crosses_origin=True)


def test_origin_crossing_gene_with_spanning_peptide() -> None:
    p = "SPANNINGORIGINPEPTIDK"[:20]
    r = run(
        [hit(p, "GSIM_00100", L - 30, L + 30, crosses_origin=True)],
        [pep(p, 0.001)],
        genes=[OGENE],
        replicons=(("chr1", L, "circular"),),
    )
    row = r.row("GSIM_00100")
    assert float(row["coverage_frac"]) == pytest.approx(60 / 300)
    assert int(row["n_accepted"]) == 1
    acc = r.accepted.iloc[0]
    assert bool(acc["crosses_origin"]) and int(acc["end"]) == L + 30


def test_origin_crossing_gene_with_wrapped_and_spanning_peptides() -> None:
    a, b = "AAAAAAAAAAAAAAAAAAAK", "CCCCCCCCCK"
    hits = [
        hit(a, "GSIM_00100", L - 30, L + 30, crosses_origin=True),
        hit(b, "GSIM_00100", 60, 90),
        hit("DDDDDDDDDK", "GSIM_00100", L - 150, L - 120),
    ]
    peps = [pep(a, 0.001), pep(b, 0.002), pep("DDDDDDDDDK", 0.003)]
    r = run(hits, peps, genes=[OGENE], replicons=(("chr1", L, "circular"),))
    assert float(r.row("GSIM_00100")["coverage_frac"]) == pytest.approx((60 + 30 + 30) / 300)


def test_origin_crossing_overlapping_parts_not_double_counted() -> None:
    a, b = "AAAAAAAAAAAAAAAAAAAK", "CCCCCCCCCK"
    hits = [hit(a, "GSIM_00100", L - 30, L + 30, crosses_origin=True), hit(b, "GSIM_00100", 0, 30)]
    r = run(
        hits, [pep(a, 0.001), pep(b, 0.002)], genes=[OGENE], replicons=(("chr1", L, "circular"),)
    )
    assert float(r.row("GSIM_00100")["coverage_frac"]) == pytest.approx(60 / 300)


def test_genes_on_two_contigs_do_not_share_coverage() -> None:
    genes = [eo.gene(G, 1000, 1300), eo.gene("GSIM_00200", 1000, 1300, contig="p1")]
    r = run(
        [hit("AAAAAAAAAK", G, 1000, 1030)],
        [pep("AAAAAAAAAK", 0.001)],
        genes=genes,
        replicons=(("chr1", 100_000, "circular"), ("p1", 5_000, "circular")),
    )
    assert float(r.row("GSIM_00200")["coverage_frac"]) == 0.0
    assert float(r.row(G)["coverage_frac"]) == pytest.approx(0.1)


# ------------------------------------------------------------------------------- empty inputs


def _assert_all_zero(r: Run) -> None:
    for tag in r.genes["locus_tag"]:
        row = r.row(tag)
        for c in (
            "n_peptides_unique",
            "n_peptides_multi_identical",
            "n_peptides_shared",
            "n_peptides_tolerant",
            "n_accepted",
        ):
            assert int(row[c]) == 0
        assert float(row["coverage_frac"]) == 0.0
        assert _is_na(row["best_peptide_q"])
        assert bool(row["detected_2pep"]) is False
        assert flags_of(row) == set()


def test_proteins_only_no_c9_header_only_c11() -> None:
    r = run([], None)
    assert r.accepted.empty
    assert set(eo.ACCEPTED_COLUMNS) <= set(r.accepted.columns)
    _assert_all_zero(r)


def test_header_only_c11_with_c9() -> None:
    _assert_all_zero(run([], [pep(P1, 0.001)]))


def test_no_genes() -> None:
    r = run([], None, genes=[])
    assert len(r.evidence) == 0
    assert set(eo.EVIDENCE_COLUMNS) <= set(r.evidence.columns)


def test_evidence_columns_present() -> None:
    r = run([EX1_HIT], [pep(P1, 0.001)])
    assert set(eo.EVIDENCE_COLUMNS) <= set(r.evidence.columns)


# ------------------------------------------------------------------------------- determinism


def _scenario() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    hits = [
        *EXAMPLES["ex2_two_peptides"].hits,
        *EXAMPLES["ex6_multi_identical"].hits[1:],
        hit("SHAREDPEPK", G, 1200, uniqueness="shared"),
        hit("SHAREDPEPK", H, 9100, uniqueness="shared"),
        hit("TOLERANTK", H, 9150, strand="-", route="tolerant", placement_q=0.002),
    ]
    peps = [pep(P1, 0.001), pep(P2, 0.004), pep("SHAREDPEPK", 0.001), pep("TOLERANTK", 0.003)]
    return hits, peps


def test_deterministic_and_row_order_free() -> None:
    hits, peps = _scenario()
    a = run(hits, peps)
    b = run(hits, peps)
    pd.testing.assert_frame_equal(a.evidence, b.evidence)
    pd.testing.assert_frame_equal(a.accepted, b.accepted)
    c = run(hits[::-1], peps[::-1])
    for tag in a.genes["locus_tag"]:
        ra, rc = a.row(tag), c.row(tag)
        assert flags_of(ra) == flags_of(rc)
        for col in eo.EVIDENCE_COLUMNS[:-1]:
            assert (_is_na(ra[col]) and _is_na(rc[col])) or ra[col] == rc[col], (tag, col)
    key = ["gene", "peptide"]
    pd.testing.assert_frame_equal(
        a.accepted.sort_values(key).reset_index(drop=True)[list(eo.ACCEPTED_COLUMNS)],
        c.accepted.sort_values(key).reset_index(drop=True)[list(eo.ACCEPTED_COLUMNS)],
        check_dtype=False,
    )


def test_inputs_not_mutated() -> None:
    ev = _ev()
    hits_rows, pep_rows = _scenario()
    hits, peps = eo.hits_frame(hits_rows), eo.peptides_frame(pep_rows)
    genes, reps = eo.genes_frame(GENES), eo.replicons_frame(REPLICONS)
    copies = [x.copy(deep=True) for x in (hits, peps, genes, reps)]
    acc = ev.accepted_peptides(hits, peps, _cfg())
    ev.gene_evidence(genes, hits, acc, reps)
    for x, y in zip((hits, peps, genes, reps), copies, strict=True):
        pd.testing.assert_frame_equal(x, y)


# ------------------------------------------------------------------------------- simulated fixtures


def _sim_inputs(n_genes: int, seed: int, origin: bool) -> tuple[Any, ...]:
    from giapsim.expression import ExpressionConfig, simulate_expression
    from giapsim.layout import LayoutConfig, simulate_layout

    layout = simulate_layout(
        LayoutConfig(
            n_genes_chromosome=n_genes, plasmids=(20,), linear_contigs=(15,), origin_crossing=origin
        ),
        np.random.default_rng(seed),
    )
    expr = simulate_expression(layout, ExpressionConfig(), np.random.default_rng(seed + 1))
    detected = dict(
        zip(expr.genes["locus_tag"], expr.genes["detected_true"].astype(bool), strict=True)
    )
    peps, hits = eo.perfect_evidence(
        layout.genes, layout.replicons, detected, np.random.default_rng(seed + 2)
    )
    return layout, detected, peps, hits


@pytest.mark.parametrize("origin", [True, False])
def test_simulated_layout_matches_oracle(origin: bool) -> None:
    layout, detected, peps, hits = _sim_inputs(400, 11, origin)
    ev = _ev()
    acc = ev.accepted_peptides(hits, peps, _cfg())
    out = ev.gene_evidence(layout.genes, hits, acc, layout.replicons)
    r = Run(acc, out, layout.genes)
    assert_matches_oracle(
        r, eo.gene_evidence(layout.genes, hits, peps, layout.replicons, eo.Rules())
    )
    for tag, det in detected.items():
        if det:
            assert int(r.row(tag)["n_accepted"]) >= 1, tag
    if origin:
        crossing = layout.genes.loc[layout.genes["crosses_origin"].astype(bool), "locus_tag"]
        assert len(crossing) >= 1


def test_simulated_4500_genes_detected_have_evidence() -> None:
    layout, detected, peps, hits = _sim_inputs(4500, 5, True)
    ev = _ev()
    acc = ev.accepted_peptides(hits, peps, _cfg())
    out = ev.gene_evidence(layout.genes, hits, acc, layout.replicons)
    assert len(out) == len(layout.genes) >= 4500
    r = Run(acc, out, layout.genes)
    missing = [t for t, d in detected.items() if d and int(r.row(t)["n_accepted"]) < 1]
    assert not missing
    # The M4-11 time budget (best of 3 < 1 s for 4,500 genes) is checked in
    # benchmarks/test_genetable_evidence_bench.py (wall clock is flaky on loaded nodes).


def test_m409_default_fixture_detected_true_have_evidence() -> None:
    """Acceptance criterion on the M4-09 default fixture (skipped until M4-09 merges)."""
    evmod = pytest.importorskip("giapsim.evidence")
    from giapsim.expression import ExpressionConfig, simulate_expression
    from giapsim.layout import LayoutConfig, simulate_layout

    layout = simulate_layout(LayoutConfig(), np.random.default_rng(1))
    expr = simulate_expression(layout, ExpressionConfig(), np.random.default_rng(2))
    sim = evmod.simulate_evidence(layout, expr, evmod.EvidenceConfig(), np.random.default_rng(3))
    ev = _ev()
    acc = ev.accepted_peptides(sim.peptide_hits, sim.peptides, _cfg())
    out = ev.gene_evidence(layout.genes, sim.peptide_hits, acc, layout.replicons)
    r = Run(acc, out, layout.genes)
    det = expr.genes.loc[expr.genes["detected_true"].astype(bool), "locus_tag"]
    assert all(int(r.row(t)["n_accepted"]) >= 1 for t in det)
    # The M4-11 time budget (< 1 s per 4,500 genes) on this fixture is checked in
    # benchmarks/test_genetable_evidence_bench.py (wall clock is flaky on loaded nodes).
