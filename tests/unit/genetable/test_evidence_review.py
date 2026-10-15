"""M4-11 review pins (docs/specs/M4-11-accepted-evidence.md, lead decisions 5-6, PR #74 review).

Decision 5: ``coverage_frac`` is the union of ALL accepted placement intervals of a gene (repeats of
one peptide inside one gene included); the one-row-per-(gene, peptide) rule of decision 4 applies to
the accepted rows and the counts only.
Decision 6: nested intervals in the union; ``cterm_extension`` not abundance-eligible; contaminant
evidence carried to the gene; representative row exact-first then leftmost; without C9 the C11
``best_q_input`` is the engine q; an internal repeat gives union coverage.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import evidence_oracle as eo
from evidence_oracle import hit, pep

G = "GSIM_00005"  # [1000, 1300) on chr1, 300 bp
GENES = [eo.gene(G, 1000, 1300), eo.gene("GSIM_00010", 5000, 5300)]
P1 = "PEPTIDEK"  # 8 aa -> 24 bp by default


def _cfg(q: float = 0.01, relation: str = "related") -> Any:
    from giap.config import GiapConfig

    return GiapConfig(
        genome_relation=relation, input={"q_value_max": q}, stats={"alpha": {"placement": 0.01}}
    )


def run(
    hit_rows: list[dict[str, Any]], pep_rows: list[dict[str, Any]] | None
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any]]:
    """Accepted rows, the evidence row of G, and the oracle's expectation for G."""
    from giap.genetable import evidence as ev

    genes = eo.genes_frame(GENES)
    hits = eo.hits_frame(hit_rows)
    peps = eo.peptides_frame(pep_rows) if pep_rows is not None else None
    reps = eo.replicons_frame()
    acc = ev.accepted_peptides(hits, peps, _cfg())
    out = ev.gene_evidence(genes, hits, acc, reps)
    row = out[out["locus_tag"] == G].iloc[0].to_dict()
    exp = eo.gene_evidence(genes, hits, peps, reps, eo.Rules())[G]
    return acc, row, exp


# ---------------------------------------------------------------- decision 5: union coverage


def test_internal_repeat_gives_union_coverage() -> None:
    # P1 placed twice inside G (internal repeat, both accepted): 24 + 24 bp.
    acc, row, exp = run([hit(P1, G, 1000), hit(P1, G, 1100)], [pep(P1, 0.001)])
    assert len(acc) == 1, "one accepted row per (gene, peptide)"
    assert int(row["n_accepted"]) == 1 and int(row["n_peptides_unique"]) == 1
    assert exp["coverage_frac"] == pytest.approx(48 / 300)
    assert float(row["coverage_frac"]) == pytest.approx(48 / 300)


def test_union_covers_tolerant_and_exact_placements_of_one_peptide() -> None:
    # The representative row is the exact one, but the tolerant placement still adds coverage.
    acc, row, exp = run(
        [hit(P1, G, 1000, route="tolerant", placement_q=0.004), hit(P1, G, 1200)],
        [pep(P1, 0.001)],
    )
    assert len(acc) == 1
    assert exp["coverage_frac"] == pytest.approx(48 / 300)
    assert float(row["coverage_frac"]) == pytest.approx(48 / 300)


def test_rejected_repeat_placement_adds_no_coverage() -> None:
    # Only accepted placements count: a second placement with placement_q 0.03 adds nothing.
    _, row, exp = run([hit(P1, G, 1000), hit(P1, G, 1100, placement_q=0.03)], [pep(P1, 0.001)])
    assert exp["coverage_frac"] == pytest.approx(24 / 300)
    assert float(row["coverage_frac"]) == pytest.approx(24 / 300)


# ---------------------------------------------------------------- decision 6


def test_nested_intervals_in_coverage_union() -> None:
    # Gene-relative [0, 90), [10, 20), [30, 40): the union is 90 bp, not 90 + 10.
    rows = [
        hit("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAK", G, 1000, 1090),
        hit("CCCK", G, 1010, 1020),
        hit("DDDK", G, 1030, 1040),
    ]
    peps = [pep(r["peptide"], 0.001) for r in rows]
    _, row, exp = run(rows, peps)
    assert exp["coverage_frac"] == pytest.approx(90 / 300)
    assert float(row["coverage_frac"]) == pytest.approx(90 / 300)


def test_cterm_extension_not_abundance_eligible() -> None:
    acc, row, _ = run([hit(P1, G, 1290, locus_class="cterm_extension")], [pep(P1, 0.001)])
    assert acc["locus_class"].tolist() == ["cterm_extension"]
    assert acc["abundance_eligible"].astype(bool).tolist() == [False]
    assert int(row["n_accepted"]) == 1


def test_contaminant_evidence_and_is_contaminant() -> None:
    acc, row, _ = run(
        [hit(P1, G, 1030)], [pep(P1, 0.001, is_contaminant=True, protein_ids=["CONT_P02769"])]
    )
    assert acc["is_contaminant"].astype(bool).tolist() == [True]
    assert bool(row["contaminant_evidence"]) is True
    assert "contaminant" not in {str(f) for f in row["evidence_flags"]}

    acc, row, _ = run([hit(P1, G, 1030)], [pep(P1, 0.001)])
    assert acc["is_contaminant"].astype(bool).tolist() == [False]
    assert bool(row["contaminant_evidence"]) is False


def test_row_choice_exact_before_leftmost_tolerant() -> None:
    acc, _, _ = run(
        [hit(P1, G, 1000, route="tolerant", placement_q=0.004), hit(P1, G, 1200)],
        [pep(P1, 0.001)],
    )
    r = acc.iloc[0]
    assert (r["route"], int(r["start"]), int(r["end"])) == ("exact", 1200, 1224)


def test_row_choice_leftmost_of_two_exact() -> None:
    acc, _, _ = run([hit(P1, G, 1200), hit(P1, G, 1100)], [pep(P1, 0.001)])
    r = acc.iloc[0]
    assert (r["route"], int(r["start"]), int(r["end"])) == ("exact", 1100, 1124)


def test_without_c9_best_q_input_above_cutoff_rejected() -> None:
    acc, row, exp = run([hit(P1, G, 1030, best_q=0.02)], None)
    assert acc.empty
    assert int(row["n_accepted"]) == 0 and exp["n_accepted"] == 0

    acc, row, _ = run([hit(P1, G, 1030, best_q=0.004)], None)
    assert float(acc.iloc[0]["q_value"]) == pytest.approx(0.004)
    assert float(row["best_peptide_q"]) == pytest.approx(0.004)
