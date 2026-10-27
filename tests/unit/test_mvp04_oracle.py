"""The MVP-04 gene-level oracle (tests/helpers/mvp04_decoys.py ``gene_evidence`` /
``gene_competition``) on a hand-built view (MVP-04 lead decision 11; statistician review of PR
#247, B3(b)): the ``GIPPAORF_nterm_extension_*`` channel and ``multi_identical`` acceptance.

Rules (as ``gippa.genetable.detection`` and M7-12 decision 28):

* a six-frame null on an ``nterm_extension`` entry gives target evidence to the one eligible gene
  the entry overlaps on its strand; its ``rev_GIPPAORF_nterm_extension_*`` decoy pairs with that
  gene only when the entry was voted (the peptide has an accepted target row);
* a protein-entry null found (I = L) in several loci gives evidence to all of them when their
  proteins are identical (``multi_identical``) and to none when they differ (``shared``);
* correct peptides: one locus, or several identical ones;
* other six-frame nulls (and their decoys) reach no gene; rows above ``q_value_max`` are no
  evidence; a tie goes to the decoy.
"""

from __future__ import annotations

from pathlib import Path

from helpers import mvp04_data as M
from helpers import mvp04_decoys as D

C_ = M.CONTIG
NT_E = f"GIPPAORF_nterm_extension_{C_}_4901_5300_p"  # overlaps gene E (+, [5000, 5600))
NT_F = f"GIPPAORF_nterm_extension_{C_}_6901_7300_p"  # overlaps gene F (+, [7000, 7600))
AS = f"GIPPAORF_antisense_{C_}_9001_9150_m"

PARALOG = "MSTLKWEEPLQGNAVRDDLK"
GENES = [
    M.Gene("A", "cds", 1000, 1300, "+", "MAAAPEPTIDEKAAAK", None, False, False),
    M.Gene("B", "cds", 2000, 2063, "+", PARALOG, None, False, False),
    M.Gene("C", "cds", 3000, 3063, "-", PARALOG, None, False, False),
    M.Gene("D", "cds", 4000, 4063, "+", "MQQQKLQGNAVRQQQ", None, False, False),
    M.Gene("E", "cds", 5000, 5600, "+", "MEEEEEEEEK", None, False, False),
    M.Gene("F", "cds", 7000, 7600, "+", "MFFFFFFFFK", None, False, False),
    M.Gene("T", "cds", 8000, 8063, "+", "MTRUNCATEDK", None, True, False),  # not eligible
]
PROTEINS = {g.locus: [g.aa] for g in GENES}


def _row(pep: str, pid: str, q: float, score: float, decoy: bool, run: str = "r1"
         ) -> dict[str, str]:  # fmt: skip
    return {"peptide": pep, "protein_ids": pid, "q_value": str(q), "pep": str(q * 10),
            "score": str(score), "is_decoy": "true" if decoy else "false", "intensity": "NA",
            "run": run}  # fmt: skip


def _view() -> D.View:
    peps = [
        _row("APEPTIDEK", "A", 0.001, 9.0, False),  # correct, unique to A
        _row("WEEPLQGNAVR", "B", 0.002, 8.0, False),  # protein null in B and C (identical)
        _row("WEEPLQGNAVR"[::-1], "rev_B", 0.003, 7.0, True, "r2"),  # its decoy
        _row("LQGNAVR", "B", 0.002, 8.5, False),  # protein null in B, C and D (shared)
        _row("NTERMNULLK", NT_E, 0.004, 6.0, False),  # nterm null on E's entry, accepted
        _row("NTERMNULLK"[::-1], "rev_" + NT_E, 0.002, 6.5, True, "r2"),  # voted -> pairs with E
        _row("NTERMDECOYONLYK"[::-1], "rev_" + NT_F, 0.001, 9.5, True),  # F: never voted
        _row("ANTISENSEK", AS, 0.001, 9.9, False),  # antisense null: no gene
        _row("ANTISENSEK"[::-1], "rev_" + AS, 0.001, 9.9, True, "r2"),
        _row("TRUNCATEDK", "T", 0.001, 9.0, False),  # correct, but T is not eligible
        _row("AAAK", "A", 0.5, 1.0, False),  # above q_value_max: no evidence
    ]
    truth = [
        {"peptide": "APEPTIDEK", "is_wrong": "false", "origin": "cds", "locus_tag": "A"},
        {"peptide": "TRUNCATEDK", "is_wrong": "false", "origin": "cds", "locus_tag": "T"},
        {"peptide": "AAAK", "is_wrong": "false", "origin": "cds", "locus_tag": "A"},
        {"peptide": "WEEPLQGNAVR", "is_wrong": "true", "origin": D.PROTEIN_NULL_ORIGIN,
         "locus_tag": "B"},
        {"peptide": "LQGNAVR", "is_wrong": "true", "origin": D.PROTEIN_NULL_ORIGIN,
         "locus_tag": "B"},
        {"peptide": "NTERMNULLK", "is_wrong": "true", "origin": "sixframe_noncoding",
         "locus_tag": "NA"},
        {"peptide": "NTERMDECOYONLYK", "is_wrong": "true", "origin": "sixframe_noncoding",
         "locus_tag": "NA"},
        {"peptide": "ANTISENSEK", "is_wrong": "true", "origin": "sixframe_noncoding",
         "locus_tag": "NA"},
    ]  # fmt: skip
    nulls: dict[str, D.Null] = {}
    for t in truth:
        if t["is_wrong"] == "true":
            kind = "sixframe" if t["origin"] == "sixframe_noncoding" else "protein"
            nulls[t["peptide"]] = D.Null(t["peptide"], kind, t["origin"], truth=[t])
    rev = {p[::-1]: p for p in nulls}
    for r in peps:
        if r["is_decoy"] == "false" and r["peptide"] in nulls:
            nulls[r["peptide"]].targets[r["run"]] = r
        elif r["is_decoy"] == "true" and r["peptide"] in rev:
            nulls[rev[r["peptide"]]].decoys[r["run"]] = r
    correct = {"APEPTIDEK", "TRUNCATEDK", "AAAK"}
    loci = {"APEPTIDEK": {"A"}, "TRUNCATEDK": {"T"}, "AAAK": {"A"}}
    return D.View(Path("."), peps, truth, {}, {}, ["r1", "r2"], nulls, correct, loci)


def _by_gene() -> dict[str, set[tuple[bool, str, str]]]:
    out: dict[str, set[tuple[bool, str, str]]] = {}
    for e in D.gene_evidence(_view(), GENES, PROTEINS):
        out.setdefault(e.gene, set()).add((e.is_decoy, e.channel, e.peptide))
    return out


def test_nterm_gene_is_the_one_overlapping_eligible_gene_on_the_strand() -> None:
    assert D.nterm_gene(NT_E, GENES) == "E"
    assert D.nterm_gene(f"GIPPAORF_nterm_extension_{C_}_4901_5300_m", GENES) is None
    assert D.nterm_gene(f"GIPPAORF_nterm_extension_{C_}_7901_8100_p", GENES) is None  # T only


def test_nterm_null_target_and_voted_decoy_reach_the_gene() -> None:
    ev = _by_gene()
    assert ev["E"] == {(False, "nterm_null", "NTERMNULLK"),
                       (True, "nterm_null", "NTERMNULLK"[::-1])}  # fmt: skip


def test_unvoted_nterm_decoy_leaves_pairing() -> None:
    assert "F" not in _by_gene()


def test_multi_identical_null_reaches_every_identical_locus_and_shared_none() -> None:
    ev = _by_gene()
    assert ev["B"] == {(False, "protein_null", "WEEPLQGNAVR"),
                       (True, "protein_null", "WEEPLQGNAVR"[::-1])}  # fmt: skip
    assert ev["C"] == {(False, "protein_null", "WEEPLQGNAVR")}
    assert "D" not in ev  # LQGNAVR is shared (B, C, D differ): no evidence anywhere


def test_correct_ineligible_antisense_and_unaccepted_rows() -> None:
    ev = _by_gene()
    assert ev["A"] == {(False, "correct", "APEPTIDEK")}
    assert "T" not in ev
    assert set(ev) == {"A", "B", "C", "E"}


def test_competition_uses_is_wrong_false_truth_and_ties_go_to_the_decoy() -> None:
    comp = {c.gene: c for c in D.gene_competition(_view(), GENES, PROTEINS)}
    assert set(comp) == {"A", "B", "C", "E"}
    assert comp["A"] == D.Competition("A", True, True)
    # B: null target q 0.002 beats its decoy q 0.003; it has a null truth row, but no true one
    assert comp["B"] == D.Competition("B", False, True)
    assert comp["C"] == D.Competition("C", False, True)
    # E: the voted decoy (q 0.002) beats the target (q 0.004)
    assert comp["E"] == D.Competition("E", False, False)
