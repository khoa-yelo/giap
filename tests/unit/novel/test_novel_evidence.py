"""Novel-ORF class evidence, RBS support, promoter/terminator evidence and ranking (issue M7-09;
docs/specs/M7-09-novel-evidence-and-ranking.md; GIAP_SPEC.md §7.6.3-§7.6.6, §9.4; ADR-0012,
ADR-0013).

Examples 1-12 and every degenerate row of the spec, its lead decisions (RBS quantile needs 20
detected genes; pseudogene junction = first premature in-frame stop; ``cds_outframe`` needs 2
frame-specific peptides; ``cterm_extension`` is labelled ``readthrough_candidate`` and never
called; only called ORFs are ranked) and the ADR-0012 consequences.

Readings pinned by the test author (see tests/helpers/novel_evidence_fixture.py for the
``apply_evidence`` signature):

* ``class_evidence(row, cfg, rbs_min)`` returns ``(class rule of rule 4, labels)``; it does **not**
  include the error-control rule 8 (``called`` = rule 8 and the class rule is ``apply_evidence``'s
  job, rule 9), so with decoys (M7-11) only rule 8 changes. ``row`` holds the ``novel_orfs``
  cells plus ``start_upstream_nt`` and ``spans_annotated_stop`` (fixture ``row``). Its labels are
  ``rbs_supported`` (any class, when the chosen start passes), ``rbs_threshold_unestimated``
  (``nterm_extension``/``intergenic`` when ``rbs_min`` is ``None``) and ``readthrough_candidate``.
* the ``cds_outframe``/``antisense`` "≥ 2 distinct peptides" is a fixed 2 (decision 3), not
  ``min_peptides_without_decoys``.
* "codons from the chosen start to the stop" = ``length_aa_best``; "a start codon exists" =
  ``start_longest`` not NA.
* ``terminator_downstream(orf, terminators, rep)`` -> ``(conf, distance)``; distance = 5' end of
  the terminator (its first base in transcription direction) minus the last base of the stop
  codon, so example 9 gives 77 and the closed window is distance 0..300. Equal best confidences
  -> the nearer call.
* ``pseudogene_junction(gene, genome, rep)`` -> the lowest 0-based coordinate (wrapped) of the first
  in-frame stop codon from the annotated start that is not the gene's last codon; ``None`` if none.
* ``junction_spanned``: some peptide interval contains the whole junction codon;
  ``homopolymer_junction``: a maximal run of >= 5 identical bases (forward strand) overlaps
  ``[junction - 5, junction + 8)``. A pseudogene without a junction (or without its gene) gets
  ``junction_unknown``.
* ``promoter_supported`` is NA when ``promoter_score`` is NA.
* ``flags`` are ``;``-joined in FLAG_ORDER (§9.4 Novel ORFs, then the schema-pattern order).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.novel_evidence_fixture import (
    FLAG_ORDER,
    GT15,
    GT30,
    ID1,
    ID6,
    PSEUDO_JUNCTION,
    THRESHOLD_30,
    TIRS_30,
    cfg,
    circ,
    ev,
    evidence_table,
    ex1_case,
    ex1_table,
    flags_of,
    gene_table,
    is_na,
    lin,
    nterm2,
    pseudo_case,
    replace_orf,
    row,
    term,
)
from helpers.novel_oracle import cds
from helpers.novel_outputs_fixture import NOVEL_ORFS_COLUMNS, orfs_of, scored, table, two_contigs

from giap.model.coords import Interval
from giap.novel.candidates import Orf

SCHEMA = Path(__file__).resolve().parents[3] / "src/giap/schemas/tables/novel_orfs.schema.json"


def _orf(strand: str, stop: int, contig: str = "c1") -> Orf:
    return Orf(f"NOVEL_{contig}_1_{'p' if strand == '+' else 'm'}", "intergenic", contig,
               strand, 1 if strand == "+" else -1, stop, (), None, False, False,  # type: ignore[arg-type]
               ("AAAAAAA", "CCCCCCC"), 0.001, None)  # fmt: skip


# ---------------------------------------------------------------- rule 2: RBS threshold


def test_rbs_threshold_spec_example_is_390() -> None:
    assert ev().rbs_threshold(GT30, 0.10) == pytest.approx(THRESHOLD_30)


def test_rbs_threshold_is_numpy_linear_quantile() -> None:
    rng = np.random.default_rng(3)
    tirs = list(rng.uniform(1, 9000, 57))
    for q in (0.0, 0.1, 0.37, 0.5, 1.0):
        want = float(np.quantile(np.array(tirs), q, method="linear"))
        assert ev().rbs_threshold(gene_table(tirs), q) == pytest.approx(want), q


def test_rbs_threshold_quantile_zero_is_the_minimum() -> None:
    """Degenerate row "rbs_min_quantile = 0"."""
    assert ev().rbs_threshold(GT30, 0.0) == pytest.approx(100.0)


def test_rbs_threshold_only_eligible_detected_scored_genes() -> None:
    """Partial, undetected, NA-detected and NA-TIR rows never enter the quantile."""
    extra = [(1.0, True, True), (2.0, False, False), (3.0, None, False), (None, True, False)]
    gt = gene_table(TIRS_30, extra=extra)
    assert ev().rbs_threshold(gt, 0.10) == pytest.approx(THRESHOLD_30)
    assert ev().rbs_threshold(gt, 0.0) == pytest.approx(100.0)


@pytest.mark.parametrize(("n", "estimated"), [(15, False), (19, False), (20, True), (30, True)])
def test_rbs_threshold_needs_20_detected_genes(n: int, estimated: bool) -> None:
    """Decision 1: fewer than 20 eligible detected genes with a score -> ``None``."""
    got = ev().rbs_threshold(gene_table(TIRS_30[:n]), 0.10)
    if estimated:
        assert got == pytest.approx(float(np.quantile(np.array(TIRS_30[:n]), 0.10)))
    else:
        assert got is None


def test_rbs_threshold_undetected_genes_do_not_count_towards_20() -> None:
    gt = gene_table(TIRS_30[:19], extra=[(5000.0, False, False)] * 10)
    assert ev().rbs_threshold(gt, 0.10) is None


def test_rbs_threshold_min_n_parameter() -> None:
    assert ev().rbs_threshold(GT15, 0.10, min_n=15) == pytest.approx(
        float(np.quantile(np.array(TIRS_30[:15]), 0.10))
    )
    assert ev().rbs_threshold(GT15, 0.10, min_n=16) is None


# ---------------------------------------------------------------- rule 4: class_evidence


def _ce(r: dict[str, Any], rbs_min: float | None = THRESHOLD_30, **c: Any) -> tuple[bool, list]:
    ok, labels = ev().class_evidence(r, cfg(**c), rbs_min)
    return bool(ok), list(labels)


def test_example_1_intergenic_passes_with_rbs_support() -> None:
    assert _ce(row()) == (True, ["rbs_supported"])


def test_example_2_intergenic_tir_below_threshold_fails() -> None:
    assert _ce(row(ostir_tir=300.0)) == (False, [])


def test_rbs_threshold_is_inclusive() -> None:
    assert _ce(row(ostir_tir=THRESHOLD_30)) == (True, ["rbs_supported"])


def test_example_3_intergenic_of_9_codons_fails() -> None:
    assert _ce(row(length_aa_best=9, length_aa_longest=15))[0] is False
    assert _ce(row(length_aa_best=10))[0] is True


def test_min_codons_is_read_from_config() -> None:
    assert _ce(row(), min_codons=21)[0] is False
    assert _ce(row(), min_codons=20)[0] is True


def test_intergenic_needs_a_start_codon() -> None:
    r = row(start_longest=None, start_best_ostir=None, start_codon_best=None, ostir_tir=None)
    assert _ce(r)[0] is False


def test_example_4_class_rule_ignores_peptide_count() -> None:
    """Example 4 fails on rule 8, not on the class rule (rule 9 composes them)."""
    assert _ce(row(n_peptides=1, ostir_tir=5000.0)) == (True, ["rbs_supported"])


def test_ostir_missing_at_chosen_start_fails_rbs_classes() -> None:
    """Degenerate row "OSTIR missing at the chosen start": ``rbs_supported`` NA -> not called."""
    r = row(start_best_ostir=None, ostir_tir=None, start_codon_best="ATG", length_aa_best=26)
    assert _ce(r) == (False, [])
    n = row(**{"class": "nterm_extension"}, start_best_ostir=None, ostir_tir=None,
            start_upstream_nt=51, overlapping_gene="G")  # fmt: skip
    assert _ce(n) == (False, [])


def test_example_5_cds_outframe_two_peptides_no_ostir_passes() -> None:
    r = row(**{"class": "cds_outframe"}, start_best_ostir=None, ostir_tir=None,
            overlapping_gene="G")  # fmt: skip
    assert _ce(r) == (True, [])
    assert _ce(r, rbs_min=None) == (True, [])  # no RBS rule -> no rbs_threshold_unestimated


@pytest.mark.parametrize("cls", ["cds_outframe", "antisense"])
def test_two_distinct_peptides_rule_is_fixed_at_two(cls: str) -> None:
    """Decision 3: the class rule is >= 2 distinct peptides whatever the error-control minimum."""
    one = row(**{"class": cls}, n_peptides=1, peptides="ALEKRDSPHW")
    assert _ce(one)[0] is False
    assert _ce(one, min_peptides_without_decoys=1)[0] is False
    assert _ce(row(**{"class": cls}, n_peptides=2))[0] is True
    assert _ce(row(**{"class": cls}, n_peptides=3), min_peptides_without_decoys=5)[0] is True


def test_example_6_antisense_one_peptide_fails() -> None:
    assert _ce(row(**{"class": "antisense"}, n_peptides=1, ostir_tir=None))[0] is False


def test_example_7_nterm_extension_upstream_ttg_passes() -> None:
    r = row(**{"class": "nterm_extension"}, start_longest=249, start_best_ostir=249,
            start_codon_best="TTG", ostir_tir=800.0, start_upstream_nt=51, overlapping_gene="G",
            length_aa_best=216, length_aa_longest=216)  # fmt: skip
    assert _ce(r) == (True, ["rbs_supported"])
    assert _ce({**r, "ostir_tir": 300.0}) == (False, [])


@pytest.mark.parametrize("upstream", [0, -9, None])
def test_nterm_extension_without_an_upstream_start_fails(upstream: int | None) -> None:
    """Degenerate row: ``start_longest`` is the annotated start itself (0), downstream of it, or
    absent."""
    r = row(**{"class": "nterm_extension"}, ostir_tir=5000.0, start_upstream_nt=upstream,
            overlapping_gene="G")  # fmt: skip
    if upstream is None:
        r.update(start_longest=None, start_best_ostir=None, start_codon_best=None)
    assert _ce(r)[0] is False


@pytest.mark.parametrize("cls", ["nterm_extension", "intergenic"])
def test_example_10_threshold_unestimated_flags_rbs_classes(cls: str) -> None:
    r = row(**{"class": cls}, start_upstream_nt=51, overlapping_gene="G", ostir_tir=5000.0)
    assert _ce(r, rbs_min=None) == (False, ["rbs_threshold_unestimated"])


@pytest.mark.parametrize("cls", ["cds_outframe", "antisense", "pseudogene", "sorf"])
def test_rbs_supported_label_on_any_class(cls: str) -> None:
    """Rule 3: ``rbs_supported`` is set for every ORF whose chosen start passes."""
    assert "rbs_supported" in _ce(row(**{"class": cls}, ostir_tir=5000.0))[1]
    assert "rbs_supported" not in _ce(row(**{"class": cls}, ostir_tir=10.0))[1]


@pytest.mark.parametrize("cls", ["pseudogene", "sorf"])
def test_pseudogene_and_sorf_have_no_extra_rule(cls: str) -> None:
    r = row(**{"class": cls}, n_peptides=1, start_longest=None, start_best_ostir=None,
            start_codon_best=None, ostir_tir=None, length_aa_best=3,
            overlapping_gene="P")  # fmt: skip
    assert _ce(r) == (True, [])
    assert _ce(r, rbs_min=None) == (True, [])


def test_cterm_extension_passes_only_when_a_peptide_spans_the_stop() -> None:
    """M7-16 rule 12 (replaces decision 4's ``ok = False``): a peptide spanning the annotated stop
    passes the class rule, labelled ``readthrough_candidate``; without one it does not pass."""
    r = row(**{"class": "cterm_extension"}, n_peptides=3, ostir_tir=9000.0, overlapping_gene="G")
    ok, labels = _ce({**r, "spans_annotated_stop": True})
    assert ok is True and "readthrough_candidate" in labels
    ok, labels = _ce({**r, "spans_annotated_stop": False})
    assert ok is False and "readthrough_candidate" not in labels


def test_class_evidence_labels_in_flag_order() -> None:
    r = row(**{"class": "cterm_extension"}, ostir_tir=9000.0, spans_annotated_stop=True)
    labels = _ce(r)[1]
    assert labels == sorted(labels, key=FLAG_ORDER.index)
    assert labels == ["rbs_supported", "readthrough_candidate"]


# ---------------------------------------------------------------- rule 7: terminators


def test_example_9_best_downstream_terminator_same_strand() -> None:
    o = _orf("+", 120)  # stop codon [120, 123): 1-based last base 123
    terms = [term("c1", 199, 230, "+", 0.8), term("c1", 449, 480, "+", 0.95),
             term("c1", 120, 150, "-", 0.99)]  # fmt: skip
    conf, dist = ev().terminator_downstream(o, terms, circ())
    assert conf == pytest.approx(0.8) and dist == 77


def test_example_9_minus_strand_mirror() -> None:
    o = _orf("-", 877)  # stop codon [877, 880): last base in transcription direction 877
    terms = [term("c1", 770, 801, "-", 0.8), term("c1", 520, 551, "-", 0.95),
             term("c1", 850, 880, "+", 0.99)]  # fmt: skip
    conf, dist = ev().terminator_downstream(o, terms, circ())
    assert conf == pytest.approx(0.8) and dist == 77


@pytest.mark.parametrize(("start", "found"), [(422, True), (423, False), (122, True), (300, True)])
def test_terminator_window_is_closed_0_to_300(start: int, found: bool) -> None:
    """Degenerate rows "exactly 300 nt after the stop" (found) and "301 nt after" (not); the 5'
    end at the stop's last base is distance 0 (ADR-0012 "0-300 nt")."""
    got = ev().terminator_downstream(_orf("+", 120), [term("c1", start, start + 25, "+", 0.7)],
                                     circ())  # fmt: skip
    assert got == ((pytest.approx(0.7), start - 122) if found else (None, None))


def test_terminator_upstream_of_the_stop_is_not_found() -> None:
    got = ev().terminator_downstream(_orf("+", 120), [term("c1", 60, 90, "+", 0.9)], circ())
    assert got == (None, None)


def test_terminator_on_the_other_strand_is_not_found() -> None:
    got = ev().terminator_downstream(_orf("+", 120), [term("c1", 199, 230, "-", 0.9)], circ())
    assert got == (None, None)


def test_terminator_no_calls_gives_none() -> None:
    assert ev().terminator_downstream(_orf("+", 120), [], circ()) == (None, None)


@pytest.mark.parametrize(
    ("strand", "stop", "iv", "dist"),
    [
        ("+", 990, (50, 80), 58),  # 5' end 50 after the origin: 1000 - 992 + 50
        ("+", 990, (995, 1010), 3),  # stored unwrapped across the origin
        ("+", 999, (10, 40), 9),  # stop codon [999, 1002) itself crosses the origin
        ("-", 5, (920, 951), 55),  # 5' end 950 (highest base): 5 + 1000 - 950
        ("-", 5, (990, 1003), 3),  # stored unwrapped, 5' end 1002 -> 2
    ],
)
def test_terminator_across_the_origin_is_found(
    strand: str, stop: int, iv: tuple[int, int], dist: int
) -> None:
    got = ev().terminator_downstream(_orf(strand, stop), [term("c1", *iv, strand, 0.6)], circ())
    assert got == (pytest.approx(0.6), dist)


def test_terminator_linear_contig_clips_without_wrap() -> None:
    """Degenerate row "linear contig, window running off the end": clipped; no wrap."""
    o = _orf("+", 990)
    assert ev().terminator_downstream(o, [term("c1", 50, 80, "+", 0.9)], lin()) == (None, None)
    got = ev().terminator_downstream(o, [term("c1", 995, 1000, "+", 0.9)], lin())
    assert got == (pytest.approx(0.9), 3)


def test_terminator_best_confidence_ties_take_the_nearer() -> None:
    terms = [term("c1", 199, 230, "+", 0.8), term("c1", 152, 180, "+", 0.8)]
    assert ev().terminator_downstream(_orf("+", 120), terms, circ()) == (pytest.approx(0.8), 30)


def test_terminator_below_threshold_is_reported_but_not_downstream() -> None:
    """``terminator_conf`` is the best confidence whatever it is; ``terminator_downstream`` needs
    ``terminatornet_yes`` (0.5)."""
    t = ex1_table(terminators=[term("c1", 199, 230, "+", 0.3)])
    assert t["terminator_conf"] == pytest.approx(0.3) and t["terminator_distance_nt"] == 77
    assert t["terminator_downstream"] == False  # noqa: E712
    assert "terminator_supported" not in flags_of(t)


def test_apply_terminator_columns_example_9() -> None:
    terms = [term("c1", 199, 230, "+", 0.8), term("c1", 449, 480, "+", 0.95),
             term("c1", 120, 150, "-", 0.99)]  # fmt: skip
    r = ex1_table(terminators=terms)
    assert r["terminator_conf"] == pytest.approx(0.8)
    assert r["terminator_distance_nt"] == 77
    assert r["terminator_downstream"] == True  # noqa: E712
    assert "terminator_supported" in flags_of(r)


def test_apply_no_terminator_call_gives_na_and_false() -> None:
    r = ex1_table(terminators=[])
    assert is_na(r["terminator_conf"]) and is_na(r["terminator_distance_nt"])
    assert r["terminator_downstream"] == False  # noqa: E712
    assert "terminator_supported" not in flags_of(r)


def test_apply_terminator_layer_missing_gives_three_na() -> None:
    r = ex1_table(terminators=None)
    for col in ("terminator_conf", "terminator_distance_nt", "terminator_downstream"):
        assert is_na(r[col]), col
    assert r["called"] == True  # noqa: E712  (never gates)


# ---------------------------------------------------------------- rule 6: promoter


@pytest.mark.parametrize(("score", "supported"), [(0.5, True), (0.93, True), (0.49, False)])
def test_promoter_supported_at_promotech_yes(score: float, supported: bool) -> None:
    r = ex1_table(promoter={ID1: score})
    assert r["promoter_supported"] == supported
    assert ("promoter_supported" in flags_of(r)) is supported
    assert r["called"] == True  # noqa: E712


def test_promoter_score_na_gives_na_promoter_supported() -> None:
    r = ex1_table()
    assert is_na(r["promoter_supported"])
    assert "promoter_supported" not in flags_of(r)


# ---------------------------------------------------------------- rule 5: pseudogene junction


def _pseudo_gene() -> Any:
    return cds("P", "c1", 540, 900, "+", "pseudogene")


def test_pseudogene_junction_first_premature_stop_plus() -> None:
    seq = pseudo_case(None).genome["c1"]
    assert ev().pseudogene_junction(_pseudo_gene(), seq, circ()) == PSEUDO_JUNCTION


def test_pseudogene_junction_takes_the_first_of_two() -> None:
    from helpers.novel_oracle import Seq

    s = Seq(1000).fill(543, 897).put(540, "ATG").put(897, "TAA").put(600, "TGA").put(750, "TAG")
    assert ev().pseudogene_junction(_pseudo_gene(), s.text(), circ()) == 600


def test_pseudogene_junction_none_without_premature_stop() -> None:
    seq = pseudo_case(None, premature=False).genome["c1"]
    assert ev().pseudogene_junction(_pseudo_gene(), seq, circ()) is None


def test_pseudogene_junction_out_of_frame_stop_is_ignored() -> None:
    from helpers.novel_oracle import Seq

    s = Seq(1000).fill(543, 897).put(540, "ATG").put(897, "TAA").put(601, "TGA")
    assert ev().pseudogene_junction(_pseudo_gene(), s.text(), circ()) is None


def test_pseudogene_junction_minus_strand() -> None:
    """Mirror of the plus case: gene - [100, 460), its start codon [457, 460), internal stop at
    [397, 400) (in frame: 460 - 400 = 60), terminal stop [100, 103)."""
    from helpers.novel_oracle import Seq, revcomp

    plus = Seq(1000).fill(543, 897).put(540, "ATG").put(897, "TAA").put(600, "TGA").text()
    seq = revcomp(plus)  # [a, b) -> [1000 - b, 1000 - a)
    gene = cds("P", "c1", 100, 460, "-", "pseudogene")
    assert ev().pseudogene_junction(gene, seq, circ()) == 1000 - 603


def test_pseudogene_junction_across_the_origin_is_wrapped() -> None:
    from helpers.novel_oracle import Seq

    s = Seq(1000).fill(953, 1100).put(950, "ATG").put(1097, "TAA").put(1040, "TGA")
    gene = cds("P", "c1", 950, 1100, "+", "pseudogene")
    assert ev().pseudogene_junction(gene, s.text(), circ()) == 40


def test_example_8_junction_labels() -> None:
    """Example 8: internal TGA [600, 603), a hand-built peptide interval [590, 620) over it and
    ``AAAAA`` at [595, 600) -> ``junction_spanned;homopolymer_junction`` in this order."""
    case = pseudo_case((595, "AAAAA"))
    (orf,) = orfs_of(case)
    t = evidence_table(case, peptide_intervals={orf.orf_id: [Interval("c1", 590, 620, "+")]})
    (r,) = t.to_dict("records")
    f = flags_of(r)
    assert f[:2] == ["junction_spanned", "homopolymer_junction"]
    assert "junction_unknown" not in f
    assert r["called"] == True  # noqa: E712  (labels only; pseudogene has no extra rule)


@pytest.mark.parametrize(
    ("iv", "spanned"),
    [((590, 620), True), ((600, 603), True), ((570, 597), False), ((590, 602), False),
     ((601, 630), False)],
)  # fmt: skip
def test_junction_spanned_needs_the_whole_codon(iv: tuple[int, int], spanned: bool) -> None:
    case = pseudo_case(None)
    (orf,) = orfs_of(case)
    t = evidence_table(case, peptide_intervals={orf.orf_id: [Interval("c1", *iv, "+")]})
    assert ("junction_spanned" in flags_of(t.iloc[0])) is spanned


@pytest.mark.parametrize(
    ("run", "label"),
    [((595, "AAAAA"), True), ((590, "AAAAA"), False), ((596, "AAAA"), False),
     ((603, "TTTTT"), True), ((608, "TTTTT"), False), (None, False)],
)  # fmt: skip
def test_homopolymer_window_is_junction_minus_5_to_plus_8(
    run: tuple[int, str] | None, label: bool
) -> None:
    case = pseudo_case(run)
    t = evidence_table(case, peptide_intervals={})
    f = flags_of(t.iloc[0])
    assert ("homopolymer_junction" in f) is label
    assert "junction_spanned" not in f


def test_pseudogene_without_premature_stop_flags_junction_unknown() -> None:
    """Degenerate row: no premature stop -> ``junction_unknown``, no junction labels."""
    case = pseudo_case((595, "AAAAA"), premature=False)
    (orf,) = orfs_of(case)
    t = evidence_table(case, peptide_intervals={orf.orf_id: [Interval("c1", 590, 620, "+")]})
    f = flags_of(t.iloc[0])
    assert "junction_unknown" in f
    assert "junction_spanned" not in f and "homopolymer_junction" not in f


def test_non_pseudogene_rows_get_no_junction_flags() -> None:
    r = ex1_table()
    assert not {"junction_spanned", "homopolymer_junction", "junction_unknown"} & set(flags_of(r))


# ---------------------------------------------------------------- rules 8-9: apply_evidence


def test_example_1_called_and_rbs_supported() -> None:
    r = ex1_table()
    assert r["rbs_supported"] == True  # noqa: E712
    assert r["called"] == True  # noqa: E712
    assert r["rank"] == 1
    assert r["class_fdr"] == "unestimated"
    assert flags_of(r) == ["class_fdr_unestimated", "rbs_supported"]


def test_example_2_tir_300_not_called() -> None:
    r = ex1_table(tir={ID1: 300.0})
    assert r["rbs_supported"] == False  # noqa: E712
    assert r["called"] == False and is_na(r["rank"])  # noqa: E712
    assert "rbs_supported" not in flags_of(r)


@pytest.mark.parametrize(("min_codons", "called"), [(20, True), (21, False)])
def test_example_3_min_codons_through_the_table(min_codons: int, called: bool) -> None:
    """Example 1's chosen GTG start gives 20 codons: ``min_codons`` 21 makes it a short ORF."""
    r = ex1_table(config=cfg(min_codons=min_codons))
    assert r["called"] == called


def test_example_4_one_peptide_not_called() -> None:
    r = ex1_table(case=ex1_case(n_hits=1), tir={ID1: 5000.0})
    assert r["n_peptides"] == 1
    assert r["rbs_supported"] == True  # noqa: E712
    assert r["called"] == False  # noqa: E712


def test_rule_8_minimum_from_config() -> None:
    r = ex1_table(case=ex1_case(n_hits=1), tir={ID1: 5000.0},
                  config=cfg(min_peptides_without_decoys=1))  # fmt: skip
    assert r["called"] == True  # noqa: E712
    r = ex1_table(config=cfg(min_peptides_without_decoys=3))
    assert r["called"] == False  # noqa: E712


def test_example_5_cds_outframe_called_without_ostir() -> None:
    r = ex1_table(case=ex1_case("cds_outframe"), best={}, tir={}, gt=GT15)
    assert is_na(r["ostir_tir"]) and is_na(r["rbs_supported"])
    assert r["called"] == True  # noqa: E712
    assert "rbs_threshold_unestimated" not in flags_of(r)


def test_cds_outframe_one_peptide_not_called_even_at_min_one() -> None:
    r = ex1_table(case=ex1_case("cds_outframe", n_hits=1),
                  config=cfg(min_peptides_without_decoys=1))  # fmt: skip
    assert r["called"] == False  # noqa: E712


def test_example_6_antisense_one_peptide_not_called() -> None:
    r = ex1_table(case=ex1_case("antisense", n_hits=1))
    assert r["called"] == False  # noqa: E712
    r = ex1_table(case=ex1_case("antisense"))
    assert r["called"] == True  # noqa: E712


def test_example_7_nterm_extension_called() -> None:
    case = nterm2()
    (orf,) = orfs_of(case)
    assert orf.orf_id == ID6 and orf.cls == "nterm_extension" and orf.start_longest == 249
    r = ex1_table(case=case, best={ID6: 249}, tir={ID6: 800.0})
    assert r["start_codon_best"] == "TTG"
    assert r["rbs_supported"] == True  # noqa: E712
    assert r["called"] == True and r["rank"] == 1  # noqa: E712


def test_nterm_extension_whose_longest_start_is_the_annotated_start_not_called() -> None:
    """Degenerate row: ``start_longest`` = the annotated start (hand-built ORF: the M7-06 walk
    never yields it)."""
    case = nterm2()
    (orf,) = orfs_of(case)
    same = replace_orf(orf, starts=(300,), start_longest=300)
    r = ex1_table(case=case, orfs=[same], best={ID6: 300}, tir={ID6: 9000.0})
    assert r["rbs_supported"] == True  # noqa: E712
    assert r["called"] == False  # noqa: E712


def test_nterm_extension_below_threshold_not_called() -> None:
    r = ex1_table(case=nterm2(), best={ID6: 249}, tir={ID6: 300.0})
    assert r["called"] == False  # noqa: E712


@pytest.mark.parametrize("gt", [GT15, gene_table(TIRS_30[:19])])
def test_example_10_few_detected_genes(gt: pd.DataFrame) -> None:
    r = ex1_table(gt=gt)
    assert is_na(r["rbs_supported"])
    assert "rbs_threshold_unestimated" in flags_of(r)
    assert "rbs_supported" not in flags_of(r)
    assert r["called"] == False and is_na(r["rank"])  # noqa: E712


def test_example_10_nterm_extension_with_few_genes() -> None:
    r = ex1_table(case=nterm2(), best={ID6: 249}, tir={ID6: 800.0}, gt=GT15)
    assert "rbs_threshold_unestimated" in flags_of(r) and r["called"] == False  # noqa: E712


def test_quantile_zero_every_scored_start_passes() -> None:
    r = ex1_table(tir={ID1: 100.0}, config=cfg(rbs_min_quantile=0.0))
    assert r["rbs_supported"] == True and r["called"] == True  # noqa: E712
    r = ex1_table(tir={ID1: 99.0}, config=cfg(rbs_min_quantile=0.0))
    assert r["called"] == False  # noqa: E712


def test_ostir_missing_intergenic_not_called() -> None:
    r = ex1_table(best={}, tir={})
    assert is_na(r["rbs_supported"])
    assert "ostir_missing" in flags_of(r)
    assert r["called"] == False  # noqa: E712


@pytest.mark.parametrize("cls", ["pseudogene", "sorf"])
def test_pseudogene_and_sorf_called_on_detection_only(cls: str) -> None:
    case = pseudo_case(None) if cls == "pseudogene" else ex1_case("sorf", gene="S")
    r = evidence_table(case, gt=GT15).to_dict("records")[0]
    assert r["called"] == True  # noqa: E712


def test_cterm_extension_called_when_a_peptide_spans_the_stop() -> None:
    """M7-16 rule 12: one of the two distinct peptides spans the stop [72, 75) and extends past it,
    so the ORF is called (M7-11 decision 2: unestimated, >= 2 distinct peptides); with no
    spanning peptide it is not."""
    case = ex1_case("cterm_extension", gene="G2")  # gene G2 + [42, 75): stop codon [72, 75)
    (orf,) = orfs_of(case)
    spans = {orf.orf_id: [Interval("c1", 72, 102, "+"), Interval("c1", 90, 111, "+")]}
    r = ex1_table(case=case, peptide_intervals=spans, tir={orf.orf_id: 9000.0},
                  best={orf.orf_id: 60})  # fmt: skip
    assert "readthrough_candidate" in flags_of(r)
    assert r["called"] == True and not is_na(r["rank"])  # noqa: E712
    inside = {orf.orf_id: [Interval("c1", 75, 102, "+"), Interval("c1", 90, 111, "+")]}
    r = ex1_table(case=case, peptide_intervals=inside, tir={orf.orf_id: 9000.0},
                  best={orf.orf_id: 60})  # fmt: skip
    assert "readthrough_candidate" not in flags_of(r)
    assert r["called"] == False  # noqa: E712


# ---------------------------------------------------------------- rule 10: ranking


def _rk(**kw: Any) -> dict[str, Any]:
    base = {"orf_id": "NOVEL_c1_1_p", "n_peptides": 2, "best_peptide_q": 0.001,
            "start_codon_best": "ATG", "ostir_tir": 500.0, "promoter_score": 0.5,
            "terminator_conf": 0.5, "length_aa_best": 50, "called": True}  # fmt: skip
    return row(**{**base, **kw})


def _order(rows: list[dict[str, Any]]) -> list[str]:
    return [r["orf_id"] for r in sorted(rows, key=ev().rank_key)]


def test_example_11_peptides_then_q_then_start_codon() -> None:
    a = _rk(orf_id="A", n_peptides=3, best_peptide_q=0.001, start_codon_best="ATG")
    b = _rk(orf_id="B", n_peptides=3, best_peptide_q=0.001, start_codon_best="GTG")
    c = _rk(orf_id="C", n_peptides=2, best_peptide_q=1e-5, start_codon_best="ATG")
    assert _order([c, b, a]) == ["A", "B", "C"]
    assert _order([b, a, c]) == ["A", "B", "C"]


def test_example_12_full_tie_lower_orf_id_first() -> None:
    rows = [_rk(orf_id=i) for i in ("NOVEL_c1_90_p", "NOVEL_c1_100_p", "NOVEL_c1_5_m")]
    assert _order(rows) == sorted(r["orf_id"] for r in rows)


@pytest.mark.parametrize(
    ("key", "better", "worse"),
    [
        ("n_peptides", 4, 3),
        ("best_peptide_q", 1e-4, 1e-3),
        ("best_peptide_q", 0.5, None),
        ("start_codon_best", "ATG", "GTG"),
        ("start_codon_best", "GTG", "TTG"),
        ("start_codon_best", "TTG", "CTG"),
        ("start_codon_best", "CTG", None),
        ("ostir_tir", 900.0, 800.0),
        ("ostir_tir", 1.0, None),
        ("promoter_score", 0.9, 0.8),
        ("promoter_score", 0.0, None),
        ("terminator_conf", 0.9, 0.8),
        ("terminator_conf", 0.1, None),
        ("length_aa_best", 60, 59),
    ],
)
def test_each_rank_key_in_isolation(key: str, better: Any, worse: Any) -> None:
    """Rule 10 direction of each key; NA sorts last in every key. ``orf_id`` is set against the
    key (the worse row has the lower ID) so a skipped key would show."""
    w = _rk(orf_id="NOVEL_c1_1_p", **{key: worse})
    b = _rk(orf_id="NOVEL_c1_2_p", **{key: better})
    assert _order([w, b]) == ["NOVEL_c1_2_p", "NOVEL_c1_1_p"]


def test_nan_counts_as_na_in_rank_key() -> None:
    w = _rk(orf_id="NOVEL_c1_1_p", ostir_tir=float("nan"))
    b = _rk(orf_id="NOVEL_c1_2_p", ostir_tir=0.0)
    assert _order([w, b]) == ["NOVEL_c1_2_p", "NOVEL_c1_1_p"]
    w = _rk(orf_id="NOVEL_c1_1_p", promoter_score=pd.NA)
    b = _rk(orf_id="NOVEL_c1_2_p", promoter_score=0.0)
    assert _order([w, b]) == ["NOVEL_c1_2_p", "NOVEL_c1_1_p"]


def test_key_priority_terminator_after_promoter_before_length() -> None:
    a = _rk(orf_id="A", promoter_score=0.9, terminator_conf=0.1, length_aa_best=10)
    b = _rk(orf_id="B", promoter_score=0.8, terminator_conf=0.9, length_aa_best=90)
    c = _rk(orf_id="C", promoter_score=0.8, terminator_conf=0.8, length_aa_best=99)
    assert _order([c, b, a]) == ["A", "B", "C"]
    d = _rk(orf_id="D", ostir_tir=501.0, promoter_score=0.0)
    assert _order([a, d]) == ["D", "A"]


def test_apply_ranks_called_orfs_only() -> None:
    """Decision 5: ranks 1..k over called rows; uncalled rows NA. ``two_contigs`` has two called
    example-1 ORFs (zz and c1) and a one-peptide ORF on c1."""
    case = two_contigs()
    orfs = orfs_of(case)
    ids = [o.orf_id for o in orfs]
    one = next(o.orf_id for o in orfs if len(o.peptides) == 1)
    tir = dict.fromkeys(ids, 500.0)
    t = evidence_table(case, tir=tir, best={o.orf_id: o.starts[-1] for o in orfs})
    by = {r["orf_id"]: r for r in t.to_dict("records")}
    assert by[one]["called"] == False and is_na(by[one]["rank"])  # noqa: E712
    called = [r for r in by.values() if r["called"]]
    assert len(called) == 2
    assert sorted(int(r["rank"]) for r in called) == [1, 2]
    want = [r["orf_id"] for r in sorted(called, key=ev().rank_key)]
    assert [r["orf_id"] for r in sorted(called, key=lambda r: int(r["rank"]))] == want


def test_apply_ranking_follows_promoter_then_terminator() -> None:
    """The two example-1 copies of ``two_contigs`` tie on peptides, q, codon and TIR."""
    case = two_contigs()
    orfs = [o for o in orfs_of(case) if len(o.peptides) == 2]
    zz, c1 = sorted(orfs, key=lambda o: o.contig != "zz")
    kw: dict[str, Any] = {"orfs": orfs, "tir": {zz.orf_id: 500.0, c1.orf_id: 500.0},
                          "best": {zz.orf_id: 60, c1.orf_id: 60}}  # fmt: skip
    t = evidence_table(case, promoter={zz.orf_id: 0.2, c1.orf_id: 0.7}, **kw)
    rank = dict(zip(t["orf_id"], t["rank"], strict=True))
    assert (rank[c1.orf_id], rank[zz.orf_id]) == (1, 2)
    terms = [term("zz", 199, 230, "+", 0.9), term("c1", 199, 230, "+", 0.6)]
    t = evidence_table(case, terminators=terms, **kw)
    rank = dict(zip(t["orf_id"], t["rank"], strict=True))
    assert (rank[zz.orf_id], rank[c1.orf_id]) == (1, 2)


# ---------------------------------------------------------------- flags and schema


def test_apply_keeps_columns_rows_and_order() -> None:
    case = two_contigs()
    orfs = orfs_of(case)
    t = evidence_table(case, tir={o.orf_id: 500.0 for o in orfs},
                       best={o.orf_id: o.starts[-1] for o in orfs})  # fmt: skip
    assert tuple(t.columns) == NOVEL_ORFS_COLUMNS
    base = table(case, orfs, scored(orfs))
    assert list(t["orf_id"]) == list(base["orf_id"])
    assert str(t["rank"].dtype) == "Int64"
    for col in ("rbs_supported", "promoter_supported", "terminator_downstream"):
        assert str(t[col].dtype) == "boolean", col
    assert str(t["called"].dtype) in ("bool", "boolean")


def test_flags_in_section_9_4_order() -> None:
    """Every flag combination this module writes is ``;``-joined in FLAG_ORDER."""
    case = pseudo_case((595, "AAAAA"))
    (orf,) = orfs_of(case)
    t = evidence_table(case, best={orf.orf_id: orf.starts[-1]}, tir={orf.orf_id: 9000.0},
                       promoter={orf.orf_id: 0.9},
                       terminators=[term("c1", 650, 680, "+", 0.9)],
                       peptide_intervals={orf.orf_id: [Interval("c1", 590, 620, "+")]})  # fmt: skip
    f = flags_of(t.iloc[0])
    assert f == ["junction_spanned", "homopolymer_junction", "class_fdr_unestimated",
                 "rbs_supported", "promoter_supported", "terminator_supported"]  # fmt: skip


def test_flag_order_matches_the_schema_pattern() -> None:
    pattern = next(f for f in json.loads(SCHEMA.read_text(encoding="utf-8"))["fields"]
                   if f["name"] == "flags")["constraints"]["pattern"]  # fmt: skip
    alt = pattern.split("(?:", 2)[1].split(")", 1)[0].split("|")
    assert tuple(alt) == FLAG_ORDER


def test_flags_mixed_with_m7_07_flags_keep_order() -> None:
    r = ex1_table(best={}, tir={}, gt=GT15, promoter={ID1: 0.9})
    f = flags_of(r)
    assert f == sorted(f, key=FLAG_ORDER.index)
    assert set(f) == {"class_fdr_unestimated", "ostir_missing", "rbs_threshold_unestimated",
                      "promoter_supported"}  # fmt: skip


def test_tsv_round_trip_validates(tmp_path: Path) -> None:
    from helpers import novel_steps_fixture as nf

    from giap.io.novel_tables import read_novel_orfs, write_novel_orfs

    case = two_contigs()
    orfs = orfs_of(case)
    t = evidence_table(case, tir={o.orf_id: 500.0 for o in orfs}, promoter={orfs[0].orf_id: 0.8},
                       best={o.orf_id: o.starts[-1] for o in orfs},
                       terminators=[term("zz", 199, 230, "+", 0.8)])  # fmt: skip
    p = tmp_path / "novel_orfs.tsv"
    write_novel_orfs(p, t)
    back = read_novel_orfs(p)
    assert list(back["rank"].astype("object").where(back["rank"].notna(), None)) == list(
        t["rank"].astype("object").where(t["rank"].notna(), None)
    )
    nf.assert_valid_novel_orfs(p)


def test_gene_table_unused_columns_not_required() -> None:
    """``apply_evidence`` reads only ``partial``, ``detected`` and ``tir_ostir``."""
    gt = GT30[["partial", "detected", "tir_ostir"]]
    assert ex1_table(gt=gt)["called"] == True  # noqa: E712
