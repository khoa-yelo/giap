"""Review round 1 of PR #180 (M7-09; docs/specs/M7-09-novel-evidence-and-ranking.md): the
decision-5 peptide-interval path, rule 7's "best confidence" and the minus-strand / origin
geometry of ``apply_evidence``.

* B1 (decision 5): ``novel_candidates.parquet`` carries each peptide's mapped interval as
  ``peptide_start``/``peptide_end`` (0-based half-open, exactly the C11 hit's ``[start, end)``),
  found by the I/L-folded peptide; ``read_peptide_intervals`` returns them per ORF (``{}`` for a
  pre-M7-09 file, rows without an interval skipped) and they drive ``junction_spanned``.
* B2 (rule 7): of two same-strand calls in the 0..300 window the higher confidence wins, even
  when it is farther.
* B3: the minus pseudogene scan excludes the gene's last codon; the minus annotated stop codon
  is ``[g.start, g.start + 3)``; ``start_upstream_nt`` folds across the origin; codon
  containment and homopolymer runs wrap on circular replicons and not on linear ones.

The step path is exercised through the functions the ``novel_candidates`` / ``novel_orfs``
steps call (``candidate_rows`` -> ``reconstruct`` -> ``_candidate_records`` ->
``write_parquet_deterministic(_CAND_COLS)`` -> ``read_peptide_intervals`` -> ``apply_evidence``).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from helpers.novel_evidence_fixture import (
    GT30,
    ID1,
    circ,
    evidence_table,
    ex1_case,
    ex1_table,
    flags_of,
    lin,
    nterm2,
    pseudo_case,
    term,
)
from helpers.novel_oracle import Case, Seq, cds, hit_row, il, mirror_case, revcomp, translate
from helpers.novel_outputs_fixture import C1, L, ex1, orfs_of

from giap.model.coords import Interval
from giap.novel.candidates import Orf, candidate_rows, peptide_loci, reconstruct
from giap.novel.evidence import (
    junction_labels,
    pseudogene_junction,
    terminator_downstream,
)

# ---------------------------------------------------------------- B1: decision-5 data path


def _candidates_parquet(case: Case, path: Path) -> list[Orf]:
    """What the ``novel_candidates`` step writes for ``case``."""
    from helpers.novel_outputs_fixture import cfg as out_cfg

    from giap.io.parquet import write_parquet_deterministic
    from giap.steps.novel import _CAND_COLS, _candidate_records

    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    orfs = reconstruct(rows, case.gm, case.genome, out_cfg())
    recs = _candidate_records(orfs, case.replicons, peptide_loci(rows))
    write_parquet_deterministic(path, _CAND_COLS, recs)
    return orfs


def test_peptide_loci_is_keyed_on_the_il_folded_peptide() -> None:
    case = ex1()
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    with_i = [h for h in case.hits if "I" in h["peptide"]]
    assert with_i, "fixture must contain an I peptide"
    loci = peptide_loci(rows)
    for h in case.hits:
        assert loci[("c1", il(h["peptide"]))] == Interval("c1", h["start"], h["end"], h["strand"])


def test_novel_candidates_peptide_columns_are_the_exact_hit_intervals(tmp_path: Path) -> None:
    """Example 1 (an I-containing peptide): ``peptide_start``/``peptide_end`` = the hit's
    0-based half-open ``[start, end)``, keyed by the folded peptide the table stores."""
    case = ex1()
    path = tmp_path / "novel_candidates.parquet"
    _candidates_parquet(case, path)
    df = pd.read_parquet(path)
    want = {il(h["peptide"]): (h["start"], h["end"]) for h in case.hits}
    got = {
        str(p): (int(a), int(b))
        for p, a, b in zip(df["peptide"], df["peptide_start"], df["peptide_end"], strict=True)
    }
    assert got == want
    assert any("I" in h["peptide"] for h in case.hits)


def test_read_peptide_intervals_round_trip(tmp_path: Path) -> None:
    from giap.steps.novel import read_peptide_intervals

    case = ex1()
    path = tmp_path / "novel_candidates.parquet"
    _candidates_parquet(case, path)
    got = read_peptide_intervals(path)
    want = sorted((h["start"], h["end"]) for h in case.hits)
    assert set(got) == {ID1}
    assert sorted((iv.start, iv.end) for iv in got[ID1]) == want
    assert all(iv.contig == "c1" and iv.strand == "+" for iv in got[ID1])


def test_read_peptide_intervals_old_file_without_columns_is_empty(tmp_path: Path) -> None:
    from giap.steps.novel import read_peptide_intervals

    path = tmp_path / "old.parquet"
    _candidates_parquet(ex1(), path)
    df = pd.read_parquet(path).drop(columns=["peptide_start", "peptide_end"])
    df.to_parquet(path)
    assert read_peptide_intervals(path) == {}


def test_read_peptide_intervals_skips_a_peptide_without_a_locus(tmp_path: Path) -> None:
    from giap.io.parquet import write_parquet_deterministic
    from giap.steps.novel import _CAND_COLS, _candidate_records, read_peptide_intervals

    case = ex1()
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    from helpers.novel_outputs_fixture import cfg as out_cfg

    orfs = reconstruct(rows, case.gm, case.genome, out_cfg())
    recs = _candidate_records(orfs, case.replicons, peptide_loci(rows))
    assert len(recs) == 2
    first, second = recs
    recs = [(*first[:-2], None, None), second]
    path = tmp_path / "novel_candidates.parquet"
    write_parquet_deterministic(path, _CAND_COLS, recs)
    got = read_peptide_intervals(path)
    assert got == {ID1: [Interval("c1", second[-2], second[-1], "+")]}


# A pseudogene P + [540, 900) (ATG [540, 543), TAA [897, 900)) whose first premature in-frame
# stop is TGA [600, 603); the peptide is mapped in frame +2 (start 592 = 540 + 52) across it,
# so its interval [592, 622) contains the junction codon: AAAELLLLLI (an I).
SHIFTED_NT = "GCTGCTGCTGAG" + "CTG" * 5 + "ATC"


def _shifted_pseudo_case() -> Case:
    seq = Seq(L).fill(543, 897).put(540, "ATG").put(897, "TAA").put(592, SHIFTED_NT).text()
    assert seq[600:603] == "TGA"
    pep = translate(seq[592:622])
    assert pep == "AAAELLLLLI"
    hits = [hit_row(pep, "c1", 592, 622, "+", "pseudogene", L, gene="P")]
    return Case({"c1": seq}, dict(C1), hits, genes=[cds("P", "c1", 540, 900, "+", "pseudogene")])


def test_junction_spanned_from_the_written_peptide_intervals(tmp_path: Path) -> None:
    """Decision 5 end to end: the step's ``novel_candidates`` intervals, read back by
    ``read_peptide_intervals``, label the pseudogene ``junction_spanned``."""
    from giap.steps.novel import read_peptide_intervals

    case = _shifted_pseudo_case()
    assert pseudogene_junction(case.genes[0], case.genome["c1"], circ()) == 600
    path = tmp_path / "novel_candidates.parquet"
    orfs = _candidates_parquet(case, path)
    (orf,) = orfs
    assert orf.cls == "pseudogene"
    ivs = read_peptide_intervals(path)
    assert ivs == {orf.orf_id: [Interval("c1", 592, 622, "+")]}
    (r,) = evidence_table(case, orfs=orfs, peptide_intervals=ivs).to_dict("records")
    f = flags_of(r)
    assert "junction_spanned" in f and "junction_unknown" not in f


def test_readthrough_candidate_from_the_written_peptide_intervals(tmp_path: Path) -> None:
    """Decision 5 for ``cterm_extension``: the written intervals reach the readthrough label
    (gene G2 + [42, 75), its stop codon [72, 75) inside the second peptide)."""
    from giap.steps.novel import read_peptide_intervals

    case = ex1_case("cterm_extension", gene="G2")
    path = tmp_path / "novel_candidates.parquet"
    orfs = _candidates_parquet(case, path)
    ivs = read_peptide_intervals(path)
    (orf,) = orfs
    assert sorted((iv.start, iv.end) for iv in ivs[orf.orf_id]) == [(72, 102), (90, 111)]
    (r,) = evidence_table(case, orfs=orfs, peptide_intervals=ivs).to_dict("records")
    assert "readthrough_candidate" in flags_of(r)


# ---------------------------------------------------------------- B2: rule 7 best confidence


def _orf(strand: str, stop: int) -> Orf:
    return Orf(f"NOVEL_c1_1_{'p' if strand == '+' else 'm'}", "intergenic", "c1",
               strand, 1 if strand == "+" else -1, stop, (), None, False, False,  # type: ignore[arg-type]
               ("AAAAAAA", "CCCCCCC"), 0.001, None)  # fmt: skip


@pytest.mark.parametrize("order", [1, -1])
def test_terminator_best_confidence_beats_the_nearer_plus(order: int) -> None:
    """Stop [120, 123) (last base 122): d = 30 at conf 0.4, d = 100 at conf 0.9 -> 0.9, 100."""
    terms = [term("c1", 152, 180, "+", 0.4), term("c1", 222, 250, "+", 0.9)][::order]
    assert terminator_downstream(_orf("+", 120), terms, circ()) == (pytest.approx(0.9), 100)


@pytest.mark.parametrize("order", [1, -1])
def test_terminator_best_confidence_beats_the_nearer_minus(order: int) -> None:
    """Minus stop [877, 880) (last base 877): 5' ends 847 (d = 30, 0.4) and 777 (d = 100, 0.9)."""
    terms = [term("c1", 820, 848, "-", 0.4), term("c1", 750, 778, "-", 0.9)][::order]
    assert terminator_downstream(_orf("-", 877), terms, circ()) == (pytest.approx(0.9), 100)


def test_apply_best_confidence_terminator_is_supported() -> None:
    """``terminatornet_yes`` 0.5: the 0.9 call at 100 nt makes the row supported."""
    r = ex1_table(terminators=[term("c1", 152, 180, "+", 0.4), term("c1", 222, 250, "+", 0.9)])
    assert r["terminator_conf"] == pytest.approx(0.9)
    assert r["terminator_distance_nt"] == 100
    assert r["terminator_downstream"] == True  # noqa: E712
    assert "terminator_supported" in flags_of(r)


# ---------------------------------------------------------------- B3: minus strand


def test_minus_pseudogene_without_premature_stop_has_no_junction() -> None:
    """Gene - [100, 460): its own terminal stop [100, 103) is not a junction."""
    plus = Seq(L).fill(543, 897).put(540, "ATG").put(897, "TAA").text()
    gene = cds("P", "c1", 100, 460, "-", "pseudogene")
    assert pseudogene_junction(gene, revcomp(plus), circ()) is None


def test_minus_pseudogene_without_premature_stop_flags_junction_unknown() -> None:
    case = mirror_case(pseudo_case(None, premature=False))
    (g,) = case.genes
    assert g.interval.strand == "-" and (g.interval.start, g.interval.end) == (100, 460)
    (orf,) = orfs_of(case)
    t = evidence_table(case, peptide_intervals={orf.orf_id: [Interval("c1", 90, 130, "-")]})
    f = flags_of(t.iloc[0])
    assert "junction_unknown" in f
    assert "junction_spanned" not in f


def test_minus_pseudogene_with_premature_stop_mirror_of_example_8() -> None:
    case = mirror_case(pseudo_case(None))
    (orf,) = orfs_of(case)
    junction = 1000 - 603
    assert pseudogene_junction(case.genes[0], case.genome["c1"], circ()) == junction
    spans = {orf.orf_id: [Interval("c1", 380, 410, "-")]}
    f = flags_of(evidence_table(case, peptide_intervals=spans).iloc[0])
    assert "junction_spanned" in f and "junction_unknown" not in f


def _minus_cterm() -> tuple[Case, Orf]:
    """Mirror of ``ex1_case("cterm_extension", gene="G2")``: gene G2 - [925, 958), its annotated
    stop codon [925, 928)."""
    case = mirror_case(ex1_case("cterm_extension", gene="G2"))
    (g,) = case.genes
    assert g.interval.strand == "-" and (g.interval.start, g.interval.end) == (925, 958)
    (orf,) = orfs_of(case)
    assert orf.cls == "cterm_extension" and orf.strand == "-"
    return case, orf


@pytest.mark.parametrize(
    ("iv", "spans"),
    [((898, 928), True), ((925, 928), False), ((898, 925), False), ((926, 940), False),
     ((928, 958), False)],
)  # fmt: skip
def test_minus_cterm_extension_readthrough_needs_the_annotated_stop(
    iv: tuple[int, int], spans: bool
) -> None:
    case, orf = _minus_cterm()
    peps = {orf.orf_id: [Interval("c1", *iv, "-")]}
    (r,) = evidence_table(case, peptide_intervals=peps).to_dict("records")
    # M7-16 rule 12: the stop codon [925, 928) must be contained AND passed on the 3' side (below
    # 925 on -); the codon alone no longer spans it, and a spanning ORF is now called
    assert ("readthrough_candidate" in flags_of(r)) is spans
    assert bool(r["called"]) is spans


# ---------------------------------------------------------------- B3: origin geometry


def _rotated_nterm2(shift: int) -> Case:
    """``nterm2`` moved by ``+shift`` on the circular contig (gene included)."""
    base = nterm2()
    seq = base.genome["c1"]
    k = shift % L
    genome = {"c1": seq[L - k :] + seq[: L - k]}
    hits = []
    for h in base.hits:
        a = (h["start"] + shift) % L
        hits.append(hit_row(h["peptide"], "c1", a, a + h["end"] - h["start"], "+",
                            h["locus_class"], L, gene=h["gene"]))  # fmt: skip
    (g,) = base.genes
    gs = (g.interval.start + shift) % L
    gene = cds(g.locus_tag, "c1", gs, gs + g.interval.end - g.interval.start, "+")
    return Case(genome, dict(C1), hits, genes=[gene])


def test_nterm_extension_upstream_start_across_the_origin_is_called() -> None:
    """Gene G + [10, 610), TTG start_longest at 959 (= 249 - 290 + 1000): 51 nt upstream through
    the origin, so the row is called exactly like example 7."""
    case = _rotated_nterm2(-290)
    (orf,) = orfs_of(case)
    assert orf.cls == "nterm_extension" and orf.start_longest == 959
    (r,) = evidence_table(case, orfs=[orf], best={orf.orf_id: 959},
                          tir={orf.orf_id: 800.0}, gt=GT30).to_dict("records")  # fmt: skip
    assert r["start_codon_best"] == "TTG"
    assert r["rbs_supported"] == True  # noqa: E712
    assert r["called"] == True and r["rank"] == 1  # noqa: E712


def test_nterm_extension_start_at_n_minus_20_gene_at_10_is_upstream_30() -> None:
    """The review's case through ``apply_evidence``: gene + at 10, ``start_longest`` = n - 20."""
    case = _rotated_nterm2(-290)
    (orf,) = orfs_of(case)
    from helpers.novel_evidence_fixture import replace_orf

    near = replace_orf(orf, starts=(L - 20,), start_longest=L - 20)
    seq = list(case.genome["c1"])
    seq[L - 20 : L - 17] = "ATG"
    case = Case({"c1": "".join(seq)}, case.replicons, case.hits, case.genes)
    (r,) = evidence_table(case, orfs=[near], best={orf.orf_id: L - 20},
                          tir={orf.orf_id: 800.0}).to_dict("records")  # fmt: skip
    assert r["called"] == True  # noqa: E712


def _genome_gct(n: int = L) -> list[str]:
    return list(("GCT" * (n // 3 + 1))[:n])


def test_junction_spanned_by_a_peptide_across_the_origin() -> None:
    """Circular: junction codon [2, 5) inside a peptide stored unwrapped as [990, 1010)."""
    g = "".join(_genome_gct())
    assert junction_labels(2, [Interval("c1", 990, 1010, "+")], g, circ()) == {"junction_spanned"}
    assert "junction_spanned" not in junction_labels(
        12, [Interval("c1", 990, 1010, "+")], g, circ()
    )


def test_pseudogene_across_the_origin_junction_spanned_through_apply() -> None:
    """Gene P + [950, 1100) (unwrapped), first premature TGA at [1040, 1043) = 40; a peptide
    [1020, 1050) (unwrapped) spans it."""
    s = Seq(L).fill(953, 1100).put(950, "ATG").put(1097, "TAA").put(1040, "TGA")
    seq = s.text()
    gene = cds("P", "c1", 950, 1100, "+", "pseudogene")
    hits = [hit_row(translate(seq[960:990]), "c1", 960, 990, "+", "pseudogene", L, gene="P")]
    case = Case({"c1": seq}, dict(C1), hits, genes=[gene])
    (orf,) = orfs_of(case)
    peps = {orf.orf_id: [Interval("c1", 1020, 1050, "+")]}
    f = flags_of(evidence_table(case, peptide_intervals=peps).iloc[0])
    assert "junction_spanned" in f


def test_homopolymer_across_the_origin() -> None:
    """AAA [997, 1000) + AA [0, 2) is one run of 5 on a circular contig; junction 3 -> window
    [-2, 11) overlaps it. On a linear contig the run is two pieces (3 and 2): no label."""
    b = _genome_gct()
    b[997:1000] = "AAA"
    b[0:2] = "AA"
    g = "".join(b)
    assert g[996] != "A" and g[2] != "A"
    assert junction_labels(3, [], g, circ()) == {"homopolymer_junction"}
    assert junction_labels(3, [], g, lin()) == set()


def test_homopolymer_across_the_origin_left_of_the_junction_window() -> None:
    """Run AAAAA [998, 1003) wrapped; junction 7 -> window [2, 15) only touches its last base."""
    b = _genome_gct()
    b[998:1000] = "AA"
    b[0:3] = "AAA"
    g = "".join(b)
    assert g[997] != "A" and g[3] != "A"
    assert junction_labels(7, [], g, circ()) == {"homopolymer_junction"}
    assert junction_labels(7, [], g, lin()) == set()


@pytest.mark.parametrize(
    ("iv", "spanned"),
    [((590, 620), True), ((600, 603), True), ((590, 602), False), ((601, 630), False),
     ((590, 601), False)],
)  # fmt: skip
def test_junction_containment_on_a_linear_replicon(iv: tuple[int, int], spanned: bool) -> None:
    g = "".join(_genome_gct())
    got = junction_labels(600, [Interval("c1", *iv, "+")], g, lin())
    assert ("junction_spanned" in got) is spanned


def test_linear_pseudogene_containment_through_apply() -> None:
    case = pseudo_case(None)
    lin_case = Case(case.genome, {"c1": lin()}, case.hits, case.genes)
    (orf,) = orfs_of(lin_case)
    out: dict[tuple[int, int], bool] = {}
    for iv in ((590, 620), (590, 602)):
        t = evidence_table(lin_case, peptide_intervals={orf.orf_id: [Interval("c1", *iv, "+")]})
        out[iv] = "junction_spanned" in flags_of(t.iloc[0])
    assert out == {(590, 620): True, (590, 602): False}
