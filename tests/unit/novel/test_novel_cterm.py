"""``cterm_extension`` calls and patch rows (M7-16b rule 12; spec examples 6 and 12; ADR-0013
items 2-3; Q3 default: a window ending on the stop codon is no ``cterm_extension`` evidence).

* M7-09 ``class_rules``: ``cterm_extension`` is called when ``spans_annotated_stop`` holds and the
  in-frame continuation exists (an M7-06 ORF: a next stop or ``end_3p_open``); label
  ``readthrough_candidate``; then the M7-11 decision-2 rule (``unestimated``: >= 2 distinct
  peptides). M7-09 decision 4's ``ok = False`` is gone.
* ``spans_annotated_stop``: some peptide interval strictly contains the gene's terminal stop codon
  and extends past it on the 3' side (lower coordinates on ``-``).
* M7-07 ``patch_features``: a called ``cterm_extension`` row is a ``CDS`` with ``ID`` =
  ``locus_tag`` = the gene, ``Note=cterm_extension``, from the gene start to the ORF's stop-codon
  end; ``partial=true`` when ``end_3p_open``.

World: toy contig ``c1`` (linear, 126 nt) of the spec: ``G1`` = CDS ``+`` ``[0, 57)`` with its TGA
at ``[54, 57)``; the next in-frame stop is the TAA at ``[90, 93)``. Peptides: example 1's
``WQLNEFKWVPAGHTL`` ``[33, 78)`` and ``KWVPAGHTLDMR`` ``[51, 87)`` (both cost 1, over the stop).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from helpers.novel_evidence_fixture import cfg as ev_cfg
from helpers.novel_evidence_fixture import evidence_table, flags_of
from helpers.novel_oracle import Case, cds, hit_row, mirror_case, mirror_interval
from helpers.novel_outputs_fixture import orfs_of, scored, table
from helpers.readthrough_toys import toy_c1

from giap.io.gff3 import write_gff3
from giap.model import Replicon
from giap.model.coords import Interval

N = 126
P1, P2 = ("WQLNEFKWVPAGHTL", 33, 78), ("KWVPAGHTLDMR", 51, 87)


def _hit(pep: str, s: int, e: int, n: int = N) -> dict[str, Any]:
    r = hit_row(pep, "c1", s, e, "+", "cterm_extension", n, gene="G1", route="tolerant")
    r["n_diff"], r["diff_types"] = 1, ["substitution"]
    return r


def case(
    peps: tuple[tuple[str, int, int], ...] = (P1, P2),
    genome: dict[str, str] | None = None,
    n: int = N,
) -> Case:
    g = genome or toy_c1()[0]
    reps = {"c1": Replicon("c1", n, "linear")}
    return Case(
        dict(g), reps, [_hit(p, s, e, n) for p, s, e in peps], genes=[cds("G1", "c1", 0, 57, "+")]
    )


def intervals(c: Case) -> dict[str, list[Interval]]:
    (orf,) = orfs_of(c)
    return {orf.orf_id: [Interval(h["contig"], h["start"], h["end"], h["strand"]) for h in c.hits]}


def row_of(c: Case, peptide_intervals: dict[str, list[Interval]] | None = None) -> dict[str, Any]:
    (r,) = evidence_table(c, peptide_intervals=peptide_intervals or intervals(c)).to_dict("records")
    return r


def patch_lines(c: Case, tmp: Path, peptide_intervals: Any = None) -> list[str]:
    from helpers.novel_evidence_fixture import GT30

    from giap.novel.evidence import apply_evidence
    from giap.novel.outputs import patch_features

    orfs = orfs_of(c)
    t = table(c, orfs, scored(orfs))
    t = apply_evidence(
        t,
        orfs,
        gene_table=GT30,
        terminators=None,
        gm=c.gm,
        genome=c.genome,
        cfg=ev_cfg(),
        peptide_intervals=peptide_intervals or intervals(c),
    )
    path = tmp / "annotation_patch.gff3"
    write_gff3(path, patch_features(t, c.gm), list(c.replicons.values()))
    return [ln for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]


# ---------------------------------------------------------------- example 6


def test_example_6_orf_is_g1_start_to_the_next_stop() -> None:
    (orf,) = orfs_of(case())
    assert (orf.cls, orf.overlapping_gene, orf.strand, orf.start_longest) == (
        "cterm_extension",
        "G1",
        "+",
        0,
    )
    assert orf.stop == 90 and not orf.end_3p_open


def test_example_6_called_with_readthrough_candidate() -> None:
    r = row_of(case())
    assert r["class"] == "cterm_extension"
    assert "readthrough_candidate" in flags_of(r)
    assert r["class_fdr"] == "unestimated"
    assert r["called"] == True  # noqa: E712
    assert r["rank"] == 1


def test_example_6_class_rule() -> None:
    from helpers.novel_evidence_fixture import GT30

    from giap.novel.evidence import class_rules

    c = case()
    orfs = orfs_of(c)
    t = table(c, orfs, scored(orfs))
    got = class_rules(
        t,
        orfs,
        gene_table=GT30,
        gm=c.gm,
        genome=c.genome,
        cfg=ev_cfg(),
        peptide_intervals=intervals(c),
    )
    assert got == {orfs[0].orf_id: True}


def test_one_spanning_peptide_is_a_candidate_but_not_called() -> None:
    """``unestimated`` needs >= 2 distinct peptides (M7-11 decision 2)."""
    r = row_of(case((P1,)))
    assert "readthrough_candidate" in flags_of(r)
    assert r["called"] == False  # noqa: E712


def test_example_6_patch_row(tmp_path: Path) -> None:
    lines = patch_lines(case(), tmp_path)
    assert lines == ["c1\tGIAP\tCDS\t1\t93\t.\t+\t0\tID=G1;locus_tag=G1;Note=cterm_extension"]


# ---------------------------------------------------------------- spans_annotated_stop


@pytest.mark.parametrize(
    ("ivs", "spans"),
    [
        ([(33, 78)], True),  # across the stop into the 3' flank
        ([(51, 87)], True),
        ([(54, 66)], True),  # starts on the codon and extends past it
        ([(27, 57)], False),  # example 12: ends on the stop (Q3), nothing past it
        ([(54, 57)], False),  # the codon itself
        ([(57, 87)], False),  # only the flank
        ([(33, 54)], False),  # only the gene
        ([(27, 57), (57, 87)], False),  # neither interval alone spans the stop
    ],
)
def test_spans_annotated_stop_needs_containment_and_the_3p_side(
    ivs: list[tuple[int, int]], spans: bool
) -> None:
    c = case()
    (orf,) = orfs_of(c)
    r = row_of(c, {orf.orf_id: [Interval("c1", s, e, "+") for s, e in ivs]})
    assert ("readthrough_candidate" in flags_of(r)) is spans
    assert bool(r["called"]) is spans


@pytest.mark.parametrize(("iv", "spans"), [((33, 78), True), ((27, 57), False), ((57, 87), False)])
def test_minus_strand_mirror(iv: tuple[int, int], spans: bool) -> None:
    """Mirror of the world: G1 ``-`` ``[69, 126)``, stop ``[69, 72)``; the 3' side is below 69."""
    m = mirror_case(case())
    (orf,) = orfs_of(m)
    a, b = mirror_interval(*iv, N)
    r = row_of(m, {orf.orf_id: [Interval("c1", a, b, "-")]})
    assert ("readthrough_candidate" in flags_of(r)) is spans


def test_minus_strand_patch_row(tmp_path: Path) -> None:
    m = mirror_case(case())
    (orf,) = orfs_of(m)
    ivs = {orf.orf_id: [Interval("c1", *mirror_interval(s, e, N), "-") for _p, s, e in (P1, P2)]}
    lines = patch_lines(m, tmp_path, ivs)
    assert lines == ["c1\tGIAP\tCDS\t34\t126\t.\t-\t0\tID=G1;locus_tag=G1;Note=cterm_extension"]


def test_open_3p_end_patch_row_is_partial(tmp_path: Path) -> None:
    """The contig ends (linear, 90 nt) before any further in-frame stop: ``end_3p_open``."""
    g = {"c1": toy_c1()[0]["c1"][:90]}
    c = case(genome=g, n=90)
    (orf,) = orfs_of(c)
    assert orf.end_3p_open
    r = row_of(c)
    assert r["called"] == True and "readthrough_candidate" in flags_of(r)  # noqa: E712
    (line,) = patch_lines(c, tmp_path)
    f = line.split("\t")
    assert (f[0], f[2], f[3], f[4], f[6]) == ("c1", "CDS", "1", "90", "+")
    kv = dict(x.split("=", 1) for x in f[8].split(";"))
    assert (kv["ID"], kv["locus_tag"], kv["Note"], kv.get("partial")) == (
        "G1",
        "G1",
        "cterm_extension",
        "true",
    )


def test_uncalled_cterm_row_writes_no_patch_row(tmp_path: Path) -> None:
    assert patch_lines(case((P1,)), tmp_path) == []
