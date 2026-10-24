"""``novel_candidate_id`` back-fill (M7-06 Q6, M7-15 item 6; review #278 B2, lead decision 11):
``giap.novel.outputs.candidate_ids`` and ``giap.steps.novel.backfill_candidate_ids``.

Kills the review's surviving mutants "locus-class filter removed" and "strand dropped from the
owner key". Pins the decision-11 nits: a hit matches an ORF only within the ORF span (positions;
``start_longest`` > ``end`` on ``-``), and the back-fill writes every row (``NA`` when a row has
no id), so a stale id from an earlier run is cleared.
"""

from __future__ import annotations

import json
from pathlib import Path

from hypothesis import given
from hypothesis import strategies as st

from giap.novel.candidates import NOVEL_CLASSES
from giap.novel.outputs import candidate_ids
from giap.steps.novel import backfill_candidate_ids

Row = dict[str, str | None]


def hit(
    peptide: str = "PEPTLDEK",
    contig: str = "chr1",
    strand: str = "+",
    locus_class: str | None = "intergenic",
    start: int = 40,
    end: int = 63,
) -> Row:
    return {
        "peptide": peptide,
        "contig": contig,
        "strand": strand,
        "locus_class": locus_class,
        "start": str(start),
        "end": str(end),
        "novel_candidate_id": None,
    }


def orf(
    orf_id: str = "NOVEL_chr1_100_p",
    peptides: str | None = "PEPTIDEK;LLAMAR",
    contig: str = "chr1",
    strand: str = "+",
    start_longest: int = 10,
    end: int = 100,
) -> Row:
    return {
        "orf_id": orf_id,
        "contig": contig,
        "strand": strand,
        "start_longest": str(start_longest),
        "end": str(end),
        "peptides": peptides,
    }


# ---------------------------------------------------------------- candidate_ids


def test_il_folded_peptide_on_contig_and_strand_matches() -> None:
    assert candidate_ids([hit("PEPTLDEK")], [orf(peptides="PEPTIDEK")]) == ["NOVEL_chr1_100_p"]
    assert candidate_ids([hit("PEPTIDEK")], [orf(peptides="PEPTLDEK")]) == ["NOVEL_chr1_100_p"]


def test_second_listed_peptide_matches() -> None:
    assert candidate_ids([hit("LLAMAR")], [orf()]) == ["NOVEL_chr1_100_p"]


def test_strand_mismatch_gives_none() -> None:
    """Kills the 'strand dropped from the owner key' mutant."""
    assert candidate_ids([hit(strand="-")], [orf(strand="+")]) == [None]


def test_contig_mismatch_gives_none() -> None:
    assert candidate_ids([hit(contig="plasmid")], [orf()]) == [None]


def test_cds_inframe_row_gives_none() -> None:
    """Kills the 'locus-class filter removed' mutant: same peptide, contig, strand and span."""
    assert candidate_ids([hit(locus_class="cds_inframe")], [orf()]) == [None]


def test_rna_overlap_and_na_locus_class_give_none() -> None:
    hits = [hit(locus_class="rna_overlap"), hit(locus_class=None)]
    assert candidate_ids(hits, [orf()]) == [None, None]


def test_every_novel_class_matches() -> None:
    hits = [hit(locus_class=c) for c in NOVEL_CLASSES]
    assert candidate_ids(hits, [orf()]) == ["NOVEL_chr1_100_p"] * len(NOVEL_CLASSES)


def test_first_orf_listing_the_peptide_wins() -> None:
    a = orf("NOVEL_chr1_100_p", start_longest=10, end=100)
    b = orf("NOVEL_chr1_120_p", start_longest=4, end=120)
    assert candidate_ids([hit()], [a, b]) == ["NOVEL_chr1_100_p"]
    assert candidate_ids([hit()], [b, a]) == ["NOVEL_chr1_120_p"]


def test_orf_without_peptides_or_id_owns_nothing() -> None:
    assert candidate_ids([hit()], [orf(peptides=None)]) == [None]
    assert candidate_ids([hit()], [{**orf(), "orf_id": None}]) == [None]


def test_unlisted_peptide_gives_none() -> None:
    assert candidate_ids([hit("QQQQQK")], [orf()]) == [None]


def test_one_value_per_hit_in_order() -> None:
    hits = [hit(), hit(strand="-"), hit("LLAMAR"), hit(locus_class="cds_inframe")]
    want = ["NOVEL_chr1_100_p", None, "NOVEL_chr1_100_p", None]
    assert candidate_ids(hits, [orf()]) == want
    assert candidate_ids([], [orf()]) == []
    assert candidate_ids(hits, []) == [None] * 4


# -- position within the ORF span (decision 11 nit; review N5): fails until positions are used


def test_hit_outside_the_orf_span_gives_none() -> None:
    """A second copy of the peptide (a repeat) on the same contig and strand, past the ORF."""
    assert candidate_ids([hit(start=5000, end=5023)], [orf(start_longest=10, end=100)]) == [None]


def test_hit_inside_and_outside_the_span_on_one_contig() -> None:
    hits = [hit(start=40, end=63), hit(start=400, end=423)]
    assert candidate_ids(hits, [orf()]) == ["NOVEL_chr1_100_p", None]


def test_minus_strand_span_runs_from_end_to_start_longest() -> None:
    # on '-' start_longest is the highest coordinate: the span is [end, start_longest]
    minus = orf("NOVEL_chr1_10_m", strand="-", start_longest=100, end=10)
    assert candidate_ids([hit(strand="-", start=40, end=63)], [minus]) == ["NOVEL_chr1_10_m"]
    assert candidate_ids([hit(strand="-", start=200, end=223)], [minus]) == [None]


def test_the_orf_covering_the_hit_wins_over_an_earlier_one_elsewhere() -> None:
    far = orf("NOVEL_chr1_9000_p", start_longest=8000, end=9000)
    near = orf("NOVEL_chr1_100_p", start_longest=10, end=100)
    assert candidate_ids([hit(start=40, end=63)], [far, near]) == ["NOVEL_chr1_100_p"]


@given(
    st.lists(
        st.tuples(
            st.sampled_from(["PEPTIDEK", "PEPTLDEK", "LLAMAR", "QQQK"]),
            st.sampled_from(["+", "-"]),
            st.sampled_from([*NOVEL_CLASSES, "cds_inframe", "rna_overlap", None]),
        ),
        max_size=8,
    )
)
def test_property_non_novel_rows_never_get_an_id(rows: list[tuple[str, str, str | None]]) -> None:
    hits = [hit(p, strand=s, locus_class=c) for p, s, c in rows]
    ids = candidate_ids(
        hits, [orf(), orf("NOVEL_chr1_10_m", strand="-", start_longest=100, end=10)]
    )
    assert len(ids) == len(hits)
    for h, oid in zip(hits, ids, strict=True):
        if h["locus_class"] not in NOVEL_CLASSES:
            assert oid is None
        if oid is not None:
            assert oid.endswith("_p" if h["strand"] == "+" else "_m")


# ---------------------------------------------------------------- backfill_candidate_ids

HITS_HEADER = [
    "peptide", "peptide_len", "source", "best_q_input", "pep_input", "score_input", "contig",
    "start", "end", "strand", "frame", "crosses_origin", "route", "n_diff", "diff_types",
    "placement_q", "length_class", "locus_class", "uniqueness", "gene", "novel_candidate_id",
    "flags",
]  # fmt: skip
ORFS_HEADER = ["orf_id", "class", "contig", "strand", "start_longest", "end", "peptides", "called"]


def _hits_line(h: Row, cand: str = "NA") -> str:
    cells = {k: "NA" for k in HITS_HEADER}
    cells.update({k: v if v is not None else "NA" for k, v in h.items()})
    cells.update(peptide_len=str(len(h["peptide"] or "")), route="exact", novel_candidate_id=cand)
    return "\t".join(cells[k] for k in HITS_HEADER)


def _orfs_line(o: Row) -> str:
    cells = {"class": "intergenic", "called": "true", **o}
    return "\t".join(str(cells[k]) if cells[k] is not None else "NA" for k in ORFS_HEADER)


def _write(out: Path, hits: list[str], orfs: list[Row] | None) -> Path:
    (out / "mapping").mkdir(parents=True)
    hp = out / "mapping" / "peptide_hits.tsv"
    hp.write_text("\n".join(["\t".join(HITS_HEADER), *hits]) + "\n", encoding="utf-8")
    if orfs is not None:
        (out / "tables").mkdir()
        lines = ["\t".join(ORFS_HEADER), *(_orfs_line(o) for o in orfs)]
        (out / "tables" / "novel_orfs.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return hp


def _cells(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]


def test_backfill_without_files_is_a_no_op(tmp_path: Path) -> None:
    assert backfill_candidate_ids(tmp_path) == 0
    hp = _write(tmp_path / "a", [_hits_line(hit())], None)
    before = hp.read_bytes()
    assert backfill_candidate_ids(tmp_path / "a") == 0
    assert hp.read_bytes() == before


def test_backfill_fills_matching_rows_and_changes_only_that_column(tmp_path: Path) -> None:
    rows = [
        _hits_line(hit()),
        _hits_line(hit(locus_class="cds_inframe")),
        _hits_line(hit("LLAMAR")),
    ]
    hp = _write(tmp_path, rows, [orf()])
    before = _cells(hp)
    assert backfill_candidate_ids(tmp_path) == 2
    after = _cells(hp)
    col = HITS_HEADER.index("novel_candidate_id")
    assert [r[col] for r in after[1:]] == ["NOVEL_chr1_100_p", "NA", "NOVEL_chr1_100_p"]
    for b, a in zip(before, after, strict=True):
        assert b[:col] + b[col + 1 :] == a[:col] + a[col + 1 :]
    assert hp.read_text(encoding="utf-8").endswith("\n")
    assert not list(hp.parent.glob(".tmp.*"))


def test_backfill_is_idempotent(tmp_path: Path) -> None:
    hp = _write(tmp_path, [_hits_line(hit()), _hits_line(hit(strand="-"))], [orf()])
    backfill_candidate_ids(tmp_path)
    once = hp.read_bytes()
    backfill_candidate_ids(tmp_path)
    assert hp.read_bytes() == once


def test_backfill_clears_a_stale_id(tmp_path: Path) -> None:
    """Decision 11 nit (review N4): every row is written, NA when it has no id, so a published
    file kept by a cached task (``-resume``) loses ids whose ORF is gone."""
    rows = [_hits_line(hit(), "NOVEL_chr1_100_p"), _hits_line(hit(strand="-"), "NOVEL_chr1_7_m")]
    hp = _write(tmp_path, rows, [orf("NOVEL_chr1_300_p", start_longest=10, end=300)])
    backfill_candidate_ids(tmp_path)
    col = HITS_HEADER.index("novel_candidate_id")
    assert [r[col] for r in _cells(hp)[1:]] == ["NOVEL_chr1_300_p", "NA"]


def test_backfill_ignores_a_table_without_the_needed_columns(tmp_path: Path) -> None:
    out = tmp_path
    (out / "mapping").mkdir()
    hp = out / "mapping" / "peptide_hits.tsv"
    hp.write_text("peptide\tcontig\nPEPTIDEK\tchr1\n", encoding="utf-8")
    (out / "tables").mkdir()
    (out / "tables" / "novel_orfs.tsv").write_text(
        "\n".join(["\t".join(ORFS_HEADER), _orfs_line(orf())]) + "\n", encoding="utf-8"
    )
    before = hp.read_bytes()
    assert backfill_candidate_ids(out) == 0
    assert hp.read_bytes() == before


def test_backfill_result_reads_back_as_a_table(tmp_path: Path) -> None:
    from giap.io.tsv import read_tsv

    hp = _write(tmp_path, [_hits_line(hit()), _hits_line(hit(strand="-"))], [orf()])
    backfill_candidate_ids(tmp_path)
    header, rows = read_tsv(hp)
    assert header == HITS_HEADER
    col = header.index("novel_candidate_id")
    assert json.dumps([r[col] for r in rows]) == json.dumps(["NOVEL_chr1_100_p", None])
