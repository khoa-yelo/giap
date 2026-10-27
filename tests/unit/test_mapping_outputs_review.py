"""M3-16 review pins (docs/specs/M3-16-mapping-outputs.md, "Lead decisions after code review", 1-3).

1. Origin-crossing frameshift placements: ``match_part`` children (and any feature) whose unwrapped
   start is >= L are shifted by -L; a child spanning the origin is handled like other
   origin-crossing features; GFF3 writing never fails on valid placements.
2. ``δ = -1``: the parent ``match`` carries ``frameshift_reread=<1-based pos>``; ``δ = +1`` carries
   ``frameshift_skipped=<1-based pos>``.
3. Mutant killers: ``match_part`` strand equals the parent's; tied sort keys in shuffled input give
   identical output; ``;``-joined ``flags``/``diff_types`` round-trip; ``[]`` is missing.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any

import pytest
from helpers.mapping_outputs_oracle import (
    bed12_core,
    bed_core_of_line,
    covered_mod,
    frameshift_base,
)
from hypothesis import given, settings
from hypothesis import strategies as st
from test_mapping_outputs import (  # type: ignore[import-not-found]
    FS_ROW,
    OUTPUT_FILES,
    attr,
    data_lines,
    gff_features,
    pep,
    prot,
    raw_tsv,
    write,
)

from gippa.model import Interval

L = 10_000
FS_PEP = "PEPTIDEKLLR"

# --------------------------------------------------- 1. origin-crossing frameshift placements

OC_FS_CASES = {
    # the lead's example: second block entirely beyond L
    "part2_beyond": (9990, 10027, "9990-10008;10009-10027"),
    # first block spans the origin, second entirely beyond L
    "part1_spans": (9999, 10037, "9999-10018;10019-10037"),
    # δ = -1 across the origin
    "reread_beyond": (9990, 10026, "9990-10008;10007-10026"),
}


def _oc_row(start: int, end: int, blocks: str) -> dict[str, Any]:
    return pep(
        start=start,
        end=end,
        route="tolerant",
        n_diff=1,
        diff_types=("frameshift",),
        blocks=blocks,
        frame=1,
    )


def _spans(blocks: str) -> list[tuple[int, int]]:
    return [(int(s), int(e)) for s, e in (p.split("-") for p in blocks.split(";"))]


@pytest.mark.parametrize("case", sorted(OC_FS_CASES))
def test_origin_crossing_frameshift_writes_all_files(tmp_path: Path, case: str) -> None:
    start, end, blocks = OC_FS_CASES[case]
    out = write(tmp_path, [_oc_row(start, end, blocks)])
    for name in OUTPUT_FILES:
        assert (out / name).is_file(), name


def test_origin_crossing_frameshift_lead_example_gff3(tmp_path: Path) -> None:
    out = write(tmp_path, [_oc_row(9990, 10027, "9990-10008;10009-10027")])
    feats = gff_features(out / "mappings.gff3")
    (parent,) = [f for f in feats if f.type == "nucleotide_to_protein_match"]
    parts = [f for f in feats if f.type == "match_part"]
    assert parent.interval == Interval("chr1", 9990, 10027, "+")  # parent keeps end > L
    # the child starting at 10009 (>= L) is shifted by -L
    assert sorted((p.interval.start, p.interval.end) for p in parts) == [(9, 27), (9990, 10008)]
    pid = attr(parent, "ID")
    assert all(attr(p, "Parent") == pid for p in parts)
    raw = [
        ln.split("\t")
        for ln in (out / "mappings.gff3").read_text(encoding="utf-8").splitlines()
        if "\tmatch_part\t" in ln
    ]
    assert sorted((c[3], c[4]) for c in raw) == [("10", "27"), ("9991", "10008")]


@pytest.mark.parametrize("case", sorted(OC_FS_CASES))
def test_origin_crossing_frameshift_children_cover_blocks(tmp_path: Path, case: str) -> None:
    start, end, blocks = OC_FS_CASES[case]
    out = write(tmp_path, [_oc_row(start, end, blocks)])
    feats = gff_features(out / "mappings.gff3")
    (parent,) = [f for f in feats if f.type == "nucleotide_to_protein_match"]
    parts = [f for f in feats if f.type == "match_part"]
    assert parts
    assert all(p.interval.start < L for p in parts)  # no child starts beyond the replicon
    assert all(attr(p, "Parent") == attr(parent, "ID") for p in parts)
    assert all(p.interval.strand == parent.interval.strand for p in parts)
    got = covered_mod([(p.interval.start, p.interval.end) for p in parts], L)
    assert got == covered_mod(_spans(blocks), L)


@pytest.mark.parametrize("case", sorted(OC_FS_CASES))
def test_origin_crossing_frameshift_bed12(tmp_path: Path, case: str) -> None:
    start, end, blocks = OC_FS_CASES[case]
    out = write(tmp_path, [_oc_row(start, end, blocks)])
    got = sorted(bed_core_of_line(f) for f in data_lines(out / "mappings.bed"))
    assert got == sorted(bed12_core("chr1", start, end, "+", FS_PEP, 500, L, blocks))
    probed = sorted(bed_core_of_line(f) for f in data_lines(out / "mappings.probed"))
    assert probed == got


# ------------------------------------------------------------ 2. frameshift odd-base attribute


@pytest.mark.parametrize(
    ("row", "key", "pos"),
    [
        # δ = +1, M3-08 example 1: blocks 2-20;21-39, 0-based base 20 skipped -> 1-based 21
        (FS_ROW, "frameshift_skipped", 21),
        # δ = -1, M3-08 example 6: blocks 2-20;19-37, 0-based base 19 read twice -> 1-based 20
        ({**FS_ROW, "end": 37, "blocks": "2-20;19-37"}, "frameshift_reread", 20),
        # δ = +1 on the minus strand (M3-08 example 5): 0-based base 21 skipped -> 1-based 22
        (
            {**FS_ROW, "start": 3, "end": 40, "strand": "-", "frame": -3, "blocks": "3-21;22-40"},
            "frameshift_skipped",
            22,
        ),
        # δ = -1 on the minus strand
        (
            {**FS_ROW, "start": 3, "end": 38, "strand": "-", "frame": -3, "blocks": "3-21;20-38"},
            "frameshift_reread",
            21,
        ),
    ],
    ids=["plus_skipped", "plus_reread", "minus_skipped", "minus_reread"],
)
def test_frameshift_odd_base_attribute(
    tmp_path: Path, row: dict[str, Any], key: str, pos: int
) -> None:
    assert frameshift_base(str(row["blocks"])) == (key, pos)  # oracle agrees with the pin
    out = write(tmp_path, [pep(**row)])
    feats = gff_features(out / "mappings.gff3")
    (parent,) = [f for f in feats if f.type == "nucleotide_to_protein_match"]
    other = "frameshift_skipped" if key == "frameshift_reread" else "frameshift_reread"
    assert attr(parent, key) == str(pos)
    assert attr(parent, other) is None
    for p in (f for f in feats if f.type == "match_part"):
        assert attr(p, key) is None


def test_no_frameshift_attribute_without_blocks(tmp_path: Path) -> None:
    out = write(tmp_path, [pep()])
    (f,) = gff_features(out / "mappings.gff3")
    assert attr(f, "frameshift_reread") is None
    assert attr(f, "frameshift_skipped") is None


# --------------------------------------------------------------------- 3. mutant killers


@pytest.mark.parametrize(
    ("strand", "frame", "start", "end", "blocks"),
    [("-", -3, 3, 40, "3-21;22-40"), ("+", 3, 2, 39, "2-20;21-39")],
)
def test_match_part_strand_equals_parent(
    tmp_path: Path, strand: str, frame: int, start: int, end: int, blocks: str
) -> None:
    row = {**FS_ROW, "strand": strand, "frame": frame, "start": start, "end": end}
    out = write(tmp_path, [pep(**{**row, "blocks": blocks})])
    feats = gff_features(out / "mappings.gff3")
    parts = [f for f in feats if f.type == "match_part"]
    assert len(parts) == 2
    assert [p.interval.strand for p in parts] == [strand, strand]


def _tied_peps() -> list[dict[str, Any]]:
    # Identical contig/start/end/strand/peptide: only the full-row tie-break orders them.
    return [
        pep(score_input=55.5, placement_q=0.002, gene="GIPPA_00010"),
        pep(score_input=1.0, placement_q=0.004, gene="GIPPA_00020", flags="low_score"),
        pep(score_input=30.0, placement_q=0.001, gene=None, uniqueness="shared"),
    ]


def _tied_prots() -> list[dict[str, Any]]:
    return [
        prot(rank=1, identity=1.0),
        prot(rank=2, identity=0.97, **{"class": "near_identical"}),
    ]


def _all_bytes(out: Path) -> dict[str, bytes]:
    return {name: (out / name).read_bytes() for name in OUTPUT_FILES}


def test_tied_sort_keys_order_independent(tmp_path: Path) -> None:
    peps, prots = _tied_peps(), _tied_prots()
    a = _all_bytes(write(tmp_path / "a", peps, prots))
    b = _all_bytes(write(tmp_path / "b", peps[::-1], prots[::-1]))
    assert a == b


@settings(max_examples=25, deadline=None)
@given(
    pep_perm=st.permutations(range(3)),
    prot_perm=st.permutations(range(2)),
)
def test_tied_sort_keys_any_permutation(pep_perm: list[int], prot_perm: list[int]) -> None:
    peps, prots = _tied_peps(), _tied_prots()
    with tempfile.TemporaryDirectory() as d:
        ref = _all_bytes(write(Path(d) / "ref", peps, prots))
        got = _all_bytes(
            write(Path(d) / "perm", [peps[i] for i in pep_perm], [prots[i] for i in prot_perm])
        )
    assert got == ref


@pytest.mark.parametrize("col", ["flags", "diff_types"])
def test_semicolon_joined_string_round_trips(tmp_path: Path, col: str) -> None:
    value = "substitution;swap" if col == "diff_types" else "low_score;near_origin"
    out = write(tmp_path, [pep(route="tolerant", n_diff=2, **{col: value})])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert r[col] == value
    (f,) = gff_features(out / "mappings.gff3")
    got = f.attributes[col]
    assert not isinstance(got, str)  # a GFF3 multi-value, not one escaped "a%3Bb" value
    assert list(got) == value.split(";")
    assert ";".join(got) == value


@pytest.mark.parametrize("empty", [(), []], ids=["tuple", "list"])
def test_empty_sequence_is_missing(tmp_path: Path, empty: Any) -> None:
    out = write(tmp_path, [pep(diff_types=empty)])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert r["diff_types"] == "NA"
    (f,) = gff_features(out / "mappings.gff3")
    assert "diff_types" not in f.attributes
    assert "diff_types=" not in (out / "mappings.gff3").read_text(encoding="utf-8")
