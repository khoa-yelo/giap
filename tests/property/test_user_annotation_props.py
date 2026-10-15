"""Hypothesis properties for ``--annotation`` GFF3/GBFF reading and the id map (M2-11).

docs/specs/M2-11-user-annotation-id-map.md; GIAP_SPEC.md §5.3, §6.5:

* GFF3 written by an "other annotator" (origin-crossing CDS as two rows sharing one ID, or one
  unwrapped row) reads back to the same 0-based intervals;
* GBFF written by Biopython (``join`` across the origin) reads back to the same intervals;
* ``build_id_map`` against genes on the same coordinates is all ``identical`` with the right locus
  tags, one row per user CDS, and invariant to input order; a shifted copy is all ``unmatched``;
* ``write_bakta_regions`` reads back (M1-05 ``read_gff3``) to the same non-crossing intervals;
* ``write_id_map`` is 1-based closed and ``NA`` exactly for unmatched rows.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from hypothesis import given, settings
from hypothesis import strategies as st

from giap.annotate.user_annotation import (
    UserCds,
    build_id_map,
    read_user_cds,
    write_bakta_regions,
    write_id_map,
)
from giap.io.gff3 import read_gff3
from giap.model import Gene, Interval, Replicon

FAST = settings(max_examples=60, deadline=None)

CHR = Replicon("chr1", 20_000, "circular", "chromosome")
PL = Replicon("plasmid_1", 6_000, "circular", "plasmid")
LIN = Replicon("linear_1", 4_000, "linear", "unknown")
REPS = {r.id: r for r in (CHR, PL, LIN)}
IDENT = {k: k for k in REPS}
ORDER = list(REPS)  # input FASTA order (build_id_map requires it)


@st.composite
def cds_sets(draw: st.DrawFn, allow_origin: bool = True) -> list[UserCds]:
    """Distinct-coordinate CDS with unique IDs; some cross the origin of circular replicons."""
    n = draw(st.integers(0, 12))
    seen: set[tuple[str, int, int, str]] = set()
    out: list[UserCds] = []
    for i in range(n):
        rep = draw(st.sampled_from(list(REPS.values())))
        length = 3 * draw(st.integers(10, 300))
        cross = allow_origin and rep.topology == "circular" and draw(st.booleans())
        if cross:
            start = draw(st.integers(rep.length - length + 1, rep.length - 1))
        else:
            start = draw(st.integers(0, rep.length - length))
        strand = draw(st.sampled_from(["+", "-"]))
        key = (rep.id, start, start + length, strand)
        if key in seen:
            continue
        seen.add(key)
        iv = Interval(rep.id, start, start + length, strand)  # type: ignore[arg-type]
        out.append(UserCds(f"SRC_{i:03d}", iv, draw(st.booleans())))
    return out


def _gff_rows(c: UserCds, split: bool) -> list[str]:
    iv = c.interval
    rep = REPS[iv.contig]
    attrs = f"ID=cds-{c.source_id};locus_tag={c.source_id}" + (";pseudo=true" if c.pseudo else "")
    parts = [(iv.start + 1, iv.end)]
    if split and iv.end > rep.length:
        parts = [(iv.start + 1, rep.length), (1, iv.end - rep.length)]
    return [
        "\t".join([iv.contig, "other", "CDS", str(a), str(b), ".", iv.strand, "0", attrs])
        for a, b in parts
    ]


def _key(c: UserCds) -> tuple[str, int, int, str, str]:
    return (c.interval.contig, c.interval.start, c.interval.end, c.interval.strand, c.source_id)


def _gene(tag: str, iv: Interval) -> Gene:
    return Gene(tag, "cds", iv, "M" * max(1, iv.length_bp // 3 - 1), None, "ATG")


@FAST
@given(cds_sets(), st.booleans(), st.randoms(use_true_random=False))
def test_gff3_read_round_trip(cds: list[UserCds], split: bool, rnd: object) -> None:
    lines = [ln for c in cds for ln in _gff_rows(c, split)]
    rnd.shuffle(lines)  # type: ignore[attr-defined]
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "user.gff3"
        p.write_text("\n".join(["##gff-version 3", *lines]) + "\n")
        got = read_user_cds(p, "gff3", REPS, IDENT)
    assert sorted(got, key=_key) == sorted(cds, key=_key)


@FAST
@given(cds_sets())
def test_gbff_read_round_trip(cds: list[UserCds]) -> None:
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqIO import write
    from Bio.SeqRecord import SeqRecord

    records = []
    for rep in REPS.values():
        rec = SeqRecord(Seq("A" * rep.length), id=rep.id, name=rep.id, description="test")
        rec.annotations = {"molecule_type": "DNA", "topology": rep.topology}
        for c in cds:
            iv = c.interval
            if iv.contig != rep.id:
                continue
            s = 1 if iv.strand == "+" else -1
            if iv.end > rep.length:
                parts = [
                    FeatureLocation(iv.start, rep.length, s),
                    FeatureLocation(0, iv.end - rep.length, s),
                ]
                loc = CompoundLocation(parts if s == 1 else parts[::-1])
            else:
                loc = FeatureLocation(iv.start, iv.end, s)
            q: dict[str, list[str]] = {"locus_tag": [c.source_id]}
            if c.pseudo:
                q["pseudo"] = [""]
            rec.features.append(SeqFeature(loc, type="CDS", qualifiers=q))
        records.append(rec)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "user.gbff"
        write(records, p, "genbank")
        got = read_user_cds(p, "gbff", REPS, IDENT)
    assert sorted(got, key=_key) == sorted(cds, key=_key)


@FAST
@given(cds_sets(), st.randoms(use_true_random=False))
def test_id_map_identical_on_same_coordinates(cds: list[UserCds], rnd: object) -> None:
    genes = [_gene(f"LT_{i:05d}", c.interval) for i, c in enumerate(cds)]
    want = {c.source_id: g.locus_tag for c, g in zip(cds, genes, strict=True)}
    rows = build_id_map(cds, genes, contig_order=ORDER)
    assert len(rows) == len(cds)
    assert all(r.match == "identical" for r in rows)
    assert {r.source_id: r.locus_tag for r in rows} == want
    for r in rows:  # coordinates are the user CDS's, 0-based
        c = next(x for x in cds if x.source_id == r.source_id)
        assert (r.contig, r.start, r.end, r.strand) == (
            c.interval.contig,
            c.interval.start,
            c.interval.end,
            c.interval.strand,
        )
    shuffled_c, shuffled_g = list(cds), list(genes)
    rnd.shuffle(shuffled_c)  # type: ignore[attr-defined]
    rnd.shuffle(shuffled_g)  # type: ignore[attr-defined]
    assert build_id_map(shuffled_c, shuffled_g, contig_order=ORDER) == rows


@FAST
@given(cds_sets(), st.integers(1, 2))
def test_id_map_shifted_start_is_unmatched(cds: list[UserCds], shift: int) -> None:
    genes = [
        _gene(f"LT_{i:05d}", Interval(c.interval.contig, c.interval.start + shift,
                                      c.interval.end, c.interval.strand))
        for i, c in enumerate(cds)
        if c.interval.length_bp > shift
    ]  # fmt: skip
    rows = build_id_map(cds, genes, contig_order=ORDER)
    assert len(rows) == len(cds)
    assert all(r.match == "unmatched" and r.locus_tag is None for r in rows)


@FAST
@given(cds_sets(allow_origin=False))
def test_regions_round_trip_m1_05_reader(cds: list[UserCds]) -> None:
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "regions.gff3"
        write_bakta_regions(p, cds, list(REPS.values()))
        text = p.read_text()
        feats = read_gff3(p, REPS)
    assert text.startswith("##gff-version 3")
    assert {f.type for f in feats} <= {"CDS"}
    got = sorted((str(f.attributes["ID"]), f.interval) for f in feats)
    assert got == sorted(((c.source_id, c.interval) for c in cds), key=lambda x: x[0])


@FAST
@given(cds_sets(), st.data())
def test_id_map_tsv_one_based_and_na(cds: list[UserCds], data: st.DataObject) -> None:
    keep = [c for c in cds if data.draw(st.booleans())]
    genes = [_gene(f"LT_{i:05d}", c.interval) for i, c in enumerate(keep)]
    rows = build_id_map(cds, genes, contig_order=ORDER)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "id_map.tsv"
        write_id_map(p, rows)
        lines = p.read_text().splitlines()
    assert lines[0].split("\t") == [
        "source_id", "locus_tag", "contig", "start", "end", "strand", "match"
    ]  # fmt: skip
    assert len(lines) == len(rows) + 1
    for r, line in zip(rows, lines[1:], strict=True):
        cells = line.split("\t")
        assert cells == [
            r.source_id,
            "NA" if r.locus_tag is None else r.locus_tag,
            r.contig,
            str(r.start + 1),
            str(r.end),
            r.strand,
            r.match,
        ]
        assert (cells[1] == "NA") == (cells[6] == "unmatched")
