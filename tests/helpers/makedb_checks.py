"""Readers and shared assertions for the `giap makedb` outputs (issue M3-23 tests)."""

from __future__ import annotations

from pathlib import Path

from makedb_oracle import classified_orfs, nterm_entries, nterm_header_options, orf_header

from giap.annotate.bakta import parse_bakta_json
from giap.io.tsv import read_tsv

ID_MAP_COLUMNS = ["db_id", "kind", "locus_tag", "class", "contig", "start", "end", "strand"]


def fasta_records(path: Path) -> list[tuple[str, str]]:
    """(header without '>', sequence) in file order; checks the line layout on the way."""
    text = path.read_text()
    assert text == "" or (text.endswith("\n") and not text.endswith("\n\n"))
    out: list[tuple[str, str]] = []
    for line in text.split("\n")[:-1]:
        if line.startswith(">"):
            out.append((line[1:], ""))
        else:
            assert out, "sequence before the first header"
            assert 0 < len(line) <= 60, f"line of {len(line)} residues (wrap at 60)"
            out[-1] = (out[-1][0], out[-1][1] + line)
    return out


def wrapped_exactly(path: Path) -> None:
    """Every sequence line but a record's last holds exactly 60 residues."""
    lines = path.read_text().split("\n")[:-1]
    for i, line in enumerate(lines):
        if line.startswith(">"):
            continue
        nxt = lines[i + 1] if i + 1 < len(lines) else ">"
        if not nxt.startswith(">"):
            assert len(line) == 60


def db_id(header: str) -> str:
    return header.split(" ", 1)[0]


def id_map(outdir: Path) -> list[dict[str, str | None]]:
    cols, rows = read_tsv(outdir / "id_map.tsv")
    assert cols == ID_MAP_COLUMNS
    return [dict(zip(cols, r, strict=True)) for r in rows]


def by_id(outdir: Path) -> dict[str, str]:
    return {db_id(h): s for h, s in fasta_records(outdir / "strain_db.faa")}


def orf_ids(outdir: Path) -> list[str]:
    return [i for i in by_id(outdir) if i.startswith("GIAPORF_")]


def assert_orfs_match_oracle(
    out: Path,
    annotation: Path,
    genome: dict[str, str],
    min_aa: int,
    skip: tuple[str, int, int] | None = None,
) -> None:
    """Non-nterm ORFs equal the oracle (ids, sequences, order); nterm entries match per locus.

    ``skip`` = (contig, start0, end0): a region with ``N`` whose ORF treatment is not asserted.
    """
    gm = parse_bakta_json(annotation).model

    def skipped(rid: str) -> bool:
        if skip is None:
            return False
        *_, contig, start1, end1, _ = rid.split("_")
        return contig == skip[0] and int(start1) - 1 < skip[2] and int(end1) > skip[1]

    want = [
        (orf_header(c, o.contig, o.start, o.end, o.strand), o.aa)
        for c, o in classified_orfs(gm, genome, min_aa)
    ]
    got = [
        (i, s)
        for i, s in by_id(out).items()
        if i.startswith("GIAPORF_") and "_nterm_extension_" not in i
    ]
    assert [x for x in got if not skipped(x[0])] == [x for x in want if not skipped(x[0])]
    nterms, undefined = nterm_entries(gm, genome)
    rows = {r["db_id"]: r for r in id_map(out)}
    got_n: dict[str | None, tuple[str, str]] = {}
    for i, s in by_id(out).items():
        if i.startswith("GIAPORF_nterm_extension_"):
            tag = rows[i]["locus_tag"]
            assert tag not in got_n, f"two nterm_extension entries for {tag}"
            got_n[tag] = (i, s)
    assert set(got_n) - undefined == {n.locus_tag for n in nterms}
    for n in nterms:
        rid, seq = got_n[n.locus_tag]
        assert rid in nterm_header_options(n, gm), (rid, nterm_header_options(n, gm))
        assert seq == n.aa, n.locus_tag
