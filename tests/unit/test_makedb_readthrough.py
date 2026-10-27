"""``gippa makedb --patch`` with a ``cterm_extension`` row (M7-16b rule 13; ADR-0058 item 5; spec
examples 7 and 7b).

For a patch CDS with ``Note=cterm_extension`` and a ``locus_tag`` naming a CDS/sORF on the same
strand, the residue at that gene's terminal stop codon translates as ``X``. The codon must lie in
frame inside the patch interval, else ``ValueError`` naming the CDS. After that substitution any
remaining internal ``*`` raises as today ("patch CDS <id> has an internal stop codon"); the
terminal ``*`` is stripped as today; non-cterm patch rows are unchanged.

World: the spec's toy contig ``c1`` (linear, 126 nt, frame +0 ``MSDKPLTEAGRWQLNEFK*VPAGHTLDMRK*
AAGLKEPSTR*``) with ``G1`` = CDS ``+`` ``[0, 57)`` (TGA ``[54, 57)``), as a Bakta-like JSON; the
minus world is its reverse complement (``G1`` ``-`` ``[69, 126)``, stop ``[69, 72)``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from makedb_checks import by_id, fasta_records, id_map
from makedb_genomes import Contig, feature, write_bakta_json
from readthrough_toys import rc, toy_c1

from gippa.makedb.build import MakedbOptions, build_strain_db

G1_AA = "MSDKPLTEAGRWQLNEFK"
PATCHED = "MSDKPLTEAGRWQLNEFKXVPAGHTLDMRK"
N = 126


def _world(tmp: Path, strand: str = "+") -> Path:
    nt = toy_c1()[0]["c1"]
    if strand == "-":
        nt = rc(nt)
    c = Contig("c1", N, "linear", "chromosome", list(nt))
    s, e = (0, 57) if strand == "+" else (N - 57, N)
    f = feature("cds", "c1", s, e, strand, locus="G1", aa=G1_AA, product="toy protein")
    return write_bakta_json(tmp / "bakta.json", [c], [f])


def _patch(tmp: Path, lines: list[str]) -> Path:
    p = tmp / "annotation_patch.gff3"
    p.write_text("\n".join(["##gff-version 3", f"##sequence-region c1 1 {N}", *lines]) + "\n")
    return p


def _cds(
    start1: int, end: int, strand: str = "+", attrs: str = "ID=G1;locus_tag=G1;Note=cterm_extension"
) -> str:
    return f"c1\tGIPPA\tCDS\t{start1}\t{end}\t.\t{strand}\t0\t{attrs}"


def _build(tmp: Path, lines: list[str], strand: str = "+") -> Path:
    out = tmp / "DB"
    build_strain_db(
        MakedbOptions(
            annotation=_world(tmp, strand), genome=None, outdir=out, patch=_patch(tmp, lines)
        )
    )
    return out


def test_example_7_patch_record_reads_x_at_the_old_stop(tmp_path: Path) -> None:
    out = _build(tmp_path, [_cds(1, 93)])
    seqs = by_id(out)
    assert seqs["G1"] == G1_AA  # the Bakta protein is still written
    assert seqs["G1_patch"] == PATCHED
    headers = [h for h, _ in fasta_records(out / "strain_db.faa")]
    assert any(h.startswith("G1_patch ") for h in headers)
    (row,) = [r for r in id_map(out) if r["kind"] == "patch"]
    assert (
        row["db_id"],
        row["locus_tag"],
        row["contig"],
        row["start"],
        row["end"],
        row["strand"],
    ) == ("G1_patch", "G1", "c1", "1", "93", "+")


def test_example_7b_a_second_internal_stop_still_raises(tmp_path: Path) -> None:
    """``[0, 126) +`` also covers the TAA at p = 30: only the annotated stop becomes ``X``."""
    with pytest.raises(ValueError, match="patch CDS G1_patch has an internal stop codon"):
        _build(tmp_path, [_cds(1, 126)])


def test_minus_strand_patch_reads_x(tmp_path: Path) -> None:
    out = _build(tmp_path, [_cds(34, 126, "-")], strand="-")
    assert by_id(out)["G1_patch"] == PATCHED
    (row,) = [r for r in id_map(out) if r["kind"] == "patch"]
    assert (row["start"], row["end"], row["strand"]) == ("34", "126", "-")


def test_non_cterm_patch_rows_are_unchanged(tmp_path: Path) -> None:
    """Without ``Note=cterm_extension`` the annotated stop is an internal ``*`` as today."""
    with pytest.raises(ValueError, match="patch CDS G1_patch has an internal stop codon"):
        _build(tmp_path, [_cds(1, 93, attrs="ID=G1;locus_tag=G1;Note=nterm_extension")])


@pytest.mark.parametrize(
    ("start1", "end"),
    [
        (1, 45),  # the stop codon [54, 57) is not inside the patch interval
        (4, 54),  # ends before the codon
    ],
)
def test_stop_outside_the_patch_interval_raises(tmp_path: Path, start1: int, end: int) -> None:
    with pytest.raises(ValueError, match="G1_patch"):
        _build(tmp_path, [_cds(start1, end)])


def test_stop_out_of_frame_in_the_patch_raises(tmp_path: Path) -> None:
    """``[1, 94)``: the codon ``[54, 57)`` lies inside but (54 - 1) mod 3 != 0."""
    with pytest.raises(ValueError, match="G1_patch"):
        _build(tmp_path, [_cds(2, 94)])
