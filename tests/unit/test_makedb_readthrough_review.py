"""``gippa makedb --patch`` ``cterm_extension`` rows: the review-335 guards and error cases (M7-16b
rule 13; code review of PR #335, must-fix items 2 and 3 and the ``d < 3`` nit).

* The residue that becomes ``X`` must be a stop: a ``cterm_extension`` patch whose gene ends on a
  sense codon (a ``not_a_stop`` partial) raises ``ValueError`` naming the patch CDS instead of
  writing ``X`` over a sense residue.
* A ``locus_tag`` naming an opposite-strand gene, an unknown tag, a pseudogene or a tRNA raises
  the "cterm_extension without a <strand> CDS" error; an out-of-frame stop and a patch that starts
  on the stop codon itself (``d < 3``) raise the "not in frame" error (matched on its own text, so
  the internal-stop error can't satisfy the test).
* On a circular replicon an origin-crossing patch puts the ``X`` at the gene's stop, both when the
  gene lies after the origin (the codon offset wraps) and when the gene crosses it.

World: the spec's toy contig ``c1`` (126 nt, frame +0 ``MSDKPLTEAGRWQLNEFK*VPAGHTLDMRK*
AAGLKEPSTR*``; TGA ``[54, 57)``, TAA ``[90, 93)``), as in ``test_makedb_readthrough.py``.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from makedb_checks import by_id
from makedb_genomes import Contig, feature, write_bakta_json
from readthrough_toys import rc, toy_c1

from gippa.makedb.build import MakedbOptions, build_strain_db

G1_AA = "MSDKPLTEAGRWQLNEFK"
PATCHED = "MSDKPLTEAGRWQLNEFKXVPAGHTLDMRK"
N = 126


def _nt(strand: str = "+", shift: int = 0) -> str:
    """Toy ``c1`` (reverse-complemented for ``-``), rotated so toy position p sits at
    ``(p - shift) % N``."""
    nt = toy_c1()[0]["c1"]
    if strand == "-":
        nt = rc(nt)
    return nt[shift:] + nt[:shift]


def _g1(strand: str = "+") -> dict[str, Any]:
    s, e = (0, 57) if strand == "+" else (N - 57, N)
    return feature("cds", "c1", s, e, strand, locus="G1", aa=G1_AA, product="toy protein")


def _build(
    tmp: Path,
    lines: list[str],
    features: list[dict[str, Any]],
    nt: str | None = None,
    topology: str = "linear",
) -> Path:
    c = Contig("c1", N, topology, "chromosome", list(nt or _nt()))
    ann = write_bakta_json(tmp / "bakta.json", [c], features)
    patch = tmp / "annotation_patch.gff3"
    patch.write_text("\n".join(["##gff-version 3", f"##sequence-region c1 1 {N}", *lines]) + "\n")
    out = tmp / "DB"
    build_strain_db(MakedbOptions(annotation=ann, genome=None, outdir=out, patch=patch))
    return out


def _cds(start1: int, end: int, strand: str = "+", tag: str = "G1", pid: str | None = None) -> str:
    attrs = f"ID={pid or tag};locus_tag={tag};Note=cterm_extension"
    return f"c1\tGIPPA\tCDS\t{start1}\t{end}\t.\t{strand}\t0\t{attrs}"


# ---------------------------------------------------------------- must-fix 2: the residue guard


@pytest.mark.parametrize(
    ("gene_end", "patch_end"),
    [
        (90, 90),  # G2 [57, 90) ends on K (p = 29); the patch ends with it
        (84, 90),  # G2 [57, 84) ends on M (p = 27); the patch runs one codon past it
    ],
)
def test_patch_whose_gene_ends_on_a_sense_codon_raises(
    tmp_path: Path, gene_end: int, patch_end: int
) -> None:
    """Today the sense residue silently becomes ``X`` (``MPAGHTLDMRX`` / ``MPAGHTLDXRK``): no
    ``*`` is left, so no other check fires. The guard must raise and name the patch CDS."""
    aa = "VPAGHTLDMRK"[: (gene_end - 57) // 3]
    g2 = feature("cds", "c1", 57, gene_end, "+", locus="G2", aa=aa, product="partial")
    with pytest.raises(ValueError, match="G2_patch"):
        _build(tmp_path, [_cds(58, patch_end, tag="G2")], [_g1(), g2])


def test_guard_does_not_touch_a_real_stop(tmp_path: Path) -> None:
    """Control for the guard: the example-7 patch over G1's TGA still reads ``X``."""
    out = _build(tmp_path, [_cds(1, 93)], [_g1()])
    assert by_id(out)["G1_patch"] == PATCHED


# ---------------------------------------------------------------- must-fix 3: the tag errors


def _no_cds(rid: str, strand: str, tag: str) -> str:
    return re.escape(f"patch CDS {rid}: cterm_extension without a {strand} CDS {tag!r}")


def test_opposite_strand_locus_tag_raises(tmp_path: Path) -> None:
    """G1 is ``+``; a ``-`` patch naming it is rejected by the strand check itself."""
    with pytest.raises(ValueError, match=_no_cds("G1_patch", "-", "G1")):
        _build(tmp_path, [_cds(1, 93, strand="-")], [_g1()])


def test_opposite_strand_locus_tag_raises_minus_world(tmp_path: Path) -> None:
    """Mirror: G1 ``-`` ``[69, 126)``; a ``+`` patch naming it."""
    with pytest.raises(ValueError, match=_no_cds("G1_patch", "+", "G1")):
        _build(tmp_path, [_cds(34, 126, strand="+")], [_g1("-")], nt=_nt("-"))


def test_unknown_locus_tag_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=_no_cds("P9", "+", "NOPE")):
        _build(tmp_path, [_cds(1, 93, tag="NOPE", pid="P9")], [_g1()])


def test_pseudogene_locus_tag_raises(tmp_path: Path) -> None:
    """A pseudogene is not a CDS/sORF, even when its interval ends on G1's stop."""
    ps = feature("cds", "c1", 0, 57, "+", locus="PS1", aa=G1_AA, pseudo=True)
    with pytest.raises(ValueError, match=_no_cds("PS1", "+", "PS1")):
        _build(tmp_path, [_cds(1, 93, tag="PS1")], [ps])


def test_trna_locus_tag_raises(tmp_path: Path) -> None:
    trna = feature("tRNA", "c1", 93, 117, "+", locus="T1", product="tRNA-Ala", amino_acid="Ala")
    with pytest.raises(ValueError, match=_no_cds("T1", "+", "T1")):
        _build(tmp_path, [_cds(1, 93, tag="T1")], [_g1(), trna])


# ---------------------------------------------------------------- must-fix 3: the frame errors


def _not_in_frame(rid: str, tag: str = "G1") -> str:
    return re.escape(f"patch CDS {rid}: the stop of {tag} is not in frame inside the CDS")


def test_out_of_frame_stop_raises_the_frame_error(tmp_path: Path) -> None:
    """``[1, 94)``: the codon ``[54, 57)`` lies inside but (54 - 1) mod 3 = 2. The shifted frame
    also has internal stops, so the test matches the frame error's own text."""
    with pytest.raises(ValueError, match=_not_in_frame("G1_patch")):
        _build(tmp_path, [_cds(2, 94)], [_g1()])


@pytest.mark.parametrize("strand", ["+", "-"])
def test_patch_starting_on_the_stop_codon_raises(tmp_path: Path, strand: str) -> None:
    """``d = 0`` (the patch's first codon is the stop): ``d < 3`` raises. ``+``: ``[54, 93)``;
    ``-``: ``[33, 72)`` with the stop ``[69, 72)`` at its 5' end."""
    line = _cds(55, 93) if strand == "+" else _cds(34, 72, strand="-")
    nt = _nt(strand)
    with pytest.raises(ValueError, match=_not_in_frame("G1_patch")):
        _build(tmp_path, [line], [_g1(strand)], nt=nt)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_patch_one_codon_before_the_stop_reads_x_second(tmp_path: Path, strand: str) -> None:
    """``d = 3``, the smallest accepted offset: ``M`` (the start codon read as M) then ``X``."""
    line = _cds(52, 93) if strand == "+" else _cds(34, 75, strand="-")
    out = _build(tmp_path, [line], [_g1(strand)], nt=_nt(strand))
    assert by_id(out)["G1_patch"] == "MXVPAGHTLDMRK"


# ---------------------------------------------------------------- must-fix 3: circular origin


def test_origin_crossing_patch_with_the_gene_after_the_origin(tmp_path: Path) -> None:
    """Toy rotated by 26 on a circular ``c1``: toy p sits at ``(p - 26) % 126``. Gene G1' =
    toy ``[27, 57)`` (``GRWQLNEFK`` + TGA) lies at ``[1, 31)``, its stop at ``[28, 31)``. The
    patch toy ``[0, 93)`` is ``[100, 193)`` unwrapped, so the codon offset ``28 - 100`` wraps
    to 54."""
    g = feature("cds", "c1", 1, 31, "+", locus="G1", aa="GRWQLNEFK", product="toy protein")
    out = _build(tmp_path, [_cds(101, 193)], [g], nt=_nt("+", 26), topology="circular")
    assert by_id(out)["G1_patch"] == PATCHED


def test_origin_crossing_patch_with_the_gene_across_the_origin(tmp_path: Path) -> None:
    """Toy rotated by 45: G1 toy ``[0, 57)`` is ``[81, 138)`` unwrapped (stop ``[9, 12)``) and the
    patch toy ``[0, 93)`` is ``[81, 174)``."""
    g = feature("cds", "c1", 81, 138, "+", locus="G1", aa=G1_AA, product="toy", length=N)
    out = _build(tmp_path, [_cds(82, 174)], [g], nt=_nt("+", 45), topology="circular")
    assert by_id(out)["G1_patch"] == PATCHED


def test_origin_crossing_patch_minus_strand(tmp_path: Path) -> None:
    """Minus world (G1 ``-`` ``[69, 126)``, stop ``[69, 72)``) rotated by 60: G1 sits at
    ``[9, 66)``, its stop at ``[9, 12)``, and the patch ``[33, 126)`` at ``[99, 192)`` unwrapped,
    so the codon offset ``189 - 9`` wraps to 54."""
    g = feature("cds", "c1", 9, 66, "-", locus="G1", aa=G1_AA, product="toy")
    out = _build(tmp_path, [_cds(100, 192, strand="-")], [g], nt=_nt("-", 60), topology="circular")
    assert by_id(out)["G1_patch"] == PATCHED
