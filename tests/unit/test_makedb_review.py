"""`gippa makedb`: lead decisions 5-9 after code review of PR #79 (docs/specs/M3-23-makedb.md).

* 5: the junction part of an ``nterm_extension`` entry (first 30 codons of the CDS) is translated
  codon by codon from the genome (GTG -> V), not copied from Bakta's ``aa``.
* 6: the borrowed private helpers become public APIs: ``LocusClassifier(gm, genome).classify(iv)``
  and ``in_frame`` in ``gippa.mapping.loci``; ``translate11_text(nt, offset)`` and the public
  ambiguity sentinel ``AMBIG`` in ``gippa.mapping.sixframe`` (M3-01); ``gippa.makedb.build`` imports
  no private name from another module.
* 7: exit 1 (ValueError) on a patch CDS whose translation has an internal ``*`` and on an
  annotation JSON without features (``{}``).
* 8: the id map is ``db_id_map.tsv`` (no ``id_map.tsv``).
* 9: one test per mutant that survived the review (mutant names of the reviewer's script in
  brackets): nterm interval pinned [nterm_span0]; ``X`` in the upstream part
  [upstream_x_continue]; ``*`` stripping for contaminants [cont_rstrip] and for the fallback
  translation [protein_rstrip]; ``--genome`` length mismatch [genome_len_check]; duplicate ORFs
  at the origin [circ_dup]; in-frame drop across the origin [inframe_noshift]; an upstream walk
  that finds no stop within one circular turn [upstream_circ_term].
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path
from typing import Any

import pytest

import gippa
from gippa.annotate.bakta import parse_bakta_json
from gippa.exitcodes import ExitCode
from gippa.mapping.loci import classify_locus
from gippa.model.coords import Interval

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from makedb_checks import by_id, id_map, orf_ids
from makedb_genomes import (
    G1,
    G3,
    Contig,
    ExampleGenome,
    build_example,
    feature,
    stretch,
    write_bakta_json,
    write_fasta,
)
from makedb_oracle import revcomp, translate

# ---------------------------------------------------------------- helpers


def _api() -> Any:
    from gippa.makedb import build

    return build


def build_db(annotation: Path, outdir: Path, **kw: Any) -> Path:
    api = _api()
    api.build_strain_db(
        api.MakedbOptions(annotation=annotation, genome=kw.pop("genome", None), outdir=outdir, **kw)
    )
    return outdir


def _main(argv: list[str]) -> int:
    from gippa.cli import main

    with pytest.raises(SystemExit) as info:
        main(argv)
    code = info.value.code
    return 0 if code is None else int(code)


@pytest.fixture(scope="module")
def example() -> ExampleGenome:
    return build_example()


@pytest.fixture
def example_json(tmp_path: Path, example: ExampleGenome) -> Path:
    return write_bakta_json(tmp_path / "bakta.json", example.contigs, example.features)


def _one_contig_db(
    tmp: Path, contig: Contig, feats: list[dict[str, Any]], **kw: Any
) -> tuple[Path, Path]:
    annotation = write_bakta_json(tmp / "bakta.json", [contig], feats)
    return annotation, build_db(annotation, tmp / "DB", **kw)


def _nterm(out: Path) -> dict[str, tuple[str, str]]:
    """locus_tag -> (db_id, sequence) of the nterm_extension entries."""
    seqs = by_id(out)
    return {
        r["locus_tag"] or "": (r["db_id"] or "", seqs[r["db_id"] or ""])
        for r in id_map(out)
        if r["class"] == "nterm_extension"
    }


# ---------------------------------------------------------------- decision 5: junction


def test_junction_is_translated_from_the_genome_gtg_stays_v(tmp_path: Path) -> None:
    c = Contig("c1", 660, "circular")
    up, body = stretch(8, 0), stretch(40, 3)
    c.put(297, "TAA" + up + "GTG" + body + "TAA")  # CDS + [324, 450), upstream [300, 324)
    aa = "M" + translate(body)  # Bakta's aa starts with M
    feats = [feature("cds", "c1", 324, 450, "+", locus="R_00005", aa=aa, start_type="GTG")]
    _, out = _one_contig_db(tmp_path, c, feats, six_frame_orfs=True)
    rid, seq = _nterm(out)["R_00005"]
    assert rid == "GIPPAORF_nterm_extension_c1_301_414_p"
    assert seq == translate(up) + "V" + translate(body)[:29]
    assert seq[8] == "V"  # not the M of Bakta's aa


def test_junction_minus_strand_ttg_stays_l(tmp_path: Path) -> None:
    c = Contig("c1", 660, "circular")
    up, body = stretch(9, 2), stretch(35, 5)
    gene = "TTG" + body + "TAA"  # 37 codons on the minus strand
    c.put_minus(200, "TAA" + up + gene)  # minus string occupies [200, 341): CDS [200, 311)
    aa = "M" + translate(body)
    feats = [feature("cds", "c1", 200, 311, "-", locus="R_00010", aa=aa, start_type="TTG")]
    _, out = _one_contig_db(tmp_path, c, feats, six_frame_orfs=True)
    rid, seq = _nterm(out)["R_00010"]
    assert rid == "GIPPAORF_nterm_extension_c1_222_338_m"  # [221, 338): 30 codons + upstream
    assert seq == translate(up) + "L" + translate(body)[:29]
    assert translate(revcomp(c.text[200:311]))[0] == "L"  # precondition


# ---------------------------------------------------------------- decision 6: public APIs


def test_locus_classifier_is_public_and_matches_classify_locus(
    example: ExampleGenome, example_json: Path
) -> None:
    from gippa.mapping.loci import LocusClassifier

    gm = parse_bakta_json(example_json).model
    genome = example.genome
    clf = LocusClassifier(gm, genome)
    ivs = [
        Interval("chr1", 1200, 1290, "+"),  # intergenic
        Interval("chr1", 300, 390, "-"),  # antisense
        Interval("chr1", 603, 700, "+"),  # cds_inframe in G1
        Interval("chr1", 604, 700, "+"),  # cds_outframe
        Interval("chr1", 564, 650, "+"),  # nterm_extension of G1
        Interval("chr1", 2490, 2580, "+"),  # rna_overlap
        Interval("chr1", 2731, 2821, "+"),  # pseudogene
        Interval("chr1", 3240, 3390, "+"),  # origin-crossing
        Interval("lin1", 585, 660, "+"),
    ]
    for iv in ivs:
        assert clf.classify(iv) == classify_locus(iv, gm, genome), iv


def test_in_frame_is_public() -> None:
    from gippa.mapping.loci import in_frame

    plus = Interval("chr1", 600, 900, "+")
    assert in_frame(603, 700, plus) and in_frame(561, 600, plus)
    assert not in_frame(604, 700, plus) and not in_frame(605, 700, plus)
    minus = Interval("chr1", 2000, 2300, "-")
    assert in_frame(2100, 2297, minus) and in_frame(2300, 2330, minus)
    assert not in_frame(2100, 2298, minus) and not in_frame(2100, 2299, minus)


def test_translate11_text_is_public() -> None:
    from gippa.mapping.sixframe import translate11_text

    nt = "ATGGTGATATAANNNAC"  # trailing partial codon dropped
    assert translate11_text(nt, 0) == "MVI*X"  # codon by codon, no I -> L, non-ACGT -> X
    assert translate11_text(nt, 1) == "W*YXX"
    assert translate11_text(nt, 2) == translate(nt[2:])
    assert translate11_text("TTGCTA", 0) == "LL"
    assert isinstance(translate11_text("", 0), str) and translate11_text("AT", 0) == ""


def test_ambiguity_sentinel_is_a_public_constant() -> None:
    from gippa.mapping import sixframe

    assert isinstance(sixframe.AMBIG, int)
    assert sixframe.decode([sixframe.AMBIG]) == "X"  # type: ignore[arg-type]


def test_makedb_uses_only_public_names_of_other_modules() -> None:
    src = Path(gippa.__file__).parent / "makedb" / "build.py"
    tree = ast.parse(src.read_text())
    private = [
        f"{node.module}.{a.name}"
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("gippa")
        for a in node.names
        if a.name.startswith("_")
    ]
    assert private == []
    assert "_AMBIG" not in {
        t.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        for t in node.targets
        if isinstance(t, ast.Name)
    }, "use the public ambiguity constant of gippa.mapping.sixframe"


# ---------------------------------------------------------------- decision 7: rejected inputs


def _patch_with_internal_stop(path: Path) -> Path:
    # chr1 + [1200, 1500): the in-frame stop of the example-4 ORF at [1290, 1293) is internal
    lines = [
        "##gff-version 3",
        "##sequence-region chr1 1 3300",
        "##sequence-region lin1 1 660",
        "chr1\tgippa\tCDS\t1201\t1500\t.\t+\t0\tID=NOVEL_chr1_1201_p;product=bad",
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_patch_with_internal_stop_is_rejected(
    example_json: Path, example: ExampleGenome, tmp_path: Path
) -> None:
    assert "*" in translate(example.genome["chr1"][1200:1497])  # precondition
    patch = _patch_with_internal_stop(tmp_path / "patch.gff3")
    with pytest.raises(ValueError, match="NOVEL_chr1_1201_p"):
        build_db(example_json, tmp_path / "DB", patch=patch)


def test_cli_patch_with_internal_stop_exits_1(example_json: Path, tmp_path: Path) -> None:
    patch = _patch_with_internal_stop(tmp_path / "patch.gff3")
    argv = ["makedb", "--annotation", str(example_json), "--outdir", str(tmp_path / "DB"),
            "--patch", str(patch), "--no-contaminants"]  # fmt: skip
    assert _main(argv) == ExitCode.INVALID_INPUT


def test_annotation_without_features_is_rejected(tmp_path: Path) -> None:
    empty = tmp_path / "empty.json"
    empty.write_text("{}\n")
    with pytest.raises(ValueError):
        build_db(empty, tmp_path / "DB")
    assert not (tmp_path / "DB" / "strain_db.faa").exists()


def test_cli_annotation_without_features_exits_1(tmp_path: Path) -> None:
    empty = tmp_path / "empty.json"
    empty.write_text("{}\n")
    argv = ["makedb", "--annotation", str(empty), "--outdir", str(tmp_path / "DB"),
            "--no-contaminants"]  # fmt: skip
    assert _main(argv) == ExitCode.INVALID_INPUT


# ---------------------------------------------------------------- decision 8: db_id_map.tsv


def test_id_map_is_written_as_db_id_map_tsv(example_json: Path, tmp_path: Path) -> None:
    api = _api()
    res = api.build_strain_db(
        api.MakedbOptions(annotation=example_json, genome=None, outdir=tmp_path / "DB")
    )
    assert res.id_map == tmp_path / "DB" / "db_id_map.tsv"
    assert res.id_map.is_file()
    assert not (tmp_path / "DB" / "id_map.tsv").exists()


# ---------------------------------------------------------------- decision 9: surviving mutants


def test_nterm_interval_is_upstream_plus_30_codons(example_json: Path, tmp_path: Path) -> None:
    """[nterm_span0] the header/id_map interval covers the upstream part and the junction."""
    out = build_db(example_json, tmp_path / "DB", six_frame_orfs=True)
    rows = {r["locus_tag"]: r for r in id_map(out) if r["class"] == "nterm_extension"}
    assert (rows[G1]["start"], rows[G1]["end"], rows[G1]["strand"]) == ("565", "690", "+")
    assert (rows[G3]["start"], rows[G3]["end"], rows[G3]["strand"]) == ("2211", "2330", "-")


def test_x_in_the_upstream_part_never_reaches_an_nterm_entry(tmp_path: Path) -> None:
    """[upstream_x_continue] the in-frame upstream walk does not run through an ``X`` codon."""
    c = Contig("c1", 660, "circular")
    gene = "ATG" + stretch(40, 3) + "TAA"
    c.put(297, "TAA" + stretch(8, 0) + "NNN" + stretch(9, 2) + gene)  # CDS + [354, 480)
    feats = [feature("cds", "c1", 354, 480, "+", locus="R_00015", aa=translate(gene[:-3]))]
    _, out = _one_contig_db(tmp_path, c, feats, six_frame_orfs=True)
    entries = _nterm(out)
    assert all("X" not in s for _, s in entries.values()), entries
    if "R_00015" in entries:  # decision 4: undefined, but never through the X
        assert int(entries["R_00015"][0].split("_")[-3]) - 1 >= 327


def test_contaminant_terminal_stop_is_stripped(example_json: Path, tmp_path: Path) -> None:
    """[cont_rstrip] no `*` is written, even when the reader leaves one (it drops only one)."""
    lib = tmp_path / "contaminants.fasta"
    lib.write_text(">P00761\nMKKLLSIVAG**\n>P02769\nMKWVTFISLL\n")
    out = build_db(example_json, tmp_path / "DB", contaminants=lib)
    seqs = by_id(out)
    assert seqs["CONT_P00761"] == "MKKLLSIVAG"
    assert seqs["CONT_P02769"] == "MKWVTFISLL"


def test_fallback_translation_has_no_terminal_stop(tmp_path: Path) -> None:
    """[protein_rstrip] a CDS without ``aa`` is translated from the genome, without the ``*``."""
    c = Contig("c1", 660, "circular")
    body = stretch(18, 1)
    c.put(297, "TAAGCA" + "GTG" + body + "TAA")  # CDS + [303, 363)
    feats = [feature("cds", "c1", 303, 363, "+", locus="R_00020", product="no aa")]
    _, out = _one_contig_db(tmp_path, c, feats)
    assert by_id(out)["R_00020"] == "M" + translate(body)


def test_genome_with_a_contig_of_the_wrong_length_is_rejected(
    example_json: Path, example: ExampleGenome, tmp_path: Path
) -> None:
    """[genome_len_check]"""
    chr1, lin1 = example.contigs
    short = Contig("lin1", lin1.length - 3, "linear", "contig", seq=list(lin1.text[:-3]))
    fa = write_fasta(tmp_path / "genome.fasta", [chr1, short])
    with pytest.raises(ValueError, match="lin1"):
        build_db(example_json, tmp_path / "DB", genome=fa)
    argv = ["makedb", "--annotation", str(example_json), "--genome", str(fa),
            "--outdir", str(tmp_path / "DB2"), "--no-contaminants"]  # fmt: skip
    assert _main(argv) == ExitCode.INVALID_INPUT


@pytest.mark.parametrize("n", [660, 661, 662])
def test_orf_after_a_stop_at_the_origin_is_written_once(tmp_path: Path, n: int) -> None:
    """[circ_dup] a stop codon at positions 0-2 of a circular contig opens one ORF, not two."""
    c = Contig("c1", n, "circular")
    c.put(0, "TAA" + stretch(25, 0) + "TAA")
    _, out = _one_contig_db(tmp_path, c, [], six_frame_orfs=True)
    ids = orf_ids(out)
    assert len(ids) == len(set(ids))
    assert ids.count("GIPPAORF_intergenic_c1_4_78_p") == 1
    rows = [r["db_id"] for r in id_map(out)]
    assert len(rows) == len(set(rows))


def test_orf_in_frame_with_an_sorf_across_the_origin_is_dropped(tmp_path: Path) -> None:
    """[inframe_noshift] the ORF [630, 732) unwrapped holds the sORF [12, 72) one turn later."""
    c = Contig("c1", 660, "circular")
    sorf = "ATG" + stretch(18, 1) + "TAA"
    c.put(627, "TAA" + stretch(14, 0) + sorf)  # wraps: upstream [630, 672) = [630, 660) + [0, 12)
    feats = [feature("sorf", "c1", 12, 72, "+", locus="R_00025", aa=translate(sorf[:-3]))]
    _, out = _one_contig_db(tmp_path, c, feats, six_frame_orfs=True)
    others = [i for i in orf_ids(out) if "_nterm_extension_" not in i]
    assert others == []
    rid, seq = _nterm(out)["R_00025"]
    assert rid == "GIPPAORF_nterm_extension_c1_631_729_p"
    assert seq == translate(stretch(14, 0)) + translate(sorf[:-3])[:19]


# n = 100 (n % 3 == 1): reading upstream of the CDS + [10, 40) in its frame meets no stop for a
# whole turn (the frame changes after the wrap); the first stop lies 49 codons upstream.
_NO_STOP_IN_ONE_TURN = (
    "TCCGGGGGGAATGCGATACCATATGTACTCTACAAATTAACGGGGACCGCACCTATCATCTAACTGCATCCGCCATGTCCACGG"
    "ACCCACGTGGACGGCA"
)


def test_no_nterm_entry_when_no_stop_within_one_circular_turn(tmp_path: Path) -> None:
    """[upstream_circ_term] the upstream part never spans more than the replicon."""
    seq = _NO_STOP_IN_ONE_TURN
    c = Contig("c1", len(seq), "circular", seq=list(seq))
    feats = [feature("cds", "c1", 10, 40, "+", locus="R_00030", aa=translate(seq[10:37]))]
    _, out = _one_contig_db(tmp_path, c, feats, six_frame_orfs=True)
    assert "R_00030" not in _nterm(out)
