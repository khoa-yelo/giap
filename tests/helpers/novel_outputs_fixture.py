"""Hand-built scenarios, scored frames and file checks for the M7-07 novel-ORF output tests.

The genomes reuse tests/helpers/novel_oracle.py (WALL background, stop-free filler) and the
layouts of docs/specs/M7-06-novel-candidates-and-orfs.md examples 1, 2, 4, 6 and 7. The ORFs come
from the merged M7-06 ``reconstruct`` so the outputs are tested on the real candidate records.

Interface reading pinned here (test author; docs/specs/M7-07-novel-orf-outputs.md):

* ``orf_table(orfs, scored, genome=..., replicons=...)``: the spec's signature cannot read
  ``start_codon_best`` or measure an origin-crossing length without the genome and the replicons,
  so both are passed by keyword (the same names as ``protein_records``);
* ``scored`` is keyed by ``orf_id`` and holds the M7-06 ``minimal_call`` columns, optionally
  ``start_best_ostir`` (0-based lowest coordinate of the chosen start codon, an element of
  ``Orf.starts``; NA when OSTIR has no score) and any other table column (``ostir_tir``,
  ``promoter_score``, the ADR-0012 columns, ``rbs_supported``, ``rank``, ``class_fdr``);
* ``write_novel_orfs(path, table)`` writes ``novel_orfs.tsv``; the two GFF3 files are written with
  the M1-05 ``write_gff3`` and the FASTA with the M1-03 ``write_protein_fasta`` from
  ``(record.header, record.sequence)``.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

from giap.config.model import GiapConfig
from giap.io.gff3 import read_gff3, write_gff3
from giap.io.protein_fasta import write_protein_fasta
from giap.io.tsv import read_tsv
from giap.model.coords import Replicon
from giap.novel.candidates import Orf, candidate_rows, minimal_call, reconstruct
from helpers.novel_oracle import Case, Seq, cds, hit_row, mirror_case, revcomp, translate, wrapped

L = 1000
C1 = {"c1": Replicon("c1", L, "circular", "chromosome")}
LIN = {"lin": Replicon("lin", 300, "linear", "plasmid")}
NOVEL_ORFS_COLUMNS: tuple[str, ...] = (
    "orf_id", "class", "contig", "strand", "frame", "start_longest", "start_best_ostir", "end",
    "length_aa_longest", "length_aa_best", "start_codon_best", "n_peptides", "peptides",
    "best_peptide_q", "ostir_tir", "promoter_score", "promoter_supported",
    "terminator_downstream", "terminator_conf", "terminator_distance_nt", "rbs_supported",
    "class_fdr", "called", "rank", "overlapping_gene", "flags",
)  # fmt: skip
PATCH_PRODUCT = "hypothetical protein (GIAP novel ORF)"
FILES = ("novel_orfs.tsv", "novel_orfs.gff3", "novel_orfs.faa", "annotation_patch.gff3")

EX1_CODONS = ["GCT", "ATT", "GAA", "AAA", "CGT", "GAT", "TCT", "CCG", "CAT", "TGG", "TTC", "AAC",
              "CAG"]  # fmt: skip
P3_CODONS = ["GAC", "TAC", "AAG", "CGC", "GTT", "GAG", "CAC", "TTT", "ACC", "ATC"]
ZZ_CODONS = EX1_CODONS[4:] + EX1_CODONS[:4]  # other peptides: multi-locus ones are excluded


def outputs() -> Any:
    import giap.novel.outputs as m

    return m


def novel_io() -> Any:
    import giap.io.novel_tables as m

    return m


def cfg() -> GiapConfig:
    return GiapConfig(seed=1)


# ---------------------------------------------------------------- scenarios


def ex1_seq(*, atg: bool = True, gtg: bool = True, codons: list[str] = EX1_CODONS) -> Seq:
    """M7-06 example 1: stop TAA [30, 33), ATG [42, 45), GTG [60, 63), peptides from 72,
    stop TGA [120, 123)."""
    s = Seq(L).fill(33, 120).put(30, "TAA").put(120, "TGA")
    if atg:
        s.put(42, "ATG")
    if gtg:
        s.put(60, "GTG")
    return s.put(72, "".join(codons))


def ex1_hits(seq: str, contig: str = "c1", cls: str = "intergenic", n: int = L) -> list[dict]:
    return [
        hit_row(translate(seq[72:102]), contig, 72, 102, "+", cls, n),
        hit_row(translate(seq[90:111]), contig, 90, 111, "+", cls, n),
    ]


def ex1(**kw: bool) -> Case:
    seq = ex1_seq(**kw).text()
    return Case({"c1": seq}, dict(C1), ex1_hits(seq))


def ex6_put(s: Seq) -> Seq:
    """M7-06 example 6 on ``s``: gene G + [300, 900), stop TAA [240, 243), TTG [249, 252),
    peptide [255, 285) classed ``nterm_extension`` (example 2's P3 codons, so it never shares a
    peptide with example 1)."""
    s.fill(243, 897).put(240, "TAA").put(249, "TTG").put(300, "ATG")
    return s.put(255, "".join(P3_CODONS)).put(897, "TAA")


def ex6_hit(seq: str) -> dict:
    return hit_row(translate(seq[255:285]), "c1", 255, 285, "+", "nterm_extension", L, gene="G")


def g_gene() -> Any:
    return cds("G", "c1", 300, 900, "+")


def ex6() -> Case:
    seq = ex6_put(Seq(L)).text()
    return Case({"c1": seq}, dict(C1), [ex6_hit(seq)], genes=[g_gene()])


def ex1_ex6() -> Case:
    """Examples 1 and 6 on one circular genome (M7-07 example 7: the makedb round trip)."""
    seq = ex6_put(ex1_seq()).text()
    return Case({"c1": seq}, dict(C1), [*ex1_hits(seq), ex6_hit(seq)], genes=[g_gene()])


def ex2_put(s: Seq) -> Seq:
    """M7-06 example 2's second ORF: frame-2 peptide [151, 181), ATG [139, 142), stop [190, 193)."""
    s.fill(133, 190).put(130, "TAA").put(139, "ATG").put(151, "".join(P3_CODONS))
    return s.put(190, "TAG")


def two_contigs() -> Case:
    """``zz`` (first in FASTA order) with example 1; ``c1`` with examples 1 and 2."""
    zz = ex1_seq(codons=ZZ_CODONS).text()
    c1 = ex2_put(ex1_seq()).text()
    reps = {"zz": Replicon("zz", L, "circular"), "c1": Replicon("c1", L, "circular")}
    hits = [*ex1_hits(zz, "zz"), *ex1_hits(c1), hit_row(translate(c1[151:181]), "c1", 151, 181,
            "+", "intergenic", L)]  # fmt: skip
    return Case({"zz": zz, "c1": c1}, reps, hits)


def plus_minus() -> Case:
    """Example 1 on ``+`` and, on the same contig, a mirrored copy (other peptides) on ``-``:
    ORFs ``NOVEL_c1_43_p`` ([42, 123)) and ``NOVEL_c1_878_m`` ([877, 958), start ATG [955, 958))."""
    plus = ex1_seq()
    other = ex1_seq(codons=ZZ_CODONS).text()
    m = revcomp(other)
    plus.put(860, m[860:])
    seq = plus.text()
    mir = mirror_case(Case({"c1": other}, dict(C1), ex1_hits(other)))
    return Case({"c1": seq}, dict(C1), [*ex1_hits(seq), *mir.hits])


def wrap_case() -> Case:
    """+ ORF across c1's origin: stop TAA [960, 963), ATG [972, 975), peptides [990, 1011) and
    [14, 35) (unwrapped [1014, 1035)), stop TGA at [50, 53) (unwrapped [1050, 1053))."""
    a = ["GCT", "ATT", "GAA", "AAA", "CGT", "GAT", "TCT"]
    b = ["CCG", "CAT", "TGG", "TTC", "AAC", "CAG", "GAC"]
    s = Seq(L).fill(963, 1050).put(960, "TAA").put(972, "ATG")
    seq = s.put(990, "".join(a)).put(1014, "".join(b)).put(1050, "TGA").text()
    hits = [
        hit_row(translate("".join(a)), "c1", 990, 1011, "+", "intergenic", L),
        hit_row(translate("".join(b)), "c1", 14, 35, "+", "intergenic", L),
    ]
    return Case({"c1": seq}, dict(C1), hits)


def open3_case() -> Case:
    """Linear 300 nt: stop TAA [198, 201), ATG [210, 213), peptides [240, 270), [270, 288); no
    stop downstream before the contig end (last whole codon [297, 300))."""
    tail = ["TTC", "AAC", "CAG", "GAC", "TAC", "AAG"]
    s = Seq(300).put(198, "TAA").fill(201, 300).put(210, "ATG")
    seq = s.put(240, "".join(EX1_CODONS[:10])).put(270, "".join(tail)).text()
    hits = [
        hit_row(translate(seq[240:270]), "lin", 240, 270, "+", "intergenic", 300),
        hit_row(translate(seq[270:288]), "lin", 270, 288, "+", "intergenic", 300),
    ]
    return Case({"lin": seq}, dict(LIN), hits)


def rotate_case(case: Case, shift: int) -> Case:
    """Move every base of the circular contigs by ``+shift`` (hits and genes follow)."""
    genome, hits = {}, []
    for c, seq in case.genome.items():
        n = len(seq)
        k = shift % n
        genome[c] = seq[n - k :] + seq[: n - k] if k else seq
    for h in case.hits:
        n = case.replicons[h["contig"]].length
        a = (h["start"] + shift) % n
        hits.append(hit_row(h["peptide"], h["contig"], a, a + h["end"] - h["start"], h["strand"],
                            h["locus_class"], n, gene=h["gene"]))  # fmt: skip
    if case.genes:
        raise ValueError("rotate_case supports gene-free cases only")
    return Case(genome, dict(case.replicons), hits)


def orfs_of(case: Case) -> list[Orf]:
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    return reconstruct(rows, case.gm, case.genome, cfg())


def scored(
    orfs: list[Orf],
    *,
    best: dict[str, int | None] | None = None,
    called: dict[str, bool] | None = None,
    **cols: dict[str, Any],
) -> pd.DataFrame:
    """``minimal_call`` plus ``start_best_ostir`` (0-based codon, NA when absent) and overrides."""
    df = minimal_call(orfs, cfg())
    ids = df["orf_id"].tolist()
    if best is not None:
        df["start_best_ostir"] = pd.array([best.get(i) for i in ids], dtype="Int64")
    if called:
        df["called"] = [bool(called.get(i, c)) for i, c in zip(ids, df["called"], strict=True)]
    for col, values in cols.items():
        df[col] = pd.Series([values.get(i) for i in ids], dtype="object")
    return df


def table(case: Case, orfs: list[Orf], sc: pd.DataFrame | None) -> pd.DataFrame:
    return outputs().orf_table(orfs, sc, genome=case.genome, replicons=case.replicons)


def proteins(case: Case, t: pd.DataFrame) -> dict[str, str]:
    recs = outputs().protein_records(t, case.genome, case.replicons)
    return {r.header: r.sequence for r in recs}


def write_all(outdir: Path, case: Case, t: pd.DataFrame) -> dict[str, Path]:
    """The four M7-07 files under ``outdir`` (§9.1 ``tables/``)."""
    out = outputs()
    outdir.mkdir(parents=True, exist_ok=True)
    reps = list(case.replicons.values())
    paths = {name: outdir / name for name in FILES}
    novel_io().write_novel_orfs(paths["novel_orfs.tsv"], t)
    write_gff3(paths["novel_orfs.gff3"], out.gff_features(t), reps)
    recs = out.protein_records(t, case.genome, case.replicons)
    write_protein_fasta(paths["novel_orfs.faa"], [(r.header, r.sequence) for r in recs])
    write_gff3(paths["annotation_patch.gff3"], out.patch_features(t, case.gm), reps)
    return paths


def tsv_rows(path: Path) -> list[dict[str, str | None]]:
    header, rows = read_tsv(path)
    assert tuple(header) == NOVEL_ORFS_COLUMNS
    return [dict(zip(header, r, strict=True)) for r in rows]


def feature_lines(path: Path) -> list[list[str]]:
    """Feature lines of a GIAP GFF3 (pragmas and ``Is_circular`` region lines dropped)."""
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if cols[2] == "region" and "Is_circular=true" in cols[8]:
            continue
        out.append(cols)
    return out


def attrs_of(col9: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in col9.split(";"))


def attr_keys(col9: str) -> list[str]:
    return [part.split("=", 1)[0] for part in col9.split(";")]


_BAD_ESCAPE = re.compile(r"%(?![0-9A-Fa-f]{2})")


def check_gff3(path: Path, replicons: dict[str, Replicon]) -> list[Any]:
    """A GFF3 validity check (§9.5): header, pragmas, 9 columns, coordinates, CDS phase, unique
    IDs, escaped attributes; then the M1-05 reader parses the file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines and lines[0] == "##gff-version 3"
    regions = [ln for ln in lines if ln.startswith("##sequence-region")]
    assert regions == [f"##sequence-region {r.id} 1 {r.length}" for r in replicons.values()]
    for r in replicons.values():
        region = f"{r.id}\tgiap\tregion\t1\t{r.length}\t.\t.\t.\tID={r.id};Is_circular=true"
        assert (region in lines) == (r.topology == "circular")
    ids: list[str] = []
    for cols in feature_lines(path):
        assert len(cols) == 9
        seqid, _src, ftype, s, e, score, strand, phase, col9 = cols
        rep = replicons[seqid]
        assert 1 <= int(s) <= int(e)
        assert int(s) <= rep.length
        assert int(e) <= (2 * rep.length if rep.topology == "circular" else rep.length)
        assert strand in ("+", "-")
        assert score == "."
        if ftype == "CDS":
            assert phase == "0"
        for part in col9.split(";"):
            key, sep, value = part.partition("=")
            assert sep and key and value
            assert not _BAD_ESCAPE.search(value) and "&" not in value
        ids.append(attrs_of(col9)["ID"])
    assert len(ids) == len(set(ids)), "duplicate GFF3 IDs"
    return read_gff3(path, replicons)


def fasta_records(path: Path) -> list[tuple[str, str]]:
    recs: list[tuple[str, list[str]]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            recs.append((line[1:], []))
        else:
            recs[-1][1].append(line)
    return [(h, "".join(parts)) for h, parts in recs]


def protein_from(seq: str, lo: int, hi: int, strand: str) -> str:
    """Table-11 protein of unwrapped ``[lo, hi)``: initiator read as M, trailing stop dropped."""
    nt = wrapped(seq, lo, hi)
    nt = nt if strand == "+" else revcomp(nt)
    aa = translate(nt)
    return ("M" + aa[1:]).removesuffix("*")


__all__ = [
    "C1",
    "EX1_CODONS",
    "FILES",
    "LIN",
    "NOVEL_ORFS_COLUMNS",
    "PATCH_PRODUCT",
    "Case",
    "L",
    "attr_keys",
    "attrs_of",
    "cfg",
    "check_gff3",
    "ex1",
    "ex1_ex6",
    "ex6",
    "fasta_records",
    "feature_lines",
    "mirror_case",
    "novel_io",
    "open3_case",
    "orfs_of",
    "outputs",
    "plus_minus",
    "protein_from",
    "proteins",
    "rotate_case",
    "scored",
    "table",
    "tsv_rows",
    "two_contigs",
    "wrap_case",
    "write_all",
]  # fmt: skip
