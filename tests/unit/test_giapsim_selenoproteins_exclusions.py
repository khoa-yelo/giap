"""giapsim ``n_selenoproteins`` exclusion rules (M7-16c rule 23; test plan C11; #346 review S4).

``test_giapsim_selenoproteins.py`` runs BASE at seed 1, which has no novel ORFs and no
N-terminal extensions, and 20 draws rarely reach a margin. These tests pin each rule that
#346's mutation run left unpinned:

- the codon choice ``10 <= i <= len(aa) - 11``, both bounds exactly;
- ``len(aa) >= 21``;
- a codon overlapping any other annotated feature (any type, either strand) is not eligible;
- a codon overlapping a planted novel ORF or the RBS window before it is not eligible;
- a codon overlapping an N-terminal extension or the RBS window before the true start is not
  eligible (mutant M33);
- the expressed protein carries ``U`` at annotated residue ``i`` after the extension;
- the draws come from their own pinned stream (a digest of ``truth/selenoproteins.tsv``).

The toy tests call ``plant_selenoproteins`` with ``n = 1`` on a layout with exactly one
eligible residue, over 40 seeds: a mutant that admits a second residue is drawn at least once
with probability ``1 - 2**-40``.
"""

from __future__ import annotations

import copy
import hashlib
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from giapsim.codons import revcomp, translate_cds
from giapsim.features import Feat
from giapsim.proteomics import ProteomicsConfig, expressed_proteins
from giapsim.selenoproteins import plant_selenoproteins
from giapsim.simulate import SimConfig, simulate

SEEDS = range(40)
N = 400  # toy contig length
OFF = 100  # first base of the toy gene A
L = 40  # residues of gene A (codons 0..39, then the stop): margins allow i in 10..29


def _gene(strand: str = "+", n_aa: int = L, off: int = OFF, locus: str = "A") -> Feat:
    return Feat(
        "cds", "c", off, off + 3 * (n_aa + 1), strand, aa="M" + "A" * (n_aa - 1), locus=locus
    )


def _block(s: int, e: int, strand: str = "-", locus: str = "blk") -> Feat:
    """A non-CDS annotated feature covering ``[s, e)``."""
    return Feat("tRNA", "c", s, e, strand, locus=locus)


def _codon(i: int, strand: str = "+", n_aa: int = L, off: int = OFF) -> tuple[int, int]:
    """``[s, e)`` of codon ``i`` of gene A."""
    cs = off + 3 * i if strand == "+" else off + 3 * (n_aa + 1) - 3 * i - 3
    return cs, cs + 3


def _draws(feats: Sequence[Feat], novel: Sequence[Sequence[Any]] = ()) -> set[tuple[str, int]]:
    """Every (locus, residue) drawn with ``n = 1`` over ``SEEDS``."""
    out: set[tuple[str, int]] = set()
    for seed in SEEDS:
        seqs = {"c": bytearray(b"A" * N)}
        rows = plant_selenoproteins(1, seed, seqs, copy.deepcopy(list(feats)), novel)
        assert len(rows) == 1
        out.add((str(rows[0][0]), int(rows[0][4])))
    return out


# ---------------------------------------------------------------- the toy layouts are sound


def test_without_exclusions_every_margin_residue_is_drawn() -> None:
    """Guard for the tests below: alone, gene A offers residues 10..29 and the draw is uniform
    enough to reach both ends; 2000 draws, each residue has p = 1/20."""
    seen: set[int] = set()
    for seed in range(2000):
        seqs = {"c": bytearray(b"A" * N)}
        seen.add(int(plant_selenoproteins(1, seed, seqs, [_gene()], [])[0][4]))
    assert seen == set(range(10, 30))


# ---------------------------------------------------------------- margins and minimum length


@pytest.mark.parametrize("strand", ["+", "-"])
def test_lower_margin_is_residue_10(strand: str) -> None:
    """Residues 0..10 free, 11..39 covered: only residue 10 (``i >= 10``, not 9 or 11)."""
    s10, e10 = _codon(10, strand)
    s11, e11 = _codon(11, strand)
    s39, e39 = _codon(39, strand)
    cover = (min(s11, s39), max(e11, e39))
    assert not (cover[0] < e10 and s10 < cover[1])
    assert _draws([_gene(strand), _block(*cover)]) == {("A", 10)}


@pytest.mark.parametrize("strand", ["+", "-"])
def test_upper_margin_is_len_minus_11(strand: str) -> None:
    """Residues 0..28 covered, 29..39 free: only residue 29 = len(aa) - 11 (not 30 or 28)."""
    s0, e0 = _codon(0, strand)
    s28, e28 = _codon(28, strand)
    assert _draws([_gene(strand), _block(min(s0, s28), max(e0, e28))]) == {("A", L - 11)}


def test_a_21_residue_gene_is_eligible_and_a_20_residue_one_is_not() -> None:
    """``len(aa) >= 21``: a 21-aa gene has the single residue 10; a 20-aa gene has none.

    (``_MIN_AA`` 21 -> 5 is an equivalent mutant: with the margins, ``range(10, len(aa) - 10)``
    is empty for every ``len(aa) <= 20``. This pins the other direction, 21 -> 22.)"""
    assert _draws([_gene(n_aa=21)]) == {("A", 10)}
    with pytest.raises(ValueError, match="only 0 eligible genes"):
        plant_selenoproteins(1, 0, {"c": bytearray(b"A" * N)}, [_gene(n_aa=20)], [])


# ---------------------------------------------------------------- other annotated features


@pytest.mark.parametrize(("strand", "other"), [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")])
def test_a_codon_overlapping_another_feature_is_not_eligible(strand: str, other: str) -> None:
    """Two features (either strand) cover every margin residue but 20: only 20 is drawn; a codon
    sharing even one base with another feature is out (coverage exactly 1)."""
    s20, e20 = _codon(20, strand)
    # the blockers end and start exactly at codon 20's bounds
    left, right = _block(OFF, s20, other, "L"), _block(e20, OFF + 3 * (L + 1), other, "R")
    assert _draws([_gene(strand), left, right]) == {("A", 20)}
    # shift the right blocker one base into codon 20: no eligible residue is left
    right_in = _block(e20 - 1, OFF + 3 * (L + 1), other, "R")
    with pytest.raises(ValueError, match="only 0 eligible genes"):
        plant_selenoproteins(1, 0, {"c": bytearray(b"A" * N)}, [_gene(strand), left, right_in], [])


def test_an_overlapping_neighbour_cds_is_not_eligible_either() -> None:
    """A neighbour CDS sharing codons 11..39 of A (opposite strand, inside A): A offers only
    residue 10, and the neighbour's own codons are all shared, so it offers none."""
    s11 = _codon(11)[0]
    nb = Feat("cds", "c", s11, OFF + 3 * (L + 1), "-", aa="M" + "A" * 28, locus="B")
    assert (nb.end - nb.start) // 3 - 1 == len(nb.aa) == 29
    assert _draws([_gene(), nb]) == {("A", 10)}


# ---------------------------------------------------------------- planted novel ORFs


def test_a_plus_novel_orf_and_its_rbs_window_are_excluded() -> None:
    """A ``+`` novel ORF on 0-based ``[153, 300)``: its span covers codons 18..29 of A and the
    20-bp window before it, ``[133, 153)``, covers codons 11..17. Only residue 10 is left."""
    assert _codon(10) == (130, 133) and _codon(11)[0] == 133 and _codon(17)[1] > 153 - 3
    novel = [("orf1", "cds_outframe", "c", 154, 300, "+", "M" * 48)]
    assert _draws([_gene()], novel) == {("A", 10)}


def test_a_minus_novel_orf_and_its_rbs_window_are_excluded() -> None:
    """A ``-`` novel ORF on ``[40, 130)``: its RBS window lies above it, ``[130, 150)``, over
    codons 10..16 of A; a tRNA covers codons 17..28. Only residue 29 is left."""
    assert _codon(10)[0] == 130 and _codon(16)[0] < 150 and _codon(29) == (187, 190)
    novel = [("orf1", "intergenic", "c", 41, 130, "-", "M" * 29)]
    assert _draws([_gene(), _block(_codon(17)[0], _codon(28)[1])], novel) == {("A", 29)}


# ---------------------------------------------------------------- N-terminal extensions (M33)


def test_a_plus_n_extension_and_its_rbs_window_are_excluded() -> None:
    """Gene B (``+``) has its true start at 0-based 153 and its annotated start 15 bp later;
    its span (extension included) covers codons 18..39 of A (18..22 by the extension only) and
    the 20-bp RBS window before the true start, ``[133, 153)``, covers codons 11..17. Only
    residue 10 is left."""
    b = Feat("cds", "c", 168, 223, "+", aa="M" * 5, locus="B", true_start=154)
    assert b.span() == (153, 223)
    assert _draws([_gene(), b]) == {("A", 10)}


def test_a_minus_n_extension_and_its_rbs_window_are_excluded() -> None:
    """Gene B (``-``) is annotated on ``[100, 115)`` with its true start at 0-based 129 (span
    ``[100, 130)``, codons 0..9 of A); its RBS window lies above the true start, ``[130,
    150)``, over codons 10..16; a tRNA covers 17..28. Only residue 29 is left."""
    b = Feat("cds", "c", 100, 115, "-", aa="M" * 3, locus="B", true_start=130)
    assert b.span() == (100, 130)
    assert _draws([_gene(), b, _block(_codon(17)[0], _codon(28)[1])]) == {("A", 29)}


# ---------------------------------------------------------------- the expressed protein


@pytest.mark.parametrize("strand", ["+", "-"])
def test_expressed_protein_has_u_after_the_extension(strand: str) -> None:
    """A gene with a 5-codon N-terminal extension and Sec at annotated residue 10: the
    expressed protein (true start) reads U at 5 + 10 = 15 and has no ``*``."""
    ext, body = "ATG" + "GCT" * 4, "ATG" + "GCT" * 9 + "TGA" + "GCT" * 14 + "TAA"
    plus = "CC" + ext + body + "CC"
    n_aa = len(body) // 3 - 1
    assert n_aa == 25
    if strand == "+":
        seq = plus
        f = Feat("cds", "c", 2 + len(ext), 2 + len(ext) + len(body), "+", true_start=3)
    else:
        seq = revcomp(plus)
        f = Feat("cds", "c", 2, 2 + len(body), "-", true_start=2 + len(ext) + len(body))
    f.aa = "M" + "A" * 9 + "U" + "A" * 14
    f.sec_residues = (10,)
    f.locus = "A"
    (p,) = expressed_proteins([f], [], {"c": seq})
    assert p.aa == "M" + "A" * 4 + f.aa
    assert p.aa[15] == "U" and "*" not in p.aa
    assert translate_cds([ext[i : i + 3] for i in range(0, len(ext), 3)] + ["TAA"]) == "MAAAA"


# ---------------------------------------------------------------- a full build with these edits

CFG = SimConfig(
    genome_length=100_000,
    n_nterm_extensions=30,
    n_novel_orfs=(("intergenic", 10), ("cds_outframe", 10), ("antisense", 10)),
    n_selenoproteins=40,
    proteomics=ProteomicsConfig(search_db="same"),
)
SEED = 1


@pytest.fixture(scope="module")
def build(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("sec_x") / "o"
    simulate(CFG, SEED, out)
    return out


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _overlaps(a: int, b: int, s: int, e: int, n: int) -> bool:
    return any(s + d < b and a < e + d for d in (-n, 0, n))


def test_full_build_keeps_sec_codons_off_novel_orfs_and_extensions(build: Path) -> None:
    """With novel ORFs (all classes) and N-terminal extensions planted, no Sec codon overlaps a
    novel ORF or the 20 bp before its start, nor another gene's extension or the 20 bp before
    its true start; at seed 1 some Sec genes carry an extension themselves."""
    sec = _tsv(build / "truth" / "selenoproteins.tsv")
    nov = _tsv(build / "truth" / "novel_orfs.tsv")
    ft = _tsv(build / "truth" / "features.tsv")
    n = 100_000  # the chromosome; BASE-like builds have one contig
    assert set(sec["contig"]) == set(ft["contig"]) == {sec["contig"].iloc[0]}
    windows: list[tuple[str, int, int]] = []
    for r in nov.itertuples(index=False):
        s, e = int(r.start) - 1, int(r.end)
        rbs = (s - 20, s) if r.strand == "+" else (e, e + 20)
        windows += [(str(r.orf_id), s, e), (str(r.orf_id), *rbs)]
    ext = ft[(ft["type"] == "cds") & (ft["true_start"] != "") & (ft["pseudo_reason"] == "")]
    for r in ext.itertuples(index=False):
        s, e, ts = int(r.start) - 1, int(r.end), int(r.true_start)
        span = (ts - 1, e) if r.strand == "+" else (s, ts)
        rbs = (span[0] - 20, span[0]) if r.strand == "+" else (span[1], span[1] + 20)
        windows += [(str(r.feature_id), *span), (str(r.feature_id), *rbs)]
    assert len(nov) == 30 and len(ext) == 30
    for x in sec.itertuples(index=False):
        a = int(x.codon_start)
        for who, s, e in windows:
            if who != x.locus_tag:
                assert not _overlaps(a, a + 3, s, e, n), (x.locus_tag, who, s, e)
    assert sec["locus_tag"].isin(set(ext["feature_id"])).sum() >= 1


def test_full_build_peptides_over_an_extended_sec_gene_carry_u(build: Path) -> None:
    """True peptides of an extended Sec gene that span the codon read U there."""
    sec = _tsv(build / "truth" / "selenoproteins.tsv")
    ft = _tsv(build / "truth" / "features.tsv").set_index("feature_id")
    pep = _tsv(build / "truth" / "peptides.tsv")
    pep = pep[(pep["is_wrong"] == "false") & (pep["origin"] == "gene")]
    spanning = 0
    for x in sec.itertuples(index=False):
        if ft.loc[x.locus_tag, "true_start"] == "":
            continue
        a = int(x.codon_start)
        for r in pep[pep["locus_tag"] == x.locus_tag].itertuples(index=False):
            s0, e = int(r.start) - 1, int(r.end)
            if not (s0 <= a and a + 3 <= e):
                continue
            k = (a - s0) // 3 if x.strand == "+" else (e - a - 3) // 3
            assert r.peptide[k] == "U", (x.locus_tag, r.peptide, k)
            spanning += 1
    assert spanning >= 1, "no true peptide spans a Sec codon of an extended gene"


# ---------------------------------------------------------------- the stream


#: sha256 of truth/selenoproteins.tsv for CFG at seed 1, recorded on 8ad7f316
PINNED_TRUTH = "58eaa53a8eb7a4f1333a59e6d44c2982b3c90e8f0e1f6306098a0d6de1c825cd"


def test_pinned_truth_digest(build: Path) -> None:
    """Pins the draws (stream key ``[seed, 12]``, gene and residue choice) of a build with every
    exclusion source present."""
    digest = hashlib.sha256((build / "truth" / "selenoproteins.tsv").read_bytes()).hexdigest()
    assert digest == PINNED_TRUTH
