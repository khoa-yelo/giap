"""gippasim read-through: the #350 review round (S1-S4, S6) for M7-16 rule 14.

The builds in ``test_sim_readthrough.py`` don't reach some boundaries and branches, so several
mutants of ``gippasim.readthrough`` / ``gippasim.proteomics`` survive there. This module covers
them:

* **S1.** ``eligible_readthrough`` on handmade sequences: the next stop at k = 4, 5, 60 and 61;
  a minus-strand scan that runs below 0 on a linear contig (rejected) and wraps on a circular
  one (accepted, ``ext_stop`` wrapped). The plus-strand wrap is covered too.
* **S2.** ``expressed_proteins(readthrough=...)`` against the oracle's ``expressed()``,
  C-terminus included, and ``len(aa) == (end - start) / 3``.
* **S3.** The property test's small config at seed 239 has a minus-strand extension across the
  origin of ``plasmid_1`` (the ``s0 < 0`` branch of ``expressed_proteins``). Truth peptides are
  range-checked: ``1 <= start <= n`` and ``end - start + 1 == 3 * len(peptide)`` (M3-19:
  1-based closed; ``end`` may pass ``n`` across the origin).
* **S4.** ``_holds`` and ``readthrough_rows`` directly, on both strands, including a peptide
  that holds a wrapped stop at ``codon_start + n``.
* **S6.** The shortfall ``ValueError`` must be gippasim's (``match="eligible"``). numpy's own
  ``rng.choice`` error must not satisfy the test.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

import pytest
from helpers.readthrough_sim_oracle import (
    RESIDUE,
    eligible,
    expressed,
    read_genome,
    readthrough_table,
    tsv,
)

import gippasim.proteomics as proteomics
from gippasim.codons import revcomp
from gippasim.features import Feat
from gippasim.proteomics import Protein, ProteomicsConfig
from gippasim.readthrough import (
    ReadThrough,
    _holds,
    eligible_readthrough,
    pick_readthrough,
    readthrough_rows,
)
from gippasim.simulate import SimConfig, simulate

# ---------------------------------------------------------------- handmade sequences (S1)

LEAD = "GCT" * 2  # 6 nt before the gene
GENE = "ATG" + "GCT" * 3  # the gene's sense codons; its stop follows
G0, G1 = len(LEAD), len(LEAD) + len(GENE) + 3  # gene [6, 21), stop at [18, 21)
CS = G1 - 3
TAIL = "GCT" * 5


def _plus(k: int, stop: str = "TAA", new: str = "TGA") -> str:
    """A + strand gene whose next in-frame stop is ``k`` codons after its own (``new``)."""
    return LEAD + GENE + stop + "GCT" * (k - 1) + new + TAIL


def _feat(contig: str, start: int, end: int, strand: str, locus: str = "G1") -> Feat:
    return Feat("cds", contig=contig, start=start, end=end, strand=strand, locus=locus)


def _minus(seq: str, start: int, end: int) -> tuple[str, int, int]:
    """The mirror image: the reverse complement and the mirrored feature interval."""
    n = len(seq)
    return revcomp(seq), n - end, n - start


def _elig(seq: str, f: Feat, circular: bool) -> list[ReadThrough]:
    return eligible_readthrough([f], {f.contig: seq}, {f.contig: circular})


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize(("k", "ok"), [(4, False), (5, True), (60, True), (61, False)])
def test_k_window_is_5_to_60_inclusive(strand: str, circular: bool, k: int, ok: bool) -> None:
    seq = _plus(k)
    n = len(seq)
    if strand == "+":
        f = _feat("c", G0, G1, "+")
        cs, ext = CS, CS + 3 * k
    else:
        seq, s, e = _minus(seq, G0, G1)
        f = _feat("c", s, e, "-")
        cs, ext = n - G1, n - (CS + 3 * k) - 3
    got = _elig(seq, f, circular)
    if not ok:
        assert got == [], (k, strand)
        return
    assert got == [ReadThrough("G1", "c", strand, cs, "TAA", k, ext, "A" * (k - 1))]
    assert got[0].residue == "Q"


def test_k_61_means_no_stop_within_60() -> None:
    """At k = 61 the scan has read 60 sense codons and stops: no gene, not a k = 61 gene."""
    seq = _plus(61)
    assert _elig(seq, _feat("c", G0, G1, "+"), True) == []


#: rotation that puts the origin after the second extension codon
ROT = G1 + 6


def _wrapped_plus(k: int) -> tuple[str, Feat, int, int]:
    """``_plus(k)`` rotated so the extension crosses the origin; (seq, feat, cs, ext_stop)."""
    seq = _plus(k, stop="TGA", new="TAG")
    n = len(seq)
    rot = seq[ROT:] + seq[:ROT]
    f = _feat("c", G0 - ROT + n, G1 - ROT + n, "+")
    return rot, f, CS - ROT + n, (CS + 3 * k - ROT) % n


@pytest.mark.parametrize("k", [5, 24, 60])
def test_plus_strand_extension_wraps_on_a_circular_contig(k: int) -> None:
    seq, f, cs, ext = _wrapped_plus(k)
    assert _elig(seq, f, False) == [], "linear: the scan leaves the contig at its end"
    got = _elig(seq, f, True)
    assert got == [ReadThrough("G1", "c", "+", cs, "TGA", k, ext, "A" * (k - 1))]
    assert ext < cs and got[0].residue == "W"


@pytest.mark.parametrize("k", [5, 24, 60])
def test_minus_strand_scan_below_zero(k: int) -> None:
    """The minus-strand gene sits at [6, 21): its stop at [6, 9), extension codons at 3, 0,
    then -3. Linear: rejected (off the contig). Circular: wraps to n - 3 and on."""
    rot, f, _, _ = _wrapped_plus(k)
    n = len(rot)
    seq, s, e = _minus(rot, f.start, f.end)
    assert (s, e) == (6, 21)
    m = _feat("c", s, e, "-")
    assert _elig(seq, m, False) == [], "linear: the minus-strand scan runs below 0"
    got = _elig(seq, m, True)
    ext = (6 - 3 * k) % n
    assert got == [ReadThrough("G1", "c", "-", 6, "TGA", k, ext, "A" * (k - 1))]
    assert ext > 6, "ext_stop is wrapped into [0, n)"


def test_minus_strand_linear_scan_stopping_exactly_at_zero_is_eligible() -> None:
    """The new stop's lowest base is 0: inside the linear contig, so eligible."""
    k = 5
    plus = _plus(k)[: CS + 3 * k + 3]  # the contig ends right after the new stop
    seq, s, e = _minus(plus, G0, G1)
    got = _elig(seq, _feat("c", s, e, "-"), False)
    assert [(r.k, r.ext_stop) for r in got] == [(k, 0)]


# ---------------------------------------------------------------- S6: our error, not numpy's


def test_shortfall_is_gippasims_error() -> None:
    seq = _plus(5)
    f = _feat("c", G0, G1, "+")
    with pytest.raises(ValueError, match="eligible"):
        pick_readthrough(2, 0, [f], {"c": seq}, {"c": False})
    with pytest.raises(ValueError, match="eligible"):
        pick_readthrough(1, 0, [], {}, {})
    assert len(pick_readthrough(1, 0, [f], {"c": seq}, {"c": False})) == 1


def test_shortfall_in_a_build_is_gippasims_error(tmp_path: Path) -> None:
    cfg = SimConfig(genome_length=12_000, proteomics=ProteomicsConfig())
    n = len(eligible(simulate(cfg, 3, tmp_path / "n0").outdir))
    with pytest.raises(ValueError, match="eligible"):
        simulate(dataclasses.replace(cfg, n_readthrough=n + 1), 3, tmp_path / "nx")


# ---------------------------------------------------------------- S4: _holds, readthrough_rows


@pytest.mark.parametrize(
    ("s0", "e0", "cs", "n", "circ", "want"),
    [
        (0, 30, 9, 100, False, True),  # in frame, inside
        (9, 12, 9, 100, False, True),  # exactly the codon
        (10, 30, 9, 100, False, False),  # starts after the codon
        (0, 11, 9, 100, False, False),  # ends inside the codon
        (1, 31, 11, 100, False, False),  # inside but out of frame
        (0, 5, 3, 100, False, False),  # x + 3 > e0 although x < e0
        (94, 106, 0, 100, True, True),  # the wrapped stop held at cs + n
        (94, 106, 0, 100, False, False),  # linear: no cs + n
        (95, 107, 0, 100, True, False),  # cs + n out of frame
        (97, 103, 0, 100, True, True),  # the codon at [100, 103) only
        (97, 102, 0, 100, True, False),  # cs + n not fully inside
    ],
)
def test_holds(s0: int, e0: int, cs: int, n: int, circ: bool, want: bool) -> None:
    assert _holds(s0, e0, cs, n, circ) is want


def _rec(
    pep: str, locus: str, s1: int, e1: int, strand: str, wrong: bool = False
) -> dict[str, str]:
    return {"peptide": pep, "is_wrong": "true" if wrong else "false", "locus_tag": locus,
            "start": str(s1), "end": str(e1), "strand": strand}  # fmt: skip


def test_readthrough_rows_count_peptides_over_a_wrapped_stop_on_both_strands() -> None:
    """Contig ``c`` (circular, 300 nt). ``P``: + strand, stop at [0, 3) past the origin; ``M``:
    - strand, stop at [3, 6), extension across the origin. Truth rows are 1-based closed with
    ``start <= n`` and ``end`` past ``n`` across the origin."""
    n = 300
    plus = ReadThrough("P", "c", "+", 0, "TGA", 8, 24, "A" * 7)
    minus = ReadThrough("M", "c", "-", 3, "TAA", 8, (3 - 24) % n, "A" * 7)
    truth = [
        _rec("AAAAWAAAA", "P", 289, 315, "+"),  # [288, 315) holds 300 = cs + n
        _rec("AAAAWAAAA", "P", 289, 315, "+"),  # the same peptide again: counted once
        _rec("AAWAAAAAA", "P", 295, 321, "+"),  # [294, 321) holds 300
        _rec("AAAAAAAAA", "P", 304, 330, "+"),  # [303, 330): after the stop
        _rec("AAAWAAAAA", "P", 292, 318, "+", wrong=True),  # wrong: not counted
        _rec("AAAQAAAAA", "M", 292, 318, "-"),  # [291, 318) holds 303 = cs + n
        _rec("AAQAAAAAA", "M", 1, 27, "-"),  # [0, 27) holds 3
        _rec("AAAAAAAAA", "M", 7, 33, "-"),  # [6, 33): the gene body only
        _rec("AAAAWAAAA", "X", 289, 315, "+"),  # another locus
    ]
    rows = readthrough_rows([plus, minus], truth, {"c": n}, {"c": True})
    assert rows == [
        ["P", "c", "+", 0, "TGA", "W", 24, 2],
        ["M", "c", "-", 3, "TAA", "Q", 279, 2],
    ]
    lin = readthrough_rows([plus, minus], truth, {"c": n}, {"c": False})
    assert [r[-1] for r in lin] == [0, 1], "linear: only the unwrapped stop counts"


# ---------------------------------------------------------------- S2/S3: the expressed protein


@dataclasses.dataclass
class Build:
    out: Path
    prots: dict[str, Protein]  # expressed proteins of the n > 0 build
    base: dict[str, Protein]  # the same build's proteins without read-through


def _spy_build(cfg: SimConfig, seed: int, out: Path) -> Build:
    seen: list[tuple[list[Protein], list[Protein]]] = []
    orig = proteomics.expressed_proteins

    def spy(
        feats: Sequence[Feat],
        novel: Sequence[Sequence[Any]],
        seqs: dict[str, str],
        readthrough: Mapping[str, tuple[str, str]] | None = None,
    ) -> list[Protein]:
        got = orig(feats, novel, seqs, readthrough)
        seen.append((got, orig(feats, novel, seqs)))
        return got

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(proteomics, "expressed_proteins", spy)
        outdir = simulate(cfg, seed, out).outdir
    assert len(seen) == 1
    got, base = seen[0]
    return Build(outdir, {p.locus: p for p in got}, {p.locus: p for p in base})


#: the property test's small config; the review found a minus-strand extension across the origin
#: of plasmid_1 at seed 239 (GSIM_00060, k = 24)
SMALL = SimConfig(
    genome_length=12_000,
    n_plasmids=1,
    plasmid_length=(3000, 3000),
    origin_crossing_gene=True,
    n_sorfs=1,
    proteomics=ProteomicsConfig(),
)
MIXED = SimConfig(
    genome_length=30_000,
    n_plasmids=1,
    plasmid_length=(5000, 5000),
    linear_contig_length=4000,
    origin_crossing_gene=True,
    n_sorfs=3,
    n_pseudogenes=2,
    proteomics=ProteomicsConfig(),
)
CASES = {"small239": (SMALL, 239), "mixed39": (MIXED, 39)}


_BUILDS: dict[str, Build] = {}


@pytest.fixture(scope="module")
def builds() -> dict[str, Build]:
    """Lazily built cases, every eligible gene planted so the draw can't skip the
    origin-crossing ones."""
    return _BUILDS


def _get(builds: dict[str, Build], case: str, tmp: pytest.TempPathFactory) -> Build:
    if case not in builds:
        cfg, seed = CASES[case]
        base = tmp.mktemp(f"rt_review_{case}")
        n = len(eligible(simulate(cfg, seed, base / "n0").outdir))
        builds[case] = _spy_build(dataclasses.replace(cfg, n_readthrough=n), seed, base / "n")
    return builds[case]


@pytest.fixture(params=sorted(CASES))
def build(
    request: pytest.FixtureRequest,
    builds: dict[str, Build],
    tmp_path_factory: pytest.TempPathFactory,
) -> Build:
    return _get(builds, request.param, tmp_path_factory)


def _planted(b: Build) -> Iterator[tuple[Any, Any]]:
    elig = eligible(b.out)
    for x in readthrough_table(b.out).itertuples(index=False):
        yield x, elig[x.locus_tag]


def test_expressed_protein_equals_the_oracle(build: Build) -> None:
    """S2: the annotated protein, the stop residue, the k - 1 extension residues, nothing more;
    the interval holds exactly those codons and ends (+) / starts (-) next to the new stop."""
    g = read_genome(build.out)
    planted = list(_planted(build))
    assert planted
    for x, e in planted:
        p, n = build.prots[x.locus_tag], g.n(x.contig)
        assert p.aa == expressed(build.out, e, g), x.locus_tag
        assert p.aa == build.base[x.locus_tag].aa + RESIDUE[e.codon] + p.aa[-(e.k - 1) :]
        assert len(p.aa) == len(build.base[x.locus_tag].aa) + e.k
        assert (p.end - p.start) == 3 * len(p.aa), x.locus_tag
        assert 0 <= p.start < n and p.end <= p.start + n, (x.locus_tag, p.start, p.end)
        if x.strand == "+":
            assert p.end % n == x.ext_stop
        else:
            assert (p.start - 3) % n == x.ext_stop
    for locus, p in build.base.items():
        if locus not in set(readthrough_table(build.out)["locus_tag"]):
            assert build.prots[locus] == p, "other proteins are untouched"


def test_small239_has_a_minus_strand_extension_across_the_origin(
    builds: dict[str, Build], tmp_path_factory: pytest.TempPathFactory
) -> None:
    """S3: the case exists (so the s0 < 0 branch runs) and folds into [0, n)."""
    b = _get(builds, "small239", tmp_path_factory)
    g = read_genome(b.out)
    crossing = [
        x.locus_tag
        for x, e in _planted(b)
        if x.strand == "-" and g.circular[x.contig] and x.codon_start - 3 * (e.k - 1) < 0
    ]
    assert "GSIM_00060" in crossing
    for t in crossing:
        p = b.prots[t]
        n = g.n(p.contig)
        assert 0 <= p.start < n < p.end, (t, p.start, p.end)


def test_truth_peptides_are_in_range(build: Build) -> None:
    """S3: 1-based closed, ``1 <= start <= n``, ``end - start + 1 == 3 * len``, for every correct
    row of a planted gene (and every other correct gene row)."""
    g = read_genome(build.out)
    t = tsv(build.out / "truth" / "peptides.tsv")
    t = t[(t["is_wrong"] == "false") & (t["origin"] == "gene")]
    planted = set(readthrough_table(build.out)["locus_tag"])
    seen = 0
    for r in t.itertuples(index=False):
        n = g.n(r.contig)
        s1, e1 = int(r.start), int(r.end)
        assert 1 <= s1 <= n, (r.locus_tag, r.peptide, s1, e1)
        assert e1 - s1 + 1 == 3 * len(r.peptide), (r.locus_tag, r.peptide, s1, e1)
        assert e1 <= n or g.circular[r.contig], (r.locus_tag, s1, e1)
        seen += r.locus_tag in planted
    assert seen >= 1
