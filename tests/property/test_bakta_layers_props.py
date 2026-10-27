"""Hypothesis properties for the Bakta layers (M2-04; GIPPA_SPEC.md §7.1.2, §11.6 TE6).

- ``anti_sd`` always returns a 9-nt DNA string; without the default flag it is the tail of a
  copy and contains ``CCTCC``; with it, it is ``ACCTCCTTA``; a unique majority tail is order-free.
- ``nearest_genes`` equals a brute-force gap computation (origin-aware on circular replicons,
  never wrapping on linear ones): distances ascending, ``min(n, #genes on the contig)`` results.
- ``split_xrefs`` reconstructs every cross-reference; ``trna_counts`` totals the tRNA count.
"""

from __future__ import annotations

from collections import Counter

from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.model import Interval, Replicon
from gippa.model.genes import Gene, NcFeature

DNA = st.text(alphabet="ACGT", min_size=0, max_size=30)


@settings(max_examples=200, deadline=None)
@given(st.lists(DNA, max_size=8))
def test_anti_sd_shape(copies: list[str]) -> None:
    from gippa.annotate.bakta_layers import DEFAULT_ANTI_SD, anti_sd

    seq, flags = anti_sd(copies)
    assert len(seq) == 9 and set(seq) <= set("ACGT")
    assert flags in ((), ("anti_sd_default",))
    if flags:
        assert seq == DEFAULT_ANTI_SD
    else:
        assert "CCTCC" in seq
        assert seq in {c[-9:] for c in copies if len(c) >= 9}


@settings(max_examples=100, deadline=None)
@given(
    st.lists(DNA, min_size=1, max_size=6),
    st.text(alphabet="ACGT", min_size=4, max_size=4),
    st.randoms(use_true_random=False),
)
def test_anti_sd_unique_majority_is_order_free(prefixes: list[str], pad: str, rnd: object) -> None:
    from gippa.annotate.bakta_layers import anti_sd

    tail = "AC" + "CCTCC" + pad[:2]
    others = [p + "T" * 9 for p in prefixes[: max(0, len(prefixes) // 2)]]
    copies = [p + tail for p in prefixes] + others
    shuffled = list(copies)
    rnd.shuffle(shuffled)  # type: ignore[attr-defined]
    assert anti_sd(copies) == anti_sd(shuffled) == (tail, ())


def _gap(a: Interval, b: Interval, rep: Replicon) -> int:
    pa = {p % rep.length for p in range(a.start, a.end)}
    pb = {p % rep.length for p in range(b.start, b.end)}
    best = rep.length
    for p in pa:
        for q in pb:
            d = abs(p - q)
            if rep.topology == "circular":
                d = min(d, rep.length - d)
            best = min(best, d)
    return max(0, best - 1)


@st.composite
def _scene(draw: st.DrawFn) -> tuple[Replicon, NcFeature, list[Gene], int]:
    length = draw(st.integers(40, 160))
    topo = draw(st.sampled_from(["circular", "linear"]))
    rep = Replicon("c1", length, topo, "chromosome")

    def iv() -> Interval:
        start = draw(st.integers(0, length - 1))
        max_end = length + (length // 3 if topo == "circular" else 0)
        end = draw(st.integers(start + 1, min(max_end, start + length // 2 + 1)))
        end = min(end, length) if topo == "linear" else end
        return Interval("c1", start, end, draw(st.sampled_from(["+", "-"])))

    feat = NcFeature(id="ori", type="oric", interval=iv())
    k = draw(st.integers(0, 6))
    genes = [Gene(f"G{i}", "cds", iv(), "M", None, None) for i in range(k)]
    genes.append(Gene("OTHER", "cds", Interval("c2", 0, 10, "+"), "M", None, None))
    return rep, feat, genes, draw(st.integers(1, 4))


@settings(max_examples=200, deadline=None)
@given(_scene())
def test_nearest_genes_matches_brute_force(
    scene: tuple[Replicon, NcFeature, list[Gene], int],
) -> None:
    from gippa.annotate.bakta_layers import nearest_genes

    rep, feat, genes, n = scene
    got = nearest_genes(feat, genes, rep, n=n)
    on_contig = {
        g.locus_tag: _gap(feat.interval, g.interval, rep)
        for g in genes
        if g.interval.contig == rep.id
    }
    assert len(got) == min(n, len(on_contig))
    assert len({loc for loc, _ in got}) == len(got)
    for loc, d in got:
        assert d == on_contig[loc], loc
    dists = [d for _, d in got]
    assert dists == sorted(dists)
    assert dists == sorted(on_contig.values())[: len(got)]


@settings(max_examples=100, deadline=None)
@given(
    st.lists(
        st.lists(
            st.tuples(
                st.sampled_from(["COG", "GO", "EC", "KEGG", "UniRef", "RefSeq", "SO"]),
                st.text(alphabet="ABCDEFKPU0123456789._:", min_size=1, max_size=12),
            ),
            max_size=4,
        ),
        max_size=5,
    )
)
def test_split_xrefs_reconstructs(per_gene: list[list[tuple[str, str]]]) -> None:
    from gippa.annotate.bakta_layers import split_xrefs

    genes = [
        Gene(
            f"L{i}",
            "cds",
            Interval("c1", 10 * i, 10 * i + 9, "+"),
            "M",
            None,
            None,
            db_xrefs=tuple(f"{db}:{x}" for db, x in xr),
        )
        for i, xr in enumerate(per_gene)
    ]
    rows = split_xrefs(genes)
    want = [(g.locus_tag, x) for g in genes for x in g.db_xrefs]
    assert len(rows) == len(want)
    for (locus, db, ident), (wl, x) in zip(rows, want, strict=True):
        assert locus == wl
        assert db == x.split(":", 1)[0]
        assert (ident if db == "GO" else f"{db}:{ident}") == x


@settings(max_examples=100, deadline=None)
@given(st.lists(st.sampled_from(["CAG", "GCC", "TTC", "AAA", "CAT"]), max_size=20))
def test_trna_counts_total(anticodons: list[str]) -> None:
    from gippa.annotate.bakta_layers import trna_counts

    feats = [
        NcFeature(
            id=f"t{i}",
            type="trna",
            interval=Interval("c1", 100 * i, 100 * i + 76, "+"),
            anticodon=ac,
            amino_acid="Xaa",
        )
        for i, ac in enumerate(anticodons)
    ]
    counts = trna_counts(feats)
    assert counts == dict(sorted(Counter(anticodons).items()))
    assert list(counts) == sorted(counts)
