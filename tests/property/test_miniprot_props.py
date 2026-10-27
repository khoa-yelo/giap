"""Hypothesis properties of the miniprot route (issue M3-12; GIPPA_SPEC.md §7.2.3, §6.5).

- GFF ``mRNA`` coordinates (1-based closed) and ``Target`` (1-based closed on the protein) are
  converted exactly once to 0-based half-open.
- ``reciprocal_overlap`` is symmetric, in ``[0, 1]``, 1 only for equal spans, and invariant
  under a whole-replicon shift of either interval on a circular replicon.
- ``similarity_class`` partitions identity x coverage as the spec defines it.
- ``assign_alignments``: every assigned row passes all three thresholds; assigned proteins and
  leftovers partition the aligned proteins; the result does not depend on input order.
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.config.model import ProteinsConfig
from gippa.model import Interval, Replicon
from gippa.model.genes import Gene, GeneModel

L = 20_000
REP = Replicon("chr1", L, "circular", "chromosome")
REPS = {"chr1": REP}
CFG = ProteinsConfig()


def _mp():  # type: ignore[no-untyped-def]
    from gippa.mapping import miniprot

    return miniprot


def _pa():  # type: ignore[no-untyped-def]
    from gippa.mapping import protein_assign

    return protein_assign


@st.composite
def _span(draw: st.DrawFn, max_len: int = 3000) -> tuple[int, int]:
    start = draw(st.integers(0, L - 1))
    length = draw(st.integers(1, max_len))
    return start, start + length


@settings(max_examples=150, deadline=None)
@given(
    span=_span(),
    strand=st.sampled_from(["+", "-"]),
    plen=st.integers(41, 2000),
    data=st.data(),
    ident=st.integers(0, 10_000),
    fs=st.integers(0, 3),
    sc=st.integers(0, 3),
)
def test_gff_roundtrip(
    span: tuple[int, int], strand: str, plen: int, data: st.DataObject, ident: int, fs: int, sc: int
) -> None:
    start, end = span
    end = min(end, L)  # miniprot reports on the linear sequence
    if end <= start:
        end = start + 1
    qs1 = data.draw(st.integers(1, plen))
    qe1 = data.draw(st.integers(qs1, plen))
    extra = ""
    if fs:
        extra += f";Frameshift={fs}"
    if sc:
        extra += f";StopCodon={sc}"
    attrs = (
        f"ID=MP000001;Rank=1;Identity={ident / 10_000:.4f};Positive=1.0000{extra};"
        f"Target=Q {qs1} {qe1}"
    )
    line = "\t".join(
        ["chr1", "miniprot", "mRNA", str(start + 1), str(end), "0", strand, ".", attrs]
    )
    with TemporaryDirectory() as d:
        p = Path(d) / "a.gff"
        p.write_text("##gff-version 3\n" + line + "\n")
        (a,) = _mp().parse_miniprot_gff(p, REPS, {"Q": plen})
    assert a.interval == Interval("chr1", start, end, strand)  # type: ignore[arg-type]
    assert (a.q_start, a.q_end) == (qs1 - 1, qe1)
    assert abs(a.identity - ident / 10_000) < 1e-9
    assert (a.frameshift, a.stop_codon) == (fs, sc)
    assert abs(a.query_coverage - (qe1 - qs1 + 1) / plen) < 1e-12


@settings(max_examples=300, deadline=None)
@given(a=_span(), b=_span())
def test_reciprocal_overlap_symmetric_bounded(a: tuple[int, int], b: tuple[int, int]) -> None:
    ia = Interval("chr1", a[0], a[1], "+")
    ib = Interval("chr1", b[0], b[1], "+")
    r1 = _pa().reciprocal_overlap(ia, ib, REP)
    r2 = _pa().reciprocal_overlap(ib, ia, REP)
    assert abs(r1 - r2) < 1e-12
    assert 0.0 <= r1 <= 1.0
    assert (r1 == 1.0) == (a == b)


@settings(max_examples=200, deadline=None)
@given(a=_span(), b=_span())
def test_reciprocal_overlap_matches_brute_force(a: tuple[int, int], b: tuple[int, int]) -> None:
    """Origin-aware brute force: positions modulo L on a circular replicon."""
    pa = {p % L for p in range(*a)}
    pb = {p % L for p in range(*b)}
    ov = len(pa & pb)
    expected = min(ov / (a[1] - a[0]), ov / (b[1] - b[0]))
    got = _pa().reciprocal_overlap(Interval("chr1", *a, "+"), Interval("chr1", *b, "+"), REP)
    assert abs(got - expected) < 1e-9


@settings(max_examples=300, deadline=None)
@given(
    identity=st.floats(0.8, 1.0, allow_nan=False),
    coverage=st.floats(0.8, 1.0, allow_nan=False),
    near=st.floats(0.81, 0.99, allow_nan=False),
)
def test_similarity_class_partition(identity: float, coverage: float, near: float) -> None:
    cls = _pa().similarity_class(identity, coverage, near)
    if identity == 1.0 and coverage == 1.0:
        assert cls == "identical"
    elif identity >= near:
        assert cls == "near"
    else:
        assert cls == "divergent"


# ---------------------------------------------------------------- assignment invariants

_GENES = (
    Gene("G1", "cds", Interval("chr1", 100, 403, "+"), "M" + "A" * 99, None, "ATG"),
    Gene("G2", "cds", Interval("chr1", 1000, 1603, "-"), "M" + "A" * 199, None, "ATG"),
    Gene("G3", "cds", Interval("chr1", 1300, 1603, "+"), "M" + "A" * 99, None, "ATG"),
    Gene("G4", "sorf", Interval("chr1", 5000, 5141, "+"), "M" + "A" * 45, None, "ATG"),
    Gene("G5", "cds", Interval("chr1", L - 90, L + 210, "-"), "M" + "A" * 98, None, "ATG"),
)
GM = GeneModel(REPS, _GENES, ())


@st.composite
def _alns(draw: st.DrawFn) -> list[object]:
    n = draw(st.integers(0, 8))
    out = []
    for k in range(n):
        pid = draw(st.sampled_from(["P1", "P2", "P3", "P4"]))
        plen = {"P1": 100, "P2": 200, "P3": 46, "P4": 99}[pid]
        if draw(st.booleans()):  # near a gene
            g = draw(st.sampled_from(_GENES))
            s = max(0, g.interval.start + draw(st.integers(-60, 60)))
            e = min(L, g.interval.end + draw(st.integers(-60, 60)))
            if g.interval.end > L:
                s, e = draw(st.sampled_from([(g.interval.start, L), (0, g.interval.end - L)]))
            strand = draw(st.sampled_from([g.interval.strand, "+", "-"]))
        else:
            s = draw(st.integers(0, L - 400))
            e = s + draw(st.integers(30, 400))
            strand = draw(st.sampled_from(["+", "-"]))
        if e <= s:
            e = s + 1
        qs = draw(st.integers(0, plen - 1))
        qe = draw(st.integers(qs + 1, plen))
        out.append(
            _mp().ProteinAlignment(
                protein_id=pid,
                protein_len_aa=plen,
                route="miniprot",
                interval=Interval("chr1", s, e, strand),
                rank=k + 1,
                identity=draw(st.sampled_from([0.5, 0.79, 0.8, 0.9, 0.95, 0.99, 1.0])),
                positive=None,
                q_start=qs,
                q_end=qe,
                frameshift=draw(st.integers(0, 1)),
                stop_codon=draw(st.integers(0, 1)),
            )
        )
    return out


@settings(max_examples=200, deadline=None)
@given(alns=_alns(), perm=st.randoms(use_true_random=False))
def test_assignment_invariants(alns: list[object], perm: object) -> None:
    df, left = _pa().assign_alignments(alns, GM, CFG)
    pids = {a.protein_id for a in alns}  # type: ignore[attr-defined]
    assigned = set(df.loc[df["gene"].notna(), "protein_id"])
    assert assigned.isdisjoint(left)
    assert assigned | left == pids
    by_tag = {g.locus_tag: g for g in _GENES}
    for _, r in df.iterrows():
        assert r["identity"] >= CFG.min_identity
        assert r["query_coverage"] >= CFG.min_coverage
        assert r["class"] in {"identical", "near", "divergent"}
        assert r["route"] == "miniprot"
        if pd.notna(r["gene"]):
            g = by_tag[r["gene"]]
            assert g.feature_type in ("cds", "sorf")
            assert g.interval.strand == r["strand"]
            assert r["reciprocal_overlap"] >= 0.80 - 1e-12  # §7.2.3 fixed threshold
            iv = Interval(r["contig"], int(r["start"]), int(r["end"]), r["strand"])
            ro = _pa().reciprocal_overlap(iv, g.interval, REP)
            assert abs(ro - r["reciprocal_overlap"]) < 1e-9
        flags = set(str(r["flags"]).split(";")) if isinstance(r["flags"], str) else set()
        assert flags <= {"fs", "st"}
    shuffled = list(alns)
    perm.shuffle(shuffled)  # type: ignore[attr-defined]
    df2, left2 = _pa().assign_alignments(shuffled, GM, CFG)
    pd.testing.assert_frame_equal(df.reset_index(drop=True), df2.reset_index(drop=True))
    assert left == left2
