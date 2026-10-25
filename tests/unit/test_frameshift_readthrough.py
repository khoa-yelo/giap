"""The frameshift route with marked stops (M7-16a rule 7; ADR-0058 item 3; statistician review
change 9 = symmetry property (e)).

- A split side never consumes a marked cell (the side DP sees code 22, which is >= 19).
- Precedence: a non-frameshift read-through alignment at cost <= the split's cost that overlaps
  the split locus suppresses the split.
- A decoy frameshift search uses **that decoy type's mirrored cells**: ``nonfs_overlap`` then
  suppresses a decoy split exactly as it suppresses the target split of the corresponding query
  (``reverse``: the decoy search of ``reversed(Q)`` is the target search of ``Q``; ``shuffle``,
  here with the target frames as the "shuffled" frames: the decoy search of ``Q`` is the target
  search of ``Q``). A deliberately broken build (decoy cells removed) shows the decoy split.

World ``H`` (homopolymer flank): frame ``+0`` of a linear contig reads a gene ``...CFHPYWCF``,
its annotated TGA stop, then ``AAAA...`` (``K`` in every frame). ``Q = CFHPYWCF E KKKKKK``:

- without marks the only alignment is the split (prefix on the gene in ``+0``, suffix
  ``GAA AAA ...`` = ``EKKKKKK`` one nt later in ``+1``): cost 1, ``frameshift``;
- with the stop marked, ``E`` over the stop is a cost-1 non-frameshift alignment overlapping the
  split, so the split is suppressed and only the tolerant read-through row remains.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from helpers.readthrough_toys import CODON, gene, model

from giap.config.model import TolerantConfig
from giap.mapping.decoys import build_decoys
from giap.mapping.frameshift import match_frameshift
from giap.mapping.sixframe import Frame, SixFrameIndex, build_sixframe, decode
from giap.mapping.tolerant import match_tolerant
from giap.model import Replicon
from giap.rng import for_step

CFG = TolerantConfig()
GENE_AA = "MCFHPYWCFHPYWCF"
Q = GENE_AA[-8:] + "E" + "K" * 6


def world_h():  # type: ignore[no-untyped-def]
    nt = (
        "".join(CODON[a] for a in GENE_AA)
        + "TGA"
        + "A" * 30
        + "".join(CODON[a] for a in "CFHPYW" * 3)
    )
    reps = {"h": Replicon("h", len(nt), "linear")}
    gm = model(reps, [gene("H1", "h", 0, 48, "+")])
    return {"h": nt}, reps, gm


def stop_cells(index: SixFrameIndex, decoys, gm, genome, *, broken: bool = False):  # type: ignore[no-untyped-def]
    import giap.mapping.readthrough as rt

    sites, _ = rt.annotated_stops(gm, genome)
    target = rt.frame_cells(index, sites)
    dec = {k: rt.mirror_cells(k, fr, target) for k, fr in decoys.items()}
    if broken:
        dec = {k: tuple(np.zeros(0, dtype=np.int64) for _ in c) for k, c in dec.items()}
    return rt.StopCells(sites=sites, target=target, decoys=dec)


def _fs(df: pd.DataFrame, pep: str, search: str) -> pd.DataFrame:
    return df[(df["peptide"] == pep) & (df["search"] == search)]


def test_world_h_split_exists_without_marks() -> None:
    genome, reps, _gm = world_h()
    index = build_sixframe(genome, reps)
    f0 = decode(index.frames[0].aa)
    assert f0[:16] == GENE_AA + "*" and f0[16:26] == "K" * 10
    rev = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 1))}
    df = match_frameshift([Q], index, rev, {Q: "reverse"}, {Q}, CFG)
    tgt = _fs(df, Q, "target")
    assert len(tgt) >= 1 and set(tgt["n_diff"]) == {1}
    assert set(tgt["diff_types"]) == {"frameshift"}
    empty = stop_cells(index, rev, model(reps, []), genome)
    df0 = match_frameshift([Q], index, rev, {Q: "reverse"}, {Q}, CFG, cells=empty)
    pd.testing.assert_frame_equal(df0.reset_index(drop=True), df.reset_index(drop=True))


def test_precedence_suppresses_the_split_under_a_read_through_alignment() -> None:
    genome, reps, gm = world_h()
    index = build_sixframe(genome, reps)
    rev = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 1))}
    cells = stop_cells(index, rev, gm, genome)
    fs = match_frameshift([Q], index, rev, {Q: "reverse"}, {Q}, CFG, cells=cells)
    assert len(_fs(fs, Q, "target")) == 0
    tol = match_tolerant([Q], index, rev, {Q: "reverse"}, {Q}, CFG, cells=cells)
    t = _fs(tol, Q, "target")
    assert len(t) == 1
    r = t.iloc[0]
    assert (r["start"], r["end"], r["strand"], r["n_diff"], r["diff_types"]) == (
        21, 66, "+", 1, "substitution",
    )  # fmt: skip


def _side_world(seed: int):  # type: ignore[no-untyped-def]
    """Gene + TGA + random flank; ``Q`` = 6 gene residues, ``W`` over the stop, 2 flank residues
    (frame +0), then 6 residues of frame +1 starting 1 nt after the prefix ends."""
    rng = np.random.default_rng(seed)
    letters = list("ACDEFGHKLMNPQRSTVWY")
    for _ in range(500):
        flank = "".join(CODON[a] for a in rng.choice(letters, size=30))
        nt = "".join(CODON[a] for a in GENE_AA) + "TGA" + flank
        reps = {"s": Replicon("s", len(nt), "linear")}
        index = build_sixframe({"s": nt}, reps)
        f0, f1 = decode(index.frames[0].aa), decode(index.frames[1].aa)
        c = 15  # the stop
        prefix = f0[c - 6 : c] + "W" + f0[c + 1 : c + 3]
        suffix = f1[c + 3 : c + 9]
        other = f0[c + 3 : c + 9]
        if (
            "*" in suffix
            or "X" in suffix
            or sum(a != b for a, b in zip(suffix, other, strict=True)) < 4
        ):
            continue
        gm = model(reps, [gene("S1", "s", 0, 48, "+")])
        return {"s": nt}, reps, gm, index, prefix + suffix
    raise AssertionError("no side world found")


def test_a_split_side_never_consumes_a_marked_cell() -> None:
    """The only split of ``Q`` would put ``W`` over the marked stop on its prefix side (cost 1 +
    shift 1 = 2 <= k): it must not be reported (rule 7, Q2 default)."""
    genome, _reps, gm, index, q = _side_world(16)
    assert len(q) == 15
    rev = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 1))}
    cells = stop_cells(index, rev, gm, genome)
    fs = match_frameshift([q], index, rev, {q: "reverse"}, {q}, CFG, cells=cells)
    assert len(_fs(fs, q, "target")) == 0
    plain = match_frameshift([q], index, rev, {q: "reverse"}, {q}, CFG)
    assert len(_fs(plain, q, "target")) == 0  # no other split (the construction is clean)


def _symmetry(kind: str, broken: bool) -> tuple[int, int, int, int]:
    """(target fs rows of Q, decoy fs rows of Q', target tolerant rows of Q, decoy tolerant rows
    of Q') with Q' = reversed(Q) for ``reverse`` and Q for ``shuffle``."""
    genome, reps, gm = world_h()
    index = build_sixframe(genome, reps)
    frames = (
        build_decoys(index, "reverse", 50, for_step("sixframe_index", 1))
        if kind == "reverse"
        else tuple(Frame(f.contig, f.strand, f.offset, f.aa.copy(), f.ext_nt) for f in index.frames)
    )  # a "shuffle" decoy equal to the target (a copy): its mirrored cells are the target cells
    decoys = {kind: frames}
    cells = stop_cells(index, decoys, gm, genome, broken=broken)
    qd = Q[::-1] if kind == "reverse" else Q
    peps = [Q, qd] if qd != Q else [Q]
    of = dict.fromkeys(peps, kind)
    fs = match_frameshift(peps, index, decoys, of, {Q}, CFG, cells=cells)  # type: ignore[arg-type]
    tol = match_tolerant(peps, index, decoys, of, {Q}, CFG, cells=cells)  # type: ignore[arg-type]
    return (len(_fs(fs, Q, "target")), len(_fs(fs, qd, "decoy")), len(_fs(tol, Q, "target")),
            len(_fs(tol, qd, "decoy")))  # fmt: skip


def test_reverse_decoy_split_is_suppressed_with_the_mirrored_cells() -> None:
    assert _symmetry("reverse", broken=False) == (0, 0, 1, 1)


def test_shuffle_decoy_split_is_suppressed_with_the_mirrored_cells() -> None:
    assert _symmetry("shuffle", broken=False) == (0, 0, 1, 1)


def test_broken_build_without_decoy_cells_shows_the_decoy_split() -> None:
    """Teeth: with the decoy cells removed the decoy search keeps its split (and loses its
    read-through row) while the target does not, so the counts differ."""
    for kind in ("reverse", "shuffle"):
        t_fs, d_fs, t_tol, d_tol = _symmetry(kind, broken=True)
        assert (t_fs, t_tol) == (0, 1), kind
        assert d_fs >= 1 and d_tol == 0, kind
