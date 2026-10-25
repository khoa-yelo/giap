"""Property: ``consumed_sites`` agrees with the cell test helper ``holds_cell`` (M7-16b rule 15,
acceptance "On the property strategies it agrees with ``holds_cell`` on every target row. Exact,
decoy and frameshift rows are never flagged").

Worlds and queries are the M7-16a strategies of ``test_readthrough_props`` (random genomes, linear
and circular, both strands, origin-near stops, queries planted over marked and unmarked stops).
For every target row of ``match_tolerant`` with the implementation's cells, ``consumed_sites(row)
!= ()`` iff ``holds_cell(row)``, and the consumed sites are exactly the sites of the cells inside
the row's frame window. Frameshift rows (``match_frameshift``) and decoy rows consume nothing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import test_readthrough_props as rtp
from helpers.readthrough_counts import holds_cell
from hypothesis import assume, given

from giap.mapping.frameshift import match_frameshift
from giap.mapping.tolerant import match_tolerant


def _consumed(rows: pd.DataFrame, w: rtp.World, cells: object) -> list[tuple[object, ...]]:
    from giap.mapping.readthrough import consumed_sites

    n = w.reps["g"].length
    return consumed_sites(
        rows,
        cells.sites,  # type: ignore[attr-defined]
        {"g": n},
        {"g": w.reps["g"].topology == "circular"},
    )


def _cell_sites(w: rtp.World, r: object, cells: object) -> set[tuple[str, str, int]]:
    """Site keys of the target cells inside the row's frame window."""
    fid, a, b = int(r.frame_id), int(r.aa_start), int(r.aa_end)  # type: ignore[attr-defined]
    n = w.index.replicons[w.index.frames[fid].contig].length
    out = set()
    for p in cells.target[fid]:  # type: ignore[attr-defined]
        if a <= int(p) < b:
            iv = w.index.to_interval(fid, int(p), int(p) + 1)
            out.add((iv.contig, iv.strand, iv.start % n))
    return out


@rtp.SETTINGS
@given(plan=rtp.plans)
@rtp._with_examples
def test_consumed_sites_equal_holds_cell_on_tolerant_rows(plan: rtp.Plan) -> None:
    w = rtp.build_world(plan)
    pep = w.query
    assume(10 <= len(pep) <= 40)
    cells = rtp.impl_cells(w)
    kind = plan.kind
    df = match_tolerant([pep], w.index, w.decoys, {pep: kind}, {pep}, rtp.CFG, cells=cells)
    got = _consumed(df, w, cells)
    assert len(got) == len(df)
    held = holds_cell(df, cells)
    for r, sites, h in zip(df.itertuples(index=False), got, held, strict=True):
        if r.search != "target":
            assert sites == (), ("a decoy row consumed a site", r)
            continue
        assert (sites != ()) == bool(h), (r, sites, h)
        keys = [(s.contig, s.strand, s.codon_start) for s in sites]
        assert keys == sorted(keys, key=lambda k: k[2])
        assert set(keys) == _cell_sites(w, r, cells), (r, sites)


@rtp.SETTINGS
@given(plan=rtp.plans)
@rtp._with_examples
def test_frameshift_rows_never_consume(plan: rtp.Plan) -> None:
    w = rtp.build_world(plan)
    pep = w.query
    assume(10 <= len(pep) <= 40)
    cells = rtp.impl_cells(w)
    kind = plan.kind
    df = match_frameshift([pep], w.index, w.decoys, {pep: kind}, {pep}, rtp.CFG, cells=cells)
    got = _consumed(df, w, cells)
    assert len(got) == len(df)
    fs = np.array([isinstance(b, str) and bool(b) for b in df["blocks"]], dtype=bool)
    for sites, is_fs in zip(got, fs, strict=True):
        if is_fs:
            assert sites == ()
