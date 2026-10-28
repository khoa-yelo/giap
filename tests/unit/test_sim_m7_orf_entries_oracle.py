"""gippasim ``OrfEntries`` against real ``gippa makedb --six-frame-orfs`` on every pool product
(M7-12 lead decision 28; code review #225 R1).

``gippasim.proteomics.OrfEntries`` restates ``gippa.makedb.build`` to name the makedb entry a
six-frame product lies in (the ``rev_GIPPAORF_…`` decoy ID of decision 21). This test checks the
restatement exhaustively: every product of the six-frame non-coding pool, in all four classes, gets
exactly the ``db_id`` of the unique entry ``m7_round3_oracle.EntryIndex.entry_of`` finds in the
database GIPPA itself builds (None when there is no entry). The config is stressed so the edge
cases of the restatement are exercised: rRNA operons and tRNAs, pseudogenes, sORFs, a linear
contig, a plasmid and an origin-crossing gene.
"""

from __future__ import annotations

import collections
import dataclasses
from pathlib import Path
from typing import Any

import pytest
from helpers import m7_round3_oracle as Q

import gippasim
import gippasim.proteomics as P
import gippasim.scenarios as S

pytestmark = pytest.mark.slow

STRESS = {
    "n_rrna_operons": 2,
    "n_trna": 10,
    "n_pseudogenes": 6,
    "n_sorfs": 5,
    "linear_contig_length": 20000,
    "n_plasmids": 1,
    "origin_crossing_gene": True,
}


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_orf_entries_match_makedb_on_a_stressed_config(
    seed: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cap: dict[str, Any] = {}
    orig_pool = P._noncoding_pool

    def pool(*a: Any, **k: Any) -> Any:
        cap["pool"] = orig_pool(*a, **k)
        return cap["pool"]

    class Capture(P.OrfEntries):
        def __init__(self, *a: Any, **k: Any) -> None:
            super().__init__(*a, **k)
            cap["oe"] = self

    monkeypatch.setattr(P, "_noncoding_pool", pool)
    monkeypatch.setattr(P, "OrfEntries", Capture)

    cfg = dataclasses.replace(S.novel_orf_power(seed, "100k"), **STRESS)
    out = tmp_path / "sim"
    gippasim.simulate(cfg, seed, out)
    assert "pool" in cap and "oe" in cap, "simulate did not build the pool / OrfEntries"
    entries = Q.makedb_entries(out, tmp_path / "db")
    lens = Q.lengths(out)
    idx = Q.EntryIndex(entries, lens)
    assert len(lens) >= 3, lens  # chromosome, linear contig, plasmid

    pl, oe = cap["pool"], cap["oe"]
    named: collections.Counter[str] = collections.Counter()
    contigs: set[str] = set()
    bad = []
    for cls in P.SIXFRAME_CLASSES:
        for i in range(len(pl.idx[cls])):
            loc = pl.locus(cls, i)
            got = oe.of(loc)
            hits = idx.entry_of(loc[1], loc[2], loc[3], loc[4])
            assert len(hits) <= 1, (cls, loc, [h.db_id for h in hits])
            want = hits[0].db_id if hits else None
            if got != want:
                bad.append((cls, loc, got, want))
            elif got is not None:
                named[cls] += 1
                contigs.add(loc[1])
    assert not bad, f"{len(bad)} mismatches, first: {bad[:8]}"
    for cls in P.SIXFRAME_CLASSES:
        assert named[cls] > 0, (cls, dict(named))
    assert len(contigs) >= 2, contigs
