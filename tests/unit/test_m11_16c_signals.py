"""M11-16c: simulator-v2 signals at genetic codes 4 and 25 (docs/specs/
M11-16-translation-tables-4-25.md item 8: ``signals`` takes the table from the config; lead
decision on #474).

An RBS written into the upstream gene's body (``gippasim.signals._body_rbs``, reached through
``plan_tu(rng, cfg, genes, length)``) must not create a stop codon *of the table*: at 4/25 a TGA
it creates is a sense codon (W/G) and stays; TAA/TAG are still replaced; the codons outside the
RBS window are untouched. End to end, every CDS of a v2 genome at 4/25 keeps a translation equal
to its ``nt`` under the table.

Interface assumed: ``plan_tu`` reads the table from ``cfg.genetic_code`` (``SimConfig``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

N_RUNS = 120
UP = ["ATG", *["TGC"] * 40, "TAA"]  # body codons TGC: an RBS at frame offset 2 makes TG|A
DOWN = ["ATG", *["GCC"] * 20, "TAA"]


def _codons(nt: str) -> list[str]:
    return [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]


def _plan(table: int, seed: int) -> tuple[Any, int]:
    """Plan a two-gene TU with no intra-TU gap (the RBS of gene 2 lies in gene 1's body);
    returns (gene 1, the RBS offset in gene 1)."""
    from gippasim.signals import plan_tu
    from gippasim.simulate import SimConfig, _Gene

    a = _Gene(list(UP), tu=1, strand="+", start=0, end=3 * len(UP))
    b = _Gene(list(DOWN), tu=1, strand="+", start=a.end, end=a.end + 3 * len(DOWN))
    cfg = SimConfig(
        genetic_code=table, rbs_strength=1.0, promoter_strength=0.0, terminator_strength=0.0
    )
    plan = plan_tu(np.random.default_rng(seed), cfg, [a, b], b.end)
    sp = dict((g.start, s) for g, s in plan.rbs)[b.start]
    return a, (b.start - sp - 6) - a.start


@pytest.mark.parametrize("table", [4, 11, 25])
def test_body_rbs_keeps_tga_only_where_it_is_sense(table: int) -> None:
    seen_offset2 = 0
    for seed in range(N_RUNS):
        a, rel = _plan(table, seed)
        nt = a.nt
        cod = _codons(nt)
        assert nt[rel : rel + 6] == "AGGAGG"
        assert cod[0] == "ATG" and cod[-1] == "TAA"
        assert not {"TAA", "TAG"} & set(cod[1:-1]), "the RBS write created TAA/TAG"
        if table == 11:
            assert "TGA" not in cod[1:-1]
        c0 = rel - rel % 3
        if rel % 3 == 2:
            seen_offset2 += 1
            if table != 11:  # TG|A is sense at 4/25: kept, not redrawn
                assert nt[c0 : c0 + 3] == "TGA", (seed, rel, nt[c0 : c0 + 3])
        # codons outside the RBS window are unchanged
        hi = -(-(rel + 6) // 3)
        for i, c in enumerate(cod):
            if i < c0 // 3 or i >= hi:
                assert c == UP[i], (seed, i)
    assert seen_offset2 > 0


def _slice(seq: str, s0: int, e0: int, strand: str) -> str:
    from gippasim.codons import revcomp

    n = len(seq)
    nt = seq[s0:e0] if e0 <= n else seq[s0:] + seq[: e0 - n]
    return revcomp(nt) if strand == "-" else nt


@pytest.mark.parametrize("table", [4, 25])
def test_v2_genome_cds_translation_intact(table: int, tmp_path: Path) -> None:
    """With every RBS planted (strength 1), every CDS's ``nt`` is the genome slice, has no
    TAA/TAG inside and its ``aa`` is the translation under the table."""
    from gippasim.codons import code
    from gippasim.simulate import SimConfig, simulate

    res = simulate(
        SimConfig(genetic_code=table, genome_length=60_000, rbs_strength=1.0, tga_sense_share=0.9),
        2,
        tmp_path,
    )
    data = json.loads(res.bakta_json.read_text())
    seqs = {s["id"]: s.get("nt") or s.get("sequence") for s in data["sequences"]}
    t = code(table)
    n_cds = n_tga = 0
    for f in data["features"]:
        if f["type"] not in ("cds", "sorf") or f.get("pseudo"):
            continue
        n_cds += 1
        cid = f.get("sequence") or f.get("contig")
        s0, e0 = f["start"] - 1, f["stop"]
        if e0 <= s0:
            e0 += len(seqs[cid])
        assert f["nt"] == _slice(seqs[cid], s0, e0, f["strand"])
        cod = _codons(f["nt"])
        assert not {"TAA", "TAG"} & set(cod[1:-1])
        assert f["aa"] == "M" + "".join(t[c] for c in cod[1:-1])
        n_tga += cod[1:-1].count("TGA")
    assert n_cds > 20 and n_tga > 0
