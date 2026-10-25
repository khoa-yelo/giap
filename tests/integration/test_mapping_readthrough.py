"""Annotated-stop read-through through the mapping steps (M7-16a; spec test plan
``test_mapping_readthrough.py``, example 10; ADR-0058 items 1-3).

- giapsim 100 kb (``parse_bakta`` on its Bakta-like JSON; sORFs, pseudogenes and an
  origin-crossing gene switched on): ``sixframe_index`` -> ``map_exact`` -> ``map_tolerant``.
  Index format 2 round trip; the ``mapping.stop_sites`` event equals the counts derived here
  without ``giap.mapping.readthrough``; five read-through peptides (each annotated stop read as
  ``W`` with 7 gene and 7 flank residues, both strands) get a cost-1 target row over their stop,
  classified ``cterm_extension`` of their gene (M3-14, unchanged); exact peptides are unchanged.

  The spec names giapsim ``n_readthrough = 5``; that simulator knob is 16b tooling (rule 14),
  so 16a plants the read-through peptides from the gene model instead.
- Example 10 (``slow``): the committed Nissle recording has 4,689 terminal sites, 0
  ``not_a_stop``, 4,689 target cells and equal decoy cell counts.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from helpers.readthrough_toys import rc

from giap.config import resolve_config
from giap.io.fasta import read_genome_fasta
from giap.mapping.loci import classify_locus
from giap.mapping.sixframe import STOP, build_sixframe, decode
from giap.model.coords import Interval
from giap.steps import load_step_modules
from giap.steps import registry as reg
from giap.steps.mapping import read_gene_model

REPO = Path(__file__).resolve().parents[2]
NISSLE_JSON = REPO / "tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz"
NISSLE_GENOME = REPO / "tests/data/genomes/nissle1917.fna.gz"
STOPS = {"TAA", "TAG", "TGA"}
SEED = 1616


def _terminal(
    seq: str, start: int, end: int, strand: str, n: int, circular: bool
) -> tuple[int, str] | None:
    """(codon_start, gene-strand codon) of a feature's last codon, None off a linear edge."""
    if strand == "+":
        cs = end - 3
        if not circular and (cs < 0 or end > n):
            return None
        codon = "".join(seq[(cs + i) % n] for i in range(3))
        return cs % n, codon
    if not circular and start + 3 > n:
        return None
    return start % n, rc("".join(seq[(start + i) % n] for i in range(3)))


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory):  # type: ignore[no-untyped-def]
    from giapsim.simulate import SimConfig, simulate

    load_step_modules()
    tmp = tmp_path_factory.mktemp("rt")
    res = simulate(
        SimConfig(genome_length=100_000, n_sorfs=4, n_pseudogenes=3, origin_crossing_gene=True),
        SEED, tmp / "sim",
    )  # fmt: skip
    cfg = resolve_config(overrides=[f"seed={SEED}"])
    files: dict[str, Path] = {"genome": Path(res.genome_fasta), "bakta_json": Path(res.bakta_json)}
    spec = reg.get("parse_bakta")
    outs = {k: tmp / f"{k}.out" for k in spec.outputs}
    reg.run_step("parse_bakta", cfg, {k: files[k] for k in spec.inputs}, outs)
    files.update(outs)
    gm = read_gene_model(files["genes"], files["replicons"], files["features"])
    genome = {r.id: r.sequence for r in read_genome_fasta(files["genome"])}
    index = build_sixframe(genome, dict(gm.replicons))
    # the expected sites, derived here
    sites: dict[tuple[str, str, int], list[str]] = {}
    not_a_stop = 0
    for g in gm.genes:
        if g.feature_type not in ("cds", "sorf"):
            continue
        rep = gm.replicons[g.interval.contig]
        t = _terminal(genome[rep.id], g.interval.start, g.interval.end, g.interval.strand,
                      rep.length, rep.topology == "circular")  # fmt: skip
        if t is None or t[1] not in STOPS:
            not_a_stop += 1
            continue
        sites.setdefault((rep.id, g.interval.strand, t[0]), []).append(g.locus_tag)
    cells: dict[int, list[int]] = {}
    for fid, f in enumerate(index.frames):
        n = gm.replicons[f.contig].length
        for p in np.nonzero(f.aa == STOP)[0]:
            iv = index.to_interval(fid, int(p), int(p) + 1)
            if (iv.contig, iv.strand, iv.start % n) in sites:
                cells.setdefault(fid, []).append(int(p))
    # five read-through peptides over clean windows, both strands, and ten exact ones
    cand: dict[str, list[tuple[str, Interval, str]]] = {"+": [], "-": []}
    for fid, ps in sorted(cells.items()):
        f = index.frames[fid]
        text = decode(f.aa)
        n = gm.replicons[f.contig].length
        for p in ps:
            w = text[p - 7 : p + 8] if p >= 7 else ""
            if len(w) != 15 or w.count("*") != 1 or "X" in w:
                continue
            iv = index.to_interval(fid, p - 7, p + 8)
            tags = sites[(iv.contig, iv.strand, index.to_interval(fid, p, p + 1).start % n)]
            if any(getattr(x, "locus_tag", None) not in tags for x in gm.overlapping(iv)):
                continue  # the flank must not run into another annotated feature
            cand[iv.strand].append((w[:7] + "W" + w[8:], iv, tags[0]))
    rt = {pep: (iv, tag) for pep, iv, tag in cand["+"][:3] + cand["-"][:2]}
    exact = []
    for g in gm.genes[:40:4]:
        if g.feature_type == "cds" and len(g.aa) > 40:
            exact.append(g.aa[10:22])
    peps = sorted(rt) + sorted(set(exact) - set(rt))
    n = len(peps)
    table = pd.DataFrame({"peptide": peps, "source": ["database"] * n, "q_value": [0.001] * n,
                          "pep": [0.01] * n, "score": [np.nan] * n, "is_decoy": [False] * n,
                          "search_db": ["same"] * n})  # fmt: skip
    files["peptides"] = tmp / "peptides.parquet"
    table.to_parquet(files["peptides"], index=False)
    logger = logging.getLogger("giap.steps.sixframe_index")
    events: list[logging.LogRecord] = []

    class _H(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            if getattr(record, "event", None) == "mapping.stop_sites":
                events.append(record)

    h = _H(level=logging.INFO)
    logger.addHandler(h)
    old = logger.level
    logger.setLevel(logging.INFO)
    try:
        for name in ("sixframe_index", "map_exact", "map_tolerant"):
            spec = reg.get(name)
            outs = {k: tmp / f"{k}.out" for k in spec.outputs}
            reg.run_step(name, cfg, {k: files[k] for k in spec.inputs}, outs)
            files.update(outs)
    finally:
        logger.removeHandler(h)
        logger.setLevel(old)
    return {"files": files, "gm": gm, "genome": genome, "sites": sites, "not_a_stop": not_a_stop,
            "cells": cells, "rt": rt, "exact": exact, "events": events}  # fmt: skip


def test_the_world_has_both_strands_and_an_origin_copy(chain) -> None:  # type: ignore[no-untyped-def]
    assert len(chain["rt"]) == 5
    assert {iv.strand for iv, _ in chain["rt"].values()} == {"+", "-"}
    assert sum(len(v) for v in chain["cells"].values()) >= len(chain["sites"])


def test_stop_sites_event_counts(chain) -> None:  # type: ignore[no-untyped-def]
    ev = chain["events"]
    assert len(ev) == 1
    n_cells = sum(len(v) for v in chain["cells"].values())
    assert ev[0].data == {  # type: ignore[attr-defined]
        "terminal": len(chain["sites"]), "recoded": 0, "not_a_stop": chain["not_a_stop"],
        "cells_target": n_cells, "cells_decoy": {"reverse": n_cells},
    }  # fmt: skip


def test_index_format_2_round_trip(chain) -> None:  # type: ignore[no-untyped-def]
    from giap.mapping.sixframe import load_stop_cells

    cells = load_stop_cells(chain["files"]["sixframe"])
    assert {(s.contig, s.strand, s.codon_start) for s in cells.sites} == set(chain["sites"])
    for s in cells.sites:
        assert list(s.locus_tags) == sorted(chain["sites"][(s.contig, s.strand, s.codon_start)])
    got = {fid: [int(x) for x in c] for fid, c in enumerate(cells.target) if len(c)}
    assert got == {k: sorted(v) for k, v in chain["cells"].items()}
    assert [len(c) for c in cells.decoys["reverse"]] == [len(c) for c in cells.target]


def test_read_through_peptides_get_a_cost_1_cterm_row(chain) -> None:  # type: ignore[no-untyped-def]
    hits = pd.read_parquet(chain["files"]["hits_tolerant"])
    for pep, (iv, tag) in chain["rt"].items():
        t = hits[(hits["peptide"] == pep) & (hits["search"] == "target")]
        assert len(t) >= 1, pep
        row = t[(t["start"] == iv.start) & (t["end"] == iv.end) & (t["strand"] == iv.strand)]
        assert len(row) == 1, (pep, t)
        r = row.iloc[0]
        assert (r["n_diff"], r["diff_types"], r["route"]) == (1, "substitution", "tolerant")
        cls, g = classify_locus(Interval(r["contig"], int(r["start"]), int(r["end"]), r["strand"]),
                                chain["gm"], chain["genome"])  # fmt: skip
        assert (cls, g) == ("cterm_extension", tag), (pep, cls, g)


def test_exact_peptides_stay_exact(chain) -> None:  # type: ignore[no-untyped-def]
    hx = pd.read_parquet(chain["files"]["hits_exact"])
    found = set(hx.loc[hx["search"] == "target", "peptide"])
    assert set(chain["exact"]) <= found


@pytest.mark.slow
def test_example_10_nissle_sites_and_cells() -> None:
    from giap.annotate.bakta import parse_bakta_json
    from giap.mapping.decoys import build_decoys
    from giap.mapping.readthrough import annotated_stops, frame_cells, mirror_cells
    from giap.rng import for_step

    genome = {r.id: r.sequence for r in read_genome_fasta(NISSLE_GENOME)}
    ann = parse_bakta_json(NISSLE_JSON, genome=genome)
    gm = ann.model
    # all four replicons are circular; 3 partial genes at plasmid edges end without a stop
    assert all(r.topology == "circular" for r in gm.replicons.values())
    sites, n_not = annotated_stops(gm, genome)
    assert (len(sites), n_not) == (4689, 3)
    index = build_sixframe(genome, dict(gm.replicons))
    cells = frame_cells(index, sites)
    assert sum(len(c) for c in cells) == 4689
    for kind in ("reverse", "shuffle"):
        frames = build_decoys(index, kind, 50, for_step("sixframe_index", 1))
        assert [len(c) for c in mirror_cells(kind, frames, cells)] == [len(c) for c in cells]
    # the fdnG Sec codons (internal TGA) are not sites under items 1-5 (example 2)
    keys = {(s.contig, s.strand, s.codon_start) for s in sites}
    assert not any(k[2] in (2578627, 4825395) for k in keys)
