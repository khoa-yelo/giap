"""Recoded Sec codons on the Nissle recording (M7-16c; spec example 29 = example 10's 16c column,
example 2's and example 30's search side, rule 26; test plan C13; ``slow``).

- ``annotated_stops(gm, genome, recoded=True)``: 4,689 terminal + 3 recoded sites, ``not_a_stop``
  3, no skips; the default call stays 4,689 (example 10 unchanged);
- the ``sixframe_index`` step: event ``mapping.stop_sites`` = terminal 4,689, recoded 3,
  not_a_stop 3, cells_target 4,692, 4,692 cells per decoy type, ``recoded_skipped`` zeros;
  the recoded sites are JMJJOG_02511 ``-`` 2578627, JMJJOG_04733 ``+`` 4825395 and JMJJOG_04490
  ``+`` 4576497 (residues 195, 195, 139), all TGA;
- ``map_tolerant``: example 2's Sec peptide and example 30's ``W`` peptide each get a cost-1
  target row over ``[2578585, 2578669) -`` (search side; the placement FDR of a two-peptide run
  is not the point here).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.readthrough_toys import write_annotation, write_genome

pytestmark = pytest.mark.slow

REPO = Path(__file__).resolve().parents[2]
NISSLE_JSON = REPO / "tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz"
NISSLE_GENOME = REPO / "tests/data/genomes/nissle1917.fna.gz"
EX2 = "SLGMLAVDNQARVUHGPTVASLAPTFGR"
EX30 = "SLGMLAVDNQARVWHGPTVASLAPTFGR"
RECODED = {
    ("JMJJOG_02511", "-", 2578627, 195),
    ("JMJJOG_04733", "+", 4825395, 195),
    ("JMJJOG_04490", "+", 4576497, 139),
}


@pytest.fixture(scope="module")
def nissle() -> tuple[dict[str, str], Any]:
    from giap.annotate.bakta import parse_bakta_json
    from giap.io.fasta import read_genome_fasta

    genome = {r.id: r.sequence for r in read_genome_fasta(NISSLE_GENOME)}
    return genome, parse_bakta_json(NISSLE_JSON, genome=genome).model


@pytest.fixture(scope="module")
def steps(tmp_path_factory: pytest.TempPathFactory, nissle: Any) -> dict[str, Any]:
    from giap.config import resolve_config
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    genome, gm = nissle
    tmp = tmp_path_factory.mktemp("nissle16c")
    ann = write_annotation(gm, tmp / "ann")
    peps = pd.DataFrame(
        {
            "peptide": [EX2, EX30, "SAMPLEWLDK"],
            "source": ["database", "database", "denovo"],
            "q_value": [0.001, 0.001, np.nan],
            "pep": [0.01, 0.01, np.nan],
            "score": [np.nan, np.nan, 0.9],
            "is_decoy": [False] * 3,
            "search_db": ["same"] * 3,
        }
    )
    peps.to_parquet(tmp / "peptides.parquet", index=False)
    files: dict[str, Path] = {
        "genome": write_genome(genome, tmp / "genome.fa"),
        "replicons": ann["replicons"],
        "genes": ann["genes"],
        "peptides": tmp / "peptides.parquet",
    }
    cfg = resolve_config(overrides=["seed=1"])
    events: list[logging.LogRecord] = []

    class _H(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            if getattr(record, "event", None) == "mapping.stop_sites":
                events.append(record)

    logger = logging.getLogger("giap.steps.sixframe_index")
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
    return {"files": files, "events": events}


def _contig_of(gm: Any, tag: str) -> str:
    return str(gm.by_locus(tag).interval.contig)


def test_example_29_pure_sites(nissle: Any) -> None:
    from giap.mapping.readthrough import annotated_stops, recoded_skips

    genome, gm = nissle
    sites, n_not = annotated_stops(gm, genome, recoded=True)
    rec = [s for s in sites if s.kind == "recoded"]
    assert (len(sites), len(rec), n_not) == (4692, 3, 3)
    assert recoded_skips(gm, genome) == {"length_mismatch": 0, "codon_mismatch": 0}
    want = set()
    for tag, strand, cs, i in RECODED:
        g = gm.by_locus(tag)
        assert g.aa[i] == "U"
        s, e = g.interval.start, g.interval.end
        assert (s + 3 * i if strand == "+" else e - 3 * i - 3) == cs
        want.add((_contig_of(gm, tag), strand, cs, (tag,), "TGA"))
    assert {(s.contig, s.strand, s.codon_start, s.locus_tags, s.codon) for s in rec} == want
    assert annotated_stops(gm, genome)[0] == tuple(s for s in sites if s.kind == "terminal")


def test_example_29_event(steps: dict[str, Any]) -> None:
    ev = steps["events"]
    assert len(ev) == 1
    assert ev[0].data == {  # type: ignore[attr-defined]
        "terminal": 4689,
        "recoded": 3,
        "not_a_stop": 3,
        "cells_target": 4692,
        "cells_decoy": {"reverse": 4692, "shuffle": 4692},
        "recoded_skipped": {"length_mismatch": 0, "codon_mismatch": 0},
    }


def test_example_29_index_cells(steps: dict[str, Any]) -> None:
    from giap.mapping.sixframe import load_stop_cells

    cells = load_stop_cells(steps["files"]["sixframe"])
    rec = {(s.locus_tags, s.strand, s.codon_start) for s in cells.sites if s.kind == "recoded"}
    assert rec == {((t,), st, cs) for t, st, cs, _ in RECODED}
    assert sum(len(c) for c in cells.target) == 4692
    for k, d in cells.decoys.items():
        assert [len(c) for c in d] == [len(c) for c in cells.target], k


@pytest.mark.parametrize("pep", [EX2, EX30], ids=["example_2_U", "example_30_W"])
def test_examples_2_and_30_get_a_cost_1_row_over_the_fdng_codon(
    steps: dict[str, Any], nissle: Any, pep: str
) -> None:
    _genome, gm = nissle
    ht = pd.read_parquet(steps["files"]["hits_tolerant"])
    t = ht[(ht["peptide"] == pep) & (ht["search"] == "target")]
    here = t[(t["start"] == 2578585) & (t["end"] == 2578669) & (t["strand"] == "-")]
    assert len(here) == 1, t.to_dict("records")
    r = here.iloc[0]
    assert r["contig"] == _contig_of(gm, "JMJJOG_02511")
    assert (int(r["n_diff"]), r["diff_types"], r["route"]) == (1, "substitution", "tolerant")
