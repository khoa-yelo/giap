"""Performance budget for the M6-04 pair table (docs/specs/M6-04-pair-table.md "Acceptance").

4,500 genes on one circular 4.6 Mb chromosome (strand runs of geometric length, mean 3), all
three evidence layers present (about 10 promoter windows per gene, 3,000 TSS rates, 2,000
terminators) plus 100 RNA features: ``pair_frame`` < 3 s on 1 CPU. Run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
from helpers import pairs_data as pdat

from gippa.annotate.operons import adjacent_pairs
from gippa.annotate.terminators import Terminator
from gippa.model.coords import Interval

N_GENES = 4500
LENGTH = 4_600_000
BUDGET_S = 3.0


def _inputs(rng: np.random.Generator) -> tuple[Any, ...]:
    genes, strand, pos = [], "+", 200
    for i in range(N_GENES):
        if rng.random() < 1 / 3:
            strand = "-" if strand == "+" else "+"
        n = int(rng.integers(300, 1500))
        genes.append(pdat.gene(f"G{i:05d}", pos, pos + n, strand))
        pos += n + int(rng.integers(-20, 300))
    reps = pdat.reps(("c1", max(LENGTH, pos + 500), "circular"))
    probs = {
        (a.locus_tag, b.locus_tag): float(rng.random()) for a, b in adjacent_pairs(genes, reps)
    }
    strands = rng.choice(["+", "-"], 10 * N_GENES)
    windows = pd.DataFrame({
        "contig": "c1", "strand": strands,
        "end3p": rng.integers(0, LENGTH, 10 * N_GENES), "score": rng.random(10 * N_GENES),
    })  # fmt: skip
    tx = pd.DataFrame({
        "contig": "c1", "strand": rng.choice(["+", "-"], 3000),
        "tss": rng.integers(0, LENGTH, 3000), "tx_rate": rng.random(3000) * 10,
    })  # fmt: skip
    terms = []
    for s in rng.integers(0, LENGTH - 50, 2000):
        st = "+" if rng.random() < 0.5 else "-"
        terms.append(
            Terminator(
                Interval("c1", int(s), int(s) + 40, st), float(rng.random()), "TerminatorNet"
            )
        )
    rna = pdat.rna_frame(
        [
            (f"r{i}", "trna", int(s), int(s) + 80, "+")
            for i, s in enumerate(rng.integers(0, LENGTH, 100))
        ]
    )
    ab = {g.locus_tag: float(rng.normal(16, 2)) for g in genes if rng.random() < 0.7}
    ev = pdat.evidence(windows, terms, tx, rna)
    return genes, reps, probs, ev, pdat.gene_tbl(genes, ab)


def test_pair_frame_4500_genes_budget(benchmark: Any) -> None:
    genes, reps, probs, ev, tbl = _inputs(np.random.default_rng(604))

    def run() -> float:
        secs, out = timed(lambda: pdat.frame(genes, reps, probs, ev=ev, tbl=tbl))
        assert set(pdat.PAIR_COLS) <= set(out.columns)
        assert (out["class"] != "SW").sum() == N_GENES
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < BUDGET_S, f"pair table took {secs:.2f} s"
