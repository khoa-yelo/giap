"""Invariants of gippasim on the real K-12 genome for any seed (MVP-04 part 1; spec steps 2-4, 8).

On the test-time 200-kb excerpt (tests/helpers/mvp04_data.py), for random seeds:

* the sample differs from K-12 by single-base substitutions only (SNVs only, Q4), so lengths
  and coordinates agree, and the K-12 input file is never modified;
* every TU of ``truth/genes.tsv`` lies on one strand and every recorded protein gene is in one;
* every novel ORF follows the rule (existing K-12 ORF, 30-150 aa, no origin crossing, not an
  annotated protein, class matches position) and every N-terminal extension is a clean in-frame
  upstream start;
* every ``alt_sites`` row reads ``ref_aa`` on K-12 and ``alt_aa`` on the sample;
* the abundance is centred to ``LOG_ABUNDANCE[0]``.

Slow: each example is a full excerpt simulation.
"""

from __future__ import annotations

import functools
import importlib
import math
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from helpers import mvp04_data as M
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippasim.expression import IslandSpec

pytestmark = pytest.mark.slow

ISLANDS = ((2, 3.0), (18, -3.0), (34, 3.0), (50, -3.0), (66, 3.0), (82, -3.0))


@functools.cache
def _excerpt() -> M.Excerpt:
    return M.make_excerpt(Path(tempfile.mkdtemp(prefix="mvp04_props_")))


@settings(
    max_examples=4, deadline=None, suppress_health_check=[HealthCheck.too_slow], derandomize=True
)
@given(seed=st.integers(min_value=0, max_value=2**31 - 1))
def test_excerpt_invariants_hold_for_any_seed(seed: int) -> None:
    from gippasim.proteomics import LOG_ABUNDANCE

    ex = _excerpt()
    base = importlib.import_module("gippasim.scenarios").k12_evidence()
    expr = replace(
        base.expression, islands=tuple(IslandSpec(M.CONTIG, f, 8, e) for f, e in ISLANDS)
    )
    cfg = replace(base, expression=expr, n_novel_per_class=15)
    before = M.sha256(ex.fasta)
    out = Path(tempfile.mkdtemp(prefix=f"mvp04_props_{seed}_"))
    importlib.import_module("gippasim.real_genome").simulate_on_genome(
        cfg, seed, out, genome_fasta=ex.fasta, bakta_json=ex.bakta_json, tu_set=None,
        check_hashes=False,
    )  # fmt: skip
    assert M.sha256(ex.fasta) == before

    sample = M.read_fasta(out / "truth" / "sample_genome.fna")[M.CONTIG]
    assert len(sample) == ex.length

    genes = ex.genes
    tus = M.tu_of(out)
    assert set(tus) == {g.locus for g in genes}
    strands: dict[str, set[str]] = {}
    for g in genes:
        strands.setdefault(tus[g.locus], set()).add(g.strand)
    assert all(len(s) == 1 for s in strands.values())

    rows = M.read_tsv(out / "truth" / "novel_orfs.tsv")
    bad = M.novel_rule_violations(rows, ex.seq, sample, genes)
    bad += M.nterm_violations(rows, ex.seq, genes)
    assert not bad, bad[:10]

    for r in M.read_tsv(out / "truth" / "alt_sites.tsv"):
        c1 = int(r["codon_start"])
        assert M.translate(M.segment(ex.seq, c1, c1 + 2, r["strand"])) == r["ref_aa"]
        assert M.translate(M.segment(sample, c1, c1 + 2, r["strand"])) == r["alt_aa"]

    ab = M.abundance_of(out)
    vals = [ab[g.locus] for g in genes if g.locus in ab]
    assert math.isclose(float(np.mean(vals)), LOG_ABUNDANCE[0], abs_tol=0.1)
