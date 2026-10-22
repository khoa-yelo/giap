"""giapsim writes the ``rbs_novel`` layer (issue M7-08 rule 6; ADR-0024; GIAP_SPEC.md §11.3).

Rule 6: ``precomputed/rbs_novel.tsv`` has one row per start codon (``annotation.start_codons``)
on both strands of every replicon; ``tir_ostir`` is drawn from the "true" distribution
(lognormal(8, 0.8), as ``rbs_row``) at planted novel-ORF starts and planted N-terminal true starts
and from the background (lognormal(6, 0.8)) elsewhere; ``config.yaml`` sets
``annotation.precomputed.rbs_novel``.

Position convention (tests/helpers/novel_steps_fixture.py): ``start`` = 1-based first base of the
codon in transcription direction (on ``-`` its highest coordinate). Codons crossing the origin of
a circular replicon may be listed; every non-crossing start codon must be.

Truth starts (transcription-first base, 0-based): planted novel ORFs (``truth/novel_orfs.tsv``,
1-based closed ``start..end``) at ``start - 1`` on ``+`` and ``end - 1`` on ``-``; planted
N-terminal extensions (``truth/features.tsv`` ``true_start``, 1-based) at ``true_start - 1``.

Calibration bounds (fixed seed, so deterministic): mean log TIR over n background rows has SE
0.8 / sqrt(n) < 0.01 for n ~ 10^4, bound +-0.05; over the ~11 truth starts SE ~ 0.24, bound
+-0.9 (~3.7 SE), and the two means differ by > 1.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import pytest
from helpers import novel_steps_fixture as nf

STARTS = ("ATG", "GTG", "TTG")
COUNTS = (("intergenic", 3), ("cds_outframe", 2), ("antisense", 2))


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    return nf.simulate(tmp_path_factory.mktemp("sim_rbs_novel") / "sim", n_novel_orfs=COUNTS,
                       n_nterm_extensions=4)  # fmt: skip


def _rows(sim: Any) -> list[dict[str, str | None]]:
    return nf.tsv_dicts(nf.rbs_novel_layer(sim))


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _truth(sim: Any) -> set[tuple[str, str, int]]:
    out = set()
    for r in _tsv(Path(sim.novel_orfs_truth)):
        first = int(r["start"]) - 1 if r["strand"] == "+" else int(r["end"]) - 1
        out.add((r["contig"], r["strand"], first))
    for r in _tsv(Path(sim.features_truth)):
        if r["type"] == "cds" and r["true_start"]:
            out.add((r["contig"], r["strand"], int(r["true_start"]) - 1))
    return out


def test_layer_file_and_config_key(sim: Any) -> None:
    path = nf.rbs_novel_layer(sim)
    assert path.resolve() == (Path(sim.outdir) / "precomputed" / "rbs_novel.tsv").resolve()
    assert path.read_text().splitlines()[0].split("\t")[: len(nf.RBS_NOVEL_COLUMNS)] == list(
        nf.RBS_NOVEL_COLUMNS
    )


def test_layer_loads_with_load_precomputed(sim: Any) -> None:
    from giap.annotate.precomputed import load_precomputed

    reps = dict(nf.annotation_of(sim).model.replicons)
    df = load_precomputed("rbs_novel", nf.rbs_novel_layer(sim), reps)
    assert len(df) > 0 and not df["tir_ostir"].isna().any()


def test_one_row_per_start_codon_on_both_strands(sim: Any) -> None:
    genome = nf.genome_of(sim)
    reps = nf.annotation_of(sim).model.replicons
    keys = [(r["contig"], r["strand"], int(r["start"] or 0) - 1) for r in _rows(sim)]
    assert len(keys) == len(set(keys)), "duplicate (contig, strand, start)"
    for r, (c, s, p) in zip(_rows(sim), keys, strict=True):
        circ = reps[c].topology == "circular"
        codon = nf.codon_at(genome[c], s, p, circular=circ)
        assert codon in STARTS and r["start_codon"] == codon, (c, s, p, codon)
    want = set()
    for c, seq in genome.items():
        for i in range(len(seq) - 2):
            if seq[i : i + 3] in STARTS:
                want.add((c, "+", i))
            if nf.codon_at(seq, "-", i + 2, circular=False) in STARTS:
                want.add((c, "-", i + 2))
    assert want <= set(keys), sorted(want - set(keys))[:5]
    assert {c for c, _, _ in keys} == set(genome)


def test_truth_starts_are_in_the_layer(sim: Any) -> None:
    keys = {(r["contig"], r["strand"], int(r["start"] or 0) - 1) for r in _rows(sim)}
    truth = _truth(sim)
    assert len(truth) >= 9
    assert truth <= keys, sorted(truth - keys)


def test_tir_distributions_true_vs_background(sim: Any) -> None:
    truth = _truth(sim)
    t, b = [], []
    for r in _rows(sim):
        key = (r["contig"], r["strand"], int(r["start"] or 0) - 1)
        (t if key in truth else b).append(math.log(float(r["tir_ostir"] or "nan")))
    mt, mb = sum(t) / len(t), sum(b) / len(b)
    assert len(b) > 5_000
    assert abs(mb - 6.0) < 0.05, mb
    assert abs(mt - 8.0) < 0.9, mt
    assert mt - mb > 1.0
    sd = math.sqrt(sum((x - mb) ** 2 for x in b) / (len(b) - 1))
    assert abs(sd - 0.8) < 0.05, sd


def test_same_seed_byte_identical_layer(sim: Any, tmp_path: Path) -> None:
    again = nf.simulate(tmp_path / "sim", n_novel_orfs=COUNTS, n_nterm_extensions=4)
    assert nf.rbs_novel_layer(again).read_bytes() == nf.rbs_novel_layer(sim).read_bytes()
