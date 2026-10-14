"""Calibration of giapsim v2's statistical knobs (issue M2-05; GIAP_SPEC.md §11.3, §8.4).

Claims checked, pooled over ``N_SEEDS`` replicate genomes (seeds ``SIM100K_SEED + i``):

* signal strength s is the per-element probability of a planted consensus: the fraction of TUs
  with ``promoter = true``, of genes with a true RBS row and the ratio of true terminators to TUs
  are within 4 binomial SDs of s;
* ``layer_false_positive_rate`` r: pooled false rows of the promoter, terminator and RBS layers
  are within 4 Poisson SDs of r times the pooled true rows;
* the prophage block's GC is 0.08 ± 0.02 below the rest of ``chr1`` in every replicate
  (``prophage_gc_shift = -0.08``);
* UniOP-like probabilities separate within-TU from between-TU adjacent pairs (AUROC >= 0.8).

Conventions as in tests/unit/test_sim_v2.py.
"""

from __future__ import annotations

import csv
import dataclasses
import itertools
import math
from pathlib import Path
from typing import Any

import pytest

import giapsim
from giapsim import simulate

N_SEEDS = 6


def _sim100k() -> tuple[Any, int]:
    for mod in (giapsim, giapsim.simulate):
        if hasattr(mod, "SIM100K"):
            return mod.SIM100K, mod.SIM100K_SEED
    raise AssertionError("giapsim exports no SIM100K / SIM100K_SEED")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _runs(tmp: Path, **kw: Any) -> list[Path]:
    cfg, seed = _sim100k()
    cfg = dataclasses.replace(cfg, **kw)
    outs = []
    for i in range(N_SEEDS):
        out = tmp / f"r{i}"
        simulate(cfg, seed + i, out)
        outs.append(out)
    return outs


def _within(k: int, n: int, p: float) -> None:
    sd = math.sqrt(n * p * (1 - p))
    assert abs(k - n * p) <= 4 * sd + 1, (k, n, p)


@pytest.mark.calibration
def test_signal_strength_is_planting_probability(tmp_path: Path) -> None:
    s = 0.5
    outs = _runs(tmp_path, rbs_strength=s, promoter_strength=s, terminator_strength=s)
    n_tu = k_prom = k_term = n_gene = k_rbs = 0
    for out in outs:
        tus = _rows(out / "truth" / "tus.tsv")
        n_tu += len(tus)
        k_prom += sum(r["promoter"] == "true" for r in tus)
        k_term += sum(r["true_element"] == "true" for r in _rows(out / "truth" / "terminators.tsv"))
        n_gene += len(_rows(out / "truth" / "genes.tsv"))
        k_rbs += sum(r["true_element"] == "true" for r in _rows(out / "truth" / "rbs.tsv"))
    assert n_tu > 100 and n_gene > 200
    _within(k_prom, n_tu, s)
    _within(k_term, n_tu, s)
    # genes.tsv excludes sORFs (which may also get an RBS), so allow their small excess
    n_sorf = N_SEEDS * _sim100k()[0].n_sorfs
    assert abs(k_rbs - n_gene * s) <= 4 * math.sqrt(n_gene * s * (1 - s)) + n_sorf + 1


@pytest.mark.calibration
def test_false_positive_rate(tmp_path: Path) -> None:
    r = 0.2
    outs = _runs(tmp_path, layer_false_positive_rate=r)
    n_true = n_false = 0
    for out in outs:
        for key in ("promoters", "terminators", "rbs"):
            rows = _rows(out / "truth" / f"{key}.tsv")
            n_true += sum(x["true_element"] == "true" for x in rows)
            n_false += sum(x["true_element"] == "false" for x in rows)
    mu = r * n_true
    assert n_true > 300
    assert abs(n_false - mu) <= 4 * math.sqrt(mu), (n_false, mu)


@pytest.mark.calibration
def test_prophage_gc_shift(tmp_path: Path) -> None:
    diffs = []
    for out in _runs(tmp_path):
        (pro,) = [
            r
            for r in _rows(out / "truth" / "regions.tsv")
            if r["region_type"] == "prophage" and r["true_element"] == "true"
        ]
        lines = (out / "genome.fasta").read_text().split(">")[1].splitlines()
        chr1 = "".join(lines[1:])
        s, e = int(pro["start"]) - 1, int(pro["end"])
        block, rest = chr1[s:e], chr1[:s] + chr1[e:]

        def gc(x: str) -> float:
            return (x.count("G") + x.count("C")) / len(x)

        diffs.append(gc(rest) - gc(block))
    assert all(0.06 <= d <= 0.10 for d in diffs), diffs


@pytest.mark.calibration
def test_operon_probabilities_separate_tus(tmp_path: Path) -> None:
    pos: list[float] = []
    neg: list[float] = []
    for out in _runs(tmp_path):
        within = set()
        for tu in _rows(out / "truth" / "tus.tsv"):
            m = tu["genes"].split(";")
            within |= {frozenset(p) for p in itertools.pairwise(m)}
        for r in _rows(out / "precomputed" / "operon_pairs.tsv"):
            key = frozenset((r["gene_a"], r["gene_b"]))
            (pos if key in within else neg).append(float(r["uniop_p"]))
    assert len(pos) > 50 and len(neg) > 50
    wins = sum((p > q) + 0.5 * (p == q) for p in pos for q in neg)
    auc = wins / (len(pos) * len(neg))
    assert auc >= 0.8, auc
