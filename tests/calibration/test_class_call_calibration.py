"""Calibration of the M7-11b call path (docs/specs/M7-11b-class-fdr-wiring.md rules 2, 3, 5, 6;
GIPPA_SPEC.md §7.6.4, §8.4 row "Novel-ORF class FDR"; ADR-0015).

The M7-11 suite (tests/calibration/test_novel_class_fdr_calibration.py) calibrates the pure
``picked_class_q`` / ``denovo_class_q`` on a given ``id_map``. This one runs the *wired* path the
step runs: ``class_fdr`` with its own ``makedb_id_map`` (no ``db_id_map.tsv``, rule 2) and
replicon lengths, ``denovo_inputs`` on a ``placement_costs`` frame plus hit tables (rule 3), and
the TSV cell ``class_fdr_cell`` (rule 5, ``lost_pick`` = 1). The statistical claim is rule 6 with
the class rule met: the ORFs whose cell is numeric and <= a have mean FDP <= a + 2 SE at a in
``levels``. R = ``replicates_fast.novel_orfs`` (full: ``replicates_full.novel_orfs``).

Generators: as the M7-11 suite (picked: three makedb classes, a false target and its decoy are
exchangeable best-of-k nulls; de novo: single-peptide ORFs, π̂ = the bin's wrong rate), smaller
(classes 600 / 300 / 300, 2000 de novo peptides) because the wired path also builds the ID map.
"""

from __future__ import annotations

import math
import os
from collections.abc import Callable

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers.class_call_fixture import (
    cc,
    config,
    costs_frame,
    decoy_hits,
    gene_model,
    target_hits,
)
from helpers.novel_fdr_fixture import (
    ContigClassifier,
    PickedWorld,
    make_orf,
    pep_name,
    pep_row,
    peptides_frame,
)

from gippa.model.coords import Interval

pytestmark = pytest.mark.calsuite("novel_orfs")

SUITE = "novel_orfs"
CLASS_SIZES = {"intergenic": 600, "cds_outframe": 300, "antisense": 300}
CLASSES = tuple(CLASS_SIZES)
P_TRUE = 0.4
P_HIT = 0.5
N_DENOVO = 2000
P_CHANCE = 0.25
CHANCE_W = np.array([0.6, 0.25, 0.15])
CORRECT_W = np.array([0.2, 0.4, 0.4])
CONTIG_LEN = 400_000
GENOME = {f"k{i}": "A" * CONTIG_LEN for i in range(len(CLASSES))}
GM = gene_model(GENOME)


def _q_in(score: float) -> float:
    return 0.01 / (1.0 + float(np.exp(score)))


def _best(rng: np.random.Generator, mu: float) -> float:
    return float(rng.normal(mu, 1.0, rng.integers(1, 4)).max())


Replicate = tuple[dict[str, float | str], set[str], set[str]]


def _cells(res: object) -> dict[str, float | str]:
    m = cc()
    return {k: m.class_fdr_cell(r) for k, r in res.orfs.items()}  # type: ignore[attr-defined]


def _picked_replicate(rng: np.random.Generator) -> Replicate:
    w = PickedWorld()
    false_ids: set[str] = set()
    for cls, n in CLASS_SIZES.items():
        contig = f"k{CLASSES.index(cls)}"
        for _ in range(n):
            true = rng.random() < P_TRUE
            has_t = true or rng.random() < P_HIT
            has_d = rng.random() < P_HIT
            if not (has_t or has_d):
                continue
            entry = w.entry(cls, contig, "+", in_map=False)
            if has_t:
                s = _best(rng, 3.0 if true else 0.0)
                orf, _ = w.target(cls, contig, 0, entry=entry, q=_q_in(s), score=s)
                if not true:
                    false_ids.add(orf.orf_id)
            if has_d:
                s = _best(rng, 0.0)
                w.decoy(entry[0], 0, q=_q_in(s), score=s)
    res = cc().class_fdr(w.orfs, w.peptides, gm=GM, genome=GENOME, cfg=config(), denovo=None)
    return _cells(res), {o.orf_id for o in w.orfs}, false_ids


def _denovo_replicate(rng: np.random.Generator) -> Replicate:
    rates = np.linspace(0.02, 0.9, 10)
    pi_true = rates[rng.integers(0, rates.size, N_DENOVO)]
    wrong = rng.random(N_DENOVO) < pi_true
    has_t = ~wrong | (rng.random(N_DENOVO) < P_CHANCE)
    has_d = rng.random(N_DENOVO) < P_CHANCE
    e_t = np.where(wrong, 10 ** rng.normal(-1.5, 1.0, N_DENOVO),
                   10 ** rng.normal(-3.5, 1.0, N_DENOVO))  # fmt: skip
    e_d = 10 ** rng.normal(-1.5, 1.0, N_DENOVO)
    t_cls = np.where(wrong, rng.choice(3, N_DENOVO, p=CHANCE_W),
                     rng.choice(3, N_DENOVO, p=CORRECT_W))  # fmt: skip
    d_cls = rng.choice(3, N_DENOVO, p=CHANCE_W)
    peps = [pep_name(k, "DNVO") for k in range(N_DENOVO)]
    costs = costs_frame(peps, e_target=np.where(has_t, e_t, np.nan),
                        e_decoy=np.where(has_d, e_d, np.nan), pi=pi_true)  # fmt: skip
    ti, di = np.flatnonzero(has_t), np.flatnonzero(has_d)
    hits = pd.concat([
        target_hits([peps[i] for i in ti], [Interval("k0", 100, 130, "-")] * len(ti)),
        decoy_hits([peps[i] for i in di],
                   [Interval(f"k{d_cls[i]}", 100, 130, "-") for i in di]),
    ], ignore_index=True)  # fmt: skip
    clf = ContigClassifier({f"k{i}": c for i, c in enumerate(CLASSES)})
    d = cc().denovo_inputs(costs, hits, clf)
    orfs = [make_orf(f"N_{i}", CLASSES[t_cls[i]], f"k{t_cls[i]}", "-", 97, [peps[i]])
            for i in ti]  # fmt: skip
    c9 = peptides_frame([pep_row(peps[i], [], 0, is_decoy=False, source="denovo", q=math.nan)
                         for i in ti])  # fmt: skip
    res = cc().class_fdr(orfs, c9, gm=GM, genome=GENOME, cfg=config(), denovo=d)
    return _cells(res), {o.orf_id for o in orfs}, {f"N_{i}" for i in ti if wrong[i]}


def _check(name: str, make: Callable[[np.random.Generator], Replicate], full: bool) -> None:
    levels = C.load_tolerances().levels
    n_acc: list[int] = []

    def fdp(rng: np.random.Generator) -> np.ndarray:
        cells, ids, false_ids = make(rng)
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = {o for o in ids if not isinstance(cells[o], str) and cells[o] <= a}
            out[j] = len(called & false_ids) / len(called) if called else 0.0
            if j == len(levels) - 1:
                n_acc.append(len(called))
        return out

    r = C.replicates(SUITE, full=full)
    stream = f"{SUITE}.class_call.{name}"
    gate = calibrate_levels(SUITE, "fdr", levels, fdp, r, stream=stream)
    gate.assert_ok()
    if not C.owns_direct_checks():  # non-triviality is judged once, in shard 0, on block 1
        C.skip_without_verdict(gate.records)
    if C.SHARD_ENV in os.environ:  # shard 0: the loosest-level calls of every block-1 replicate
        n_acc.clear()
        for i in range(r):
            fdp(C.replicate_rng(stream, i))
    assert n_acc and np.mean(n_acc) > 0  # non-trivial: ORFs are called at the loosest level


@pytest.mark.calibration
def test_wired_picked_call_controls_fdp_fast() -> None:
    _check("picked", _picked_replicate, full=False)


@pytest.mark.calibration
def test_wired_denovo_call_controls_fdp_fast() -> None:
    _check("denovo", _denovo_replicate, full=False)


@pytest.mark.calibration_full
def test_wired_picked_call_controls_fdp_full() -> None:
    _check("picked", _picked_replicate, full=True)


@pytest.mark.calibration_full
def test_wired_denovo_call_controls_fdp_full() -> None:
    _check("denovo", _denovo_replicate, full=True)
