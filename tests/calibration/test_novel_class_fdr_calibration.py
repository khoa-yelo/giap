"""Statistics-only calibration of the M7-11 novel-ORF class FDR (GIAP_SPEC.md §7.6.4, §8.4 row
"Novel-ORF class FDR"; docs/specs/M7-11-novel-class-fdr.md; ADR-0009, ADR-0010, ADR-0015).

A fast check on synthetic lists, one per regime; the giapsim suite with planted ORFs, false
peptides and de novo errors (and the two new attacks) is M7-13's. Criterion: mean FDP over R
replicates <= q + 2 SE at q in ``levels``; R = ``replicates_fast.novel_orfs`` (full:
``replicates_full.novel_orfs``), suite ``novel_orfs``, two-block confirmation via
``_gate.calibrate_levels``.

Picked generator (three makedb classes, sizes 1500 / 600 / 600): an entry is a true novel ORF with
probability 0.4. A true ORF always has its target, best of 1-3 peptides with scores
N(3, 1); a false entry has a target with probability ``P_HIT``. Every entry has a decoy with
probability ``P_HIT``. False targets and all decoys draw best-of-k null scores N(0, 1) with the
same k law, so a false target and its decoy are exchangeable. Input q = 0.01 / (1 + e^score)
(all <= ``input.q_value_max``; ordering = score). FDP = false ORFs among ORFs with class q <= a,
pooled over classes.

De novo generator (single-peptide ORFs so ORF level = peptide level; decision 3's ORF-level FDP
is M7-13's): N peptides in 10 score bins with wrong rates 0.02..0.9 and π̂ = the bin's rate
(ADR-0010 scenario). A correct peptide has its target hit (E ~ 10^N(-3.5, 1)); a wrong one has a
chance target hit with probability ``P_CHANCE``; every peptide has a chance decoy hit with
probability ``P_CHANCE``. Chance hits (E ~ 10^N(-1.5, 1)) fall into the three classes with the
same weights for targets and decoys; correct peptides are spread with other weights.
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels
from helpers.novel_fdr_fixture import (
    ContigClassifier,
    PickedWorld,
    hit_rows,
    make_orf,
    pep_name,
    q_of,
)

from giap.model.coords import Interval

pytestmark = pytest.mark.calsuite("novel_orfs")

SUITE = "novel_orfs"
CLASS_SIZES = {"intergenic": 1500, "cds_outframe": 600, "antisense": 600}
P_TRUE = 0.4
P_HIT = 0.5
N_DENOVO = 4000
P_CHANCE = 0.25
CHANCE_W = np.array([0.6, 0.25, 0.15])
CORRECT_W = np.array([0.2, 0.4, 0.4])
CLASSES = tuple(CLASS_SIZES)


def _fdr():  # type: ignore[no-untyped-def]
    import giap.novel.fdr as m  # fails (ModuleNotFoundError) until M7-11 lands

    return m


def _q_in(score: float) -> float:
    return 0.01 / (1.0 + float(np.exp(score)))


def _best(rng: np.random.Generator, mu: float) -> float:
    return float(rng.normal(mu, 1.0, rng.integers(1, 4)).max())


def _picked_replicate(rng: np.random.Generator) -> tuple[dict[str, float], set[str], set[str]]:
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
            entry = w.entry(cls, contig, "+")
            if has_t:
                s = _best(rng, 3.0 if true else 0.0)
                orf, _ = w.target(cls, contig, 0, entry=entry, q=_q_in(s), score=s)
                if not true:
                    false_ids.add(orf.orf_id)
            if has_d:
                s = _best(rng, 0.0)
                w.decoy(entry[0], 0, q=_q_in(s), score=s)
    qs = _fdr().picked_class_q(w.orfs, w.peptides, w.id_map)
    return qs, {o.orf_id for o in w.orfs}, false_ids


def _denovo_replicate(rng: np.random.Generator) -> tuple[dict[str, float], set[str], set[str]]:
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
    ti, di = np.flatnonzero(has_t), np.flatnonzero(has_d)
    iv = [Interval(CLASSES[t_cls[i]], 100, 130, "-") for i in ti]
    targets = hit_rows([peps[i] for i in ti], iv, e_t[ti], search="target")
    dv = [Interval(CLASSES[d_cls[i]], 100, 130, "-") for i in di]
    decoys = hit_rows([peps[i] for i in di], dv, e_d[di], search="decoy")
    orfs = [make_orf(f"N_{i}", CLASSES[t_cls[i]], CLASSES[t_cls[i]], "-", 97, [peps[i]])
            for i in ti]  # fmt: skip
    pi = {p: float(v) for p, v in zip(peps, pi_true, strict=True)}
    qs = _fdr().denovo_class_q(orfs, targets, decoys, ContigClassifier({c: c for c in CLASSES}), pi)
    return qs, {o.orf_id for o in orfs}, {f"N_{i}" for i in ti if wrong[i]}


def _check(
    name: str,
    make: Callable[[np.random.Generator], tuple[dict[str, float], set[str], set[str]]],
    full: bool,
) -> None:
    levels = C.load_tolerances().levels
    n_acc: list[int] = []

    def fdp(rng: np.random.Generator) -> np.ndarray:
        qs, ids, false_ids = make(rng)
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = {o for o in ids if q_of(qs, o) <= a}  # NaN never passes
            out[j] = len(called & false_ids) / len(called) if called else 0.0
            if j == len(levels) - 1:
                n_acc.append(len(called))
        return out

    r = C.replicates(SUITE, full=full)
    calibrate_levels(SUITE, "fdr", levels, fdp, r, stream=f"{SUITE}.class_fdr.{name}").assert_ok()
    assert not n_acc or np.mean(n_acc) > 0  # non-trivial: ORFs are called at the loosest level


@pytest.mark.calibration
def test_picked_class_fdr_controls_fdp_fast() -> None:
    _check("picked", _picked_replicate, full=False)


@pytest.mark.calibration
def test_denovo_subgroup_fdr_controls_fdp_fast() -> None:
    _check("denovo", _denovo_replicate, full=False)


@pytest.mark.calibration_full
def test_picked_class_fdr_controls_fdp_full() -> None:
    _check("picked", _picked_replicate, full=True)


@pytest.mark.calibration_full
def test_denovo_subgroup_fdr_controls_fdp_full() -> None:
    _check("denovo", _denovo_replicate, full=True)
