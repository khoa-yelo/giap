"""Attack nulls for the M7-11 novel-ORF class FDR (GIAP_SPEC.md §7.6.4, §8.1, §8.5; ADR-0009,
ADR-0010; docs/specs/M7-11-novel-class-fdr.md; statistician, PR #237 review).

Statistics-only generators on the frame layouts of ``tests/helpers/novel_fdr_fixture.py``; the
giapsim attacks of the spec (reversed stretches concentrated in intergenic space, GC-poor false
peptides) stay with M7-13. Criterion ``fdr`` (mean FDP <= q + k SE, ADR-0015 two-block
confirmation), suite ``novel_orfs``, q in ``levels``, streams ``novel_orfs.attack.<name>``.

Picked regime (makedb decoys; a false entry's target and its decoy are exchangeable, best of 1-3
N(0, 1); a true ORF's target is best of 1-3 N(3, 1); input q = 0.01 / (1 + e^score)):

- ``low_prior_class``: intergenic (2000 entries, 3 % true) next to cds_outframe (600, 80 % true).
  FDP is measured **within the intergenic class**: "within each class" (§7.6.4) is the claim, and
  the fast suite pools FDP over classes with identical priors, so pooling would pass it.
  Teeth: the pooled list (every entry relabelled one class) is a confirmed ``fail``.
- ``decoy_loss``: 10 % / 20 % of decoy rows re-pointed to class-tagged IDs absent from
  ``id_map`` (the ADR-0009 attack, novel analogue): the 1/f rescaling (or the 80 % guard) must
  keep the level.
- ``unscored_decoys``: an all-null world where decoy rows carry a score but no input q (C9 allows a
  missing ``q_value``). Today they rank below every scored target, so D = 0 and every target gets
  q = 1/T (the gene-detection hazard, same shape). Remedy-agnostic claim: at most one ORF reaches
  q <= 0.05 when every ORF is false (raising, returning unestimated or ranking on the score all
  pass). ``xfail(strict=True)``: PR #237 block item 1.

De novo regime (ADR-0010 weights; 10 score bins, π̂ = the bin's wrong rate):

- ``low_prior_class``: chance hits fall 80 % intergenic, correct peptides 5 % intergenic; FDP
  within intergenic.
- ``nan_decoy_e``: decoy hits with a missing ``e_decoy`` must not vanish from D. Claim: the call
  raises ``ValueError``, or every q is >= the q obtained with those decoys counted at the best E.
  ``xfail(strict=True)``: PR #237 block item 1.
- ``multi_peptide_orfs`` (decision 3, recorded for M7-13 G4): false ORFs of 1-5 wrong peptides,
  ORF q = the minimum peptide q. Expected anti-conservative (about 2 q at q = 0.01 in review);
  ``xfail(strict=True)`` so that the ADR that fixes decision 3 must remove the marker.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import _criteria as C
import numpy as np
import pandas as pd
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
CLASSES = ("intergenic", "cds_outframe", "antisense")
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
Replicate = tuple[dict[str, float], list[str], set[str]]


def _fdr():  # type: ignore[no-untyped-def]
    import giap.novel.fdr as m

    return m


def _q_in(score: float) -> float:
    return 0.01 / (1.0 + math.exp(score))


def _best(rng: np.random.Generator, mu: float) -> float:
    return float(rng.normal(mu, 1.0, rng.integers(1, 4)).max())


# ---------------------------------------------------------------- picked generator


def _picked_world(
    rng: np.random.Generator,
    sizes: dict[str, int],
    p_true: dict[str, float],
    *,
    p_hit: float = 0.5,
    loss: float = 0.0,
) -> tuple[PickedWorld, set[str], dict[str, str]]:
    w = PickedWorld()
    false_ids: set[str] = set()
    cls_of: dict[str, str] = {}
    n_lost = 0
    for ci, (cls, n) in enumerate(sizes.items()):
        contig = f"k{ci}"
        for _ in range(n):
            true = rng.random() < p_true[cls]
            has_t = true or rng.random() < p_hit
            has_d = rng.random() < p_hit
            if not (has_t or has_d):
                continue
            entry = w.entry(cls, contig, "+")
            if has_t:
                s = _best(rng, 3.0 if true else 0.0)
                orf, _ = w.target(cls, contig, 0, entry=entry, q=_q_in(s), score=s)
                cls_of[orf.orf_id] = cls
                if not true:
                    false_ids.add(orf.orf_id)
            if has_d:
                s = _best(rng, 0.0)
                did = entry[0]
                if rng.random() < loss:  # tagged, parseable, but not an id_map entry
                    n_lost += 1
                    lo = 9_000_001 + 600 * n_lost
                    did = f"GIAPORF_{cls}_{contig}_{lo}_{lo + 299}_p"
                w.decoy(did, 0, q=_q_in(s), score=s)
    return w, false_ids, cls_of


def _picked_q(w: PickedWorld, pooled: bool = False) -> dict[str, float]:
    id_map = w.id_map
    if pooled:  # rejected rule: one list over all classes
        id_map = id_map.assign(**{"class": "intergenic"})
    return dict(_fdr().picked_class_q(w.orfs, w.peptides, id_map))


HETERO_SIZES = {"intergenic": 2000, "cds_outframe": 600}
HETERO_TRUE = {"intergenic": 0.03, "cds_outframe": 0.8}


def _low_prior_picked(pooled: bool) -> Callable[[np.random.Generator], Replicate]:
    def make(rng: np.random.Generator) -> Replicate:
        w, false_ids, cls_of = _picked_world(rng, HETERO_SIZES, HETERO_TRUE)
        ids = [o for o, c in cls_of.items() if c == "intergenic"]
        return _picked_q(w, pooled), ids, false_ids

    return make


def _decoy_loss(loss: float) -> Callable[[np.random.Generator], Replicate]:
    def make(rng: np.random.Generator) -> Replicate:
        w, false_ids, cls_of = _picked_world(rng, {"intergenic": 1500}, {"intergenic": 0.4},
                                             loss=loss)  # fmt: skip
        return _picked_q(w), list(cls_of), false_ids

    return make


# ---------------------------------------------------------------- de novo generator

RATES = np.linspace(0.02, 0.9, 10)
N_DENOVO = 4000
P_CHANCE = 0.25


def _denovo(
    rng: np.random.Generator,
    *,
    chance_w: tuple[float, float, float],
    correct_w: tuple[float, float, float],
    kmax: int = 1,
) -> tuple[list, pd.DataFrame, pd.DataFrame, dict[str, float], set[str]]:  # type: ignore[type-arg]
    """Searched peptides with π̂ = true wrong rate; correct peptides hit their target (low E);
    wrong peptides get a chance target hit, every peptide a chance decoy hit, both w.p.
    ``P_CHANCE`` with the same E law and class weights ``chance_w``. Correct (wrong) peptides
    are grouped into true (false) ORFs of 1..``kmax`` peptides of one class."""
    names = [pep_name(k, "DNVO") for k in range(N_DENOVO)]
    pis = RATES[rng.integers(0, RATES.size, N_DENOVO)]
    wrong = rng.random(N_DENOVO) < pis
    e_t = np.where(wrong, 10 ** rng.normal(-1.5, 1.0, N_DENOVO),
                   10 ** rng.normal(-3.5, 1.0, N_DENOVO))  # fmt: skip
    has_d = rng.random(N_DENOVO) < P_CHANCE
    e_d = 10 ** rng.normal(-1.5, 1.0, N_DENOVO)
    d_cls = rng.choice(3, N_DENOVO, p=chance_w)
    chance_t = wrong & (rng.random(N_DENOVO) < P_CHANCE)
    orfs = []
    false_ids: set[str] = set()
    tp: list[str] = []
    tiv: list[Interval] = []
    te: list[float] = []
    for members, is_false, wts in (
        (rng.permutation(np.flatnonzero(~wrong)), False, correct_w),
        (rng.permutation(np.flatnonzero(chance_t)), True, chance_w),
    ):
        j = 0
        while j < members.size:
            k = int(rng.integers(1, kmax + 1))
            mem = members[j : j + k]
            j += k
            cls = CLASSES[int(rng.choice(3, p=wts))]
            oid = f"{'F' if is_false else 'T'}_{len(orfs)}"
            orfs.append(make_orf(oid, cls, cls, "-", 97 + 1000 * len(orfs),
                                 [names[m] for m in mem]))  # fmt: skip
            if is_false:
                false_ids.add(oid)
            for m in mem:
                tp.append(names[m])
                tiv.append(Interval(cls, 100, 130, "-"))
                te.append(float(e_t[m]))
    di = np.flatnonzero(has_d)
    targets = hit_rows(tp, tiv, te, search="target")
    decoys = hit_rows([names[i] for i in di], [Interval(CLASSES[d_cls[i]], 100, 130, "-")
                                               for i in di], e_d[di], search="decoy")  # fmt: skip
    pi = {p: float(v) for p, v in zip(names, pis, strict=True)}
    return orfs, targets, decoys, pi, false_ids


def _denovo_q(orfs, targets, decoys, pi) -> dict[str, float]:  # type: ignore[no-untyped-def]
    clf = ContigClassifier({c: c for c in CLASSES})
    return dict(_fdr().denovo_class_q(orfs, targets, decoys, clf, pi))


def _low_prior_denovo(rng: np.random.Generator) -> Replicate:
    orfs, t, d, pi, false_ids = _denovo(rng, chance_w=(0.8, 0.1, 0.1),
                                        correct_w=(0.05, 0.5, 0.45))  # fmt: skip
    ids = [o.orf_id for o in orfs if o.cls == "intergenic"]
    return _denovo_q(orfs, t, d, pi), ids, false_ids


def _multi_peptide(rng: np.random.Generator) -> Replicate:
    orfs, t, d, pi, false_ids = _denovo(rng, chance_w=(0.6, 0.25, 0.15),
                                        correct_w=(0.2, 0.4, 0.4), kmax=5)  # fmt: skip
    return _denovo_q(orfs, t, d, pi), [o.orf_id for o in orfs], false_ids


# ---------------------------------------------------------------- gate


def _fdp_fn(make: Callable[[np.random.Generator], Replicate]):  # type: ignore[no-untyped-def]
    levels = C.load_tolerances().levels

    def fdp(rng: np.random.Generator) -> np.ndarray:
        qs, ids, false_ids = make(rng)
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = {o for o in ids if q_of(qs, o) <= a}
            out[j] = len(called & false_ids) / len(called) if called else 0.0
        return out

    return levels, fdp


def _gate(name: str, make: Callable[[np.random.Generator], Replicate], full: bool):  # type: ignore[no-untyped-def]
    levels, fdp = _fdp_fn(make)
    stream = f"{SUITE}.attack.{name}" + (".full" if full else "")
    return calibrate_levels(SUITE, "fdr", levels, fdp, C.replicates(SUITE, full=full),
                            stream=stream)  # fmt: skip


@pytest.mark.parametrize("full", MODES)
def test_picked_fdr_holds_within_a_low_prior_class(full: bool) -> None:
    _gate("picked.low_prior_class", _low_prior_picked(pooled=False), full).assert_ok()


@pytest.mark.parametrize("full", MODES)
def test_low_prior_attack_has_power_against_pooled_classes(full: bool) -> None:
    """Teeth: one list over all classes is a confirmed ``fail`` at q = 0.05 and 0.10 (review:
    FDP 0.29 / 0.46 within intergenic)."""
    gate = _gate("picked.low_prior_class.rejected.pooled", _low_prior_picked(pooled=True), full)
    by_level = {r["level"]: r["status"] for r in gate.records}
    assert by_level.get(0.05) == "fail" and by_level.get(0.1) == "fail", gate.detail()


@pytest.mark.parametrize("loss", [10, 20])
@pytest.mark.parametrize("full", MODES)
def test_picked_fdr_under_random_decoy_loss(full: bool, loss: int) -> None:
    _gate(f"picked.decoy_loss.L{loss}", _decoy_loss(loss / 100), full).assert_ok()


@pytest.mark.parametrize("full", MODES)
def test_denovo_fdr_holds_within_a_low_prior_class(full: bool) -> None:
    _gate("denovo.low_prior_class", _low_prior_denovo, full).assert_ok()


@pytest.mark.xfail(strict=True, reason="M7-11 decision 3: ORF q = min peptide q (M7-13 G4 / ADR)")
@pytest.mark.parametrize("full", MODES)
def test_denovo_orf_fdr_with_multi_peptide_false_orfs(full: bool) -> None:
    _gate("denovo.multi_peptide_orfs", _multi_peptide, full).assert_ok()


# ---------------------------------------------------------------- unscored decoys (hazard)


@pytest.mark.calibration
def test_picked_decoys_without_input_q_do_not_vanish() -> None:
    rng = C.replicate_rng(f"{SUITE}.attack.picked.unscored_decoys", 0)
    w = PickedWorld()
    for _ in range(400):
        entry = w.entry("intergenic", "k0", "+")
        s_t, s_d = _best(rng, 0.0), _best(rng, 0.0)  # all null, exchangeable
        w.target("intergenic", "k0", 0, entry=entry, q=_q_in(s_t), score=s_t)
        w.decoy(entry[0], 0, q=math.nan, score=s_d)
    try:
        qs = _picked_q(w)
    except ValueError:
        return
    n_called = sum(q_of(qs, o.orf_id) <= 0.05 for o in w.orfs)
    assert n_called <= 1, f"{n_called} of {len(w.orfs)} null ORFs at q <= 0.05"


@pytest.mark.calibration
def test_denovo_decoy_hits_without_e_do_not_vanish() -> None:
    rng = C.replicate_rng(f"{SUITE}.attack.denovo.nan_decoy_e", 0)
    orfs, t, d, pi, _ = _denovo(rng, chance_w=(0.6, 0.25, 0.15), correct_w=(0.2, 0.4, 0.4))
    lost = rng.random(len(d)) < 0.3
    d_nan = d.assign(e_decoy=np.where(lost, np.nan, d["e_decoy"].to_numpy()))
    d_best = d.assign(e_decoy=np.where(lost, 0.0, d["e_decoy"].to_numpy()))
    try:
        q_nan = _denovo_q(orfs, t, d_nan, pi)
    except ValueError:
        return
    q_ref = _denovo_q(orfs, t, d_best, pi)
    worse = [o for o in q_ref if q_of(q_nan, o) < q_of(q_ref, o) - 1e-12]
    assert not worse, f"{len(worse)} ORFs gain from dropping decoys with missing E"
