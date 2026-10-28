"""Attack (GIPPA_SPEC.md §8.5, statistician): the free-entrapment estimate of M10-07
(docs/specs/M10-07-free-entrapment.md, "Statistician decisions" S3, S4 item 3, S6).

Null world (one universe, built once from the seed): a K-12-like database of k = 1,000 K-12-only
(entrapment, E) and s = 4,000 shared (target, T) peptides, so r = 0.25 and pi = k/(k+s) = 0.2.
Peptides are grouped into proteins whose sizes are geometric (p = 0.2) truncated to 1-30:

- K-12-specific proteins carry 600 E peptides (every peptide of the protein is E);
- homolog proteins carry 4,400 peptides, of which 400 drawn at random are E (variant peptides)
  and the rest are shared;
- whole homolog proteins, in random order, are "present" until their shared peptides reach
  2,000 (the last protein is trimmed). The other 2,000 shared peptides are absent.

Per replicate, each present peptide is a true discovery with probability 0.9, and M false
discoveries (distinct peptides, never true) are added. E peptides are false by construction. The
true FDP is M / N with N = N_T + N_E. The scenarios differ only in where the false discoveries
land:

- ``A`` edge equal chance: M ~ Poisson(30); each false discovery is in E with probability pi,
  otherwise on an absent shared peptide. This is the model ``ucb95`` inverts.
- ``B`` uniform over the null elements E and absent shared (P(E) = 1/3 > pi).
- ``C`` homolog attraction: as A with the E weight multiplied by lambda in {1, 2, 4}.
- ``D`` attraction to T: lambda = 0.5.
- ``E`` clustered: false discoveries arrive as whole proteins; Poisson(30 / mean size) clusters,
  each with probability pi a K-12-specific protein (all its E peptides), otherwise an absent
  homolog protein (all its shared peptides). Proteins are drawn without replacement.
- ``F`` a controlled search: as A with M ~ Poisson(9), so the true FDR is about 0.5 alpha. In
  A-E the true FDP is about 1.6 % (> alpha), so only 20-30 replicates per scenario have
  FDP <= alpha; F characterises the labels where the search does control.

Gated (fast R = 2,000; full R = 20,000, ``calibration_full``, Slurm), k = ``se_multiplier``:

- A, B, C(lambda = 1), E: mean(``fdp_combined`` - FDP) >= -k SE (the paper's averaged bound);
- C: the mean excess is strictly increasing in lambda, each step by more than k SE;
- D: mean(``fdp_combined`` - FDP) + k SE < 0 (the low-direction row of the bias table);
- every scenario and replicate: ``fdp_lower`` <= FDP (exact, in Fractions), and no replicate with
  FDP <= alpha is labelled ``evidence_of_failure``;
- A, B: P(FDP <= ``ucb95``) >= 0.95 - k sqrt(0.95 0.05 / R).

Recorded, not gated (``compare_attack_free_entrapment_recorded.jsonl``): ``ucb95`` coverage in
C, D, E and F (E is the reason its ``assumes`` string excludes clustered errors), and the label
rates, including F's (how often a controlled search reads ``consistent_with_control`` or
``inconclusive``). F gates only the deterministic per-replicate checks.

Teeth: in A, ``fdp_lower`` used as an estimate fails the mean check, and a median bound
(``ucb95(level=0.5)``) fails the coverage check. Without them, a check that cannot fail would
look like a pass.

Calls the implementation's pure API (``combined``, ``lower``, ``ucb95``, ``assess`` in
``validation/pilot/entrapment.py``) unchanged; the test fails, not errors, while it does not
exist. Randomness only through ``gippa.rng.for_step``, seeded ``seeds.calibration_base + 707``
from ``tolerances.yaml``.
"""

from __future__ import annotations

import functools
import json
import math
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import _criteria as C
import numpy as np
import pytest

from gippa.rng import for_step

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from helpers.entrapment_m10_07 import load_entrapment

SUITE = "compare"
pytestmark = [pytest.mark.calibration, pytest.mark.calsuite(SUITE)]
MODES = [pytest.param(False, id="fast"),
         pytest.param(True, marks=pytest.mark.calibration_full, id="full")]  # fmt: skip

R_FAST = 2_000
R_FULL = 20_000
SEED_OFFSET = 707
ALPHA = Fraction(1, 100)

N_K12_ONLY = 1_000
N_SHARED = 4_000
N_E_SPECIFIC = 600  # E peptides on K-12-specific proteins; the other 400 are variant peptides
N_PRESENT = 2_000
P_TRUE = 0.9
MEAN_FALSE = 30.0
MEAN_FALSE_CONTROLLED = 9.0  # F: true FDR about 0.5 alpha
GEOM_P = 0.2
MAX_SIZE = 30
PI = Fraction(N_K12_ONLY, N_K12_ONLY + N_SHARED)

#: (scenario, lambda) pairs; lambda multiplies the E weight of A
RUNS = {"A": 1.0, "B": None, "C1": 1.0, "C2": 2.0, "C4": 4.0, "D": 0.5, "E": None, "F": 1.0}


def _seed() -> int:
    return C.load_tolerances().calibration_base + SEED_OFFSET


def _k() -> float:
    return C.load_tolerances().se_multiplier


def _R(full: bool) -> int:
    return R_FULL if full else R_FAST


# ------------------------------------------------------------------ the universe


def _sizes(rng: np.random.Generator, total: int) -> list[int]:
    """Truncated-geometric protein sizes summing to exactly ``total`` (the last one trimmed)."""
    out: list[int] = []
    while sum(out) < total:
        size = int(rng.geometric(GEOM_P))
        if size <= MAX_SIZE:
            out.append(min(size, total - sum(out)))
    return out


@dataclass(frozen=True)
class Universe:
    e_proteins: tuple[np.ndarray, ...]  # K-12-specific proteins: their E peptide ids
    absent_proteins: tuple[np.ndarray, ...]  # absent homolog proteins: their shared peptide ids
    e_ids: np.ndarray  # every E peptide id
    absent_shared: np.ndarray  # every absent shared peptide id
    present: np.ndarray  # the N_PRESENT truly present shared peptide ids
    mean_size: float


@functools.cache
def universe() -> Universe:
    rng = for_step("attack.free_entrapment.universe", _seed())
    next_id = 0

    def ids(n: int) -> np.ndarray:
        nonlocal next_id
        a = np.arange(next_id, next_id + n)
        next_id += n
        return a

    e_prot = [ids(n) for n in _sizes(rng, N_E_SPECIFIC)]
    hom = [ids(n) for n in _sizes(rng, N_SHARED + N_K12_ONLY - N_E_SPECIFIC)]
    hom_all = np.concatenate(hom)
    variant = set(rng.choice(hom_all, N_K12_ONLY - N_E_SPECIFIC, replace=False).tolist())
    hom_shared = [np.array([p for p in prot if p not in variant], dtype=int) for prot in hom]
    present: list[int] = []
    absent_prot: list[np.ndarray] = []
    for j in rng.permutation(len(hom)):
        sh = hom_shared[j]
        if len(present) >= N_PRESENT:
            if sh.size:
                absent_prot.append(sh)
        else:
            take = min(sh.size, N_PRESENT - len(present))
            present.extend(sh[:take].tolist())  # a trimmed protein is not an absent protein
    e_ids = np.concatenate([*e_prot, np.array(sorted(variant), dtype=int)])
    absent = np.setdiff1d(np.concatenate(hom_shared), np.array(present))
    assert e_ids.size == N_K12_ONLY and len(present) == N_PRESENT
    assert absent.size == N_SHARED - N_PRESENT
    sizes = [p.size for p in e_prot + hom]
    return Universe(tuple(e_prot), tuple(absent_prot), e_ids, absent, np.array(present),
                    float(np.mean(sizes)))  # fmt: skip


# ------------------------------------------------------------------ one replicate


def _false_counts(u: Universe, scenario: str, lam: float | None,
                  rng: np.random.Generator) -> tuple[int, int]:  # fmt: skip
    """(false discoveries in E, false discoveries in T), distinct peptides."""
    if scenario == "E":
        n_cl = int(rng.poisson(MEAN_FALSE / u.mean_size))
        is_e = rng.random(n_cl) < float(PI)
        n_ce = int(is_e.sum())
        e_pick = rng.choice(len(u.e_proteins), min(n_ce, len(u.e_proteins)), replace=False)
        t_pick = rng.choice(len(u.absent_proteins), min(n_cl - n_ce, len(u.absent_proteins)),
                            replace=False)  # fmt: skip
        return (int(sum(u.e_proteins[j].size for j in e_pick)),
                int(sum(u.absent_proteins[j].size for j in t_pick)))  # fmt: skip
    m = int(rng.poisson(MEAN_FALSE_CONTROLLED if scenario == "F" else MEAN_FALSE))
    if scenario == "B":
        pool = u.e_ids.size + u.absent_shared.size
        hit = rng.choice(pool, m, replace=False)
        m_e = int(np.count_nonzero(hit < u.e_ids.size))
        return m_e, m - m_e
    assert lam is not None
    w = lam * float(PI)
    p_e = w / (w + 1.0 - float(PI))
    m_e = int(rng.binomial(m, p_e))
    m_e = min(m_e, u.e_ids.size)
    return m_e, min(m - m_e, u.absent_shared.size)


@dataclass(frozen=True)
class Run:
    fdp: np.ndarray
    combined: np.ndarray
    lower: np.ndarray
    ucb95: np.ndarray
    ucb50: np.ndarray
    labels: np.ndarray
    lower_ok: bool
    failure_at_control: int
    n_at_control: int
    seconds: float


@functools.cache
def run(name: str, R: int) -> Run:
    E = load_entrapment()
    u = universe()
    lam = RUNS[name]
    scenario = name[0]
    rng = for_step(f"attack.free_entrapment.{name}", _seed())
    t0 = time.perf_counter()
    fdp, comb, low, ub, ub50 = (np.empty(R) for _ in range(5))
    labels: list[str] = []
    lower_ok, fail_ctl, n_ctl = True, 0, 0
    for i in range(R):
        n_true = int(rng.binomial(N_PRESENT, P_TRUE))
        m_e, m_t = _false_counts(u, scenario, lam, rng)
        n_t, n_e = n_true + m_t, m_e
        n = n_t + n_e
        true_fdp = Fraction(m_e + m_t, n)
        c = E.combined(n_t, n_e, N_K12_ONLY, N_SHARED)
        lo = E.lower(n_t, n_e)
        label = E.assess(c, lo)
        lower_ok &= Fraction(lo) <= true_fdp
        if true_fdp <= ALPHA:
            n_ctl += 1
            fail_ctl += label == "evidence_of_failure"
        fdp[i], comb[i], low[i] = float(true_fdp), float(c), float(lo)
        ub[i] = E.ucb95(n_t, n_e, N_K12_ONLY, N_SHARED)
        ub50[i] = E.ucb95(n_t, n_e, N_K12_ONLY, N_SHARED, level=0.5) if scenario == "A" else 0.0
        labels.append(label)
    return Run(fdp, comb, low, ub, ub50, np.array(labels), bool(lower_ok), fail_ctl, n_ctl,
               time.perf_counter() - t0)  # fmt: skip


# ------------------------------------------------------------------ criteria


def _mean_se(d: np.ndarray) -> tuple[float, float]:
    return float(d.mean()), float(d.std(ddof=1) / math.sqrt(d.size))


def _coverage(r: Run, ub: np.ndarray) -> float:
    return float(np.mean(r.fdp <= ub + 1e-12))


def _coverage_bound(R: int) -> float:
    return 0.95 - _k() * math.sqrt(0.95 * 0.05 / R)


def _log_recorded(record: dict[str, object]) -> None:
    path = C.log_dir() / f"{SUITE}_attack_free_entrapment_recorded.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(record, sort_keys=True) + "\n")


def _gate(name: str, ok: bool, observed: float, bound: float, R: int, detail: str,
          r: Run, failures: list[str]) -> None:  # fmt: skip
    res = C.CheckResult(ok, observed, bound, R, detail)
    C.log_result(SUITE, f"free_entrapment_attack.{name}", float(ALPHA), res, _seed(), r.seconds)
    if not ok:
        failures.append(f"{name}: {detail}")


def _per_replicate_checks(name: str, r: Run, R: int, failures: list[str]) -> None:
    _gate(f"{name}.lower_le_fdp", r.lower_ok, float(r.lower_ok), 1.0, R,
          "fdp_lower <= true FDP in every replicate", r, failures)  # fmt: skip
    _gate(f"{name}.no_failure_label_under_control", r.failure_at_control == 0,
          float(r.failure_at_control), 0.0, R,
          f"{r.failure_at_control} of {r.n_at_control} replicates with FDP <= alpha labelled "
          "evidence_of_failure", r, failures)  # fmt: skip
    values, counts = np.unique(r.labels, return_counts=True)
    _log_recorded({"run": name, "R": R, "mean_fdp": float(r.fdp.mean()),
                   "mean_combined": float(r.combined.mean()),
                   "ucb95_coverage": _coverage(r, r.ucb95), "n_fdp_le_alpha": r.n_at_control,
                   "labels": dict(zip(map(str, values), map(int, counts), strict=True)),
                   })  # fmt: skip


@pytest.mark.parametrize("name", ["A", "B", "C1", "E"])
@pytest.mark.parametrize("full", MODES)
def test_combined_is_an_upper_bound_in_mean(full: bool, name: str) -> None:
    """A, B, C(lambda=1), E: mean(fdp_combined) >= mean(FDP) - k SE; plus the per-replicate
    checks."""
    R = _R(full)
    r = run(name, R)
    failures: list[str] = []
    mean, se = _mean_se(r.combined - r.fdp)
    _gate(f"{name}.mean_bias", mean >= -_k() * se, mean, -_k() * se, R,
          f"mean(fdp_combined - FDP) = {mean:.3g}, SE {se:.3g}", r, failures)  # fmt: skip
    _per_replicate_checks(name, r, R, failures)
    assert not failures, failures


@pytest.mark.parametrize("name", ["A", "B"])
@pytest.mark.parametrize("full", MODES)
def test_ucb95_coverage_under_equal_chance(full: bool, name: str) -> None:
    R = _R(full)
    r = run(name, R)
    cov, bound = _coverage(r, r.ucb95), _coverage_bound(R)
    failures: list[str] = []
    _gate(f"{name}.ucb95_coverage", cov >= bound, cov, bound, R,
          f"P(FDP <= ucb95) = {cov:.4f} vs {bound:.4f}", r, failures)  # fmt: skip
    assert not failures, failures


@pytest.mark.parametrize("full", MODES)
def test_homolog_attraction_excess_grows_with_lambda(full: bool) -> None:
    """C: mean(fdp_combined - FDP) strictly increasing over lambda 1, 2, 4."""
    R = _R(full)
    names = ["C1", "C2", "C4"]
    stats = [_mean_se(run(n, R).combined - run(n, R).fdp) for n in names]
    failures: list[str] = []
    for (a, sa), (b, sb), n in zip(stats, stats[1:], names[1:], strict=False):
        step, bound = b - a, _k() * math.hypot(sa, sb)
        _gate(f"{n}.excess_step", step > bound, step, bound, R,
              f"excess step to {n}: {step:.3g} vs {bound:.3g}", run(n, R), failures)  # fmt: skip
    for n in names[1:]:
        _per_replicate_checks(n, run(n, R), R, failures)
    assert not failures, failures


@pytest.mark.parametrize("full", MODES)
def test_attraction_to_targets_biases_low(full: bool) -> None:
    """D (lambda = 0.5): mean(fdp_combined) < mean(FDP), by more than k SE."""
    R = _R(full)
    r = run("D", R)
    mean, se = _mean_se(r.combined - r.fdp)
    failures: list[str] = []
    _gate("D.mean_bias_low", mean + _k() * se < 0, mean, -_k() * se, R,
          f"mean(fdp_combined - FDP) = {mean:.3g}, SE {se:.3g}", r, failures)  # fmt: skip
    _per_replicate_checks("D", r, R, failures)
    assert not failures, failures


@pytest.mark.parametrize("full", MODES)
def test_controlled_search_is_recorded(full: bool) -> None:
    """F (true FDR about 0.5 alpha): label rates and ucb95 coverage recorded, not gated; the
    deterministic per-replicate checks gate, now on a population mostly at FDP <= alpha."""
    R = _R(full)
    r = run("F", R)
    failures: list[str] = []
    _per_replicate_checks("F", r, R, failures)
    assert r.n_at_control > R // 2, f"F: only {r.n_at_control} of {R} replicates at FDP <= alpha"
    assert not failures, failures


# ------------------------------------------------------------------ teeth


@pytest.mark.parametrize("full", MODES)
def test_teeth_lower_bound_fails_the_mean_check(full: bool) -> None:
    R = _R(full)
    r = run("A", R)
    mean, se = _mean_se(r.lower - r.fdp)
    assert mean < -_k() * se, f"fdp_lower passes the mean check ({mean:.3g}, SE {se:.3g})"


@pytest.mark.parametrize("full", MODES)
def test_teeth_median_bound_fails_the_coverage_check(full: bool) -> None:
    R = _R(full)
    r = run("A", R)
    cov = _coverage(r, r.ucb50)
    assert cov < _coverage_bound(R), f"ucb95(level=0.5) covers {cov:.4f}: the check cannot bite"


def test_universe_structure() -> None:
    """The universe keeps the protein grouping and the counts the spec fixes."""
    u = universe()
    assert u.e_ids.size == N_K12_ONLY and np.unique(u.e_ids).size == N_K12_ONLY
    assert u.present.size == N_PRESENT and u.absent_shared.size == N_SHARED - N_PRESENT
    assert not set(u.present.tolist()) & set(u.absent_shared.tolist())
    assert not set(u.e_ids.tolist()) & set(u.absent_shared.tolist())
    sizes = [p.size for p in u.e_proteins]
    assert sum(sizes) == N_E_SPECIFIC and max(sizes) <= MAX_SIZE
    assert all(p.size <= MAX_SIZE for p in u.absent_proteins)
    assert sum(p.size for p in u.absent_proteins) <= u.absent_shared.size
    assert 3.0 < u.mean_size < 6.0
