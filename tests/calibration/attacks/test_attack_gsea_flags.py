"""Calibration and attack arms of the M8-05b GSEA row flags (M8-05 spec, lead decision 9;
ADR-0025, ADR-0026; GIPPA_SPEC.md §8.4, §8.5).

The strict xfails on *all* rows (one run of 15 under AR(1) phi 0.5; detectability confounded with
membership) stay in test_attack_gsea_structure.py. This file gates what the flags claim, in the
worlds of the statistician's probe (local/stat/m805flags/flags_probe.py; Slurm 1866662):

* **Unflagged rows hold level** (gated, levels 0.01 / 0.05 / 0.10). A replicate draws worlds until
  the set's row exists and does not carry the arm's flag, and returns that row's p: the replicate
  values are then a sample of p among unflagged rows. Probe rates of the unflagged rows:
  - ``genomic_cluster``: phi 0.5 mix 0.0104 / 0.0519 / 0.0961; phi 0.2 mix 0.0088 / 0.0492 /
    0.0980; phi 0 mix 0.0102 / 0.0451 / 0.1013; phi 0.5, 15 scattered 0.0110 / 0.0513 / 0.1034.
  - ``detection_confounded``: the mix (s ~ U(0, 1.5), beta in {1, 2}, 15-120 TUs) 0.0106 / 0.0544
    / 0.1025 (sigma inf); null s 0, beta 2 0.0089 / 0.0533 / 0.0971; null s 0, beta 1, 15-120 TUs
    0.0123 / 0.0544 / 0.1080; s 1.0, beta 1 0.0028 / 0.0100 / 0.0213.
    Statistician #174 (Slurm 1867332, R 4,000): the detection flag is a screen, so its arms are
    gated at rejection <= 1.5 alpha (criterion ``fwer`` on the rejection indicator), with teeth:
    all rows of the mix reject at >= 3 alpha.
* **Teeth** (gated): one run of 15 and five runs of 3 are flagged ``genomic_cluster`` in >= 99 %
  of replicates; ``detection_confounded`` covers >= 99 % of s 1.5 rows and >= 75 % of s 1.0 rows.
  A coverage claim c passes when the observed share is >= c - 2 sqrt(c (1 - c) / R) (the §8.4
  bound mirrored), so fast mode (R 50) does not fail on Monte Carlo noise alone.
* **Flag rates under the null.** ``genomic_cluster`` on random placement: <= 0.002 on a
  2,816-TU ring (all E. coli K-12 operons ranked) and on 1,710 ranked TUs (60 %), for every
  k in {10, 15, 30, 60, 100, 200, 300, 500}; pooled level check of 200 placements per replicate.
  ``detection_confounded`` in the unconfounded null: <= 0.25 (probe 0.133-0.160).

Not gated, recorded in the spec's Limits: phi 0.7 (unflagged rows 0.0165 at alpha 0.01, the
cluster flag is calibrated for phi <= 0.5 only) and the marginal sigma-0 mix at alpha 0.01
(0.0147 vs 0.0142). The covariate sigma of the probe enters neither the ranking nor the flag
(``p_detect`` only decides eligibility), so one mix arm stands for sigma 0 / 1 / inf; ``p_detect``
is the overall detected fraction here (an uninformative fit, sigma inf).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from functools import cache
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers.gsea_data import as_list, gene_sets

from gippa.config.model import GseaConfig
from gippa.rng import for_step

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("gsea")]

SUITE = "gsea"
CFG = GseaConfig()
SIZE = (10, 500)
ALPHA = 0.05
MAX_DRAWS = 400
REPLICONS = {"chr": "circular"}

FAST_FULL = [
    pytest.param(False, id="fast"),
    pytest.param(True, id="full", marks=pytest.mark.calibration_full),
]
#: the null flag-rate checks pool 200 placements per replicate (k up to 500): nightly in fast mode
NULL_FAST_FULL = [
    pytest.param(False, id="fast", marks=pytest.mark.calibration_nightly),
    pytest.param(True, id="full", marks=pytest.mark.calibration_full),
]


def _g() -> Any:
    from gippa.function import gsea

    return gsea


# --------------------------------------------------------------------------- worlds


def _world(
    rng: np.random.Generator, n: int, rho: float, phi: float
) -> tuple[np.ndarray, np.ndarray]:
    """E. coli-like operon lengths (geometric, ~half singletons), TU effect AR(1) phi along the TU
    order, within-TU correlation rho; gene -> TU index and gene log2 abundance."""
    e = rng.standard_normal(n)
    t = e.copy()
    for i in range(1, n):
        t[i] = phi * t[i - 1] + math.sqrt(1 - phi**2) * e[i]
    tu = np.repeat(np.arange(n), np.minimum(rng.geometric(0.48, n), 14))
    x = math.sqrt(rho) * t[tu] + math.sqrt(1 - rho) * rng.standard_normal(tu.size) + 20.0
    return tu, x


def _frame(tu: np.ndarray, x: np.ndarray, det: np.ndarray | None = None) -> pd.DataFrame:
    """Gene frame on one contig ``chr``: gene g<i+1> in TU_<tu+1>, positions in gene order."""
    d = np.ones(tu.size, dtype=bool) if det is None else det.astype(bool)
    pos = pd.Series(tu).groupby(tu).cumcount().to_numpy() + 1
    pdet = float(d.mean())
    return pd.DataFrame(
        {
            "locus_tag": [f"g{i + 1}" for i in range(tu.size)],
            "tu_id": [f"TU_{int(t) + 1}" for t in tu],
            "tu_position": pos.astype(np.int64),
            "contig": "chr",
            "detected": d,
            "abundance_log2": np.where(d, x, np.nan).astype(np.float64),
            "abundance_imputed": False,
            "abundance_scale": "molar",
            "partial": False,
            "p_detect": pdet,
        }
    )


def _row(genes: pd.DataFrame, members: list[str], rng: np.random.Generator) -> pd.Series | None:
    seed = int(rng.integers(0, 2**31 - 1))
    out = _g().gsea_table(genes, [gene_sets("kegg", {"S": members})], SIZE, CFG, ALPHA,
                          for_step("function.gsea", seed), replicons=REPLICONS)  # fmt: skip
    if len(out) == 0:
        return None
    assert len(out) == 1, out
    return out.iloc[0]


def _p(r: pd.Series) -> float:
    p = float(r["p_value"])
    return 1.0 if math.isnan(p) else p  # an unestimated row does not reject


def _first_genes(tu: np.ndarray, pick: np.ndarray) -> list[str]:
    first = {int(t): g for g, t in reversed(list(enumerate(tu)))}
    return [f"g{first[int(t)] + 1}" for t in pick]


# ------------------------------------------------------------------- shapes (flags_probe.py)


def _shape(rng: np.random.Generator, n: int, kind: str, k: int) -> np.ndarray:
    if kind == "scatter":
        return rng.choice(n, k, replace=False)
    if kind == "run":
        return (int(rng.integers(0, n)) + np.arange(k)) % n
    if kind.startswith("runs"):  # runsL:r -> r runs of length L, the rest scattered
        L, r = (int(v) for v in kind[4:].split(":"))
        starts = rng.choice(n // (L + 1), r, replace=False) * (L + 1)
        pos = set(np.concatenate([s + np.arange(L) for s in starts]).tolist())
        rest = [i for i in rng.permutation(n).tolist() if i not in pos]
        out = list(pos) + rest[: max(0, k - len(pos))]
        return np.array(out[:k])
    raise ValueError(kind)


def _block_plus_scatter(rng: np.random.Generator, n: int, k: int, step: int) -> np.ndarray:
    m = int(rng.integers(2, k + 1))
    blk = (int(rng.integers(0, n)) + step * np.arange(m)) % n
    taken = set(blk.tolist())
    rest = [i for i in rng.permutation(n).tolist() if i not in taken][: k - m]
    return np.concatenate([blk, np.array(rest, dtype=int)])


def _mix(rng: np.random.Generator, n: int) -> np.ndarray:
    k = int(rng.integers(10, 41))
    t = int(rng.integers(0, 6))
    if t == 0:
        return _shape(rng, n, "scatter", k)
    if t == 1:
        return _shape(rng, n, f"runs2:{int(rng.integers(1, k // 2 + 1))}", k)
    if t == 2:
        L = int(rng.integers(3, 6))
        return _shape(rng, n, f"runs{L}:{int(rng.integers(1, k // L + 1))}", k)
    if t == 3:
        return _block_plus_scatter(rng, n, k, 2)  # a spaced block plus scattered TUs
    if t == 4:
        return _block_plus_scatter(rng, n, k, 1)  # one run plus scattered TUs
    return _shape(rng, n, "run", k)


def _positions(rng: np.random.Generator, n: int, kind: str) -> np.ndarray:
    if kind == "mix":
        return np.unique(_mix(rng, n))
    if kind == "run15":
        return np.unique(_shape(rng, n, "run", 15))
    if kind == "5x3":
        return np.unique(_shape(rng, n, "runs3:5", 15))
    return np.unique(_shape(rng, n, "scatter", 15))


def _delta(pos: np.ndarray, n: int) -> float:
    return float(_g().cluster_delta({"chr": sorted(pos.tolist())}, {"chr": n}, {"chr": True}))


# --------------------------------------------------------------------------- genomic_cluster


N_AR = 1500


def _cluster_unflagged(phi: float, kind: str) -> Callable[[np.random.Generator], float]:
    """p of an unflagged row: shapes are drawn until delta < tau (the world is independent of the
    shape, so drawing it afterwards leaves the conditional law unchanged)."""

    def rep(rng: np.random.Generator) -> float:
        g = _g()
        for _ in range(MAX_DRAWS):
            pos = _positions(rng, N_AR, kind)
            if g.is_genomic_cluster(_delta(pos, N_AR)):
                continue
            tu, x = _world(rng, N_AR, 0.8, phi)
            r = _row(_frame(tu, x), _first_genes(tu, pos), rng)
            assert r is not None
            assert "genomic_cluster" not in as_list(r["flags"]), (pos, r["flags"])
            return _p(r)
        raise AssertionError(f"no unflagged row in {MAX_DRAWS} draws")

    return rep


CLUSTER_LEVEL = {
    "phi0.5_mix": (0.5, "mix"),
    "phi0.2_mix": (0.2, "mix"),
    "phi0_mix": (0.0, "mix"),
    "phi0.5_scattered15": (0.5, "scatter"),
}


@pytest.mark.parametrize("full", FAST_FULL)
@pytest.mark.parametrize("arm", sorted(CLUSTER_LEVEL))
def test_genomic_cluster_unflagged_rows_at_level(arm: str, full: bool) -> None:
    levels = C.load_tolerances().levels
    reps = C.replicates(SUITE, full=full)
    rec = calibrate_levels(SUITE, "level", levels, _cluster_unflagged(*CLUSTER_LEVEL[arm]), reps,
                           stream=f"{SUITE}.attack.flags.cluster.{arm}")  # fmt: skip
    rec.assert_ok()


def _coverage_ok(hits: list[bool], claim: float) -> tuple[bool, str]:
    R = len(hits)
    share = float(np.mean(hits))
    bound = claim - 2.0 * math.sqrt(claim * (1 - claim) / R)
    return share >= bound, f"share {share:.4f} vs >= {bound:.4f} (claim {claim}, R {R})"


@pytest.mark.parametrize("full", FAST_FULL)
@pytest.mark.parametrize("kind", ["run15", "5x3"])
def test_genomic_cluster_teeth_clustered_sets_are_flagged(kind: str, full: bool) -> None:
    def rep(rng: np.random.Generator) -> bool:
        pos = _positions(rng, N_AR, kind)
        tu, x = _world(rng, N_AR, 0.8, 0.5)
        r = _row(_frame(tu, x), _first_genes(tu, pos), rng)
        assert r is not None
        return "genomic_cluster" in as_list(r["flags"])

    hits = C.run_replicates(f"{SUITE}.attack.flags.cluster_teeth.{kind}", rep,
                            C.replicates(SUITE, full=full))  # fmt: skip
    ok, msg = _coverage_ok(hits, 0.99)
    assert ok, msg


@cache
def _ring_rank(n: int) -> tuple[Any, pd.DataFrame]:
    tu = np.arange(n)
    x = 20.0 + 0.001 * np.arange(n)
    genes = _frame(tu, x)
    return _g().tu_ranking(genes), genes


PLACEMENTS = 200
NULL_RATE = 0.002


def _null_flag_share(n: int, k: int, through_frame: bool) -> Callable[[np.random.Generator], float]:
    def rep(rng: np.random.Generator) -> float:
        g = _g()
        if through_frame:
            rank, genes = _ring_rank(n)
            sets = {f"S{j}": frozenset(rank.tu_ids[int(i)] for i in rng.choice(n, k, replace=False))
                    for j in range(PLACEMENTS)}  # fmt: skip
            deltas = list(g.cluster_deltas(rank, genes, REPLICONS, sets).values())
        else:
            deltas = [_delta(rng.choice(n, k, replace=False), n) for _ in range(PLACEMENTS)]
        assert len(deltas) == PLACEMENTS
        return sum(bool(g.is_genomic_cluster(float(d))) for d in deltas) / PLACEMENTS

    return rep


@pytest.mark.parametrize("full", NULL_FAST_FULL)
@pytest.mark.parametrize("k", [10, 15, 30, 60, 100, 200, 300, 500])
@pytest.mark.parametrize(("n", "through_frame"), [(2816, True), (1710, False)],
                         ids=["ring2816", "ranked1710"])  # fmt: skip
def test_genomic_cluster_null_flag_rate(n: int, through_frame: bool, k: int, full: bool) -> None:
    rec = C.calibrate(SUITE, "level", NULL_RATE, _null_flag_share(n, k, through_frame),
                      C.replicates(SUITE, full=full), m=PLACEMENTS,
                      stream=f"{SUITE}.attack.flags.cluster_null.n{n}.k{k}")  # fmt: skip
    assert rec["status"] != "fail", rec


# --------------------------------------------------------------------- detection_confounded


N_DET = 2000


def _det_world(
    rng: np.random.Generator, shift: float, beta: float, k: int
) -> tuple[pd.DataFrame, list[str]]:
    """Detection ~ logistic(beta (a - shift [in set] + 0.5)), a = abundance offset; only detected
    genes are ranked; the gene-level set holds every gene of its k TUs."""
    tu, x = _world(rng, N_DET, 0.8, 0.0)
    st = rng.choice(N_DET, k, replace=False)
    ins = np.isin(tu, st)
    a = x - 20.0
    pdet = 1.0 / (1.0 + np.exp(-beta * (a - np.where(ins, shift, 0.0) + 0.5)))
    det = rng.random(tu.size) < pdet
    return _frame(tu, x, det), [f"g{int(i) + 1}" for i in np.flatnonzero(ins)]


Draw = Callable[[np.random.Generator], tuple[float, float, int]]


def _det_params(kind: str) -> Draw:
    def draw(rng: np.random.Generator) -> tuple[float, float, int]:
        if kind == "mix":
            return (
                float(rng.uniform(0, 1.5)),
                float(rng.choice([1.0, 2.0])),
                int(rng.integers(15, 121)),
            )
        if kind == "null_b2":
            return 0.0, 2.0, 60
        if kind == "null_b1_mix":
            return 0.0, 1.0, int(rng.integers(15, 121))
        if kind == "s1_b1":
            return 1.0, 1.0, 60
        if kind == "s1.5_b2":
            return 1.5, 2.0, 60
        raise ValueError(kind)

    return draw


def _det_row(rng: np.random.Generator, kind: str) -> pd.Series:
    """The first replicate world whose set has a row (size_tus in [10, 500])."""
    draw = _det_params(kind)
    for _ in range(MAX_DRAWS):
        s, b, k = draw(rng)
        genes, members = _det_world(rng, s, b, k)
        r = _row(genes, members, rng)
        if r is not None:
            return r
    raise AssertionError(f"no row in {MAX_DRAWS} draws")


def _det_unflagged(kind: str) -> Callable[[np.random.Generator], float]:
    def rep(rng: np.random.Generator) -> float:
        for _ in range(MAX_DRAWS):
            r = _det_row(rng, kind)
            if "detection_confounded" not in as_list(r["flags"]):
                return _p(r)
        raise AssertionError(f"no unflagged row in {MAX_DRAWS} draws")

    return rep


DET_LEVEL = ("mix", "null_b2", "null_b1_mix", "s1_b1")
SCREEN = 1.5  # statistician #174: the detection flag is a screen, unflagged rows <= 1.5 alpha


@pytest.mark.parametrize("full", FAST_FULL)
@pytest.mark.parametrize("arm", DET_LEVEL)
def test_detection_confounded_unflagged_rows_within_the_screen(arm: str, full: bool) -> None:
    """Unflagged rows reject at <= 1.5 alpha (the §8.4 bound around 1.5 alpha), not at alpha: at
    R 4,000 the mix gives 0.0138 / 0.0592 / 0.1125 and the pooled mix (R about 12,600) 0.0128 /
    0.0567 / 0.1073, above the level bounds (Slurm 1867332)."""
    levels = C.load_tolerances().levels
    reps = C.replicates(SUITE, full=full)
    unflagged = _det_unflagged(arm)

    def rejects(rng: np.random.Generator) -> np.ndarray:
        p = unflagged(rng)
        return np.array([float(p <= a) for a in levels])

    rec = calibrate_levels(SUITE, "fwer", [SCREEN * a for a in levels], rejects, reps,
                           stream=f"{SUITE}.attack.flags.detection.{arm}")  # fmt: skip
    rec.assert_ok()


@pytest.mark.parametrize("full", FAST_FULL)
def test_detection_screen_teeth_all_rows_of_the_mix_reject_at_3_alpha(full: bool) -> None:
    """Without the screen, rows of the confounded mix reject at about 0.25 at alpha 0.05: >= 3 alpha
    (coverage convention above), so the screen visibly removes rejections."""
    alpha = 0.05

    def rep(rng: np.random.Generator) -> bool:
        return _p(_det_row(rng, "mix")) <= alpha

    hits = C.run_replicates(f"{SUITE}.attack.flags.detection_teeth.all_rows_mix", rep,
                            C.replicates(SUITE, full=full))  # fmt: skip
    ok, msg = _coverage_ok(hits, 3 * alpha)
    assert ok, msg


@pytest.mark.parametrize("full", FAST_FULL)
@pytest.mark.parametrize(("kind", "claim"), [("s1.5_b2", 0.99), ("s1_b1", 0.75)])
def test_detection_confounded_teeth_coverage(kind: str, claim: float, full: bool) -> None:
    def rep(rng: np.random.Generator) -> bool:
        return "detection_confounded" in as_list(_det_row(rng, kind)["flags"])

    hits = C.run_replicates(f"{SUITE}.attack.flags.detection_teeth.{kind}", rep,
                            C.replicates(SUITE, full=full))  # fmt: skip
    ok, msg = _coverage_ok(hits, claim)
    assert ok, msg


@pytest.mark.parametrize("full", FAST_FULL)
@pytest.mark.parametrize("arm", ["null_b2", "null_b1_mix"])
def test_detection_confounded_null_flag_rate(arm: str, full: bool) -> None:
    """Unconfounded null: the share of flagged rows is <= 0.25 (probe 0.133-0.160)."""

    def rep(rng: np.random.Generator) -> float:
        return float("detection_confounded" in as_list(_det_row(rng, arm)["flags"]))

    rec = C.calibrate(SUITE, "fwer", 0.25, rep, C.replicates(SUITE, full=full),
                      stream=f"{SUITE}.attack.flags.detection_null_rate.{arm}")  # fmt: skip
    assert rec["status"] != "fail", rec
