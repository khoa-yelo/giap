"""M8-05b: the GSEA row flags ``genomic_cluster`` and ``detection_confounded``
(docs/specs/M8-05-gsea-tu-collapsed.md, lead decision 9; ADR-0025, ADR-0026).

Cases G1-G15 and D1-D14 of decision 9, the flag order of the ``flags`` cell, the new
``gsea_table`` argument ``replicons`` and the degenerate cases.

Interface pinned here (decision 9 names the statistics, not the functions):

* ``giap.function.gsea.CLUSTER_TAU`` = 0.15 and ``DETECTION_Z0`` = 1.0.
* ``cluster_delta(positions, n_ranked, circular) -> float``: the decision-9a kernel. ``positions``
  maps a replicon to its member positions (0-based indices in the ranked order of that replicon);
  ``n_ranked`` maps *every* replicon with ranked TUs to its count n_r (a replicon without members
  still enters n and P_w); ``circular`` maps a replicon to True when circular (missing = linear).
* ``cluster_deltas(rank, genes, replicons, sets) -> {term: delta}``: the same statistic from a
  ``TuRanking``, the gene frame (``contig`` per gene) and per-term member TUs (``tu_sets``).
* ``detection_z(genes, members) -> float``: the decision-9b z for one gene-level set.
* ``is_genomic_cluster(delta) -> bool`` (delta >= tau) and ``is_detection_confounded(z, es)``.
* ``gsea_table(genes, colls, size, cfg, alpha, rng, replicons=None)``: ``replicons`` is a
  ``contig -> topology`` mapping or the C1 frame (``contig``, ``topology``). Its default ``None``
  keeps the M8-05 six-argument calls valid and makes every contig linear; a contig missing from
  ``replicons`` is linear (decision 9a: any topology other than "circular" is linear); a gene frame
  without a ``contig`` column is one linear replicon.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.gsea_data import Gene, as_list, frame, gene_sets
from helpers.gsea_flags_oracle import (
    FLAG_ORDER,
    delta_bruteforce,
    flags_frame,
    replicons_frame,
    ring_frame,
    z_counts,
)

from giap.config.model import GseaConfig
from giap.rng import for_step
from giap.stats.fgsea import SimpleResult

DELTA_TOL = 5e-7  # "delta to 6 decimals"
Z_TOL = 5e-5  # "z to 4 decimals"
STREAM = "function.gsea"


def _g() -> Any:
    from giap.function import gsea

    return gsea


def _kernel(
    positions: Mapping[str, Sequence[int]],
    n_ranked: Mapping[str, int],
    circular: Mapping[str, bool],
) -> float:
    return float(_g().cluster_delta(positions, n_ranked, circular))


def _one(pos: Sequence[int], n: int, circular: bool = True) -> float:
    return _kernel({"r": list(pos)}, {"r": n}, {"r": circular})


def _run(
    genes: pd.DataFrame,
    colls: Sequence[Any],
    replicons: Any = None,
    size: tuple[int, int] = (10, 500),
    n_perm: int = 200,
    seed: int = 1,
    **kw: Any,
) -> pd.DataFrame:
    cfg = GseaConfig(n_perm_simple=n_perm)
    return _g().gsea_table(
        genes, colls, size, cfg, 0.05, for_step(STREAM, seed), replicons=replicons, **kw
    )


def _row(out: pd.DataFrame, term: str) -> pd.Series:
    rows = out[out["term"] == term]
    assert len(rows) == 1, f"expected one row for {term}, got {len(rows)}"
    return rows.iloc[0]


def _flags(r: pd.Series) -> list[str]:
    return as_list(r["flags"])


# --------------------------------------------------------------------------- constants


def test_thresholds_are_the_decision_9_values() -> None:
    g = _g()
    assert g.CLUSTER_TAU == 0.15
    assert g.DETECTION_Z0 == 1.0


# --------------------------------------------------------------------------- genomic_cluster


G_CASES = [
    # id, positions, n_ranked, circular, delta, flag
    ("G1", {"r": range(15)}, {"r": 100}, {"r": True}, 1.706342, True),
    ("G2", {"r": range(0, 89, 11)}, {"r": 100}, {"r": True}, -0.177427, False),
    ("G3", {"r": [97, 98, 99, 0, 1, 2]}, {"r": 100}, {"r": True}, 1.322169, True),
    ("G4", {"r": [97, 98, 99, 0, 1, 2]}, {"r": 100}, {"r": False}, 0.781309, True),
    ("G5", {"A": [59, 20], "B": [0, 20]}, {"A": 60, "B": 40}, {"A": True, "B": True}, -0.063070,
     False),
    ("G6", {"r": [59, 60, 20, 80]}, {"r": 100}, {"r": True}, 0.197347, True),
    ("G7", {"r": range(0, 29, 2)}, {"r": 100}, {"r": True}, 0.381740, True),
    ("G8", {"r": [s + i for s in (0, 300, 600, 900, 1200) for i in range(3)]}, {"r": 1500},
     {"r": True}, 0.822901, True),
    ("G9", {"r": [0, 1, 100, 101, 200, 201, *range(400, 1201, 100)]}, {"r": 1500}, {"r": True},
     0.183171, True),
    ("G10", {"r": [0, 1, 100, 101, *range(300, 1301, 100)]}, {"r": 1500}, {"r": True}, 0.115831,
     False),
    ("G11", {"r": range(0, 1401, 100)}, {"r": 1500}, {"r": True}, -0.018849, False),
    ("G12", {"r": [0, 1]}, {"r": 2, "s": 1}, {"r": True, "s": True}, 1.0, True),
]  # fmt: skip


@pytest.mark.parametrize(
    ("case", "positions", "n_ranked", "circular", "delta", "flag"),
    G_CASES,
    ids=[c[0] for c in G_CASES],
)
def test_cluster_delta_decision_9_cases(
    case: str,
    positions: Mapping[str, Sequence[int]],
    n_ranked: Mapping[str, int],
    circular: Mapping[str, bool],
    delta: float,
    flag: bool,
) -> None:
    pos = {r: list(p) for r, p in positions.items()}
    got = _kernel(pos, n_ranked, circular)
    assert got == pytest.approx(delta, abs=DELTA_TOL), case
    assert got == pytest.approx(delta_bruteforce(pos, n_ranked, circular), abs=1e-9)
    assert _g().is_genomic_cluster(got) is flag


def test_g4_the_seam_pair_does_not_count_on_a_linear_replicon() -> None:
    circ, lin = _one([97, 98, 99, 0, 1, 2], 100, True), _one([97, 98, 99, 0, 1, 2], 100, False)
    assert circ > lin
    assert lin == pytest.approx(_one([0, 1, 2, 50, 51, 52], 100, False), abs=1e-12)


def test_g5_members_on_different_replicons_are_never_paired() -> None:
    # A59 and B0 would be adjacent if the replicons were concatenated
    split = _kernel({"A": [59, 20], "B": [0, 20]}, {"A": 60, "B": 40}, {"A": True, "B": True})
    joined = _one([59, 20, 60, 80], 100, True)  # G6
    assert split < 0 < joined


@pytest.mark.parametrize(
    ("positions", "n"),
    [([], 10), ([3], 10), ([0], 1), (list(range(10)), 10), ([0, 1], 2)],
    ids=["k0", "k1", "k1_n1", "k_eq_n", "k_eq_n_2"],
)
def test_g13_delta_is_zero_when_k_below_2_or_k_at_least_n(positions: list[int], n: int) -> None:
    for circular in (True, False):
        assert _one(positions, n, circular) == 0.0


def test_g13_k_equal_n_over_several_replicons() -> None:
    got = _kernel({"a": [0, 1], "b": [0, 1, 2]}, {"a": 2, "b": 3}, {"a": True, "b": False})
    assert got == 0.0


def test_replicon_without_members_still_counts_in_n() -> None:
    alone = _one([0, 1, 2], 50, True)
    with_plasmid = _kernel({"r": [0, 1, 2]}, {"r": 50, "p": 30}, {"r": True, "p": True})
    assert with_plasmid != pytest.approx(alone)
    assert with_plasmid == pytest.approx(
        delta_bruteforce({"r": [0, 1, 2]}, {"r": 50, "p": 30}, {"r": True, "p": True}), abs=1e-12
    )


def test_missing_topology_is_linear() -> None:
    assert _kernel({"r": [0, 99]}, {"r": 100}, {}) == pytest.approx(_one([0, 99], 100, False))


@pytest.mark.parametrize("n", [3, 4, 5, 20, 21])
def test_circular_closed_form_at_half_ring(n: int) -> None:
    """c_d = n_r/2 when 2d = n_r (even rings), checked by brute force on small rings."""
    pos = [0, n // 2]
    assert _one(pos, n, True) == pytest.approx(
        delta_bruteforce({"r": pos}, {"r": n}, {"r": True}), abs=1e-12
    )


def test_weights_stop_beyond_lag_10() -> None:
    assert _one([0, 10], 200, False) > _one([0, 11], 200, False)
    assert _one([0, 11], 200, False) == pytest.approx(_one([0, 50], 200, False), abs=1e-12)


@pytest.mark.parametrize(
    ("delta", "flag"), [(0.15, True), (0.1499999, False), (2.0, True), (-0.3, False), (0.0, False)]
)
def test_cluster_threshold_is_inclusive(delta: float, flag: bool) -> None:
    assert _g().is_genomic_cluster(delta) is flag


# ------------------------------------------------------------- cluster_deltas (TuRanking level)


def _deltas(
    genes: pd.DataFrame, replicons: Any, sets: Mapping[str, Sequence[str]]
) -> dict[str, float]:
    g = _g()
    rank = g.tu_ranking(genes)
    out = g.cluster_deltas(rank, genes, replicons, {t: frozenset(s) for t, s in sets.items()})
    return {t: float(v) for t, v in out.items()}


def test_cluster_deltas_match_the_kernel_on_a_ring() -> None:
    genes = ring_frame(100)
    sets = {"run": [f"TU_{i + 1}" for i in range(15)], "seam": ["TU_98", "TU_99", "TU_100",
            "TU_1", "TU_2", "TU_3"]}  # fmt: skip
    got = _deltas(genes, {"chr": "circular"}, sets)
    assert got["run"] == pytest.approx(1.706342, abs=DELTA_TOL)  # G1
    assert got["seam"] == pytest.approx(1.322169, abs=DELTA_TOL)  # G3
    lin = _deltas(genes, {"chr": "linear"}, sets)
    assert lin["seam"] == pytest.approx(0.781309, abs=DELTA_TOL)  # G4


def test_g14_unranked_tus_are_skipped_in_the_order() -> None:
    """TU_1..TU_12, only the odd ones ranked: {TU_1, TU_3, TU_5} is a run (d = 1)."""
    tus = {}
    for i in range(1, 13):
        if i % 2:
            tus[f"TU_{i}"] = [Gene(f"g{i}", 20.0 + 0.1 * i)]
        else:
            tus[f"TU_{i}"] = [Gene(f"g{i}", None, detected=False)]
    genes = flags_frame(tus)
    got = _deltas(genes, {"chr": "linear"}, {"S": ["TU_1", "TU_3", "TU_5"]})["S"]
    assert got == pytest.approx(_one([0, 1, 2], 6, False), abs=1e-12)
    assert got == pytest.approx(0.591667, abs=1e-6)
    assert got != pytest.approx(_one([0, 2, 4], 12, False), abs=1e-3)  # not the unranked order


def test_g14_a_tu_without_a_scoring_gene_is_skipped() -> None:
    """A TU with only a partial or NA-abundance gene is not ranked, so it is not a gap either."""
    tus = {
        "TU_1": [Gene("g1", 20.1)],
        "TU_2": [Gene("g2", 20.2, partial=True)],
        "TU_3": [Gene("g3", None)],
        "TU_4": [Gene("g4", 20.4)],
        **{f"TU_{i}": [Gene(f"g{i}", 20.0 + 0.1 * i)] for i in range(5, 31)},
    }
    got = _deltas(flags_frame(tus), {"chr": "linear"}, {"S": ["TU_1", "TU_4"]})["S"]
    assert got == pytest.approx(_one([0, 1], 28, False), abs=1e-12)


def test_g15_natural_sort_of_tu_ids() -> None:
    genes = ring_frame(12)
    d = _deltas(genes, {"chr": "linear"}, {"a": ["TU_9", "TU_10"], "b": ["TU_1", "TU_10"]})
    assert d["a"] == pytest.approx(_one([8, 9], 12, False), abs=1e-12)  # adjacent
    assert d["b"] == pytest.approx(_one([0, 9], 12, False), abs=1e-12)  # d = 9, not 1
    assert d["a"] > d["b"]


def test_ranked_order_does_not_follow_gene_row_order() -> None:
    genes = ring_frame(40)
    shuffled = genes.sample(frac=1.0, random_state=3).reset_index(drop=True)
    sets = {"S": [f"TU_{i}" for i in (5, 6, 7, 30)]}
    assert _deltas(genes, {"chr": "circular"}, sets) == _deltas(shuffled, {"chr": "circular"}, sets)


def test_origin_crossing_tu_is_adjacent_to_the_first_tu_on_a_circular_replicon() -> None:
    """build_tus numbers an origin-crossing TU last on its replicon (TU_20 here); the circular
    distance makes it adjacent to TU_1."""
    tus = {f"TU_{i}": [Gene(f"g{i}", 20.0 + 0.1 * i)] for i in range(1, 20)}
    tus["TU_20"] = [Gene("g20a", 22.0), Gene("g20b", 22.1)]  # spans the origin
    genes = flags_frame(tus)
    circ = _deltas(genes, {"chr": "circular"}, {"S": ["TU_20", "TU_1", "TU_2"]})["S"]
    lin = _deltas(genes, {"chr": "linear"}, {"S": ["TU_20", "TU_1", "TU_2"]})["S"]
    assert circ == pytest.approx(_one([19, 0, 1], 20, True), abs=1e-12)
    assert lin == pytest.approx(_one([19, 0, 1], 20, False), abs=1e-12)
    assert circ > lin


def test_multiple_replicons_split_by_contig_never_paired_across() -> None:
    """TU_60 (last on chr) and TU_61 (first on the plasmid) are consecutive IDs, not neighbours."""
    tus = {f"TU_{i}": [Gene(f"g{i}", 20.0 + 0.01 * i)] for i in range(1, 101)}
    contig = {f"TU_{i}": ("chr" if i <= 60 else "pl") for i in range(1, 101)}
    genes = flags_frame(tus, contig)
    reps = {"chr": "circular", "pl": "circular"}
    d = _deltas(genes, reps, {"x": ["TU_60", "TU_61"], "g5": ["TU_60", "TU_21", "TU_61", "TU_81"]})
    assert d["x"] == pytest.approx(
        _kernel({"chr": [59], "pl": [0]}, {"chr": 60, "pl": 40}, {"chr": True, "pl": True}),
        abs=1e-12,
    )
    assert d["x"] < 0.0
    assert d["g5"] == pytest.approx(-0.063070, abs=DELTA_TOL)  # G5 through the frame


def test_replicons_mapping_and_c1_frame_agree() -> None:
    tus = {f"TU_{i}": [Gene(f"g{i}", 20.0 + 0.01 * i)] for i in range(1, 101)}
    contig = {f"TU_{i}": ("chr" if i <= 60 else "pl") for i in range(1, 101)}
    genes = flags_frame(tus, contig)
    sets = {"a": ["TU_60", "TU_1", "TU_2"], "b": ["TU_100", "TU_61"]}
    m = _deltas(genes, {"chr": "circular", "pl": "linear"}, sets)
    f = _deltas(genes, replicons_frame({"chr": "circular", "pl": "linear"}), sets)
    assert m == f
    assert m["b"] == pytest.approx(
        _kernel({"pl": [39, 0]}, {"chr": 60, "pl": 40}, {"chr": True, "pl": False}), abs=1e-12
    )


def test_unknown_topology_or_missing_contig_is_linear() -> None:
    genes = ring_frame(30)
    sets = {"S": ["TU_30", "TU_1"]}
    lin = _deltas(genes, {"chr": "linear"}, sets)
    assert _deltas(genes, {"chr": "unknown"}, sets) == lin
    assert _deltas(genes, {}, sets) == lin
    assert _deltas(genes, None, sets) == lin


def test_gene_frame_without_contig_is_one_linear_replicon() -> None:
    genes = ring_frame(30).drop(columns="contig")
    sets = {"S": ["TU_30", "TU_1", "TU_2"]}
    assert _deltas(genes, None, sets)["S"] == pytest.approx(_one([29, 0, 1], 30, False), abs=1e-12)


# --------------------------------------------------------------------------- detection z


def _count_frame(M: int, N: int, n: int, x: int) -> tuple[pd.DataFrame, list[str]]:
    """M eligible genes (one per TU): the first N score, the rest are undetected. The set holds
    x scoring genes and n - x undetected ones."""
    assert 0 <= x <= min(n, N) and n - x <= M - N
    tus = {}
    for i in range(M):
        if i < N:
            tus[f"TU_{i + 1}"] = [Gene(f"g{i + 1}", 20.0 + 0.001 * i)]
        else:
            tus[f"TU_{i + 1}"] = [Gene(f"g{i + 1}", None, detected=False)]
    members = [f"g{i + 1}" for i in range(x)] + [f"g{N + i + 1}" for i in range(n - x)]
    return flags_frame(tus), members


D_CASES = [
    # id, M, N, n, x, ES, z, flag
    ("D1", 1000, 500, 40, 12, 0.4, -2.5807, True),
    ("D2", 1000, 500, 40, 12, -0.4, -2.5807, False),
    ("D3", 1000, 500, 40, 18, 0.4, -0.6452, False),
    ("D4", 1000, 500, 40, 27, -0.3, 2.2581, True),
    ("D5", 1000, 500, 40, 27, 0.3, 2.2581, False),
    ("D6", 1000, 600, 30, 15, 0.5, -1.1346, True),
    ("D7", 1000, 600, 30, 16, 0.5, -0.7564, False),
    ("D8", 1000, 500, 40, 12, math.nan, -2.5807, True),
    ("D9", 1000, 1000, 40, 40, 0.4, 0.0, False),
    ("D10", 1000, 500, 1000, 500, 0.4, 0.0, False),
]


@pytest.mark.parametrize(
    ("case", "M", "N", "n", "x", "es", "z", "flag"), D_CASES, ids=[c[0] for c in D_CASES]
)
def test_detection_z_decision_9_cases(
    case: str, M: int, N: int, n: int, x: int, es: float, z: float, flag: bool
) -> None:
    genes, members = _count_frame(M, N, n, x)
    got = float(_g().detection_z(genes, members))
    assert got == pytest.approx(z, abs=Z_TOL), case
    assert got == pytest.approx(z_counts(M, N, n, x), abs=1e-12)
    assert _g().is_detection_confounded(got, es) is flag


@pytest.mark.parametrize(
    ("z", "es", "flag"),
    [
        (-1.0, 0.2, True),  # |z| >= 1 is inclusive
        (1.0, -0.2, True),
        (-0.9999, 0.2, False),
        (-3.0, -0.2, False),  # the bias works against the observed sign
        (3.0, 0.2, False),
        (-1.0, math.nan, True),  # ES NaN: |z| rule
        (1.0, math.nan, True),
        (0.9999, math.nan, False),
        (2.0, 0.0, True),  # ES 0: |z| rule
        (-2.0, 0.0, True),
        (-0.5, 0.0, False),
        (0.0, 0.0, False),
        (0.0, math.nan, False),
    ],
)
def test_detection_rule(z: float, es: float, flag: bool) -> None:
    assert _g().is_detection_confounded(z, es) is flag


def test_zero_variance_gives_z_zero() -> None:
    g = _g()
    # f = 1 (every eligible gene scores)
    genes, members = _count_frame(50, 50, 10, 10)
    assert g.detection_z(genes, members) == 0.0
    # n = 0: no eligible gene of the set is in the frame
    assert g.detection_z(genes, ["nope1", "nope2"]) == 0.0
    # M = 1: one eligible gene (no 0/0)
    one = flags_frame({"TU_1": [Gene("g1", 20.0)], "TU_2": [Gene("g2", 21.0, partial=True)]})
    assert g.detection_z(one, ["g1"]) == 0.0


def test_d11_partial_gene_counts_in_neither_n_nor_m() -> None:
    genes, members = _count_frame(1000, 500, 40, 12)
    extra = flags_frame({"TU_2001": [Gene("p1", 25.0, partial=True)],
                         "TU_2002": [Gene("p2", None, detected=False, partial=True)]})  # fmt: skip
    more = pd.concat([genes, extra], ignore_index=True)
    assert float(_g().detection_z(more, [*members, "p1", "p2"])) == pytest.approx(
        z_counts(1000, 500, 40, 12), abs=1e-12
    )


def test_d11_non_finite_p_detect_is_ineligible() -> None:
    genes, members = _count_frame(1000, 500, 40, 12)
    extra = flags_frame(
        {
            "TU_2001": [Gene("q1", 25.0, p_detect=math.nan)],
            "TU_2002": [Gene("q2", None, detected=False, p_detect=math.inf)],
        }
    )
    more = pd.concat([genes, extra], ignore_index=True)
    assert float(_g().detection_z(more, [*members, "q1", "q2"])) == pytest.approx(
        z_counts(1000, 500, 40, 12), abs=1e-12
    )


def test_d12_undetected_member_and_na_abundance_member_count_in_n_not_x() -> None:
    genes, members = _count_frame(1000, 500, 40, 12)
    extra = flags_frame({"TU_2001": [Gene("u1", None, detected=False)],
                         "TU_2002": [Gene("na1", None, detected=True)]})  # fmt: skip
    more = pd.concat([genes, extra], ignore_index=True)
    got = float(_g().detection_z(more, [*members, "u1", "na1"]))
    assert got == pytest.approx(z_counts(1002, 500, 42, 12), abs=1e-12)


def test_d13_an_undetected_gene_in_no_set_changes_m() -> None:
    genes, members = _count_frame(1000, 500, 40, 12)
    extra = flags_frame({"TU_2001": [Gene("bg", None, detected=False)]})
    more = pd.concat([genes, extra], ignore_index=True)
    got = float(_g().detection_z(more, members))
    assert got == pytest.approx(z_counts(1001, 500, 40, 12), abs=1e-12)
    assert got != pytest.approx(z_counts(1000, 500, 40, 12), abs=1e-9)


def test_detection_z_without_partial_or_p_detect_columns() -> None:
    genes, members = _count_frame(200, 100, 20, 5)
    bare = genes.drop(columns=["partial", "p_detect"])
    assert float(_g().detection_z(bare, members)) == pytest.approx(
        z_counts(200, 100, 20, 5), abs=1e-12
    )


def test_detection_z_counts_genes_not_tus() -> None:
    """A two-gene TU with one detected gene: both genes count in M and n; one in N and x."""
    tus = {f"TU_{i}": [Gene(f"g{i}", 20.0 + 0.01 * i)] for i in range(1, 51)}
    tus["TU_51"] = [Gene("a", 21.0), Gene("b", None, detected=False)]
    genes = flags_frame(tus)
    got = float(_g().detection_z(genes, ["a", "b", "g1", "g2"]))
    assert got == pytest.approx(z_counts(52, 51, 4, 3), abs=1e-12)


# --------------------------------------------------------------------------- gsea_table


class _Fake:
    """``fgsea_simple`` returning (p, ES) keyed by set size in TUs (NaN p = unestimated)."""

    def __init__(self, by_size: Mapping[int, tuple[float, float]]) -> None:
        self.by_size = dict(by_size)

    def __call__(self, prep: Any, n_perm: int, sample_size: int, rng: Any) -> SimpleResult:
        sizes = [int(s) for s in np.asarray(prep.sizes)]
        k = len(sizes)
        p = np.array([self.by_size.get(s, (0.5, 0.5))[0] for s in sizes], dtype=float)
        es = np.array([self.by_size.get(s, (0.5, 0.5))[1] for s in sizes], dtype=float)
        bad = np.isnan(p)
        edges = tuple(
            np.asarray(prep.order)[np.asarray(m)[:1] - 1].astype(np.int64) for m in prep.members
        )
        return SimpleResult(
            es,
            np.where(bad, np.nan, 1.5 * np.sign(np.nan_to_num(es))),
            p,
            np.where(bad, np.nan, 0.1),
            np.zeros(k, dtype=np.int64),
            np.where(bad, 3, 500).astype(np.int64),
            np.zeros(k, dtype=bool),
            edges,
        )


@pytest.fixture
def fake(monkeypatch: pytest.MonkeyPatch) -> Callable[[Mapping[int, tuple[float, float]]], _Fake]:
    def install(by_size: Mapping[int, tuple[float, float]]) -> _Fake:
        import giap.stats.fgsea as fg

        eng = _Fake(by_size)
        monkeypatch.setattr(fg, "fgsea_simple", eng)
        mod = _g()
        if hasattr(mod, "fgsea_simple"):
            monkeypatch.setattr(mod, "fgsea_simple", eng)
        return eng

    return install


def _ring_with_undetected(n: int = 200, n_undet: int = 20) -> pd.DataFrame:
    """``n`` singleton scoring TUs on a ring; TU_(n-n_undet+1)..TU_n also hold an undetected gene
    u1..u<n_undet> at tu_position 2 (eligible, not scoring)."""
    tus: dict[str, list[Gene]] = {}
    for i in range(1, n + 1):
        tus[f"TU_{i}"] = [Gene(f"g{i}", 20.0 + 0.01 * ((7 * i) % n))]
        j = i - (n - n_undet)
        if j > 0:
            tus[f"TU_{i}"].append(Gene(f"u{j}", None, detected=False))
    return flags_frame(tus)


UNDET = [f"u{j}" for j in range(1, 21)]


def test_flags_order_unestimated_row(fake: Any) -> None:
    """ES NaN (unestimated): the |z| rule; order ci_unestimated, unestimated, genomic_cluster,
    detection_confounded."""
    fake({15: (math.nan, math.nan)})
    genes = _ring_with_undetected()
    run15 = [f"g{i}" for i in range(101, 116)]
    out = _run(genes, [gene_sets("kegg", {"U": run15 + UNDET})], {"chr": "circular"})
    assert _flags(_row(out, "U")) == [
        "ci_unestimated",
        "unestimated",
        "genomic_cluster",
        "detection_confounded",
    ]


def test_flags_order_size_filtered_row_and_flags_from_the_main_run(fake: Any) -> None:
    """GO parent P = child C (TU_1..TU_10 plus 20 undetected genes) + TU_100 + TU_150. The
    re-test without C keeps 2 TUs (size_filtered), scattered and fully detected; P's new flags come
    from its main run (D14), so it still carries both."""
    fake({12: (0.001, 0.5), 10: (0.001, 0.5), 2: (0.9, 0.5)})
    genes = _ring_with_undetected()
    child = [f"g{i}" for i in range(1, 11)] + UNDET
    parent = [*child, "g100", "g150"]
    coll = gene_sets("go_bp", {"P": parent, "C": child}, parents={"C": ["P"]})
    out = _run(genes, [coll], {"chr": "circular"})
    assert _flags(_row(out, "P")) == [
        "ci_unestimated",
        "size_filtered",
        "genomic_cluster",
        "detection_confounded",
    ]
    assert _flags(_row(out, "C")) == ["ci_unestimated", "genomic_cluster", "detection_confounded"]
    assert len(out) == 2  # re-tests are not rows


def test_flags_follow_the_decision_9_order_in_every_row(fake: Any) -> None:
    fake({15: (math.nan, math.nan), 12: (0.001, 0.5), 10: (0.001, 0.5), 2: (0.9, 0.5),
          13: (0.2, -0.5)})  # fmt: skip
    genes = _ring_with_undetected()
    child = [f"g{i}" for i in range(1, 11)] + UNDET
    sets = {
        "P": [*child, "g100", "g150"],
        "C": child,
        "U": [f"g{i}" for i in range(101, 116)] + UNDET,
        "Z": [f"g{i}" for i in range(5, 200, 15)],
    }
    out = _run(genes, [gene_sets("go_bp", sets, parents={"C": ["P"]})], {"chr": "circular"})
    for _, r in out.iterrows():
        f = _flags(r)
        assert f == [x for x in FLAG_ORDER if x in f], (r["term"], f)
        assert len(f) == len(set(f))
        assert set(f) <= set(FLAG_ORDER)


def test_no_new_flags_for_a_scattered_fully_detected_set(fake: Any) -> None:
    fake({13: (0.2, 0.5)})
    genes = _ring_with_undetected()
    out = _run(genes, [gene_sets("kegg", {"Z": [f"g{i}" for i in range(5, 200, 15)]})],
               {"chr": "circular"})  # fmt: skip
    assert _flags(_row(out, "Z")) == ["ci_unestimated"]


@pytest.mark.parametrize(("es", "flag"), [(0.0, True), (0.5, False), (-0.5, True)])
def test_detection_flag_uses_the_row_es(fake: Any, es: float, flag: bool) -> None:
    """A set detected *more* often than the background (z > +1): flagged for ES < 0 and ES = 0."""
    fake({12: (0.2, es)})
    tus: dict[str, list[Gene]] = {}
    for i in range(1, 201):
        det = i % 2 == 1 or i > 180  # TUs 181..200 all scoring; half of the rest undetected
        tus[f"TU_{i}"] = [Gene(f"g{i}", 20.0 + 0.01 * i if det else None, detected=det)]
    genes = flags_frame(tus)
    members = [f"g{i}" for i in range(181, 200, 2)] + ["g1", "g101"]  # 12 TUs, all scoring
    z = z_counts(200, 110, 12, 12)
    assert z >= 1.0
    out = _run(genes, [gene_sets("kegg", {"S": members})], {"chr": "linear"})
    assert ("detection_confounded" in _flags(_row(out, "S"))) is flag


def test_detection_flag_with_the_real_engine_d1_d2() -> None:
    """D1 / D2 through fgsea: x = 12 scoring members at the top (ES > 0, flagged) or the bottom
    (ES < 0, not flagged) of the ranking, with 28 undetected members (z = -2.5807)."""
    M, N = 1000, 500
    pos = list(range(0, 480, 40))  # 12 scattered ranked positions: no genomic_cluster
    rest = [i for i in range(N) if i not in pos]

    def world(top: bool) -> tuple[pd.DataFrame, list[str]]:
        vals = np.empty(N)
        hi = 30.0 + 0.01 * np.arange(12)
        lo = 10.0 + 0.01 * np.arange(12)
        vals[pos] = hi if top else lo
        vals[rest] = 20.0 + 0.001 * np.arange(len(rest))
        tus: dict[str, list[Gene]] = {}
        for i in range(M):
            if i < N:
                tus[f"TU_{i + 1}"] = [Gene(f"g{i + 1}", float(vals[i]))]
            else:
                tus[f"TU_{i + 1}"] = [Gene(f"g{i + 1}", None, detected=False)]
        members = [f"g{i + 1}" for i in pos] + [f"g{N + j + 1}" for j in range(28)]
        return flags_frame(tus), members

    for top, flagged in ((True, True), (False, False)):
        genes, members = world(top)
        assert float(_g().detection_z(genes, members)) == pytest.approx(-2.5807, abs=Z_TOL)
        out = _run(genes, [gene_sets("kegg", {"S": members})], {"chr": "circular"})
        r = _row(out, "S")
        assert bool(r["es"] > 0) is top
        f = _flags(r)
        assert ("detection_confounded" in f) is flagged
        assert "genomic_cluster" not in f


@pytest.mark.parametrize(
    ("replicons", "flag"),
    [
        ({"chr": "circular"}, True),
        ({"chr": "linear"}, False),
        (None, False),
        ({}, False),
    ],
    ids=["circular", "linear", "default_none", "missing_contig"],
)
def test_gsea_table_replicons_argument(fake: Any, replicons: Any, flag: bool) -> None:
    """{TU_200, TU_1} on a 200-TU ring: delta 0.49 when circular, -0.01 when linear."""
    fake({2: (0.2, 0.5)})
    out = _run(
        ring_frame(200), [gene_sets("kegg", {"S": ["g200", "g1"]})], replicons, size=(2, 500)
    )
    assert ("genomic_cluster" in _flags(_row(out, "S"))) is flag


def test_gsea_table_replicons_as_c1_frame(fake: Any) -> None:
    fake({2: (0.2, 0.5)})
    sets = [gene_sets("kegg", {"S": ["g200", "g1"]})]
    a = _run(ring_frame(200), sets, replicons_frame({"chr": "circular"}), size=(2, 500))
    b = _run(ring_frame(200), sets, {"chr": "circular"}, size=(2, 500))
    assert _flags(_row(a, "S")) == _flags(_row(b, "S"))
    assert "genomic_cluster" in _flags(_row(a, "S"))


def test_gsea_table_two_replicons_are_not_joined(fake: Any) -> None:
    fake({2: (0.2, 0.5)})
    tus = {f"TU_{i}": [Gene(f"g{i}", 20.0 + 0.01 * i)] for i in range(1, 101)}
    contig = {f"TU_{i}": ("chr" if i <= 60 else "pl") for i in range(1, 101)}
    genes = flags_frame(tus, contig)
    reps = {"chr": "circular", "pl": "circular"}
    out = _run(genes, [gene_sets("kegg", {"X": ["g60", "g61"], "Y": ["g59", "g60"]})], reps,
               size=(2, 500))  # fmt: skip
    assert "genomic_cluster" not in _flags(_row(out, "X"))
    assert "genomic_cluster" in _flags(_row(out, "Y"))


def test_gsea_table_six_argument_call_still_works() -> None:
    """The M8-05 call without ``replicons`` (and a frame without ``contig``) stays valid."""
    rng = np.random.default_rng(0)
    vals = 20.0 + 0.01 * rng.permutation(60)
    genes = frame({f"TU_{i + 1}": [Gene(f"g{i + 1}", float(vals[i]))] for i in range(60)})
    cfg = GseaConfig(n_perm_simple=100)
    out = _g().gsea_table(genes, [gene_sets("kegg", {"S": [f"g{i}" for i in range(1, 13)]})],
                          (10, 500), cfg, 0.05, for_step(STREAM, 1))  # fmt: skip
    assert len(out) == 1
    assert "genomic_cluster" in _flags(out.iloc[0])  # a run of 12 on one linear replicon


def test_gsea_table_flags_need_no_rng(fake: Any) -> None:
    """Both flags are pure functions of the inputs: the same rows for any seed."""
    fake({15: (0.3, 0.4), 13: (0.3, -0.4)})
    genes = _ring_with_undetected()
    sets = [gene_sets("kegg", {"U": [f"g{i}" for i in range(101, 116)] + UNDET,
                               "Z": [f"g{i}" for i in range(5, 200, 15)]})]  # fmt: skip
    a = _run(genes, sets, {"chr": "circular"}, seed=1)
    b = _run(genes, sets, {"chr": "circular"}, seed=999)
    assert [_flags(r) for _, r in a.iterrows()] == [_flags(r) for _, r in b.iterrows()]
