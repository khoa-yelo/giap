"""M8-05b properties of the GSEA row flags (M8-05 spec, decision 9 "Property tests").

* Permuting gene rows, or relabelling the TU scores, leaves delta (and the detection z) unchanged.
* delta is invariant to rotating all positions on a circular replicon.
* For a random placement the mean of C equals E[C], i.e. delta has mean 0: exactly, by enumerating
  every placement on small genomes, and by Monte Carlo on a 300-TU ring.
* The kernel matches the brute-force reference (every pair of ranked TUs) on random genomes.
* The detection z is the hypergeometric z of the counts, and the rule is sign-symmetric.

Interface as pinned in tests/unit/function/test_gsea_flags.py.
"""

from __future__ import annotations

import itertools
import math
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.gsea_data import Gene
from helpers.gsea_flags_oracle import delta_bruteforce, flags_frame, z_counts
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

SETTINGS = settings(max_examples=40, deadline=None, suppress_health_check=list(HealthCheck))


def _g() -> Any:
    from giap.function import gsea

    return gsea


# --------------------------------------------------------------------------- strategies


@st.composite
def genomes(draw: st.DrawFn, max_n: int = 60) -> dict[str, Any]:
    """1-3 replicons with ranked counts and topologies, and members on each."""
    n_rep = draw(st.integers(1, 3))
    n_ranked, circular, positions = {}, {}, {}
    for r in range(n_rep):
        name = f"r{r}"
        n = draw(st.integers(1, max_n))
        n_ranked[name] = n
        circular[name] = draw(st.booleans())
        positions[name] = sorted(draw(st.sets(st.integers(0, n - 1), max_size=n)))
    return {"positions": positions, "n_ranked": n_ranked, "circular": circular}


@st.composite
def gene_worlds(draw: st.DrawFn) -> tuple[pd.DataFrame, dict[str, str], dict[str, list[str]]]:
    """A gene frame of 10-50 TUs (1-3 genes, some undetected or partial) on 1-2 contigs, the
    contig topologies and 1-3 gene sets."""
    seed = draw(st.integers(0, 2**31 - 1))
    rng = np.random.default_rng(seed)
    n_tu = int(rng.integers(10, 51))
    n_chr = int(rng.integers(5, n_tu + 1)) if rng.random() < 0.5 else n_tu
    tus: dict[str, list[Gene]] = {}
    contig: dict[str, str] = {}
    k = 0
    for i in range(n_tu):
        genes = []
        for _ in range(int(rng.integers(1, 4))):
            k += 1
            det = bool(rng.random() < 0.7)
            genes.append(
                Gene(
                    f"g{k}",
                    float(rng.normal(20, 2)) + 1e-6 * k if det else None,
                    detected=det,
                    partial=bool(rng.random() < 0.05),
                )
            )
        tus[f"TU_{i + 1}"] = genes
        contig[f"TU_{i + 1}"] = "chr" if i < n_chr else "pl"
    df = flags_frame(tus, contig)
    topo = {"chr": str(rng.choice(["circular", "linear"])), "pl": "circular"}
    tags = list(df["locus_tag"])
    sets = {
        f"S{j}": [str(tags[int(x)]) for x in rng.choice(len(tags), int(rng.integers(1, len(tags))),
                                                         replace=False)]
        for j in range(int(rng.integers(1, 4)))
    }  # fmt: skip
    return df, topo, sets


def _deltas(
    genes: pd.DataFrame, topo: dict[str, str], sets: dict[str, list[str]]
) -> dict[str, float]:
    g = _g()
    rank = g.tu_ranking(genes)
    tsets = g.tu_sets(rank, _coll(sets))
    return {t: float(v) for t, v in g.cluster_deltas(rank, genes, topo, tsets).items()}


def _coll(sets: dict[str, list[str]]) -> Any:
    from helpers.gsea_data import gene_sets

    return gene_sets("kegg", sets)


# --------------------------------------------------------------------------- genomic_cluster


@SETTINGS
@given(genomes())
def test_kernel_matches_bruteforce(gen: dict[str, Any]) -> None:
    got = float(_g().cluster_delta(gen["positions"], gen["n_ranked"], gen["circular"]))
    want = delta_bruteforce(gen["positions"], gen["n_ranked"], gen["circular"])
    assert got == pytest.approx(want, rel=1e-9, abs=1e-12)


@SETTINGS
@given(gene_worlds(), st.integers(0, 2**31 - 1))
def test_permuting_gene_rows_leaves_delta_and_z_unchanged(world: Any, perm_seed: int) -> None:
    genes, topo, sets = world
    shuffled = genes.sample(frac=1.0, random_state=perm_seed % (2**32)).reset_index(drop=True)
    assert _deltas(genes, topo, sets) == _deltas(shuffled, topo, sets)
    g = _g()
    for members in sets.values():
        assert float(g.detection_z(genes, members)) == float(g.detection_z(shuffled, members))


@SETTINGS
@given(gene_worlds(), st.integers(0, 2**31 - 1))
def test_relabelling_tu_scores_leaves_delta_and_z_unchanged(world: Any, perm_seed: int) -> None:
    """Shuffle the abundance values among the scoring genes: the scores move between TUs, the
    set of ranked TUs and their order do not, so delta and z are unchanged."""
    genes, topo, sets = world
    g = _g()
    rank = g.tu_ranking(genes)
    mask = genes["locus_tag"].isin([x for v in rank.tu_genes.values() for x in v]).to_numpy()
    vals = genes.loc[mask, "abundance_log2"].to_numpy(float)
    moved = genes.copy()
    moved.loc[mask, "abundance_log2"] = np.random.default_rng(perm_seed).permutation(vals)
    assert _deltas(genes, topo, sets) == pytest.approx(_deltas(moved, topo, sets), abs=1e-12)
    for members in sets.values():
        assert float(g.detection_z(genes, members)) == float(g.detection_z(moved, members))


@SETTINGS
@given(genomes(max_n=80), st.integers(0, 10_000))
def test_delta_is_invariant_to_rotation_on_a_circular_replicon(gen: dict[str, Any], s: int) -> None:
    circular = {r: True for r in gen["n_ranked"]}
    rotated = {
        r: sorted((p + s) % gen["n_ranked"][r] for p in pos) for r, pos in gen["positions"].items()
    }
    a = float(_g().cluster_delta(gen["positions"], gen["n_ranked"], circular))
    b = float(_g().cluster_delta(rotated, gen["n_ranked"], circular))
    assert a == pytest.approx(b, rel=1e-12, abs=1e-12)


@settings(max_examples=30, deadline=None, suppress_health_check=list(HealthCheck))
@given(
    st.lists(st.tuples(st.integers(1, 6), st.booleans()), min_size=1, max_size=2),
    st.integers(2, 6),
)
def test_mean_delta_over_all_placements_is_exactly_zero(
    reps: list[tuple[int, bool]], k: int
) -> None:
    """Enumerate every placement of k members among the n ranked TUs: mean C = E[C] exactly."""
    n_ranked = {f"r{i}": n for i, (n, _) in enumerate(reps)}
    circular = {f"r{i}": c for i, (_, c) in enumerate(reps)}
    slots = [(r, i) for r, n in n_ranked.items() for i in range(n)]
    n = len(slots)
    if not 2 <= k < n:
        return
    total, count = 0.0, 0
    for combo in itertools.combinations(slots, k):
        pos: dict[str, list[int]] = {}
        for r, i in combo:
            pos.setdefault(r, []).append(i)
        total += float(_g().cluster_delta(pos, n_ranked, circular))
        count += 1
    assert total / count == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize(("k", "circular"), [(5, True), (30, True), (30, False), (150, True)])
def test_mean_c_equals_e_c_monte_carlo(k: int, circular: bool) -> None:
    """Random placements on a 300-TU replicon: the mean of delta is 0 within 4 standard errors."""
    n, draws = 300, 3000
    rng = np.random.default_rng([20261021, k, int(circular)])
    d = np.array([
        float(_g().cluster_delta({"r": sorted(rng.choice(n, k, replace=False).tolist())},
                                 {"r": n}, {"r": circular}))
        for _ in range(draws)
    ])  # fmt: skip
    se = d.std(ddof=1) / math.sqrt(draws)
    assert abs(d.mean()) <= 4 * se + 1e-12, (d.mean(), se)


# --------------------------------------------------------------------------- detection z


@SETTINGS
@given(gene_worlds())
def test_detection_z_is_the_hypergeometric_z_of_the_counts(world: Any) -> None:
    genes, _, sets = world
    elig = ~genes["partial"].astype(bool) & np.isfinite(genes["p_detect"].to_numpy(float))
    score = elig & genes["detected"].astype(bool) & genes["abundance_log2"].notna()
    M, N = int(elig.sum()), int(score.sum())
    for members in sets.values():
        ins = genes["locus_tag"].isin(members)
        n, x = int((elig & ins).sum()), int((score & ins).sum())
        assert float(_g().detection_z(genes, members)) == pytest.approx(
            z_counts(M, N, n, x), rel=1e-9, abs=1e-12
        )


@SETTINGS
@given(
    st.floats(-6, 6, allow_nan=False),
    st.one_of(st.floats(-1, 1, allow_nan=False), st.just(math.nan)),
)
def test_detection_rule_is_sign_symmetric(z: float, es: float) -> None:
    g = _g()
    assert g.is_detection_confounded(z, es) is g.is_detection_confounded(-z, -es)
    if abs(z) < 1.0:
        assert g.is_detection_confounded(z, es) is False
