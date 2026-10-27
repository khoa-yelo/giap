"""Calibration of gippasim layout knobs (issue M4-07; GIPPA_SPEC.md §11.3, §8.4).

Every distributional claim of docs/specs/M4-07-sim-layout-tables.md, pooled over ``N``
replicate 1,000-gene layouts (seeds ``calibration_base + i`` from tolerances.yaml). Proportions
must lie within 4 binomial SDs; distributions must pass KS / chi-square at ``ks_min_p``:

* TU sizes geometric with mean ``mean_tu_size`` (1.8); TU strand kept with ``p_same_strand`` (0.6);
* CDS lengths log-normal (median 280 aa, sigma 0.45);
* ``frac_sorf`` 0.02, ``frac_pseudogene`` 0.01, ``frac_hypothetical`` 0.25 of all genes;
* a pseudogene has ``aa = NA`` with probability 0.5;
* intra-TU pair p ~ Beta(8, 2), inter-TU adjacent same-strand pair p ~ Beta(2, 8) (and the knob
  is honoured for other parameters);
* prophage subsets with probabilities (0.05, 0.05, 0.4, 0.1, 0.4); ``core_biosynthetic`` for 30 %
  of BGC genes; xrefs for 70 % of genes;
* without ``short_genes_in`` prophage genes are not shorter than the rest.
"""

from __future__ import annotations

import itertools
import math
import sys
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import layout_oracle as lo

TOL = yaml.safe_load((Path(__file__).parent / "tolerances.yaml").read_text())
BASE = int(TOL["seeds"]["calibration_base"])
KS_MIN_P = float(TOL["ks_min_p"])
N = 20

pytestmark = pytest.mark.calibration


def _layout(seed: int, **kw: Any) -> Any:
    from gippasim.layout import LayoutConfig, simulate_layout

    return simulate_layout(LayoutConfig(**kw), np.random.default_rng(seed))


@cache
def _pool() -> tuple[Any, ...]:
    return tuple(_layout(BASE + i) for i in range(N))


def _genes(layout: Any) -> pd.DataFrame:
    return layout.genes


def _within(k: float, n: int, p: float) -> None:
    assert lo.within_sd(k, n, p), (k, n, p, k / n)


def _pair_ps(layout: Any) -> tuple[list[float], list[float]]:
    tmap = lo.tu_of(layout.tus)
    pairs = layout.operon_pairs
    col = lo.pair_prob_column(pairs)
    intra, inter = [], []
    for a, b, p in zip(pairs["gene_a"], pairs["gene_b"], pairs[col], strict=True):
        (intra if tmap[a] == tmap[b] else inter).append(float(p))
    return intra, inter


def test_tu_size_geometric_mean() -> None:
    sizes = np.array([len(lo.as_list(c)) for lay in _pool() for c in lay.tus["genes"]])
    m, p = 1.8, 1 / 1.8
    se = math.sqrt((1 - p) / p**2 / len(sizes))
    assert abs(sizes.mean() - m) <= 4 * se, (sizes.mean(), se)
    # geometric shape: P(size = 1) = p
    _within(int((sizes == 1).sum()), len(sizes), p)


def test_tu_strand_kept_with_p_same_strand() -> None:
    k = n = 0
    for lay in _pool():
        tus = lay.tus.sort_values(["contig", "start"], kind="mergesort")
        for _, sub in tus.groupby("contig"):
            s = sub["strand"].tolist()
            n += len(s) - 1
            k += sum(a == b for a, b in itertools.pairwise(s))
    _within(k, n, 0.6)


def test_cds_lengths_lognormal() -> None:
    g = pd.concat([_genes(lay) for lay in _pool()])
    cds = g[g["feature_type"] == "cds"]
    len_aa = (cds["end"] - cds["start"]).to_numpy() / 3 - 1
    logs = np.log(len_aa)
    p = stats.kstest(logs, stats.norm(math.log(280.0), 0.45).cdf).pvalue
    assert p >= KS_MIN_P, p
    se = 1.2533 * 0.45 / math.sqrt(len(logs))
    assert abs(np.median(logs) - math.log(280.0)) <= 4 * se


@pytest.mark.parametrize(
    ("what", "p"), [("sorf", 0.02), ("pseudogene", 0.01), ("hypothetical", 0.25)]
)
def test_gene_class_fractions(what: str, p: float) -> None:
    g = pd.concat([_genes(lay) for lay in _pool()])
    k = (
        int(g["hypothetical"].astype(bool).sum())
        if what == "hypothetical"
        else int((g["feature_type"] == what).sum())
    )
    _within(k, len(g), p)


def test_sorfs_short() -> None:
    g = pd.concat([_genes(lay) for lay in _pool()])
    s = g[g["feature_type"] == "sorf"]
    assert len(s) > 0
    assert ((s["end"] - s["start"]) // 3 - 1 < 30).all()


def test_pseudogene_aa_na_half() -> None:
    g = pd.concat([_genes(lay) for lay in _pool()])
    ps = g[g["feature_type"] == "pseudogene"]
    assert len(ps) >= 100
    _within(int(ps["aa"].isna().sum()), len(ps), 0.5)


def test_uniop_beta_in_and_out() -> None:
    intra: list[float] = []
    inter: list[float] = []
    for lay in _pool():
        a, b = _pair_ps(lay)
        intra += a
        inter += b
    assert stats.kstest(intra, stats.beta(8, 2).cdf).pvalue >= KS_MIN_P
    assert stats.kstest(inter, stats.beta(2, 8).cdf).pvalue >= KS_MIN_P


def test_uniop_beta_knobs_honoured() -> None:
    intra: list[float] = []
    inter: list[float] = []
    for i in range(8):
        a, b = _pair_ps(
            _layout(BASE + 100 + i, uniop_beta_in=(3, 3), uniop_beta_out=(1, 4), regions=())
        )
        intra += a
        inter += b
    assert stats.kstest(intra, stats.beta(3, 3).cdf).pvalue >= KS_MIN_P
    assert stats.kstest(inter, stats.beta(1, 4).cdf).pvalue >= KS_MIN_P


def test_prophage_subset_probabilities() -> None:
    counts = dict.fromkeys(lo.SUBSETS, 0)
    for lay in _pool():
        rg = lay.region_genes
        for s in rg.loc[rg["region_id"] == "prophage_1", "subset"]:
            counts[s] += 1
    p = lo.chi2_gof_p([counts[s] for s in lo.SUBSETS], lo.SUBSET_P)
    assert p >= KS_MIN_P, (counts, p)


def test_core_biosynthetic_thirty_percent() -> None:
    k = n = 0
    for lay in _pool():
        rg = lay.region_genes
        bgc = rg[rg["region_id"] == "bgc_1"]
        n += len(bgc)
        k += int(bgc["core_biosynthetic"].astype(bool).sum())
    assert n == 25 * N
    _within(k, n, 0.3)


def test_xrefs_for_seventy_percent(tmp_path: Path) -> None:
    from gippasim.contract_tables import write_contract_tables

    k = n = 0
    for i, lay in enumerate(_pool()[:4]):
        work = tmp_path / f"w{i}"
        write_contract_tables(lay, work)
        x = pd.read_parquet(work / "xrefs.parquet")
        k += x["locus_tag"].nunique()
        n += len(lay.genes)
    _within(k, n, 0.7)


def test_prophage_not_short_without_attack() -> None:
    ratios = []
    for lay in _pool():
        g = lay.genes.set_index("locus_tag")
        span = (g["end"] - g["start"]).astype(float)
        pro = span.index.isin(lo.members(lay.region_genes, "prophage_1"))
        ratios.append(span[pro].median() / span[~pro].median())
    assert 0.8 <= float(np.median(ratios)) <= 1.25, ratios


def test_short_genes_attack_every_replicate() -> None:
    for i in range(10):
        lay = _layout(BASE + 200 + i, short_genes_in="prophage")
        g = lay.genes.set_index("locus_tag")
        span = (g["end"] - g["start"]).astype(float)
        pro = span.index.isin(lo.members(lay.region_genes, "prophage_1"))
        assert span[pro].median() <= 0.5 * span[~pro].median(), i
