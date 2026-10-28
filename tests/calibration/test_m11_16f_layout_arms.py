"""Genetic-code arms 4 and 25 of the gippasim layout calibration (M11-16f; docs/specs/
M11-16-translation-tables-4-25.md item 9 row ``test_sim_layout_calibration.py``: "sense usage and
stop weights at 4/25", arms 4 and 25 in both tiers; item 8).

Arm 11 is ``test_sim_layout_calibration.py``, unchanged. The same pool (``N`` = 20 layouts of
1,000 genes, seeds ``calibration_base + i``) and the same criteria (4 binomial SDs, KS /
chi-square at ``ks_min_p``) at ``LayoutConfig(genetic_code=t)``: the M4-07 knob distributions do
not move when the table changes, and the table's own claims hold on that pool (pseudogene
internal stops TAA/TAG with the renormalised weights 0.871 / 0.129; TGA's share of the W (4) / G
(25) codons is ``sense_usage(t)``'s). ``tests/unit/test_m11_16c_layout.py`` pins the same table
claims at 3 SE on 5 x 400 genes; these arms run on the calibration pool. The suite has no
``calibration_full`` tests, so the fast arms are the only tier, as for arm 11.
"""

from __future__ import annotations

import math
import sys
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import yaml
from helpers.genetic_code_arms import SENSE_TABLES, TGA_AA
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import layout_oracle as lo

TOL = yaml.safe_load((Path(__file__).parent / "tolerances.yaml").read_text())
BASE = int(TOL["seeds"]["calibration_base"])
KS_MIN_P = float(TOL["ks_min_p"])
N = 20

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("missed_operons")]
TABLE = pytest.mark.parametrize("table", SENSE_TABLES)


@cache
def _pool(table: int) -> tuple[Any, ...]:
    from gippasim.layout import LayoutConfig, simulate_layout

    cfg = LayoutConfig(genetic_code=table)
    return tuple(simulate_layout(cfg, np.random.default_rng(BASE + i)) for i in range(N))


def _genes(table: int) -> pd.DataFrame:
    return pd.concat([lay.genes for lay in _pool(table)], ignore_index=True)


def _within(k: float, n: int, p: float) -> None:
    assert lo.within_sd(k, n, p), (k, n, p, k / n)


def _codons(nt: str) -> list[str]:
    return [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]


@TABLE
def test_tu_size_geometric_mean(table: int) -> None:
    sizes = np.array([len(lo.as_list(c)) for lay in _pool(table) for c in lay.tus["genes"]])
    m, p = 1.8, 1 / 1.8
    se = math.sqrt((1 - p) / p**2 / len(sizes))
    assert abs(sizes.mean() - m) <= 4 * se, (table, sizes.mean(), se)
    _within(int((sizes == 1).sum()), len(sizes), p)


@TABLE
def test_cds_lengths_lognormal(table: int) -> None:
    g = _genes(table)
    cds = g[g["feature_type"] == "cds"]
    logs = np.log((cds["end"] - cds["start"]).to_numpy() / 3 - 1)
    assert stats.kstest(logs, stats.norm(math.log(280.0), 0.45).cdf).pvalue >= KS_MIN_P
    se = 1.2533 * 0.45 / math.sqrt(len(logs))
    assert abs(np.median(logs) - math.log(280.0)) <= 4 * se


@TABLE
@pytest.mark.parametrize(
    ("what", "p"), [("sorf", 0.02), ("pseudogene", 0.01), ("hypothetical", 0.25)]
)
def test_gene_class_fractions(table: int, what: str, p: float) -> None:
    g = _genes(table)
    k = (
        int(g["hypothetical"].astype(bool).sum())
        if what == "hypothetical"
        else int((g["feature_type"] == what).sum())
    )
    _within(k, len(g), p)


@TABLE
def test_pseudogene_aa_na_half(table: int) -> None:
    ps = _genes(table).query("feature_type == 'pseudogene'")
    assert len(ps) >= 100
    _within(int(ps["aa"].isna().sum()), len(ps), 0.5)


@TABLE
def test_uniop_beta_in_and_out(table: int) -> None:
    intra: list[float] = []
    inter: list[float] = []
    for lay in _pool(table):
        tmap = lo.tu_of(lay.tus)
        pairs = lay.operon_pairs
        col = lo.pair_prob_column(pairs)
        for a, b, p in zip(pairs["gene_a"], pairs["gene_b"], pairs[col], strict=True):
            (intra if tmap[a] == tmap[b] else inter).append(float(p))
    assert stats.kstest(intra, stats.beta(8, 2).cdf).pvalue >= KS_MIN_P
    assert stats.kstest(inter, stats.beta(2, 8).cdf).pvalue >= KS_MIN_P


@TABLE
def test_terminal_stops_taa_tag_with_renormalised_weights(table: int) -> None:
    last = [_codons(nt)[-1] for nt in _genes(table)["nt"]]
    assert set(last) <= {"TAA", "TAG"}
    _within(last.count("TAA"), len(last), 0.61 / 0.70)


@TABLE
def test_pseudogene_internal_stops_are_taa_tag(table: int) -> None:
    k = n = 0
    for ftype, nt in zip(*(_genes(table)[c] for c in ("feature_type", "nt")), strict=True):
        if ftype != "pseudogene":
            continue
        stops_ = [c for c in _codons(nt)[1:-1] if c in {"TAA", "TAG"}]
        assert stops_, nt  # TGA is sense: the planted stop is TAA or TAG
        k += stops_.count("TAA")
        n += len(stops_)
    assert n >= 100
    _within(k, n, 0.61 / 0.70)


@TABLE
def test_tga_share_of_family_on_the_pool(table: int) -> None:
    from gippasim.codons import code, sense_usage

    t, u = code(table), sense_usage(table)
    fam = [c for c in u if t[c] == TGA_AA[table]]
    share = u["TGA"] / sum(u[c] for c in fam)
    k = n = 0
    for ftype, nt in zip(*(_genes(table)[c] for c in ("feature_type", "nt")), strict=True):
        if ftype == "pseudogene":
            continue
        body = _codons(nt)[1:-1]
        assert "TAA" not in body and "TAG" not in body
        k += body.count("TGA")
        n += sum(c in fam for c in body)
    _within(k, n, share)
