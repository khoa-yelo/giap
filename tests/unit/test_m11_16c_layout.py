"""M11-16c: the sequence-free layout at genetic codes 4 and 25 (docs/specs/
M11-16-translation-tables-4-25.md item 8: "Every module that reads them (``features``,
``layout``, ...) takes the table from the config"; item 9 gives ``test_sim_layout_calibration``
arms 4 and 25; lead decision on #474).

Interface assumed: ``LayoutConfig.genetic_code: int = 11`` in {4, 11, 25} (else ``ValueError``),
used by ``simulate_layout``: the back-translation draws synonymous codons from
``gippasim.codons.sense_usage(table)`` (TGA joins W at 4, G at 25, at the table's default
``tga_sense_share``), terminal and pseudogene stops come from the table's stops with the
default weights renormalised (TAA 0.61 / TAG 0.09 -> 0.871 / 0.129). Code 11 is pinned in
``test_m11_16c_identity.py``.

Shares are pooled over ``N_SEEDS`` layouts (``calibration_base + i``) and must lie within 3
binomial SEs (spec AC 14).
"""

from __future__ import annotations

import hashlib
import math
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import yaml

TOL = Path(__file__).resolve().parents[1] / "calibration" / "tolerances.yaml"
BASE = int(yaml.safe_load(TOL.read_text())["seeds"]["calibration_base"])
N_SEEDS = 5
N_GENES = 400
FAMILY = {4: "W", 25: "G"}
STOPS_4_25 = frozenset({"TAA", "TAG"})


def _layout(table: int, seed: int, **kw: Any) -> Any:
    from gippasim.layout import LayoutConfig, simulate_layout

    cfg = LayoutConfig(n_genes_chromosome=N_GENES, genetic_code=table, **kw)
    return simulate_layout(cfg, np.random.default_rng(seed))


def _codons(nt: str) -> list[str]:
    return [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]


def _within_3se(k: int, n: int, p: float) -> None:
    se = math.sqrt(p * (1 - p) / n)
    assert abs(k / n - p) <= 3 * se, f"{k}/{n} = {k / n:.4f}, expected {p:.4f} +- 3 SE ({se:.4f})"


@cache
def _genes(table: int) -> pd.DataFrame:
    frames = [_layout(table, BASE + i).genes for i in range(N_SEEDS)]
    return pd.concat(frames, ignore_index=True)


@pytest.mark.parametrize("table", [0, 1, 12, 26])
def test_layout_refuses_other_tables(table: int) -> None:
    from gippasim.layout import LayoutConfig

    with pytest.raises(ValueError):
        LayoutConfig(genetic_code=table)


@pytest.mark.parametrize("table", [4, 25])
def test_layout_stops_only_from_table(table: int) -> None:
    """Every gene ends in TAA/TAG; intact genes have no internal TAA/TAG; each pseudogene's
    planted internal stop is TAA/TAG (TGA is sense at 4/25)."""
    g = _genes(table)
    assert (g["feature_type"] == "pseudogene").sum() > 0
    for ftype, nt in zip(g["feature_type"], g["nt"], strict=True):
        cod = _codons(nt)
        assert cod[-1] in STOPS_4_25, (ftype, cod[-1])
        inner = STOPS_4_25.intersection(cod[1:-1])
        if ftype == "pseudogene":
            assert inner, "pseudogene without a TAA/TAG internal stop"
        else:
            assert not inner


@pytest.mark.parametrize("table", [4, 25])
def test_layout_aa_is_translation_under_table(table: int) -> None:
    from gippasim.codons import code

    t = code(table)
    g = _genes(table)
    intact = g[g["feature_type"] != "pseudogene"]
    assert len(intact) > 0
    for nt, aa in zip(intact["nt"], intact["aa"], strict=True):
        cod = _codons(nt)
        assert aa == "M" + "".join(t[c] for c in cod[1:-1])
        assert "*" not in aa


@pytest.mark.calibration
@pytest.mark.parametrize("table", [4, 25])
def test_layout_tga_share_of_family(table: int) -> None:
    """TGA's share of the W (4) / G (25) codons in gene bodies is the table's usage share."""
    from gippasim.codons import code, sense_usage

    t = code(table)
    u = sense_usage(table)
    fam = [c for c in u if t[c] == FAMILY[table]]
    share = u["TGA"] / sum(u[c] for c in fam)
    g = _genes(table)
    k = n = 0
    for ftype, nt in zip(g["feature_type"], g["nt"], strict=True):
        if ftype == "pseudogene":
            continue
        body = _codons(nt)[1:-1]
        k += body.count("TGA")
        n += sum(c in fam for c in body)
    _within_3se(k, n, share)


@pytest.mark.calibration
@pytest.mark.parametrize("table", [4, 25])
def test_layout_terminal_stop_weights(table: int) -> None:
    g = _genes(table)
    last = [_codons(nt)[-1] for nt in g["nt"]]
    _within_3se(last.count("TAA"), len(last), 0.61 / 0.70)


def _digest(layout: Any) -> str:
    h = hashlib.sha256()
    for n in ("replicons", "genes", "tus", "operon_pairs", "regions", "region_genes", "truth"):
        h.update(getattr(layout, n).to_csv().encode())
    return h.hexdigest()


@pytest.mark.parametrize("table", [4, 25])
def test_layout_deterministic_per_seed(table: int) -> None:
    a, b, c = _layout(table, 3), _layout(table, 3), _layout(table, 4)
    assert _digest(a) == _digest(b)
    assert _digest(a) != _digest(c)
