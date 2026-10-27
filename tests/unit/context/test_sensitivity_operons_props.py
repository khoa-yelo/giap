"""M6-10 properties (docs/specs/M6-10-sensitivity-operon-flags.md acceptance criteria 2-3):

* the family is pure (inputs untouched, same set on the same inputs) and order-independent
  (row order of the gene table, ``pairs_raw`` and C5 ``tus`` does not matter);
* the set at any UniOP threshold is a set of ``pair_id`` values of that run's pair frame, equal to
  the chain's ``operon_flags`` IDs; its ``jaccard`` against the primary lies in [0, 1] or is NA;
* no recomputation drift: the primary set equals the ``pair_id`` column of the
  ``operon_flags.tsv`` that the registered ``operons`` step writes on the same inputs.

Inputs: the SIM100K producer chain and the gene tables of tests/helpers/sensitivity_operons_data.py.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import sensitivity_operons_data as so
from helpers.layers_chain import Chain, build_chain
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

SETTINGS = settings(max_examples=12, deadline=None, derandomize=True,
                    suppress_health_check=[HealthCheck.function_scoped_fixture,
                                           HealthCheck.too_slow])  # fmt: skip


def _S() -> Any:
    from gippa.context import sensitivity as S

    return S


def _rng() -> np.random.Generator:
    from gippa.rng import for_step

    return for_step(so.STREAM, 7)


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(Path(tmp_path_factory.mktemp("sens_operons_props")))


@pytest.fixture(scope="module")
def tabs(chain: Chain) -> dict[str, pd.DataFrame]:
    return so.tables(chain)


@pytest.fixture(scope="module")
def gt(tabs: dict[str, pd.DataFrame], tmp_path_factory: pytest.TempPathFactory) -> pd.DataFrame:
    g = so.gene_table(tabs["genes"], imputed=(1, 2, 3, 10, 11))
    return so.roundtrip(g, tmp_path_factory.mktemp("gt") / "g.parquet")


def _perm(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    order = np.random.default_rng(seed).permutation(len(df))
    return df.iloc[order].reset_index(drop=True)


@SETTINGS
@given(seed=st.integers(0, 2**32 - 1), t=st.sampled_from([0.3, 0.5, 0.7]))
def test_pure_and_order_independent(
    seed: int, t: float, gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]
) -> None:
    cfg = so.with_threshold(so.config(), t)
    before = {k: v.copy(deep=True) for k, v in tabs.items()}
    gt0 = gt.copy(deep=True)
    want = so.call_family(gt, cfg, _rng(), tabs)
    pd.testing.assert_frame_equal(gt, gt0)
    for k, v in tabs.items():
        pd.testing.assert_frame_equal(v, before[k], obj=k)
    shuffled = dict(tabs)
    for i, key in enumerate(("pairs_raw", "tus")):
        shuffled[key] = _perm(tabs[key], seed + i)
    got = so.call_family(_perm(gt, seed + 7), cfg, _rng(), shuffled)
    assert got == want


@SETTINGS
@given(t=st.floats(0.05, 0.95, allow_nan=False))
def test_any_threshold_gives_pair_ids_and_a_bounded_jaccard(
    t: float, gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]
) -> None:
    S = _S()
    cfg = so.config()
    vcfg = so.with_threshold(cfg, t)
    var = so.call_family(gt, vcfg, _rng(), tabs)
    merged, _ = so.oracle_frames(gt, vcfg, tabs)
    assert var <= set(merged["pair_id"].astype(str))
    assert var == so.oracle_flagged(gt, vcfg, tabs)
    prim = so.call_family(gt, cfg, _rng(), tabs)
    j = S.jaccard(prim, var)
    if prim or var:
        assert j is not None and 0.0 <= j <= 1.0
    else:
        assert j is None, "lead decision 2: J of two empty sets is NA"


@SETTINGS
@given(
    a=st.frozensets(st.integers(1, 60).map(lambda i: f"PAIR_{i}"), max_size=40),
    b=st.frozensets(st.integers(1, 60).map(lambda i: f"PAIR_{i}"), max_size=40),
)
def test_rows_jaccard_in_unit_interval_or_na(a: frozenset[str], b: frozenset[str]) -> None:
    df = _S().sensitivity_rows({so.FAMILY: a}, {(so.FAMILY, "pair_threshold=0.3"): b})
    r = df.iloc[0]
    assert (r["n_primary"], r["n_variant"], r["n_shared"]) == (len(a), len(b), len(a & b))
    j = r["jaccard"]
    if not (a or b):
        assert j is pd.NA or j is None or (isinstance(j, float) and math.isnan(j))
    else:
        assert 0.0 <= float(j) <= 1.0
        assert float(j) == pytest.approx(len(a & b) / len(a | b))


# --------------------------------------------------------------------------- no drift


def _operons_step(inputs: dict[str, Path], out: Path, cfg: Any) -> pd.DataFrame:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    names = {"pairs": "pairs.tsv", "operons": "operons.tsv", "operon_flags": "operon_flags.tsv"}
    outputs = {k: out / v for k, v in names.items()}
    keys = ("gene_table", "genes", "replicons", "tus", "pairs_raw", "features", *so.LAYER_KEYS)
    reg.run_step("operons", cfg, {k: inputs[k] for k in keys if k in inputs}, outputs)
    return pd.read_csv(outputs["operon_flags"], sep="\t", dtype=str, keep_default_na=False)


@pytest.mark.parametrize(
    ("kind", "t", "drop"),
    [("plain", 0.5, ()), ("plain", 0.3, ()), ("imputed", 0.5, ()),
     ("plain", 0.5, so.LAYER_KEYS)],
)  # fmt: skip
def test_primary_set_equals_the_main_run_operon_flags(
    kind: str, t: float, drop: tuple[str, ...], chain: Chain, tabs: dict[str, pd.DataFrame],
    tmp_path: Path,
) -> None:  # fmt: skip
    """The family on the primary inputs reproduces ``tables/operon_flags.tsv`` of the main run
    (written by the registered ``operons`` step), also at the rerun B = 199."""
    imputed = (1, 2, 3, 10, 11) if kind == "imputed" else ()
    g = so.gene_table(tabs["genes"], imputed=imputed)
    ins = so.write_inputs(chain, tmp_path / "in", g, drop=drop)
    cfg = so.with_threshold(so.config(), t)
    main = _operons_step(ins, tmp_path / "out", cfg)
    want = frozenset(main["pair_id"])
    assert want or t != 0.5 or drop, "fixture: the main run flags no pair"
    from gippa.io.genes_tsv import read_genes_parquet

    gt = read_genes_parquet(ins["gene_table"])
    t_in = {k: pd.read_parquet(ins[k]) for k in so.FAMILY_KW if k in ins}
    assert so.call_family(gt, cfg, _rng(), t_in) == want
    b199 = cfg.model_copy(
        update={"stats": cfg.stats.model_copy(update={"permutations": so.SENS_B})}
    )
    assert so.call_family(gt, b199, _rng(), t_in) == want
