"""M8-16 item 2: the ``determinants`` sensitivity family (docs/specs/M8-16-sensitivity-function-
determinants.md, examples 3-4, degenerate rows, lead decision 2; GIAP_SPEC.md §8.6 "the
determinants CV R² (without its permutation null)", "plus ``agreement`` for the dominant
determinant block").

Interface assumed: ``tests/helpers/sensitivity_function_data.py`` (module docstring). The oracle is
the M8-10 chain (``det_oracle``): model rows of the features frame with the variant's detection /
abundance, ``tu_folds`` from the family's generator, ``cross_validate`` (seed ``cfg.seed``), SHAP
and ``attribution_table``; never ``null_p``.
"""

from __future__ import annotations

import inspect
import math
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import sensitivity_function_data as sf
from hypothesis import given, settings
from hypothesis import strategies as st


def _S() -> Any:
    from giap.context import sensitivity as S

    return S


def _rng() -> np.random.Generator:
    from giap.rng import for_step

    return for_step(f"sensitivity/{sf.DET}", sf.SEED)


@pytest.fixture(scope="module")
def world() -> tuple[pd.DataFrame, pd.DataFrame]:
    return sf.planted_world()


@pytest.fixture
def no_null(monkeypatch: pytest.MonkeyPatch) -> None:
    """§8.6 "without its permutation null": the null is never computed."""
    import giap.determinants.null as N

    def boom(*_a: Any, **_k: Any) -> Any:
        raise AssertionError("the sensitivity rerun computed the permutation null")

    for name in ("null_p", "permute_response"):
        monkeypatch.setattr(N, name, boom)
    import giap.determinants.sensitivity as DS

    for name in ("null_p", "permute_response"):
        monkeypatch.setattr(DS, name, boom, raising=False)


# --------------------------------------------------------------------------- interface


def test_registered_in_families() -> None:
    from giap.determinants.sensitivity import determinants_family

    assert _S().FAMILIES[sf.DET] is determinants_family


def test_signature_binds_features_by_keyword() -> None:
    params = inspect.signature(sf.det_family()).parameters
    assert list(params)[:3] == ["genes", "cfg", "rng"]
    assert params["features"].kind is inspect.Parameter.KEYWORD_ONLY


# --------------------------------------------------------------------------- the rerun


def test_matches_the_m8_chain_without_null(
    world: tuple[pd.DataFrame, pd.DataFrame], no_null: None
) -> None:
    """CV R², the ``abundance_not_explained_by_sequence`` set and the dominant blocks of the
    explained genes (M8-16c) equal the M8-10 chain's on the same folds; ``null_p`` is never
    called."""
    genes, feats = world
    cfg = sf.base_config()
    res = sf.det_result(genes, feats, cfg, _rng())
    want = sf.det_oracle(genes, feats, cfg, _rng())
    assert math.isfinite(want.cv_r2) and want.cv_r2 > 0.3, "fixture: no planted signal"
    assert want.flagged and want.dominant, "fixture: empty sets"
    assert float(res.cv_r2) == pytest.approx(want.cv_r2, rel=1e-9, abs=1e-12)
    assert frozenset(res.flagged) == want.flagged
    assert sf.non_none(res.dominant) == want.dominant
    assert set(res.dominant) == want.explained


def test_uses_the_variant_gene_table(world: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    """The variant's ``abundance_log2`` / ``detected`` / ``abundance_imputed`` come from
    ``genes``, not from the primary features frame (a top3-like change of abundance, and 15
    genes losing their abundance as under ``exclude``)."""
    genes, feats = world
    g = genes.copy()
    r = np.random.default_rng(11)
    g["abundance_log2"] = g["abundance_log2"].to_numpy(float) + r.normal(0, 0.7, len(g))
    g.loc[g.index[:15], "abundance_log2"] = math.nan
    cfg = sf.base_config()
    res = sf.det_result(g, feats, cfg, _rng())
    want = sf.det_oracle(g, feats, cfg, _rng())
    prim = sf.det_oracle(genes, feats, cfg, _rng())
    assert want.cv_r2 != pytest.approx(prim.cv_r2), "fixture: variant does not change R²"
    assert float(res.cv_r2) == pytest.approx(want.cv_r2, rel=1e-9, abs=1e-12)
    assert frozenset(res.flagged) == want.flagged
    assert sf.non_none(res.dominant) == want.dominant
    assert not set(res.dominant) & set(sf.tags(range(15)))


def test_min_and_exclude_leave_the_result_unchanged(
    world: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    """M8-08 rule 1 keeps imputed genes out of the models, so under common random numbers the
    ``min`` and ``exclude`` reruns equal the primary (agreement 1, jaccard 1)."""
    S = _S()
    genes, feats = world
    g = genes.copy()
    hit = g.index[::9]
    med = float(np.median(g["abundance_log2"].drop(index=hit)))
    g.loc[hit, "abundance_log2"] = med
    g.loc[hit, "abundance_imputed"] = True
    cfg = sf.base_config()
    prim = sf.det_result(g, feats, cfg, _rng())
    reruns = {}
    for p in ("min", "exclude"):
        v = S.Variant(f"missing={p}", {"abundance.missing": p})
        reruns[(sf.DET, v.name)] = sf.det_result(S.variant_genes(g, v, cfg), feats, cfg, _rng())
    for r in reruns.values():
        assert float(r.cv_r2) == pytest.approx(float(prim.cv_r2), rel=1e-12)
        assert frozenset(r.flagged) == frozenset(prim.flagged)
    rows = S.sensitivity_rows({sf.DET: prim}, reruns)
    assert (rows["agreement"].astype(float) == 1.0).all()
    same = (rows["n_primary"] == rows["n_variant"]) & (rows["n_shared"] == rows["n_primary"])
    assert same.all()


def test_example_4_same_seed_identical(world: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genes, feats = world
    cfg = sf.base_config()
    a = sf.det_result(genes, feats, cfg, _rng())
    b = sf.det_result(genes, feats, cfg, _rng())
    assert float(a.cv_r2) == float(b.cv_r2)
    assert frozenset(a.flagged) == frozenset(b.flagged)
    assert dict(a.dominant) == dict(b.dominant)


# --------------------------------------------------------------------------- degenerate


def _empty_result(res: Any) -> bool:
    return res is None or (math.isnan(float(res.cv_r2)) and not res.flagged and not res.dominant)


def test_below_the_100_row_minimum_fits_nothing(no_null: None) -> None:
    genes, feats = sf.planted_world(n=60)
    assert _empty_result(sf.det_result(genes, feats, sf.base_config(), _rng()))


def test_no_abundance_fits_nothing(world: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    genes, feats = world
    g = genes.copy()
    g["abundance_log2"] = math.nan
    assert _empty_result(sf.det_result(g, feats, sf.base_config(), _rng()))


# --------------------------------------------------------------------------- rows, agreement


def test_example_3_agreement_two_thirds(tmp_path: Any) -> None:
    """Example 3: dominant blocks (codon, codon, transcription) vs ``exclude`` (codon, protein,
    transcription) on the same 3 top-decile genes -> agreement 2/3; n_* / jaccard on the
    ``abundance_not_explained_by_sequence`` sets (lead decision 2), written to 6 digits."""
    from giap.io.sensitivity_tsv import write_sensitivity_tsv

    S = _S()
    prim = sf.det_sens(["a", "b"], {"g1": "codon", "g2": "codon", "g3": "transcription"}, 0.41)
    var = sf.det_sens(["a"], {"g1": "codon", "g2": "protein", "g3": "transcription"}, 0.38)
    rows = S.sensitivity_rows({sf.DET: prim}, {(sf.DET, "missing=exclude"): var})
    r = rows.iloc[0]
    assert (r["family"], r["variant"]) == (sf.DET, "missing=exclude")
    assert (r["n_primary"], r["n_variant"], r["n_shared"]) == (2, 1, 1)
    assert float(r["jaccard"]) == pytest.approx(0.5)
    assert float(r["agreement"]) == pytest.approx(2 / 3)
    out = tmp_path / "sensitivity.tsv"
    write_sensitivity_tsv(out, rows)
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == list(S.COLUMNS)
    assert lines[1].split("\t")[-1] == "0.666667"


def test_agreement_counts_genes_in_both_top_deciles_only() -> None:
    S = _S()
    prim = sf.det_sens([], {"g1": "codon", "g2": "codon", "g3": "protein", "g4": "codon"}, 0.4)
    var = sf.det_sens([], {"g1": "codon", "g2": "protein", "g3": "protein", "g9": "codon"}, 0.4)
    r = S.sensitivity_rows({sf.DET: prim}, {(sf.DET, "top3"): var}).iloc[0]
    assert float(r["agreement"]) == pytest.approx(2 / 3)
    assert pd.isna(r["jaccard"])


@pytest.mark.parametrize("which", ["primary", "variant", "disjoint"])
def test_empty_top_decile_gives_agreement_na(which: str) -> None:
    S = _S()
    full = {"g1": "codon", "g2": "protein"}
    p = {} if which == "primary" else full
    v = {} if which == "variant" else ({"g7": "codon"} if which == "disjoint" else full)
    rows = S.sensitivity_rows({sf.DET: sf.det_sens(["a"], p, 0.3)},
                              {(sf.DET, "missing=min"): sf.det_sens(["a"], v, 0.3)})  # fmt: skip
    assert pd.isna(rows.iloc[0]["agreement"])
    assert float(rows.iloc[0]["jaccard"]) == 1.0


def test_set_families_keep_agreement_na_beside_determinants() -> None:
    S = _S()
    rows = S.sensitivity_rows(
        {"regions": frozenset({"x"}), sf.DET: sf.det_sens([], {"g1": "codon"}, 0.3)},
        {("regions", "missing=min"): frozenset({"x"}),
         (sf.DET, "missing=min"): sf.det_sens([], {"g1": "codon"}, 0.2)},
    )  # fmt: skip
    by = dict(zip(rows["family"], rows["agreement"], strict=True))
    assert pd.isna(by["regions"])
    assert float(by[sf.DET]) == 1.0


_BLOCKS = st.sampled_from(["transcription", "translation", "codon", "protein", "context"])
_DOM = st.dictionaries(st.sampled_from([f"g{i}" for i in range(12)]), _BLOCKS, max_size=10)


@settings(max_examples=60, deadline=None)
@given(p=_DOM, v=_DOM)
def test_agreement_property(p: dict[str, str], v: dict[str, str]) -> None:
    """Property: agreement = |same block on shared top genes| / |shared top genes| (NA when no
    gene is shared), symmetric in primary / variant, in [0, 1]."""
    S = _S()

    def agree(a: dict[str, str], b: dict[str, str]) -> Any:
        rows = S.sensitivity_rows({sf.DET: sf.det_sens([], a, 0.1)},
                                  {(sf.DET, "v"): sf.det_sens([], b, 0.1)})  # fmt: skip
        return rows.iloc[0]["agreement"]

    shared = set(p) & set(v)
    ab, ba = agree(p, v), agree(v, p)
    if not shared:
        assert pd.isna(ab) and pd.isna(ba)
        return
    want = sum(p[g] == v[g] for g in shared) / len(shared)
    assert float(ab) == pytest.approx(want)
    assert float(ba) == pytest.approx(want)
    assert 0.0 <= float(ab) <= 1.0
