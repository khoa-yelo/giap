"""M6-10 ``operon_flags`` sensitivity family: registration, signature, the ``config_patch``
variant kind (lead decision 1), the flagged-set rule and spec examples 1-5 and 7 plus the
degenerate cases on the giapsim SIM100K producer chain
(docs/specs/M6-10-sensitivity-operon-flags.md; GIAP_SPEC.md §8.6, G9; M4-19 decisions 2, 5).

The interface assumed is spelled out in ``tests/helpers/sensitivity_operons_data.py``. The
reference flagged set is the ``pair_id`` column of the merged M6-08 chain's ``operon_flags``
frame (``run_operons``) on the same tables and config. Step-level checks (example 6, the patch
reaching the step, CRN, determinism) live in
tests/integration/test_sensitivity_operons_step.py; the properties in
test_sensitivity_operons_props.py.
"""

from __future__ import annotations

import dataclasses
import inspect
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import sensitivity_operons_data as so
from helpers.layers_chain import Chain, build_chain

PRIMARY_T = 0.5  # GiapConfig default annotation.uniop.pair_threshold
OP_FLAGS = frozenset({"internal_promoter_supported", "internal_terminator_supported",
                      "unexplained_discordance", "conflict_terminator_in_operon"})  # fmt: skip
NO_FLAGS = frozenset({"novel_operon_candidate", "coabundant_despite_terminator"})
IMPUTED = (1, 2, 3, 10, 11, 12, 20, 21)
SINGLE = (4, 5, 13, 14, 22, 23, 31, 32)


def _S() -> Any:
    from giap.context import sensitivity as S

    return S


def _rng(seed: int = 20261022) -> np.random.Generator:
    from giap.rng import for_step

    return for_step(so.STREAM, seed)


def _is_na(x: Any) -> bool:
    return x is None or x is pd.NA or (isinstance(x, float) and math.isnan(x))


def _names(vs: list[Any]) -> list[str]:
    return [v.name for v in vs]


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(Path(tmp_path_factory.mktemp("sens_operons_chain")))


@pytest.fixture(scope="module")
def tabs(chain: Chain) -> dict[str, pd.DataFrame]:
    return so.tables(chain)


@pytest.fixture(scope="module")
def gt(chain: Chain, tabs: dict[str, pd.DataFrame], tmp_path_factory: pytest.TempPathFactory
       ) -> pd.DataFrame:  # fmt: skip
    return so.roundtrip(so.gene_table(tabs["genes"]), tmp_path_factory.mktemp("gt") / "g.parquet")


# --------------------------------------------------------------------------- interface


def test_operon_flags_family_is_registered() -> None:
    S = _S()
    assert S.FAMILIES["operon_flags"] is S.operon_flags_family
    assert "operon_flags_family" in S.__all__
    assert {"regions", "islands"} <= set(S.FAMILIES), "M4-19 / M5-10 families stay registered"


def test_signature_binds_tables_by_keyword() -> None:
    """Behaviour 1: ``(genes, cfg, rng)`` positional, the tables keyword-only (bound by the step
    as M4-19 decision 5). C1 / C5 / ``pairs_raw`` / the C2 gene model / C3 features are
    required; the three evidence layers default to ``None`` (missing layer)."""
    params = inspect.signature(_S().operon_flags_family).parameters
    assert list(params)[:3] == ["genes", "cfg", "rng"]
    for key in ("pairs_raw", "tus", "replicons", "genes", "features"):
        p = params[so.FAMILY_KW[key]]
        assert p.kind is inspect.Parameter.KEYWORD_ONLY, key
        assert p.default is inspect.Parameter.empty, key
    for key in so.LAYER_KEYS:
        p = params[so.FAMILY_KW[key]]
        assert p.kind is inspect.Parameter.KEYWORD_ONLY, key
        assert p.default is None, key
    assert {"regions", "region_genes"}.isdisjoint(params), "operons do not read C7 / C8"


@pytest.mark.parametrize("missing", ["pairs_raw", "tus", "replicons", "genes"])
def test_family_requires_its_tables(
    missing: str, gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]
) -> None:
    t = {k: v for k, v in tabs.items() if k != missing}
    with pytest.raises(TypeError):
        so.call_family(gt, so.config(), _rng(), t)


# --------------------------------------------------------------------------- Variant.config_patch


def test_variant_has_an_optional_config_patch() -> None:
    """Lead decision 1: ``Variant`` gains an optional config patch; a variant built as before
    (three fields) has an empty one and is otherwise unchanged (frozen, slotted)."""
    S = _S()
    names = [f.name for f in dataclasses.fields(S.Variant)]
    assert names[:3] == ["name", "overrides", "single_peptide_detected"]
    assert "config_patch" in names
    v = S.Variant("missing=min", {"abundance.missing": "min"})
    assert not v.config_patch
    assert v.single_peptide_detected is False
    p = S.Variant("pair_threshold=0.3", {}, config_patch={so.PATCH_KEY: 0.3})
    assert dict(p.config_patch) == {so.PATCH_KEY: 0.3}
    assert S.Variant.__dataclass_params__.frozen
    assert hasattr(S.Variant, "__slots__")
    with pytest.raises(dataclasses.FrozenInstanceError):
        p.config_patch = {}  # type: ignore[misc]


@pytest.mark.parametrize("tier", ["mayu", "two_peptide"])
def test_basic_variants_unchanged_without_patch(tier: str) -> None:
    """Lead decision 1, "existing variants are byte-identical": under ``basic`` every variant
    equals the three-field ``Variant`` M4-19 builds, and none carries a patch."""
    S = _S()
    vs = S.variants(so.config(), tier)
    want = ["missing=exclude"] + (["single_peptide_detected"] if tier == "two_peptide" else [])
    assert _names(vs) == want
    for v in vs:
        assert not v.config_patch, v.name
        assert v == S.Variant(v.name, v.overrides, v.single_peptide_detected)
        assert not any("uniop" in str(k) for k in v.overrides), v.name


def test_full_adds_the_two_uniop_threshold_variants() -> None:
    """Behaviour 3: ``full`` adds ``pair_threshold`` 0.3 and 0.7 after ``top3``; they patch the
    config (classes) and change no abundance (empty ``overrides``, no single-peptide switch)."""
    S = _S()
    vs = S.variants(so.config(output={"sensitivity": "full"}), "two_peptide")
    assert _names(vs) == ["missing=exclude", "single_peptide_detected", "top3",
                          "pair_threshold=0.3", "pair_threshold=0.7"]  # fmt: skip
    for v, t in zip(vs[-2:], so.UNIOP_THRESHOLDS, strict=True):
        assert dict(v.config_patch) == {so.PATCH_KEY: t}
        assert not v.overrides
        assert v.single_peptide_detected is False
    for v in vs[:-2]:
        assert not v.config_patch, v.name


@pytest.mark.parametrize("primary", so.UNIOP_THRESHOLDS)
def test_primary_threshold_is_not_rerun(primary: float) -> None:
    """As M4-19 example 5: a threshold equal to the primary is not a variant."""
    S = _S()
    cfg = so.config(
        output={"sensitivity": "full"}, annotation={"uniop": {"pair_threshold": primary}}
    )
    names = _names(S.variants(cfg, "mayu"))
    other = [t for t in so.UNIOP_THRESHOLDS if t != primary]
    assert [n for n in names if n.startswith("pair_threshold=")] == [
        so.variant_name(t) for t in other
    ]


def test_basic_has_no_uniop_variant() -> None:
    names = _names(_S().variants(so.config(), "mayu"))
    assert not any(n.startswith("pair_threshold=") for n in names)


def test_variant_genes_ignores_the_patch(gt: pd.DataFrame) -> None:
    """A patch-only variant re-derives classes, not abundances: ``variant_genes`` gives the
    primary gene table back (abundance, imputation, detection and flags)."""
    S = _S()
    v = S.Variant("pair_threshold=0.3", {}, config_patch={so.PATCH_KEY: 0.3})
    got = S.variant_genes(gt, v, so.config())
    for col in ("abundance_log2", "abundance_imputed", "detected", "detection_status", "flags"):
        assert got[col].tolist() == gt[col].tolist() or (
            col == "abundance_log2"
            and np.allclose(got[col].astype(float), gt[col].astype(float), equal_nan=True)
        ), col


# --------------------------------------------------------------------------- flagged-set rule


@pytest.mark.parametrize("t", [0.3, PRIMARY_T, 0.7])
def test_flagged_set_is_the_operon_flags_pair_ids(
    t: float, gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]
) -> None:
    """Behaviour 2: the set is the ``pair_id`` column of ``operon_flags`` for that config."""
    cfg = so.with_threshold(so.config(), t)
    got = so.call_family(gt, cfg, _rng(), tabs)
    assert isinstance(got, frozenset)
    assert all(isinstance(x, str) and x.startswith("PAIR_") for x in got)
    assert got == so.oracle_flagged(gt, cfg, tabs)


def test_fixture_primary_has_flags(gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]) -> None:
    assert so.oracle_flagged(gt, so.config(), tabs), "fixture: the primary run flags no pair"


def test_missing_layers_are_missing_not_absent_evidence(
    gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]
) -> None:
    """A layer the family is not given is a missing layer, as in the ``operons`` step."""
    t = {k: v for k, v in tabs.items() if k not in so.LAYER_KEYS}
    cfg = so.config()
    assert so.call_family(gt, cfg, _rng(), t) == so.oracle_flagged(gt, cfg, t)


def test_family_ignores_the_generator(gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]) -> None:
    """Behaviour 4: Module 4 draws no random number, so the set does not depend on the stream."""
    cfg = so.config()
    assert so.call_family(gt, cfg, _rng(1), tabs) == so.call_family(gt, cfg, _rng(2), tabs)


def test_b_199_does_not_change_the_set(gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]) -> None:
    """Behaviour 4: reruns set ``stats.permutations = 199``; Module 4 has no permutation test."""
    cfg = so.config()
    b = so.config(stats={"permutations": so.SENS_B})
    assert so.call_family(gt, b, _rng(), tabs) == so.call_family(gt, cfg, _rng(), tabs)


# --------------------------------------------------------------------------- examples


def _row(primary: frozenset[str], variant: frozenset[str], name: str = "missing=exclude") -> Any:
    df = _S().sensitivity_rows({so.FAMILY: primary}, {(so.FAMILY, name): variant})
    assert len(df) == 1
    return df.iloc[0]


def test_example_1_jaccard_36_of_42() -> None:
    p = frozenset(f"PAIR_{i}" for i in range(1, 41))
    v = frozenset(f"PAIR_{i}" for i in range(5, 41)) | {"PAIR_41", "PAIR_42"}
    r = _row(p, v)
    assert (r["family"], r["variant"]) == (so.FAMILY, "missing=exclude")
    assert (r["n_primary"], r["n_variant"], r["n_shared"]) == (40, 38, 36)
    assert r["jaccard"] == pytest.approx(36 / 42)
    assert round(float(r["jaccard"]), 3) == 0.857
    assert _is_na(r["agreement"]), "agreement belongs to the determinants family (Behaviour 6)"


def test_example_2_exclude_drops_imputed_pairs(
    chain: Chain, tabs: dict[str, pd.DataFrame], tmp_path: Path
) -> None:
    """``missing = exclude``: a pair with an imputed gene is never flagged in the variant; the
    variant set is the chain's set on ``variant_genes``; n_variant <= n_primary, J in [0, 1]."""
    S = _S()
    cfg = so.config()
    gt = so.roundtrip(so.gene_table(tabs["genes"], imputed=IMPUTED), tmp_path / "g.parquet")
    (v,) = [v for v in S.variants(cfg, "mayu") if v.name == "missing=exclude"]
    vg = S.variant_genes(gt, v, cfg)
    vcfg = so.config(abundance={"missing": "exclude"})
    prim = so.call_family(gt, cfg, _rng(), tabs)
    var = so.call_family(vg, vcfg, _rng(), tabs)
    assert var == so.oracle_flagged(vg, vcfg, tabs)
    merged, _ = so.oracle_frames(vg, vcfg, tabs)
    imp_tags = set(gt["locus_tag"].iloc[list(IMPUTED)])
    pairs = zip(merged["pair_id"], merged["gene_a"], merged["gene_b"], strict=True)
    with_imp = {str(p) for p, a, b in pairs if a in imp_tags or b in imp_tags}
    assert with_imp, "fixture: no pair touches an imputed gene"
    assert var.isdisjoint(with_imp)
    r = _row(prim, var)
    assert r["n_variant"] <= r["n_primary"]
    assert _is_na(r["jaccard"]) or 0.0 <= float(r["jaccard"]) <= 1.0


def _pair_genes(gt: pd.DataFrame, t: float, tabs: dict[str, pd.DataFrame]) -> dict[str, tuple]:
    merged, _ = so.oracle_frames(gt, so.with_threshold(so.config(), t), tabs)
    return {str(p): (a, b) for p, a, b in zip(merged["pair_id"], merged["gene_a"],
                                              merged["gene_b"], strict=True)}  # fmt: skip


@pytest.mark.parametrize("t", so.UNIOP_THRESHOLDS)
def test_examples_3_4_threshold_variant(
    t: float, gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]
) -> None:
    """0.3 (more OP pairs) / 0.7 (more NO pairs): the variant set is the chain's set at that
    threshold, on the *primary* gene table (no abundance re-derived), and every ``pair_id``
    names the same gene pair as in the primary run (M6-04 decision 3), so J is meaningful."""
    cfg = so.config()
    vcfg = so.with_threshold(cfg, t)
    var = so.call_family(gt, vcfg, _rng(), tabs)
    assert var == so.oracle_flagged(gt, vcfg, tabs)
    assert _pair_genes(gt, t, tabs) == _pair_genes(gt, PRIMARY_T, tabs)
    m_p, _ = so.oracle_frames(gt, cfg, tabs)
    m_v, f_v = so.oracle_frames(gt, vcfg, tabs)
    n_op_p, n_op_v = int((m_p["class"] == "OP").sum()), int((m_v["class"] == "OP").sum())
    assert (n_op_v > n_op_p) if t < PRIMARY_T else (n_op_v < n_op_p), "fixture: classes move"
    prim = so.call_family(gt, cfg, _rng(), tabs)
    r = _row(prim, var, so.variant_name(t))
    assert _is_na(r["jaccard"]) or 0.0 <= float(r["jaccard"]) <= 1.0
    if t < PRIMARY_T:  # example 3 on this fixture: discordance flags rise, NO flags do not
        _, f_p = so.oracle_frames(gt, cfg, tabs)
        assert f_v["flag"].isin(OP_FLAGS).sum() >= f_p["flag"].isin(OP_FLAGS).sum()
        assert f_v["flag"].isin(NO_FLAGS).sum() <= f_p["flag"].isin(NO_FLAGS).sum()


def test_degenerate_07_leaves_too_few_op_pairs(
    gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]
) -> None:
    """At 0.7 fewer than 20 OP pairs are usable: discordance is ``unestimated`` and the variant
    set holds only pairs that are NO in that run (missed-operon flags)."""
    vcfg = so.with_threshold(so.config(), 0.7)
    merged, _ = so.oracle_frames(gt, vcfg, tabs)
    op = merged[merged["class"] == "OP"]
    assert len(op) and all("unestimated" in list(f) for f in op["flags"]), "fixture: OP tested"
    var = so.call_family(gt, vcfg, _rng(), tabs)
    no_ids = set(merged.loc[merged["class"] == "NO", "pair_id"].astype(str))
    assert var <= no_ids
    assert var == so.oracle_flagged(gt, vcfg, tabs)


def test_example_5_two_peptide_variant(chain: Chain, tabs: dict[str, pd.DataFrame],
                                       tmp_path: Path) -> None:  # fmt: skip
    """``single_peptide_detected``: more detected genes -> more both-detected tests; the variant
    set is the chain's set on ``variant_genes``."""
    S = _S()
    cfg = so.config()
    gt = so.roundtrip(so.gene_table(tabs["genes"], single_peptide=SINGLE, tier="two_peptide"),
                      tmp_path / "g.parquet")  # fmt: skip
    (v,) = [v for v in S.variants(cfg, "two_peptide") if v.single_peptide_detected]
    vg = S.variant_genes(gt, v, cfg)
    var = so.call_family(vg, cfg, _rng(), tabs)
    assert var == so.oracle_flagged(vg, cfg, tabs)
    m_p, _ = so.oracle_frames(gt, cfg, tabs)
    m_v, _ = so.oracle_frames(vg, cfg, tabs)
    n_p = int((m_p["test"] == "both_detected_t").sum())
    n_v = int((m_v["test"] == "both_detected_t").sum())
    assert n_v > n_p


def test_no_abundance_gives_the_empty_set(tabs: dict[str, pd.DataFrame], tmp_path: Path) -> None:
    """Degenerate case 1 (sim100k without abundance): no test, no flag; J of two empty sets is NA
    (lead decision 2)."""
    gt = so.roundtrip(so.gene_table(tabs["genes"], abundance=False), tmp_path / "g.parquet")
    for t in (PRIMARY_T, *so.UNIOP_THRESHOLDS):
        assert so.call_family(gt, so.with_threshold(so.config(), t), _rng(), tabs) == frozenset()
    r = _row(frozenset(), frozenset())
    assert (r["n_primary"], r["n_variant"], r["n_shared"]) == (0, 0, 0)
    assert _is_na(r["jaccard"])


def test_relative_scale_gives_the_empty_set(tabs: dict[str, pd.DataFrame], tmp_path: Path) -> None:
    """Degenerate case 2: on a relative scale neither test runs, so the set is empty."""
    gt = so.roundtrip(so.gene_table(tabs["genes"], scale="relative"), tmp_path / "g.parquet")
    assert so.call_family(gt, so.config(), _rng(), tabs) == frozenset()


def test_example_7_same_inputs_same_set(gt: pd.DataFrame, tabs: dict[str, pd.DataFrame]) -> None:
    cfg = so.config()
    assert so.call_family(gt, cfg, _rng(), tabs) == so.call_family(gt, cfg, _rng(), tabs)
