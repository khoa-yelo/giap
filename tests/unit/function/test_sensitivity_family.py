"""M8-16 item 1: the ``enriched_categories`` sensitivity family (docs/specs/M8-16-sensitivity-
function-determinants.md, examples 1, 2 and 4, degenerate rows, lead decision 1; GIAP_SPEC.md §8.6
"enriched categories (q ≤ 0.05)", "Sensitivity reruns use B = 199").

Interface assumed: ``tests/helpers/sensitivity_function_data.py`` (module docstring). GSEA runs
through an injected engine (``FakeEngine`` / ``SIMPLE_ENGINE``), so no R is needed.
"""

from __future__ import annotations

import inspect
import math
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import sensitivity_function_data as sf
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st


def _S() -> Any:
    from giap.context import sensitivity as S

    return S


def _rng() -> np.random.Generator:
    from giap.rng import for_step

    return for_step(f"sensitivity/{sf.EC}", sf.SEED)


def _variant(name: str) -> Any:
    S = _S()
    if name == "single_peptide_detected":
        return S.Variant(name, {}, single_peptide_detected=True)
    return S.Variant(name, {"abundance.missing": name.split("=")[1]})


# --------------------------------------------------------------------------- interface


def test_registered_in_families() -> None:
    from giap.function.sensitivity import enriched_categories_family

    assert _S().FAMILIES[sf.EC] is enriched_categories_family
    assert {"regions", "islands", "operon_flags"} <= set(_S().FAMILIES)


def test_signature_binds_tables_by_keyword() -> None:
    """M4-19 decisions 5, 13: ``(genes, cfg, rng)`` positional, the tables keyword-only; the
    engine is keyword-only with a default (the step binds R fgsea, tests a fake)."""
    params = inspect.signature(sf.ec_family()).parameters
    names = list(params)
    assert names[:3] == ["genes", "cfg", "rng"]
    for kw in ("collections", "tus", "engine"):
        assert params[kw].kind is inspect.Parameter.KEYWORD_ONLY, kw
    assert params["engine"].default is not inspect.Parameter.empty


# --------------------------------------------------------------------------- IDs and alpha


def test_ids_are_collection_term_kind() -> None:
    genes = sf.null_genes()
    colls = [sf.gene_sets("cog", {"K": sf.tags(range(12))}),
             sf.gene_sets("go_bp", {"GO:0000001": sf.tags(range(12, 24))})]  # fmt: skip
    got = sf.call_ec(genes, sf.base_config(), _rng(), colls, sf.constant_engine(0.001))
    assert got == {"cog:K:gsea", "go_bp:GO:0000001:gsea"}


@pytest.mark.parametrize(("p", "alpha", "flagged"), [
    (0.05, None, True), (0.0500001, None, False), (0.08, 0.1, True), (0.08, None, False),
])  # fmt: skip
def test_q_at_most_alpha_enrichment(p: float, alpha: float | None, flagged: bool) -> None:
    """One term per collection, so q = p: q <= ``stats.alpha.enrichment`` (inclusive)."""
    over = {} if alpha is None else {"stats": {"alpha": {"enrichment": alpha}}}
    cfg = sf.base_config(**over)
    colls = [sf.gene_sets("kegg", {"map90001": sf.tags(range(12))})]
    got = sf.call_ec(sf.null_genes(), cfg, _rng(), colls, sf.constant_engine(p))
    assert got == ({"kegg:map90001:gsea"} if flagged else set())


def test_nan_p_is_never_flagged() -> None:
    colls = [sf.gene_sets("kegg", {"map90001": sf.tags(range(12))})]
    got = sf.call_ec(sf.null_genes(), sf.base_config(), _rng(), colls, sf.constant_engine(math.nan))
    assert got == frozenset()


@settings(max_examples=30, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(p=st.lists(st.sampled_from([1e-4, 0.004, 0.01, 0.02, 0.03, 0.049, 0.05, 0.2, 0.9]),
                  min_size=1, max_size=6))  # fmt: skip
def test_gsea_ids_are_bh_within_collection(p: list[float]) -> None:
    """Property: the GSEA IDs are the terms with BH(p) <= alpha within their collection (the
    fixture's ORA flags nothing)."""
    cfg = sf.base_config()
    genes = sf.null_genes()
    colls = sf.null_collections("kegg", len(p))
    by_term = {f"T{j + 1}": v for j, v in enumerate(p)}
    assert sf.ora_ids(genes, colls, cfg) == frozenset()
    got = sf.call_ec(genes, cfg, _rng(), colls, sf.table_engine(by_term))
    assert got == sf.gsea_ids(by_term, "kegg", cfg.stats.alpha.enrichment)


# --------------------------------------------------------------------------- B = 199


@pytest.mark.parametrize("n_perm", [1000, 5000])
def test_gsea_rerun_at_b_199_and_sample_size_101(n_perm: int) -> None:
    """§8.6 B = 199: every engine call runs ``n_perm_simple = 199``, multilevel
    ``sample_size = 101``, whatever ``function.gsea.n_perm_simple`` says."""
    cfg = sf.base_config(function={"gsea": {"n_perm_simple": n_perm}})
    eng = sf.constant_engine(0.5)
    sf.call_ec(sf.example1_genes(), cfg, _rng(), sf.example1_collections(), eng)
    assert eng.calls, "GSEA was not run"
    assert {c.params.n_perm_simple for c in eng.calls} == {sf.SENS_B}
    assert {c.params.sample_size for c in eng.calls} == {sf.SAMPLE_SIZE}


# --------------------------------------------------------------------------- examples


def test_example_1_min_variant_row() -> None:
    """Example 1: primary {cog:K:gsea, go_bp:X:gsea}, ``missing=min`` {cog:K:gsea} -> n 2, 1,
    1, jaccard 0.5 (X's genes are imputed; under ``min`` they drop to the lowest abundance)."""
    S = _S()
    cfg = sf.base_config()
    genes, colls = sf.example1_genes(), sf.example1_collections()
    prim = sf.call_ec(genes, cfg, _rng(), colls, sf.example1_engine())
    vgen = S.variant_genes(genes, _variant("missing=min"), cfg)
    var = sf.call_ec(vgen, cfg, _rng(), colls, sf.example1_engine())
    assert prim == {"cog:K:gsea", "go_bp:X:gsea"}
    assert var == {"cog:K:gsea"}
    row = S.sensitivity_rows({sf.EC: prim}, {(sf.EC, "missing=min"): var}).iloc[0]
    assert (row["n_primary"], row["n_variant"], row["n_shared"]) == (2, 1, 1)
    assert row["jaccard"] == pytest.approx(0.5)
    assert pd.isna(row["agreement"])


def test_example_1_exclude_drops_the_imputed_set() -> None:
    """Under ``exclude`` X's genes have no abundance, so X has no ranked TU and is not tested."""
    S = _S()
    cfg = sf.base_config()
    genes, colls = sf.example1_genes(), sf.example1_collections()
    vgen = S.variant_genes(genes, _variant("missing=exclude"), cfg)
    assert sf.call_ec(vgen, cfg, _rng(), colls, sf.example1_engine()) == {"cog:K:gsea"}


def test_example_2_single_peptide_variant_reruns_ora_and_gsea() -> None:
    """Example 2 (``two_peptide``): the single-peptide variant reruns the ORA (cog:J moves from
    under- to over-detected) and GSEA (the engine is called on the variant table too)."""
    S = _S()
    cfg = sf.base_config()
    genes, colls = sf.example2_genes(), sf.example2_collections()
    vgen = S.variant_genes(genes, _variant("single_peptide_detected"), cfg)
    e1, e2 = sf.constant_engine(0.9), sf.constant_engine(0.9)
    prim = sf.call_ec(genes, cfg, _rng(), colls, e1)
    var = sf.call_ec(vgen, cfg, _rng(), colls, e2)
    assert prim == sf.ora_ids(genes, colls, cfg)
    assert var == sf.ora_ids(vgen, colls, cfg)
    assert "cog:J:ora_under" in prim and "cog:J:ora_over" in var
    assert prim != var
    assert e1.calls and e2.calls, "GSEA not rerun"


def test_example_4_same_seed_identical() -> None:
    """Example 4 with the M8-02 port: the same stream twice gives the same set (non-empty)."""
    from giap.function.gsea import SIMPLE_ENGINE

    cfg = sf.base_config()
    genes, colls = sf.example1_genes(), sf.example1_collections()
    a = sf.call_ec(genes, cfg, _rng(), colls, SIMPLE_ENGINE)
    b = sf.call_ec(genes, cfg, _rng(), colls, SIMPLE_ENGINE)
    assert a == b
    assert "cog:K:gsea" in a


def test_pure_inputs_unchanged() -> None:
    cfg = sf.base_config()
    genes, colls = sf.example1_genes(), sf.example1_collections()
    before = genes.copy()
    sf.call_ec(genes, cfg, _rng(), colls, sf.example1_engine())
    assert genes.equals(before)


# --------------------------------------------------------------------------- ORA invariance


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(data=st.data())
def test_ora_ids_unchanged_under_min_and_exclude(data: st.DataObject) -> None:
    """Lead decision 1: detection does not change under ``min`` / ``exclude``, so with a GSEA that
    flags nothing the family's set is the primary ORA set for both variants."""
    S = _S()
    cfg = sf.base_config()
    n = 90
    det = data.draw(st.lists(st.booleans(), min_size=n, max_size=n))
    p = data.draw(st.lists(st.sampled_from([0.2, 0.5, 0.8]), min_size=n, max_size=n))
    imp_draw = data.draw(st.lists(st.booleans(), min_size=n, max_size=n))
    cog = data.draw(st.lists(st.sampled_from(["J", "K", "L"]), min_size=n, max_size=n))
    imp = [d and i for d, i in zip(det, imp_draw, strict=True)]
    ab = [math.nan if (not d or i) else float(k % 9) for k, (d, i) in
          enumerate(zip(det, imp, strict=True))]  # fmt: skip
    obs = [a for a in ab if math.isfinite(a)]
    med = float(np.median(obs)) if obs else math.nan
    ab = [med if i else a for a, i in zip(ab, imp, strict=True)]
    imp = [i and bool(obs) for i in imp]
    genes = sf.gene_table(ab, det, p, imputed=imp, cog=cog)
    colls = [sf.gene_sets("cog", {c: [t for t, x in zip(genes["locus_tag"], cog, strict=True)
                                      if x == c] for c in sorted(set(cog))})]  # fmt: skip
    want = sf.ora_ids(genes, colls, cfg)
    prim = sf.call_ec(genes, cfg, _rng(), colls, sf.constant_engine(0.9))
    assert prim == want
    for name in ("missing=min", "missing=exclude"):
        g = S.variant_genes(genes, _variant(name), cfg)
        assert sf.call_ec(g, cfg, _rng(), colls, sf.constant_engine(0.9)) == want, name


# --------------------------------------------------------------------------- degenerate rows


def test_nothing_significant_gives_n_zero_and_jaccard_na() -> None:
    S = _S()
    cfg = sf.base_config()
    genes, colls = sf.null_genes(), sf.null_collections("cog", 3)
    prim = sf.call_ec(genes, cfg, _rng(), colls, sf.constant_engine(0.9))
    var = sf.call_ec(genes, cfg, _rng(), colls, sf.constant_engine(0.9))
    assert prim == var == frozenset()
    row = S.sensitivity_rows({sf.EC: prim}, {(sf.EC, "missing=min"): var}).iloc[0]
    assert (row["n_primary"], row["n_variant"], row["n_shared"]) == (0, 0, 0)
    assert pd.isna(row["jaccard"])


def test_no_collection_gives_empty_set() -> None:
    eng = sf.constant_engine(0.001)
    assert sf.call_ec(sf.null_genes(), sf.base_config(), _rng(), [], eng) == frozenset()


def test_no_abundance_keeps_ora_ids_without_gsea() -> None:
    """Without abundance GSEA has no row (``layer_missing:abundance``); the ORA still flags."""
    genes = sf.example2_genes()
    genes["abundance_log2"] = math.nan
    genes["abundance_scale"] = None
    cfg = sf.base_config()
    colls = sf.example2_collections()
    got = sf.call_ec(genes, cfg, _rng(), colls, sf.constant_engine(0.001))
    assert got == sf.ora_ids(genes, colls, cfg)
    assert not any(i.endswith(":gsea") for i in got)
