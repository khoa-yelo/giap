"""M4-19 sensitivity framework and the ``regions`` family: spec examples 1-5 and 7, edge cases and
lead decisions 2-4 (docs/specs/M4-19-sensitivity-regions.md; GIAP_SPEC.md §8.6, §5.2, §7.3.2).

Examples 6 (``off``: no file) and the step-level CRN checks live in
tests/integration/test_sensitivity_step.py. Inputs: hand-built §9.3 gene tables
(tests/helpers/sensitivity_data.py); the reported flagged set comes from merged M4-17
``test_all_regions`` with the lead-decision-4 rule applied in the helper (``flagged_from``).
"""

from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import sensitivity_data as sd

NAN = math.nan


def _S() -> Any:
    from giap.context import sensitivity as S

    return S


def _names(vs: list[Any]) -> list[str]:
    return [v.name for v in vs]


def _is_na(x: Any) -> bool:
    return x is None or x is pd.NA or (isinstance(x, float) and math.isnan(x))


# --------------------------------------------------------------------------- interface


def test_variant_is_a_frozen_slotted_dataclass() -> None:
    S = _S()
    assert dataclasses.is_dataclass(S.Variant)
    assert [f.name for f in dataclasses.fields(S.Variant)] == [
        "name", "overrides", "single_peptide_detected", "config_patch",  # + M6-10 decision 1
        "families",  # + M6-10 decision 8 (review B2)
    ]  # fmt: skip
    assert S.Variant.__dataclass_params__.frozen
    assert hasattr(S.Variant, "__slots__")
    v = S.Variant("missing=min", {"abundance.missing": "min"})
    assert v.single_peptide_detected is False
    with pytest.raises(dataclasses.FrozenInstanceError):
        v.name = "x"  # type: ignore[misc]


def test_families_registers_regions_only_in_m4() -> None:
    """M4-index decision 7: M4 registers ``regions``; M5-M8 add theirs (M5-10: ``islands``, M6-10:
    ``operon_flags``, M8-16: ``enriched_categories`` and ``determinants``)."""
    S = _S()
    known = {"regions", "islands", "operon_flags", "enriched_categories", "determinants"}
    assert known >= set(S.FAMILIES) >= {"regions"}
    assert callable(S.FAMILIES["regions"])


# --------------------------------------------------------------------------- variants


def test_example_3_basic_mayu_median() -> None:
    S = _S()
    vs = S.variants(sd.config(), "mayu")
    assert _names(vs) == ["missing=min", "missing=exclude"]
    assert not any(v.single_peptide_detected for v in vs)


def test_missing_variant_overrides_the_missing_policy() -> None:
    S = _S()
    for v in S.variants(sd.config(), "mayu"):
        policy = v.name.removeprefix("missing=")
        keys = [k for k, val in v.overrides.items() if str(k).endswith("missing") and val == policy]
        assert keys, (v.name, dict(v.overrides))


@pytest.mark.parametrize("tier", ["decoy_picked", "mayu", "protein_list"])
def test_no_single_peptide_variant_outside_two_peptide(tier: str) -> None:
    S = _S()
    assert "single_peptide_detected" not in _names(S.variants(sd.config(), tier))


def test_example_4_two_peptide_adds_single_peptide_detected() -> None:
    S = _S()
    vs = S.variants(sd.config(), "two_peptide")
    assert _names(vs) == ["missing=min", "missing=exclude", "single_peptide_detected"]
    sp = vs[-1]
    assert sp.single_peptide_detected is True
    assert not any(v.single_peptide_detected for v in vs[:-1])


def test_example_5_primary_policy_is_not_rerun() -> None:
    S = _S()
    vs = S.variants(sd.config(abundance={"missing": "min"}), "mayu")
    assert _names(vs) == ["missing=exclude"]


def test_variants_follow_abundance_sensitivity_order() -> None:
    S = _S()
    cfg = sd.config(abundance={"sensitivity": ["exclude", "zero", "min"]})
    assert _names(S.variants(cfg, "mayu")) == ["missing=exclude", "missing=zero", "missing=min"]


def test_off_gives_no_variants() -> None:
    S = _S()
    for tier in ("mayu", "two_peptide"):
        assert S.variants(sd.config(output={"sensitivity": "off"}), tier) == []


def test_full_adds_top3_and_the_uniop_threshold_variants() -> None:
    """``full`` adds ``top3`` and (M6-10 lead decision 1) the UniOP thresholds 0.3 / 0.7 as
    config-patch variants."""
    S = _S()
    vs = S.variants(sd.config(output={"sensitivity": "full"}), "two_peptide")
    assert _names(vs) == ["missing=min", "missing=exclude", "single_peptide_detected", "top3",
                          "pair_threshold=0.3", "pair_threshold=0.7"]  # fmt: skip
    top3 = vs[3]
    assert any(str(k).endswith("method") and v == "top3" for k, v in top3.overrides.items())
    assert all(not v.config_patch for v in vs[:4])


def test_basic_has_no_top3() -> None:
    S = _S()
    assert "top3" not in _names(S.variants(sd.config(), "mayu"))


def test_full_with_top3_primary_does_not_rerun_top3() -> None:
    """As example 5: a variant equal to the primary choice is not rerun."""
    S = _S()
    cfg = sd.config(output={"sensitivity": "full"}, abundance={"method": "top3"})
    assert "top3" not in _names(S.variants(cfg, "mayu"))


def test_variant_names_are_unique() -> None:
    S = _S()
    cfg = sd.config(output={"sensitivity": "full"})
    names = _names(S.variants(cfg, "two_peptide"))
    assert len(names) == len(set(names))


# --------------------------------------------------------------------------- jaccard


def test_jaccard_both_empty_is_none() -> None:
    assert _S().jaccard(frozenset(), frozenset()) is None


@pytest.mark.parametrize(
    ("a", "b", "want"),
    [
        ({"x", "y"}, {"x"}, 0.5),
        ({"x"}, set(), 0.0),
        (set(), {"x"}, 0.0),
        ({"x", "y"}, {"y", "x"}, 1.0),
        ({"a", "b", "c"}, {"b", "c", "d"}, 0.5),
        ({"a"}, {"b"}, 0.0),
    ],
)
def test_jaccard_values(a: set[str], b: set[str], want: float) -> None:
    got = _S().jaccard(frozenset(a), frozenset(b))
    assert isinstance(got, float)
    assert got == pytest.approx(want, abs=1e-15)


# --------------------------------------------------------------------------- sensitivity_rows


def test_example_1_one_row() -> None:
    S = _S()
    a = sd.family_id("prophage_1", "chr1", "o_minus_e")
    b = sd.family_id("bgc_1", "chr1", "hl_shift")
    df = S.sensitivity_rows({"regions": frozenset({a, b})},
                            {("regions", "missing=min"): frozenset({a})})  # fmt: skip
    assert isinstance(df, pd.DataFrame)
    assert tuple(df.columns[:6]) == sd.COLUMNS
    assert len(df) == 1
    r = df.iloc[0]
    assert (r["family"], r["variant"]) == ("regions", "missing=min")
    assert (int(r["n_primary"]), int(r["n_variant"]), int(r["n_shared"])) == (2, 1, 1)
    assert float(r["jaccard"]) == pytest.approx(0.5)


def test_example_2_both_empty_gives_na_and_zero_counts() -> None:
    S = _S()
    df = S.sensitivity_rows({"regions": frozenset()},
                            {("regions", "missing=min"): frozenset()})  # fmt: skip
    r = df.iloc[0]
    assert (int(r["n_primary"]), int(r["n_variant"]), int(r["n_shared"])) == (0, 0, 0)
    assert _is_na(r["jaccard"])


def test_rows_one_per_rerun_and_agreement_na_for_regions() -> None:
    S = _S()
    x = frozenset({"r:c:o_minus_e"})
    reruns = {("regions", "missing=min"): x, ("regions", "missing=exclude"): frozenset()}
    df = S.sensitivity_rows({"regions": x}, reruns)
    assert sorted(zip(df["family"], df["variant"], strict=True)) == sorted(reruns)
    got = dict(zip(df["variant"], df["jaccard"], strict=True))
    assert float(got["missing=min"]) == 1.0
    assert float(got["missing=exclude"]) == 0.0
    if "agreement" in df.columns:  # §8.6: only for the dominant determinant block
        assert list(df.columns) == [*sd.COLUMNS, "agreement"]
        assert df["agreement"].map(_is_na).all()
    else:
        assert list(df.columns) == list(sd.COLUMNS)
    assert pd.api.types.is_integer_dtype(df["n_primary"])
    assert pd.api.types.is_integer_dtype(df["n_shared"])


def test_rows_are_deterministic() -> None:
    S = _S()
    p = {"regions": frozenset({"a:c:o_minus_e", "b:c:hl_shift"})}
    r = {("regions", "missing=min"): frozenset({"a:c:o_minus_e"}),
         ("regions", "missing=exclude"): frozenset({"b:c:hl_shift", "z:c:hl_shift"})}  # fmt: skip
    pd.testing.assert_frame_equal(S.sensitivity_rows(p, r), S.sensitivity_rows(p, r))


# --------------------------------------------------------------------------- regions family:
# flagged-set rule (lead decisions 3 and 4) on a stubbed M4-17 result


def _stub_rows() -> list[tuple[str, str, str, str, float, float, str]]:
    cs = "circular_shift"
    return [
        ("prophage_1", "prophage", "chr1", "o_minus_e", 0.004, 0.01, cs),  # flagged
        ("prophage_1", "prophage", "chr1", "hl_shift", 0.04, 0.05, cs),  # q == alpha: flagged
        ("bgc_1", "bgc", "chr1", "o_minus_e", 0.03, 0.0500001, cs),  # just above
        ("bgc_1", "bgc", "chr1", "hl_shift", 0.5, 0.6, cs),
        ("amr_1", "amr", "chr1", "o_minus_e", 0.001, 0.002, cs),  # split region, part 1
        ("amr_1", "amr", "p1", "o_minus_e", 0.001, 0.002, cs),  # split region, part 2
        ("amr_1", "amr", "chr1", "hl_shift", NAN, NAN, cs),  # p NA
        ("amr_1", "amr", "p1", "hl_shift", NAN, 0.001, cs),  # p NA, q set: never flagged
        ("crispr_1", "crispr", "chr1", "o_minus_e", 0.001, 0.001, "none"),  # not a shift test
        ("crispr_1", "crispr", "chr1", "hl_shift", 0.2, 0.2, cs),
    ]


def test_flagged_set_rule_and_id_format(monkeypatch: pytest.MonkeyPatch) -> None:
    _S()
    sd.patch_region_tests(monkeypatch, sd.fake_tests(_stub_rows()))
    got = sd.call_regions_family(sd.family_tables(), sd.config(), np.random.default_rng(1))
    assert isinstance(got, frozenset)
    assert got == {
        "prophage_1:chr1:o_minus_e",
        "prophage_1:chr1:hl_shift",
        "amr_1:chr1:o_minus_e",
        "amr_1:p1:o_minus_e",
    }


def test_split_contig_parts_do_not_collide(monkeypatch: pytest.MonkeyPatch) -> None:
    """Lead decision 3: one ID per (region, contig, stat); ``<region_id>:<stat>`` would merge
    the two parts of a ``split_contigs`` region into one."""
    _S()
    cs = "circular_shift"
    rows = [("amr_1", "amr", "chr1", "o_minus_e", 0.001, 0.002, cs),
            ("amr_1", "amr", "p1", "o_minus_e", 0.001, 0.002, cs)]  # fmt: skip
    sd.patch_region_tests(monkeypatch, sd.fake_tests(rows))
    got = sd.call_regions_family(sd.family_tables(), sd.config(), np.random.default_rng(1))
    assert len(got) == 2
    assert got == {"amr_1:chr1:o_minus_e", "amr_1:p1:o_minus_e"}


def test_p_na_rows_are_never_flagged(monkeypatch: pytest.MonkeyPatch) -> None:
    _S()
    cs = "circular_shift"
    rows = [("bgc_1", "bgc", "chr1", "hl_shift", NAN, 0.0, cs),
            ("bgc_1", "bgc", "chr1", "o_minus_e", NAN, NAN, cs)]  # fmt: skip
    sd.patch_region_tests(monkeypatch, sd.fake_tests(rows))
    got = sd.call_regions_family(sd.family_tables(), sd.config(), np.random.default_rng(1))
    assert got == frozenset()


def test_alpha_context_is_the_threshold(monkeypatch: pytest.MonkeyPatch) -> None:
    _S()
    sd.patch_region_tests(monkeypatch, sd.fake_tests(_stub_rows()))
    t = sd.family_tables()
    loose = sd.call_regions_family(t, sd.config(stats={"alpha": {"context": 0.10}}),
                                   np.random.default_rng(1))  # fmt: skip
    assert "bgc_1:chr1:o_minus_e" in loose
    strict = sd.call_regions_family(t, sd.config(stats={"alpha": {"context": 0.005}}),
                                    np.random.default_rng(1))  # fmt: skip
    assert strict == {"amr_1:chr1:o_minus_e", "amr_1:p1:o_minus_e"}


def test_family_runs_the_region_tests_with_the_given_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The step sets ``stats.permutations = 199``; the family must not override it."""
    _S()
    calls = sd.patch_region_tests(monkeypatch, sd.fake_tests(_stub_rows()))
    cfg = sd.config(stats={"permutations": sd.SENS_B})
    sd.call_regions_family(sd.family_tables(), cfg, np.random.default_rng(1))
    assert len(calls) == 1
    passed = [a for a in (*calls[0]["args"], *calls[0]["kwargs"].values())
              if type(a).__name__ == "GiapConfig"]  # fmt: skip
    assert passed and all(c.stats.permutations == sd.SENS_B for c in passed)


# --------------------------------------------------------------------------- regions family:
# real M4-17 tests on the hand-built tables


def test_recomputed_primary_equals_reported() -> None:
    """Lead decision 2: with ``exhaustive_shift_max_genes = 20000`` every replicon here uses all
    shifts, so the primary recomputed at B = 199 is the reported (B = 1,000) flagged set."""
    t = sd.family_tables()
    reported = sd.reported_flagged(t, sd.config(), np.random.default_rng(11))
    assert reported == {
        "amr_1:chr1:o_minus_e", "amr_1:p1:o_minus_e",
        "prophage_1:chr1:o_minus_e", "prophage_1:chr1:hl_shift",
    }  # fmt: skip
    cfg = sd.config(stats={"permutations": sd.SENS_B})
    for seed in (1, 2, 3):
        assert sd.call_regions_family(t, cfg, np.random.default_rng(seed)) == reported


def test_example_7_exclude_keeps_o_minus_e_flags() -> None:
    """Example 7: ``missing=exclude`` changes only abundances: the O - E IDs equal the
    primary's; here the prophage HL flag goes (5 of its 10 abundances become NA)."""
    S = _S()
    t = sd.family_tables()
    cfg = sd.config(stats={"permutations": sd.SENS_B})
    prim = sd.call_regions_family(t, cfg, np.random.default_rng(5))
    excl_genes = sd.with_missing_abundance(t.genes, sd.tags_of(t.genes, "chr1", [0, 2, 4, 6, 8]))
    excl = sd.call_regions_family(dataclasses.replace(t, genes=excl_genes), cfg,
                                  np.random.default_rng(5))  # fmt: skip

    def oe(s: frozenset[str]) -> set[str]:
        return {x for x in s if x.endswith(":o_minus_e")}

    assert oe(prim) == oe(excl) == {"amr_1:chr1:o_minus_e", "amr_1:p1:o_minus_e",
                                    "prophage_1:chr1:o_minus_e"}  # fmt: skip
    assert "prophage_1:chr1:hl_shift" in prim
    assert "prophage_1:chr1:hl_shift" not in excl
    assert S.jaccard(prim, excl) == pytest.approx(3 / 4)


def test_exclude_with_nothing_imputed_is_identical() -> None:
    """Edge case: ``exclude`` with every detected gene observed equals the primary."""
    S = _S()
    t = sd.family_tables()
    cfg = sd.config(stats={"permutations": sd.SENS_B})
    prim = sd.call_regions_family(t, cfg, np.random.default_rng(5))
    same = sd.call_regions_family(dataclasses.replace(t, genes=sd.with_missing_abundance(
        t.genes, [])), cfg, np.random.default_rng(5))  # fmt: skip
    assert prim == same
    assert S.jaccard(prim, same) == 1.0


def test_sampled_shifts_same_stream_same_set() -> None:
    """Above ``exhaustive_shift_max_genes`` the shifts are drawn: the same generator state gives
    the same flagged set (what common random numbers rely on)."""
    t = sd.family_tables()
    cfg = sd.config(stats={"permutations": 60}, context={"exhaustive_shift_max_genes": 50})
    a = sd.call_regions_family(t, cfg, np.random.default_rng(99))
    b = sd.call_regions_family(t, cfg, np.random.default_rng(99))
    assert a == b


def test_region_types_off_give_no_ids() -> None:
    t = sd.family_tables()
    cfg = sd.config(stats={"permutations": sd.SENS_B}, context={"region_types": ["bgc", "crispr"]})
    got = sd.call_regions_family(t, cfg, np.random.default_rng(1))
    assert not any(x.startswith(("prophage_", "amr_")) for x in got)


# --------------------------------------------------------------------------- lead decisions 5, 13:
# the family uses the bound C7 / C8 / C1, never tables derived from ``region_ids``


@pytest.mark.parametrize(
    "given", [(), ("regions",), ("regions", "region_genes"), ("region_genes", "replicons")]
)
def test_regions_family_requires_bound_tables(given: tuple[str, ...]) -> None:
    """Decision 13: ``regions``, ``region_genes`` and ``replicons`` are keyword-only and required
    (no fallback to tables derived from the gene table)."""
    S = _S()
    t = sd.family_tables()
    tables = {"regions": t.regions, "region_genes": t.region_genes, "replicons": t.replicons}
    with pytest.raises(TypeError):
        S.regions_family(t.genes, sd.config(), np.random.default_rng(1),
                         **{k: tables[k] for k in given})  # fmt: skip


def test_regions_family_tables_are_keyword_only() -> None:
    import inspect

    params = inspect.signature(_S().regions_family).parameters
    for name in ("regions", "region_genes", "replicons"):
        assert params[name].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert params[name].default is inspect.Parameter.empty, name


@pytest.mark.parametrize("topology", ["circular", "linear"])
def test_bound_tables_and_topology_are_used(topology: str) -> None:
    """M27: C7 / C8 hold ``bgc_2`` that no gene's ``region_ids`` lists, and C1 makes p1
    ``topology``: the flagged set is the one of the M4-17 tests on exactly those tables (on the
    linear p1 ``bgc_2`` has 5 valid shifts and is never flagged; circular: p = 0.01)."""
    t = sd.bound_tables(topology)
    cfg = sd.config(stats={"permutations": sd.SENS_B})
    want = sd.reported_flagged(t, sd.config(), np.random.default_rng(3))
    got = sd.call_regions_family(t, cfg, np.random.default_rng(4))
    assert got == want
    assert ("bgc_2:p1:o_minus_e" in got) is (topology == "circular")
    assert "prophage_1:chr1:o_minus_e" in got  # region_ids are empty: C7 / C8 were used
