"""M4-14 review fixes: lead decisions 1-3 and 5 after code review (PR #82, changes requested).

docs/specs/M4-14-abundance.md "Lead decisions after code review":
1. Protein-route scale: ``gene_abundance(..., c10_scale=...)``; explicit ``abundance.scale`` (not
   ``auto``) first, then the preflight scale (``M4Inputs.abundance_scale``, passed as
   ``c10_scale``), else ValueError; mixed scales in C10 raise. C10 has NO ``scale`` column.
2. ``cross_strain_protein`` comes from the C12 protein-hit flags (M3-16 ``protein_hits``, passed as
   ``gene_abundance(..., protein_hits=...)``), not from the protein -> gene map (protein_id, gene,
   method only); the gene flag needs the flag on a C12 hit of a protein mapped by ``protein_hit``.
3. Kill the surviving review mutants (local/scratch/review-m4-14/mutants.py): M11 explicit scale
   wins, M24 mixed C10 scales raise, M12 the ``protein_hit`` condition, M17 +-inf never survives
   the log2 route, M14 ``no_observed_abundance`` logged at WARNING.
5. Detected partial genes with a missing value ARE filled; partial genes stay out of the fill
   statistic.

Resolved ambiguities (pinned here):
- "Mixed scales in C10": an (optional, non-contract) ``scale`` column carrying more than one value
  raises ValueError when ``abundance.scale = auto``, whatever ``c10_scale`` says.
- The scale is needed only on the protein route: peptide intensities, or no C10 at all, never
  raise for a missing scale.
- ``protein_hits`` flags are lists (M4Inputs parses the ``;`` TSV lists); a protein carries the
  flag when any of its C12 rows does; ``protein_hits`` omitted or empty means no such flag.
- "Logged at WARNING": a ``logging`` record at level WARNING whose message contains
  ``no_observed_abundance`` (a Python ``warnings.warn`` alone does not count).
"""

from __future__ import annotations

import logging
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import abundance_oracle as ao
from abundance_oracle import (
    acc,
    accepted_frame,
    c10,
    covariates_frame,
    detection,
    pmap_frame,
    protein_hits_frame,
)
from test_abundance_gene import TAGS, G, H, K, L, M, cfg, flags, genes_frame, na, row, run

SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])


def _protein_run(values: dict[str, float], owners: dict[str, str], **kw: Any) -> Any:
    pm = pmap_frame([(p, owners[p], "protein_hit") for p in values])
    kw.setdefault("detected", tuple(sorted(set(owners.values()))))
    return run(None, abundance=c10(values), pmap=pm, type="ibaq", **kw)


# ------------------------------------------------------------------ decision 1: scale precedence


def test_c10_fixture_has_no_scale_column() -> None:
    assert list(c10({"P1": 1.0}).columns) == ["protein_id", "value"]
    assert list(pmap_frame([("P1", G, "protein_hit")]).columns) == ["protein_id", "gene", "method"]


def test_preflight_scale_used_when_config_auto() -> None:
    res = _protein_run({"P1": 1024.0}, {"P1": G}, c10_scale="linear")
    assert row(res, G)["abundance_log2"] == pytest.approx(10.0)
    assert row(res, G)["abundance_raw"] == pytest.approx(1024.0)


@pytest.mark.parametrize(
    ("explicit", "preflight", "value", "log2"),
    [("linear", "log2", 1024.0, 10.0), ("log2", "linear", 20.0, 20.0)],
)
def test_explicit_abundance_scale_wins_over_preflight(
    explicit: str, preflight: str, value: float, log2: float
) -> None:
    """Mutant M11: ignoring an explicit ``abundance.scale`` reads the preflight scale instead."""
    res = _protein_run({"P1": value}, {"P1": G}, c10_scale=preflight, scale=explicit)
    assert row(res, G)["abundance_log2"] == pytest.approx(log2)
    assert row(res, G)["abundance_raw"] == pytest.approx(value)


def test_explicit_scale_without_preflight() -> None:
    res = _protein_run({"P1": 1024.0}, {"P1": G}, c10_scale=None, scale="linear")
    assert row(res, G)["abundance_log2"] == pytest.approx(10.0)


def test_no_scale_anywhere_raises() -> None:
    with pytest.raises(ValueError):
        _protein_run({"P1": 20.0}, {"P1": G}, c10_scale=None)


@pytest.mark.parametrize("preflight", ["log2", "linear", None])
def test_mixed_c10_scales_raise(preflight: str | None) -> None:
    """Mutant M24: a C10 frame mixing two measurement scales is accepted."""
    ab = c10({"P1": 20.0, "P2": 1024.0}).assign(scale=["log2", "linear"])
    pm = pmap_frame([("P1", G, "protein_hit"), ("P2", H, "protein_hit")])
    with pytest.raises(ValueError):
        run(None, detected=(G, H), abundance=ab, pmap=pm, type="ibaq", c10_scale=preflight)


def test_scale_not_needed_on_peptide_route() -> None:
    rows = [acc(G, "PEPGK", 1024.0)]
    res = run(accepted_frame(rows), n_obs={G: 1}, detected=(G,), abundance=c10({"P1": 5.0}),
              pmap=pmap_frame([("P1", G, "protein_hit")]), c10_scale=None)  # fmt: skip
    assert res.source == "peptide_intensity"
    assert row(res, G)["abundance_log2"] == pytest.approx(10.0)


def test_scale_not_needed_without_c10() -> None:
    res = run(None, detected=(G,), c10_scale=None)
    assert res.source == "none"


# ------------------------------------------------------- decision 2: cross_strain_protein from C12


def test_cross_strain_flag_needs_protein_hit_mapping_and_c12_flag() -> None:
    """Mutant M12: the flag is set whatever the mapping method."""
    ab = c10({"Q1": 21.0, "Q2": 20.0, "Q3": 19.0, "P4": 18.0})
    pm = pmap_frame(
        [
            ("Q1", G, "protein_hit"),  # protein_hit + C12 flag -> flag
            ("Q2", H, "peptide_vote"),  # C12 flag but mapped by peptide vote -> no flag
            ("Q3", K, "locus_tag"),  # C12 flag but mapped by locus tag -> no flag
            ("P4", L, "protein_hit"),  # protein_hit, C12 hit without the flag -> no flag
        ]
    )
    hits = protein_hits_frame(
        [
            ("Q1", G, ["cross_strain_protein"]),
            ("Q2", M, ["cross_strain_protein"]),
            ("Q3", K, ["cross_strain_protein", "two_proteins_one_gene"]),
            ("P4", L, []),
        ]
    )
    res = run(None, detected=(G, H, K, L), abundance=ab, pmap=pm, type="ibaq", c10_scale="log2",
              protein_hits=hits)  # fmt: skip
    assert "cross_strain_protein" in flags(row(res, G))
    for t in (H, K, L, M):
        assert "cross_strain_protein" not in flags(row(res, t)), t
    assert row(res, G)["abundance_log2"] == pytest.approx(21.0)
    assert row(res, H)["abundance_log2"] == pytest.approx(20.0)


def test_cross_strain_flag_on_a_lower_rank_c12_row_counts() -> None:
    ab = c10({"Q1": 21.0})
    pm = pmap_frame([("Q1", G, "protein_hit")])
    hits = protein_hits_frame([("Q1", G, [], 1), ("Q1", H, ["cross_strain_protein"], 2)])
    res = run(None, detected=(G,), abundance=ab, pmap=pm, type="ibaq", c10_scale="log2",
              protein_hits=hits)  # fmt: skip
    assert "cross_strain_protein" in flags(row(res, G))


def test_no_protein_hits_no_cross_strain_flag() -> None:
    res = _protein_run({"Q1": 21.0}, {"Q1": G}, c10_scale="log2", protein_hits=None)
    assert "cross_strain_protein" not in flags(row(res, G))
    assert row(res, G)["abundance_log2"] == pytest.approx(21.0)


def test_pmap_flags_column_is_not_read() -> None:
    """A stray ``flags`` column on the map (pre-decision fixture shape) gives no gene flag."""
    ab = c10({"Q1": 21.0})
    pm = pmap_frame([("Q1", G, "protein_hit")]).assign(flags=[["cross_strain_protein"]])
    res = run(None, detected=(G,), abundance=ab, pmap=pm, type="ibaq", c10_scale="log2")
    assert "cross_strain_protein" not in flags(row(res, G))


# ---------------------------------------------------------------- decision 3: +-inf, WARNING


@pytest.mark.parametrize("scale", ["log2", "linear"])
def test_infinite_c10_values_never_survive(scale: str) -> None:
    """Mutant M17: +-inf kept on the log2 route (and then poisons the fill statistic)."""
    good = 10.0 if scale == "log2" else 1024.0
    good2 = 14.0 if scale == "log2" else 16384.0
    values = {"P1": math.inf, "P2": -math.inf, "P3": good, "P4": good2, "P5": math.inf}
    owners = {"P1": G, "P2": H, "P3": K, "P4": L, "P5": M}
    res = _protein_run(values, owners, detected=(G, H, K, L), c10_scale=scale)
    log2 = res.table["abundance_log2"].astype(float)
    assert not np.isinf(log2).any(), res.table
    assert res.fill_value == pytest.approx(12.0)
    for t in (G, H):
        assert row(res, t)["abundance_log2"] == pytest.approx(12.0)
        assert bool(row(res, t)["abundance_imputed"])
    assert na(row(res, M)["abundance_log2"]), "undetected inf -> NA, not imputed"


def test_infinite_peptide_intensity_skipped() -> None:
    rows = [acc(G, "AAAAAAAK", math.inf), acc(G, "CCCCCCCK", 256.0, start=30)]
    res = run(accepted_frame(rows), n_obs={G: 1}, detected=(G,))
    assert row(res, G)["abundance_log2"] == pytest.approx(8.0)
    assert row(res, G)["top3_log2"] == pytest.approx(8.0)


@SETTINGS
@given(
    st.lists(
        st.one_of(
            st.floats(min_value=-1e6, max_value=1e12, allow_nan=False, allow_infinity=False),
            st.sampled_from([math.inf, -math.inf, math.nan, 0.0]),
        ),
        min_size=len(TAGS),
        max_size=len(TAGS),
    ),
    st.sampled_from(["log2", "linear"]),
    st.sampled_from(ao.POLICIES),
)
def test_protein_route_log2_always_finite_or_na(
    values: list[float], scale: str, policy: str
) -> None:
    vals = {f"P{i}": v for i, v in enumerate(values)}
    owners = {f"P{i}": t for i, t in enumerate(TAGS)}
    res = _protein_run(vals, owners, c10_scale=scale, policy=policy)
    log2 = res.table["abundance_log2"].astype(float).to_numpy()
    assert not np.isinf(log2).any()
    assert res.fill_value is None or math.isfinite(res.fill_value)
    for i, t in enumerate(TAGS):
        exp = ao.to_log2(values[i], scale)
        got = row(res, t)["abundance_log2"]
        if not bool(row(res, t)["abundance_imputed"]):
            assert (na(got) and na(exp)) or math.isclose(float(got), exp, abs_tol=1e-9)


def test_no_observed_abundance_logged_at_warning(caplog: pytest.LogCaptureFixture) -> None:
    """Mutant M14: the record is emitted below WARNING."""
    rows = [acc(H, "PEPHK", 4096.0)]
    with caplog.at_level(logging.WARNING):
        res = run(accepted_frame(rows), n_obs={H: 1}, detected=(G, K))
    assert res.fill_value is None
    hits = [
        r
        for r in caplog.records
        if r.levelno == logging.WARNING and "no_observed_abundance" in r.getMessage()
    ]
    assert hits, [(r.levelname, r.getMessage()) for r in caplog.records]


def test_no_warning_when_a_detected_gene_is_observed(caplog: pytest.LogCaptureFixture) -> None:
    rows = [acc(G, "PEPGK", 4096.0)]
    with caplog.at_level(logging.DEBUG):
        run(accepted_frame(rows), n_obs={G: 1}, detected=(G, K))
    assert not any("no_observed_abundance" in r.getMessage() for r in caplog.records)


# ----------------------------------------------------- decision 5: detected partial genes filled


def _partial_run(**kw: Any) -> Any:
    rows = [acc(t, f"PEP{t}K", 2.0**v) for t, v in ((G, 10), (H, 12), (L, 30))]
    genes = genes_frame(TAGS, partial=(L, M))
    det = detection(TAGS, (G, H, K, L, M))  # a stub marking the partial L and M detected
    from gippa.genetable.abundance import gene_abundance

    return gene_abundance(genes, accepted_frame(rows), covariates_frame({t: 1 for t in TAGS}, TAGS),
                          det, None, pmap_frame([]), cfg(**kw), protein_hits=protein_hits_frame([]),
                          c10_scale=None)  # fmt: skip


def test_detected_partial_gene_with_missing_value_is_filled() -> None:
    res = _partial_run()
    assert res.fill_value == pytest.approx(11.0), "median(10, 12): the partial L (30) is excluded"
    for t in (K, M):  # K eligible, M partial: both detected and missing -> filled
        r = row(res, t)
        assert r["abundance_log2"] == pytest.approx(11.0), t
        assert bool(r["abundance_imputed"]), t
        assert "imputed_abundance" in flags(r), t
    r = row(res, L)
    assert r["abundance_log2"] == pytest.approx(30.0)
    assert not bool(r["abundance_imputed"])


def test_detected_partial_gene_not_filled_under_exclude() -> None:
    res = _partial_run(missing="exclude")
    assert na(row(res, M)["abundance_log2"])
    assert not res.table["abundance_imputed"].astype(bool).any()


@st.composite
def _partial_scenarios(draw: st.DrawFn) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for t in TAGS:
        for j in range(draw(st.integers(min_value=0, max_value=3))):
            v = draw(st.floats(min_value=1.0, max_value=1e9, allow_nan=False))
            rows.append(acc(t, f"P{t}X{j}K", v, start=30 * j))
    return {
        "rows": rows,
        "detected": tuple(t for t in TAGS if draw(st.booleans())),
        "partial": tuple(t for t in TAGS if draw(st.booleans())),
        "policy": draw(st.sampled_from(ao.POLICIES)),
    }


@SETTINGS
@given(_partial_scenarios())
def test_partial_genes_filled_but_not_in_statistic_matches_oracle(sc: dict[str, Any]) -> None:
    from gippa.genetable.abundance import gene_abundance

    if not sc["rows"]:
        return  # source "none": covered elsewhere
    res = gene_abundance(
        genes_frame(TAGS, sc["partial"]),
        accepted_frame(sc["rows"]),
        covariates_frame({t: 1 for t in TAGS}, TAGS),
        detection(TAGS, sc["detected"]),  # the stub may mark partial genes detected
        None,
        pmap_frame([]),
        cfg(),
        policy=sc["policy"],
        c10_scale=None,
        protein_hits=protein_hits_frame([]),
    )
    vals = ao.peptide_values(accepted_frame(sc["rows"]))
    log2 = {t: (math.log2(ao.ibaq(vals[t], 1)) if t in vals else math.nan) for t in TAGS}
    det = {t: t in sc["detected"] for t in TAGS}
    elig = {t: t not in sc["partial"] for t in TAGS}
    filled, mask, fill = ao.impute(log2, det, sc["policy"], elig)
    assert (res.fill_value is None) == (fill is None)
    if fill is not None:
        assert math.isclose(res.fill_value, fill, abs_tol=1e-9)
    t = res.table.set_index("locus_tag")
    for g in TAGS:
        got, exp = t.at[g, "abundance_log2"], filled[g]
        assert (na(got) and na(exp)) or math.isclose(float(got), exp, abs_tol=1e-9), (g, got, exp)
        assert bool(t.at[g, "abundance_imputed"]) == mask[g], g
    assert pd.Series(mask).sum() == int(t["abundance_imputed"].astype(bool).sum())
