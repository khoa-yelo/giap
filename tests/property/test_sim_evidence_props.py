"""Hypothesis properties of giapsim evidence (issue M4-09; GIAP_SPEC.md §11.3, §7.2.6, §7.3.2).

For any small layout (plasmids, linear contigs, origin crossing), any expression seed and any valid
``EvidenceConfig`` (pi_wrong, false source, decoy count, emit subset, protein IDs, uniqueness
knobs, proteins-only, abundance table): the C9/C11/C12 frames pass ``validate_table`` and the
written work dir loads through ``load_m4_inputs`` without errors; C9 has one row per peptide;
decoys are never placed; false-target and matched-decoy counts follow the spec formulas; PEP is
the mixture posterior and q the mean target PEP at or above the score (so mean PEP of rows with
q <= t is <= t for every t); ``gene_truth`` agrees with the truth; placements obey §6.5 coordinate
rules; the frames are a pure function of (layout, expression, config, seed).
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import evidence_oracle as eo

SETTINGS = settings(max_examples=25, deadline=None,
                    suppress_health_check=[HealthCheck.too_slow])  # fmt: skip
EMIT_FIELDS = ("q_value", "pep", "score", "is_decoy", "intensity", "protein_ids")


def _mods() -> tuple[Any, Any, Any]:
    import giapsim.evidence as ev
    import giapsim.expression as ex
    import giapsim.layout as la

    return ev, ex, la


@st.composite
def scenarios(draw: st.DrawFn) -> tuple[Any, Any, Any, int]:
    ev, ex, la = _mods()
    lay = la.simulate_layout(
        la.LayoutConfig(
            n_genes_chromosome=draw(st.integers(80, 300)),
            plasmids=tuple(draw(st.lists(st.integers(1, 30), max_size=1))),
            linear_contigs=tuple(draw(st.lists(st.integers(1, 12), max_size=1))),
            frac_sorf=draw(st.sampled_from([0.02, 0.3])),
            origin_crossing=draw(st.booleans()),
        ),
        np.random.default_rng(draw(st.integers(0, 2**31))),
    )
    expr = ex.simulate_expression(lay, ex.ExpressionConfig(),
                                  np.random.default_rng(draw(st.integers(0, 2**31))))  # fmt: skip
    kw: dict[str, Any] = {
        "pi_wrong": draw(st.sampled_from([0.0, 0.01, 0.05, 0.2])),
        "false_source": draw(st.sampled_from(["random_gene", "unplaceable"])),
        "n_decoys": draw(st.one_of(st.just("matched"), st.integers(0, 30))),
        "mixture": draw(st.sampled_from([(0.0, 3.0), (-1.0, 2.0), (0.5, 4.0)])),
        "emit": frozenset(draw(st.sets(st.sampled_from(EMIT_FIELDS)))),
        "protein_ids": draw(st.sampled_from(["locus_tag", "foreign"])),
        "multi_identical_pairs": draw(st.integers(0, 2)),
        "shared_frac": draw(st.sampled_from([0.0, 0.1])),
        "tolerant_frac": draw(st.sampled_from([0.0, 0.1])),
        "proteins_only": draw(st.booleans()),
        "abundance_table": draw(st.sampled_from(["none", "ibaq", "maxlfq"])),
    }
    if draw(st.booleans()):  # half the scenarios use the full emit set (decoys, scores)
        kw["emit"] = frozenset(EMIT_FIELDS)
    return lay, expr, ev.EvidenceConfig(**kw), draw(st.integers(0, 2**31))


def _run(lay: Any, expr: Any, cfg: Any, seed: int) -> Any:
    ev, _, _ = _mods()
    return ev.simulate_evidence(lay, expr, cfg, np.random.default_rng(seed))


@SETTINGS
@given(scenarios())
def test_contract_valid_and_loadable(sc: tuple[Any, Any, Any, int]) -> None:
    from giap.genetable.contract import validate_inputs, validate_table
    from giapsim.contract_tables import write_contract_tables

    lay, expr, cfg, seed = sc
    ev = _run(lay, expr, cfg, seed)
    assert eo.errors(validate_table("peptides", ev.peptides)) == []
    assert eo.errors(validate_table("peptide_hits", ev.peptide_hits)) == []
    if ev.protein_hits is not None:
        assert eo.errors(validate_table("protein_hits", ev.protein_hits)) == []
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        write_contract_tables(lay, work)
        _mods()[0].write_evidence_tables(ev, work)
        inputs, issues = eo.load(work)
        assert eo.errors(issues) == []
        assert validate_inputs(inputs) == []


@SETTINGS
@given(scenarios())
def test_counts_uniqueness_and_decoys(sc: tuple[Any, Any, Any, int]) -> None:
    lay, expr, cfg, seed = sc
    ev = _run(lay, expr, cfg, seed)
    p, h = ev.peptides, ev.peptide_hits
    tags = set(lay.genes["locus_tag"])
    assert p["peptide"].is_unique
    assert set(h["gene"].astype(str)) <= tags
    if cfg.proteins_only:
        assert len(p) == 0 and len(h) == 0
        return
    dec = eo.decoys(p)
    assert not set(dec["peptide"]) & set(h["peptide"])
    if "is_decoy" not in cfg.emit:
        assert len(dec) == 0
    tgt = eo.targets(p)
    bad = eo.false_set(ev.truth) & set(tgt["peptide"])
    n_true = len(tgt) - len(bad)
    assert len(bad) == round(cfg.pi_wrong / (1 - cfg.pi_wrong) * n_true)
    if "is_decoy" in cfg.emit:
        want = len(bad) if cfg.n_decoys == "matched" else cfg.n_decoys
        assert len(dec) == want
    if cfg.false_source == "unplaceable":
        assert not bad & set(h["peptide"])
    else:
        assert bad <= set(h["peptide"])
    for col in ("q_value", "pep", "score", "intensity"):
        if col not in cfg.emit:
            assert p[col].isna().all(), col
    for col in ("q_value", "pep"):
        if col in cfg.emit:
            v = p[col].to_numpy(float)
            assert ((v >= 0) & (v <= 1)).all(), col


@SETTINGS
@given(scenarios(), st.floats(0.0, 0.3))
def test_pep_and_q_rules(sc: tuple[Any, Any, Any, int], t: float) -> None:
    lay, expr, cfg, seed = sc
    if cfg.proteins_only or not {"q_value", "pep", "score"} <= cfg.emit:
        return
    ev = _run(lay, expr, cfg, seed)
    p = ev.peptides
    if len(p) == 0:
        return
    mu0, mu1 = cfg.mixture
    np.testing.assert_allclose(p["pep"].to_numpy(float),
                               eo.pep(p["score"].to_numpy(float), cfg.pi_wrong, mu0, mu1),
                               atol=1e-6)  # fmt: skip
    tg = eo.targets(p)
    q = eo.q_from_targets(tg["score"], tg["score"], tg["pep"])
    np.testing.assert_allclose(tg["q_value"].to_numpy(float), q, atol=1e-6)
    acc = tg[tg["q_value"].to_numpy(float) <= t]
    if len(acc):
        assert float(acc["pep"].astype(float).mean()) <= t + 1e-9


@SETTINGS
@given(scenarios())
def test_gene_truth_and_detection(sc: tuple[Any, Any, Any, int]) -> None:
    lay, expr, cfg, seed = sc
    ev = _run(lay, expr, cfg, seed)
    gt = ev.gene_truth.set_index("locus_tag")
    assert sorted(gt.index) == sorted(lay.genes["locus_tag"])
    aa = dict(zip(lay.genes["locus_tag"], lay.genes["aa"], strict=True))
    det = dict(zip(expr.genes["locus_tag"], expr.genes["detected_true"], strict=True))
    for g, row in gt.iterrows():
        n = int(row["n_true_peptides"])
        assert bool(row["present"]) == (n >= 1)
        if not cfg.proteins_only:
            assert bool(row["present"]) == (bool(det[g]) and eo.digestible(aa[g])), g


@SETTINGS
@given(scenarios())
def test_placement_coordinates_follow_section_6_5(sc: tuple[Any, Any, Any, int]) -> None:
    lay, expr, cfg, seed = sc
    ev = _run(lay, expr, cfg, seed)
    length = dict(zip(lay.replicons["contig"], lay.replicons["length"], strict=True))
    for r in ev.peptide_hits.to_dict("records"):
        n, s, e = int(length[r["contig"]]), int(r["start"]), int(r["end"])
        assert 0 <= s < e
        assert e - s == 3 * len(r["peptide"])
        if bool(r["crosses_origin"]):
            assert s < n < e < 2 * n
        else:
            assert e <= n


@SETTINGS
@given(scenarios())
def test_pure_function_of_seed(sc: tuple[Any, Any, Any, int]) -> None:
    lay, expr, cfg, seed = sc
    a, b = _run(lay, expr, cfg, seed), _run(lay, expr, cfg, seed)
    for name in ("peptides", "peptide_hits", "truth", "gene_truth"):
        pd.testing.assert_frame_equal(getattr(a, name), getattr(b, name), obj=name)
