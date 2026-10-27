"""M2-19 Hypothesis properties (docs/specs/M2-19-regions.md acceptance: "``overlap_frac`` ∈
[0, 1]; equals brute-force base counting on random intervals of both topologies").

The spec names this file ``tests/property/test_region_membership.py``; it carries the
``regions_`` prefix used by every M2-19 test module.
"""

from __future__ import annotations

import importlib

import pytest
from helpers.regions_oracle import brute_members, brute_overlap_frac, oracle_subset
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.annotate.genomad import GenomadResult, Provirus
from gippa.config import resolve_config
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene

regions = importlib.import_module("gippa.annotate.regions")  # M2-19

SETTINGS = settings(max_examples=150, deadline=None)


@st.composite
def replicon_and_interval(
    draw: st.DrawFn, rep: Replicon | None = None
) -> tuple[Replicon, Interval]:
    if rep is None:
        length = draw(st.integers(1, 400))
        topo = draw(st.sampled_from(["circular", "linear"]))
        rep = Replicon("c", length, topo)  # type: ignore[arg-type]
    return rep, draw(interval_on(rep))


@st.composite
def interval_on(draw: st.DrawFn, rep: Replicon) -> Interval:
    start = draw(st.integers(0, rep.length - 1))
    max_end = rep.length if rep.topology == "linear" else start + rep.length
    end = draw(st.integers(start + 1, max_end))
    return Interval(rep.id, start, end, draw(st.sampled_from(["+", "-", "."])))


@st.composite
def rep_gene_region(draw: st.DrawFn) -> tuple[Replicon, Interval, Interval]:
    rep, g = draw(replicon_and_interval())
    return rep, g, draw(interval_on(rep))


@SETTINGS
@given(rep_gene_region())
def test_overlap_frac_in_unit_interval_and_equals_base_counting(
    case: tuple[Replicon, Interval, Interval],
) -> None:
    rep, g, r = case
    f = regions.overlap_frac(g, r, rep)
    assert 0.0 <= f <= 1.0
    assert f == pytest.approx(brute_overlap_frac(g, r, rep), abs=1e-12)


@SETTINGS
@given(replicon_and_interval())
def test_whole_replicon_contains_every_gene(case: tuple[Replicon, Interval]) -> None:
    rep, g = case
    assert regions.overlap_frac(g, Interval(rep.id, 0, rep.length, "."), rep) == pytest.approx(1)


@st.composite
def genes_region_frac(draw: st.DrawFn) -> tuple[Replicon, list[Gene], Interval, float]:
    rep, region = draw(replicon_and_interval())
    ivs = draw(st.lists(interval_on(rep), max_size=12))
    genes = [Gene(f"g{k}", "cds", v, "M", None, "ATG") for k, v in enumerate(ivs)]
    frac = draw(st.sampled_from([0.0001, 0.25, 0.5, 0.75, 1.0]))
    return rep, genes, region, frac


@SETTINGS
@given(genes_region_frac())
def test_members_are_exactly_the_genes_at_or_above_min_frac(
    case: tuple[Replicon, list[Gene], Interval, float],
) -> None:
    rep, genes, region, frac = case
    got = dict(regions.members(genes, region, rep, min_frac=frac))
    want = brute_members(genes, region, rep, frac)
    assert set(got) == set(want)
    for tag, f in got.items():
        assert f == pytest.approx(want[tag], abs=1e-12)


_WORDS = st.sampled_from(
    [
        "phage", "protein", "integrase", "Cro", "cro", "cI", "CI", "repressor", "tail", "fiber",
        "holin", "lysozyme", "Rz", "micro", "circ", "DNA", "methylase", "capsid", "portal",
        "excisionase", "hypothetical", "family", "domain-containing", "tape", "measure", "head",
        "-like", "/", "scaffold", "terminase", "spanin", "endolysin", "antirepressor",
    ]
)  # fmt: skip
_TEXT = st.lists(_WORDS, max_size=6).map(" ".join)


@SETTINGS
@given(_TEXT, st.one_of(st.none(), _TEXT))
def test_subset_matches_the_independent_keyword_oracle(product: str, marker: str | None) -> None:
    got = regions.prophage_subset(product, marker)
    assert got in {"integrase", "repressor", "structural", "lysis", "cargo"}
    assert got == oracle_subset(product, marker)


@SETTINGS
@given(_TEXT, st.one_of(st.none(), _TEXT))
def test_integrase_wins_over_every_later_rule(product: str, marker: str | None) -> None:
    assert regions.prophage_subset(f"{product} integrase", marker) == "integrase"


@SETTINGS
@given(
    st.lists(st.tuples(st.integers(0, 99_000), st.integers(1, 5000)), min_size=1, max_size=6),
    st.randoms(use_true_random=False),
)
def test_prophage_ids_follow_start_order_for_any_input_order(spans, rnd) -> None:  # type: ignore[no-untyped-def]
    rep = Replicon("chr1", 100_000, "circular", "chromosome")
    pvs = [Provirus("chr1", Interval("chr1", s, s + n, "."), 0.9, None, 1, "x") for s, n in spans]
    shuffled = list(pvs)
    rnd.shuffle(shuffled)
    cfg = resolve_config(overrides=["context.region_types=[prophage]"])
    regs, _ = regions.build_regions([], [], {"chr1": rep}, GenomadResult(tuple(shuffled), (), ()),
                                    None, cfg)  # fmt: skip
    assert sorted(r.region_id for r in regs) == sorted(
        f"prophage_{k}" for k in range(1, len(pvs) + 1)
    )
    starts = [r.interval.start for r in sorted(regs, key=lambda r: int(r.region_id.split("_")[1]))]
    assert starts == sorted(starts)
