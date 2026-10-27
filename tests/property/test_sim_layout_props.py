"""Hypothesis properties of gippasim's sequence-free layout (issue M4-07; GIPPA_SPEC.md §11.3).

For any replicon mix (chromosome of 2-250 genes, 0-2 plasmids, 0-2 linear contigs, with or without
the default regions) and seed: the written tables pass the M4-06 contract with no error; genes,
TUs, gaps, locus tags and layers satisfy the spec (tests/helpers/layout_oracle.py); exactly the
first and last gene of each linear contig are partial; every adjacent same-strand pair has a
UniOP-like p; the layout is a pure function of (config, seed).
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
import layout_oracle as lo

SETTINGS = settings(
    max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow]
)  # fmt: skip


@st.composite
def configs(draw: st.DrawFn) -> dict[str, Any]:
    with_regions = draw(st.booleans())
    n = draw(st.integers(150, 250) if with_regions else st.integers(2, 120))
    kw: dict[str, Any] = {
        "n_genes_chromosome": n,
        "plasmids": tuple(draw(st.lists(st.integers(1, 15), max_size=2))),
        "linear_contigs": tuple(draw(st.lists(st.integers(1, 15), max_size=2))),
        "mean_tu_size": draw(st.sampled_from([1.0, 1.8, 3.0])),
        "p_same_strand": draw(st.sampled_from([0.0, 0.6, 1.0])),
    }
    if not with_regions:
        kw["regions"] = ()
    return kw


def _run(kw: dict[str, Any], seed: int, root: str) -> tuple[Any, dict[str, pd.DataFrame], Path]:
    from gippasim.contract_tables import write_contract_tables
    from gippasim.layout import LayoutConfig, simulate_layout

    layout = simulate_layout(LayoutConfig(**kw), np.random.default_rng(seed))
    work = Path(tempfile.mkdtemp(dir=root)) / "work"
    write_contract_tables(layout, work)
    return layout, dict(lo.read_parquets(work)), work


@SETTINGS
@given(kw=configs(), seed=st.integers(0, 2**32 - 1))
def test_written_tables_pass_contract(kw: dict[str, Any], seed: int) -> None:
    with tempfile.TemporaryDirectory() as root:
        _, t, work = _run(kw, seed, root)
        _, issues = lo.load_work(work)
    assert lo.errors(issues) == []
    n = kw["n_genes_chromosome"] + sum(kw["plasmids"]) + sum(kw["linear_contigs"])
    assert len(t["genes"]) == n


@SETTINGS
@given(kw=configs(), seed=st.integers(0, 2**32 - 1))
def test_structure_matches_spec(kw: dict[str, Any], seed: int) -> None:
    with tempfile.TemporaryDirectory() as root:
        _, t, _ = _run(kw, seed, root)
    g, rep, tus = t["genes"], t["replicons"], t["tus"]
    assert lo.check_gene_coordinates(g, rep) == []
    assert lo.check_gaps(g, tus) == []
    assert lo.check_tus(g, tus, rep) == []
    assert lo.check_locus_tags(g, rep) == []
    # replicons: chr1, p1.., ctg1.. in that order with the configured gene counts
    names = ["chr1", *(f"p{i + 1}" for i in range(len(kw["plasmids"])))]
    names += [f"ctg{i + 1}" for i in range(len(kw["linear_contigs"]))]
    assert rep["contig"].tolist() == names
    counts = g.groupby("contig").size().to_dict()
    sizes = [kw["n_genes_chromosome"], *kw["plasmids"], *kw["linear_contigs"]]
    assert [counts.get(c, 0) for c in names] == sizes
    # partial: exactly the first and last gene of each linear contig
    for c, gg in lo.genes_by_contig(g).items():
        part = gg["partial"].astype(bool).tolist()
        if c.startswith("ctg"):
            want = [True] * len(gg) if len(gg) <= 2 else [True, *[False] * (len(gg) - 2), True]
        else:
            want = [False] * len(gg)
        assert part == want, c
    # one plasmid region per plasmid with all its genes
    for i in range(len(kw["plasmids"])):
        mem = lo.members(t["region_genes"], f"plasmid_{i + 1}")
        assert sorted(mem) == sorted(g.loc[g["contig"] == f"p{i + 1}", "locus_tag"])


@SETTINGS
@given(kw=configs(), seed=st.integers(0, 2**32 - 1))
def test_pairs_and_layers_agree(kw: dict[str, Any], seed: int) -> None:
    with tempfile.TemporaryDirectory() as root:
        layout, t, _ = _run(kw, seed, root)
    pairs = layout.operon_pairs
    got = {frozenset(p) for p in zip(pairs["gene_a"], pairs["gene_b"], strict=True)}
    need = {frozenset(p) for p in lo.adjacent_same_strand(t["genes"], t["replicons"])}
    allowed = {frozenset(p) for p in lo.adjacent_same_strand(t["genes"], t["replicons"], True)}
    assert need <= got <= allowed
    want = lo.expected_uniop(t["genes"], t["replicons"], pairs)
    lay = t["layers"]
    for tag, prev, nxt in zip(
        lay["locus_tag"], lay["uniop_p_prev"], lay["uniop_p_next"], strict=True
    ):
        for v, e in ((prev, want[tag][0]), (nxt, want[tag][1])):
            assert pd.isna(v) if e is None else abs(float(v) - e) < 1e-12
    tmap = lo.tu_of(t["tus"])
    assert dict(zip(lay["locus_tag"], lay["tu_id"], strict=True)) == tmap
    if kw["mean_tu_size"] == 1.0:
        assert set(lay["tu_size"]) == {1}
    if kw["p_same_strand"] == 1.0:
        for c, gg in lo.genes_by_contig(t["genes"]).items():
            assert gg["strand"].nunique() == 1, c


@settings(max_examples=10, deadline=None)
@given(kw=configs(), seed=st.integers(0, 2**32 - 1))
def test_pure_function_of_config_and_seed(kw: dict[str, Any], seed: int) -> None:
    with tempfile.TemporaryDirectory() as root:
        a, ta, _ = _run(kw, seed, root)
        b, tb, _ = _run(kw, seed, root)
    for key in lo.WRITTEN_TABLES:
        pd.testing.assert_frame_equal(ta[key], tb[key])
    pd.testing.assert_frame_equal(
        a.operon_pairs.reset_index(drop=True), b.operon_pairs.reset_index(drop=True)
    )


@settings(max_examples=50, deadline=None)
@given(
    idx=st.lists(st.integers(0, 29), min_size=1, max_size=30, unique=True),
    n=st.just(30),
)
def test_oracle_circular_consecutive(idx: list[int], n: int) -> None:
    """The oracle itself: one run on the circle iff a rotation makes the set an interval."""
    s = set(idx)
    brute = any(s == {(r + k) % n for k in range(len(s))} for r in range(n))
    assert lo.circular_consecutive(idx, n) == brute


@settings(max_examples=50, deadline=None)
@given(aa=st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=40))
def test_oracle_translate_inverts_back_translation(aa: str) -> None:
    rev: dict[str, str] = {}
    for c, a in lo.CODE_11.items():
        rev.setdefault(a, c)
    nt = "".join(rev[a] for a in aa) + "TAA"
    assert lo.translate(nt) == aa + "*"
