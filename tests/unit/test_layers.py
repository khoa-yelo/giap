"""M2-20 ``giap.annotate.layers``: ``LAYER_COLUMNS``, ``gene_layers``, ``layer_flags``
(docs/specs/M2-20-context-layers.md examples 2-6, edge cases; GIAP_SPEC.md §9.3, §9.4, §6.3).

Inputs are the merged producers' outputs on giapsim SIM100K (``helpers.layers_chain``): every
``gene_layers`` argument is the DataFrame read with ``pd.read_parquet`` from the producer's file
(``genes``, ``tus``, ``pairs_raw``, ``promoter_windows``, ``tx_rates``, ``terminators``, ``rbs``,
``codon``, ``region_genes``, ``replicons``), ``None`` for a missing optional layer.

Pinned where the spec is silent (test author):

* One row per gene of ``genes.parquet``, in that order; columns exactly ``LAYER_COLUMNS``.
* ``region_ids`` and ``flags`` are native list cells (C6 contract, ``giap.genetable.contract``);
  ``flags`` is exactly the set {partial, crosses_origin, pseudogene, cai_reference} that applies.
* Promoter / tx / terminator columns equal the merged per-gene functions (M2-13
  ``gene_promoter_features`` with ``promotech_yes``; M2-14 ``tu_tx_features`` of the gene's TU;
  M2-15 ``terminator_conf_downstream``, checked where the next same-strand gene is a pair
  partner); rbs / codon columns equal the ``rbs`` / ``codon`` rows of the gene (NA without one).
* ``layer_flags(cfg, present)``: ``present`` is the set of ``context_layers`` input keys given;
  result sorted, unique, every item ``layer_missing:<layer>``; ``genomad`` / ``antismash`` follow
  ``cfg.modules`` (``mobile`` / ``bgc``), ``tx_rates`` follows the key (and
  ``promoter_calculator.enabled``).
* ``gc_local_dev`` (lead decision 4; GIAP_SPEC §7.8.1 / §9.3 "10-kb window GC - genome GC"):
  a float ``LAYER_COLUMNS`` column right after ``region_ids`` (its §9.3 context block);
  ``gene_layers(..., genome=<contig -> sequence>)`` fills it, ``genome=None`` gives NA. The pure
  function ``gc_local_dev(genes, replicons, genome, window=10_000)`` returns one float per gene
  row (window rules in ``helpers.layers_oracle.gc_local_dev``). "Genome GC" is read literally as
  the GC of the whole genome (all replicons): §7.8.1 / §9.3 say "genome GC", not "replicon GC".
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest
from helpers import layers_oracle as orc
from helpers.layers_chain import (
    NO_RBS_EVERY,
    PARTIAL_LOCUS,
    Chain,
    build_chain,
    call_gene_layers,
    fresh_steps,
    genome_seqs,
    precomputed_operon_pairs,
    read_layers_frames,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig

SPEC_COLUMNS = (
    "locus_tag", "tu_id", "tu_position", "tu_size", "uniop_p_prev", "uniop_p_next", "region_ids",
    "gc_local_dev", "promoter_score", "promoter_score_tu", "n_promoter_windows_tu",
    "tx_rate_max_tu", "tx_rate_sum_tu", "tss_distance_bp", "terminator_conf_downstream",
    "tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg", "spacer_bp", "start_codon",
    "cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10", "flags",
)  # fmt: skip
RBS_COLS = ("tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg", "spacer_bp", "start_codon")
CODON_COLS = ("cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10")
TX_COLS = ("tx_rate_max_tu", "tx_rate_sum_tu", "tss_distance_bp")
ALL_KEYS = {
    "genes", "features", "replicons", "pairs_raw", "tus", "promoter_windows", "tx_rates",
    "terminators", "terminators_crosscheck", "rbs", "codon", "regions", "region_genes", "genome",
}  # fmt: skip


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(tmp_path_factory.mktemp("layers_chain"))


@pytest.fixture(scope="module")
def layers(chain: Chain) -> pd.DataFrame:
    return call_gene_layers(read_layers_frames(chain), genome_seqs(chain))


def _row(df: pd.DataFrame, tag: str) -> dict[str, Any]:
    rows = df[df["locus_tag"] == tag]
    assert len(rows) == 1, f"{tag}: {len(rows)} rows"
    return rows.iloc[0].to_dict()


def _list(cell: Any) -> list[str]:
    assert cell is not None and not isinstance(cell, str), f"expected a list cell, got {cell!r}"
    return [str(x) for x in cell]


def _gene_runs(chain: Chain, strand: str) -> list[tuple[str, str, str]]:
    """(A, B, C) consecutive same-strand genes in genome order linked by pairs (A,B), (B,C)."""
    p = chain.frame("pairs_raw")
    g = chain.frame("genes").set_index("locus_tag")
    nxt = dict(zip(p["gene_a"], p["gene_b"], strict=True))
    return [
        (a, b, nxt[b]) for a, b in nxt.items()
        if b in nxt and g.loc[a, "strand"] == strand and g.loc[a, "contig"] == "chr1"
        and g.loc[a, "start"] < g.loc[b, "start"] < g.loc[nxt[b], "start"]
    ]  # fmt: skip


# ----------------------------------------------------------------------------- columns / shape


def test_layer_columns_are_the_spec_tuple() -> None:
    from giap.annotate.layers import LAYER_COLUMNS

    assert tuple(LAYER_COLUMNS) == SPEC_COLUMNS


def test_layer_columns_equal_the_c6_contract_columns() -> None:
    """M4-06 C6 consumes layers.parquet with exactly these columns."""
    from giap.annotate.layers import LAYER_COLUMNS
    from giap.genetable.contract import CONTRACT_COLUMNS

    assert tuple(LAYER_COLUMNS) == CONTRACT_COLUMNS["layers"]


def test_one_row_per_gene_in_gene_table_order(chain: Chain, layers: pd.DataFrame) -> None:
    assert tuple(layers.columns) == SPEC_COLUMNS
    assert layers["locus_tag"].tolist() == chain.frame("genes")["locus_tag"].tolist()


def test_layers_frame_passes_the_c6_contract(chain: Chain, layers: pd.DataFrame) -> None:
    from giap.genetable.contract import validate_inputs, validate_table

    errors = [i for i in validate_table("layers", layers) if i.severity == "error"]
    assert not errors, errors
    empty = pd.DataFrame()
    tables = SimpleNamespace(
        replicons=chain.frame("replicons"), genes=chain.frame("genes"), tus=chain.frame("tus"),
        layers=layers, regions=chain.frame("regions"), region_genes=chain.frame("region_genes"),
        peptides=None, peptide_hits=empty, protein_hits=empty,
    )  # fmt: skip
    bad = [i for i in validate_inputs(tables) if i.code.startswith("contract.layers")]
    assert not bad, bad


# ----------------------------------------------------------------------------- TU context


def test_tu_columns_follow_tus_parquet(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.tu_positions(chain.frame("tus"))
    for r in layers.itertuples(index=False):
        assert (r.tu_id, int(r.tu_position), int(r.tu_size)) == want[r.locus_tag], r.locus_tag


def test_every_gene_has_a_tu(layers: pd.DataFrame) -> None:
    """Edge case: genes without a TU cannot exist (singletons are TUs)."""
    assert layers["tu_id"].notna().all()
    assert (layers["tu_size"] >= 1).all() and (layers["tu_position"] >= 1).all()


def test_uniop_p_prev_next_in_transcription_order(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.prev_next(chain.frame("genes"), chain.frame("pairs_raw"))
    for r in layers.itertuples(index=False):
        wp, wn = want[r.locus_tag]
        assert orc.same(r.uniop_p_prev, wp), (r.locus_tag, r.uniop_p_prev, wp)
        assert orc.same(r.uniop_p_next, wn), (r.locus_tag, r.uniop_p_next, wn)


def test_example_2_plus_tu_middle_gene(chain: Chain, tmp_path: Path) -> None:
    """TU (A, B, C) on +, p(A,B)=0.9, p(B,C)=0.7 -> B: position 2, size 3, prev 0.9, next 0.7."""
    runs = _gene_runs(chain, "+")
    assert runs, "SIM100K has three linked + genes"
    a, b, c = runs[0]
    p = chain.frame("pairs_raw")
    over = {(a, b): 0.9, (b, c): 0.7}
    over |= {(x, a): 0.1 for x in p.loc[p["gene_b"] == a, "gene_a"]}
    over |= {(c, y): 0.1 for y in p.loc[p["gene_a"] == c, "gene_b"]}
    pairs, tus = precomputed_operon_pairs(chain, tmp_path, over)
    tu = [list(g) for g in pd.read_parquet(tus)["genes"] if a in list(g)]
    assert tu == [[a, b, c]], "the precomputed pairs make TU (A, B, C)"
    df = call_gene_layers(read_layers_frames(chain, pairs_raw=pairs, tus=tus))
    rb = _row(df, b)
    assert (int(rb["tu_position"]), int(rb["tu_size"])) == (2, 3)
    assert rb["uniop_p_prev"] == pytest.approx(0.9) and rb["uniop_p_next"] == pytest.approx(0.7)
    ra = _row(df, a)
    assert int(ra["tu_position"]) == 1 and ra["uniop_p_next"] == pytest.approx(0.9)
    if (p["gene_b"] == a).any():  # lead Q2: the first gene keeps p to its upstream neighbour
        assert ra["uniop_p_prev"] == pytest.approx(0.1)


def test_example_3_minus_tu_last_gene(chain: Chain, tmp_path: Path) -> None:
    """- TU (C, B, A) (genome order A, B, C), the gene before A in genome order is on + ->
    A: position 3, uniop_p_prev = p(A,B), uniop_p_next = NA."""
    g = chain.frame("genes")
    order = g["locus_tag"].tolist()
    strand = dict(zip(g["locus_tag"], g["strand"], strict=True))
    runs = [r for r in _gene_runs(chain, "-")
            if order.index(r[0]) > 0 and strand[order[order.index(r[0]) - 1]] == "+"]  # fmt: skip
    assert runs, "SIM100K has a - run whose genome predecessor is on +"
    a, b, c = runs[0]
    p = chain.frame("pairs_raw")
    over = {(a, b): 0.8, (b, c): 0.85}
    over |= {(c, y): 0.1 for y in p.loc[p["gene_a"] == c, "gene_b"]}
    pairs, tus = precomputed_operon_pairs(chain, tmp_path, over)
    tu = [list(x) for x in pd.read_parquet(tus)["genes"] if a in list(x)]
    assert tu == [[c, b, a]], "transcription order of a - TU"
    df = call_gene_layers(read_layers_frames(chain, pairs_raw=pairs, tus=tus))
    ra = _row(df, a)
    assert (int(ra["tu_position"]), int(ra["tu_size"])) == (3, 3)
    assert ra["uniop_p_prev"] == pytest.approx(0.8)
    assert orc.isna(ra["uniop_p_next"])
    rc = _row(df, c)
    assert int(rc["tu_position"]) == 1 and rc["uniop_p_next"] == pytest.approx(0.85)
    if (p["gene_a"] == c).any():  # upstream same-strand neighbour below the threshold (Q2)
        assert rc["uniop_p_prev"] == pytest.approx(0.1)


@settings(
    max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(data=st.data())
def test_property_tu_positions_and_shared_pair_probabilities(
    chain: Chain, tmp_path_factory: pytest.TempPathFactory, data: st.DataObject
) -> None:
    """For any pair probabilities (wrap pairs unlinked): positions are 1..size in each TU;
    consecutive genes of a TU share their pair p (next of k == prev of k+1, >= threshold); a
    TU's first gene has prev NA or < threshold and its last gene next NA or < threshold."""
    p = chain.frame("pairs_raw")
    probs = data.draw(st.lists(st.sampled_from([0.05, 0.3, 0.49, 0.5, 0.51, 0.8, 0.99]),
                               min_size=len(p), max_size=len(p)))  # fmt: skip
    start = dict(zip(chain.frame("genes")["locus_tag"], chain.frame("genes")["start"], strict=True))
    # wrap pairs (last, first) stay unlinked: a fully linked ring breaks at the origin (M2-12)
    over = {(a, b): (x if start[a] < start[b] else 0.05)
            for a, b, x in zip(p["gene_a"], p["gene_b"], probs, strict=True)}  # fmt: skip
    pairs, tus = precomputed_operon_pairs(chain, tmp_path_factory.mktemp("prop"), over)
    df = call_gene_layers(read_layers_frames(chain, pairs_raw=pairs, tus=tus)).set_index(
        "locus_tag"
    )
    thr = GiapConfig().annotation.uniop.pair_threshold
    for gs in pd.read_parquet(tus)["genes"]:
        gs = list(gs)
        assert [int(df.loc[x, "tu_position"]) for x in gs] == list(range(1, len(gs) + 1))
        assert {int(df.loc[x, "tu_size"]) for x in gs} == {len(gs)}
        for x, y in itertools.pairwise(gs):
            assert df.loc[x, "uniop_p_next"] == df.loc[y, "uniop_p_prev"] >= thr
        first, last = df.loc[gs[0], "uniop_p_prev"], df.loc[gs[-1], "uniop_p_next"]
        assert orc.isna(first) or first < thr
        assert orc.isna(last) or last < thr


# ----------------------------------------------------------------------------- gc_local_dev
# Lead decision 4: GIAP_SPEC §7.8.1 / §9.3 "gc_local_dev (10-kb window GC - genome GC)".
#
# Small genome (hand-computed, window = 10 so the half-window is 5):
#   c1 circular, 20 bp: GGGGG AAAAAAAAAA CCCCC   (0-4 G, 5-14 A, 15-19 C)  ACGT 20, GC 10
#   l1 linear,   12 bp: ATATAT GCGC NN           (N at 10, 11)              ACGT 10, GC 4
#   genome GC = (10 + 4) / (20 + 10) = 14/30 (whole genome, N excluded; with N: 14/32;
#   per replicon it would be 10/20 on c1 and 4/10 on l1).
GC_SEQS = {"c1": "GGGGGAAAAAAAAAACCCCC", "l1": "ATATATGCGCNN"}
GC_GENOME = 14 / 30
GC_CASES = [
    # (locus_tag, contig, start, end, strand, window GC by hand)
    # m = 2, window [-3, 7) wraps: 17 18 19 | 0..6 -> CCC GGGGG AA -> 8/10
    ("g1", "c1", 0, 4, "+", 8 / 10),
    # m = 10, window [5, 15): all A -> 0
    ("g2", "c1", 8, 12, "-", 0.0),
    # origin-crossing gene stored unwrapped (end 23 > 20): m = 20 % 20 = 0, window [-5, 5):
    # 15..19 C + 0..4 G -> 10/10
    ("g3", "c1", 18, 23, "+", 1.0),
    # linear: m = 1, window [-4, 6) truncated to [0, 6) = ATATAT -> 0/6 (wrapping would add GCNN)
    ("g4", "l1", 0, 2, "+", 0.0),
    # linear: m = 10, window [5, 15) truncated to [5, 12) = T GCGC NN -> 4/5 (N excluded; 4/7 if
    # N counted)
    ("g5", "l1", 9, 12, "-", 4 / 5),
]


def _gc_frames(
    cases: list[tuple[str, str, int, int, str, float]] | None = None,
    seqs: dict[str, str] | None = None,
    topology: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cases = GC_CASES if cases is None else cases
    seqs = GC_SEQS if seqs is None else seqs
    topology = topology or {"c1": "circular", "l1": "linear"}
    genes = pd.DataFrame(
        [(tag, c, s, e, st_) for tag, c, s, e, st_, _ in cases],
        columns=["locus_tag", "contig", "start", "end", "strand"],
    )
    reps = pd.DataFrame(
        [(c, len(s), topology[c], "unknown") for c, s in seqs.items()],
        columns=["contig", "length", "topology", "kind"],
    )
    return genes, reps


def _gc(genes: pd.DataFrame, reps: pd.DataFrame, seqs: dict[str, str], **kw: int) -> list[float]:
    from giap.annotate.layers import gc_local_dev

    got = gc_local_dev(genes, reps, seqs, **kw)
    assert len(got) == len(genes)
    return [float(x) for x in list(got)]


def test_layer_columns_include_gc_local_dev() -> None:
    from giap.annotate.layers import LAYER_COLUMNS

    assert "gc_local_dev" in LAYER_COLUMNS
    assert LAYER_COLUMNS.index("gc_local_dev") == LAYER_COLUMNS.index("region_ids") + 1


@pytest.mark.parametrize("case", GC_CASES, ids=[c[0] for c in GC_CASES])
def test_gc_local_dev_hand_computed(case: tuple[str, str, int, int, str, float]) -> None:
    """Window centred on the gene midpoint; circular wrap (also for an origin-crossing gene);
    truncation at linear ends; N excluded from numerator and denominator; whole-genome GC."""
    genes, reps = _gc_frames([case])
    (got,) = _gc(genes, reps, GC_SEQS, window=10)
    assert got == pytest.approx(case[5] - GC_GENOME, abs=1e-12)


def test_gc_local_dev_all_rows_match_the_oracle() -> None:
    genes, reps = _gc_frames()
    got = _gc(genes, reps, GC_SEQS, window=10)
    want = orc.gc_local_dev(genes, reps, GC_SEQS, window=10)
    assert got == pytest.approx(want, abs=1e-12)
    assert got == pytest.approx([c[5] - GC_GENOME for c in GC_CASES], abs=1e-12)


def test_gc_local_dev_genome_gc_is_whole_genome_not_per_replicon() -> None:
    """Gene on c1 whose window is all A: dev = 0 - 14/30, not 0 - 10/20 (c1's own GC)."""
    genes, reps = _gc_frames([GC_CASES[1]])
    (got,) = _gc(genes, reps, GC_SEQS, window=10)
    assert got == pytest.approx(-14 / 30, abs=1e-12)
    assert got != pytest.approx(-10 / 20, abs=1e-6)


def test_gc_local_dev_n_excluded_from_genome_gc() -> None:
    """Adding an all-N replicon changes nothing (N is in neither numerator nor denominator)."""
    seqs = {**GC_SEQS, "n1": "N" * 40}
    genes, reps = _gc_frames(seqs=seqs, topology={"c1": "circular", "l1": "linear",
                                                  "n1": "linear"})  # fmt: skip
    base_genes, base_reps = _gc_frames()
    assert _gc(genes, reps, seqs, window=10) == pytest.approx(
        _gc(base_genes, base_reps, GC_SEQS, window=10), abs=1e-12
    )


def test_gc_local_dev_window_without_acgt_is_na() -> None:
    seqs = {"c1": GC_SEQS["c1"], "n1": "N" * 30}
    cases = [("g1", "n1", 10, 14, "+", float("nan")), ("g2", "c1", 8, 12, "+", 0.0)]
    genes, reps = _gc_frames(cases, seqs, {"c1": "circular", "n1": "linear"})
    got = _gc(genes, reps, seqs, window=10)
    assert math.isnan(got[0])
    assert got[1] == pytest.approx(0.0 - 10 / 20, abs=1e-12)


def test_gc_local_dev_circular_replicon_shorter_than_window_counts_each_base_once() -> None:
    """c2 = GGGAAA (circular, 6 bp) with a 10-bp window: the whole replicon (3/6), not a
    10-position wrap with repeats (4/10). Single replicon, so dev = 0."""
    seqs = {"c2": "GGGAAA"}
    genes, reps = _gc_frames([("g1", "c2", 0, 2, "+", 0.5)], seqs, {"c2": "circular"})
    (got,) = _gc(genes, reps, seqs, window=10)
    assert got == pytest.approx(0.0, abs=1e-12)


def test_gc_local_dev_default_window_is_10_kb() -> None:
    """30-kb circular replicon, G on [0, 15 000), A after; gene midpoint 12 000. Default window
    [7 000, 17 000) -> 8 000 / 10 000 = 0.8; genome GC 0.5 -> dev 0.3 (a 9-kb window gives
    7 500 / 9 000, a 20-kb one 0.75)."""
    seqs = {"c3": "G" * 15_000 + "A" * 15_000}
    genes, reps = _gc_frames([("g1", "c3", 11_900, 12_100, "+", 0.8)], seqs, {"c3": "circular"})
    (got,) = _gc(genes, reps, seqs)
    assert got == pytest.approx(0.3, abs=1e-12)


def test_gc_local_dev_case_insensitive() -> None:
    """Soft-masked (lower-case) bases count like upper-case ones."""
    genes, reps = _gc_frames()
    lower = {c: s.lower() for c, s in GC_SEQS.items()}
    assert _gc(genes, reps, lower, window=10) == pytest.approx(
        _gc(genes, reps, GC_SEQS, window=10), abs=1e-12
    )


@st.composite
def _gc_genomes(draw: st.DrawFn) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str], int]:
    n = draw(st.integers(1, 3))
    seqs, topo, cases = {}, {}, []
    for i in range(n):
        c = f"r{i}"
        seqs[c] = draw(st.text(alphabet="ACGTNacgtnRY", min_size=1, max_size=60))
        topo[c] = draw(st.sampled_from(["circular", "linear"]))
        length = len(seqs[c])
        for k in range(draw(st.integers(0, 4))):
            s = draw(st.integers(0, length - 1))
            hi = s + length if topo[c] == "circular" else length
            e = draw(st.integers(s + 1, hi))
            cases.append((f"{c}_{k}", c, s, e, draw(st.sampled_from("+-")), 0.0))
    if not cases:
        cases.append(("r0_0", "r0", 0, 1, "+", 0.0))
    genes, reps = _gc_frames(cases, seqs, topo)
    return genes, reps, seqs, draw(st.integers(1, 80))


@settings(max_examples=150, deadline=None)
@given(case=_gc_genomes())
def test_property_gc_local_dev_equals_the_brute_force_oracle(
    case: tuple[pd.DataFrame, pd.DataFrame, dict[str, str], int],
) -> None:
    """Any genome (mixed case, N / IUPAC letters), any topology, any gene, any window: the
    implementation equals the brute-force window oracle; values lie in [-1, 1] or are NA."""
    genes, reps, seqs, window = case
    got = _gc(genes, reps, seqs, window=window)
    want = orc.gc_local_dev(genes, reps, seqs, window=window)
    for g, w in zip(got, want, strict=True):
        assert orc.same(g, w) or math.isclose(g, w, abs_tol=1e-12), (g, w)
        assert math.isnan(g) or -1.0 <= g <= 1.0


@settings(max_examples=60, deadline=None)
@given(case=_gc_genomes())
def test_property_gc_local_dev_strand_invariant(
    case: tuple[pd.DataFrame, pd.DataFrame, dict[str, str], int],
) -> None:
    """Strand plays no role (the window is centred on the midpoint, not on the start)."""
    genes, reps, seqs, window = case
    base = _gc(genes, reps, seqs, window=window)
    flipped = genes.assign(strand=genes["strand"].map({"+": "-", "-": "+"}))
    assert _gc(flipped, reps, seqs, window=window) == pytest.approx(base, abs=1e-12, nan_ok=True)


def test_gene_layers_gc_local_dev_on_sim100k(chain: Chain, layers: pd.DataFrame) -> None:
    """``gene_layers(..., genome=...)`` fills the column with the per-gene values (float)."""
    want = orc.gc_local_dev(chain.frame("genes"), chain.frame("replicons"), genome_seqs(chain))
    assert pd.api.types.is_float_dtype(layers["gc_local_dev"])
    got = layers["gc_local_dev"].tolist()
    for tag, g, w in zip(layers["locus_tag"], got, want, strict=True):
        assert orc.same(float(g), w) or math.isclose(g, w, abs_tol=1e-12), (tag, g, w)
    assert layers["gc_local_dev"].notna().all()
    assert layers["gc_local_dev"].std() > 0


def test_gene_layers_without_genome_gives_na(chain: Chain) -> None:
    df = call_gene_layers(read_layers_frames(chain))
    assert pd.api.types.is_float_dtype(df["gc_local_dev"])
    assert df["gc_local_dev"].isna().all()


# ----------------------------------------------------------------------------- regions


def test_region_ids_sorted_from_region_genes(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.region_ids(chain.frame("region_genes"))
    for r in layers.itertuples(index=False):
        assert _list(r.region_ids) == want.get(r.locus_tag, []), r.locus_tag


def test_example_4_gene_in_prophage_and_plasmid(chain: Chain, layers: pd.DataFrame) -> None:
    rg = chain.frame("region_genes")
    both = set(rg.loc[rg["region_id"] == "prophage_1", "locus_tag"]) & set(
        rg.loc[rg["region_id"] == "plasmid_1", "locus_tag"]
    )
    assert both, "the chain plants a provirus on plasmid_1"
    for tag in sorted(both):
        assert _list(_row(layers, tag)["region_ids"]) == ["plasmid_1", "prophage_1"]


def test_example_5_core_only_has_no_prophage_or_bgc_ids(chain: Chain, tmp_path: Path) -> None:
    """modules=[core]: region_layer runs without geNOMAD / antiSMASH (annotate.nf), so no
    prophage / bgc ids; the run flags name both missing tools."""
    from giap.annotate.layers import layer_flags

    outs = {"regions": tmp_path / "regions.parquet", "region_genes": tmp_path / "rg.parquet"}
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        reg.run_step(
            "region_layer", GiapConfig(seed=1, modules=["core"]),
            {k: chain.inputs[k] for k in ("genes", "features", "replicons")}, outs,
        )  # fmt: skip
    df = call_gene_layers(read_layers_frames(chain, region_genes=outs["region_genes"]))
    ids = {x for cell in df["region_ids"] for x in _list(cell)}
    assert not {x for x in ids if x.startswith(("prophage_", "bgc_"))}
    assert "plasmid_1" in ids  # the Bakta replicon table still gives the plasmid region
    flags = layer_flags(GiapConfig(seed=1, modules=["core"]), ALL_KEYS)
    assert {"layer_missing:genomad", "layer_missing:antismash"} <= set(flags)


# ----------------------------------------------------------------------- tx, promoters, terminators


def test_example_6_no_tx_rates_gives_na_and_flag(chain: Chain) -> None:
    from giap.annotate.layers import layer_flags

    df = call_gene_layers(read_layers_frames(chain, tx_rates=None))
    for c in TX_COLS:
        assert df[c].isna().all(), c
    cfg = GiapConfig.model_validate({"annotation": {"promoter_calculator": {"enabled": False}}})
    assert "layer_missing:tx_rates" in layer_flags(cfg, ALL_KEYS - {"tx_rates"})


def test_tx_columns_are_the_tu_features(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.tx_features(read_layers_frames(chain))
    assert (want["tx_rate_max"] > 0).any(), "some TUs have Promoter Calculator hits"
    for r in layers.itertuples(index=False):
        w = want.loc[r.tu_id]
        assert orc.same(r.tx_rate_max_tu, float(w["tx_rate_max"])), r.locus_tag
        assert orc.same(r.tx_rate_sum_tu, float(w["tx_rate_sum"])), r.locus_tag
        wd = w["tss_distance_bp"]
        assert orc.same(r.tss_distance_bp, None if pd.isna(wd) else int(wd)), r.locus_tag


def test_promoter_columns_are_the_gene_promoter_features(
    chain: Chain, layers: pd.DataFrame
) -> None:
    yes = GiapConfig().annotation.promoters.promotech_yes
    want = orc.promoter_features(read_layers_frames(chain), yes)
    for r in layers.itertuples(index=False):
        w = want.loc[r.locus_tag]
        for c in ("promoter_score", "promoter_score_tu", "n_promoter_windows_tu"):
            assert orc.same(getattr(r, c), w[c]), (r.locus_tag, c)


def test_no_promoter_layer_gives_na(chain: Chain) -> None:
    df = call_gene_layers(read_layers_frames(chain, promoter_windows=None))
    for c in ("promoter_score", "promoter_score_tu", "n_promoter_windows_tu"):
        assert df[c].isna().all(), c


def test_terminator_conf_downstream(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.terminator_downstream_with_partner(read_layers_frames(chain))
    assert any(not math.isnan(v) for v in want.values()), "recorded TerminatorNet calls hit"
    got = dict(zip(layers["locus_tag"], layers["terminator_conf_downstream"], strict=True))
    for tag, v in want.items():
        assert orc.same(got[tag], v), tag
    conf = set(chain.frame("terminators")["confidence"].astype(float))
    assert all(orc.isna(v) or float(v) in conf for v in got.values())


def test_no_terminator_layer_gives_na(chain: Chain) -> None:
    df = call_gene_layers(read_layers_frames(chain, terminators=None))
    assert df["terminator_conf_downstream"].isna().all()


# ----------------------------------------------------------------------------- rbs / codon / flags


def test_rbs_columns_and_missing_rbs_rows_are_na(chain: Chain, layers: pd.DataFrame) -> None:
    rbs = chain.frame("rbs").set_index("locus_tag")
    n_missing = 0
    for r in layers.itertuples(index=False):
        for c in RBS_COLS:
            w = rbs.loc[r.locus_tag, c] if r.locus_tag in rbs.index else None
            assert orc.same(getattr(r, c), w), (r.locus_tag, c)
        n_missing += r.locus_tag not in rbs.index
    assert n_missing >= len(layers) // NO_RBS_EVERY  # edge case exercised


def test_codon_columns_and_gene_missing_from_codon_is_na(chain: Chain) -> None:
    codon = chain.frame("codon")
    dropped = set(codon["locus_tag"].iloc[::5])
    frames = read_layers_frames(chain)
    frames["codon"] = codon[~codon["locus_tag"].isin(dropped)].reset_index(drop=True)
    df = call_gene_layers(frames)
    full = codon.set_index("locus_tag")
    for r in df.itertuples(index=False):
        for c in CODON_COLS:
            w = None if r.locus_tag in dropped else full.loc[r.locus_tag, c]
            assert orc.same(getattr(r, c), w), (r.locus_tag, c)


def test_no_rbs_or_codon_layer_gives_na(chain: Chain) -> None:
    df = call_gene_layers(read_layers_frames(chain, rbs=None, codon=None))
    for c in (*RBS_COLS, *CODON_COLS):
        assert df[c].isna().all(), c


def test_gene_flags(chain: Chain, layers: pd.DataFrame) -> None:
    want = orc.gene_flags(chain.frame("genes"), chain.frame("codon"))
    assert any("partial" in f for f in want.values())
    assert any("pseudogene" in f for f in want.values())
    assert any("cai_reference" in f for f in want.values())
    for r in layers.itertuples(index=False):
        assert set(_list(r.flags)) == want[r.locus_tag], r.locus_tag


def test_partial_gene_keeps_its_layer_values(chain: Chain, layers: pd.DataFrame) -> None:
    """Edge case: partial genes keep layer values (statistics exclude them later, §6.5)."""
    r = _row(layers, PARTIAL_LOCUS)
    assert "partial" in _list(r["flags"])
    assert not orc.isna(r["promoter_score"])
    assert not orc.isna(r["tu_id"]) and not orc.isna(r["at_2_10"])


# ----------------------------------------------------------------------------- layer_flags


def test_layer_flags_all_present_all_modules_is_empty() -> None:
    from giap.annotate.layers import layer_flags

    assert tuple(layer_flags(GiapConfig(), ALL_KEYS)) == ()


@pytest.mark.parametrize(
    ("modules", "present", "want"),
    [
        (["core"], ALL_KEYS, {"layer_missing:genomad", "layer_missing:antismash"}),
        (["core", "mobile"], ALL_KEYS, {"layer_missing:antismash"}),
        (["core", "bgc"], ALL_KEYS, {"layer_missing:genomad"}),
        (["core", "mobile", "bgc"], ALL_KEYS - {"tx_rates"}, {"layer_missing:tx_rates"}),
    ],
)
def test_layer_flags(modules: list[str], present: set[str], want: set[str]) -> None:
    from giap.annotate.layers import layer_flags

    got = tuple(layer_flags(GiapConfig(modules=modules), set(present)))
    assert want <= set(got)
    assert list(got) == sorted(set(got))
    assert all(f.startswith("layer_missing:") for f in got)
    if modules == ["core", "mobile"]:
        assert "layer_missing:genomad" not in got
    if modules == ["core", "bgc"]:
        assert "layer_missing:antismash" not in got
