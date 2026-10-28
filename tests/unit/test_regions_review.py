"""M2-19 lead decisions after code review (PR #91), items 3, 5 and 6
(docs/specs/M2-19-regions.md, "Lead decisions after code review").

* 3. The CRISPR cas-gene search is origin-aware on circular replicons: the gap is the circular
  distance between the array and the gene (also for origin-crossing arrays); the region extent
  stays in one frame (unwrapped, ``end > length`` per GIPPA_SPEC.md §6.5). The two review examples
  (chr1 circular, 100,000 bp): array ``[99500, 99900)`` with cas1 at ``[200, 1100)`` (gap 300) and
  the origin-crossing array ``[99800, 100300)`` with cas2 at ``[800, 1400)`` (gap 500). Linear
  replicons do not wrap.
* 5. Contract C7 (``regions``) has optional ``genes`` and ``notes`` columns: no ``extra_column``
  warning on a real ``regions.parquet``.
* 6. Tests for the review's surviving mutants: the ``region_layer`` step with an AMR gene (the
  ``amr`` column of ``genes.parquet`` → ``expert``); AMR rows in FASTA contig order with contig
  names that sort differently; ``cas`` is a gene-name *prefix* (``re.match``).

Item 4 (three Nissle plasmid regions on the circular recording) is in
``tests/contract/test_regions_contract.py``.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

import gippasim
from gippa.config import GippaConfig, resolve_config
from gippa.genetable.contract import validate_table
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene, NcFeature

regions = importlib.import_module("gippa.annotate.regions")  # M2-19

CHR1 = Replicon("chr1", 100_000, "circular", "chromosome")
LIN = Replicon("lin", 100_000, "linear", "chromosome")


def iv(start: int, end: int, contig: str = "chr1", strand: str = "+") -> Interval:
    return Interval(contig, start, end, strand)  # type: ignore[arg-type]


def gene(
    tag: str,
    start: int,
    end: int,
    contig: str = "chr1",
    *,
    name: str | None = None,
    product: str = "hypothetical protein",
    expert: tuple[str, ...] = (),
) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=iv(start, end, contig),
        aa="M",
        nt=None,
        start_codon="ATG",
        gene=name,
        product=product,
        expert=expert,
    )


def crispr(start: int, end: int, contig: str = "chr1") -> NcFeature:
    return NcFeature(id="arr", type="crispr", interval=iv(start, end, contig, "."))


def build(genes, features=(), reps=None):
    reps = reps if reps is not None else {"chr1": CHR1}
    return regions.build_regions(list(genes), list(features), reps, None, None, resolve_config())


def one_crispr(genes, array: NcFeature, reps=None):
    regs, rgs = build(genes, [array], reps)
    (r,) = [x for x in regs if x.region_type == "crispr"]
    return r, {x.locus_tag: x for x in rgs if x.region_id == r.region_id}


# ------------------------------------------------------------- 3. origin-aware cas search


def test_review_example_a_cas_gene_past_the_origin() -> None:
    """Array ``[99500, 99900)``, cas1 ``[200, 1100)``: circular gap 100 + 200 = 300 ≤ 1,000."""
    r, rows = one_crispr([gene("cas1", 200, 1100, name="cas1")], crispr(99_500, 99_900))
    assert r.genes == ("cas1",)
    assert (r.contig, r.interval.start, r.interval.end) == ("chr1", 99_500, 101_100)
    assert rows["cas1"].overlap_frac == pytest.approx(1.0)


def test_review_example_b_origin_crossing_array() -> None:
    """Array ``[99800, 100300)`` (crosses the origin), cas2 ``[800, 1400)``: gap 800 - 300 = 500."""
    r, _ = one_crispr([gene("cas2", 800, 1400, name="cas2")], crispr(99_800, 100_300))
    assert r.genes == ("cas2",)
    assert (r.interval.start, r.interval.end) == (99_800, 101_400)


def test_cas_gene_before_the_origin_with_array_after_it() -> None:
    """Array ``[100, 500)``, cas ``[99000, 99800)`` upstream across the origin: gap 200 + 100."""
    r, _ = one_crispr([gene("casU", 99_000, 99_800, name="cas3")], crispr(100, 500))
    assert r.genes == ("casU",)
    assert (r.interval.start, r.interval.end) == (99_000, 100_500)


def test_cas_genes_on_both_sides_of_the_origin_share_one_frame() -> None:
    """Upstream cas ``[98000, 99000)`` (gap 500) and downstream cas ``[200, 1100)`` (gap 300)."""
    genes = [gene("up", 98_000, 99_000, name="cas1"), gene("down", 200, 1100, name="cas2")]
    r, _ = one_crispr(genes, crispr(99_500, 99_900))
    assert set(r.genes) == {"up", "down"}
    assert (r.interval.start, r.interval.end) == (98_000, 101_100)
    assert r.interval.start < CHR1.length < r.interval.end


def test_circular_gap_boundary_1000_in_1001_out() -> None:
    """Across the origin: gap exactly 1,000 is in, 1,001 is out (array ``[99500, 99900)``)."""
    r_in, _ = one_crispr([gene("c", 900, 1500, name="cas1")], crispr(99_500, 99_900))
    assert r_in.genes == ("c",)
    r_out, _ = one_crispr([gene("c", 901, 1500, name="cas1")], crispr(99_500, 99_900))
    assert r_out.genes == ()
    assert (r_out.interval.start, r_out.interval.end) == (99_500, 99_900)


def test_linear_replicon_does_not_wrap() -> None:
    reps = {"lin": LIN}
    genes = [gene("c", 200, 1100, "lin", name="cas1")]
    r, _ = one_crispr(genes, crispr(99_500, 99_900, "lin"), reps)
    assert r.genes == ()
    assert (r.interval.start, r.interval.end) == (99_500, 99_900)


# ------------------------------------------------------------- 6. surviving mutants (unit)


@pytest.mark.parametrize("name", ["pcaS", "ycaS", "Lucas", "dcas9"])
def test_cas_is_a_gene_name_prefix(name: str) -> None:
    """``cas`` inside a gene name (not at its start) does not make a cas gene."""
    r, _ = one_crispr([gene("g", 61_000, 61_500, name=name)], crispr(60_000, 60_800))
    assert r.genes == ()
    assert (r.interval.start, r.interval.end) == (60_000, 60_800)


@pytest.mark.parametrize("name", ["cas1", "Cas9", "CASA", "casB2"])
def test_cas_prefix_names(name: str) -> None:
    r, _ = one_crispr([gene("g", 61_000, 61_500, name=name)], crispr(60_000, 60_800))
    assert r.genes == ("g",)


def test_amr_rows_follow_fasta_order_not_contig_names() -> None:
    """Replicons in FASTA order ``zeta``, ``alpha``, ``mid``: ``amr_1`` is ``zeta``."""
    reps = {
        "zeta": Replicon("zeta", 50_000, "circular", "chromosome"),
        "alpha": Replicon("alpha", 9_000, "circular", "plasmid"),
        "mid": Replicon("mid", 8_000, "circular", "plasmid"),
    }
    genes = [
        gene("m1", 100, 400, "mid", expert=("amrfinder",)),
        gene("a1", 100, 400, "alpha", expert=("amrfinder",)),
        gene("z1", 100, 400, "zeta", expert=("amrfinder",)),
    ]
    regs, _ = build(genes, reps=reps)
    amr = [(r.region_id, r.contig, r.genes) for r in regs if r.region_type == "amr"]
    assert amr == [
        ("amr_1", "zeta", ("z1",)),
        ("amr_2", "alpha", ("a1",)),
        ("amr_3", "mid", ("m1",)),
    ]


# ------------------------------------------------------------- 5. contract C7 columns


def _regions_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "region_id": ["prophage_1", "amr_1"],
            "region_type": ["prophage", "amr"],
            "contig": ["chr1", "chr1"],
            "start": pd.array([100, None], dtype="Int64"),
            "end": pd.array([5000, None], dtype="Int64"),
            "source_tool": ["geNOMAD", "AMRFinderPlus"],
            "score": [0.97, None],
            "genes": [["g1", "g2"], ["g9"]],
            "notes": ["Viruses;;", ""],
        }
    )


def test_c7_accepts_genes_and_notes_without_extra_column() -> None:
    issues = validate_table("regions", _regions_frame())
    assert [i for i in issues if "extra_column" in i.code] == []
    assert [i for i in issues if i.severity == "error"] == []


def test_c7_genes_and_notes_are_optional() -> None:
    issues = validate_table("regions", _regions_frame().drop(columns=["genes", "notes"]))
    assert [i for i in issues if i.severity == "error"] == []


def test_c7_other_extra_columns_still_warn() -> None:
    issues = validate_table("regions", _regions_frame().assign(my_extra=1))
    assert any("extra_column" in i.code for i in issues)


# ------------------------------------------------------------- 6. the step with an AMR gene

PARSE_OUT = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}


@pytest.fixture(scope="module")
def sim_core(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    out = tmp_path_factory.mktemp("regions_review_sim")
    res = gippasim.simulate(gippasim.SIM100K, gippasim.SIM100K_SEED, out / "sim")
    load_step_modules()
    work = out / "work"
    reg.run_step(
        "parse_bakta",
        GippaConfig(seed=1),
        {"bakta_json": res.bakta_json, "genome": res.genome_fasta},
        {k: work / v for k, v in PARSE_OUT.items()},
    )
    return {k: work / PARSE_OUT[k] for k in ("genes", "features", "replicons")}


def test_step_reads_the_amr_column_into_amr_regions(
    sim_core: dict[str, Path], tmp_path: Path
) -> None:
    from gippa.steps import registry as reg

    core = sim_core
    table = pq.read_table(core["genes"])
    tags = table.column("locus_tag").to_pylist()
    contigs = dict(zip(tags, table.column("contig").to_pylist(), strict=True))
    chr1 = [t for t in tags if contigs[t] == "chr1"]
    pl = [t for t in tags if contigs[t] != "chr1"]
    assert len(chr1) >= 2 and pl
    amr_tags = {chr1[0], chr1[-1], pl[0]}
    flags = [t in amr_tags for t in tags]
    col = pa.array(flags, type=pa.bool_())
    if "amr" in table.column_names:
        table = table.set_column(table.column_names.index("amr"), "amr", col)
    else:
        table = table.append_column("amr", col)
    gpath = tmp_path / "genes.parquet"
    pq.write_table(table, gpath)

    outs = {"regions": tmp_path / "regions.parquet", "region_genes": tmp_path / "rg.parquet"}
    reg.run_step("region_layer", GippaConfig(seed=1), {**core, "genes": gpath}, outs)
    rows = pq.read_table(outs["regions"]).to_pylist()
    amr = [r for r in rows if r["region_type"] == "amr"]
    assert [(r["region_id"], r["contig"]) for r in amr] == [
        ("amr_1", "chr1"),
        ("amr_2", contigs[pl[0]]),
    ]
    for r in amr:
        assert (r["start"], r["end"], r["source_tool"]) == (None, None, "AMRFinderPlus")
    cell = amr[0]["genes"]
    got = cell.split(";") if isinstance(cell, str) else list(cell)
    assert sorted(got) == sorted([chr1[0], chr1[-1]])
    rg = pq.read_table(outs["region_genes"]).to_pylist()
    assert {x["locus_tag"] for x in rg if x["region_id"].startswith("amr_")} == amr_tags


def test_step_regions_parquet_passes_c7_without_extra_column(
    sim_core: dict[str, Path], tmp_path: Path
) -> None:
    from gippa.steps import registry as reg

    outs = {"regions": tmp_path / "regions.parquet", "region_genes": tmp_path / "rg.parquet"}
    reg.run_step("region_layer", GippaConfig(seed=1), sim_core, outs)
    df = pq.read_table(outs["regions"]).to_pandas()
    assert {"genes", "notes"} <= set(df.columns) and len(df) >= 1
    issues = validate_table("regions", df)
    assert [i for i in issues if i.severity == "error"] == []
    assert [i for i in issues if "extra_column" in i.code] == []
