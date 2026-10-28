"""V1 ``circos``: where in the genome is the proteome? (M9-08; GIPPA_SPEC.md §10.2 V1, §6.5).

Tests the data TSV as written (``plots/data/V01_circos.tsv``, 1-based closed), the drawn elements
per track, colours, the title, the degenerate cases, binning, determinism and purity. No ``--mpl``
baselines (M9-25). The oriC/ter and GC-skew tracks are tested only through their empty states
(lead decision 4: M9-03 support files are not on main).

Interface assumed where the spec leaves it open:

- The plot lives in ``gippa.viz.circos``; importing it registers ``V01`` (name ``circos``,
  ``needs`` = {``preflight``, ``genes_tsv``}, module ``core``, no ``layers``; ``optional`` ⊇ the
  other inputs). Called through ``registered()["V01"].fn(view, cfg)``, it returns exactly one
  ``Instance`` with ``suffix=None``.
- Every artist of a data track has ``get_gid()`` = the track's enum value (see
  ``tests/helpers/viz_tracks.py`` for what one drawn element is). An origin-crossing feature is
  two elements; any other feature (or bin) is one.
- Data TSV ``category``: genes = ``detection_status``; ``region`` = ``region_type``;
  ``peptide_density`` = ``placements``; ``island`` = ``direction``; ``alt_site`` = ``confidence``.
  ``abundance`` rows of measured genes carry ``value`` = ``abundance_log2``. ``novel_orf`` rows
  span the ORF's lowest to highest coordinate (``novel_orfs.tsv`` writes ``-``-strand ORFs with
  ``start > end``). ``alt_site`` ``start`` = ``codon_start``. Promoter/terminator rows are the BED
  record in 1-based closed form. ``tu_start`` rows are pinned for ``+``-strand TUs only.
- Peptide-density bins are ``[1 + 10,000 k, min(10,000 (k + 1), length)]``; rows of bins with
  zero placements are optional (and then have ``value`` 0).
- Omitted tracks are "noted" by a visible figure text that mentions them ("peptide", "GC skew",
  "oriC", the missing layer's tool name; case-insensitive).
"""

from __future__ import annotations

import inspect
import json
import os
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers.viz_out_tree import gene_row, make_out, write_events, write_schema_table
from helpers.viz_tracks import elements, face_hexes, tagged
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib._pylab_helpers import Gcf
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text

from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.io.tsv import write_tsv
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
VAL_OUT = Path(os.environ.get("GIPPA_ROOT", REPO / "local")) / "runs" / "val-mvp-2026-10-22" / "OUT"
COLUMNS = ["track", "contig", "start", "end", "value", "category", "binned"]
TRACKS = [
    "sector", "region", "gene_plus", "gene_minus", "abundance", "peptide_density", "tu_start",
    "promoter", "terminator", "island", "gc_skew", "alt_site", "novel_orf", "oric", "ter",
]  # fmt: skip
REGION_COLUMNS = [
    "region_id", "region_type", "contig", "start", "end", "source_tool", "score", "genes", "notes",
]  # fmt: skip
TUS_COLUMNS = [
    "tu_id", "contig", "start", "end", "strand", "genes", "upstream_start", "upstream_end",
]  # fmt: skip
TWO = [("chr1", 100_000, "circular"), ("plasmid_1", 8_000, "circular")]
BLACK, WHITE = "#000000", "#FFFFFF"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


# --- builders -------------------------------------------------------------------------------


def _preflight(out: Path, replicons: Sequence[tuple[str, int, str]]) -> None:
    reps = [
        {"id": i, "kind": "plasmid" if "plasmid" in i else "chromosome", "length": n, "topology": t}
        for i, n, t in replicons
    ]
    path = out / "preflight" / "preflight.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "replicons": reps,
        "contig_map": [{"id": r["id"], "original_id": r["id"]} for r in reps],
    }
    path.write_text(json.dumps(data), encoding="utf-8")


def gene(
    tag: str,
    contig: str,
    start: int,
    end: int,
    strand: str = "+",
    status: str = "detected",
    log2: float | None = None,
    *,
    imputed: bool = False,
    length: int = 100_000,
) -> dict[str, object]:
    return gene_row(
        locus_tag=tag,
        contig=contig,
        start=start,
        end=end,
        strand=strand,
        feature_type="cds",
        detection_status=status,
        detected=status == "detected",
        partial=status == "partial_excluded",
        crosses_origin=end > length,
        abundance_log2=log2,
        abundance_imputed=imputed,
        abundance_type="none" if log2 is None else "ibaq",
    )


def pep(contig: str, start: int, end: int | None = None) -> dict[str, object]:
    return {
        "peptide": "PEPTIDEK",
        "peptide_len": 8,
        "source": "database",
        "contig": contig,
        "start": start,
        "end": start + 23 if end is None else end,
        "strand": "+",
        "frame": 1,
        "crosses_origin": False,
        "route": "exact",
        "n_diff": 0,
        "placement_q": 0.0,
        "locus_class": "cds_inframe",
        "uniqueness": "unique",
    }


def island(iid: str, direction: str, contig: str, start: int, end: int, p: float) -> dict:
    return {
        "island_id": iid,
        "direction": direction,
        "contig": contig,
        "start": start,
        "end": end,
        "tu_first": "TU_1",
        "tu_last": "TU_2",
        "n_tus": 2,
        "n_genes": 3,
        "genes": "g1;g2;g3",
        "scan_s": 3.0,
        "block_len": 3,
        "stat": "scan_max",
        "effect": 1.0,
        "p_value": p,
        "method": "scan_max_stat",
        "null_model": "circular_block_bootstrap",
        "n_perm": 999,
    }


def variant(vid: str, contig: str, codon_start: int, call: str, confidence: str) -> dict:
    return {
        "variant_id": vid,
        "contig": contig,
        "codon_start": codon_start,
        "strand": "+",
        "gene": "g2",
        "aa_pos": 10,
        "ref_codon": "GCT",
        "ref_aa": "A",
        "alt_aa": "V",
        "alt_type": "substitution",
        "call": call,
        "confidence": confidence,
        "n_ref": 1,
        "n_alt": 3,
        "nt_changes": 1,
    }


def novel(oid: str, contig: str, strand: str, start: int, end: int, called: bool) -> dict:
    return {
        "orf_id": oid,
        "class": "intergenic",
        "contig": contig,
        "strand": strand,
        "frame": 1,
        "start_longest": start,
        "start_best_ostir": start,
        "end": end,
        "length_aa_longest": 99,
        "length_aa_best": 99,
        "n_peptides": 2,
        "called": called,
        "rank": 1,
    }


def region(rid: str, rtype: str, contig: str, start: int | None, end: int | None) -> list:
    return [rid, rtype, contig, start, end, "Bakta", None, ["g1"], ""]


def _bed(path: Path, rows: Sequence[tuple[str, int, int, str, int, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows), encoding="utf-8")


def make_run(
    tmp: Path,
    *,
    replicons: Sequence[tuple[str, int, str]] = TWO,
    genes: Sequence[Mapping[str, object]] = (),
    peptides: Sequence[Mapping[str, object]] | None = None,
    regions: Sequence[list] | None = None,
    tus: Sequence[list] | None = None,
    promoters: Sequence[tuple[str, int, int, str, int, str]] | None = None,
    terminators: Sequence[tuple[str, int, int, str, int, str]] | None = None,
    islands: Sequence[Mapping[str, object]] | None = None,
    variants: Sequence[Mapping[str, object]] | None = None,
    novel_orfs: Sequence[Mapping[str, object]] | None = None,
    flags: Sequence[str] = (),
    name: str = "OUT",
) -> Path:
    out = make_out(tmp, name)
    _preflight(out, replicons)
    write_schema_table(out / "tables" / "genes.tsv", "genes", genes)
    if peptides is not None:
        write_schema_table(out / "mapping" / "peptide_hits.tsv", "peptide_hits", peptides)
    if regions is not None:
        p = out / "annotation" / "regions.tsv"
        p.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(p, REGION_COLUMNS, regions)
    if tus is not None:
        p = out / "annotation" / "tus.tsv"
        p.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(p, TUS_COLUMNS, tus)
    if promoters is not None:
        _bed(out / "annotation" / "promoters.bed", promoters)
    if terminators is not None:
        _bed(out / "annotation" / "terminators.bed", terminators)
    if islands is not None:
        write_schema_table(out / "tables" / "islands.tsv", "islands", islands)
    if variants is not None:
        write_schema_table(out / "tables" / "variants.tsv", "variants", variants)
    if novel_orfs is not None:
        write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", novel_orfs)
    if flags:
        records = [{"event": "layer_missing", "data": {"flag": f}} for f in flags]
        write_events(out / "logs" / "events.jsonl", records)
    return out


#: the full fixture: every track that has an input on main
FULL_GENES = [
    gene("g1", "chr1", 99_801, 100_300, "+", "detected", 20.0),  # crosses the origin
    gene("g2", "chr1", 1_001, 2_000, "+", "not_detected"),
    gene("g3", "chr1", 3_001, 4_000, "-", "single_peptide_unestimated", 15.0, imputed=True),
    gene("g4", "chr1", 5_001, 6_000, "-", "detected", 22.0),
    gene("p1", "plasmid_1", 101, 700, "+", "detected", 18.0, length=8_000),
]
FULL_PEPTIDES = [pep("chr1", s) for s in (5, 9_999, 10_001, 25_000)]  # example 3
FULL_REGIONS = [
    region("prophage_1", "prophage", "chr1", 10_001, 20_000),
    region("plasmid_1", "plasmid", "plasmid_1", 1, 8_000),
    region("bgc_1", "bgc", "chr1", 30_001, 40_000),
    region("crispr_1", "crispr", "chr1", 62_388, 62_799),
    region("amr_1", "amr", "chr1", 70_001, 71_000),
]
FULL_TUS = [
    ["TU_1", "chr1", 1_001, 2_000, "+", ["g2"], 801, 1_000],
    ["TU_2", "chr1", 3_001, 6_000, "-", ["g4", "g3"], 6_001, 6_200],
]
FULL_PROMOTERS = [("chr1", 900, 950, "promoter_1", 500, "+")]
FULL_TERMINATORS = [("chr1", 2_000, 2_030, "TerminatorNet_1", 800, "+")]
FULL_ISLANDS = [
    island("ISL_1", "high", "chr1", 1_001, 6_000, 0.01),
    island("ISL_2", "low", "chr1", 50_001, 60_000, 0.2),
    island("ISL_3", "low", "chr1", 95_001, 103_000, 0.05),  # crosses the origin; p = alpha
]
FULL_VARIANTS = [
    variant("v1", "chr1", 1_501, "ALT", "high"),
    variant("v2", "chr1", 5_101, "mixed", "candidate"),
]
FULL_NOVEL = [
    novel("n1", "chr1", "+", 80_001, 80_300, True),
    novel("n2", "chr1", "-", 84_300, 84_001, True),  # '-': start_longest > end
    novel("n3", "chr1", "+", 86_001, 86_300, False),
]
FULL_TITLE = "3 of 5 genes detected (60 %); 2 islands at FWER ≤ 0.05; 2 ALT sites, 2 novel ORFs"


def full_run(tmp: Path, name: str = "OUT", **over: Any) -> Path:
    kw: dict[str, Any] = {
        "genes": FULL_GENES,
        "peptides": FULL_PEPTIDES,
        "regions": FULL_REGIONS,
        "tus": FULL_TUS,
        "promoters": FULL_PROMOTERS,
        "terminators": FULL_TERMINATORS,
        "islands": FULL_ISLANDS,
        "variants": FULL_VARIANTS,
        "novel_orfs": FULL_NOVEL,
        "name": name,
    }
    kw.update(over)
    return make_run(tmp, **kw)


def _v01() -> base.PlotSpec:
    import gippa.viz.circos  # noqa: F401  (registers V01)

    return base.registered()["V01"]


def render(out: Path, cfg: GippaConfig) -> Instance:
    result = _v01().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    assert inst.suffix is None
    return inst


def tsv(inst: Instance, tmp: Path) -> pd.DataFrame:
    """The data TSV exactly as written (strings; ``NA`` kept as text)."""
    paths = write_instance(inst, "V01", "circos", tmp / "plots_tsv", [], 100)
    (path,) = [p for p in paths if p.suffix == ".tsv"]
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def rows(data: pd.DataFrame, track: str) -> list[tuple[str, int, int, str]]:
    d = data[data["track"] == track]
    return [
        (c, int(s), int(e), cat)
        for c, s, e, cat in zip(d["contig"], d["start"], d["end"], d["category"], strict=True)
    ]


def spans(data: pd.DataFrame, track: str) -> list[tuple[str, int, int]]:
    return [(c, s, e) for c, s, e, _ in rows(data, track)]


def texts(fig: Figure) -> list[str]:
    canvas = fig.canvas if isinstance(fig.canvas, FigureCanvasAgg) else FigureCanvasAgg(fig)
    canvas.draw()
    return [t.get_text() for t in fig.findobj(Text) if t.get_text() and t.get_visible()]


def mentions(fig: Figure, word: str) -> bool:
    return any(word.lower() in t.lower() for t in texts(fig))


def viridis(x: float) -> str:
    return mcolors.to_hex(matplotlib.colormaps["viridis"](x)).upper()


# --- registration ---------------------------------------------------------------------------


def test_v01_is_registered_as_circos() -> None:
    spec = _v01()
    assert (spec.pid, spec.name, spec.module, spec.layers) == ("V01", "circos", "core", ())
    assert set(spec.needs) == {"preflight", "genes_tsv"}
    assert set(spec.optional) >= {
        "peptide_hits", "regions_tsv", "tus_tsv", "promoters_bed", "terminators_bed", "islands",
        "variants", "novel_orfs",
    }  # fmt: skip


def test_no_genes_tsv_is_skipped(tmp_path: Path, cfg: GippaConfig) -> None:
    out = make_out(tmp_path)
    _preflight(out, TWO)
    assert _v01().fn(RunView.from_outdir(out), cfg) == Skip("input_absent:genes_tsv")


def test_layer_flags_do_not_skip_the_plot(tmp_path: Path, cfg: GippaConfig) -> None:
    out = full_run(tmp_path, flags=["layer_missing:genomad", "layer_missing:antismash"])
    render(out, cfg)


# --- data TSV: shape, order, examples -------------------------------------------------------


def test_data_tsv_columns_and_enum_order(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(full_run(tmp_path), cfg), tmp_path)
    assert list(data.columns) == COLUMNS
    assert set(data["track"]) <= set(TRACKS)
    assert set(data["binned"]) == {"false"}
    fasta = {"chr1": 0, "plasmid_1": 1}
    keys = [
        (TRACKS.index(t), fasta[c], int(s), int(e), cat)
        for t, c, s, e, cat in zip(
            data["track"], data["contig"], data["start"], data["end"], data["category"], strict=True
        )
    ]
    assert keys == sorted(keys)


def test_example1_one_sector_row_per_replicon(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(make_run(tmp_path), cfg), tmp_path)
    assert spans(data, "sector") == [("chr1", 1, 100_000), ("plasmid_1", 1, 8_000)]


def test_sectors_follow_fasta_order_not_names(tmp_path: Path, cfg: GippaConfig) -> None:
    reps = [("plasmid_1", 8_000, "circular"), ("chr1", 100_000, "circular")]
    data = tsv(render(make_run(tmp_path, replicons=reps), cfg), tmp_path)
    assert spans(data, "sector") == [("plasmid_1", 1, 8_000), ("chr1", 1, 100_000)]


def test_example2_origin_crossing_gene_one_row_two_pieces(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = [gene("g1", "chr1", 99_801, 100_300, "+", "detected")]
    inst = render(make_run(tmp_path, genes=genes), cfg)
    assert rows(tsv(inst, tmp_path), "gene_plus") == [("chr1", 99_801, 100_300, "detected")]
    assert len(elements(inst.figure, "gene_plus")) == 2
    assert elements(inst.figure, "gene_minus") == []


def test_gene_rows_by_strand_with_detection_status(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(full_run(tmp_path), cfg), tmp_path)
    assert rows(data, "gene_plus") == [
        ("chr1", 1_001, 2_000, "not_detected"),
        ("chr1", 99_801, 100_300, "detected"),
        ("plasmid_1", 101, 700, "detected"),
    ]
    assert rows(data, "gene_minus") == [
        ("chr1", 3_001, 4_000, "single_peptide_unestimated"),
        ("chr1", 5_001, 6_000, "detected"),
    ]


def test_example3_peptide_density_10kb_bins_by_start(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(full_run(tmp_path), cfg), tmp_path)
    d = data[data["track"] == "peptide_density"]
    assert set(d["category"]) == {"placements"}
    got = {
        (c, int(s), int(e)): float(v)
        for c, s, e, v in zip(d.contig, d.start, d.end, d.value, strict=True)
    }
    nonzero = {k: v for k, v in got.items() if v != 0}
    assert nonzero == {
        ("chr1", 1, 10_000): 2.0,
        ("chr1", 10_001, 20_000): 1.0,
        ("chr1", 20_001, 30_000): 1.0,
    }


def test_peptide_density_last_bin_stops_at_the_sector_end(tmp_path: Path, cfg: GippaConfig) -> None:
    peps = [pep("plasmid_1", 7_990, 8_013), pep("chr1", 99_999, 100_022)]
    data = tsv(render(make_run(tmp_path, peptides=peps), cfg), tmp_path)
    d = data[(data["track"] == "peptide_density") & (data["value"].astype(float) > 0)]
    assert spans(d, "peptide_density") == [("chr1", 90_001, 100_000), ("plasmid_1", 1, 8_000)]
    every = data[data["track"] == "peptide_density"]
    ends = {"chr1": 100_000, "plasmid_1": 8_000}
    assert all(int(e) <= ends[c] for c, e in zip(every["contig"], every["end"], strict=True))


def test_peptide_density_counts_rows_of_a_multi_locus_peptide(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    peps = [pep("chr1", 101), pep("chr1", 5_001), pep("chr1", 15_001)]  # one peptide, 3 loci
    data = tsv(render(make_run(tmp_path, peptides=peps), cfg), tmp_path)
    d = data[data["track"] == "peptide_density"]
    assert d["value"].astype(float).sum() == 3.0


def test_example4_only_claimed_islands(tmp_path: Path, cfg: GippaConfig) -> None:
    isl = FULL_ISLANDS[:2]
    data = tsv(render(make_run(tmp_path, islands=isl), cfg), tmp_path)
    assert rows(data, "island") == [("chr1", 1_001, 6_000, "high")]


def test_islands_at_alpha_included_origin_crossing_keeps_end(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    inst = render(full_run(tmp_path), cfg)
    assert rows(tsv(inst, tmp_path), "island") == [
        ("chr1", 1_001, 6_000, "high"),
        ("chr1", 95_001, 103_000, "low"),
    ]
    assert len(elements(inst.figure, "island")) == 3  # ISL_3 is two pieces


def test_island_alpha_comes_from_config(tmp_path: Path) -> None:
    cfg = resolve_config(overrides=["stats.alpha.context=0.25"])
    inst = render(make_run(tmp_path, islands=FULL_ISLANDS), cfg)
    assert len(rows(tsv(inst, tmp_path), "island")) == 3
    assert "3 islands at FWER ≤ 0.25" in inst.title


def test_example5_alt_sites(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(full_run(tmp_path), cfg), tmp_path)
    got = [(c, s, cat) for c, s, _, cat in rows(data, "alt_site")]
    assert got == [("chr1", 1_501, "high"), ("chr1", 5_101, "candidate")]


def test_novel_orfs_called_only_lowest_to_highest(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(full_run(tmp_path), cfg), tmp_path)
    assert spans(data, "novel_orf") == [("chr1", 80_001, 80_300), ("chr1", 84_001, 84_300)]


def test_regions_all_five_types(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(full_run(tmp_path), cfg), tmp_path)
    assert rows(data, "region") == [
        ("chr1", 10_001, 20_000, "prophage"),
        ("chr1", 30_001, 40_000, "bgc"),
        ("chr1", 62_388, 62_799, "crispr"),
        ("chr1", 70_001, 71_000, "amr"),
        ("plasmid_1", 1, 8_000, "plasmid"),
    ]


def test_tu_promoter_terminator_rows(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(full_run(tmp_path), cfg), tmp_path)
    assert ("chr1", 1_001) in [(c, s) for c, s, _ in spans(data, "tu_start")]
    assert spans(data, "promoter") == [("chr1", 901, 950)]  # BED 900-950 -> 1-based 901-950
    assert spans(data, "terminator") == [("chr1", 2_001, 2_030)]


def test_abundance_values_of_measured_genes(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(full_run(tmp_path), cfg), tmp_path)
    d = data[(data["track"] == "abundance") & (data["value"] != "NA")]
    got = {(c, int(s)): float(v) for c, s, v in zip(d.contig, d.start, d.value, strict=True)}
    assert got == {("chr1", 5_001): 22.0, ("chr1", 99_801): 20.0, ("plasmid_1", 101): 18.0}


# --- drawn elements and colours -------------------------------------------------------------


def test_gene_colours_by_detection_status(tmp_path: Path, cfg: GippaConfig) -> None:
    fig = render(full_run(tmp_path), cfg).figure
    det, nd, un = (style.STATUS_COLOUR[k].upper() for k in
                   ("detected", "not_detected", "single_peptide_unestimated"))  # fmt: skip
    assert face_hexes(elements(fig, "gene_plus")) == sorted([det, det, det, nd])  # g1 two pieces
    assert face_hexes(elements(fig, "gene_minus")) == sorted([un, det])


def test_partial_gene_is_grey_and_hatched(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = [gene("gp", "chr1", 1_001, 2_000, "+", "partial_excluded")]
    (el,) = elements(render(make_run(tmp_path, genes=genes), cfg).figure, "gene_plus")
    assert el.face_hex == style.STATUS_COLOUR["partial_excluded"].upper()
    assert el.hatch


def test_region_colours(tmp_path: Path, cfg: GippaConfig) -> None:
    fig = render(full_run(tmp_path), cfg).figure
    got = Counter(e.face_hex for e in elements(fig, "region"))
    assert got == Counter(v.upper() for v in style.REGION.values())


def test_island_colours_and_hatches(tmp_path: Path, cfg: GippaConfig) -> None:
    els = elements(render(full_run(tmp_path), cfg).figure, "island")
    by_hatch = Counter(e.hatch for e in els)
    assert by_hatch == Counter({style.ISLAND["high"][1]: 1, style.ISLAND["low"][1]: 2})
    for e in els:
        direction = "high" if e.hatch == style.ISLAND["high"][1] else "low"
        assert style.ISLAND[direction][0].upper() in e.colours()


def test_tick_colours(tmp_path: Path, cfg: GippaConfig) -> None:
    fig = render(full_run(tmp_path), cfg).figure
    for track, colour in (("tu_start", BLACK), ("promoter", "#009E73"), ("terminator", "#D55E00")):
        els = elements(fig, track)
        assert els, track
        assert all(colour in e.colours() for e in els), track


def test_marker_styles(tmp_path: Path, cfg: GippaConfig) -> None:
    fig = render(full_run(tmp_path), cfg).figure
    alt = elements(fig, "alt_site")
    assert len(alt) == 2
    filled = [e for e in alt if not e.hollow]
    hollow = [e for e in alt if e.hollow]
    assert len(filled) == 1 and filled[0].face_hex == BLACK
    assert len(hollow) == 1 and BLACK in hollow[0].colours()
    orfs = elements(fig, "novel_orf")
    assert len(orfs) == 2
    assert all("#CC79A7" in e.colours() for e in orfs)


def test_abundance_viridis_2nd_to_98th_percentile_na_and_imputed_white(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    genes = [gene(f"m{v}", "chr1", 1 + 900 * v, 500 + 900 * v, "+", "detected", float(v))
             for v in range(101)]  # fmt: skip
    genes += [
        gene("na", "chr1", 95_001, 95_500, "+", "not_detected"),
        gene("imp", "chr1", 96_001, 96_500, "+", "detected", 200.0, imputed=True),
    ]
    inst = render(make_run(tmp_path, genes=genes), cfg)
    els = elements(inst.figure, "abundance")
    assert len(els) == len(genes)  # one cell per gene
    faces = Counter(e.face_hex for e in els)
    assert faces[WHITE] == 2
    assert faces[viridis(0.0)] >= 3  # 0, 1, 2 lie at or below the 2nd percentile
    assert faces[viridis(1.0)] >= 3  # 98, 99, 100 lie at or above the 98th
    assert faces[viridis(0.5)] >= 1  # 50 is the midpoint of [2, 98]


def test_track_patch_counts_match_features(tmp_path: Path, cfg: GippaConfig) -> None:
    fig = render(full_run(tmp_path), cfg).figure
    assert len(elements(fig, "region")) == 5
    assert len(elements(fig, "gene_plus")) == 4
    assert len(elements(fig, "gene_minus")) == 2
    assert len(elements(fig, "promoter")) == 1
    assert len(elements(fig, "terminator")) == 1


# --- title ----------------------------------------------------------------------------------


def test_title_full(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = render(full_run(tmp_path), cfg)
    assert inst.title == FULL_TITLE
    assert FULL_TITLE in texts(inst.figure)


def test_title_without_marker_tables(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = render(full_run(tmp_path, variants=None, novel_orfs=None), cfg)
    assert inst.title == "3 of 5 genes detected (60 %); 2 islands at FWER ≤ 0.05"


# --- degenerate cases -----------------------------------------------------------------------


def test_one_linear_contig_does_not_wrap(tmp_path: Path, cfg: GippaConfig) -> None:
    reps = [("chr1", 100_000, "linear")]
    genes = [gene("g1", "chr1", 99_501, 100_000, "+", "detected")]
    inst = render(make_run(tmp_path, replicons=reps, genes=genes), cfg)
    data = tsv(inst, tmp_path)
    assert spans(data, "sector") == [("chr1", 1, 100_000)]
    assert rows(data, "gene_plus") == [("chr1", 99_501, 100_000, "detected")]
    assert len(elements(inst.figure, "gene_plus")) == 1


def test_contig_without_genes_keeps_its_sector(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = [gene("g1", "chr1", 1_001, 2_000)]
    data = tsv(render(make_run(tmp_path, genes=genes), cfg), tmp_path)
    assert [c for c, _, _ in spans(data, "sector")] == ["chr1", "plasmid_1"]
    gene_rows = data[data["track"].isin(["gene_plus", "gene_minus", "abundance"])]
    assert "plasmid_1" not in set(gene_rows["contig"])


def test_no_peptide_hits_omits_density_and_notes_it(tmp_path: Path, cfg: GippaConfig) -> None:
    with_hits = render(full_run(tmp_path, name="A"), cfg)
    without = render(full_run(tmp_path, name="B", peptides=None), cfg)
    assert (tsv(without, tmp_path)["track"] == "peptide_density").sum() == 0
    assert elements(without.figure, "peptide_density") == []
    extra = set(texts(without.figure)) - set(texts(with_hits.figure))
    assert any("peptide" in t.lower() for t in extra), extra


def test_marker_tables_absent_omit_markers(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = render(full_run(tmp_path, variants=None, novel_orfs=None), cfg)
    data = tsv(inst, tmp_path)
    assert not data["track"].isin(["alt_site", "novel_orf"]).any()
    assert elements(inst.figure, "alt_site") == elements(inst.figure, "novel_orf") == []


def test_only_variants_absent_keeps_novel_orfs(tmp_path: Path, cfg: GippaConfig) -> None:
    data = tsv(render(full_run(tmp_path, variants=None), cfg), tmp_path)
    assert (data["track"] == "alt_site").sum() == 0
    assert (data["track"] == "novel_orf").sum() == 2


def test_no_support_files_omit_gc_skew_and_landmarks_noted(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    inst = render(full_run(tmp_path), cfg)
    data = tsv(inst, tmp_path)
    assert not data["track"].isin(["gc_skew", "oric", "ter"]).any()
    for track in ("gc_skew", "oric", "ter"):
        assert elements(inst.figure, track) == []
    assert mentions(inst.figure, "GC skew")
    assert mentions(inst.figure, "oriC")


def test_missing_layers_bakta_arcs_still_drawn_and_noted(tmp_path: Path, cfg: GippaConfig) -> None:
    regs = [FULL_REGIONS[1], FULL_REGIONS[3]]  # plasmid and CRISPR (Bakta)
    out = full_run(
        tmp_path, regions=regs, flags=["layer_missing:genomad", "layer_missing:antismash"]
    )
    inst = render(out, cfg)
    assert [cat for *_, cat in rows(tsv(inst, tmp_path), "region")] == ["crispr", "plasmid"]
    assert len(elements(inst.figure, "region")) == 2
    assert mentions(inst.figure, "genomad")
    assert mentions(inst.figure, "antismash")


def test_optional_annotation_files_absent_render(tmp_path: Path, cfg: GippaConfig) -> None:
    out = full_run(tmp_path, regions=None, tus=None, promoters=None, terminators=None, islands=None)
    data = tsv(render(out, cfg), tmp_path)
    absent = ["region", "tu_start", "promoter", "terminator", "island"]
    assert not data["track"].isin(absent).any()


def test_amr_gene_set_without_coordinates_renders(tmp_path: Path, cfg: GippaConfig) -> None:
    # annotation/regions.tsv writes AMR gene sets with NA start/end (gippa.annotate.regions)
    regs = [*FULL_REGIONS[:4], region("amr_1", "amr", "chr1", None, None)]
    render(full_run(tmp_path, regions=regs), cfg)


def test_small_replicon_among_large_ones_keeps_its_sector(tmp_path: Path, cfg: GippaConfig) -> None:
    reps = [("chr1", 5_000_000, "circular"), ("chr2", 5_000_000, "circular"),
            ("plasmid_1", 20_000, "circular")]  # fmt: skip
    inst = render(make_run(tmp_path, replicons=reps), cfg)
    assert spans(tsv(inst, tmp_path), "sector")[-1] == ("plasmid_1", 1, 20_000)
    assert "plasmid_1" in texts(inst.figure)
    assert any(re.search(r"\dMb|\d Mb", t) for t in texts(inst.figure))


def test_kb_tick_labels_on_a_100kb_genome(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = render(make_run(tmp_path), cfg)
    assert any(re.search(r"\d ?kb", t) for t in texts(inst.figure))


# --- binning (example 6) --------------------------------------------------------------------


def test_example6_25000_genes_are_binned(tmp_path: Path, cfg: GippaConfig) -> None:
    length = 5_000_000
    genes = [
        gene(f"g{i}", "chr1", 1 + 200 * i, 150 + 200 * i, "+", "detected", 20.0, length=length)
        for i in range(25_000)
    ]
    inst = render(make_run(tmp_path, replicons=[("chr1", length, "circular")], genes=genes), cfg)
    data = tsv(inst, tmp_path)
    plus = data[data["track"] == "gene_plus"]
    assert 0 < len(plus) <= 20_000
    assert set(plus["binned"]) == {"true"}
    assert set(plus["value"].astype(float)) == {1.0}  # every gene detected: share 1
    ab = data[data["track"] == "abundance"]
    assert 0 < len(ab) <= 20_000
    assert set(ab["binned"]) == {"true"}
    assert set(ab["value"].astype(float)) == {20.0}  # bin medians
    assert (data["track"] == "gene_minus").sum() == 0
    for track, arts in tagged(inst.figure).items():
        assert len(elements(inst.figure, track)) <= 20_000, (track, len(arts))
    assert mentions(inst.figure, "detected share")


# --- figure conventions, determinism, purity ------------------------------------------------


def test_double_width_and_minimum_font_size(tmp_path: Path, cfg: GippaConfig) -> None:
    fig = render(full_run(tmp_path), cfg).figure
    assert isinstance(fig, Figure)
    assert round(float(fig.get_size_inches()[0]), 6) == style.DOUBLE_IN
    texts(fig)
    sizes = [t.get_fontsize() for t in fig.findobj(Text) if t.get_text() and t.get_visible()]
    assert sizes and min(sizes) >= 6


def test_no_pyplot_figure_is_created(tmp_path: Path, cfg: GippaConfig) -> None:
    before = Gcf.get_num_fig_managers()
    render(full_run(tmp_path), cfg)
    assert Gcf.get_num_fig_managers() == before


def test_viz_module_is_pure() -> None:
    import gippa.viz.circos as mod

    src = inspect.getsource(mod)
    for banned in ("pyplot", "open(", "read_text", "read_csv", "to_csv", "savefig"):
        assert banned not in src, banned


def test_rerenders_are_byte_identical(tmp_path: Path, cfg: GippaConfig) -> None:
    out = full_run(tmp_path)
    written = []
    for name in ("r1", "r2"):
        paths = write_instance(
            render(out, cfg), "V01", "circos", tmp_path / name, ["svg", "png"], 100
        )
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert written[0].keys() == written[1].keys()
    for key in written[0]:
        assert written[0][key] == written[1][key], key


def test_run_view_from_inputs_gives_the_same_data(tmp_path: Path, cfg: GippaConfig) -> None:
    out = full_run(tmp_path)
    inputs = {
        "preflight": out / "preflight" / "preflight.json",
        "genes_tsv": out / "tables" / "genes.tsv",
        "islands": out / "tables" / "islands.tsv",
        "variants": out / "tables" / "variants.tsv",
        "novel_orfs": out / "tables" / "novel_orfs.tsv",
        "annotation": out / "annotation",
        "mapping": out / "mapping",
    }
    a = tsv(render(out, cfg), tmp_path / "a")
    result = _v01().fn(RunView.from_inputs(inputs), cfg)
    assert not isinstance(result, Skip)
    pd.testing.assert_frame_equal(a, tsv(result[0], tmp_path / "b"))


# --- property: peptide-density bins ---------------------------------------------------------


@settings(
    max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    starts=st.lists(st.integers(1, 100_000), min_size=1, max_size=30),
    pstarts=st.lists(st.integers(1, 8_000), max_size=10),
)
def test_density_bins_tile_from_1_and_conserve_rows(
    tmp_path_factory: pytest.TempPathFactory, starts: list[int], pstarts: list[int]
) -> None:
    tmp = tmp_path_factory.mktemp("prop")
    peps = [pep("chr1", s, s + 8) for s in starts] + [pep("plasmid_1", s, s + 8) for s in pstarts]
    data = tsv(render(make_run(tmp, peptides=peps), resolve_config()), tmp)
    d = data[data["track"] == "peptide_density"]
    length = {"chr1": 100_000, "plasmid_1": 8_000}
    for contig, n in (("chr1", len(starts)), ("plasmid_1", len(pstarts))):
        c = d[d["contig"] == contig]
        assert c["value"].astype(float).sum() == n
        for s, e, v in zip(c["start"].astype(int), c["end"].astype(int), c["value"], strict=True):
            assert (s - 1) % 10_000 == 0
            assert e == min(s + 9_999, length[contig])
            expect = sum(s <= x <= e for x in (starts if contig == "chr1" else pstarts))
            assert float(v) == expect


@settings(
    max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(start=st.integers(1, 99_999), width=st.integers(1, 5_000))
def test_gene_pieces_two_iff_crossing_origin(
    tmp_path_factory: pytest.TempPathFactory, start: int, width: int
) -> None:
    tmp = tmp_path_factory.mktemp("prop")
    end = start + width - 1
    genes = [gene("g1", "chr1", start, end, "+", "detected")]
    inst = render(make_run(tmp, genes=genes), resolve_config())
    assert rows(tsv(inst, tmp), "gene_plus") == [("chr1", start, end, "detected")]
    assert len(elements(inst.figure, "gene_plus")) == (2 if end > 100_000 else 1)


# --- example 8: the validation run ----------------------------------------------------------


@pytest.mark.slow
@pytest.mark.skipif(not VAL_OUT.is_dir(), reason="validation run OUT/ not present")
def test_example8_validation_run(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = render(VAL_OUT, cfg)
    assert inst.title == "0 of 87 genes detected (0 %); 0 islands at FWER ≤ 0.05"
    data = tsv(inst, tmp_path)
    assert not data["track"].isin(["island", "alt_site", "novel_orf"]).any()
    assert spans(data, "sector") == [("chr1", 1, 100_000), ("plasmid_1", 1, 8_000)]
