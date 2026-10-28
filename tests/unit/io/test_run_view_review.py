"""PR #191 review round (M9-01): the four ``RunView`` / writer blockers, as the lead decided them.

Pinned here (docs/specs/M9-01-style-and-plot-framework.md, index decision 11; GIPPA_SPEC.md §6.5):
- B1. ``RunView.table`` returns 0-based half-open coordinates for every key. A TSV whose schema
  declares ``"coordinates": "1-based-closed"`` is converted on read (``start - 1`` for each start
  column, ``end`` unchanged), like ``gippa.io.genes_tsv``, ``gippa.io.work_tables`` and
  ``gippa.io.context_tsv``; ``novel_orfs`` has two start columns (``start_longest``,
  ``start_best_ostir``, ``gippa.novel.outputs.POSITION_COLUMNS``). BED and Parquet are unchanged.
  ``write_instance`` writes the data TSV 1-based: ``start + 1`` for each column the Instance
  names in ``Instance.coords`` (a tuple of column names; empty by default).
- B3. ``from_inputs`` takes the local runner's input map: ``preflight`` may be the directory
  ``work/preflight/``, values may be lists of Paths (``InputPath``), and the published annotation
  files are keyed ``tus_tsv``, ``regions_tsv``, ``promoters_bed`` and ``terminators_bed`` so the
  work keys ``tus``, ``regions``, ``terminators`` (Parquet) never shadow them.
- B4. An empty BED file (``write_bed`` with no records) reads as a zero-row frame with the BED6
  columns.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest
from helpers.viz_out_tree import ISLAND_ROW, make_out, write_schema_table

from gippa.io.bed import BedRecord, write_bed
from gippa.io.genes_tsv import write_genes_parquet, write_genes_tsv
from gippa.io.plot_files import write_instance
from gippa.io.run_view import PATHS, RunView
from gippa.model import Interval, Replicon
from gippa.viz import style
from gippa.viz.base import Instance

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import gene_table_fixture as gt

BED_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture(scope="module")
def genes_frame() -> pd.DataFrame:
    return gt.build_world().assemble()


def _data_tsv(paths: list[Path]) -> list[list[str]]:
    tsv = next(p for p in paths if p.suffix == ".tsv")
    return [line.split("\t") for line in tsv.read_text(encoding="utf-8").splitlines()]


def _column(rows: list[list[str]], name: str) -> list[str]:
    i = rows[0].index(name)
    return [r[i] for r in rows[1:]]


# ---------------------------------------------------------------- B1 coordinates


def _coords(frame: pd.DataFrame) -> dict[str, tuple[int, int]]:
    rows = zip(frame["locus_tag"], frame["start"], frame["end"], strict=True)
    return {str(t): (int(s), int(e)) for t, s, e in rows}


def test_gene_from_parquet_and_from_genes_tsv_has_equal_coordinates(
    tmp_path: Path, genes_frame: pd.DataFrame
) -> None:
    out = make_out(tmp_path)
    (out / "genetable").mkdir()
    write_genes_parquet(out / "genetable" / "genes.parquet", genes_frame)
    write_genes_tsv(out / "tables" / "genes.tsv", genes_frame)
    view = RunView.from_outdir(out)
    pq, tsv = view.table("gene_table"), view.table("genes_tsv")
    assert pq is not None and tsv is not None
    a, b, internal = (_coords(f) for f in (pq, tsv, genes_frame))
    assert a == b
    assert a == internal  # both are the internal (0-based half-open) coordinates


def test_one_based_schema_table_start_is_converted_on_read(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    write_schema_table(out / "tables" / "islands.tsv", "islands", [ISLAND_ROW])
    isl = RunView.from_outdir(out).table("islands")
    assert isl is not None
    assert isl["start"].tolist() == [1000]  # 1001 in the file
    assert isl["end"].tolist() == [5400]  # unchanged


def test_novel_orfs_both_start_columns_are_converted_on_read(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    row = {
        "orf_id": "NORF0001",
        "class": "intergenic",
        "contig": "chr1",
        "strand": "+",
        "frame": 1,
        "start_longest": 1001,
        "start_best_ostir": 1010,
        "end": 2000,
    }
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", [row])
    orfs = RunView.from_outdir(out).table("novel_orfs")
    assert orfs is not None
    assert orfs["start_longest"].tolist() == [1000]
    assert orfs["start_best_ostir"].tolist() == [1009]
    assert orfs["end"].tolist() == [2000]


def test_bed_coordinates_are_read_unchanged(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    (out / "annotation").mkdir()
    reps = {"chr1": Replicon("chr1", 10_000, "circular")}
    rec = BedRecord(Interval("chr1", 100, 250, "+"), "p1", 500)
    write_bed(out / "annotation" / "promoters.bed", [rec], reps)
    bed = RunView.from_outdir(out).table("promoters_bed")
    assert bed is not None
    assert bed["start"].tolist() == [100]
    assert bed["end"].tolist() == [250]


def _fig_instance(data: pd.DataFrame, coords: tuple[str, ...]) -> Instance:
    fig = style.new_figure("single", 2.0)
    fig.add_subplot().plot([0, 1], [0, 1])
    return Instance(suffix=None, figure=fig, data=data, title="t", coords=coords)


def test_write_instance_writes_declared_coordinate_columns_one_based(tmp_path: Path) -> None:
    data = pd.DataFrame(
        {
            "locus_tag": ["GSIM_00005"],
            "start": pd.array([287], dtype="int64"),
            "end": pd.array([600], dtype="int64"),
            "value_start": pd.array([7], dtype="int64"),  # not declared: written as is
        }
    )
    paths = write_instance(_fig_instance(data, ("start",)), "V99", "t", tmp_path, ["svg"], 72)
    rows = _data_tsv(paths)
    assert _column(rows, "start") == ["288"]
    assert _column(rows, "end") == ["600"]
    assert _column(rows, "value_start") == ["7"]


def test_instance_coords_default_to_none_declared() -> None:
    fig = style.new_figure("single", 2.0)
    inst = Instance(suffix=None, figure=fig, data=pd.DataFrame(), title="t")
    assert tuple(inst.coords) == ()


def test_genes_tsv_round_trips_through_run_view_and_write_instance(
    tmp_path: Path, genes_frame: pd.DataFrame
) -> None:
    out = make_out(tmp_path / "run")
    write_genes_tsv(out / "tables" / "genes.tsv", genes_frame)
    genes = RunView.from_outdir(out).table("genes_tsv")
    assert genes is not None
    data = genes[["locus_tag", "start", "end"]].reset_index(drop=True)
    paths = write_instance(_fig_instance(data, ("start",)), "V99", "t", tmp_path, ["svg"], 72)
    written = _data_tsv(paths)
    published = [
        line.split("\t") for line in (out / "tables" / "genes.tsv").read_text().splitlines()
    ]
    assert _column(written, "start") == _column(published, "start")
    assert _column(written, "end") == _column(published, "end")


# ---------------------------------------------------------------- B3 from_inputs


def test_published_annotation_keys_are_renamed() -> None:
    assert PATHS["tus_tsv"] == "annotation/tus.tsv"
    assert PATHS["regions_tsv"] == "annotation/regions.tsv"
    assert PATHS["promoters_bed"] == "annotation/promoters.bed"
    assert PATHS["terminators_bed"] == "annotation/terminators.bed"
    for work_key in ("tus", "regions", "promoters", "terminators"):
        assert not PATHS.get(work_key, "").startswith("annotation/"), work_key


def test_from_inputs_preflight_directory(tmp_path: Path) -> None:
    pre = tmp_path / "work" / "preflight"
    pre.mkdir(parents=True)
    (pre / "preflight.json").write_text(json.dumps({"replicons": [{"id": "c1"}]}))
    view = RunView.from_inputs({"preflight": pre})
    assert view.preflight() == {"replicons": [{"id": "c1"}]}


def test_from_inputs_preflight_file_still_works(tmp_path: Path) -> None:
    pre = tmp_path / "preflight.json"
    pre.write_text(json.dumps({"replicons": []}))
    assert RunView.from_inputs({"preflight": pre}).preflight() == {"replicons": []}


def test_work_parquet_keys_never_shadow_published_annotation_files(tmp_path: Path) -> None:
    work = tmp_path / "work"
    ann = work / "annotation"
    ann.mkdir(parents=True)
    (ann / "tus.tsv").write_text("tu_id\tcontig\tstart\tend\tstrand\nTU0001\tchr1\t101\t400\t+\n")
    (ann / "regions.tsv").write_text("region_id\tcontig\tstart\tend\nR0001\tchr1\t11\t90\n")
    reps = {"chr1": Replicon("chr1", 10_000, "circular")}
    write_bed(ann / "promoters.bed", [BedRecord(Interval("chr1", 10, 50, "+"), "p1")], reps)
    write_bed(ann / "terminators.bed", [BedRecord(Interval("chr1", 60, 90, "-"), "t1")], reps)
    parquet = pd.DataFrame({"x": [1]})
    inputs: dict[str, Path | list[Path]] = {"annotation": ann}
    for key in ("tus", "regions", "promoters", "terminators"):
        parquet.to_parquet(work / key)  # work files carry no extension
        inputs[key] = work / key
    view = RunView.from_inputs(inputs)
    tus = view.table("tus_tsv")
    assert tus is not None and tus["tu_id"].tolist() == ["TU0001"]
    regions = view.table("regions_tsv")
    assert regions is not None and regions["region_id"].tolist() == ["R0001"]
    prom = view.table("promoters_bed")
    assert prom is not None and prom["name"].tolist() == ["p1"]
    term = view.table("terminators_bed")
    assert term is not None and term["name"].tolist() == ["t1"]


def test_from_inputs_accepts_lists_of_paths(tmp_path: Path) -> None:
    work = tmp_path / "work"
    write_schema_table(work / "islands", "islands", [ISLAND_ROW])
    write_schema_table(work / "mapping" / "peptide_hits.tsv", "peptide_hits")
    pre = work / "preflight"
    pre.mkdir()
    (pre / "preflight.json").write_text("{}")
    inputs: dict[str, Path | list[Path]] = {
        "islands": [work / "islands"],
        "mapping": [work / "mapping"],
        "preflight": [pre],
    }
    view = RunView.from_inputs(inputs)
    isl = view.table("islands")
    assert isl is not None and isl["island_id"].tolist() == ["ISL0001"]
    hits = view.table("peptide_hits")
    assert hits is not None and len(hits) == 0
    assert view.preflight() == {}


# ---------------------------------------------------------------- B4 empty BED


@pytest.mark.parametrize("key", ["promoters_bed", "terminators_bed"])
def test_empty_bed_reads_as_zero_rows_with_bed6_columns(tmp_path: Path, key: str) -> None:
    out = make_out(tmp_path)
    (out / "annotation").mkdir()
    path = out / PATHS[key]
    write_bed(path, [], {})
    assert path.is_file() and path.stat().st_size == 0
    frame = RunView.from_outdir(out).table(key)
    assert frame is not None
    assert len(frame) == 0
    assert list(frame.columns) == BED_COLUMNS
