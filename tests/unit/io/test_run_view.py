"""``giap.io.run_view.RunView``: every file read of the plot framework (M9-01 behaviour 3-6,
examples 7 and 10, degenerate cases).

Interface assumed beyond the spec text:
- ``giap.io.run_view.PATHS: Mapping[str, str]`` is the one key -> ``OUT/``-relative-path table
  (behaviour 3). It extends ``PUBLISHED_FILES`` with, at least, ``preflight``, the four mapping
  tables (``peptide_hits``, ``protein_hits``, ``unmapped``, ``mapping_fdr``) and the annotation
  files (``tus``, ``regions``, ``promoters``, ``terminators``).
- ``from_inputs`` takes the step input map: file keys point at files (work files carry no
  extension), the ``PUBLISHED`` directory keys (``mapping``, ``annotation``) point at directories.
- ``NotAGiapRun`` lives in ``giap.io.run_view``; ``table`` of a key not in ``PATHS`` raises
  ``KeyError``.
"""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest
from helpers.viz_out_tree import (
    ISLAND_ROW,
    gene_row,
    make_out,
    schema_columns,
    write_events,
    write_resolved_config,
    write_schema_table,
)

from giap.io.run_view import PATHS, NotAGiapRun, RunView
from giap.runner.local import PUBLISHED, PUBLISHED_FILES, publish

ROOT = Path(__file__).resolve().parents[3]
VAL_RUN = Path("local/runs/val-mvp-2026-10-22/OUT")


# ---------------------------------------------------------------- keys


def test_keys_are_pinned_against_published_files() -> None:
    for key, rel in PUBLISHED_FILES.items():
        assert PATHS[key] == rel, key


def test_keys_cover_the_files_published_has_no_entry_for() -> None:
    assert PATHS["preflight"] == "preflight/preflight.json"
    for name in ("peptide_hits", "protein_hits", "unmapped", "mapping_fdr"):
        assert PATHS[name] == f"mapping/{name}.tsv"
    assert PATHS["tus_tsv"] == "annotation/tus.tsv"
    assert PATHS["regions_tsv"] == "annotation/regions.tsv"
    assert PATHS["promoters_bed"] == "annotation/promoters.bed"
    assert PATHS["terminators_bed"] == "annotation/terminators.bed"
    # every directory key's files live under its published directory
    dirs = set(PUBLISHED.values())
    for key in ("peptide_hits", "tus_tsv"):
        assert PATHS[key].split("/")[0] in dirs


def test_unknown_key_raises_key_error(tmp_path: Path) -> None:
    view = RunView.from_outdir(make_out(tmp_path))
    with pytest.raises(KeyError):
        view.table("no_such_table")


# ---------------------------------------------------------------- what counts as a run


def test_not_a_giap_run(tmp_path: Path) -> None:
    (tmp_path / "OUT").mkdir()
    (tmp_path / "OUT" / "notes.txt").write_text("hi\n")
    with pytest.raises(NotAGiapRun):
        RunView.from_outdir(tmp_path / "OUT")


def test_preflight_alone_is_a_run(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    (out / "preflight").mkdir(parents=True)
    (out / "preflight" / "preflight.json").write_text(json.dumps({"replicons": [{"id": "c1"}]}))
    view = RunView.from_outdir(out)
    assert view.preflight() == {"replicons": [{"id": "c1"}]}
    assert view.table("islands") is None


def test_preflight_absent_is_none(tmp_path: Path) -> None:
    assert RunView.from_outdir(make_out(tmp_path)).preflight() is None


# ---------------------------------------------------------------- tables


def test_absent_table_is_none(tmp_path: Path) -> None:
    view = RunView.from_outdir(make_out(tmp_path))
    assert view.table("islands") is None
    assert view.table("peptide_hits") is None


def test_header_only_table_has_zero_rows_and_the_schema_columns(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    write_schema_table(out / "tables" / "islands.tsv", "islands")
    write_schema_table(out / "mapping" / "peptide_hits.tsv", "peptide_hits")
    view = RunView.from_outdir(out)
    islands = view.table("islands")
    assert islands is not None and len(islands) == 0
    assert list(islands.columns) == schema_columns("islands")
    hits = view.table("peptide_hits")
    assert hits is not None and len(hits) == 0
    assert list(hits.columns) == schema_columns("peptide_hits")


def test_missing_schema_column_raises_naming_table_and_column(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    write_schema_table(out / "tables" / "islands.tsv", "islands", [ISLAND_ROW], drop=["q_value"])
    view = RunView.from_outdir(out)
    with pytest.raises(ValueError, match="islands") as err:
        view.table("islands")
    assert "q_value" in str(err.value)


def test_tables_are_read_with_schema_types(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    write_schema_table(out / "tables" / "islands.tsv", "islands", [ISLAND_ROW])
    write_schema_table(
        out / "tables" / "genes.tsv",
        "genes",
        [
            gene_row(),
            gene_row(locus_tag="g2", start=501, end=900, detected=False, cog_category=None),
        ],
    )
    view = RunView.from_outdir(out)
    isl = view.table("islands")
    assert isl is not None
    assert pd.api.types.is_integer_dtype(isl["start"])
    assert isl["start"].tolist() == [1000]  # #191 B1: 0-based on read
    assert pd.api.types.is_float_dtype(isl["effect"])
    assert isl["effect"].tolist() == [1.5]
    assert pd.isna(isl.loc[0, "q_value"])
    assert pd.isna(isl.loc[0, "flags"])
    assert isl.loc[0, "island_id"] == "ISL0001"
    genes = view.table("genes_tsv")
    assert genes is not None
    assert pd.api.types.is_bool_dtype(genes["detected"])
    assert genes["detected"].tolist() == [True, False]
    assert genes["locus_tag"].tolist() == ["g1", "g2"]
    assert pd.isna(genes.loc[1, "cog_category"])  # NA = missing


def test_support_files(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    support = out / "plots" / "support"
    support.mkdir(parents=True)
    (support / "gc_skew.tsv").write_text("contig\tstart\tend\tgc_skew\nc1\t1\t10000\t0.0125\n")
    with gzip.open(support / "landmarks.tsv.gz", "wt") as fh:
        fh.write("contig\tlandmark\tposition\nc1\toric\t1\n")
    view = RunView.from_outdir(out)
    skew = view.support("gc_skew")
    assert skew is not None and skew["gc_skew"].tolist() == [0.0125]
    marks = view.support("landmarks")
    assert marks is not None and marks["landmark"].tolist() == ["oric"]
    assert view.support("null_draws") is None


# ---------------------------------------------------------------- run flags and modules


def _example_7(out: Path) -> None:
    write_events(
        out / "logs" / "steps" / "context_layers.events.jsonl",
        [
            {
                "event": "context_layers.layer_missing",
                "data": {"flag": "layer_missing:genomad"},
            },
            {"event": "context_layers.done", "data": {"n": 3}},
        ],
    )
    write_events(
        out / "logs" / "events.jsonl",
        [
            {"event": "function.run_flags", "data": {"flags": ["layer_missing:kegg"]}},
            # operons.missed_operons writes -Infinity (backlog M6-06); json accepts it
            '{"event": "operons.missed_operons", "data": {"bin_edges": [-Infinity, 0.0]}}',
            {"event": "context_layers.layer_missing", "data": {"flag": "layer_missing:genomad"}},
        ],
    )


def test_example_7_run_flags(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    _example_7(out)
    flags = RunView.from_outdir(out).run_flags()
    assert isinstance(flags, frozenset)
    assert flags == {"layer_missing:genomad", "layer_missing:kegg"}


def test_run_flags_without_events_is_empty(tmp_path: Path) -> None:
    assert RunView.from_outdir(make_out(tmp_path)).run_flags() == frozenset()


@pytest.mark.parametrize("where", ["resolved_config.yaml", "work/resolved_config.yaml"])
def test_modules_off_from_resolved_config(tmp_path: Path, where: str) -> None:
    out = make_out(tmp_path)
    write_resolved_config(out / where, ["core", "bgc"])
    assert RunView.from_outdir(out).modules_off() == {"mobile", "function"}


def test_modules_off_without_resolved_config_is_empty(tmp_path: Path) -> None:
    assert RunView.from_outdir(make_out(tmp_path)).modules_off() == frozenset()


# ---------------------------------------------------------------- from_inputs == from_outdir


def test_from_inputs_and_from_outdir_give_equal_frames(tmp_path: Path) -> None:
    work = tmp_path / "work"
    work.mkdir()
    write_schema_table(work / "islands", "islands", [ISLAND_ROW])
    write_schema_table(work / "genes_tsv", "genes", [gene_row()])
    write_schema_table(work / "regions_stats", "regions_stats")
    write_schema_table(work / "mapping" / "peptide_hits.tsv", "peptide_hits")
    out = tmp_path / "OUT"
    publish(work, out)
    assert (out / "tables" / "islands.tsv").is_file()
    inputs = {
        "islands": work / "islands",
        "genes_tsv": work / "genes_tsv",
        "regions_stats": work / "regions_stats",
        "mapping": work / "mapping",
    }
    a, b = RunView.from_inputs(inputs), RunView.from_outdir(out)
    for key in ("islands", "genes_tsv", "regions_stats", "peptide_hits"):
        fa, fb = a.table(key), b.table(key)
        assert fa is not None and fb is not None, key
        pd.testing.assert_frame_equal(fa, fb)
    assert a.table("pairs") is None and b.table("pairs") is None


def _val_run() -> Path | None:
    for base in (ROOT, *ROOT.parents):
        candidate = base / VAL_RUN
        if candidate.is_dir():
            return candidate
    return None


def test_example_10_validation_run_islands(tmp_path: Path) -> None:
    run = _val_run()
    if run is None:
        pytest.skip(f"{VAL_RUN} is not on this machine (FarmShare only)")
    # copy the tables we read so the shared run directory is never touched
    out = make_out(tmp_path)
    shutil.copyfile(run / "tables" / "islands.tsv", out / "tables" / "islands.tsv")
    for view in (RunView.from_outdir(run), RunView.from_outdir(out)):
        islands = view.table("islands")
        assert islands is not None
        assert len(islands) == 0
        assert list(islands.columns) == schema_columns("islands")
        assert len(islands.columns) == 24
