"""Precomputed annotation layers (M2-02 examples 9-10; GIPPA_SPEC.md §5.2, §6.4, §9.2).

Files follow the ``gippa.io.tsv`` conventions (1-based closed coordinates, ``NA``, ``;`` lists);
``load_precomputed`` returns a validated frame with 0-based half-open intervals.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from gippa.io.tsv import write_tsv
from gippa.model.coords import Replicon

COLUMNS: dict[str, tuple[str, ...]] = {
    "operon_pairs": ("gene_a", "gene_b", "uniop_p", "source"),
    "promoters": ("contig", "strand", "end3p", "score"),
    "tx_rates": ("contig", "strand", "tss", "tx_rate"),
    "terminators": ("contig", "start", "end", "strand", "confidence", "source"),
    "rbs": (
        "locus_tag",
        "tir_ostir",
        "sd_dg",
        "mrna_dg",
        "spacing_dg",
        "standby_dg",
        "spacer_bp",
        "start_codon",
    ),
    "regions": (
        "region_id",
        "region_type",
        "contig",
        "start",
        "end",
        "source_tool",
        "score",
        "genes",
        "notes",
    ),
}

REPLICONS = {
    "chr1": Replicon("chr1", 1000, "circular", "chromosome"),
    "lin1": Replicon("lin1", 500, "linear", "unknown"),
}


def _load(key: str, path: Path) -> pd.DataFrame:
    from gippa.annotate.precomputed import load_precomputed

    return load_precomputed(key, path, REPLICONS)


def _write(tmp_path: Path, key: str, rows: list[list[Any]], cols: Any = None) -> Path:
    path = tmp_path / f"{key}.tsv"
    write_tsv(path, list(cols or COLUMNS[key]), rows)
    return path


def test_precomputed_columns_match_spec_table() -> None:
    from gippa.annotate.precomputed import PRECOMPUTED_COLUMNS

    for key, cols in COLUMNS.items():
        assert tuple(PRECOMPUTED_COLUMNS[key]) == cols, key


def test_precomputed_keys_are_config_keys() -> None:
    from gippa.annotate.precomputed import PRECOMPUTED_COLUMNS
    from gippa.config.model import PrecomputedConfig

    assert set(COLUMNS) <= set(PRECOMPUTED_COLUMNS) <= set(PrecomputedConfig.model_fields)


# ---------------------------------------------------------------- one test per key


def test_terminators_example_9(tmp_path: Path) -> None:
    path = _write(tmp_path, "terminators", [["chr1", 100, 160, "-", 0.82, "TerminatorNet"]])
    df = _load("terminators", path)
    assert list(df.columns) == list(COLUMNS["terminators"])
    assert len(df) == 1
    row = df.iloc[0]
    assert row["contig"] == "chr1"
    assert int(row["start"]) == 99
    assert int(row["end"]) == 160
    assert row["strand"] == "-"
    assert float(row["confidence"]) == pytest.approx(0.82)
    assert row["source"] == "TerminatorNet"


def test_operon_pairs(tmp_path: Path) -> None:
    rows = [["g1", "g2", 0.93, "uniop"], ["g2", "g3", None, "distance_model"]]
    df = _load("operon_pairs", _write(tmp_path, "operon_pairs", rows))
    assert list(df.columns) == list(COLUMNS["operon_pairs"])
    assert df["gene_a"].tolist() == ["g1", "g2"]
    assert df["gene_b"].tolist() == ["g2", "g3"]
    assert float(df["uniop_p"].iloc[0]) == pytest.approx(0.93)
    assert pd.isna(df["uniop_p"].iloc[1])
    assert df["source"].tolist() == ["uniop", "distance_model"]


def test_promoters(tmp_path: Path) -> None:
    rows = [["chr1", "+", 250, 0.71], ["lin1", "-", 40, 0.05]]
    df = _load("promoters", _write(tmp_path, "promoters", rows))
    assert list(df.columns) == list(COLUMNS["promoters"])
    assert df["contig"].tolist() == ["chr1", "lin1"]
    assert df["strand"].tolist() == ["+", "-"]
    assert pd.api.types.is_integer_dtype(df["end3p"])
    assert df["score"].tolist() == pytest.approx([0.71, 0.05])


def test_tx_rates(tmp_path: Path) -> None:
    rows = [["chr1", "+", 300, 1234.5], ["chr1", "-", 900, 12.0]]
    df = _load("tx_rates", _write(tmp_path, "tx_rates", rows))
    assert list(df.columns) == list(COLUMNS["tx_rates"])
    assert pd.api.types.is_integer_dtype(df["tss"])
    assert df["tx_rate"].tolist() == pytest.approx([1234.5, 12.0])


def test_rbs(tmp_path: Path) -> None:
    rows = [["G_0001", 5321.2, -7.1, -3.2, 0.4, 0.0, 7, "ATG"]]
    df = _load("rbs", _write(tmp_path, "rbs", rows))
    assert list(df.columns) == list(COLUMNS["rbs"])
    row = df.iloc[0]
    assert row["locus_tag"] == "G_0001"
    assert float(row["tir_ostir"]) == pytest.approx(5321.2)
    assert float(row["sd_dg"]) == pytest.approx(-7.1)
    assert int(row["spacer_bp"]) == 7
    assert row["start_codon"] == "ATG"


def test_regions(tmp_path: Path) -> None:
    rows = [
        ["R0001", "prophage", "chr1", 11, 400, "genomad", 0.97, ["G_0001", "G_0002"], None],
        ["R0002", "bgc", "lin1", 1, 500, "antismash", None, ["G_0100"], "NRPS"],
    ]
    df = _load("regions", _write(tmp_path, "regions", rows))
    assert list(df.columns) == list(COLUMNS["regions"])
    assert df["start"].astype(int).tolist() == [10, 0]
    assert df["end"].astype(int).tolist() == [400, 500]
    assert df["region_type"].tolist() == ["prophage", "bgc"]
    assert pd.isna(df["score"].iloc[1])
    assert df["notes"].iloc[1] == "NRPS"


# ---------------------------------------------------------------- validation (example 10)


def test_rbs_missing_column_example_10(tmp_path: Path) -> None:
    cols = [c for c in COLUMNS["rbs"] if c != "sd_dg"]
    path = _write(tmp_path, "rbs", [["G_0001", 5321.2, -3.2, 0.4, 0.0, 7, "ATG"]], cols)
    with pytest.raises(ValueError, match="precomputed rbs: missing column sd_dg"):
        _load("rbs", path)


@pytest.mark.parametrize("key", sorted(COLUMNS))
def test_every_key_reports_a_missing_column(key: str, tmp_path: Path) -> None:
    dropped = COLUMNS[key][-1]
    cols = list(COLUMNS[key][:-1])
    path = tmp_path / f"{key}.tsv"
    write_tsv(path, cols, [])
    with pytest.raises(ValueError, match=f"precomputed {key}: missing column {dropped}"):
        _load(key, path)


def test_unknown_key_rejected(tmp_path: Path) -> None:
    path = _write(tmp_path, "terminators", [])
    with pytest.raises((ValueError, KeyError)):
        _load("bogus_layer", path)


def test_unknown_contig_rejected(tmp_path: Path) -> None:
    path = _write(tmp_path, "terminators", [["chrX", 100, 160, "+", 0.5, "TransTermHP"]])
    with pytest.raises(ValueError, match="chrX"):
        _load("terminators", path)


def test_origin_crossing_allowed_on_circular(tmp_path: Path) -> None:
    path = _write(tmp_path, "terminators", [["chr1", 990, 1010, "+", 0.5, "TransTermHP"]])
    df = _load("terminators", path)
    assert int(df["start"].iloc[0]) == 989
    assert int(df["end"].iloc[0]) == 1010  # kept unwrapped (end > length)


def test_end_beyond_linear_contig_rejected(tmp_path: Path) -> None:
    path = _write(
        tmp_path, "regions", [["R1", "bgc", "lin1", 450, 510, "antismash", 1.0, [], None]]
    )
    with pytest.raises(ValueError):
        _load("regions", path)


def test_start_beyond_replicon_rejected(tmp_path: Path) -> None:
    path = _write(tmp_path, "terminators", [["chr1", 1001, 1020, "+", 0.5, "TransTermHP"]])
    with pytest.raises(ValueError):
        _load("terminators", path)


# ---------------------------------------------------------------- round trip and determinism


def test_terminators_round_trip_through_tsv_writer(tmp_path: Path) -> None:
    rows = [
        ["chr1", 100, 160, "-", 0.82, "TerminatorNet"],
        ["chr1", 5, 40, "+", 0.333333333, "TransTermHP"],
        ["lin1", 1, 500, "+", None, "BacTermFinder"],
    ]
    df = _load("terminators", _write(tmp_path, "terminators", rows))
    back = [
        [r.contig, int(r.start) + 1, int(r.end), r.strand, r.confidence, r.source]
        for r in df.itertuples(index=False)
    ]
    out = tmp_path / "back.tsv"
    write_tsv(out, list(COLUMNS["terminators"]), back)
    assert out.read_bytes() == (tmp_path / "terminators.tsv").read_bytes()
    assert math.isnan(float(df["confidence"].iloc[2]))


def test_load_is_deterministic(tmp_path: Path) -> None:
    rows = [["chr1", 100, 160, "-", 0.82, "TerminatorNet"], ["chr1", 5, 40, "+", 0.1, "X"]]
    path = _write(tmp_path, "terminators", rows)
    pd.testing.assert_frame_equal(_load("terminators", path), _load("terminators", path))
