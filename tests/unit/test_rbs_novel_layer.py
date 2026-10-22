"""``annotation.precomputed.rbs_novel`` (issue M7-08 lead decision 1, ADR-0024) and
``parse_ostir_starts`` (M7-08 public interface; GIAP_SPEC.md §6.4, §7.1.1, §7.6.2).

* config: the key exists in the model and in defaults.yaml, null by default, takes a path, and is
  in the generated config JSON schema;
* ``load_precomputed('rbs_novel', ...)``: columns ``contig, strand, start, tir_ostir, sd_dg,
  mrna_dg, spacing_dg, standby_dg, spacer_bp, start_codon``; ``start`` is a single 1-based
  position in the file (the start codon's first base in transcription direction) and 0-based in
  the frame, validated like every ``pos`` column;
* ``parse_ostir_starts(path, requested)``: every requested start of an id (1-based
  ``start_position``), unrequested positions and ids dropped, the first row of a repeated
  ``(name, start_position)`` kept, ids without a row absent.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import yaml

from giap.io.tsv import write_tsv
from giap.model.coords import Replicon

ROOT = Path(__file__).resolve().parents[2]
COLS = ("contig", "strand", "start", "tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg",
        "spacer_bp", "start_codon")  # fmt: skip
REPS = {
    "chr1": Replicon("chr1", 1000, "circular", "chromosome"),
    "lin1": Replicon("lin1", 500, "linear", "unknown"),
}
OSTIR_HEADER = ("name,start_codon,start_position,expression,RBS_distance_bp,dG_total,"
                "dG_rRNA:mRNA,dG_mRNA,dG_spacing,dG_standby,dG_start_codon")  # fmt: skip


# ---------------------------------------------------------------- config


def test_model_has_rbs_novel_null_by_default() -> None:
    from giap.config import GiapConfig

    assert GiapConfig().annotation.precomputed.rbs_novel is None  # type: ignore[attr-defined]


def test_defaults_yaml_lists_rbs_novel_null() -> None:
    data = yaml.safe_load((ROOT / "src/giap/config/defaults.yaml").read_text(encoding="utf-8"))
    pre = data["annotation"]["precomputed"]
    assert "rbs_novel" in pre and pre["rbs_novel"] is None


def test_rbs_novel_takes_a_path(tmp_path: Path) -> None:
    from giap.config import resolve_config

    p = tmp_path / "rbs_novel.tsv"
    cfg = resolve_config(None, [f"annotation.precomputed.rbs_novel={p}"])
    assert Path(cfg.annotation.precomputed.rbs_novel) == p  # type: ignore[attr-defined]


def test_config_schema_lists_rbs_novel() -> None:
    schema = json.loads((ROOT / "src/giap/schemas/config.schema.json").read_text())
    assert "rbs_novel" in json.dumps(schema)


# ---------------------------------------------------------------- load_precomputed


def _write(tmp_path: Path, rows: list[list[Any]], cols: tuple[str, ...] = COLS) -> Path:
    path = tmp_path / "rbs_novel.tsv"
    write_tsv(path, list(cols), rows)
    return path


def _load(path: Path) -> Any:
    from giap.annotate.precomputed import load_precomputed

    return load_precomputed("rbs_novel", path, REPS)


def test_precomputed_columns() -> None:
    from giap.annotate.precomputed import PRECOMPUTED_COLUMNS

    assert PRECOMPUTED_COLUMNS["rbs_novel"] == COLS


def test_load_converts_start_to_zero_based(tmp_path: Path) -> None:
    rows = [
        ["chr1", "+", 1, 2500.5, -9.1, -12.0, 0.3, -0.2, 7, "ATG"],
        ["chr1", "-", 1000, 12.25, -3.0, -8.5, 0.0, 0.0, 5, "GTG"],
        ["lin1", "+", 498, None, None, None, None, None, None, "TTG"],
    ]
    df = _load(_write(tmp_path, rows))
    assert list(df.columns) == list(COLS)
    assert df["start"].tolist() == [0, 999, 497]
    assert str(df["start"].dtype) == "int64"
    assert df["tir_ostir"].iloc[0] == 2500.5
    assert df["tir_ostir"].isna().iloc[2]
    assert df["spacer_bp"].iloc[0] == 7 and df["spacer_bp"].isna().iloc[2]
    assert df["start_codon"].tolist() == ["ATG", "GTG", "TTG"]


@pytest.mark.parametrize(
    ("row", "col"),
    [
        (["chr1", "+", 1001, 1.0, -1, -1, 0, 0, 5, "ATG"], "start"),  # beyond the contig
        (["chr1", "+", 0, 1.0, -1, -1, 0, 0, 5, "ATG"], "start"),
        (["chr1", "+", None, 1.0, -1, -1, 0, 0, 5, "ATG"], "start"),
        (["chrX", "+", 5, 1.0, -1, -1, 0, 0, 5, "ATG"], "contig"),
        (["chr1", "x", 5, 1.0, -1, -1, 0, 0, 5, "ATG"], "strand"),
        (["chr1", "+", 5, "high", -1, -1, 0, 0, 5, "ATG"], "tir_ostir"),
    ],
)
def test_load_rejects_bad_rows_naming_key_column_and_row(
    tmp_path: Path, row: list[Any], col: str
) -> None:
    good = ["chr1", "+", 10, 1.0, -1, -1, 0, 0, 5, "ATG"]
    with pytest.raises(ValueError, match=rf"rbs_novel.*row 2.*{col}"):
        _load(_write(tmp_path, [good, row]))


def test_load_missing_column_names_it(tmp_path: Path) -> None:
    cols = tuple(c for c in COLS if c != "start")
    with pytest.raises(ValueError, match=r"rbs_novel.*start"):
        _load(_write(tmp_path, [["chr1", "+", 1.0, -1, -1, 0, 0, 5, "ATG"]], cols))


# ---------------------------------------------------------------- parse_ostir_starts


def _ostir(tmp_path: Path, rows: list[list[Any]]) -> Path:
    path = tmp_path / "ostir.csv"
    lines = [OSTIR_HEADER] + [",".join(str(v) for v in r) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _row(name: str, codon: str, pos: int, expr: float, dist: int = 6) -> list[Any]:
    return [name, codon, pos, expr, dist, 1.1, -7.5, -12.3, 0.4, -0.1, -1.194]


def test_parse_ostir_starts_keeps_every_requested_start(tmp_path: Path) -> None:
    from giap.annotate.ostir import RbsFeatures, parse_ostir_starts

    path = _ostir(tmp_path, [
        _row("NOVEL_c1_10_p", "ATG", 51, 1500.25),
        _row("NOVEL_c1_10_p", "GTG", 54, 80.5, 9),   # requested
        _row("NOVEL_c1_10_p", "TTG", 66, 3.0),       # requested
        _row("NOVEL_c1_10_p", "ATG", 52, 99.0),      # not requested: dropped
        _row("NOVEL_c1_10_p", "GTG", 54, 1.0),       # repeat of (name, 54): first kept
        _row("NOVEL_c1_90_m", "ATG", 51, 7.0),
        _row("NOVEL_c1_500_p", "ATG", 51, 8.0),      # id not requested
    ])  # fmt: skip
    got = parse_ostir_starts(path, {"NOVEL_c1_10_p": {51, 54, 66}, "NOVEL_c1_90_m": {51, 57},
                                    "NOVEL_none_1_p": {51}})  # fmt: skip
    assert set(got) == {"NOVEL_c1_10_p", "NOVEL_c1_90_m"}
    assert set(got["NOVEL_c1_10_p"]) == {51, 54, 66}
    assert set(got["NOVEL_c1_90_m"]) == {51}  # no row at 57: absent (NA downstream)
    f = got["NOVEL_c1_10_p"][54]
    assert isinstance(f, RbsFeatures)
    assert (f.tir_ostir, f.spacer_bp, f.start_codon) == (80.5, 9, "GTG")
    assert (f.sd_dg, f.mrna_dg, f.spacing_dg, f.standby_dg) == (-7.5, -12.3, 0.4, -0.1)
    assert got["NOVEL_c1_10_p"][51].tir_ostir == 1500.25


def test_parse_ostir_starts_agrees_with_parse_ostir_on_single_starts(tmp_path: Path) -> None:
    from giap.annotate.ostir import parse_ostir, parse_ostir_starts

    path = _ostir(tmp_path, [_row("G1", "ATG", 51, 10.0), _row("G2", "ATG", 51, 20.0),
                             _row("G2", "ATG", 45, 1.0)])  # fmt: skip
    one = parse_ostir(path, {"G1": 51, "G2": 51})
    many = parse_ostir_starts(path, {"G1": {51}, "G2": {51}})
    assert {k: {51: v} for k, v in one.items()} == many


def test_parse_ostir_starts_header_only_is_empty(tmp_path: Path) -> None:
    from giap.annotate.ostir import parse_ostir_starts

    assert parse_ostir_starts(_ostir(tmp_path, []), {"NOVEL_c1_10_p": {51}}) == {}
