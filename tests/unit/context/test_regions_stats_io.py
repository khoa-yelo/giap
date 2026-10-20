"""M4-18 I/O: ``read_replicon_depth``, ``write_regions_stats`` and the ``regions_stats`` schema.

docs/specs/M4-18-regions-stats-table.md "Public interface", edge cases (depth <= 0 is an error
issue; header-only table; origin-crossing regions keep ``end > length`` in the 1-based TSV) and
acceptance (schema; flags are a controlled vocabulary "defined in the schema", GIAP_SPEC.md §9.2,
§9.4, §9.5). §5.3: the depth TSV has columns ``contig``, ``mean_depth``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import regions_stats_data as sd

TABLE_DIR = Path(__file__).resolve().parents[3] / "src" / "giap" / "schemas" / "tables"


def _errors(issues: list[Any]) -> list[Any]:
    return [i for i in issues if i.severity == "error"]


# --------------------------------------------------------------------------- depth reader


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_read_replicon_depth_valid(tmp_path: Path) -> None:
    from giap.io.depth import read_replicon_depth

    p = _write(tmp_path / "depth.tsv", "contig\tmean_depth\nchr1\t50\np1\t400.5\n")
    df, issues = read_replicon_depth(p)
    assert not _errors(issues)
    assert df is not None
    assert list(df["contig"].astype(str)) == ["chr1", "p1"]
    np.testing.assert_allclose(df["mean_depth"].to_numpy(float), [50.0, 400.5])


@pytest.mark.parametrize("value", ["0", "-3.5"])
def test_read_replicon_depth_nonpositive_is_error(tmp_path: Path, value: str) -> None:
    """Edge case: depth <= 0 -> error issue (log2 undefined)."""
    from giap.io.depth import read_replicon_depth

    p = _write(tmp_path / "depth.tsv", f"contig\tmean_depth\nchr1\t50\np1\t{value}\n")
    _, issues = read_replicon_depth(p)
    assert _errors(issues)


def test_read_replicon_depth_missing_column_is_error(tmp_path: Path) -> None:
    from giap.io.depth import read_replicon_depth

    p = _write(tmp_path / "depth.tsv", "contig\tdepth\nchr1\t50\n")
    df, issues = read_replicon_depth(p)
    assert _errors(issues)
    assert df is None


def test_read_replicon_depth_non_numeric_is_error(tmp_path: Path) -> None:
    from giap.io.depth import read_replicon_depth

    p = _write(tmp_path / "depth.tsv", "contig\tmean_depth\nchr1\thigh\n")
    _, issues = read_replicon_depth(p)
    assert _errors(issues)


@pytest.mark.parametrize("value", ["inf", "-inf", "nan", "1e400"])
def test_read_replicon_depth_non_finite_is_error(tmp_path: Path, value: str) -> None:
    """Decision 19: non-finite depths are errors and give no frame."""
    from giap.io.depth import read_replicon_depth

    p = _write(tmp_path / "depth.tsv", f"contig\tmean_depth\nchr1\t50\np1\t{value}\n")
    df, issues = read_replicon_depth(p)
    assert _errors(issues) and df is None


@pytest.mark.parametrize(
    "text",
    ["contig\tmean_depth\nchr1\t50\np1\t400\nchr1\t60\n",
     "contig\tmean_depth\nchr1\t50\nNA\t400\n",
     "contig\tmean_depth\nchr1\t50\n\t400\n"],
    ids=["duplicate", "na_contig", "empty_contig"],
)  # fmt: skip
def test_read_replicon_depth_bad_contig_is_error(tmp_path: Path, text: str) -> None:
    from giap.io.depth import read_replicon_depth

    df, issues = read_replicon_depth(_write(tmp_path / "depth.tsv", text))
    assert _errors(issues) and df is None


@pytest.mark.parametrize(
    "text", ["", "contig\tmean_depth\nchr1\t50\tEXTRA\n", "chr1\n50\t3\t4\n"],
    ids=["empty", "ragged", "ragged_header"],
)  # fmt: skip
def test_read_replicon_depth_unreadable_returns_issues(tmp_path: Path, text: str) -> None:
    """Decision 19: empty / ragged files give (None, issues), never an exception."""
    from giap.io.depth import read_replicon_depth

    df, issues = read_replicon_depth(_write(tmp_path / "depth.tsv", text))
    assert df is None and _errors(issues)


def test_read_replicon_depth_missing_file_returns_issues(tmp_path: Path) -> None:
    from giap.io.depth import read_replicon_depth

    df, issues = read_replicon_depth(tmp_path / "absent.tsv")
    assert df is None and _errors(issues)


# --------------------------------------------------------------------------- writer


def _origin_genome() -> tuple[sd.Genome, int]:
    g = sd.mixed_genome()
    length = int(g.replicons.set_index("contig").loc["chr1", "length"])
    j = g.regions.index[g.regions["region_id"] == "crispr_1"][0]
    g.regions.loc[j, "start"] = length - 300
    g.regions.loc[j, "end"] = length + 200  # crosses the origin: end > length (§6.5)
    return g, length


def _written(tmp_path: Path, df: pd.DataFrame) -> tuple[list[str], list[dict[str, str | None]]]:
    from giap.io.regions_tsv import write_regions_stats
    from giap.io.tsv import read_tsv

    out = tmp_path / "regions_stats.tsv"
    write_regions_stats(out, df)
    header, rows = read_tsv(out)
    return header, [dict(zip(header, r, strict=True)) for r in rows]


def test_writer_header_and_one_based_coordinates(tmp_path: Path) -> None:
    g, length = _origin_genome()
    df = sd.table(g)
    header, rows = _written(tmp_path, df)
    assert tuple(header) == sd.REGIONS_STATS_COLUMNS
    assert len(rows) == len(df)
    c7 = g.regions.set_index("region_id")
    for r in rows:
        rid = str(r["region_id"])
        if rid == "amr_1":
            assert r["start"] is None and r["end"] is None  # NA
            continue
        assert int(str(r["start"])) == int(c7.loc[rid, "start"]) + 1
        assert int(str(r["end"])) == int(c7.loc[rid, "end"])
    cr = [r for r in rows if r["region_id"] == "crispr_1"]
    assert cr and all(int(str(r["end"])) == length + 200 > length for r in cr)


def test_writer_na_and_flags(tmp_path: Path) -> None:
    df = sd.table(sd.mixed_genome(), sd.depth_frame({"chr1": 40.0, "p1": 160.0}))
    _, rows = _written(tmp_path, df)
    pl = [r for r in rows if r["region_type"] == "plasmid"]
    assert len(pl) == 3
    for r in pl:
        assert r["p_value"] is None and r["q_value"] is None
        assert r["n_perm"] == "0" and r["null_model"] == "none"
        assert (
            "dosage_plus_regulation" in str(r["flags"]).split(";")
            or r["label"] == "dosage_adjusted"
        )
    for r, (_, src) in zip(rows, df.iterrows(), strict=True):
        assert set(filter(None, str(r["flags"] or "").split(";"))) == sd.flags(src["flags"])
        for c in ("effect", "ci_low", "ci_high", "observed", "expected"):
            if sd.na(src[c]):
                assert r[c] is None, c
            else:
                assert float(str(r[c])) == pytest.approx(float(src[c]), rel=1e-12, abs=1e-12), c


def test_writer_empty_table_is_header_only(tmp_path: Path) -> None:
    from giap.io.regions_tsv import write_regions_stats

    T = sd.module()
    g = sd.genome([sd.chr_contig()])
    df = T.regions_stats(sd.tests(g), sd.plasmid_rows(g), T.subset_rows(g.genes, g.regions,
                         g.region_genes), g.regions)  # fmt: skip
    out = tmp_path / "regions_stats.tsv"
    write_regions_stats(out, df)
    assert out.read_text(encoding="utf-8") == "\t".join(sd.REGIONS_STATS_COLUMNS) + "\n"


# --------------------------------------------------------------------------- schema


def _descriptor() -> dict[str, Any]:
    found = [p for p in (TABLE_DIR / "regions_stats.json",
                         TABLE_DIR / "regions_stats.schema.json") if p.exists()]  # fmt: skip
    assert len(found) == 1, found
    return json.loads(found[0].read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def _field(d: dict[str, Any], name: str) -> dict[str, Any]:
    return next(f for f in d["fields"] if f["name"] == name)


def test_schema_fields_match_columns() -> None:
    d = _descriptor()
    assert d["name"] == "regions_stats"
    assert d.get("coordinates") == "1-based-closed"
    assert tuple(f["name"] for f in d["fields"]) == sd.REGIONS_STATS_COLUMNS
    for f in d["fields"]:
        assert f.get("type") and f.get("description"), f["name"]


def test_schema_enums() -> None:
    d = _descriptor()
    assert set(_field(d, "stat")["constraints"]["enum"]) == {"o_minus_e", "hl_shift"}
    assert {"none", "circular_shift"} <= set(_field(d, "null_model")["constraints"]["enum"])
    assert set(_field(d, "subset")["constraints"]["enum"]) == {"all", *sd.SUBSETS}
    assert set(_field(d, "region_type")["constraints"]["enum"]) == set(sd.TYPE_ORDER)


def _validate(path: Path) -> Any:
    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in _descriptor().items() if k != "coordinates"}
    resource = frictionless.Resource(path=path.name, basepath=str(path.parent),
                                     schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    return resource.validate()


def test_written_table_validates(tmp_path: Path) -> None:
    from giap.io.regions_tsv import write_regions_stats

    g, _ = _origin_genome()
    out = tmp_path / "regions_stats.tsv"
    write_regions_stats(out, sd.table(g, sd.depth_frame({"chr1": 40.0, "p1": 160.0})))
    report = _validate(out)
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def _with_flags(tmp_path: Path, value: str) -> Path:
    from giap.io.regions_tsv import write_regions_stats
    from giap.io.tsv import read_tsv, write_tsv

    src = tmp_path / "src.tsv"
    write_regions_stats(src, sd.table(sd.mixed_genome()))
    header, rows = read_tsv(src)
    j = header.index("flags")
    rows[0][j] = value
    out = tmp_path / f"flags_{abs(hash(value))}.tsv"
    write_tsv(out, header, rows)
    return out


def test_unknown_flag_fails_validation(tmp_path: Path) -> None:
    """§9.2 / §9.4: flags are a controlled vocabulary defined in the schema."""
    assert not _validate(_with_flags(tmp_path, "bogus_flag")).valid
    assert not _validate(_with_flags(tmp_path, "small_region;bogus_flag")).valid


@pytest.mark.parametrize("flag", sd.REGION_FLAGS)
def test_vocabulary_flags_validate(tmp_path: Path, flag: str) -> None:
    report = _validate(_with_flags(tmp_path, flag))
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_flag_lists_validate(tmp_path: Path) -> None:
    report = _validate(_with_flags(tmp_path, "small_region;few_shifts;split_contigs"))
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])
