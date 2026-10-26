"""Sage ``results.sage.tsv`` -> generic peptide table (SPEC §9)."""

from __future__ import annotations

import csv
import gzip
import math
import tempfile
from pathlib import Path

import pandas as pd
import pytest
from giap_campaign.sage import convert_sage
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config.model import InputConfig
from giap.io.converters.base import ConvertResult
from giap.io.peptides import PEPTIDE_COLUMNS, validate_peptides

SAMPLE = Path(__file__).parent / "data" / "converters" / "results.sage.tsv"
SCAN = "controllerType=0 controllerNumber=1 scan="

# Real Sage 0.15 column order (results.sage.tsv).
SAGE_HEADER = [
    "peptide",
    "proteins",
    "num_proteins",
    "filename",
    "scannr",
    "rank",
    "label",
    "expmass",
    "calcmass",
    "charge",
    "peptide_len",
    "missed_cleavages",
    "semi_enzymatic",
    "isotope_error",
    "precursor_ppm",
    "fragment_ppm",
    "hyperscore",
    "delta_next",
    "delta_best",
    "rt",
    "aligned_rt",
    "predicted_rt",
    "delta_rt_model",
    "ion_mobility",
    "predicted_mobility",
    "delta_mobility",
    "matched_peaks",
    "longest_b",
    "longest_y",
    "longest_y_pct",
    "matched_intensity_pct",
    "scored_candidates",
    "poisson",
    "sage_discriminant_score",
    "posterior_error",
    "spectrum_q",
    "peptide_q",
    "protein_q",
    "ms2_intensity",
]

EXPECTED_COLUMNS = {
    "peptide",
    "modified_peptide",
    "protein_ids",
    "q_value",
    "pep",
    "score",
    "is_decoy",
    "run",
    "charge",
    "spectrum_id",
    "source",
    "search_db",
}


def _rows() -> list[dict[str, str]]:
    with SAMPLE.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write(path: Path, rows: list[dict[str, str]], header: list[str] | None = None) -> Path:
    header = header if header is not None else list(rows[0])
    lines = ["\t".join(header)] + ["\t".join(r[c] for c in header) for r in rows]
    text = "\n".join(lines) + "\n"
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def _without(path: Path, *columns: str) -> Path:
    rows = _rows()
    header = [c for c in SAGE_HEADER if c not in columns]
    return _write(path, rows, header)


def _by_sid(df: pd.DataFrame) -> dict[str, dict[str, object]]:
    return {r["spectrum_id"]: r for r in df.to_dict("records")}


def _errors(result: ConvertResult) -> list:
    return [i for i in result.issues if i.severity == "error"]


# ---------------------------------------------------------------- the fixture itself


def test_fixture_has_the_real_sage_header() -> None:
    with SAMPLE.open(encoding="utf-8") as fh:
        assert fh.readline().rstrip("\n").split("\t") == SAGE_HEADER


# ---------------------------------------------------------------- mapping


def test_convert_sage_mapping() -> None:
    result = convert_sage([SAMPLE])
    assert isinstance(result, ConvertResult)
    assert result.abundance is None
    assert _errors(result) == []
    df = result.peptides
    assert len(df) == 5  # the rank-2 row is dropped
    rows = _by_sid(df)

    r1 = rows[SCAN + "1001"]
    assert r1["peptide"] == "LVNELTEFAK"
    assert r1["modified_peptide"] == "LVNELTEFAK"
    assert r1["protein_ids"] == "b0014;GIAPORF_intergenic_chr_120_300_p"
    assert r1["q_value"] == pytest.approx(0.0005)  # peptide_q, not spectrum_q (0.001)
    assert r1["pep"] == pytest.approx(10**-3.2)
    assert r1["score"] == pytest.approx(1.234)  # sage_discriminant_score, not hyperscore
    assert not bool(r1["is_decoy"])
    assert r1["run"] == "run_A"
    assert r1["charge"] == 2
    assert r1["source"] == "database"
    assert r1["search_db"] == "same"

    r2 = rows[SCAN + "1002"]
    assert r2["peptide"] == "PEPTMIDEK"
    assert r2["modified_peptide"] == "PEPTM[+15.9949]IDEK"
    assert r2["protein_ids"] == "GIAPORF_chr_1_5000_5300_m"
    assert r2["q_value"] == pytest.approx(0.003)
    assert r2["pep"] == pytest.approx(0.1)
    assert r2["score"] == pytest.approx(0.876)
    assert r2["charge"] == 3
    assert r2["run"] == "run_A"

    decoy = rows[SCAN + "1003"]
    assert decoy["peptide"] == "EDIMTPEPK"
    assert decoy["protein_ids"] == "rev_GIAPORF_intergenic_chr_10_99_p"
    assert bool(decoy["is_decoy"]) is True
    assert decoy["pep"] == pytest.approx(1.0)  # 10 ** 0.0
    assert decoy["q_value"] == pytest.approx(0.45)
    assert decoy["score"] == pytest.approx(-0.45)

    r5 = rows[SCAN + "2001"]
    assert r5["peptide"] == "MSDKIIHLTDDSFDTDVLK"
    assert r5["modified_peptide"] == "[+42.0106]-MSDKIIHLTDDSFDTDVLK"
    assert r5["run"] == "run_B"  # directory and .mzML removed
    assert r5["pep"] == pytest.approx(1.0)  # 10 ** 0.3 clipped to 1
    assert r5["charge"] == 3

    r6 = rows[SCAN + "2002"]
    assert r6["peptide"] == "ACDEFGHIK"
    assert r6["modified_peptide"] == "AC[+57.0215]DEFGHIK"
    assert r6["protein_ids"] == "b0005;b0006"
    assert r6["run"] == "run_B"
    assert r6["pep"] == pytest.approx(10**-2.5)


def test_convert_sage_decoy_flags() -> None:
    df = convert_sage([SAMPLE]).peptides
    decoys = {sid: bool(d) for sid, d in zip(df["spectrum_id"], df["is_decoy"], strict=True)}
    assert decoys == {
        SCAN + "1001": False,
        SCAN + "1002": False,
        SCAN + "1003": True,
        SCAN + "2001": False,
        SCAN + "2002": False,
    }


def test_convert_sage_only_rank_1() -> None:
    df = convert_sage([SAMPLE]).peptides
    assert "AGHLLEVK" not in set(df["peptide"])
    assert list(df["spectrum_id"]).count(SCAN + "1002") == 1


def test_convert_sage_without_rank_column_keeps_every_row(tmp_path: Path) -> None:
    df = convert_sage([_without(tmp_path / "r.tsv", "rank")]).peptides
    assert len(df) == 6
    assert "AGHLLEVK" in set(df["peptide"])


def test_convert_sage_columns_and_order() -> None:
    df = convert_sage([SAMPLE]).peptides
    cols = list(df.columns)
    assert set(cols) == EXPECTED_COLUMNS
    assert cols == [c for c in PEPTIDE_COLUMNS if c in cols]


def test_convert_sage_row_order_and_types() -> None:
    df = convert_sage([SAMPLE]).peptides
    assert list(zip(df["run"], df["peptide"], strict=True)) == [
        ("run_A", "EDIMTPEPK"),
        ("run_A", "LVNELTEFAK"),
        ("run_A", "PEPTMIDEK"),
        ("run_B", "ACDEFGHIK"),
        ("run_B", "MSDKIIHLTDDSFDTDVLK"),
    ]
    assert pd.api.types.is_integer_dtype(df["charge"])
    assert pd.api.types.is_float_dtype(df["q_value"])
    assert pd.api.types.is_float_dtype(df["pep"])
    assert pd.api.types.is_float_dtype(df["score"])
    assert pd.api.types.is_bool_dtype(df["is_decoy"])
    assert set(df["source"]) == {"database"}
    assert set(df["search_db"]) == {"same"}
    assert all(isinstance(s, str) for s in df["spectrum_id"])


def test_convert_sage_q_value_falls_back_to_spectrum_q(tmp_path: Path) -> None:
    df = convert_sage([_without(tmp_path / "r.tsv", "peptide_q")]).peptides
    q = dict(zip(df["spectrum_id"], df["q_value"], strict=True))
    assert q[SCAN + "1001"] == pytest.approx(0.001)
    assert q[SCAN + "1002"] == pytest.approx(0.004)
    assert q[SCAN + "1003"] == pytest.approx(0.5)


def test_convert_sage_score_falls_back_to_hyperscore(tmp_path: Path) -> None:
    df = convert_sage([_without(tmp_path / "r.tsv", "sage_discriminant_score")]).peptides
    score = dict(zip(df["spectrum_id"], df["score"], strict=True))
    assert score == pytest.approx(
        {
            SCAN + "1001": 35.21,
            SCAN + "1002": 28.7,
            SCAN + "1003": 12.3,
            SCAN + "2001": 40.1,
            SCAN + "2002": 22.05,
        }
    )


def test_convert_sage_placeholder_posterior_error_leaves_pep_empty(tmp_path: Path) -> None:
    """Without LDA rescoring Sage writes posterior_error = 1.0 everywhere: pep is left empty."""
    rows = _rows()
    for r in rows:
        r["posterior_error"] = "1.0"
    result = convert_sage([_write(tmp_path / "r.tsv", rows)])
    df = result.peptides
    assert len(df) == 5
    assert "pep" in df.columns
    assert df["pep"].isna().all()
    assert _errors(result) == []


def test_convert_sage_mixed_values_are_clipped_not_placeholder(tmp_path: Path) -> None:
    rows = _rows()
    values = ["1.0", "-2.0", "-2.0", "1.5", "-4.0", "1.0"]
    for r, v in zip(rows, values, strict=True):
        r["posterior_error"] = v
    df = convert_sage([_write(tmp_path / "r.tsv", rows)]).peptides
    pep = dict(zip(df["spectrum_id"], df["pep"], strict=True))
    assert pep == pytest.approx(
        {
            SCAN + "1001": 1.0,
            SCAN + "1002": 0.01,
            SCAN + "1003": 1.0,
            SCAN + "2001": 1e-4,
            SCAN + "2002": 1.0,
        }
    )


@settings(max_examples=50, deadline=None)
@given(
    values=st.lists(
        st.floats(min_value=-12.0, max_value=3.0, allow_nan=False), min_size=6, max_size=6
    ).filter(lambda v: min(v[:2] + v[3:]) < 1.0)  # v[2] is the rank-2 row
)
def test_convert_sage_pep_is_clipped_power_of_ten(values: list[float]) -> None:
    rows = _rows()
    for r, v in zip(rows, values, strict=True):
        r["posterior_error"] = repr(v)
    with tempfile.TemporaryDirectory() as d:
        df = convert_sage([_write(Path(d) / "r.tsv", rows)]).peptides
    kept = [v for r, v in zip(rows, values, strict=True) if r["rank"] == "1"]
    sids = [r["scannr"] for r in rows if r["rank"] == "1"]
    pep = dict(zip(df["spectrum_id"], df["pep"], strict=True))
    for sid, v in zip(sids, kept, strict=True):
        assert 0.0 <= pep[sid] <= 1.0
        assert pep[sid] == pytest.approx(min(1.0, 10.0**v), rel=1e-9, abs=1e-300)


def test_convert_sage_gzip_input(tmp_path: Path) -> None:
    gz = _write(tmp_path / "results.sage.tsv.gz", _rows())
    plain = convert_sage([SAMPLE]).peptides
    pd.testing.assert_frame_equal(convert_sage([gz]).peptides, plain)


def test_convert_sage_several_inputs_are_concatenated(tmp_path: Path) -> None:
    rows = _rows()
    for r in rows:
        r["filename"] = r["filename"].replace("run_A", "run_C").replace("run_B", "run_D")
    second = _write(tmp_path / "second.tsv", rows)
    df = convert_sage([SAMPLE, second]).peptides
    assert len(df) == 10
    assert sorted(set(df["run"])) == ["run_A", "run_B", "run_C", "run_D"]
    assert list(df["run"]) == sorted(df["run"])


def test_convert_sage_accepts_str_paths() -> None:
    assert len(convert_sage([str(SAMPLE)]).peptides) == 5


def test_convert_sage_mzml_gz_and_plain_mzml_runs(tmp_path: Path) -> None:
    rows = _rows()[:2]
    rows[0]["filename"] = "/data/x/Ecoli_WT_1.mzML.gz"
    rows[1]["filename"] = "Ecoli_WT_2.mzML"
    df = convert_sage([_write(tmp_path / "r.tsv", rows)]).peptides
    assert sorted(df["run"]) == ["Ecoli_WT_1", "Ecoli_WT_2"]


# ---------------------------------------------------------------- issues


@pytest.mark.parametrize(
    "missing",
    [("peptide",), ("proteins",), ("filename",), ("scannr",), ("label",), ("spectrum_q",)],
)
def test_convert_sage_missing_required_column(tmp_path: Path, missing: tuple[str, ...]) -> None:
    result = convert_sage([_without(tmp_path / "r.tsv", *missing)])
    errors = _errors(result)
    assert len(errors) == 1
    assert errors[0].code.startswith("missing_column")
    assert missing[0] in errors[0].message
    assert len(result.peptides) == 0


def test_convert_sage_several_missing_columns_give_one_error(tmp_path: Path) -> None:
    result = convert_sage([_without(tmp_path / "r.tsv", "label", "spectrum_q")])
    errors = _errors(result)
    assert len(errors) == 1
    assert errors[0].code.startswith("missing_column")
    assert len(result.peptides) == 0


def test_convert_sage_optional_columns_are_not_required(tmp_path: Path) -> None:
    path = _without(
        tmp_path / "r.tsv", "peptide_q", "sage_discriminant_score", "posterior_error", "rank"
    )
    result = convert_sage([path])
    assert _errors(result) == []
    assert len(result.peptides) == 6


def test_convert_sage_empty_peptide_is_dropped_with_a_warning(tmp_path: Path) -> None:
    rows = _rows()
    rows[0]["peptide"] = "[+42.0106]-"
    rows[5]["peptide"] = "[+57.0215]"
    result = convert_sage([_write(tmp_path / "r.tsv", rows)])
    df = result.peptides
    assert len(df) == 3
    assert "" not in set(df["peptide"])
    warns = [i for i in result.issues if i.code == "empty_peptide"]
    assert len(warns) == 1
    assert warns[0].severity == "warn"
    assert "2" in warns[0].message
    assert _errors(result) == []


# ---------------------------------------------------------------- contract


@pytest.mark.pending_sample
def test_contract_hand_written_sample_is_a_valid_giap_peptide_table() -> None:
    """Hand-written sample until a recorded Sage 0.15.0 output exists (GIAP §5.4 practice)."""
    result = convert_sage([SAMPLE])
    assert _errors(result) == []
    df = result.peptides
    valid, issues = validate_peptides(df, InputConfig())
    assert [i for i in issues if i.severity == "error"] == []
    assert len(valid) == len(df)
    assert not any(math.isnan(q) for q in valid["q_value"])
