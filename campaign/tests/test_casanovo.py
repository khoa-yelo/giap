"""Casanovo mzTab -> generic peptide table (SPEC §9)."""

from __future__ import annotations

import gzip
import math
from pathlib import Path

import pandas as pd
import pytest
from giap_campaign.casanovo import convert_casanovo

from giap.config.model import InputConfig
from giap.io.converters.base import ConvertResult
from giap.io.peptides import PEPTIDE_COLUMNS, validate_peptides

SAMPLE = Path(__file__).parent / "data" / "converters" / "casanovo.mztab"

EXPECTED_COLUMNS = {
    "peptide",
    "modified_peptide",
    "score",
    "run",
    "charge",
    "spectrum_id",
    "source",
    "search_db",
}

PSH = (
    "PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\t"
    "search_engine_score[1]\tmodifications\tretention_time\tcharge\texp_mass_to_charge\t"
    "calc_mass_to_charge\tspectra_ref\tpre\tpost\tstart\tend\topt_ms_run[1]_aa_scores"
)
ENGINE = "[MS, MS:1003281, Casanovo, 5.2.1]"


def _psm(seq: str, score: str, charge: str, ref: str, psm_id: int = 1) -> str:
    cells = [
        "PSM",
        seq,
        str(psm_id),
        "null",
        "null",
        "null",
        "null",
        ENGINE,
        score,
        "null",
        "null",
        charge,
        "500.0",
        "500.0",
        ref,
        "null",
        "null",
        "null",
        "null",
        "null",
    ]
    return "\t".join(cells)


def _mztab(locations: list[str], psms: list[str], psh: str = PSH) -> str:
    mtd = [
        "MTD\tmzTab-version\t1.0.0",
        "MTD\tmzTab-mode\tSummary",
        "MTD\tmzTab-type\tIdentification",
        f"MTD\tsoftware[1]\t{ENGINE}",
        "MTD\tpsm_search_engine_score[1]\t"
        "[MS, MS:1001143, search engine specific score for PSMs, ]",
    ]
    mtd += [f"MTD\tms_run[{i}]-location\t{loc}" for i, loc in enumerate(locations, start=1)]
    return "\n".join([*mtd, psh, *psms]) + "\n"


def _write(path: Path, text: str) -> Path:
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def _errors(result: ConvertResult) -> list:
    return [i for i in result.issues if i.severity == "error"]


def _key(df: pd.DataFrame) -> dict[tuple[str, str], dict[str, object]]:
    return {(r["run"], r["spectrum_id"]): r for r in df.to_dict("records")}


# ---------------------------------------------------------------- mapping


def test_convert_casanovo_mapping() -> None:
    result = convert_casanovo([SAMPLE])
    assert isinstance(result, ConvertResult)
    assert result.abundance is None
    assert _errors(result) == []
    df = result.peptides
    assert len(df) == 5
    rows = _key(df)
    assert set(rows) == {
        ("run_A", "index=0"),
        ("run_A", "index=4"),
        ("run_B", "index=3"),
        ("run_B", "index=7"),
        ("run_B", "index=9"),
    }

    plain = rows[("run_A", "index=0")]
    assert plain["peptide"] == "PEPTIDEK"
    assert plain["modified_peptide"] == "PEPTIDEK"
    assert plain["score"] == pytest.approx(0.9312)
    assert plain["charge"] == 2
    assert plain["source"] == "denovo"
    assert plain["search_db"] == "same"

    acetyl = rows[("run_A", "index=4")]  # N-terminal mass shift
    assert acetyl["peptide"] == "PEPTIDEK"
    assert acetyl["modified_peptide"] == "+42.011PEPTIDEK"
    assert acetyl["score"] == pytest.approx(0.8120)

    ox = rows[("run_B", "index=3")]  # .mzML.gz location: .gz then .mzML removed
    assert ox["peptide"] == "MLVNELTEFAK"
    assert ox["modified_peptide"] == "M[Oxidation]LVNELTEFAK"
    assert ox["charge"] == 3
    assert ox["score"] == pytest.approx(0.7755)

    cam = rows[("run_B", "index=7")]  # score null -> missing
    assert cam["peptide"] == "ACDEFGHIK"
    assert cam["modified_peptide"] == "AC[Carbamidomethyl]DEFGHIK"
    assert cam["score"] is None or pd.isna(cam["score"])

    neg = rows[("run_B", "index=9")]  # Casanovo 4 notation; negative score kept
    assert neg["peptide"] == "VLSEAMK"
    assert neg["modified_peptide"] == "VLSEAM+15.995K"
    assert neg["score"] == pytest.approx(-0.41)


def test_convert_casanovo_columns_order_and_types() -> None:
    df = convert_casanovo([SAMPLE]).peptides
    cols = list(df.columns)
    assert set(cols) == EXPECTED_COLUMNS
    assert cols == [c for c in PEPTIDE_COLUMNS if c in cols]
    assert pd.api.types.is_integer_dtype(df["charge"])
    assert pd.api.types.is_float_dtype(df["score"])
    assert set(df["source"]) == {"denovo"}
    assert set(df["search_db"]) == {"same"}
    # sorted by run, peptide, modified_peptide, charge
    assert list(zip(df["run"], df["modified_peptide"], strict=True)) == [
        ("run_A", "+42.011PEPTIDEK"),
        ("run_A", "PEPTIDEK"),
        ("run_B", "AC[Carbamidomethyl]DEFGHIK"),
        ("run_B", "M[Oxidation]LVNELTEFAK"),
        ("run_B", "VLSEAM+15.995K"),
    ]


def test_convert_casanovo_null_cells_are_missing(tmp_path: Path) -> None:
    text = _mztab(
        ["file:///data/run_A.mzML"],
        [
            _psm("PEPTIDEK", "null", "2", "ms_run[1]:index=0", 1),
            _psm("LVNELTEFAK", "0.5", "null", "ms_run[1]:index=1", 2),
        ],
    )
    df = convert_casanovo([_write(tmp_path / "n.mztab", text)]).peptides
    rows = _key(df)
    assert pd.isna(rows[("run_A", "index=0")]["score"])
    assert pd.isna(rows[("run_A", "index=1")]["charge"])
    assert rows[("run_A", "index=1")]["score"] == pytest.approx(0.5)
    assert "null" not in set(df["modified_peptide"].astype(str))


def test_convert_casanovo_run_names_from_locations(tmp_path: Path) -> None:
    text = _mztab(
        [
            "file:///data/raw/Ecoli_WT_1.mzML",
            "file:///scratch/x/Ecoli_WT_2.mzML.gz",
            "Ecoli_WT_3.mzml",
        ],
        [
            _psm("PEPTIDEK", "0.9", "2", "ms_run[1]:index=0", 1),
            _psm("PEPTIDEK", "0.9", "2", "ms_run[2]:index=0", 2),
            _psm("PEPTIDEK", "0.9", "2", "ms_run[3]:index=0", 3),
        ],
    )
    df = convert_casanovo([_write(tmp_path / "r.mztab", text)]).peptides
    assert list(df["run"]) == ["Ecoli_WT_1", "Ecoli_WT_2", "Ecoli_WT_3"]
    assert set(df["spectrum_id"]) == {"index=0"}


def test_convert_casanovo_scan_style_spectra_ref(tmp_path: Path) -> None:
    text = _mztab(
        ["file:///data/run_A.mzML"],
        [_psm("PEPTIDEK", "0.9", "2", "ms_run[1]:controllerType=0 controllerNumber=1 scan=17")],
    )
    df = convert_casanovo([_write(tmp_path / "r.mztab", text)]).peptides
    assert list(df["spectrum_id"]) == ["controllerType=0 controllerNumber=1 scan=17"]
    assert list(df["run"]) == ["run_A"]


def test_convert_casanovo_gzip_input(tmp_path: Path) -> None:
    gz = _write(tmp_path / "denovo.mztab.gz", SAMPLE.read_text(encoding="utf-8"))
    pd.testing.assert_frame_equal(
        convert_casanovo([gz]).peptides, convert_casanovo([SAMPLE]).peptides
    )


def test_convert_casanovo_several_inputs_are_concatenated(tmp_path: Path) -> None:
    # ms_run[1] of the second file is a different run: indices are per file
    second = _mztab(
        ["file:///data/run_C.mzML"],
        [
            _psm("SAMPLEK", "0.66", "2", "ms_run[1]:index=5", 1),
            _psm("M[Oxidation]SAMPLER", "0.61", "3", "ms_run[1]:index=6", 2),
        ],
    )
    path = _write(tmp_path / "second.mztab", second)
    result = convert_casanovo([SAMPLE, path])
    assert _errors(result) == []
    df = result.peptides
    assert len(df) == 7
    rows = _key(df)
    assert rows[("run_C", "index=5")]["peptide"] == "SAMPLEK"
    assert rows[("run_C", "index=6")]["peptide"] == "MSAMPLER"
    assert ("run_A", "index=5") not in rows
    assert sorted(set(df["run"])) == ["run_A", "run_B", "run_C"]


def test_convert_casanovo_accepts_str_paths() -> None:
    assert len(convert_casanovo([str(SAMPLE)]).peptides) == 5


# ---------------------------------------------------------------- issues


@pytest.mark.parametrize("missing", ["sequence", "spectra_ref"])
def test_convert_casanovo_missing_required_column(tmp_path: Path, missing: str) -> None:
    header = PSH.split("\t")
    j = header.index(missing)
    psh = "\t".join(header[:j] + header[j + 1 :])
    psm = _psm("PEPTIDEK", "0.9", "2", "ms_run[1]:index=0").split("\t")
    psm = psm[:j] + psm[j + 1 :]
    text = _mztab(["file:///data/run_A.mzML"], ["\t".join(psm)], psh=psh)
    result = convert_casanovo([_write(tmp_path / "m.mztab", text)])
    errors = _errors(result)
    assert len(errors) == 1
    assert errors[0].code.startswith("missing_column")
    assert missing in errors[0].message
    assert len(result.peptides) == 0


def test_convert_casanovo_empty_peptide_is_dropped_with_a_warning(tmp_path: Path) -> None:
    text = _mztab(
        ["file:///data/run_A.mzML"],
        [
            _psm("PEPTIDEK", "0.9", "2", "ms_run[1]:index=0", 1),
            _psm("+42.011", "0.2", "2", "ms_run[1]:index=1", 2),
            _psm("[Acetyl]-", "0.1", "2", "ms_run[1]:index=2", 3),
        ],
    )
    result = convert_casanovo([_write(tmp_path / "e.mztab", text)])
    assert list(result.peptides["peptide"]) == ["PEPTIDEK"]
    warns = [i for i in result.issues if i.code == "empty_peptide"]
    assert len(warns) == 1
    assert warns[0].severity == "warn"
    assert "2" in warns[0].message
    assert _errors(result) == []


def test_convert_casanovo_header_only_gives_no_rows(tmp_path: Path) -> None:
    text = _mztab(["file:///data/run_A.mzML"], [])
    result = convert_casanovo([_write(tmp_path / "h.mztab", text)])
    assert len(result.peptides) == 0
    assert not [i for i in _errors(result) if i.code.startswith("missing_column")]


# ---------------------------------------------------------------- contract


def test_contract_sample_is_a_valid_giap_peptide_table() -> None:
    df = convert_casanovo([SAMPLE]).peptides
    valid, issues = validate_peptides(df, InputConfig())
    assert [i for i in issues if i.severity == "error"] == []
    assert len(valid) == len(df)
    assert all(isinstance(p, str) and p.isupper() and p.isalpha() for p in valid["peptide"])
    assert not any(isinstance(s, float) and math.isnan(s) for s in valid["run"])
