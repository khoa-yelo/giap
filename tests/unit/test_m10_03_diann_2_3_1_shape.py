"""DIA-NN 2.3.1-shaped parquet reports through ``convert_diann`` (M10-03 criteria 2-4).

Spec: docs/specs/M10-03-diann-convert.md, examples 1-7 and facts F1-F10. Every fixture is
hand-made (``helpers.diann_convert_m10``) from the DIA-NN 2.x documentation: float32 scores and
intensities, int8 or bool ``Decoy``, dictionary-encoded text, a few unmapped extra columns,
``Cont_`` IDs and A2 locus tags. No value, row or column list of the real report is used (Q-K2).

Conditions the real sample may show later (a renamed column, ``sp|ACC|NAME`` IDs) get their own
CI test once seen (spec "Test plan", step 2); none is assumed here.
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from click.testing import CliRunner
from helpers.diann_convert_m10 import row, table, write_parquet
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.cli import cli
from gippa.config.model import InputConfig
from gippa.io.converters.base import ConvertResult
from gippa.io.converters.diann import convert_diann
from gippa.io.peptides import PEPTIDE_COLUMNS

# Example 1 (A1 target) and example 2 (its decoy), in the 2.3.1 shape.
EX1 = row(
    Run="WT_1",
    Stripped_Sequence="PEPTMIDEK",
    Modified_Sequence="PEPTM(UniMod:35)IDEK",
    Protein_Ids="P0A7V0",
    Protein_Group="P0A7V0",
    Global_Q_Value=0.0021,
    Q_Value=0.0009,
    PEP=0.003,
    Decoy=0,
    Precursor_Charge=2,
    Precursor_Normalised=15234.5,
    Precursor_Quantity=14000.0,
    PG_MaxLFQ=100.0,
)
EX2 = {
    **EX1,
    "Stripped.Sequence": "PLPTMIDDK",
    "Modified.Sequence": "PLPTM(UniMod:35)IDDK",
    "Decoy": 1,
    "Global.Q.Value": 0.4,
    "PEP": 0.97,
}
EXTRA = {"Precursor.Id": "PEPTM(UniMod:35)IDEK2", "RT": 12.5, "IM": 0.91}


def _convert(
    tmp_path: Path,
    rows: list[dict[str, Any]],
    search_db: str | None = "other",
    **types: pa.DataType,
) -> ConvertResult:
    path = write_parquet(tmp_path / "report.parquet", table(rows, types))
    return convert_diann([path], search_db=search_db, cfg=InputConfig())


def _bad(res: ConvertResult) -> list[str]:
    return [f"{i.severity}:{i.code}" for i in res.issues if i.severity in ("error", "warn")]


def _f32(x: float) -> float:
    return float(np.float32(x))


def _na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


# ------------------------------------------------------------------ example 1: the preference


def test_example1_global_q_and_normalised_intensity_win(tmp_path: Path) -> None:
    """F1: ``Global.Q.Value`` (not ``Q.Value``) and ``Precursor.Normalised`` (not ``.Quantity``)."""
    res = _convert(tmp_path, [{**EX1, **EXTRA}])
    assert not _bad(res), _bad(res)
    r = res.peptides.iloc[0]
    assert r["q_value"] == _f32(0.0021)
    assert r["intensity"] == _f32(15234.5)
    assert r["pep"] == _f32(0.003)
    assert (r["peptide"], r["modified_peptide"], r["protein_ids"]) == (
        "PEPTMIDEK",
        "PEPTM(UniMod:35)IDEK",
        "P0A7V0",
    )
    assert (r["run"], int(r["charge"]), bool(r["is_decoy"])) == ("WT_1", 2, False)
    assert (r["source"], r["search_db"]) == ("database", "other")


def test_example1_written_line_and_column_order(tmp_path: Path) -> None:
    """Criterion 3: §5.3 column order, ``%.6g`` floats (``0.0021``, not ``0.0020999999``)."""
    rep = write_parquet(tmp_path / "report.parquet", table([{**EX1, **EXTRA}]))
    out = tmp_path / "a.tsv"
    res = CliRunner().invoke(
        cli, ["convert", "--from", "diann", str(rep), "--out", str(out), "--search-db", "other"]
    )
    assert res.exit_code == 0, res.output
    header, line = out.read_text().splitlines()
    cols = header.split("\t")
    assert cols == [c for c in PEPTIDE_COLUMNS if c in cols]
    assert "score" not in cols and "spectrum_id" not in cols
    assert line.split("\t") == [
        "PEPTMIDEK",
        "PEPTM(UniMod:35)IDEK",
        "P0A7V0",
        "0.0021",
        "0.003",
        "false",
        "15234.5",
        "WT_1",
        "2",
        "database",
        "other",
    ]


def test_parquet_dtypes_of_the_generic_table(tmp_path: Path) -> None:
    """F2: float32 widened to float64, int8 charge to Int64, int8 ``Decoy`` to boolean."""
    df = _convert(tmp_path, [EX1, EX2]).peptides
    for col in ("q_value", "pep", "intensity"):
        assert df[col].dtype == np.float64, (col, df[col].dtype)
    assert str(df["charge"].dtype) == "Int64"
    assert str(df["is_decoy"].dtype) == "boolean"
    for col in ("peptide", "modified_peptide", "protein_ids", "run"):
        assert isinstance(df[col].dtype, pd.StringDtype), (col, df[col].dtype)


# ------------------------------------------------------------------ example 2: decoys kept


@pytest.mark.parametrize("decoy_type", [pa.int8(), pa.bool_(), pa.int64()], ids=str)
def test_example2_decoy_rows_kept_for_every_decoy_type(
    tmp_path: Path, decoy_type: pa.DataType
) -> None:
    """F2/F4: ``Decoy`` as int8 or ``bool`` gives ``is_decoy``; q filtering is preflight's."""
    rows = [EX1, EX2]
    if decoy_type == pa.bool_():
        rows = [{**r, "Decoy": bool(r["Decoy"])} for r in rows]
    res = _convert(tmp_path, rows, Decoy=decoy_type)
    assert not _bad(res), _bad(res)
    df = res.peptides.set_index("peptide")
    assert len(df) == 2
    assert bool(df.loc["PLPTMIDDK", "is_decoy"]) is True
    assert bool(df.loc["PEPTMIDEK", "is_decoy"]) is False
    assert df.loc["PLPTMIDDK", "protein_ids"] == "P0A7V0"
    assert df.loc["PLPTMIDDK", "q_value"] == _f32(0.4)


# ------------------------------------------------------------------ example 3: NA q on a decoy


def test_example3_decoy_with_na_global_q_has_na_q_value(tmp_path: Path) -> None:
    """Criterion 4 (Q2): no per-row ``Q.Value`` fallback; targets keep ``Global.Q.Value``."""
    decoy = {**EX2, "Global.Q.Value": None, "Q.Value": 0.3}
    res = _convert(tmp_path, [EX1, decoy])
    assert not _bad(res), _bad(res)
    df = res.peptides.set_index("peptide")
    assert _na(df.loc["PLPTMIDDK", "q_value"])
    assert df.loc["PLPTMIDDK", "pep"] == _f32(0.97)
    assert df.loc["PEPTMIDEK", "q_value"] == _f32(0.0021)


# ------------------------------------------------------------------ examples 4-6: IDs, mods


def test_example4_and_5_cont_prefixes_survive(tmp_path: Path) -> None:
    """F8: ``Cont_P00761`` and ``Cont_P00722;P00722`` reach ``protein_ids`` unchanged."""
    rows = [
        {**EX1, "Protein.Ids": "Cont_P00761", "Protein.Group": "Cont_P00761"},
        {
            **EX1,
            "Stripped.Sequence": "LLGNVLVR",
            "Modified.Sequence": "LLGNVLVR",
            "Protein.Ids": "Cont_P00722;P00722",
            "Protein.Group": "Cont_P00722;P00722",
        },
    ]
    df = _convert(tmp_path, rows).peptides.set_index("peptide")
    assert df.loc["PEPTMIDEK", "protein_ids"] == "Cont_P00761"
    assert df.loc["LLGNVLVR", "protein_ids"] == "Cont_P00722;P00722"


def test_example6_a2_locus_tag_fixed_mod_and_zero_intensity(tmp_path: Path) -> None:
    a2 = {
        **EX1,
        "Protein.Ids": "JMJJOG_01234",
        "Protein.Group": "JMJJOG_01234",
        "Modified.Sequence": "AC(UniMod:4)DLM(UniMod:35)TGK",
        "Stripped.Sequence": "ACDLMTGK",
        "Precursor.Normalised": 0.0,
    }
    res = _convert(tmp_path, [a2], search_db="same")
    assert not _bad(res), _bad(res)
    r = res.peptides.iloc[0]
    assert (r["peptide"], r["modified_peptide"]) == ("ACDLMTGK", "AC(UniMod:4)DLM(UniMod:35)TGK")
    assert r["protein_ids"] == "JMJJOG_01234"
    assert _na(r["intensity"])
    assert r["search_db"] == "same"


def test_dictionary_encoded_text_columns(tmp_path: Path) -> None:
    """F2: ``Run`` and ``Protein.Ids`` dictionary-encoded (pandas ``category``) convert as text."""
    rows = [EX1, {**EX1, "Run": "WT_2", "Protein.Ids": "P0A7V1"}, EX2]
    res = _convert(tmp_path, rows)
    assert not _bad(res), _bad(res)
    df = res.peptides
    assert sorted(df["run"]) == ["WT_1", "WT_1", "WT_2"]
    assert sorted(df["protein_ids"]) == ["P0A7V0", "P0A7V0", "P0A7V1"]
    assert isinstance(df["run"].dtype, pd.StringDtype)


# ------------------------------------------------------------------ example 7: MaxLFQ


def test_example7_maxlfq_median_over_runs(tmp_path: Path) -> None:
    """F9: one value per run, 0 excluded, decoys excluded, median of 100 and 300 = 200."""
    g = {"Protein.Ids": "JMJJOG_01234", "Protein.Group": "JMJJOG_01234"}
    peps = ["ACDLMTGK", "LLGNVLVR", "IGGIGTVPVGR"]
    rows = [
        {**EX1, **g, "Stripped.Sequence": p, "Modified.Sequence": p, "PG.MaxLFQ": 100.0}
        for p in peps
    ]
    rows += [
        {**EX1, **g, "Run": "WT_2", "PG.MaxLFQ": 300.0},
        {**EX1, **g, "Run": "WT_3", "PG.MaxLFQ": 0.0},
        {**EX2, **g, "Run": "WT_3", "PG.MaxLFQ": 900.0},
    ]
    res = _convert(tmp_path, rows, search_db="same")
    ab = res.abundance
    assert ab is not None
    assert list(ab.columns) == ["protein_id", "maxlfq"]
    assert ab["protein_id"].tolist() == ["JMJJOG_01234"]
    assert ab["maxlfq"].tolist() == [200.0]


def test_maxlfq_with_bool_decoy_excludes_decoys(tmp_path: Path) -> None:
    """A ``bool`` ``Decoy`` column must still exclude decoy rows from the abundance table."""
    rows = [
        {**EX1, "Decoy": False, "PG.MaxLFQ": 100.0},
        {**EX2, "Decoy": True, "Protein.Group": "P0A7V9", "PG.MaxLFQ": 900.0},
    ]
    ab = _convert(tmp_path, rows, Decoy=pa.bool_()).abundance
    assert ab is not None
    assert ab["protein_id"].tolist() == ["P0A7V0"]


# ------------------------------------------------------------------ F10: one row per report row

_PEP = st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=7, max_size=12)
_ROW = st.fixed_dictionaries(
    {
        "Run": st.sampled_from(["WT_1", "WT_2", "WT_3"]),
        "Stripped.Sequence": _PEP,
        "Protein.Ids": st.sampled_from(["P0A7V0", "Cont_P00761", "JMJJOG_00001;JMJJOG_00002"]),
        "Global.Q.Value": st.one_of(st.none(), st.floats(0, 1, width=32)),
        "Decoy": st.integers(0, 1),
        "Precursor.Charge": st.integers(1, 4),
        "Precursor.Normalised": st.one_of(st.none(), st.floats(0, 1e9, width=32)),
    }
)


@settings(
    max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(rows=st.lists(_ROW, min_size=1, max_size=25))
def test_property_one_output_row_per_report_row(tmp_path: Path, rows: list[dict[str, Any]]) -> None:
    """F10 and F3: same multiset of precursor keys; ``q_value == float64(Global.Q.Value)``."""
    rows = [{**r, "Modified.Sequence": r["Stripped.Sequence"]} for r in rows]
    res = _convert(tmp_path, rows)
    assert not [i for i in res.issues if i.severity == "error"], res.issues
    df = res.peptides
    assert len(df) == len(rows)

    def key(run: Any, pep: Any, charge: Any, decoy: Any, q: Any, ids: Any) -> tuple[Any, ...]:
        qq = None if _na(q) else float(q)
        return (str(run), str(pep), int(charge), bool(decoy), qq, str(ids))

    want = Counter(
        key(
            r["Run"],
            r["Stripped.Sequence"],
            r["Precursor.Charge"],
            r["Decoy"],
            None if r["Global.Q.Value"] is None else _f32(r["Global.Q.Value"]),
            r["Protein.Ids"],
        )
        for r in rows
    )
    got = Counter(
        key(*t)
        for t in df[["run", "peptide", "charge", "is_decoy", "q_value", "protein_ids"]].itertuples(
            index=False, name=None
        )
    )
    assert got == want
