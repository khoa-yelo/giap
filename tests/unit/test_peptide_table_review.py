"""Lead decisions from the M1-06 code review (generic peptide table; GIAP_SPEC.md §5.3; ADR-0006).

1. rows with an error-level issue are dropped from the evidence (aggregate and ``peptide_runs``);
2. ``nan``/``inf`` text in q_value/pep/intensity is an error (``NA``/empty stay plain missing);
3. duplicate spectrum rows (key uses the stripped peptide) are de-duplicated before summing
   intensity (``duplicate_spectrum_row`` warning);
4. ``protein_ids`` entries are stripped before prefix checks;
5. aggregated ``search_db`` same > other > auto, ``source`` database > denovo > auto;
6. mixed explicit/derived ``is_decoy`` follows the all-rows rule;
7. 200,000 rows read + validate + aggregate in < 10 s.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from giap.config.model import InputConfig
from giap.io.peptides import Issue, aggregate_peptides, read_peptide_table, validate_peptides


def _tsv(*lines: str) -> str:
    return "\n".join(line.replace("|", "\t") for line in lines) + "\n"


def _load(
    tmp_path: Path, text: str, cfg: InputConfig | None = None
) -> tuple[pd.DataFrame, list[Issue]]:
    p = tmp_path / "peptides.tsv"
    p.write_text(text, encoding="utf-8")
    df_raw, issues_read = read_peptide_table(p)
    df, issues_val = validate_peptides(df_raw, cfg or InputConfig())
    return df, [*issues_read, *issues_val]


def _pipeline(
    tmp_path: Path, text: str, cfg: InputConfig | None = None
) -> tuple[pd.DataFrame, pd.DataFrame, list[Issue]]:
    df, issues = _load(tmp_path, text, cfg)
    agg, runs = aggregate_peptides(df)
    return agg, runs, issues


def _codes(issues: list[Issue]) -> set[str]:
    return {i.code for i in issues}


def _row(agg: pd.DataFrame, peptide: str) -> pd.Series:
    sel = agg[agg["peptide"] == peptide]
    assert len(sel) == 1, f"expected one row for {peptide}, got {len(sel)}"
    return sel.iloc[0]


def _error_codes_for(issues: list[Issue], col: str) -> set[str]:
    return {
        i.code
        for i in issues
        if i.severity == "error" and i.code in (f"bad_type:{col}", f"out_of_range:{col}")
    }


# --- 1. error rows are dropped from the evidence ------------------------------------------------

ERROR_CASES = [
    ("q_value", "-0.5", "out_of_range:q_value"),
    ("q_value", "abc", "bad_type:q_value"),
    ("pep", "1.5", "out_of_range:pep"),
    ("intensity", "-3", "out_of_range:intensity"),
]


@pytest.mark.parametrize(("col", "value", "code"), ERROR_CASES)
def test_error_only_row_peptide_absent(tmp_path: Path, col: str, value: str, code: str) -> None:
    text = _tsv(f"peptide|run|{col}", f"PEPTIDEK|r1|{value}", "AAAAKPEP|r1|0.5")
    agg, runs, issues = _pipeline(tmp_path, text)
    iss = [i for i in issues if i.code == code]
    assert iss and iss[0].severity == "error" and iss[0].examples == ("1",)
    assert list(agg["peptide"]) == ["AAAAKPEP"]
    assert list(runs["peptide"]) == ["AAAAKPEP"]


@pytest.mark.parametrize(("col", "value", "code"), ERROR_CASES)
def test_good_and_bad_row_aggregates_good_only(
    tmp_path: Path, col: str, value: str, code: str
) -> None:
    # the bad value lands in ``col`` of the r2 row (better q/pep); the r1 row is fully valid
    cols = ["peptide", "run", "q_value", "pep", "intensity"]
    good = {"peptide": "PEPTIDEK", "run": "r1", "q_value": "0.02", "pep": "0.1", "intensity": "100"}
    bad = {"peptide": "PEPTIDEK", "run": "r2", "q_value": "0.001", "pep": "0.01", "intensity": "7"}
    bad[col] = value
    text = _tsv("|".join(cols), "|".join(good[c] for c in cols), "|".join(bad[c] for c in cols))
    agg, runs, issues = _pipeline(tmp_path, text)
    assert code in _codes(issues)
    r = _row(agg, "PEPTIDEK")
    assert r["n_rows"] == 1
    assert r["n_runs"] == 1
    assert r["q_value"] == pytest.approx(0.02)
    assert r["pep"] == pytest.approx(0.1)
    assert r["intensity"] == pytest.approx(100.0)
    assert list(runs["run"]) == ["r1"]
    assert list(runs["intensity"]) == pytest.approx([100.0])


def test_bad_row_dropped_from_validated_frame_or_aggregate(tmp_path: Path) -> None:
    text = _tsv("peptide|q_value|run", "PEPTIDEK|0.01|r1", "PEPTIDEK|abc|r1")
    agg, runs, _ = _pipeline(tmp_path, text)
    assert _row(agg, "PEPTIDEK")["n_rows"] == 1
    assert list(runs["n_rows"]) == [1]


# --- 2. nan / inf text is an error ---------------------------------------------------------------


@pytest.mark.parametrize("col", ["q_value", "pep", "intensity"])
@pytest.mark.parametrize("value", ["nan", "NaN", "inf", "-inf"])
def test_nan_inf_text_is_error_and_row_dropped(tmp_path: Path, col: str, value: str) -> None:
    text = _tsv(f"peptide|run|{col}", f"PEPTIDEK|r1|{value}", "AAAAKPEP|r1|0.5")
    agg, runs, issues = _pipeline(tmp_path, text)
    codes = _error_codes_for(issues, col)
    assert codes, f"no error for {col}={value!r}: {[i.code for i in issues]}"
    assert list(agg["peptide"]) == ["AAAAKPEP"]
    assert list(runs["peptide"]) == ["AAAAKPEP"]


@pytest.mark.parametrize("col", ["q_value", "pep", "intensity"])
@pytest.mark.parametrize("value", ["NA", ""])
def test_na_and_empty_are_missing_not_error(tmp_path: Path, col: str, value: str) -> None:
    agg, _, issues = _pipeline(tmp_path, _tsv(f"peptide|{col}", f"PEPTIDEK|{value}"))
    assert [i for i in issues if i.severity == "error"] == []
    assert pd.isna(_row(agg, "PEPTIDEK")[col])


# --- 3. duplicate spectrum rows ------------------------------------------------------------------


def test_identical_spectrum_rows_deduplicated(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|run|spectrum_id|charge|intensity",
        "PEPTIDEK|r1|s1|2|100",
        "PEPTIDEK|r1|s1|2|100",
    )
    agg, runs, issues = _pipeline(tmp_path, text)
    iss = [i for i in issues if i.code == "duplicate_spectrum_row"]
    assert iss and iss[0].severity == "warn"
    assert "duplicate_spectrum_id" not in _codes(issues)
    assert _row(agg, "PEPTIDEK")["intensity"] == pytest.approx(100.0)
    assert list(runs["intensity"]) == pytest.approx([100.0])


def test_duplicate_key_uses_stripped_peptide(tmp_path: Path) -> None:
    """Two modified forms of one stripped peptide on one spectrum: one spectrum's evidence.

    Key is (run, spectrum_id, stripped peptide, charge); the first row is kept, intensity counted
    once, warning only. Modified forms on different spectra of a run are still summed (§5.3).
    """
    text = _tsv(
        "peptide|modified_peptide|run|spectrum_id|charge|intensity",
        "PEPTIDEK|PEPT(UniMod:21)IDEK|r1|s1|2|100",
        "PEPTIDEK|PEPTIDEK|r1|s1|2|40",
        "PEPTIDEK|PEPTIDEK|r1|s2|2|30",
    )
    agg, runs, issues = _pipeline(tmp_path, text)
    iss = [i for i in issues if i.code == "duplicate_spectrum_row"]
    assert iss and iss[0].severity == "warn"
    assert "duplicate_spectrum_id" not in _codes(issues)
    assert [i for i in issues if i.severity == "error"] == []
    assert list(runs["intensity"]) == pytest.approx([130.0])  # first s1 row (100) + s2 (30)
    assert _row(agg, "PEPTIDEK")["intensity"] == pytest.approx(130.0)


def test_duplicate_rows_with_modified_form_warn_not_error(tmp_path: Path) -> None:
    """(a) identical rows on the de-dup key: warning only, never ``duplicate_spectrum_id``."""
    text = _tsv(
        "peptide|modified_peptide|run|spectrum_id|charge|intensity",
        "PEPTIDEK|PEPT(UniMod:21)IDEK|r1|s1|2|100",
        "PEPTIDEK|PEPT(UniMod:21)IDEK|r1|s1|2|100",
        "PEPTIDEK|PEPT(UniMod:21)IDEK|r1|s1|2|100",
    )
    agg, runs, issues = _pipeline(tmp_path, text)
    assert "duplicate_spectrum_row" in _codes(issues)
    assert "duplicate_spectrum_id" not in _codes(issues)
    assert [i for i in issues if i.severity == "error"] == []
    assert _row(agg, "PEPTIDEK")["intensity"] == pytest.approx(100.0)
    assert list(runs["intensity"]) == pytest.approx([100.0])


def test_conflicting_spectrum_id_drops_all_rows_of_that_key(tmp_path: Path) -> None:
    """(b) one (run, spectrum_id) with different peptides: error, every row of that key dropped."""
    text = _tsv(
        "peptide|run|spectrum_id|charge|intensity",
        "PEPTIDEK|r1|s1|2|100",
        "AAAAKPEP|r1|s1|2|50",
        "AAAAKPEP|r2|s1|2|7",
        "PEPTIDEK|r1|s2|2|30",
    )
    agg, runs, issues = _pipeline(tmp_path, text)
    iss = [i for i in issues if i.code == "duplicate_spectrum_id"]
    assert iss and iss[0].severity == "error"
    assert set(iss[0].examples) == {"1", "2"}
    assert "duplicate_spectrum_row" not in _codes(issues)
    assert list(agg["peptide"]) == ["AAAAKPEP", "PEPTIDEK"]
    a, p = _row(agg, "AAAAKPEP"), _row(agg, "PEPTIDEK")
    assert (a["n_rows"], a["intensity"]) == (1, pytest.approx(7.0))
    assert (p["n_rows"], p["intensity"]) == (1, pytest.approx(30.0))
    assert list(zip(runs["peptide"], runs["run"], strict=True)) == [
        ("AAAAKPEP", "r2"),
        ("PEPTIDEK", "r1"),
    ]


def test_conflicting_spectrum_id_only_rows_peptides_absent(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|run|spectrum_id|intensity",
        "PEPTIDEK|r1|s1|100",
        "AAAAKPEP|r1|s1|50",
        "LLLLIIIK|r1|s2|5",
    )
    agg, runs, issues = _pipeline(tmp_path, text)
    assert "duplicate_spectrum_id" in _codes(issues)
    assert list(agg["peptide"]) == ["LLLLIIIK"]
    assert list(runs["peptide"]) == ["LLLLIIIK"]


def test_different_charge_or_run_not_duplicates(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|run|spectrum_id|charge|intensity",
        "PEPTIDEK|r1|s1|2|100",
        "PEPTIDEK|r1|s1|3|100",
        "PEPTIDEK|r2|s1|2|100",
    )
    _, runs, issues = _pipeline(tmp_path, text)
    assert "duplicate_spectrum_row" not in _codes(issues)
    assert list(runs["intensity"]) == pytest.approx([200.0, 100.0])


def test_rows_without_spectrum_id_not_deduplicated(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|run|spectrum_id|charge|intensity",
        "PEPTIDEK|r1|NA|2|100",
        "PEPTIDEK|r1|NA|2|100",
        "AAAAKPEP|r1|7|2|5",
    )
    agg, _, issues = _pipeline(tmp_path, text)
    assert "duplicate_spectrum_row" not in _codes(issues)
    assert _row(agg, "PEPTIDEK")["intensity"] == pytest.approx(200.0)


def test_no_spectrum_id_column_not_deduplicated(tmp_path: Path) -> None:
    text = _tsv("peptide|run|charge|intensity", "PEPTIDEK|r1|2|100", "PEPTIDEK|r1|2|100")
    agg, _, issues = _pipeline(tmp_path, text)
    assert "duplicate_spectrum_row" not in _codes(issues)
    assert _row(agg, "PEPTIDEK")["intensity"] == pytest.approx(200.0)


# --- 4. protein_ids entries stripped -------------------------------------------------------------


def test_protein_ids_entries_stripped(tmp_path: Path) -> None:
    agg, _, _ = _pipeline(tmp_path, _tsv("peptide|protein_ids", "PEPTIDEK|P1; rev_P2"))
    r = _row(agg, "PEPTIDEK")
    assert list(r["protein_ids"]) == ["P1", "rev_P2"]
    assert not bool(r["is_decoy"])


def test_spaced_decoy_entries_recognised(tmp_path: Path) -> None:
    df, _ = _load(tmp_path, _tsv("peptide|protein_ids", "PEPTIDEK|rev_P1; rev_P2"))
    assert bool(df["is_decoy"].iloc[0])
    agg, _ = aggregate_peptides(df)
    r = _row(agg, "PEPTIDEK")
    assert bool(r["is_decoy"])
    assert list(r["protein_ids"]) == ["rev_P1", "rev_P2"]


def test_spaced_contaminant_entries_recognised(tmp_path: Path) -> None:
    df, _ = _load(tmp_path, _tsv("peptide|protein_ids", "PEPTIDEK|CON__A ;  CON__B"))
    assert bool(df["is_contaminant"].iloc[0])


# --- 5. search_db / source precedence ------------------------------------------------------------


@pytest.mark.parametrize(
    ("values", "expected"),
    [(("other", "same"), "same"), (("same", "other"), "same"), (("other", "NA"), "other")],
)
def test_search_db_precedence(tmp_path: Path, values: tuple[str, str], expected: str) -> None:
    text = _tsv("peptide|search_db", *(f"PEPTIDEK|{v}" for v in values))
    agg, _, _ = _pipeline(tmp_path, text, InputConfig(search_db="auto"))
    assert _row(agg, "PEPTIDEK")["search_db"] == expected


@pytest.mark.parametrize(
    ("values", "expected"),
    [(("denovo", "database"), "database"), (("NA", "denovo"), "denovo"), (("NA",), "auto")],
)
def test_source_precedence(tmp_path: Path, values: tuple[str, ...], expected: str) -> None:
    text = _tsv("peptide|source", *(f"PEPTIDEK|{v}" for v in values))
    agg, _, _ = _pipeline(tmp_path, text, InputConfig(peptide_source="auto"))
    assert _row(agg, "PEPTIDEK")["source"] == expected


# --- 6. mixed explicit / derived is_decoy --------------------------------------------------------


def test_mixed_explicit_false_and_prefix_decoy_is_not_decoy(tmp_path: Path) -> None:
    text = _tsv("peptide|protein_ids|is_decoy", "PEPTIDEK|P1|false", "PEPTIDEK|rev_P1;rev_P2|")
    df, _ = _load(tmp_path, text)
    assert [bool(x) for x in df["is_decoy"]] == [False, True]
    agg, _ = aggregate_peptides(df)
    assert not bool(_row(agg, "PEPTIDEK")["is_decoy"])


# --- 7. performance ------------------------------------------------------------------------------

_AA = np.array(list("ACDEFGHKLMNPQRSTVWY"))


def _synthetic_table(path: Path, n_rows: int = 200_000, n_peptides: int = 40_000) -> None:
    rng = np.random.default_rng(0)
    # unique 9-mers: base-19 digits of the peptide index (4 positions cover 130,321), then fixed
    codes = np.arange(n_peptides)
    digits = np.stack([(codes // 19**k) % 19 for k in range(4)], axis=1)
    stems = ["".join(r) for r in _AA[digits]]
    peptides = np.array([f"PEP{s}K" for s in stems])
    idx = rng.integers(0, n_peptides, n_rows)
    run = rng.integers(1, 4, n_rows)
    q = rng.random(n_rows)
    pep = rng.random(n_rows)
    inten = rng.lognormal(10, 2, n_rows)
    charge = rng.integers(1, 5, n_rows)
    prot = np.where(rng.random(n_rows) < 0.1, "rev_P", "P")
    lines = ["peptide\tprotein_ids\tq_value\tpep\tintensity\trun\tcharge\tspectrum_id"]
    lines += [
        f"{peptides[i]}\t{pr}{i}\t{qq:.6g}\t{pp:.6g}\t{ii:.6g}\tr{rr}\t{cc}\ts{k}"
        for k, (i, pr, qq, pp, ii, rr, cc) in enumerate(
            zip(idx, prot, q, pep, inten, run, charge, strict=True)
        )
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_200k_rows_under_10s(tmp_path: Path) -> None:
    path = tmp_path / "big.tsv"
    _synthetic_table(path)
    t0 = time.perf_counter()
    df_raw, _ = read_peptide_table(path)
    df, issues = validate_peptides(df_raw, InputConfig())
    agg, runs = aggregate_peptides(df)
    elapsed = time.perf_counter() - t0
    assert [i for i in issues if i.severity == "error"] == []
    assert len(df) == 200_000
    assert 39_000 <= len(agg) <= 40_000
    assert int(runs["n_rows"].sum()) == 200_000
    # loose 10x guard; the 10 s budget is enforced in benchmarks/test_peptide_table.py
    assert elapsed < 100.0, f"read+validate+aggregate took {elapsed:.1f} s"
