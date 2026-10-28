"""Lead decisions from the M1-07 code review (GIPPA_SPEC.md §5.3, §5.4).

1. Stripping removes ()/[]/{}/<> groups repeatedly, so nested groups go too.
2. Flanking residues (``K.PEPTIDE.R``) are removed before stripping.
3. DIA-NN MaxLFQ abundance excludes decoys and protein groups with no positive value.
4. Gaps found by mutation testing (decoy codes, flag filling, sort order, MaxLFQ details,
   generic validation, reproducible gzip).
5. Performance: a 500k-row DIA-NN parquet converts in < 15 s.
"""

from __future__ import annotations

import gzip
import math
import time
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from gippa.cli import cli
from gippa.config.model import InputConfig
from gippa.exitcodes import ExitCode
from gippa.io.converters.base import ConvertResult
from gippa.io.converters.diann import convert_diann
from gippa.io.converters.generic import convert_generic
from gippa.io.converters.strip import strip_modifications
from gippa.io.issues import Issue

ROW: dict[str, object] = {
    "Run": "r1",
    "Protein.Group": "P1",
    "Protein.Ids": "P1",
    "Stripped.Sequence": "PEPTIDEK",
    "Modified.Sequence": "PEPTIDEK",
    "Precursor.Charge": 2,
    "Q.Value": 0.002,
    "Global.Q.Value": 0.004,
    "PEP": 0.01,
    "Precursor.Quantity": 1400.0,
    "Precursor.Normalised": 1500.0,
    "PG.MaxLFQ": 10.0,
    "Decoy": 0,
}


def _tsv(tmp_path: Path, rows: Sequence[dict[str, object]], name: str = "report.tsv") -> Path:
    p = tmp_path / name
    pd.DataFrame(list(rows)).to_csv(p, sep="\t", index=False, lineterminator="\n")
    return p


def _parquet(tmp_path: Path, rows: Sequence[dict[str, object]]) -> Path:
    p = tmp_path / "report.parquet"
    pd.DataFrame(list(rows)).to_parquet(p, engine="pyarrow", index=False)
    return p


def _diann(paths: Sequence[Path]) -> ConvertResult:
    return convert_diann(list(paths), search_db=None, cfg=InputConfig())


def _generic(
    tmp_path: Path, text: str, source: str | None = None, search_db: str | None = None
) -> ConvertResult:
    p = tmp_path / "in.tsv"
    p.write_text(text, encoding="utf-8")
    return convert_generic([p], source=source, search_db=search_db, cfg=InputConfig())


def _errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


def _abundance(res: ConvertResult) -> dict[str, float]:
    assert res.abundance is not None
    return dict(zip(res.abundance["protein_id"], res.abundance["maxlfq"], strict=True))


# --- 1. repeated / nested group stripping -------------------------------------------------


@pytest.mark.parametrize(
    ("seq", "expected"),
    [
        ("_(Acetyl (Protein N-term))PEPM(Oxidation (M))IDE_", "PEPMIDE"),
        ("PEPS[Phospho[UNIMOD:21]]TIDE", "PEPSTIDE"),
        ("{Hex}PEPTIDE", "PEPTIDE"),
        ("<13C>PEPTIDE", "PEPTIDE"),
        ("PEP{a{b}c}TIDE", "PEPTIDE"),
        ("PEP<x<y>>TIDE", "PEPTIDE"),
        ("PEP[A(B{C<D>})]TIDE", "PEPTIDE"),
    ],
)
def test_nested_groups_removed(seq: str, expected: str) -> None:
    assert strip_modifications(seq) == expected


# --- 2. flanking residues -----------------------------------------------------------------


@pytest.mark.parametrize(
    ("seq", "expected"),
    [
        ("K.PEPTIDE.R", "PEPTIDE"),
        ("-.PEPTIDE.A", "PEPTIDE"),
        ("K.PEPTIDE.-", "PEPTIDE"),
        ("K.PEP(UniMod:35)TIDE.R", "PEPTIDE"),
        ("R.M[+15.995]PEPTIDE.G", "MPEPTIDE"),
    ],
)
def test_flanking_residues_removed(seq: str, expected: str) -> None:
    assert strip_modifications(seq) == expected


@pytest.mark.parametrize("seq", ["KPEPTIDER", "APEPTIDEA", "PEPTIDE"])
def test_no_dots_unaffected(seq: str) -> None:
    assert strip_modifications(seq) == seq


def test_generic_converter_strips_flanks(tmp_path: Path) -> None:
    res = _generic(tmp_path, "peptide\nK.PEPTIDEK.R\n")
    assert res.peptides["peptide"].tolist() == ["PEPTIDEK"]


# --- 3. MaxLFQ: decoys and all-0/NA groups ------------------------------------------------


def test_maxlfq_excludes_decoy_rows(tmp_path: Path) -> None:
    rows = [
        {**ROW, "Run": "r1", "PG.MaxLFQ": 10.0},
        {**ROW, "Run": "r2", "PG.MaxLFQ": 30.0},
        {**ROW, "Run": "r1", "Protein.Group": "REV_P9", "Protein.Ids": "REV_P9",
         "Stripped.Sequence": "KEDITPEP", "Modified.Sequence": "KEDITPEP", "Decoy": 1,
         "PG.MaxLFQ": 99.0},
    ]  # fmt: skip
    for path in (_tsv(tmp_path, rows), _parquet(tmp_path, rows)):
        vals = _abundance(_diann([path]))
        assert set(vals) == {"P1"}
        assert vals["P1"] == pytest.approx(20)


def test_maxlfq_decoy_row_does_not_affect_target_group(tmp_path: Path) -> None:
    # A decoy row sharing the group name (in a run listed first) must not supply that run's value.
    rows = [
        {**ROW, "Run": "r1", "Decoy": 1, "PG.MaxLFQ": 1000.0},
        {**ROW, "Run": "r1", "PG.MaxLFQ": 10.0},
        {**ROW, "Run": "r2", "PG.MaxLFQ": 30.0},
    ]
    vals = _abundance(_diann([_tsv(tmp_path, rows)]))
    assert vals == {"P1": pytest.approx(20)}


def test_maxlfq_drops_groups_with_only_zero_or_na(tmp_path: Path) -> None:
    rows = [
        {**ROW, "Run": "r1", "PG.MaxLFQ": 10.0},
        {**ROW, "Run": "r1", "Protein.Group": "Z0", "Protein.Ids": "Z0",
         "Stripped.Sequence": "AAAK", "Modified.Sequence": "AAAK", "PG.MaxLFQ": 0.0},
        {**ROW, "Run": "r2", "Protein.Group": "Z0", "Protein.Ids": "Z0",
         "Stripped.Sequence": "AAAK", "Modified.Sequence": "AAAK", "PG.MaxLFQ": float("nan")},
        {**ROW, "Run": "r1", "Protein.Group": "ZN", "Protein.Ids": "ZN",
         "Stripped.Sequence": "CCCK", "Modified.Sequence": "CCCK", "PG.MaxLFQ": float("nan")},
    ]  # fmt: skip
    res = _diann([_tsv(tmp_path, rows)])
    assert res.abundance is not None
    assert res.abundance["protein_id"].tolist() == ["P1"]
    assert not res.abundance["maxlfq"].isna().any()


# --- 4. mutation gaps ---------------------------------------------------------------------


@pytest.mark.parametrize("decoy", [2, -1, 0.5, "yes"])
def test_decoy_other_than_0_1_is_error(tmp_path: Path, decoy: object) -> None:
    rows = [ROW, {**ROW, "Run": "r2", "Decoy": decoy}]
    res = _diann([_tsv(tmp_path, rows)])
    errs = _errors(res.issues)
    assert errs, [i.code for i in res.issues]
    assert any("Decoy" in i.code for i in errs)


def test_source_flag_fills_only_empty_cells(tmp_path: Path) -> None:
    text = "peptide\trun\tsource\tsearch_db\nPEPTIDEK\tr1\tdatabase\tsame\nAAAAAAAK\tr2\t\t\n"
    res = _generic(tmp_path, text, source="denovo", search_db="other")
    assert not _errors(res.issues)
    df = res.peptides.set_index("peptide")
    assert df.loc["PEPTIDEK", "source"] == "database"
    assert df.loc["PEPTIDEK", "search_db"] == "same"
    assert df.loc["AAAAAAAK", "source"] == "denovo"
    assert df.loc["AAAAAAAK", "search_db"] == "other"


def test_sort_na_last_and_ties_stable(tmp_path: Path) -> None:
    text = (
        "peptide\trun\tcharge\tq_value\n"
        "PEPTIDEK\tr2\t2\t0.01\n"
        "AAAAAAAK\tNA\t2\t0.02\n"
        "PEPTIDEK\tr1\tNA\t0.03\n"
        "PEPTIDEK\tr1\t2\t0.04\n"
        "PEPTIDEK\tr1\t2\t0.05\n"
        "AAAAAAAK\tr1\t2\t0.06\n"
        "PEPTIDEK\tr1\t2\t0.07\n"
    )
    res = _generic(tmp_path, text)
    assert not _errors(res.issues)
    assert res.peptides["q_value"].tolist() == pytest.approx(
        [0.06, 0.04, 0.05, 0.07, 0.03, 0.01, 0.02]
    )
    runs = res.peptides["run"].tolist()
    assert runs[:-1] == ["r1", "r1", "r1", "r1", "r1", "r2"]
    assert runs[-1] is pd.NA or runs[-1] is None


def test_diann_ties_keep_input_order(tmp_path: Path) -> None:
    rows = [
        {**ROW, "Run": "r2", "PEP": 0.1},
        {**ROW, "Run": "r1", "PEP": 0.2},
        {**ROW, "Run": "r1", "PEP": 0.3},
        {**ROW, "Run": "r1", "PEP": 0.4},
    ]
    res = _diann([_tsv(tmp_path, rows)])
    assert res.peptides["pep"].tolist() == pytest.approx([0.2, 0.3, 0.4, 0.1])


def test_maxlfq_excludes_na_protein_group(tmp_path: Path) -> None:
    rows = [
        {**ROW, "Run": "r1", "PG.MaxLFQ": 10.0},
        {**ROW, "Run": "r1", "Protein.Group": None, "Stripped.Sequence": "AAAK",
         "Modified.Sequence": "AAAK", "PG.MaxLFQ": 55.0},
    ]  # fmt: skip
    for path in (_tsv(tmp_path, rows), _parquet(tmp_path, rows)):
        res = _diann([path])
        assert res.abundance is not None
        assert res.abundance["protein_id"].tolist() == ["P1"]
        assert res.abundance["maxlfq"].tolist() == pytest.approx([10.0])


def test_maxlfq_conflicting_duplicate_uses_first_not_mean(tmp_path: Path) -> None:
    rows = [
        {**ROW, "Run": "r1", "PG.MaxLFQ": 10.0},
        {**ROW, "Run": "r1", "Stripped.Sequence": "AAAK", "Modified.Sequence": "AAAK",
         "PG.MaxLFQ": 50.0},
        {**ROW, "Run": "r2", "PG.MaxLFQ": 30.0},
    ]  # fmt: skip
    vals = _abundance(_diann([_tsv(tmp_path, rows)]))
    # first() → r1=10, median(10, 30) = 20; averaging r1 would give 30 (or 30 over all rows).
    assert vals["P1"] == pytest.approx(20)


@pytest.mark.parametrize("q", ["1.5", "-0.1"])
def test_generic_reports_validation_issues(tmp_path: Path, q: str) -> None:
    res = _generic(tmp_path, f"peptide\tq_value\nPEPTIDEK\t{q}\n")
    assert "out_of_range:q_value" in {i.code for i in _errors(res.issues)}


def test_generic_validation_error_exits_1_via_cli(tmp_path: Path) -> None:
    src = tmp_path / "in.tsv"
    src.write_text("peptide\tq_value\nPEPTIDEK\t1.5\n", encoding="utf-8")
    out = tmp_path / "peptides.tsv"
    result = CliRunner().invoke(cli, ["convert", "--from", "generic", str(src), "--out", str(out)])
    assert result.exit_code == ExitCode.INVALID_INPUT
    assert not out.exists()
    assert "out_of_range:q_value" in result.stderr


def test_gz_output_is_reproducible(tmp_path: Path) -> None:
    src = _tsv(tmp_path, [ROW, {**ROW, "Run": "r2"}])
    blobs = []
    for name in ("a.tsv.gz", "other_name.tsv.gz"):
        out = tmp_path / name
        argv = ["convert", "--from", "diann", str(src), "--out", str(out)]
        result = CliRunner().invoke(cli, argv)
        assert result.exit_code == ExitCode.SUCCESS, result.output
        blob = out.read_bytes()
        assert blob[:2] == b"\x1f\x8b"
        assert blob[3] & 0x08 == 0  # no FNAME field
        assert blob[4:8] == b"\x00\x00\x00\x00"  # mtime 0
        assert gzip.decompress(blob).decode().startswith("peptide\t")
        blobs.append(blob)
        time.sleep(1.1)
    assert blobs[0] == blobs[1]


# --- 5. performance -----------------------------------------------------------------------


def _big_report(path: Path, n: int = 500_000, n_filler: int = 50) -> None:
    rng = np.random.default_rng(0)
    peps = np.array(["PEPTIDEK", "AAAAAAAK", "LLGNVLVCVLAHHFGK", "VVDLMAHMASK", "SAMPLER"])
    mods = np.array(["PEPTIDEK", "AAAAAAAK", "LLGNVLVC(UniMod:4)VLAHHFGK", "VVDLM(UniMod:35)AHMASK",
                     "SAMPLER"])  # fmt: skip
    pick = rng.integers(0, len(peps), n)
    groups = np.array([f"P{i:05d}" for i in range(5000)])
    gpick = rng.integers(0, len(groups), n)
    runs = np.array([f"run{i:02d}" for i in range(20)])
    cols: dict[str, pa.Array] = {
        "Run": pa.array(runs[rng.integers(0, len(runs), n)]),
        "Protein.Group": pa.array(groups[gpick]),
        "Protein.Ids": pa.array(groups[gpick]),
        "Stripped.Sequence": pa.array(peps[pick]),
        "Modified.Sequence": pa.array(mods[pick]),
        "Precursor.Charge": pa.array(rng.integers(1, 5, n).astype(np.int64)),
        "Q.Value": pa.array(rng.uniform(0, 0.01, n)),
        "Global.Q.Value": pa.array(rng.uniform(0, 0.01, n)),
        "PEP": pa.array(rng.uniform(0, 0.1, n)),
        "Precursor.Quantity": pa.array(rng.uniform(0, 1e6, n)),
        "Precursor.Normalised": pa.array(rng.uniform(0, 1e6, n)),
        "PG.MaxLFQ": pa.array(rng.uniform(0, 1e6, n)),
        "Decoy": pa.array((rng.uniform(size=n) < 0.05).astype(np.int64)),
    }
    for i in range(n_filler):
        if i % 2:
            cols[f"Filler.Num{i}"] = pa.array(rng.normal(size=n))
        else:
            cols[f"Filler.Str{i}"] = pa.array(runs[rng.integers(0, len(runs), n)])
    pq.write_table(pa.table(cols), path)


def test_large_parquet_converts_quickly(tmp_path: Path) -> None:
    path = tmp_path / "report.parquet"
    _big_report(path)
    assert pq.read_metadata(path).num_columns >= 60
    start = time.perf_counter()
    res = _diann([path])
    elapsed = time.perf_counter() - start
    assert not _errors(res.issues), [i.code for i in res.issues]
    assert len(res.peptides) == 500_000
    assert res.abundance is not None and len(res.abundance) > 0
    assert not any(math.isnan(v) for v in res.abundance["maxlfq"])
    # loose 10x guard; the 15 s budget is enforced in benchmarks/test_converters.py
    assert elapsed < 150, f"conversion took {elapsed:.1f} s"
