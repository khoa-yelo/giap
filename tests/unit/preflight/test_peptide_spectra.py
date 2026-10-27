"""M7-10 lead decision 1: preflight writes ``work/peptide_spectra.parquet`` (C17) for the de novo
and reference rows, an optional contract table the ``variants`` step reads
(docs/specs/M7-10-class-error-rates.md "Lead decisions (2026-10-22)" 1; GIPPA_SPEC.md §5.3
``--reference-peptides``, G2, §7.5.4).

Interface pinned here (test-author reading; see the report to the lead):

* ``gippa.preflight.run.PreflightInputs`` gains the keyword field ``reference_peptides: Path |
  None = None`` (the ``--reference-peptides`` file; ``gippa run`` already has the option).
* With de novo rows that carry a ``spectrum_id`` and a reference table that does too,
  ``write_outputs`` writes ``OUT/work/preflight/peptide_spectra.parquet`` (next to
  ``peptides.norm.parquet`` and ``peptide_runs.parquet``). Every de novo row of ``--peptides``
  is one row with ``source = denovo`` and its stripped ``peptide``; every confident
  (``q_value <= 0.01``) target reference row is one row whose ``source`` is not ``denovo``.
  Database rows of ``--peptides`` are *not* reference rows (§7.5.4 pairs with
  ``--reference-peptides`` only). Decoy and low-confidence (q > 0.01) reference rows are dropped
  here: preflight is the only decoy filter on the pipeline path (M7-10 lead decision 5).
* C17 is registered in ``gippa.genetable.contract.TABLES["peptide_spectra"]`` with the required
  columns ``run``, ``spectrum_id``, ``peptide`` and ``source`` (more columns allowed), in
  ``gippa.io.work_tables.EXTRA_FILE_NAMES`` as ``peptide_spectra.parquet`` (like C16), and
  preflight's own output validates without errors.
* ``gippa run --executor local`` (``cli._local_roots``) stages it as the root key
  ``peptide_spectra`` when ``--reference-peptides`` is given.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from gippa.config import resolve_config
from gippa.exitcodes import ExitCode

LOW = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
HEADER = ["peptide", "protein_ids", "q_value", "score", "is_decoy", "run", "charge",
          "spectrum_id", "source"]  # fmt: skip
# de novo rows of --peptides: (peptide, run, spectrum_id)
DENOVO = [("AAASAAKR", "run1", "run1.1"), ("GGCPEPTK", "run1", "run1.2"),
          ("PEPTIDEKR", "run2", "run2.1"), ("MSLNWVAK", "run2", "run2.7")]  # fmt: skip
# a database row of --peptides (not a reference row)
DATABASE = [("QQRSTVWK", "run1", "run1.9", "0.001")]
# reference rows: (peptide, run, spectrum_id, q_value, is_decoy)
REFERENCE = [("AAAAAAKR", "run1", "run1.1", "0.001", "false"),
             ("NCPEPTKA", "run1", "run1.2", "0.004", "false"),
             ("PEPTLDEKR", "run2", "run2.1", "0.02", "false"),
             ("KAVWNLSM", "run2", "run2.7", "0.001", "true")]  # fmt: skip


@pytest.fixture(autouse=True)
def _isolated_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for name in ("GIPPA_DB_DIR", "GIPPA_REPO", "GIPPA_ROOT"):
        d = tmp_path / f"env_{name.lower()}"
        d.mkdir()
        monkeypatch.setenv(name, str(d))


def _tsv(path: Path, rows: list[list[str]]) -> Path:
    path.write_text("\n".join("\t".join(r) for r in [HEADER, *rows]) + "\n", encoding="utf-8")
    return path


def _genome(root: Path) -> Path:
    rng = np.random.default_rng(5)
    seq = "".join(np.array(list("ACGT"))[rng.integers(0, 4, 20_000)])
    g = root / "genome.fa"
    g.write_text(f">chr1 circular=true\n{seq}\n", encoding="utf-8")
    return g


def _peptides(root: Path, with_sid: bool = True) -> Path:
    rows = [[p, "NA", "NA", "1.5", "false", r, "2", s if with_sid else "NA", "denovo"]
            for p, r, s in DENOVO]  # fmt: skip
    rows += [[p, "P1", q, "3.0", "false", r, "2", s, "database"] for p, r, s, q in DATABASE]
    return _tsv(root / "peptides.tsv", rows)


def _reference(root: Path, with_sid: bool = True) -> Path:
    rows = [[p, "P1", q, "3.0", d, r, "2", s if with_sid else "NA", "database"]
            for p, r, s, q, d in REFERENCE]  # fmt: skip
    return _tsv(root / "reference.tsv", rows)


def _fake_probe() -> Any:
    from gippa.preflight.environment import GIB, EnvProbe, requirements

    found = {r.name: r.min_version or "1.0"
             for r in requirements(["core", "mobile", "bgc", "function"], "genome")}  # fmt: skip
    return EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)


def run_preflight_with_reference(root: Path, *, ref_sid: bool = True) -> tuple[Any, Path]:
    from gippa.preflight.run import PreflightInputs, exit_code, run_preflight, write_outputs

    root.mkdir(parents=True, exist_ok=True)
    inp = PreflightInputs(_genome(root), None, _peptides(root), None, None,
                          reference_peptides=_reference(root, with_sid=ref_sid))  # fmt: skip
    out = root / "out"
    report, norm = run_preflight(inp, resolve_config(overrides=LOW), _fake_probe(), "preflight",
                                 out)  # fmt: skip
    assert exit_code(report, False) == ExitCode.SUCCESS, [
        (c.rule_id, c.message) for c in report.checks if c.level == "error"
    ]
    write_outputs(out, report, norm)
    return report, out / "work" / "preflight"


def test_preflight_inputs_has_reference_peptides_field() -> None:
    import dataclasses

    from gippa.preflight.run import PreflightInputs

    fields = {f.name: f for f in dataclasses.fields(PreflightInputs)}
    assert "reference_peptides" in fields
    assert fields["reference_peptides"].default is None


def test_preflight_writes_peptide_spectra(tmp_path: Path) -> None:
    _report, work = run_preflight_with_reference(tmp_path)
    path = work / "peptide_spectra.parquet"
    assert path.is_file(), sorted(p.name for p in work.iterdir())
    df = pd.read_parquet(path)
    assert {"run", "spectrum_id", "peptide", "source"} <= set(df.columns)
    rows = {(str(r.run), str(r.spectrum_id), str(r.peptide), str(r.source))
            for r in df.itertuples(index=False)}  # fmt: skip
    for p, r, s in DENOVO:
        assert (r, s, p, "denovo") in rows, (r, s, p)
    assert sum(src == "denovo" for *_x, src in rows) == len(DENOVO)
    ref = {(r, s, p) for r, s, p, src in rows if src != "denovo"}
    # the two confident target reference rows are there
    assert ("run1", "run1.1", "AAAAAAKR") in ref
    assert ("run1", "run1.2", "NCPEPTKA") in ref
    # the database row of --peptides is not a reference row
    assert not any(s == "run1.9" for _r, s, _p in ref)
    # decision 5: the confident decoy and the q 0.02 target are dropped (review #239 B2)
    assert ("run2", "run2.7", "KAVWNLSM") not in ref
    assert ("run2", "run2.1", "PEPTLDEKR") not in ref
    assert len(ref) == 2, ref


def test_peptide_spectra_is_contract_table_c17(tmp_path: Path) -> None:
    from gippa.genetable.contract import TABLES, validate_table
    from gippa.io.work_tables import EXTRA_FILE_NAMES

    assert "peptide_spectra" in TABLES
    required = {c.name for c in TABLES["peptide_spectra"] if not c.optional}
    assert {"run", "spectrum_id", "peptide", "source"} <= required
    assert EXTRA_FILE_NAMES.get("peptide_spectra") == "peptide_spectra.parquet"
    _report, work = run_preflight_with_reference(tmp_path)
    df = pd.read_parquet(work / "peptide_spectra.parquet")
    errors = [i for i in validate_table("peptide_spectra", df) if i.severity == "error"]
    assert errors == []
    for col in ("run", "spectrum_id", "peptide", "source"):
        bad = validate_table("peptide_spectra", df.drop(columns=[col]))
        assert any(i.code == "contract.peptide_spectra.missing_column" for i in bad), col


def test_load_tables_reads_c17_and_a_missing_file_is_a_warning(tmp_path: Path) -> None:
    from gippa.io.work_tables import load_tables

    _report, work = run_preflight_with_reference(tmp_path / "a")
    frames, issues = load_tables({"peptide_spectra": work / "peptide_spectra.parquet"},
                                 ["peptide_spectra"])  # fmt: skip
    assert [i for i in issues if i.severity == "error"] == []
    assert len(frames["peptide_spectra"]) >= len(DENOVO)
    frames, issues = load_tables({"peptide_spectra": tmp_path / "absent.parquet"},
                                 ["peptide_spectra"])  # fmt: skip
    assert [i for i in issues if i.severity == "error"] == []
    assert any(i.severity == "warn" for i in issues)


def test_reference_without_spectrum_id_gives_no_reference_rows(tmp_path: Path) -> None:
    """Degenerate row "reference rows without spectrum_id": nothing to pair (the warning
    ``no_spectrum_pairing`` is pinned at step level, tests/unit/steps/
    test_variants_step_class_error.py)."""
    _report, work = run_preflight_with_reference(tmp_path, ref_sid=False)
    path = work / "peptide_spectra.parquet"
    if path.is_file():
        df = pd.read_parquet(path)
        ref = df[df["source"].astype(str) != "denovo"]
        assert ref["spectrum_id"].isna().all() or len(ref) == 0


def test_local_executor_stages_peptide_spectra(tmp_path: Path) -> None:
    from helpers.mapping_steps import write_test_db

    from gippa.cli import _local_roots
    from gippa.config import GippaConfig

    db = write_test_db(tmp_path / "db")
    root = tmp_path / "in"
    root.mkdir()
    inputs: dict[str, Path | None] = {
        "proteins": None, "peptides": _peptides(root), "abundance": None, "annotation": None,
        "reference_peptides": _reference(root), "db_dir": db,
    }  # fmt: skip
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIPPA_DB_DIR", str(db))
        cfg = GippaConfig.model_validate(
            {"preflight": {"min_n50_warn": 1, "min_largest_contig_warn": 1}}
        )
        roots = _local_roots(cfg, tmp_path / "OUT", _genome(root), inputs)
    assert "peptide_spectra" in roots, sorted(roots)
    assert Path(roots["peptide_spectra"]).is_file()
    assert Path(roots["peptide_spectra"]).name == "peptide_spectra.parquet"
