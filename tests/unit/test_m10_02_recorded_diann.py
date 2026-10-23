"""The DIA-NN pin and the licence/size guard on the committed tree (M10-02 criteria 9-10).

Spec: docs/specs/M10-02-diann-searches.md deliverables 1 and 6, "The committed samples",
"Licence and size rules"; ADR-0051 "Image handling"; R5, R6, R11.

Khoa's decision of 2026-10-23 replaces the committed samples: **no DIA-NN report rows are ever
committed**. The stratified sample stays on FarmShare under
``$GIAP_ROOT/pilot/recorded/diann/2.3.1/<arm>/``. Under ``tests/data/recorded/diann/`` only these
may be committed:

* ``2.3.1/PROVENANCE.md``;
* ``2.3.1/<arm>/schema.tsv``: the column list of the full report (``column``, ``arrow_type``,
  ``nullable``; one line per column, no data). This is the test author's reading of "the column
  list of the full report"; it is checked to hold column metadata only.

The guard tests apply to whatever is tracked and pass vacuously until something lands; the
pin test fails until deliverable 1 exists. Tests that need real sample rows read ``$GIAP_ROOT``
and skip with a reason when the sample is absent (as on CI).
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest
from helpers.pilot_kit import DIANN_SIF_SHA256
from helpers.pilot_kit_m10 import PIN_FILE, PIN_NAME, README, RECORDED_DIANN, REPO

R11_MAX_BYTES = 5 * 1024 * 1024
TESTS_DATA_MAX = 50 * 1024 * 1024
ALLOWED = re.compile(r"^2\.3\.1/(PROVENANCE\.md|A[12]/schema\.tsv)$")
ROWS_SUFFIXES = {".parquet", ".tsv", ".csv", ".gz", ".zst", ".feather", ".arrow", ".json"}
SCHEMA_HEADER = ["column", "arrow_type", "nullable"]
ARROW_TYPE = re.compile(
    r"^(null|bool|u?int(8|16|32|64)|halffloat|float|double|string|large_string|binary|"
    r"large_binary|date32\[day\]|date64\[ms\]|timestamp\[\w+(, tz=[^\]]+)?\]|"
    r"(large_)?list<.*>|dictionary<.*>|struct<.*>|decimal\d*\(\d+, \d+\))$"
)


def _tracked(*paths: str) -> list[str]:
    res = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z", "--", *paths],
                         capture_output=True, text=True, check=True)  # fmt: skip
    return [p for p in res.stdout.split("\0") if p]


# ------------------------------------------------------------ criterion 9: one pin


def _pin_file_hash() -> str:
    assert PIN_FILE.is_file(), f"{PIN_FILE.relative_to(REPO)} does not exist (M10-02 deliverable 1)"
    lines = PIN_FILE.read_text().splitlines()
    assert len(lines) == 1, f"the pin file holds exactly one line: {lines}"
    m = re.fullmatch(r"([0-9a-f]{64})  (\S+)", lines[0])
    assert m, f"not sha256sum format '<64 hex>  {PIN_NAME}': {lines[0]!r}"
    assert m.group(2) == PIN_NAME, m.group(2)
    return m.group(1)


def test_the_pin_file_equals_the_test_constant() -> None:
    assert _pin_file_hash() == DIANN_SIF_SHA256


def test_the_pin_file_equals_the_readme() -> None:
    pin = _pin_file_hash()
    pins = set(re.findall(r"\b[0-9a-f]{64}\b", README.read_text()))
    assert pin in pins, "the README's image sha256 is not the pin file's"


def test_the_readme_pin_line_names_the_pin_file() -> None:
    lines = [ln for ln in README.read_text().splitlines() if DIANN_SIF_SHA256 in ln]
    assert any("diann_image.sha256" in ln for ln in lines), (
        "deliverable 7: the README's pin line names validation/pilot/diann_image.sha256"
    )


def test_the_pin_file_is_tracked() -> None:
    assert _tracked(str(PIN_FILE.relative_to(REPO))), "the pin file must be committed"


# ------------------------------------------------------------ criterion 10 (Khoa 2026-10-23)


def test_only_provenance_and_the_column_list_are_committed_under_recorded_diann() -> None:
    rel = RECORDED_DIANN.relative_to(REPO)
    bad = [p for p in _tracked(str(rel))
           if not ALLOWED.match(str(Path(p).relative_to(rel)))]  # fmt: skip
    assert not bad, f"not allowed under {rel} (Khoa 2026-10-23: no DIA-NN rows committed): {bad}"


def test_no_dia_nn_derived_table_with_rows_is_committed() -> None:
    rel = RECORDED_DIANN.relative_to(REPO)
    for p in _tracked(str(rel)):
        path = REPO / p
        if path.name == "schema.tsv":
            continue  # the column list, checked below
        assert path.suffix.lower() not in ROWS_SUFFIXES, f"DIA-NN table committed: {p}"
        assert not path.read_bytes().startswith(b"PAR1"), f"parquet committed: {p}"


def test_a_committed_column_list_holds_column_metadata_only() -> None:
    rel = RECORDED_DIANN.relative_to(REPO)
    for p in _tracked(str(rel)):
        if Path(p).name != "schema.tsv":
            continue
        rows = [ln.split("\t") for ln in (REPO / p).read_text().splitlines()]
        assert rows and rows[0] == SCHEMA_HEADER, f"{p}: header {rows[:1]}"
        names = [r[0] for r in rows[1:]]
        assert len(names) == len(set(names)), f"{p}: repeated column names"
        for r in rows[1:]:
            assert len(r) == 3, f"{p}: {r}"
            assert ARROW_TYPE.match(r[1]), f"{p}: {r[1]!r} is not an Arrow type"
            assert r[2].lower() in {"true", "false"}, f"{p}: {r}"
        assert "Run" in names and "Decoy" in names, f"{p}: not a DIA-NN report's columns"


def test_no_container_library_or_model_or_elf_file_is_tracked() -> None:
    for p in _tracked("tests", "validation"):
        path = REPO / p
        assert path.suffix.lower() not in {".sif", ".speclib", ".pt", ".so", ".dll"}, p
        assert not p.endswith(("_empirical.parquet", "_report.parquet")), p
        if path.is_file() and not path.is_symlink():
            with path.open("rb") as fh:
                assert fh.read(4) != b"\x7fELF", f"ELF binary committed: {p}"


def test_every_recorded_diann_file_is_small() -> None:
    for p in _tracked(str(RECORDED_DIANN.relative_to(REPO))):
        assert (REPO / p).stat().st_size <= R11_MAX_BYTES, p


def test_tests_data_stays_under_50_mb() -> None:
    total = sum((REPO / p).stat().st_size for p in _tracked("tests/data") if (REPO / p).is_file())
    assert total < TESTS_DATA_MAX, f"tests/data/ is {total / 2**20:.1f} MB (R11 < 50 MB)"


# ------------------------------------------------------------ the real samples (FarmShare only)


def _sample_dir(arm: str) -> Path:
    root = os.environ.get("GIAP_ROOT")
    if not root:
        pytest.skip(
            "GIAP_ROOT is not set: the DIA-NN samples live on FarmShare only (Khoa 2026-10-23)"
        )
    d = Path(root) / "pilot" / "recorded" / "diann" / "2.3.1" / arm
    if not (d / "report.sample.parquet").is_file():
        pytest.skip(f"no DIA-NN sample at {d} yet (M10-02 report job not run)")
    return d


@pytest.mark.parametrize("arm", ["A1", "A2"])
def test_the_real_sample_is_at_most_1000_rows_and_matches_its_column_list(arm: str) -> None:
    pq = pytest.importorskip("pyarrow.parquet")
    d = _sample_dir(arm)
    assert not d.resolve().is_relative_to(REPO.resolve()), "the sample must not be in the checkout"
    meta = pq.read_metadata(d / "report.sample.parquet")
    assert 0 < meta.num_rows <= 1000, meta.num_rows
    names = pq.read_schema(d / "report.sample.parquet").names
    committed = RECORDED_DIANN / "2.3.1" / arm / "schema.tsv"
    listing = committed if committed.is_file() else d / "schema.tsv"
    rows = [ln.split("\t") for ln in listing.read_text().splitlines()]
    if rows and rows[0] == SCHEMA_HEADER:
        rows = rows[1:]
    assert [r[0] for r in rows] == names, f"{listing} differs from the sample's columns"
