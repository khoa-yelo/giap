"""Contract test: DIA-NN 2.x converter against a recorded sample (GIAP_SPEC.md §5.4, §11.4).

Marked ``pending_sample`` until the pilot (§15) supplies a real ``report.parquet``/``report.tsv``
under ``tests/data/recorded/diann/``. If an ``expected_peptides.tsv`` (and optionally
``expected_abundance.tsv``) sits next to the report, the converter output must match it byte for
byte.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from giap.config.model import InputConfig
from giap.io.peptides import PEPTIDE_COLUMNS
from giap.io.tsv import write_tsv

RECORDED = Path(__file__).resolve().parents[1] / "data" / "recorded" / "diann"
REPORTS = (
    sorted(p for p in RECORDED.iterdir() if p.name in ("report.parquet", "report.tsv"))
    if RECORDED.is_dir()
    else []
)

pytestmark = [
    pytest.mark.pending_sample,
    pytest.mark.skipif(
        not REPORTS, reason=f"pending_sample: no recorded DIA-NN report in {RECORDED}"
    ),
]


@pytest.mark.parametrize("report", REPORTS, ids=[p.name for p in REPORTS])
def test_recorded_report_converts(report: Path, tmp_path: Path) -> None:
    from giap.io.converters.diann import convert_diann

    res = convert_diann([report], search_db=None, cfg=InputConfig())
    assert not [i for i in res.issues if i.severity == "error"], res.issues
    df = res.peptides
    assert len(df) > 0
    cols = list(df.columns)
    assert cols == [c for c in PEPTIDE_COLUMNS if c in cols]
    for required in ("peptide", "modified_peptide", "protein_ids", "q_value", "run", "source"):
        assert required in cols
    assert df["peptide"].map(lambda s: bool(re.fullmatch(r"[A-Z]+", s))).all()
    assert set(df["source"]) == {"database"}
    assert ((df["q_value"] >= 0) & (df["q_value"] <= 1)).all()

    expected = RECORDED / "expected_peptides.tsv"
    if expected.exists():
        out = tmp_path / "peptides.tsv"
        write_tsv(out, cols, df.itertuples(index=False, name=None))
        assert out.read_bytes() == expected.read_bytes()
    expected_ab = RECORDED / "expected_abundance.tsv"
    if expected_ab.exists():
        assert res.abundance is not None
        out = tmp_path / "abundance.tsv"
        write_tsv(
            out, list(res.abundance.columns), res.abundance.itertuples(index=False, name=None)
        )
        assert out.read_bytes() == expected_ab.read_bytes()
