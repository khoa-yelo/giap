"""Contract test: blastp (BLAST+ 2.16.0) against the contaminant library on a recorded run (M3-15).

Marked ``pending_sample`` until ``tests/data/recorded/blast/2.16.0/contaminants/`` exists. The
directory holds ``blastp.tsv`` (tabular, ``gippa.mapping.tblastn.OUTFMT``, written by the argument
vector of ``blastp_cmd``), the query ``proteins.faa`` and ``PROVENANCE.md`` (command, image digest,
input sha256, date, host). The library is ``tests/data/contaminants_mini.faa``. Query IDs encode
the truth:

- ``HIT_<ACC>_<k>``: a library entry with random substitutions to 93-99 % identity, full length
  (expected match ``<ACC>``);
- ``LOWCOV_<ACC>_<k>``: 40 % of a library entry fused to 60 % random sequence (not a contaminant:
  coverage < 0.5, example 10);
- ``RAND_<k>``: random sequences (no passing alignment).
"""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

from gippa.io.protein_fasta import read_protein_fasta
from gippa.mapping.tblastn import OUTFMT

DATA = Path(__file__).resolve().parents[1] / "data"
RECORDED = DATA / "recorded" / "blast" / "2.16.0" / "contaminants"
MINI = DATA / "contaminants_mini.faa"

pytestmark = [
    pytest.mark.pending_sample,
    pytest.mark.skipif(
        not (RECORDED / "blastp.tsv").is_file(),
        reason=f"pending_sample: no recorded blastp run in {RECORDED}",
    ),
]


def test_provenance_present() -> None:
    text = (RECORDED / "PROVENANCE.md").read_text()
    assert "blastp" in text and OUTFMT in text
    assert "sha256" in text


def test_recorded_rows_have_outfmt_columns() -> None:
    n = len(OUTFMT.split()) - 1
    lines = [x for x in (RECORDED / "blastp.tsv").read_text().splitlines() if x.strip()]
    assert lines
    assert all(len(x.split("\t")) == n for x in lines)


def test_recorded_matches_truth() -> None:
    c = importlib.import_module("gippa.mapping.contaminants")
    lib_ids = {r.id for r in read_protein_fasta(MINI)}
    queries = [r.id for r in read_protein_fasta(RECORDED / "proteins.faa")]
    got = c.match_contaminant_proteins(RECORDED / "blastp.tsv")
    hits = [q for q in queries if q.startswith("HIT_")]
    assert hits and any(q.startswith("LOWCOV_") for q in queries)
    for q in hits:
        acc = q.split("_")[1]
        assert acc in lib_ids
        assert got.get(q) == acc, q
    for q in queries:
        if not q.startswith("HIT_"):
            assert q not in got, q
    assert set(got.values()) <= lib_ids


def test_cmd_matches_recorded_outfmt(tmp_path: Path) -> None:
    c = importlib.import_module("gippa.mapping.contaminants")
    cmd = c.blastp_cmd(tmp_path / "q.faa", MINI, 2)
    assert cmd[cmd.index("-outfmt") + 1] == OUTFMT
