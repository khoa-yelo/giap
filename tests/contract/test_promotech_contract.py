"""Contract test: Promotech targeted (``-s``) output parser on recorded runs (M2-13 example 10).

Recordings live in ``tests/data/recorded/promotech/56251ad9/<input>/`` (pinned commit
``56251ad9b883ef831b4753fc623d5ec970fe65e0``) with the exact input chunk(s) ``windows_*.fa``, the
tool's ``sequences_predictions.csv`` per chunk (``<chunk stem>.sequences_predictions.csv`` or a
single ``sequences_predictions.csv``) and ``PROVENANCE.md``.

- ``annotate-real-1847900`` (committed): a real ``-s -m RF-HOT`` run from the reviewer's ANNOTATE
  job 1847900 on a stub-written chunk (locus-tag headers). Always runs.
- ``k12_mg1655_excerpt`` (committed, M2-21): the first 100 records of the GIAP-written chunk
  ``windows_000.fa`` (``w0 …``) and the last 100 of ``windows_009.fa`` (``… w950227``) of the real
  K-12 recording (``record_tools.sbatch`` job 1849513_6), with their prediction rows. Always runs.
  (The ``sim100k`` recording under ``$GIAP_ROOT`` is a genome-wide ``-g`` run, not GIAP's ``-s``
  mode; see the excerpt's ``PROVENANCE.md``.)

Real format (pinned ``sequences/process_sequences.py``): tab-separated despite the ``.csv`` name,
header ``CHROM<TAB>SEQ<TAB>PRED``, one row per input record in input order, CHROM = the FASTA id.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

from helpers.promoter_oracle import api

RECORDED = Path(__file__).resolve().parents[1] / "data" / "recorded" / "promotech" / "56251ad9"
REAL = RECORDED / "annotate-real-1847900"
K12 = RECORDED / "k12_mg1655_excerpt"
HEADER = ["CHROM", "SEQ", "PRED"]


def _fasta(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for ln in path.read_text().splitlines():
        if ln.startswith(">"):
            out.append((ln[1:].split()[0], ""))
        elif ln.strip():
            out[-1] = (out[-1][0], out[-1][1] + ln.strip())
    return out


def _pairs(sample: Path) -> list[tuple[Path, Path]]:
    fas = sorted(sample.glob("windows_*.fa"))
    single = sample / "sequences_predictions.csv"
    if len(fas) == 1 and single.is_file():
        return [(fas[0], single)]
    return [(fa, sample / f"{fa.stem}.sequences_predictions.csv") for fa in fas]


def _check_sample(sample: Path) -> list[str]:
    prov = (sample / "PROVENANCE.md").read_text()
    assert "56251ad9b883ef831b4753fc623d5ec970fe65e0" in prov
    pairs = _pairs(sample)
    assert pairs, f"{sample}: no windows_*.fa"
    headers: list[str] = []
    for fa, csv in pairs:
        assert hashlib.sha256(fa.read_bytes()).hexdigest() in prov
        records = _fasta(fa)
        lines = csv.read_text().splitlines()
        assert lines[0].split("\t") == HEADER
        rows = [ln.split("\t") for ln in lines[1:] if ln]
        assert [(r[0], r[1]) for r in rows] == [(h, s.upper()) for h, s in records]
        assert all(len(s) == 40 for _, s in records)
        scores = api("parse_promotech")(csv)
        assert list(scores) == [h for h, _ in records]
        for r in rows:
            assert scores[r[0]] == float(r[2])
            assert 0.0 <= scores[r[0]] <= 1.0
        headers += list(scores)
    return headers


def test_real_annotate_run_parses() -> None:
    headers = _check_sample(REAL)
    assert len(headers) == 86


def test_k12_excerpt_has_provenance_for_every_file() -> None:
    prov = (K12 / "PROVENANCE.md").read_text()
    files = sorted(K12.glob("windows_*"))
    assert len(files) == 4
    for f in files:
        assert hashlib.sha256(f.read_bytes()).hexdigest() in prov, f.name
    assert "1849513" in prov and "-s -m RF-HOT" in prov


def test_k12_giap_w_headers_round_trip_unmangled() -> None:
    """M2-13 example 10 on GIAP-written windows: ``w<index>`` keys survive Promotech unchanged."""
    headers = _check_sample(K12)
    first = [f"w{i}" for i in range(100)]
    last = [f"w{i}" for i in range(950128, 950228)]
    assert headers == first + last
