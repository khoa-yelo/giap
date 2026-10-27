"""Contract test: TerminatorNet parser on recorded runs (M2-15; GIPPA_SPEC.md §5.4, §7.1.1, G12).

Recorded with the M2-10 harness (commit 72180309, ``-s``; PROVENANCE in
``tests/data/recorded/terminators/``): sim100k is committed whole; K-12 MG1655 and Nissle 1917
are read from ``$GIPPA_ROOT/recordings/terminatornet/<sample>/`` (``terminatornet.tsv`` +
``genome.fa``) and skipped when that directory is missing (CI).

Contract: every 11-column line becomes one ``Terminator``; ``NO TERMINATORS FOUND`` lines are
skipped; the interval is ``[start - 1, stop)`` on the forward strand for both strands, so the
recorded hairpin equals the genome slice (reverse-complemented for ``-``);
confidence = ``pred`` in [0, 1].
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from gippa.annotate.terminators import parse_terminatornet
from gippa.io.fasta import read_genome_fasta
from gippa.model.coords import Replicon, check_interval

DATA = Path(__file__).resolve().parents[1] / "data" / "recorded" / "terminators" / "terminatornet"
ROOT = Path(os.environ.get("GIPPA_ROOT", "/nonexistent")) / "recordings" / "terminatornet"
COMP = str.maketrans("ACGTN", "TGCAN")
SAMPLES = {"k12_mg1655": 8294, "nissle1917": 8644, "sim100k": 133}


def _genome(sample: str) -> dict[str, str] | None:
    fa = ROOT / sample / "genome.fa"
    if not fa.exists():
        return None
    return {r.id: r.sequence.upper() for r in read_genome_fasta(fa)}


def _check(tsv: Path, genome: dict[str, str] | None, reps: dict[str, Replicon], n: int) -> None:
    raw = [ln.split("\t") for ln in tsv.read_text().splitlines() if ln.strip()]
    calls = [f for f in raw if len(f) == 11]
    assert all(len(f) == 11 or f[-1].strip() == "NO TERMINATORS FOUND" for f in raw)
    got = parse_terminatornet(tsv, reps)
    assert len(got) == len(calls) == n
    assert {t.interval.strand for t in got} == {"+", "-"}
    hairpins = {(f[0], int(f[8]) - 1, int(f[9]), f[10]): f[3].upper() for f in calls}
    for t in got:
        iv = t.interval
        check_interval(iv, reps[iv.contig])
        assert t.source == "TerminatorNet" and 0.0 <= t.confidence <= 1.0
        assert t.hairpin_dg is not None and t.hairpin_dg < 0
        key = (iv.contig, iv.start, iv.end, iv.strand)
        assert key in hairpins
        if genome is not None:
            seq = genome[iv.contig][iv.start : iv.end]
            assert (seq if iv.strand == "+" else seq.translate(COMP)[::-1]) == hairpins[key]


def test_sim100k_committed_recording() -> None:
    reps = {
        "chr1": Replicon("chr1", 100_000, "circular"),
        "plasmid_1": Replicon("plasmid_1", 8_000, "circular"),
    }
    _check(DATA / "sim100k" / "terminatornet.tsv", _genome("sim100k"), reps, 133)


@pytest.mark.parametrize("sample", ["k12_mg1655", "nissle1917"])
def test_full_recordings_under_gippa_root(sample: str) -> None:
    tsv = ROOT / sample / "terminatornet.tsv"
    genome = _genome(sample)
    if not tsv.exists() or genome is None:
        pytest.skip(f"no recording at {ROOT / sample} ($GIPPA_ROOT/recordings missing)")
    reps = {c: Replicon(c, len(s), "circular") for c, s in genome.items()}
    _check(tsv, genome, reps, SAMPLES[sample])
