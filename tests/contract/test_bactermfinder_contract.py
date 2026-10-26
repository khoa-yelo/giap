"""Contract test: BacTermFinder parser (M2-15; GIAP_SPEC.md §5.4, §7.1.1, §7.1.5).

Active since fix-bactermfinder-parser: every sample directory under
``tests/data/recorded/bactermfinder/024870447f62/`` (pinned commit; genome FASTA, the tool's
``*_mean.csv`` and ``PROVENANCE.md``) is parsed. The output format is fixed by that fixture
(spec: "format fixed by the recorded fixture"); the generic contract checked here: every call
becomes a ``Terminator`` with ``source == "BacTermFinder"``, a valid interval on its replicon and
a confidence in [0, 1]
(Example 10 checks that such calls never count for ``terminator_between``).
"""

from __future__ import annotations

from pathlib import Path

import pytest

RECORDED = (
    Path(__file__).resolve().parents[1] / "data" / "recorded" / "bactermfinder" / "024870447f62"
)
SAMPLES = sorted(p for p in RECORDED.glob("*") if p.is_dir()) if RECORDED.is_dir() else []


def test_recordings_exist() -> None:
    assert len(SAMPLES) >= 2, RECORDED


@pytest.mark.parametrize("sample", SAMPLES, ids=[p.name for p in SAMPLES])
def test_recorded_bactermfinder_output_parses(sample: Path) -> None:
    from giap.annotate.terminators import parse_bactermfinder
    from giap.io.fasta import read_genome_fasta
    from giap.model.coords import Replicon, check_interval

    fastas = sorted(sample.glob("*.fa")) + sorted(sample.glob("*.fasta"))
    assert fastas, f"{sample}: the recording needs its genome FASTA"
    reps = {r.id: Replicon(r.id, len(r.sequence), "circular") for r in read_genome_fasta(fastas[0])}
    outs = sorted(p for p in sample.iterdir() if p.suffix not in {".fa", ".fasta", ".md", ".txt"})
    assert outs, f"{sample}: no BacTermFinder output file"
    calls = parse_bactermfinder(outs[0], reps)
    assert calls
    for t in calls:
        check_interval(t.interval, reps[t.interval.contig])
        assert t.source == "BacTermFinder" and 0.0 <= t.confidence <= 1.0
