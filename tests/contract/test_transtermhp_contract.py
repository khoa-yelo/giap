"""Contract test: TransTermHP 2.09 ``.tt`` parser and ``.coords`` writer on recorded runs
(M2-15; GIAP_SPEC.md §5.4, §7.1.1, §7.1.5).

sim100k is committed in ``tests/data/recorded/terminators/transtermhp/sim100k/``; K-12 MG1655 and
Nissle 1917 are read from ``$GIAP_ROOT/recordings/transtermhp/<sample>/`` (``transtermhp.tt``,
``genes.coords``, ``genome.fa``) with the Bakta JSON of
``$GIAP_ROOT/recordings/bakta_full/<sample>/``, skipped when missing (CI).

Contract: one ``Terminator`` per ``  TERM`` line, contig from the preceding ``SEQUENCE`` line,
interval ``[min(a, b) - 1, max(a, b))``, strand from the line, confidence = CONF / 100; and
``write_transterm_coords`` on the Bakta gene model reproduces, byte for byte, the ``genes.coords``
TransTermHP was run on.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from giap.annotate.bakta import parse_bakta_json
from giap.annotate.terminators import parse_transterm_tt, write_transterm_coords
from giap.io.fasta import read_genome_fasta
from giap.model.coords import Replicon, check_interval

DATA = Path(__file__).resolve().parents[1] / "data" / "recorded" / "terminators" / "transtermhp"
REC = Path(os.environ.get("GIAP_ROOT", "/nonexistent")) / "recordings"
SAMPLES = {"k12_mg1655": 1847, "nissle1917": 1815}


def _raw(tt: Path) -> list[tuple[str, int, int, str, int]]:
    out, contig = [], ""
    for ln in tt.read_text().splitlines():
        if ln.startswith("SEQUENCE "):
            contig = ln.split()[1]
        elif ln.startswith("  TERM "):
            f = ln.split()
            out.append((contig, int(f[2]), int(f[4]), f[5], int(f[7])))
    return out


def _check(tt: Path, reps: dict[str, Replicon], n: int) -> None:
    raw = _raw(tt)
    got = parse_transterm_tt(tt, reps)
    assert len(got) == len(raw) == n
    want = sorted((c, min(a, b) - 1, max(a, b), s, conf) for c, a, b, s, conf in raw)
    seen = sorted(
        (t.interval.contig, t.interval.start, t.interval.end, t.interval.strand,
         round(t.confidence * 100))
        for t in got
    )  # fmt: skip
    assert seen == want
    for t in got:
        check_interval(t.interval, reps[t.interval.contig])
        assert t.source == "TransTermHP" and 0.0 <= t.confidence <= 1.0


def test_sim100k_committed_recording() -> None:
    reps = {
        "chr1": Replicon("chr1", 100_000, "circular"),
        "plasmid_1": Replicon("plasmid_1", 8_000, "circular"),
    }
    _check(DATA / "sim100k" / "transtermhp.tt", reps, 30)


def _sample(sample: str) -> Path:
    d = REC / "transtermhp" / sample
    if not (d / "transtermhp.tt").exists():
        pytest.skip(f"no recording at {d} ($GIAP_ROOT/recordings missing)")
    return d


@pytest.mark.parametrize("sample", sorted(SAMPLES))
def test_full_tt_recordings_under_giap_root(sample: str) -> None:
    d = _sample(sample)
    reps = {
        r.id: Replicon(r.id, len(r.sequence), "circular")
        for r in read_genome_fasta(d / "genome.fa")
    }
    _check(d / "transtermhp.tt", reps, SAMPLES[sample])


@pytest.mark.parametrize("sample", [*sorted(SAMPLES), "sim100k"])
def test_coords_writer_reproduces_recorded_genes_coords(sample: str, tmp_path: Path) -> None:
    d = _sample(sample)
    jsons = sorted((REC / "bakta_full" / sample / "bakta").glob("*.json"))
    if not jsons:
        pytest.skip(f"no Bakta JSON under {REC / 'bakta_full' / sample}")
    model = parse_bakta_json(jsons[0]).model
    out = tmp_path / "genes.coords"
    write_transterm_coords(out, list(model.genes))
    assert out.read_bytes() == (d / "genes.coords").read_bytes()
