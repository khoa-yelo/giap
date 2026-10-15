"""Promotech targeted command line and output parser (GIAP_SPEC.md §7.1.1; M2-13).

Promotech is GPL-3.0: it is never imported, only run as a separate process from its own
environment (``python promotech.py -s -m RF-HOT -f windows.fa -o out/``). Its output
``sequences_predictions.csv`` is tab-separated despite the name, with header
``CHROM<TAB>SEQ<TAB>PRED`` and one row per input record (CHROM = the FASTA id).
"""

from __future__ import annotations

from pathlib import Path

PROMOTECH_MODEL = "RF-HOT"
PROMOTECH_HEADER = ["CHROM", "SEQ", "PRED"]


def promotech_command(chunk: Path, outdir: Path, promotech_home: Path) -> list[str]:
    """``python <home>/promotech.py -s -m RF-HOT -f <chunk> -o <outdir>`` (targeted, 40-nt)."""
    return [
        "python",
        str(promotech_home / "promotech.py"),
        "-s",
        "-m",
        PROMOTECH_MODEL,
        "-f",
        str(chunk),
        "-o",
        str(outdir),
    ]


def parse_promotech(path: Path) -> dict[str, float]:
    """``sequences_predictions.csv`` → ``{header: score}`` (file order); ValueError if foreign."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].split("\t") != PROMOTECH_HEADER:
        got = lines[0] if lines else "<empty>"
        raise ValueError(f"{path}: expected header CHROM<TAB>SEQ<TAB>PRED, got {got!r}")
    scores: dict[str, float] = {}
    for n, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError(f"{path}:{n}: expected 3 tab-separated fields, got {len(fields)}")
        scores[fields[0]] = float(fields[2])
    return scores
