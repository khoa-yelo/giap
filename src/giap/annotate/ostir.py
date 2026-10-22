"""OSTIR input rows, input writer and output parser (GIAP_SPEC.md §7.1.1, §7.1.6; M2-16).

Each gene gets the window ``-up .. +down`` nt around its first start-codon base, strand-aware:
``+`` genes read ``genome[S - up, S + down)``, ``-`` genes ``revcomp(genome[E - down, E + up))``.
Circular replicons wrap (origin-crossing genes are stored unwrapped, ``end > length``); linear
contigs clip and ``start = end = 1 + (bases kept upstream)`` (OSTIR positions are 1-based).
The parser keeps only the rows at the requested starts of each id; ids without one are absent.
"""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from giap.model.coords import Replicon
from giap.model.genes import Gene

INPUT_HEADER = ("id", "seq", "start", "end", "anti-Shine-Dalgarno")
_COMP = str.maketrans("ACGTN", "TGCAN")


@dataclass(frozen=True, slots=True)
class OstirRow:
    id: str
    seq: str
    start: int
    end: int
    anti_sd: str


@dataclass(frozen=True, slots=True)
class RbsFeatures:
    tir_ostir: float
    sd_dg: float
    mrna_dg: float
    spacing_dg: float
    standby_dg: float
    spacer_bp: int
    start_codon: str


def _window(genome: str, lo: int, hi: int, circular: bool) -> str:
    """``genome[lo, hi)``, wrapping on circular replicons and clipping on linear ones."""
    n = len(genome)
    if not circular:
        return genome[max(lo, 0) : min(hi, n)]
    return "".join(genome[p % n] for p in range(lo, hi))


def ostir_row(
    gene: Gene, replicon: Replicon, genome_seq: str, anti_sd: str, up: int = 50, down: int = 50
) -> OstirRow:
    """The OSTIR input row of ``gene``: ``up`` nt upstream and ``down`` nt from its start."""
    iv, n = gene.interval, replicon.length
    circular = replicon.topology == "circular"
    genome = genome_seq.upper()
    if iv.strand == "+":
        lo, hi = iv.start - up, iv.start + down
        kept_up = up if circular else iv.start - max(lo, 0)
        seq = _window(genome, lo, hi, circular)
    else:
        lo, hi = iv.end - down, iv.end + up
        kept_up = up if circular else min(hi, n) - iv.end
        seq = _window(genome, lo, hi, circular).translate(_COMP)[::-1]
    return OstirRow(gene.locus_tag, seq, kept_up + 1, kept_up + 1, anti_sd)


def write_ostir_input(path: Path, rows: Sequence[OstirRow]) -> None:
    """``id,seq,start,end,anti-Shine-Dalgarno`` CSV, rows in the given order, LF line ends."""
    lines = [",".join(INPUT_HEADER)]
    lines += [f"{r.id},{r.seq},{r.start},{r.end},{r.anti_sd}" for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def read_ostir_input(path: Path) -> list[OstirRow]:
    """Rows of an ``ostir_input.csv`` written by :func:`write_ostir_input`, in file order."""
    with path.open(newline="", encoding="utf-8") as fh:
        return [
            OstirRow(r["id"], r["seq"], int(r["start"]), int(r["end"]), r["anti-Shine-Dalgarno"])
            for r in csv.DictReader(fh)
        ]


def _features(r: Mapping[str, str]) -> RbsFeatures:
    return RbsFeatures(
        tir_ostir=float(r["expression"]),
        sd_dg=float(r["dG_rRNA:mRNA"]),
        mrna_dg=float(r["dG_mRNA"]),
        spacing_dg=float(r["dG_spacing"]),
        standby_dg=float(r["dG_standby"]),
        spacer_bp=int(r["RBS_distance_bp"]),
        start_codon=r["start_codon"],
    )


def parse_ostir(path: Path, requested: Mapping[str, int]) -> dict[str, RbsFeatures]:
    """``ostir.csv`` → features per id, keeping the row whose ``start_position`` was requested."""
    many = parse_ostir_starts(path, {k: {v} for k, v in requested.items()})
    return {k: v[requested[k]] for k, v in many.items()}


def parse_ostir_starts(
    path: Path, requested: Mapping[str, set[int]]
) -> dict[str, dict[int, RbsFeatures]]:
    """``ostir.csv`` → features per id and requested 1-based ``start_position`` (M7-08): every
    requested start of an id with a row, the first row of a repeated ``(name, position)``."""
    out: dict[str, dict[int, RbsFeatures]] = {}
    with path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            name, pos = r["name"], int(r["start_position"])
            if pos not in requested.get(name, ()) or pos in out.get(name, {}):
                continue
            out.setdefault(name, {})[pos] = _features(r)
    return out
