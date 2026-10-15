"""Route 2a: proteins > 40 aa aligned with miniprot (GIAP_SPEC.md §7.2.3, §7.1.1; issue M3-12).

miniprot (MIT) is run as a separate program; this module only builds its argument vectors and
parses its GFF. §7.1.1 gotcha: ``-T 11`` is an indexing option, so the genetic code goes to the
index command and never to the alignment command.

Coordinates: ``mRNA`` columns 4-5 (1-based closed on the genome) and ``Target=ID qs qe`` (1-based
closed on the protein) are converted once to 0-based half-open (spec §6.5).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from giap.model.coords import Interval, Replicon, Strand, from_one_based_closed

Route = Literal["miniprot", "tblastn"]


def miniprot_index_cmd(genome: Path, mpi: Path, threads: int) -> list[str]:
    """``miniprot -t N -T 11 -d genome.mpi genome.fa`` (§7.1.1)."""
    return ["miniprot", "-t", str(threads), "-T", "11", "-d", str(mpi), str(genome)]


def miniprot_align_cmd(mpi: Path, proteins: Path, threads: int, args: str) -> list[str]:
    """``miniprot -t N <args> --gff genome.mpi proteins.faa``; ``args`` is whitespace-split."""
    return ["miniprot", "-t", str(threads), *args.split(), "--gff", str(mpi), str(proteins)]


@dataclass(frozen=True, slots=True)
class ProteinAlignment:
    """One protein alignment (one miniprot ``mRNA`` line or one merged tblastn hit)."""

    protein_id: str
    protein_len_aa: int
    route: Route
    interval: Interval
    rank: int
    identity: float
    positive: float | None
    q_start: int  # 0-based half-open on the protein
    q_end: int
    frameshift: int  # counts; > 0 -> flag fs
    stop_codon: int  # counts; > 0 -> flag st

    @property
    def query_coverage(self) -> float:
        return (self.q_end - self.q_start) / self.protein_len_aa


def _attributes(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in field.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def _parse_mrna(
    cols: list[str],
    where: str,
    replicons: Mapping[str, Replicon],
    protein_len: Mapping[str, int],
) -> ProteinAlignment:
    attrs = _attributes(cols[8])
    for key in ("Identity", "Target", "Rank"):
        if key not in attrs:
            raise ValueError(f"{where}: miniprot mRNA line has no {key} attribute")
    target = attrs["Target"].split()
    if len(target) < 3:
        raise ValueError(f"{where}: malformed Target={attrs['Target']!r}")
    pid = target[0]
    if pid not in protein_len:
        raise ValueError(f"{where}: protein {pid!r} has no known length")
    contig = cols[0]
    if contig not in replicons:
        raise ValueError(f"{where}: unknown contig {contig!r}")
    strand = cols[6]
    if strand not in ("+", "-"):
        raise ValueError(f"{where}: invalid strand {strand!r}")
    try:
        start1, end1 = int(cols[3]), int(cols[4])
        q_start1, q_end1 = int(target[1]), int(target[2])
        identity = float(attrs["Identity"])
        positive = float(attrs["Positive"]) if "Positive" in attrs else None
        rank = int(attrs["Rank"])
        frameshift = int(attrs.get("Frameshift", "0"))
        stop_codon = int(attrs.get("StopCodon", "0"))
    except ValueError as exc:
        raise ValueError(f"{where}: {exc}") from exc
    if start1 > end1:
        raise ValueError(f"{where}: start {start1} > end {end1}")
    interval = from_one_based_closed(contig, start1, end1, _strand(strand), replicons[contig])
    plen = protein_len[pid]
    if not 1 <= q_start1 <= q_end1 <= plen:
        raise ValueError(f"{where}: Target {q_start1}-{q_end1} outside protein length {plen}")
    return ProteinAlignment(
        protein_id=pid,
        protein_len_aa=plen,
        route="miniprot",
        interval=interval,
        rank=rank,
        identity=identity,
        positive=positive,
        q_start=q_start1 - 1,
        q_end=q_end1,
        frameshift=frameshift,
        stop_codon=stop_codon,
    )


def _strand(value: str) -> Strand:
    return "+" if value == "+" else "-"


def parse_miniprot_gff(
    path: Path, replicons: Mapping[str, Replicon], protein_len: Mapping[str, int]
) -> list[ProteinAlignment]:
    """One ``ProteinAlignment`` per ``mRNA`` line, in file order; ``##PAF`` lines are ignored.

    Raises ``ValueError`` naming the line when ``Identity`` (or ``Target``/``Rank``) is missing,
    the contig is unknown or the protein has no length (contract drift).
    """
    out: list[ProteinAlignment] = []
    with path.open(encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "mRNA":
                continue
            where = f"{path.name} line {lineno}"
            out.append(_parse_mrna(cols, where, replicons, protein_len))
    return out
