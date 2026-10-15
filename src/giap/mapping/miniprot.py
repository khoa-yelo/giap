"""Route 2a: proteins > 40 aa aligned with miniprot (GIAP_SPEC.md §7.2.3, §7.1.1; issue M3-12).

miniprot (MIT) is run as a separate program; this module only builds its argument vectors and
parses its GFF. §7.1.1 gotcha: ``-T 11`` is an indexing option, so the genetic code goes to the
index command and never to the alignment command.

Coordinates: ``mRNA`` columns 4-5 (1-based closed on the genome) and ``Target=ID qs qe`` (1-based
closed on the protein) are converted once to 0-based half-open (spec §6.5).

Alternative starts (M3-12 review decision 1): miniprot scores a GTG/TTG start codon translated as
``M`` as a substitution (``cs:Z:*gtgM:258``). The ``cs`` tag of the ``##PAF`` line preceding an
``mRNA`` line is read; a lone ``*[acgt]tgM`` at query position 0 counts as a match for identity
and sets ``ProteinAlignment.alt_start``.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, replace
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
    alt_start: bool = False  # GTG/TTG start read as M, counted as a match in ``identity``
    coverage: float | None = None  # explicit query coverage (tblastn: union of HSP ranges)

    @property
    def query_coverage(self) -> float:
        """``coverage`` if set (M3-13 lead decision 1), else the span ``(q_end - q_start) / n``."""
        if self.coverage is not None:
            return self.coverage
        return (self.q_end - self.q_start) / self.protein_len_aa


_CS_OP = re.compile(
    r":(\d+)|\*([a-z]{3})([A-Z*])|\+([A-Z*]+)|-([a-z]+)|~[a-z]{2}\d+[a-z]{2}|([a-z]+)"
)
_ALT_START = re.compile(r"\*[acgt]tgM")
_UNSAFE_ID = (";", "=")


@dataclass(frozen=True, slots=True)
class _Paf:
    """The fields of a miniprot ``##PAF`` line needed for the following ``mRNA`` line."""

    query: str
    q_start: int
    cs: str | None


def _parse_paf(line: str, where: str) -> _Paf:
    cols = line.split("\t")
    if len(cols) < 13:
        raise ValueError(f"{where}: malformed ##PAF line")
    query = cols[1]
    _check_id(query, where)
    try:
        q_start = int(cols[3])
    except ValueError as exc:
        raise ValueError(f"{where}: {exc}") from exc
    cs = next((c[5:] for c in cols[13:] if c.startswith("cs:Z:")), None)
    return _Paf(query, q_start, cs)


def _check_id(pid: str, where: str) -> None:
    if any(ch in pid for ch in _UNSAFE_ID):
        raise ValueError(f"{where}: protein ID {pid!r} contains ';' or '=' (unsafe in GFF)")


def _matches_columns(cs: str) -> tuple[int, int]:
    """Identical residues and alignment columns (residues, codon-sized gaps) of a ``cs`` string."""
    matches = columns = 0
    for m in _CS_OP.finditer(cs):
        if m.group(1) is not None:
            matches += int(m.group(1))
            columns += int(m.group(1))
        elif m.group(2) is not None:
            columns += 1
        elif m.group(4) is not None:
            columns += len(m.group(4))
        elif m.group(5) is not None:
            columns += len(m.group(5)) // 3
    return matches, columns


def _alt_start_identity(paf: _Paf | None, pid: str, q_start: int) -> float | None:
    """``(matches + 1) / columns`` when the alignment opens with a lone ``*[acgt]tgM`` at 0."""
    if paf is None or paf.cs is None or paf.query != pid or paf.q_start != 0 or q_start != 0:
        return None
    if not _ALT_START.match(paf.cs):
        return None
    matches, columns = _matches_columns(paf.cs)
    return None if columns == 0 else min(1.0, (matches + 1) / columns)


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
    _check_id(pid, where)
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
    """One ``ProteinAlignment`` per ``mRNA`` line, in file order.

    The ``##PAF`` line preceding an ``mRNA`` line applies to it only: a lone ``*[acgt]tgM`` at
    query position 0 in its ``cs`` tag is an alternative start, counted as a match
    (identity ``(matches + 1) / columns`` from ``cs``) and recorded in ``alt_start``.

    Raises ``ValueError`` naming the line when ``Identity`` (or ``Target``/``Rank``) is missing,
    the contig is unknown, the protein has no length (contract drift) or a protein ID contains
    ``;`` or ``=``.
    """
    out: list[ProteinAlignment] = []
    paf: _Paf | None = None
    with path.open(encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n\r")
            where = f"{path.name} line {lineno}"
            if line.startswith("##PAF\t"):
                paf = _parse_paf(line, where)
                continue
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "mRNA":
                continue
            aln = _parse_mrna(cols, where, replicons, protein_len)
            identity = _alt_start_identity(paf, aln.protein_id, aln.q_start)
            paf = None
            if identity is not None:
                aln = replace(aln, identity=identity, alt_start=True)
            out.append(aln)
    return out
