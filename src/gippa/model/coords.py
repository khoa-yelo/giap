"""0-based half-open coordinates, replicons and origin handling (GIPPA_SPEC.md §2.1 d16, §6.5).

Code works 0-based half-open internally; TSV, GFF3 and VCF writers convert to 1-based closed.
On a circular replicon an origin-crossing interval is stored unwrapped (``end > length``);
``crosses_origin`` is computed from the replicon, never stored. Linear contigs never wrap.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

Strand = Literal["+", "-", "."]
Topology = Literal["circular", "linear"]
RepliconKind = Literal["chromosome", "plasmid", "unknown"]

_STRAND_RANK: dict[str, int] = {"+": 0, "-": 1, ".": 2}


@dataclass(frozen=True, slots=True)
class Replicon:
    """A contig of the input genome with its length and topology."""

    id: str
    length: int
    topology: Topology
    kind: RepliconKind = "unknown"

    def __post_init__(self) -> None:
        if self.length < 1:
            raise ValueError(f"replicon {self.id!r}: length must be >= 1, got {self.length}")


@dataclass(frozen=True, slots=True)
class Interval:
    """A 0-based half-open interval ``[start, end)``; ``end`` may exceed the replicon length."""

    contig: str
    start: int
    end: int
    strand: Strand = "."

    def __post_init__(self) -> None:
        if type(self.start) is not int or type(self.end) is not int:
            raise TypeError(
                f"interval {self.contig}: start/end must be int, "
                f"got {type(self.start).__name__}/{type(self.end).__name__}"
            )
        if self.strand not in _STRAND_RANK:
            raise ValueError(f"interval {self.contig}: invalid strand {self.strand!r}")
        if not 0 <= self.start < self.end:
            raise ValueError(
                f"invalid interval {self.contig}:{self.start}-{self.end}: need 0 <= start < end"
            )

    @property
    def length_bp(self) -> int:
        return self.end - self.start

    def crosses_origin(self, replicon: Replicon) -> bool:
        return self.end > replicon.length


def check_interval(iv: Interval, replicon: Replicon) -> None:
    """Raise ``ValueError`` unless ``iv`` is a valid interval on ``replicon``."""
    where = f"{iv.contig}:{iv.start}-{iv.end}"
    if iv.contig != replicon.id:
        raise ValueError(f"interval {where} is not on replicon {replicon.id!r}")
    if iv.start >= replicon.length:
        raise ValueError(f"interval {where} starts beyond replicon length {replicon.length}")
    if replicon.topology == "linear":
        if iv.end > replicon.length:
            raise ValueError(f"interval {where} exceeds linear contig (linear contigs never wrap)")
    elif iv.length_bp > replicon.length:
        raise ValueError(f"interval {where} is longer than circular replicon ({replicon.length})")


def split_at_origin(iv: Interval, replicon: Replicon) -> tuple[Interval, ...]:
    """Split an origin-crossing interval into ``[start, length)`` and ``[0, end - length)``."""
    if not iv.crosses_origin(replicon):
        return (iv,)
    return (
        Interval(iv.contig, iv.start, replicon.length, iv.strand),
        Interval(iv.contig, 0, iv.end - replicon.length, iv.strand),
    )


def to_one_based_closed(iv: Interval) -> tuple[int, int]:
    """Return ``(start + 1, end)``; origin-crossing intervals keep ``end > length``."""
    return iv.start + 1, iv.end


def from_one_based_closed(
    contig: str, start1: int, end1: int, strand: Strand, replicon: Replicon
) -> Interval:
    """Convert 1-based closed coordinates; a wrapped ``start1 > end1`` is unwrapped if circular."""
    if start1 < 1 or end1 < 1:
        raise ValueError(f"{contig}:{start1}-{end1}: 1-based coordinates must be >= 1")
    end = end1
    if start1 > end1:
        if replicon.topology != "circular":
            raise ValueError(
                f"{contig}:{start1}-{end1}: start > end on a linear contig (never wraps)"
            )
        end = end1 + replicon.length
    iv = Interval(contig, start1 - 1, end, strand)
    check_interval(iv, replicon)
    return iv


def wrap_position(pos0: int, replicon: Replicon) -> int:
    """Map a 0-based (possibly unwrapped) position into ``[0, length)``."""
    if pos0 < 0:
        raise ValueError(f"negative position {pos0}")
    if pos0 < replicon.length:
        return pos0
    if replicon.topology != "circular":
        raise ValueError(f"position {pos0} beyond linear contig {replicon.id!r}")
    return pos0 % replicon.length


def interval_sort_key(
    iv: Interval, contig_order: Mapping[str, int], id_: str = ""
) -> tuple[int, int, int, int, str]:
    """Sort key: contig (input FASTA order), start, end, strand (``+ < - < .``), ID."""
    if iv.contig not in contig_order:
        raise ValueError(f"unknown contig {iv.contig!r} (not in contig order)")
    return (contig_order[iv.contig], iv.start, iv.end, _STRAND_RANK[iv.strand], id_)
