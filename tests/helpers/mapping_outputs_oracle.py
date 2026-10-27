"""Independent oracle for the Module 2 output writers (docs/specs/M3-16-mapping-outputs.md).

Written from GIPPA_SPEC.md §6.5/§9.2/§9.3 and the M3-16 spec, without importing gippa.io: input
rows are plain dicts with 0-based half-open, unwrapped ``start``/``end``; outputs are the expected
TSV coordinates, BED12 core fields (chrom..strand, blockCount, blockSizes, blockStarts) and sort
keys.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

STRAND_RANK = {"+": 0, "-": 1, ".": 2}

PROBED_UNIQUENESS = {
    "unique": "unique",
    "multi_identical": "not-unique[same-set]",
    "shared": "not-unique[conflict]",
}


def tsv_coords(start0: int, end0: int, length: int) -> tuple[str, str, str]:
    """``(start, end, crosses_origin)`` as written in a TSV: 1-based closed, ``end > L`` kept."""
    return str(start0 + 1), str(end0), "true" if end0 > length else "false"


def parse_blocks(blocks: str | None, start: int, end: int) -> list[tuple[int, int]]:
    """Absolute ``[s, e)`` blocks of a hit; ``None`` means one block spanning the hit.

    For ``δ = -1`` (second block starts before the first ends) the first block is shortened by
    1 nt so BED12 blocks do not overlap (M3-16 spec, M3-08 edge case).
    """
    if blocks is None:
        return [(start, end)]
    out = []
    for part in blocks.split(";"):
        s, e = part.split("-")
        out.append((int(s), int(e)))
    out.sort()
    if len(out) == 2 and out[1][0] < out[0][1]:
        out[0] = (out[0][0], out[1][0])
    return out


def bed12_core(
    contig: str,
    start: int,
    end: int,
    strand: str,
    name: str,
    score: int,
    length: int,
    blocks: str | None = None,
) -> list[tuple[str, ...]]:
    """BED12 core fields per output line (split at the origin into ``_part1``/``_part2``).

    Each tuple is ``(chrom, chromStart, chromEnd, name, score, strand, blockCount, blockSizes,
    blockStarts)``; block lists without a trailing comma.
    """
    absb = parse_blocks(blocks, start, end)
    if end <= length:
        parts = [(start, end, name, 0)]
    else:
        parts = [(start, length, f"{name}_part1", 0), (length, end, f"{name}_part2", length)]
    out = []
    for lo, hi, nm, off in parts:
        bl = [(max(s, lo), min(e, hi)) for s, e in absb if max(s, lo) < min(e, hi)]
        b0 = bl[0][0]
        out.append(
            (
                contig,
                str(b0 - off),
                str(bl[-1][1] - off),
                nm,
                str(score),
                strand,
                str(len(bl)),
                ",".join(str(e - s) for s, e in bl),
                ",".join(str(s - b0) for s, _ in bl),
            )
        )
    return out


def bed_core_of_line(fields: Sequence[str]) -> tuple[str, ...]:
    """The oracle's core tuple from a written BED12/proBed line (trailing commas stripped)."""
    return (
        *fields[0:6],
        fields[9],
        fields[10].rstrip(","),
        fields[11].rstrip(","),
    )


def bed_sort_key(core: tuple[str, ...], order: Mapping[str, int]) -> tuple[int, int, int, int, str]:
    return (order[core[0]], int(core[1]), int(core[2]), STRAND_RANK[core[5]], core[3])


def hit_sort_key(
    row: Mapping[str, object], order: Mapping[str, int], id_col: str
) -> tuple[int, int, int, int, str]:
    """§6.5 / M3-16 sort: contig (FASTA order), start, end, strand, then peptide/protein_id."""
    return (
        order[str(row["contig"])],
        int(row["start"]),  # type: ignore[call-overload]
        int(row["end"]),  # type: ignore[call-overload]
        STRAND_RANK[str(row["strand"])],
        str(row[id_col]),
    )


def fmt6g(x: float) -> str:
    """§9.2 float format (``%.6g``; ``0`` for zero)."""
    return "0" if x == 0 else f"{x:.6g}"


def frameshift_base(blocks: str) -> tuple[str, int]:
    """The GFF3 attribute recording a frameshift's odd base (M3-16 lead decision 2).

    ``δ = -1`` (blocks overlap by 1 nt) -> ``("frameshift_reread", pos)``; ``δ = +1`` (a 1-nt gap)
    -> ``("frameshift_skipped", pos)``; ``pos`` is the 1-based genome position of that base.
    """
    (_, e1), (s2, _) = sorted(tuple(int(x) for x in p.split("-")) for p in blocks.split(";"))
    if s2 == e1 - 1:
        return "frameshift_reread", s2 + 1
    if s2 == e1 + 1:
        return "frameshift_skipped", e1 + 1
    raise ValueError(f"not a +-1 frameshift: {blocks!r}")


def covered_mod(spans: Sequence[tuple[int, int]], length: int) -> list[int]:
    """Sorted 0-based genome positions of unwrapped ``[s, e)`` spans, folded mod ``length``."""
    return sorted({p % length for s, e in spans for p in range(s, e)})
