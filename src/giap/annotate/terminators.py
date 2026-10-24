"""Terminator layer: TerminatorNet, TransTermHP, BacTermFinder, cross-check and
``terminator_between`` (GIAP_SPEC.md §7.1.1, §7.1.5, G12; M2-15).

Parsers return 0-based half-open ``Terminator`` intervals validated against the replicons.
Between/downstream windows are evaluated in unwrapped coordinates on circular replicons, so a
window crossing the origin matches calls stored on either side of it. BacTermFinder calls are
source-tagged and never count for ``terminator_between`` / ``terminator_conf_downstream``.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

from giap.io._open import open_text
from giap.io.bed import BedRecord
from giap.model.coords import Interval, Replicon, Strand, check_interval
from giap.model.genes import Gene

Source = Literal["TerminatorNet", "TransTermHP", "BacTermFinder"]
NOT_PRIMARY: Source = "BacTermFinder"
DOWNSTREAM_NO_NEXT_BP = 300  # M2-15 Q3


@dataclass(frozen=True, slots=True)
class Terminator:
    interval: Interval
    confidence: float
    source: Source
    hairpin_dg: float | None = None


@dataclass(frozen=True, slots=True)
class CrossRow:
    """One row of ``annotation/terminators_crosscheck.tsv`` (M2-15 Q2)."""

    contig: str
    n_terminatornet: int
    n_transtermhp: int
    n_matched: int
    share_matched: float


def _validated(iv: Interval, replicons: Mapping[str, Replicon], where: str) -> Interval:
    rep = replicons.get(iv.contig)
    if rep is None:
        raise ValueError(f"{where}: unknown contig {iv.contig!r}")
    check_interval(iv, rep)
    return iv


def _strand(s: str, where: str) -> Strand:
    if s not in ("+", "-"):
        raise ValueError(f"{where}: invalid strand {s!r}")
    return cast(Strand, s)


_TN_HEADER = (
    "name", "TERM_i", "A_tract", "hairpin", "U_tract", "structure", "GH", "pred", "start", "stop",
    "strand",
)  # fmt: skip


def _is_tn_header(fields: Sequence[str]) -> bool:
    """Exactly the documented column names (with or without ``strand``), never a prefix match."""
    got = tuple(c.strip() for c in fields)
    return got in (_TN_HEADER, _TN_HEADER[:-1])


def parse_terminatornet(path: Path, replicons: Mapping[str, Replicon]) -> list[Terminator]:
    """TerminatorNet ``-s`` output (TAB, 11 columns, 1-based closed forward coordinates). Only
    blank, ``NO TERMINATORS FOUND`` and exact header lines are skipped (M2-15 review item 4).
    Plain or gzip (``open_text``; the K-12 recording is ``terms.tsv.gz``)."""
    out: list[Terminator] = []
    with open_text(path) as fh:
        text = fh.read()
    for n, line in enumerate(text.splitlines(), start=1):
        f = line.rstrip("\r\n").split("\t")
        if not line.strip() or "NO TERMINATORS FOUND" in line or _is_tn_header(f):
            continue
        where = f"{path.name}:{n}"
        if len(f) != 11:
            raise ValueError(f"{where}: expected 11 columns (run with -s), got {len(f)}")
        iv = Interval(f[0], int(f[8]) - 1, int(f[9]), _strand(f[10], where))
        out.append(
            Terminator(_validated(iv, replicons, where), float(f[7]), "TerminatorNet", float(f[6]))
        )
    return out


def write_transterm_coords(path: Path, genes: Sequence[Gene]) -> None:
    """``<locus_tag> <start1> <end1> <contig>``; minus-strand genes written with start > end."""
    lines = []
    for g in genes:
        iv = g.interval
        a, b = (iv.start + 1, iv.end) if iv.strand != "-" else (iv.end, iv.start + 1)
        lines.append(f"{g.locus_tag} {a} {b} {iv.contig}\n")
    path.write_text("".join(lines))


def parse_transterm_tt(path: Path, replicons: Mapping[str, Replicon]) -> list[Terminator]:
    """TransTermHP ``.tt``: ``  TERM <n> <a> - <b> <strand> <loc> <conf> …``; conf / 100."""
    out: list[Terminator] = []
    contig: str | None = None
    for n, line in enumerate(path.read_text().splitlines(), start=1):
        if line.startswith("SEQUENCE "):
            contig = line.split()[1]
            continue
        f = line.split()
        if not line[:1].isspace() or not f or f[0] != "TERM":
            continue
        where = f"{path.name}:{n}"
        if contig is None:
            raise ValueError(f"{where}: TERM line before any SEQUENCE line")
        a, b = int(f[2]), int(f[4])
        iv = Interval(contig, min(a, b) - 1, max(a, b), _strand(f[5], where))
        out.append(Terminator(_validated(iv, replicons, where), int(f[7]) / 100, "TransTermHP"))
    return out


_BTF_COLS = {
    "contig": ("contig", "chrom", "seq_id", "seqid", "chromosome"),
    "start": ("start",),
    "end": ("end", "stop"),
    "strand": ("strand",),
    "confidence": ("probability", "prob", "confidence", "score"),
}


def parse_bactermfinder(path: Path, replicons: Mapping[str, Replicon]) -> list[Terminator]:
    """BacTermFinder calls (provisional until a run is recorded: a CSV/TSV with a header naming
    contig, 1-based closed start/end, strand and probability columns)."""
    text = path.read_text().splitlines()
    if not text:
        return []
    delim = "\t" if "\t" in text[0] else ","
    rows = list(csv.reader(text, delimiter=delim))
    header = [h.strip().lower() for h in rows[0]]
    idx: dict[str, int] = {}
    for key, names in _BTF_COLS.items():
        hit = [header.index(c) for c in names if c in header]
        if not hit:
            raise ValueError(f"{path.name}: no {key} column in header {rows[0]}")
        idx[key] = hit[0]
    out = []
    for n, r in enumerate(rows[1:], start=2):
        if not any(c.strip() for c in r):
            continue
        where = f"{path.name}:{n}"
        iv = Interval(
            r[idx["contig"]],
            int(r[idx["start"]]) - 1,
            int(r[idx["end"]]),
            _strand(r[idx["strand"]].strip(), where),
        )
        conf = float(r[idx["confidence"]])
        out.append(Terminator(_validated(iv, replicons, where), conf, "BacTermFinder"))
    return out


def _gap(a: Interval, b: Interval) -> int:
    return max(0, max(a.start, b.start) - min(a.end, b.end))


def _row(contig: str, prim: list[Terminator], oth: list[Terminator], tol: int) -> CrossRow:
    m = sum(
        any(
            o.interval.strand == p.interval.strand and _gap(o.interval, p.interval) <= tol
            for o in oth
        )
        for p in prim
    )
    share = m / len(prim) if prim else math.nan
    return CrossRow(contig, len(prim), len(oth), m, share)


def crosscheck(
    primary: Sequence[Terminator],
    other: Sequence[Terminator],
    tol_bp: int = 20,
    contig_order: Sequence[str] | None = None,
) -> list[CrossRow]:
    """Share of primary calls within ``tol_bp`` of an ``other`` call (same contig and strand);
    one row per contig and a last ``all`` row. Contigs follow ``contig_order`` (FASTA order);
    contigs absent from it (all of them without an order) follow, sorted by name."""
    seen = {t.interval.contig for t in (*primary, *other)}
    first = [c for c in dict.fromkeys(contig_order or ()) if c in seen]
    contigs = first + sorted(seen.difference(first))
    rows = []
    for c in contigs:
        p = [t for t in primary if t.interval.contig == c]
        o = [t for t in other if t.interval.contig == c]
        rows.append(_row(c, p, o, tol_bp))
    n_p = sum(r.n_terminatornet for r in rows)
    n_m = sum(r.n_matched for r in rows)
    n_o = sum(r.n_transtermhp for r in rows)
    rows.append(CrossRow("all", n_p, n_o, n_m, n_m / n_p if n_p else math.nan))
    return rows


def _max_in(
    terms: Iterable[Terminator],
    contig: str,
    strand: str,
    lo: int,
    hi: int,
    replicon: Replicon,
    yes: float | None,
) -> float | None:
    """Max confidence of primary calls on ``contig``/``strand`` inside ``[lo, hi)``."""
    if lo >= hi:
        return None
    shifts = (0, -replicon.length, replicon.length) if replicon.topology == "circular" else (0,)
    best: float | None = None
    for t in terms:
        iv = t.interval
        if t.source == NOT_PRIMARY or iv.contig != contig or iv.strand != strand:
            continue
        if yes is not None and t.confidence < yes:
            continue
        if any(lo <= iv.start + s and iv.end + s <= hi for s in shifts):
            best = t.confidence if best is None else max(best, t.confidence)
    return best


def _window(up: Gene, down: Gene, replicon: Replicon, slack_bp: int) -> tuple[int, int]:
    """Window between same-strand ``up`` → ``down`` in transcription direction (unwrapped)."""
    u, d = up.interval, down.interval
    circ = replicon.topology == "circular"
    if u.strand == "-":
        s_up = u.start + replicon.length if circ and u.start < d.start else u.start
        return d.end, s_up + slack_bp
    s_down = d.start + replicon.length if circ and d.start < u.start else d.start
    return u.end - slack_bp, s_down


def terminator_between(
    up: Gene,
    down: Gene,
    terms: Sequence[Terminator],
    replicon: Replicon,
    slack_bp: int = 50,
    yes: float = 0.5,
) -> float | None:
    """Max confidence ≥ ``yes`` of a primary call between two same-strand adjacent genes.

    Raises ValueError unless ``up`` and ``down`` share contig and strand (M2-15 review item 6)."""
    u, d = up.interval, down.interval
    if u.contig != d.contig or u.strand != d.strand:
        raise ValueError(
            f"terminator_between: {up.locus_tag} ({u.contig} {u.strand}) and {down.locus_tag} "
            f"({d.contig} {d.strand}) must share contig and strand"
        )
    lo, hi = _window(up, down, replicon, slack_bp)
    iv = up.interval
    return _max_in(terms, iv.contig, iv.strand, lo, hi, replicon, yes)


def terminator_conf_downstream(
    gene: Gene,
    next_gene: Gene | None,
    terms: Sequence[Terminator],
    replicon: Replicon,
    slack_bp: int = 50,
) -> float | None:
    """Max confidence of primary calls after ``gene`` up to the next same-strand gene start, or
    300 bp past the stop when there is none (Q3); no threshold."""
    iv = gene.interval
    if next_gene is not None:
        lo, hi = _window(gene, next_gene, replicon, slack_bp)
    elif iv.strand == "-":
        lo, hi = iv.start - DOWNSTREAM_NO_NEXT_BP, iv.start + slack_bp
    else:
        lo, hi = iv.end - slack_bp, iv.end + DOWNSTREAM_NO_NEXT_BP
    return _max_in(terms, iv.contig, iv.strand, lo, hi, replicon, None)


def terminators_bed(
    terms: Sequence[Terminator], replicons: Mapping[str, Replicon]
) -> list[BedRecord]:
    """BED6 records named ``<source>_<n>`` (numbered per source in genome order), score
    ``round(1000 * confidence)``; origin-crossing intervals stay unwrapped."""
    rank = {c: k for k, c in enumerate(replicons)}

    def key(t: Terminator) -> tuple[int, int, int, str, str, float]:
        iv = t.interval
        return (rank[iv.contig], iv.start, iv.end, iv.strand, t.source, t.confidence)

    counts: dict[str, int] = {}
    out = []
    for t in sorted(terms, key=key):
        counts[t.source] = counts.get(t.source, 0) + 1
        score = min(1000, max(0, round(1000 * t.confidence)))
        out.append(BedRecord(t.interval, f"{t.source}_{counts[t.source]}", score))
    return out
