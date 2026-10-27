"""Genome preflight rules T-PRE-01..10 (GIPPA_SPEC.md §7.0, §2.1 decision 12, §7.1.2).

Pure functions over parsed ``GenomeRecord`` objects and the resolved config; no file I/O.
Summaries hold 0-based half-open positions; ``Check.value`` reports 1-based closed positions.
"""

from __future__ import annotations

import re
import statistics
from collections import Counter
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from itertools import islice
from typing import Literal

from gippa.config import GippaConfig
from gippa.io.fasta import IUPAC_DNA, GenomeRecord
from gippa.model import Replicon, RepliconKind, Topology
from gippa.preflight.checks import LIST_CAP, Check, JsonValue, Level, capped

__all__ = [
    "LIST_CAP",
    "GenomeSummary",
    "classify_replicons",
    "genome_checks",
    "n50",
    "resolve_topology",
    "summarize_genome",
]

_AMBIGUOUS = "".join(sorted(IUPAC_DNA - set("ACGTN")))
_AMBIGUOUS_RUN_RE = re.compile(f"([{_AMBIGUOUS}])\\1*")
_DELETE_IUPAC = str.maketrans("", "", "".join(IUPAC_DNA))


@dataclass(frozen=True, slots=True)
class GenomeSummary:
    """Genome facts; ``n_runs``/``ambiguous`` hold at most ``LIST_CAP`` runs each.

    ``n_runs`` entries are ``(contig, start0, end0)`` and ``ambiguous`` entries are
    ``(contig, start0, end0, code)`` (0-based half-open, adjacent identical codes merged);
    ``n_run_total``/``ambiguous_total`` count every run, listed or not.
    """

    ids: tuple[str, ...]
    original_ids: tuple[str, ...]
    lengths: tuple[int, ...]
    circular_tags: tuple[bool, ...]
    gc: int
    acgt: int
    n_runs: tuple[tuple[str, int, int], ...]
    ambiguous: tuple[tuple[str, int, int, str], ...]
    invalid: Mapping[str, int]
    empty_records: tuple[str, ...]
    n_run_min_bp: int
    n_run_total: int
    ambiguous_total: int


def _first_runs(rx: re.Pattern[str], seq: str, room: int) -> tuple[Iterator[re.Match[str]], int]:
    """Up to ``room`` leading matches plus the total match count (counted in C via ``subn``)."""
    total = rx.subn("", seq)[1]
    return islice(rx.finditer(seq), max(room, 0)), total


def summarize_genome(records: Sequence[GenomeRecord], n_run_min_bp: int) -> GenomeSummary:
    """Collect the facts the genome rules need from records in FASTA order."""
    if n_run_min_bp < 1:
        raise ValueError(f"n_run_min_bp must be >= 1, got {n_run_min_bp}")
    n_run_re = re.compile(f"N{{{n_run_min_bp},}}")
    gc = acgt = n_total = amb_total = 0
    n_runs: list[tuple[str, int, int]] = []
    ambiguous: list[tuple[str, int, int, str]] = []
    invalid: Counter[str] = Counter()
    empty: list[str] = []
    for r in records:
        seq = r.sequence
        if not seq:
            empty.append(r.id)
        g_c = seq.count("G") + seq.count("C")
        gc += g_c
        acgt += g_c + seq.count("A") + seq.count("T")
        matches, count = _first_runs(n_run_re, seq, LIST_CAP - len(n_runs))
        n_runs.extend((r.id, m.start(), m.end()) for m in matches)
        n_total += count
        matches, count = _first_runs(_AMBIGUOUS_RUN_RE, seq, LIST_CAP - len(ambiguous))
        ambiguous.extend((r.id, m.start(), m.end(), m.group(1)) for m in matches)
        amb_total += count
        invalid.update(seq.translate(_DELETE_IUPAC))
    return GenomeSummary(
        ids=tuple(r.id for r in records),
        original_ids=tuple(r.original_id for r in records),
        lengths=tuple(len(r.sequence) for r in records),
        circular_tags=tuple(r.circular_tag for r in records),
        gc=gc,
        acgt=acgt,
        n_runs=tuple(n_runs),
        ambiguous=tuple(ambiguous),
        invalid=dict(invalid),
        empty_records=tuple(empty),
        n_run_min_bp=n_run_min_bp,
        n_run_total=n_total,
        ambiguous_total=amb_total,
    )


def n50(lengths: Sequence[int]) -> int:
    """Length at which the descending cumulative sum first reaches half the total."""
    if not lengths:
        raise ValueError("n50 of an empty sequence of lengths")
    ordered = sorted(lengths, reverse=True)
    total = sum(ordered)
    cum = 0
    for length in ordered:
        cum += length
        if 2 * cum >= total:
            return length
    return ordered[-1]  # pragma: no cover - unreachable for non-negative lengths


def resolve_topology(s: GenomeSummary, complete: bool | Literal["auto"]) -> tuple[Topology, ...]:
    """Circular when tagged, or every contig when ``complete`` is ``True`` (§5.3)."""
    if complete is True:
        return tuple("circular" for _ in s.ids)
    return tuple("circular" if tag else "linear" for tag in s.circular_tags)


def classify_replicons(
    s: GenomeSummary,
    topo: Sequence[Topology],
    small_bp: int,
    chromosome_min_bp: int = 1_000_000,
) -> list[Replicon]:
    """§7.1.2 replicon kinds; empty records are skipped, ties for largest go to the first."""
    largest = max(range(len(s.lengths)), key=lambda i: s.lengths[i], default=-1)
    out: list[Replicon] = []
    for i, (rid, length, t) in enumerate(zip(s.ids, s.lengths, topo, strict=True)):
        if length < 1:
            continue
        kind: RepliconKind = "unknown"
        if i == largest and (t == "circular" or length >= chromosome_min_bp):
            kind = "chromosome"
        elif t == "circular" and length < small_bp:
            kind = "plasmid"
        out.append(Replicon(id=rid, length=length, topology=t, kind=kind))
    return out


def _median(lengths: Sequence[int]) -> int | float:
    m = statistics.median(lengths)
    return int(m) if float(m).is_integer() else float(m)


def _check(rule_id: str, level: Level, rule: str, value: JsonValue, message: str = "") -> Check:
    return Check(
        rule_id=rule_id, input="genome", level=level, rule=rule, value=value, message=message
    )


def _parse_check(s: GenomeSummary, parse_error: str | None) -> Check:
    rule = "≥ 1 record, IUPAC nucleotides only, no empty sequences → else error"
    if parse_error is not None:
        return _check("genome.parse", "error", rule, parse_error, parse_error)
    if not s.ids:
        return _check("genome.parse", "error", rule, 0, "no FASTA records")
    if s.invalid:
        letters = dict(sorted(s.invalid.items()))
        msg = "non-IUPAC characters: " + ", ".join(f"{k!r} x{v}" for k, v in letters.items())
        return _check("genome.parse", "error", rule, dict(letters), msg)
    if s.empty_records:
        empty: list[JsonValue] = list(s.empty_records)
        return _check("genome.parse", "error", rule, empty, "records with an empty sequence")
    return _check("genome.parse", "ok", rule, len(s.ids))


def genome_checks(s: GenomeSummary, cfg: GippaConfig, parse_error: str | None = None) -> list[Check]:
    """Evaluate every genome rule; only ``genome.parse`` when there are no records."""
    p = cfg.preflight
    out = [_parse_check(s, parse_error)]
    if not s.ids:
        return out

    by_safe: dict[str, list[str]] = {}
    for safe, orig in zip(s.ids, s.original_ids, strict=True):
        by_safe.setdefault(safe, []).append(orig)
    dups: list[JsonValue] = [
        {"id": safe, "original_ids": list(origs)}
        for safe, origs in by_safe.items()
        if len(origs) > 1
    ]
    out.append(
        _check(
            "genome.ids_unique",
            "error" if dups else "ok",
            "first header tokens unique after safe-ID mapping → else error",
            dups,
            "duplicate contig IDs" if dups else "",
        )
    )

    renamed: list[JsonValue] = [
        {"from": orig, "to": safe}
        for safe, orig in zip(s.ids, s.original_ids, strict=True)
        if safe != orig
    ]
    renamed_v, msg = capped(renamed)
    out.append(
        _check(
            "genome.ids_renamed",
            "warn" if renamed else "ok",
            "IDs with characters outside [A-Za-z0-9._-] mapped to '_' → warn",
            renamed_v,
            msg,
        )
    )

    count = len(s.ids)
    out.append(
        _check(
            "genome.contig_count",
            "warn" if count > p.max_contigs_warn else "ok",
            f"contig count > {p.max_contigs_warn} → warn",
            count,
        )
    )

    total = sum(s.lengths)
    out.append(
        _check(
            "genome.total_length",
            "warn" if total < p.min_total_bp_warn else "ok",
            f"total length < {p.min_total_bp_warn} bp → warn",
            total,
        )
    )

    nonempty = [x for x in s.lengths if x > 0]
    if nonempty:
        n50_v, largest = n50(nonempty), max(nonempty)
        bad = n50_v < p.min_n50_warn or largest < p.min_largest_contig_warn
        contiguity: JsonValue = {"n50": n50_v, "largest": largest, "median": _median(nonempty)}
    else:
        bad, contiguity = True, {"n50": 0, "largest": 0, "median": 0}
    out.append(
        _check(
            "genome.contiguity",
            "warn" if bad else "ok",
            f"N50 < {p.min_n50_warn} bp or largest contig < {p.min_largest_contig_warn} bp → warn;"
            " median recorded",
            contiguity,
        )
    )

    gaps: list[JsonValue] = [{"contig": c, "start": a + 1, "end": b} for c, a, b in s.n_runs]
    gaps.extend(
        {"contig": c, "start": a + 1, "end": b, "code": code} for c, a, b, code in s.ambiguous
    )
    msg = ""
    if s.n_run_total > len(s.n_runs) or s.ambiguous_total > len(s.ambiguous):
        msg = (
            f"{s.n_run_total} N runs and {s.ambiguous_total} ambiguity-code runs in total;"
            f" first {LIST_CAP} of each listed"
        )
    out.append(
        _check(
            "genome.gaps",
            "warn" if gaps else "ok",
            f"runs of ≥ {s.n_run_min_bp} N or other IUPAC ambiguity codes → warn",
            gaps,
            msg,
        )
    )

    lo, hi = p.gc_range_warn
    gc_rule = f"GC outside {lo:g}-{hi:g} → warn"
    if s.acgt:
        gc = s.gc / s.acgt
        out.append(_check("genome.gc", "ok" if lo <= gc <= hi else "warn", gc_rule, gc))
    else:
        out.append(_check("genome.gc", "warn", gc_rule, None, "no A/C/G/T bases"))

    topo = resolve_topology(s, cfg.annotation.bakta.complete)
    reps = classify_replicons(
        s, topo, small_bp=p.small_replicon_bp, chromosome_min_bp=p.min_largest_contig_warn
    )
    table: list[JsonValue] = [
        {"id": r.id, "length": r.length, "topology": r.topology, "kind": r.kind} for r in reps
    ]
    table_v, msg = capped(table)
    out.append(
        _check(
            "genome.topology",
            "info",
            "circular=true or [topology=circular] (or annotation.bakta.complete) → circular;"
            " Bakta replicon table",
            table_v,
            msg,
        )
    )

    small: list[JsonValue] = [
        r.id for r in reps if r.length < p.small_replicon_bp and r.kind != "chromosome"
    ]
    small_v, msg = capped(small)
    out.append(
        _check(
            "genome.small_replicons",
            "info" if small else "ok",
            f"contigs < {p.small_replicon_bp} bp listed as putative plasmids or fragments → info",
            small_v,
            msg,
        )
    )
    return out
