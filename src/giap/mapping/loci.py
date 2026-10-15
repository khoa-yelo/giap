"""Locus classes, uniqueness and gene assignment of placements (issue M3-14).

GIAP_SPEC.md §7.2.6 (locus class, uniqueness, peptides -> genes, proteins -> genes, consistency
check), §7.2.7 (genome relation labels) and §9.4 (placement and gene flags), following
docs/specs/M3-14-loci-and-assignment.md and its lead decisions:

* classes in priority order ``rna_overlap`` > ``pseudogene`` > ``sorf`` > ``cds_inframe`` >
  ``nterm_extension`` > ``cterm_extension`` > ``cds_outframe`` > ``antisense`` > ``intergenic``;
* ``cterm_extension`` (ADR-0013): same strand, in frame with a CDS/sORF, overlapping it and
  running past its annotated stop (3' end), also across the origin; assigned to that gene;
  ``cds_outframe`` only when no in-frame relation exists;
* ``multi_identical`` only when every gene has the same non-empty ``aa``;
* an unknown contig raises ``ValueError`` naming it;
* all comparisons on unwrapped intervals: on a circular replicon the hit is shifted by +-L when it
  and the feature lie on opposite sides of the origin;
* ``multi`` is a placement flag on every row of a peptide whose loci have different classes;
* several loci inside one gene are ``unique``; two unassigned loci are ``shared``;
* pseudogenes are assigned in any frame on the same strand;
* output rows sorted by (peptide, contig, start, strand); flags ``;``-joined alphabetically;
* ``novel_candidate_id`` stays NA (novel-locus clustering is Module 6).

Frameshift hits: the frame comes from the prefix (N-terminal) block, whose 5' end is the span's
5' end on either strand (``start`` on plus, ``end`` on minus), so the span gives the frame.

Pure functions; no file I/O.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, get_args

import pandas as pd

from giap.config.model import InputConfig
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel, NcFeature

__all__ = [
    "LOCUS_CLASSES",
    "LocusClass",
    "Uniqueness",
    "annotate_peptide_hits",
    "annotate_protein_hits",
    "classify_locus",
]

LocusClass = Literal[
    "rna_overlap",
    "pseudogene",
    "sorf",
    "cds_inframe",
    "nterm_extension",
    "cterm_extension",
    "cds_outframe",
    "antisense",
    "intergenic",
]
Uniqueness = Literal["unique", "multi_identical", "shared"]
GenomeRelation = Literal["related", "same"]

LOCUS_CLASSES: tuple[str, ...] = get_args(LocusClass)
_RNA_TYPES = frozenset({"rrna", "trna", "tmrna", "ncrna"})
_PLUS_STOPS = frozenset({"TAA", "TAG", "TGA"})  # table 11
_MINUS_STOPS = frozenset({"TTA", "CTA", "TCA"})  # reverse complements, read on the plus string
_ASSIGNING = frozenset({"cds_inframe", "sorf", "nterm_extension", "cterm_extension", "pseudogene"})
FLAG_SEP = ";"
_DEFAULT_DECOY_PREFIXES: tuple[str, ...] = tuple(InputConfig().decoy_prefixes)
_DEFAULT_CONTAMINANT_PREFIXES: tuple[str, ...] = tuple(InputConfig().contaminant_prefixes)


# ---------------------------------------------------------------- interval helpers


def _shifts(rep: Replicon) -> tuple[int, ...]:
    return (0, rep.length, -rep.length) if rep.topology == "circular" else (0,)


def _cterm(s: int, e: int, g: Interval, shifts: tuple[int, ...]) -> bool:
    """True if ``[s, e)`` (shifted) overlaps ``g`` in frame and runs past its 3' end."""
    fs, fe = g.start, g.end
    for d in shifts:
        a, b = s + d, e + d
        if max(a, fs) >= min(b, fe) or not _in_frame(a, b, g):
            continue
        if (g.strand == "+" and b > fe) or (g.strand == "-" and a < fs):
            return True
    return False


def _inside_shift(s: int, e: int, fs: int, fe: int, shifts: tuple[int, ...]) -> int | None:
    for d in shifts:
        if fs <= s + d and e + d <= fe:
            return d
    return None


def _in_frame(s: int, e: int, g: Interval) -> bool:
    if g.strand == "+":
        return (s - g.start) % 3 == 0
    return (g.end - e) % 3 == 0


# ---------------------------------------------------------------- N-terminal extension index


@dataclass(frozen=True, slots=True)
class _Upstream:
    """Same-strand CDSs of one (contig, strand), sorted by their 5' end (plus: start mod L;
    minus: end mod L on circular replicons)."""

    keys: list[int]
    genes: list[Gene]


@dataclass(frozen=True, slots=True)
class _Index:
    upstream: dict[tuple[str, str], _Upstream]


def _build_index(gm: GeneModel) -> _Index:
    groups: dict[tuple[str, str], list[tuple[int, str, Gene]]] = {}
    for g in gm.genes:
        if g.feature_type != "cds" or g.interval.strand not in ("+", "-"):
            continue
        rep = gm.replicons[g.interval.contig]
        circ = rep.topology == "circular"
        key = g.interval.start if g.interval.strand == "+" else g.interval.end
        if circ:
            key %= rep.length
        groups.setdefault((g.interval.contig, g.interval.strand), []).append((key, g.locus_tag, g))
    up = {}
    for k, items in groups.items():
        items.sort(key=lambda t: (t[0], t[1]))
        up[k] = _Upstream([t[0] for t in items], [t[2] for t in items])
    return _Index(up)


def _codon(seq: str, p: int, n: int) -> str:
    q = p % n
    if q + 3 <= n:
        return seq[q : q + 3]
    return (seq[q:] + seq[: q + 3 - n])[:3]


def _nterm_plus(s: int, e: int, rep: Replicon, seq: str, up: _Upstream) -> str | None:
    n = rep.length
    circ = rep.topology == "circular"
    keys, genes = up.keys, up.genes
    m = len(keys)
    if m == 0:
        return None
    # candidates in increasing distance D = gs - s (0 < D < n)
    i0 = bisect_right(keys, s % n if circ else s)
    order: Iterable[int] = (
        (i0 + j) % m for j in range(m)
    ) if circ else range(i0, m)  # fmt: skip
    nxt = s + 3 * -(-(e - s) // 3)  # first in-phase codon at or after the hit end
    blocked_at: int | None = None
    best: tuple[int, str] | None = None
    for i in order:
        g = genes[i]
        gs, ge = g.interval.start, g.interval.end
        d = ((gs - s - 1) % n) + 1 if circ else gs - s
        if d <= 0 or d >= n:
            continue
        if best is not None and d > best[0]:
            break
        target = s + d
        # scan in-phase codons below target for a stop
        while blocked_at is None and nxt < target:
            if not circ and nxt + 3 > n:
                break
            q = nxt % n
            cod = seq[q : q + 3] if q + 3 <= n else _codon(seq, q, n)
            if cod in _PLUS_STOPS:
                blocked_at = nxt
                break
            nxt += 3
        if blocked_at is not None and blocked_at < target:
            break
        if d % 3 != 0 or e > target + (ge - gs):
            continue
        if best is None or (d, g.locus_tag) < best:
            best = (d, g.locus_tag)
    return None if best is None else best[1]


def _nterm_minus(s: int, e: int, rep: Replicon, seq: str, up: _Upstream) -> str | None:
    n = rep.length
    circ = rep.topology == "circular"
    keys, genes = up.keys, up.genes
    m = len(keys)
    if m == 0:
        return None
    # candidates in increasing distance D = e - ge (0 < D < n): descending keys below e
    i0 = bisect_left(keys, e % n if circ else e) - 1
    order: Iterable[int] = (
        (i0 - j) % m for j in range(m)
    ) if circ else range(i0, -1, -1)  # fmt: skip
    nxt = e - 3 * -(-(e - s) // 3)  # end of the first in-phase codon at or below the hit start
    blocked_at: int | None = None
    best: tuple[int, str] | None = None
    for i in order:
        g = genes[i]
        gs, ge = g.interval.start, g.interval.end
        d = ((e - ge - 1) % n) + 1 if circ else e - ge
        if d <= 0 or d >= n:
            continue
        if best is not None and d > best[0]:
            break
        target = e - d
        while blocked_at is None and nxt > target:
            if not circ and nxt - 3 < 0:
                break
            q = (nxt - 3) % n
            cod = seq[q : q + 3] if q + 3 <= n else _codon(seq, q, n)
            if cod in _MINUS_STOPS:
                blocked_at = nxt
                break
            nxt -= 3
        if blocked_at is not None and blocked_at > target:
            break
        if d % 3 != 0 or s < target - (ge - gs):
            continue
        if best is None or (d, g.locus_tag) < best:
            best = (d, g.locus_tag)
    return None if best is None else best[1]


# ---------------------------------------------------------------- classification


def _classify(
    contig: str,
    s: int,
    e: int,
    strand: str,
    gm: GeneModel,
    genome: Mapping[str, str],
    index: _Index,
) -> tuple[LocusClass, str | None]:
    rep = gm.replicons.get(contig)
    if rep is None:
        raise ValueError(f"unknown contig {contig!r}: not a replicon of the gene model")
    shifts = _shifts(rep)
    items = gm.overlapping(Interval(contig, s, e, strand))  # type: ignore[arg-type]
    pseudo: list[str] = []
    sorf: list[str] = []
    inframe: list[str] = []
    cterm: list[str] = []
    same_strand = False
    antisense = False
    for x in items:
        if isinstance(x, NcFeature):
            if x.type in _RNA_TYPES:
                return "rna_overlap", None
            continue
        fi = x.interval
        if fi.strand != strand:
            if x.feature_type != "pseudogene" and fi.strand in ("+", "-"):
                antisense = True
            continue
        if x.feature_type == "pseudogene":
            pseudo.append(x.locus_tag)
            continue
        same_strand = True
        d = _inside_shift(s, e, fi.start, fi.end, shifts)
        if d is not None and _in_frame(s + d, e + d, fi):
            (sorf if x.feature_type == "sorf" else inframe).append(x.locus_tag)
        elif _cterm(s, e, fi, shifts):
            cterm.append(x.locus_tag)
    if pseudo:
        return "pseudogene", min(pseudo)
    if sorf:
        return "sorf", min(sorf)
    if inframe:
        return "cds_inframe", min(inframe)
    if strand in ("+", "-"):
        up = index.upstream.get((contig, strand))
        if up is not None:
            seq = genome[contig]
            fn = _nterm_plus if strand == "+" else _nterm_minus
            tag = fn(s, e, rep, seq, up)
            if tag is not None:
                return "nterm_extension", tag
    if cterm:
        return "cterm_extension", min(cterm)
    if same_strand:
        return "cds_outframe", None
    if antisense:
        return "antisense", None
    return "intergenic", None


def classify_locus(
    iv: Interval, gm: GeneModel, genome: Mapping[str, str]
) -> tuple[LocusClass, str | None]:
    """Locus class of ``iv`` and the gene it assigns to (None unless an assigning class)."""
    return _classify(iv.contig, iv.start, iv.end, iv.strand, gm, genome, _build_index(gm))


# ---------------------------------------------------------------- peptide hits


def _split_flags(x: object) -> set[str]:
    if x is None or (not isinstance(x, str) and pd.isna(x)):
        return set()
    return {f for f in str(x).split(FLAG_SEP) if f}


def _join_flags(flags: Iterable[str]) -> str:
    return FLAG_SEP.join(sorted(set(flags)))


def _engine_ids(
    x: object, decoy_prefixes: Sequence[str], contaminant_prefixes: Sequence[str]
) -> frozenset[str]:
    prefixes = (*decoy_prefixes, *contaminant_prefixes)
    return frozenset(p for p in _split_flags(x) if not p.startswith(prefixes))


def _identical_aa(gm: GeneModel, tags: Sequence[str | None]) -> bool:
    """All genes have the same, non-empty protein sequence."""
    aas = {gm.by_locus(t).aa for t in tags if t is not None}
    return len(aas) == 1 and "" not in aas


_OUT_COLUMNS = ("locus_class", "gene", "uniqueness", "flags", "novel_candidate_id")


def annotate_peptide_hits(
    accepted: pd.DataFrame,
    gm: GeneModel,
    genome: Mapping[str, str],
    peptides: pd.DataFrame,
    genome_relation: GenomeRelation,
    *,
    decoy_prefixes: Sequence[str] = _DEFAULT_DECOY_PREFIXES,
    contaminant_prefixes: Sequence[str] = _DEFAULT_CONTAMINANT_PREFIXES,
) -> pd.DataFrame:
    """Add ``locus_class``, ``gene``, ``uniqueness``, ``flags`` and ``novel_candidate_id``.

    One output row per input row, sorted by (peptide, contig, start, strand); input columns are
    kept and input flags (``same_db``, ``low_complexity``) are merged into ``flags``.
    """
    if genome_relation not in ("related", "same"):
        raise ValueError(f"genome_relation must be 'related' or 'same', got {genome_relation!r}")
    df = accepted.copy()
    n = len(df)
    contigs = df["contig"].astype(str).tolist() if n else []
    starts = [int(v) for v in df["start"].tolist()] if n else []
    ends = [int(v) for v in df["end"].tolist()] if n else []
    strands = df["strand"].astype(str).tolist() if n else []
    peps = df["peptide"].astype(str).tolist() if n else []
    routes = df["route"].tolist() if n else []
    in_flags = df["flags"].tolist() if "flags" in df.columns else [None] * n

    index = _build_index(gm)
    cache: dict[tuple[str, int, int, str], tuple[LocusClass, str | None]] = {}
    classes: list[LocusClass] = []
    genes: list[str | None] = []
    for c, s, e, st in zip(contigs, starts, ends, strands, strict=True):
        key = (c, s, e, st)
        r = cache.get(key)
        if r is None:
            r = _classify(c, s, e, st, gm, genome, index)
            cache[key] = r
        classes.append(r[0])
        genes.append(r[1])

    # per-peptide uniqueness and `multi`
    rows_of: dict[str, list[int]] = {}
    for i, p in enumerate(peps):
        rows_of.setdefault(p, []).append(i)
    uniq: list[Uniqueness] = ["unique"] * n
    multi = [False] * n
    for idx in rows_of.values():
        if len(idx) == 1:
            continue
        units = {
            ("g", genes[i]) if genes[i] is not None else ("l", contigs[i], starts[i], ends[i],
                                                          strands[i])
            for i in idx
        }  # fmt: skip
        u: Uniqueness
        if len(units) == 1:
            u = "unique"
        elif all(genes[i] is not None for i in idx) and _identical_aa(gm, [genes[i] for i in idx]):
            u = "multi_identical"
        else:
            u = "shared"
        is_multi = len({classes[i] for i in idx}) > 1
        for i in idx:
            uniq[i] = u
            multi[i] = is_multi

    # engine protein consistency (same-strain input only)
    same_db: dict[str, frozenset[str]] = {}
    if n and len(peptides):
        for p, sdb, ids in zip(
            peptides["peptide"].astype(str).tolist(),
            peptides["search_db"].tolist(),
            peptides["protein_ids"].tolist(),
            strict=True,
        ):
            if sdb == "same":
                same_db[p] = _engine_ids(ids, decoy_prefixes, contaminant_prefixes)

    related = genome_relation == "related"
    lengths = {c: r.length for c, r in gm.replicons.items()}
    flags_out: list[str] = []
    for i in range(n):
        fl = _split_flags(in_flags[i])
        if multi[i]:
            fl.add("multi")
        if uniq[i] == "multi_identical":
            fl.add("multi_identical")
        if classes[i] == "nterm_extension":
            fl.add("nterm_extension_evidence")
        if ends[i] > lengths[contigs[i]]:
            fl.add("origin_spanning")
        if related and routes[i] == "tolerant":
            fl.add("strain_difference")
        g = genes[i]
        if g is not None and peps[i] in same_db and g not in same_db[peps[i]]:
            fl.add("engine_protein_mismatch")
        flags_out.append(_join_flags(fl))

    df["locus_class"] = pd.Series(classes, index=df.index, dtype=object)
    df["gene"] = pd.Series(genes, index=df.index, dtype=object)
    df["uniqueness"] = pd.Series(uniq, index=df.index, dtype=object)
    df["flags"] = pd.Series(flags_out, index=df.index, dtype=object)
    df["novel_candidate_id"] = pd.Series([None] * n, index=df.index, dtype=object)
    sort_cols = ["peptide", "contig", "start", "strand", "end"]
    extra = [c for c in ("route", "search", "placement_q") if c in df.columns]
    if n:
        df = df.sort_values([*sort_cols, *extra], kind="mergesort")
    return df.reset_index(drop=True)


# ---------------------------------------------------------------- protein rows


def annotate_protein_hits(rows: pd.DataFrame, gm: GeneModel) -> pd.DataFrame:
    """Add gene flags ``cross_strain_protein`` and ``two_proteins_one_gene`` to protein rows.

    ``cross_strain_protein``: ``protein_id`` is not a locus tag of the model;
    ``two_proteins_one_gene``: every row of a gene assigned >= 2 distinct input proteins.
    """
    df = rows.copy()
    n = len(df)
    pids = df["protein_id"].astype(str).tolist() if n else []
    genes = df["gene"].tolist() if n else []
    in_flags = df["flags"].tolist() if "flags" in df.columns else [None] * n
    tags = {g.locus_tag for g in gm.genes}
    per_gene: dict[str, set[str]] = {}
    for p, g in zip(pids, genes, strict=True):
        if g is not None and not pd.isna(g):
            per_gene.setdefault(str(g), set()).add(p)
    out: list[str] = []
    for p, g, f in zip(pids, genes, in_flags, strict=True):
        fl = _split_flags(f)
        if p not in tags:
            fl.add("cross_strain_protein")
        if g is not None and not pd.isna(g) and len(per_gene[str(g)]) >= 2:
            fl.add("two_proteins_one_gene")
        out.append(_join_flags(fl))
    df["flags"] = pd.Series(out, index=df.index, dtype=object)
    return df
