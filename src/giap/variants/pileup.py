"""Pileup, calls and confidence at every covered codon (GIAP_SPEC.md §7.5.2 steps 2-4; M7-02).

Counting is over distinct peptide sequences after I/L folding (decision 1). A codon is keyed
``(contig, strand, codon)`` with ``codon`` its lowest 0-based coordinate wrapped into
``[0, length)`` on circular replicons. Per codon a folded peptide counts once, in one category:
REF (``match``), IL (``il``, or a substitution between I and L), or one ALT key (a substitution
residue, ``J`` for I or L; or the ``ref>alt`` token of a swap or block at its first codon only).

:func:`pileup` returns one row per covered codon in §6.5 order with the §9.3 columns, followed by
the internal column ``_alts`` (one ``(token, alt_type, alt_aa, n, source, ref_codon, peptides)``
tuple per ALT key) that :func:`sites` reads; writers select the §9.3 columns by name.
Frameshift events are never pileup rows: they travel in ``attrs["frameshifts"]`` (named tuples),
with the contig order in ``attrs["contig_order"]``. :func:`sites` puts the supporting peptides
of each row (a frozenset, aligned with the rows) in ``attrs["alt_peptides"]`` for M7-05 flags.

M7-18: an op whose peptide side holds one of ``NONSTANDARD`` (B J X Z ``ambiguous``, O U ``rare``;
the raw letter, before I/L folding) gives no category at its codon but still spans ``cover``;
the skipped ``(codon, folded peptide)`` pairs are counted in ``attrs["nonstandard"]``.

``ref_codon`` (M7-03) is the gene-strand bases of every codon the REF allele spans: 3 nt for a
substitution or a one-codon block, 6 nt for a swap or a two-codon block (first codon in peptide
order first). A frameshift row has ``alt_aa = fs+1`` or ``fs-1`` (the signed shift) and
``ref_codon`` the codon of the last prefix residue.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Literal, NamedTuple

import pandas as pd

from giap.config.model import GiapConfig
from giap.mapping.loci import in_frame
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel
from giap.variants.align import HitAlignment, ResidueOp, translate_locus

__all__ = [
    "NONSTANDARD",
    "PILEUP_COLUMNS",
    "SITE_COLUMNS",
    "AltType",
    "CallKind",
    "Support",
    "eligible",
    "pileup",
    "sites",
]

CallKind = Literal["REF", "ALT", "mixed", "no_call"]
AltType = Literal["substitution", "swap", "block", "frameshift"]

PILEUP_COLUMNS: tuple[str, ...] = (
    "contig", "codon_start", "strand", "gene", "aa_pos", "ref_aa", "n_ref", "n_il",
    "alt_counts", "call",
)  # fmt: skip
SITE_COLUMNS: tuple[str, ...] = (
    "variant_id", "contig", "codon_start", "strand", "gene", "aa_pos", "ref_codon", "ref_aa",
    "alt_aa", "alt_type", "call", "confidence", "n_ref", "n_alt", "source",
)  # fmt: skip
_INTERNAL = ("_alts",)
_DIFF = frozenset({"substitution", "swap", "block", "near_block"})
_IL = frozenset("IL")
_STRAND_RANK = {"+": 0, "-": 1}
NONSTANDARD: frozenset[str] = frozenset("BJOUXZ")  # peptide letters that are never ALT
_AMBIGUOUS = frozenset("BJXZ")  # the rest of NONSTANDARD (O, U) is ``rare``
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

_Key = tuple[str, str, int]  # (contig, strand, wrapped codon)
_Cat = tuple[str, str, str, str]  # (category, token, alt_type, alt_aa)
_REF: _Cat = ("REF", "", "", "")
_ILC: _Cat = ("IL", "", "", "")


class _Site(NamedTuple):
    """One ``sites`` row before IDs, call and confidence (frameshifts travel as these)."""

    contig: str
    codon_start: int
    strand: str
    gene: str | None
    aa_pos: int | None
    ref_codon: str
    ref_aa: str
    alt_aa: str
    alt_type: str
    n_ref: int
    n_alt: int
    source: str
    peptides: frozenset[str] = frozenset()  # the ALT's supporting peptides (attrs only)


@dataclass(frozen=True, slots=True)
class Support:
    """Per peptide, joined from C9 ``peptides.norm`` and ``peptide_hits``."""

    peptide: str
    source: Literal["database", "denovo"]
    search_db: Literal["same", "other"]
    placement_q: float
    uniqueness: str


def _fold(peptide: str) -> str:
    return peptide.replace("I", "L")


def eligible(
    aligns: Sequence[HitAlignment], support: Mapping[str, Support], cfg: GiapConfig
) -> list[HitAlignment]:
    """Alignments that may support anything (rule 1), in input order.

    Drops peptides without a support row, with placement q above ``stats.alpha.variants`` or
    ``shared``. Under ``genome_relation = same`` a ``database`` peptide searched against another
    strain keeps only its ``match``/``il`` ops and loses its frameshift event (§7.5.1 row 2).
    """
    alpha = cfg.stats.alpha.variants
    out: list[HitAlignment] = []
    for a in aligns:
        s = support.get(a.peptide)
        if s is None or s.placement_q > alpha or s.uniqueness == "shared":
            continue
        if cfg.genome_relation == "same" and s.source == "database" and s.search_db == "other":
            a = replace(a, ops=tuple(o for o in a.ops if o.kind not in _DIFF), frameshift=None)
        out.append(a)
    return out


def _category(op: ResidueOp) -> _Cat:
    if op.kind == "match":
        return _REF
    if op.kind == "il" or (op.kind == "substitution" and op.ref in _IL and op.alt in _IL):
        return _ILC
    if op.kind == "substitution":
        alt = "J" if op.alt in _IL else op.alt
        return ("ALT", alt, "substitution", alt)
    kind = "swap" if op.kind == "swap" else "block"
    return ("ALT", f"{op.ref}>{op.alt}", kind, op.alt)


def _skipped(op: ResidueOp) -> bool:
    return not NONSTANDARD.isdisjoint(op.alt)


def _nonstandard(skipped: Mapping[tuple[_Key, str], str], examples: set[str]) -> dict[str, object]:
    """``attrs["nonstandard"]`` (rule 6): JSON-native, keys sorted."""
    letters = list(skipped.values())
    rare = sum(1 for x in letters if _AMBIGUOUS.isdisjoint(x))
    return {
        "by_residue": {x: letters.count(x) for x in sorted(set(letters))},
        "by_rule": {"ambiguous": len(letters) - rare, "rare": rare},
        "examples": sorted(examples)[:10],
        "skipped_codons": len({k for k, _ in skipped}),
        "skipped_peptides": len({f for _, f in skipped}),
        "skipped_sites": len(letters),
    }


def _wrap(codon: int, rep: Replicon) -> int:
    return codon % rep.length if rep.topology == "circular" else codon


def _codon_nt(seq: str, rep: Replicon, codon: int, strand: str) -> str:
    n = rep.length
    nt = seq[codon : codon + 3] if codon + 3 <= n else seq[codon:] + seq[: codon + 3 - n]
    return nt.translate(_COMP)[::-1] if strand == "-" else nt


def _ref_nt(seq: str, rep: Replicon, key: _Key, token: str) -> str:
    """Gene-strand bases of the codons the REF allele of ``token`` spans (1 or 2 codons)."""
    _, strand, codon = key
    nt = _codon_nt(seq, rep, codon, strand)
    if len(token.partition(">")[0]) == 2:
        nt += _codon_nt(seq, rep, _wrap(codon + (3 if strand == "+" else -3), rep), strand)
    return nt


def _gene_at(gm: GeneModel, rep: Replicon, codon: int, strand: str) -> tuple[str, int] | None:
    """Smallest locus tag (with its 1-based ``aa_pos``) of a same-strand gene framing the codon."""
    shifts = (0, rep.length, -rep.length) if rep.topology == "circular" else (0,)
    best: tuple[str, int] | None = None
    for g in gm.overlapping(Interval(rep.id, codon, codon + 3, strand)):  # type: ignore[arg-type]
        if not isinstance(g, Gene) or g.interval.strand != strand:
            continue
        iv = g.interval
        for d in shifts:
            c = codon + d
            if iv.start <= c and c + 3 <= iv.end and in_frame(c, c + 3, iv):
                pos = (c - iv.start) // 3 + 1 if strand == "+" else (iv.end - c - 3) // 3 + 1
                if best is None or g.locus_tag < best[0]:
                    best = (g.locus_tag, pos)
                break
    return best


def _call(n_ref: int, n_alt: int) -> CallKind:
    if n_alt:
        return "mixed" if n_ref else "ALT"
    return "REF" if n_ref else "no_call"


def _sources(peps: set[str], support: Mapping[str, Support]) -> str:
    return ";".join(sorted({support[p].source for p in peps if p in support}))


def _locate(
    key: _Key, gm: GeneModel, genome: Mapping[str, str]
) -> tuple[str | None, int | None, str, str]:
    """``(gene, aa_pos, ref_codon, ref_aa)`` of a codon key."""
    contig, strand, codon = key
    rep, seq = gm.replicons[contig], genome[contig]
    hit = _gene_at(gm, rep, codon, strand)
    ref_aa = translate_locus(seq, rep, codon, codon + 3, strand)
    gene, aa_pos = (None, None) if hit is None else hit
    return gene, aa_pos, _codon_nt(seq, rep, codon, strand), ref_aa


def pileup(
    aligns: Sequence[HitAlignment],
    gm: GeneModel,
    genome: Mapping[str, str],
    cfg: GiapConfig,
    *,
    support: Mapping[str, Support],
) -> pd.DataFrame:
    """One row per covered codon of the eligible alignments ``aligns`` (rules 2-6, 8-10)."""
    reps = gm.replicons
    cats: dict[_Key, dict[str, set[_Cat]]] = {}
    peps: dict[tuple[_Key, str], set[str]] = {}  # (key, folded) -> original peptides
    cover: dict[_Key, set[str]] = {}
    events: dict[tuple[_Key, int], set[str]] = {}
    skipped: dict[tuple[_Key, str], str] = {}  # (key, folded) -> the skipped op's peptide side
    examples: set[str] = set()
    for a in aligns:
        if a.contig not in reps:
            raise ValueError(f"peptide {a.peptide}: unknown contig {a.contig!r}")
        rep, f = reps[a.contig], _fold(a.peptide)
        for op in a.ops:
            ks = [(a.contig, a.strand, _wrap(c, rep)) for c in op.codons]
            for k in ks:
                cover.setdefault(k, set()).add(f)
            if _skipped(op):
                skipped.setdefault((ks[0], f), op.alt)
                examples.add(a.peptide)
                continue
            cats.setdefault(ks[0], {}).setdefault(f, set()).add(_category(op))
            peps.setdefault((ks[0], f), set()).add(a.peptide)
        if a.frameshift is not None:
            ev = (a.contig, a.strand, _wrap(a.frameshift.last_prefix_codon, rep))
            events.setdefault((ev, a.frameshift.shift), set()).add(a.peptide)
    order = {c: i for i, c in enumerate(reps)}

    def sort_key(k: _Key) -> tuple[int, int, int]:
        return (order[k[0]], k[2], _STRAND_RANK[k[1]])

    rows: list[dict[str, object]] = []
    for key in sorted(cats, key=sort_key):
        n_ref = n_il = 0
        alts: dict[tuple[str, str, str], set[str]] = {}
        for f, c in sorted(cats[key].items()):
            if c == {_REF}:
                n_ref += 1
            elif c <= {_REF, _ILC}:
                n_il += 1
            else:
                _, tok, kind, alt = min(x for x in c if x[0] == "ALT")
                alts.setdefault((tok, kind, alt), set()).update(peps[(key, f)])
        counts = {k: len({_fold(p) for p in v}) for k, v in alts.items()}
        ranked = sorted(alts, key=lambda k: (-counts[k], k))
        gene, aa_pos, _, ref_aa = _locate(key, gm, genome)
        rows.append({
            "contig": key[0], "codon_start": key[2], "strand": key[1], "gene": gene,
            "aa_pos": aa_pos, "ref_aa": ref_aa, "n_ref": n_ref, "n_il": n_il,
            "alt_counts": ";".join(f"{k[0]}:{counts[k]}" for k in ranked) or None,
            "call": _call(n_ref, len(alts)),
            "_alts": tuple((*k, counts[k], _sources(alts[k], support),
                            _ref_nt(genome[key[0]], reps[key[0]], key, k[0]), frozenset(alts[k]))
                           for k in ranked),
        })  # fmt: skip
    fs: list[_Site] = []
    for (key, shift), ps in sorted(events.items(), key=lambda e: (sort_key(e[0][0]), e[0][1])):
        alt_f = {_fold(p) for p in ps}
        contig, strand, codon = key
        nxt = (contig, strand, _wrap(codon + (3 if strand == "+" else -3), reps[contig]))
        n_ref = len((cover.get(key, set()) & cover.get(nxt, set())) - alt_f)
        gene, aa_pos, ref_codon, ref_aa = _locate(key, gm, genome)
        fs.append(_Site(contig, codon, strand, gene, aa_pos, ref_codon, ref_aa, f"fs{shift:+d}",
                        "frameshift", n_ref, len(alt_f), _sources(ps, support),
                        frozenset(ps)))  # fmt: skip
    pile = _frame(rows, [*PILEUP_COLUMNS, *_INTERNAL])
    pile.attrs["frameshifts"] = tuple(fs)
    pile.attrs["contig_order"] = tuple(reps)
    pile.attrs["nonstandard"] = _nonstandard(skipped, examples)
    return pile


def _frame(rows: list[dict[str, object]], columns: Sequence[str]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=list(columns))
    for col in ("codon_start", "n_ref", "n_il", "n_alt"):
        if col in df:
            df[col] = df[col].astype("int64")
    if "aa_pos" in df:
        df["aa_pos"] = df["aa_pos"].astype("Int64")
    for col in ("gene", "alt_counts"):
        if col in df:
            df[col] = df[col].astype(object).where(df[col].notna(), pd.NA)
    return df


def sites(pile: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame:
    """The core rows of ``variants.tsv`` (rules 7 and 11): one per (codon, ALT key) of an
    ``ALT``/``mixed`` codon plus one per frameshift event, in §6.5 order with ``VAR_<n>`` IDs."""
    k = cfg.variants.high_min_peptides
    if "frameshifts" not in pile.attrs or "contig_order" not in pile.attrs:
        # a frame that lost its attrs (rebuilt, round-tripped) would drop every frameshift site
        raise ValueError("sites() needs the frame returned by pileup(); its attrs are missing")
    recs: list[_Site] = list(pile.attrs["frameshifts"])
    for r in pile.to_dict("records"):
        if r["call"] not in ("ALT", "mixed"):
            continue
        gene = None if pd.isna(r["gene"]) else str(r["gene"])
        aa_pos = None if pd.isna(r["aa_pos"]) else int(r["aa_pos"])
        for _tok, kind, alt, n, src, ref_nt, peps in r["_alts"]:
            recs.append(_Site(str(r["contig"]), int(r["codon_start"]), str(r["strand"]), gene,
                              aa_pos, ref_nt, str(r["ref_aa"]), alt, kind,
                              int(r["n_ref"]), n, src, peps))  # fmt: skip
    order = {c: i for i, c in enumerate(pile.attrs["contig_order"])}
    for c in pile["contig"]:
        order.setdefault(c, len(order))
    recs.sort(key=lambda x: (order[x.contig], x.codon_start, _STRAND_RANK[x.strand], x.alt_type,
                             x.alt_aa))  # fmt: skip
    rows: list[dict[str, object]] = [
        {"variant_id": f"VAR_{i}", **x._asdict(), "call": _call(x.n_ref, x.n_alt),
         "confidence": "high" if x.n_alt >= k else "candidate"}
        for i, x in enumerate(recs, 1)
    ]  # fmt: skip
    out = _frame(rows, SITE_COLUMNS)
    out.attrs["alt_peptides"] = tuple(x.peptides for x in recs)
    return out
