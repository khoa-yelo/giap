"""Novel-ORF outputs: the ``novel_orfs`` table, its GFF3 features, the patch features and the
protein records (issue M7-07; GIAP_SPEC.md §7.6.6, §9.3, §6.5; ADR-0012 decision 3).

Pure. The table holds the §9.3 cells; ``start_longest``, ``start_best_ostir`` and ``end`` are
0-based *positions* in transcription direction (start: the start codon's first base, on ``-``
its highest coordinate; end: the stop codon's last base), so on ``-`` start > end (lead
decision 1). Writers add 1. Intervals for GFF3 and translation are rebuilt from a start position
and an amino-acid length, which keeps origin-crossing ORFs at ``end > length`` (M1-05).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any, Literal

import pandas as pd

from giap.annotate.bakta import interval_nt
from giap.io.gff3 import GffFeature
from giap.io.protein_fasta import LetterCounts, ProteinRecord
from giap.io.tsv import format_value
from giap.mapping.sixframe import translate11_text
from giap.model.coords import Interval, Replicon, check_interval
from giap.model.genes import GeneModel
from giap.novel.candidates import NOVEL_CLASSES, Orf, _il, _stops, _strand_view, _walk_up

Kind = Literal["str", "int", "float", "bool"]
NOVEL_ORFS_TYPES: dict[str, Kind] = {
    "orf_id": "str", "class": "str", "contig": "str", "strand": "str", "frame": "int",
    "start_longest": "int", "start_best_ostir": "int", "end": "int", "length_aa_longest": "int",
    "length_aa_best": "int", "start_codon_best": "str", "n_peptides": "int", "peptides": "str",
    "best_peptide_q": "float", "ostir_tir": "float", "promoter_score": "float",
    "promoter_supported": "bool", "terminator_downstream": "bool", "terminator_conf": "float",
    "terminator_distance_nt": "int", "rbs_supported": "bool", "class_fdr": "str",
    "called": "bool", "rank": "int", "overlapping_gene": "str", "flags": "str",
}  # fmt: skip
NOVEL_ORFS_COLUMNS: tuple[str, ...] = tuple(NOVEL_ORFS_TYPES)
POSITION_COLUMNS = ("start_longest", "start_best_ostir", "end")
PATCH_PRODUCT = "hypothetical protein (GIAP novel ORF)"
NEW_ORF_CLASSES = frozenset({"intergenic", "cds_outframe", "antisense"})
# Second-round and evidence columns copied from ``scored`` when present (NA otherwise).
_SCORED = ("ostir_tir", "promoter_score", "promoter_supported", "terminator_downstream",
           "terminator_conf", "terminator_distance_nt", "rbs_supported", "rank")  # fmt: skip
_GFF_ATTRS = ("rank", "n_peptides")
_ADR12 = ("promoter_supported", "terminator_downstream", "terminator_conf",
          "terminator_distance_nt")  # fmt: skip
_DTYPE = {"str": "object", "int": "Int64", "float": "float64", "bool": "boolean"}
# The novel_orfs ``flags`` vocabulary in schema order: §9.4 "Novel ORFs", then the milestone
# additions (M7-06/07, M7-09, ADR-0012, ADR-0013). Flags are joined in this order (M7-09 d4).
FLAG_ORDER: tuple[str, ...] = (
    "junction_spanned", "homopolymer_junction", "class_fdr_unestimated", "rbs_supported",
    "no_start_codon", "partial", "ostir_missing", "rbs_threshold_unestimated", "junction_unknown",
    "promoter_supported", "terminator_supported", "readthrough_candidate",
    "denovo_only_unsupported",
)  # fmt: skip
_FLAG_RANK = {f: i for i, f in enumerate(FLAG_ORDER)}


def is_na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def order_flags(flags: set[str]) -> str | None:
    """``;``-joined flags in ``FLAG_ORDER`` (unknown flags last, by name); ``None`` if empty."""
    key = sorted(flags, key=lambda f: (_FLAG_RANK.get(f, len(FLAG_ORDER)), f))
    return ";".join(key) or None


def _na_none(v: object) -> object:
    return None if is_na(v) else v


def novel_frame(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    """The typed table (nullable ``Int64``/``boolean``, ``float64``, object strings)."""
    data = {
        c: pd.Series([None if is_na(r[c]) else r[c] for r in rows], dtype=_DTYPE[k])
        for c, k in NOVEL_ORFS_TYPES.items()
    }
    return pd.DataFrame(data, columns=list(NOVEL_ORFS_COLUMNS))


def _unique_ids(ids: Sequence[str]) -> None:
    if len(set(ids)) != len(ids):
        raise ValueError("two novel ORFs share one orf_id")


def _scored_by_id(scored: pd.DataFrame | None) -> dict[str, dict[str, Any]]:
    if scored is None:
        return {}
    recs = scored.to_dict("records")
    _unique_ids([str(r["orf_id"]) for r in recs])
    return {str(r["orf_id"]): r for r in recs}


def _class_fdr(v: Any) -> str:
    return "unestimated" if is_na(v) else v if isinstance(v, str) else format_value(float(v))


def _orf_row(
    o: Orf, s: Mapping[str, Any], genome: Mapping[str, str], rep: Replicon
) -> dict[str, Any]:
    sv = _strand_view(o.contig, o.strand, rep, genome)
    n, plus = sv.n, o.strand == "+"

    def t(x: int) -> int:  # transcription coordinate of the codon at genome x
        return x if plus else n - x - 3

    t_end = t(o.stop) + (3 if o.end_3p_open else 0)  # exclusive end of the coding codons

    def n_aa(x_t: int) -> int:
        d = t_end - x_t
        return (d % n if sv.circular else d) // 3

    best = s.get("start_best_ostir")
    best = None if is_na(best) or not o.starts else int(best)  # type: ignore[arg-type]
    if best is not None and best not in o.starts:
        raise ValueError(f"{o.orf_id}: start_best_ostir {best} is not one of its start codons")
    if o.start_longest is not None:
        len_longest = n_aa(t(o.start_longest))
    elif o.first_peptide is not None:  # from the first peptide's mapped position
        len_longest = n_aa(t(o.first_peptide))
    else:  # no position: find the first peptide (I/L folded) in the translated ORF body
        q, _ = _walk_up(sv, t_end, _stops(11))
        body = "".join(sv.codon(p) or "NNN" for p in range(q + 3, t_end, 3))
        aa = _il(translate11_text(body))
        hits = [i for i in (aa.find(p) for p in o.peptides) if i >= 0]
        if not hits:
            raise ValueError(f"{o.orf_id}: no peptide found in the ORF translation")
        len_longest = len(aa) - min(hits)
    used = best if best is not None else o.start_longest
    fdr = _class_fdr(s.get("class_fdr"))
    flags = set(str(s.get("flags") or "").split(";")) - {"", "class_fdr_unestimated"}
    flags |= {"class_fdr_unestimated"} if fdr == "unestimated" else set()
    flags |= {"no_start_codon"} if not o.starts else set()
    flags |= {"partial"} if o.end_5p_open or o.end_3p_open else set()
    flags |= {"ostir_missing"} if o.starts and best is None else set()

    def pos(x: int | None) -> int | None:  # first base of the codon in transcription direction
        return None if x is None else x if plus else (x + 2) % n

    called = s.get("called")
    row = {
        "orf_id": o.orf_id, "class": o.cls, "contig": o.contig, "strand": o.strand,
        "frame": o.frame, "start_longest": pos(o.start_longest), "start_best_ostir": pos(best),
        "end": (o.stop + 2) % n if plus else o.stop, "length_aa_longest": len_longest,
        "length_aa_best": len_longest if used is None else n_aa(t(used)),
        "start_codon_best": None if used is None else sv.codon(t(used)),
        "n_peptides": s.get("n_peptides", len(o.peptides)),
        "peptides": ";".join(o.peptides) or None, "best_peptide_q": o.best_q,
        "class_fdr": fdr, "called": False if is_na(called) else bool(called),
        "overlapping_gene": o.overlapping_gene, "flags": order_flags(flags),
    }  # fmt: skip
    row.update({c: s.get(c) for c in _SCORED})
    return row


def orf_table(
    orfs: Sequence[Orf],
    scored: pd.DataFrame | None,
    *,
    genome: Mapping[str, str],
    replicons: Mapping[str, Replicon],
) -> pd.DataFrame:
    """The ``novel_orfs`` table in §6.5 order (contig FASTA order, lowest coordinate, end, strand,
    ID). ``scored`` is keyed by ``orf_id`` (M7-06 ``minimal_call`` columns, optionally
    ``start_best_ostir`` as a 0-based codon of ``Orf.starts`` and any later table column)."""
    _unique_ids([o.orf_id for o in orfs])
    sc = _scored_by_id(scored)
    order = {c: i for i, c in enumerate(replicons)}
    rows = [_orf_row(o, sc.get(o.orf_id, {}), genome, replicons[o.contig]) for o in orfs]
    rows.sort(key=lambda r: (order[r["contig"]], _id_pos(r["orf_id"]), r["end"], r["strand"],
                             r["orf_id"]))  # fmt: skip
    return novel_frame(rows)


def _id_pos(orf_id: str) -> int:
    return int(orf_id.rsplit("_", 2)[1])


def _opt(v: object) -> int | None:
    return None if is_na(v) else int(v)  # type: ignore[call-overload]


def _span(r: Mapping[str, Any], start: int | None, n_aa: int) -> Interval:
    """The interval from ``start`` (a position, or the first peptide) over ``n_aa`` codons and
    the stop codon (none when the ORF is open at a linear 3' end)."""
    end = int(r["end"])
    if r["strand"] == "+":
        lo = start if start is not None else _id_pos(r["orf_id"]) - 1
        open3 = end + 1 == lo + 3 * n_aa
    else:
        lo = end
        flags = str(r["flags"] or "").split(";")
        # without a start, only an ORF ending at the contig's first codon can be 3'-open
        open3 = start + 1 == lo + 3 * n_aa if start is not None else "partial" in flags and lo < 3
    return Interval(str(r["contig"]), lo, lo + 3 * n_aa + (0 if open3 else 3), r["strand"])


def _records(table: pd.DataFrame) -> list[dict[str, Any]]:
    return [dict(r) for r in table.to_dict("records")]


def _best(r: Mapping[str, Any]) -> int | None:
    b = _opt(r["start_best_ostir"])
    return b if b is not None else _opt(r["start_longest"])


def gff_features(table: pd.DataFrame) -> list[GffFeature]:
    """One ``CDS`` per row (called or not): best start (or longest, or first peptide) to stop."""
    out = []
    for r in _records(table):
        attrs = {"ID": str(r["orf_id"]), "class": str(r["class"]),
                 "called": format_value(bool(r["called"]))}  # fmt: skip
        attrs.update({c: format_value(_na_none(r[c])) for c in _GFF_ATTRS})
        start = _opt(r["start_longest"])
        attrs["start_longest"] = format_value(None if start is None else start + 1)
        attrs["start_codon"] = format_value(_na_none(r["start_codon_best"]))
        attrs.update({c: format_value(_na_none(r[c])) for c in _ADR12})
        attrs["Note"] = str(r["class"])
        iv = _span(r, _best(r), int(r["length_aa_best"]))
        if "partial" in str(r["flags"] or "").split(";"):
            attrs["partial"] = "true"
        out.append(GffFeature(iv, "GIAP", "CDS", None, 0, attrs))
    return out


def patch_features(table: pd.DataFrame, gm: GeneModel) -> list[GffFeature]:
    """``annotation_patch.gff3``: called new ORFs as ``CDS`` and called N-terminal extensions as
    the gene with the new start (rule 7); no ``cterm_extension`` row until M7-16."""
    rows = _records(table)
    _unique_ids([str(r["orf_id"]) for r in rows])
    out = []
    for r in rows:
        start = _best(r)
        if not r["called"] or start is None:
            continue
        cls, contig, strand = str(r["class"]), str(r["contig"]), str(r["strand"])
        attrs: dict[str, str]
        if cls in NEW_ORF_CLASSES:
            iv = _span(r, start, int(r["length_aa_best"]))
            attrs = {"ID": str(r["orf_id"]), "product": PATCH_PRODUCT}
        elif cls == "nterm_extension":
            tag = r["overlapping_gene"]
            if is_na(tag):
                raise ValueError(f"{r['orf_id']}: nterm_extension without its gene")
            g, n = gm.by_locus(str(tag)).interval, gm.replicons[contig].length
            lo, hi = (start, g.end) if strand == "+" else (g.start, start + 1)
            iv = Interval(contig, lo, hi + n if hi <= lo else hi, strand)  # type: ignore[arg-type]
            attrs = {"ID": str(tag), "locus_tag": str(tag), "Note": "nterm_extension"}
        else:
            continue
        out.append(GffFeature(iv, "GIAP", "CDS", None, 0, attrs))
    return out


def protein_records(
    table: pd.DataFrame, genome: Mapping[str, str], replicons: Mapping[str, Replicon]
) -> list[ProteinRecord]:
    """``|longest`` and ``|best_ostir`` proteins (merged when equal) of every ORF with a start:
    table 11, the initiator read as ``M``, no ``*``."""
    out = []
    for r in _records(table):
        longest, best = _opt(r["start_longest"]), _opt(r["start_best_ostir"])
        if longest is None:
            continue
        oid = str(r["orf_id"])
        variants = [("longest", longest, int(r["length_aa_longest"]))]
        if best == longest:
            variants = [("longest,best_ostir", longest, variants[0][2])]
        elif best is not None:
            variants.append(("best_ostir", best, int(r["length_aa_best"])))
        for suffix, start, n_aa in variants:
            iv = _span(r, start, n_aa)
            check_interval(iv, replicons[iv.contig])
            nt = interval_nt(genome[iv.contig].upper(), iv)
            aa = ("M" + translate11_text(nt[3:])).removesuffix("*")
            header = f"{oid}|{suffix}"
            out.append(ProteinRecord(header, (header,), header, aa, aa.count("*"),
                                     LetterCounts(()), LetterCounts(())))  # fmt: skip
    return out


def _cell_int(v: str | None) -> int | None:
    return int(v) if v is not None and v not in ("", "NA") else None


def _orf_span(o: Mapping[str, str | None]) -> tuple[int, int] | None:
    """The ORF's 1-based closed span in transcription order, (start, end): ``start_longest``
    (``start_best_ostir``, then ``start``, when NA) to ``end``. A ``no_start_codon`` ORF has no
    start column; its first codon is ``length_aa_longest`` codons plus the stop before ``end``.
    ``None`` without an end or a way to place the start."""
    end = _cell_int(o.get("end"))
    if end is None:
        return None
    start = next(
        (s for k in ("start_longest", "start_best_ostir", "start") if (s := _cell_int(o.get(k)))),
        None,
    )
    if start is None and (n_aa := _cell_int(o.get("length_aa_longest"))) is not None:
        nt = 3 * (n_aa + 1) - 1  # the codons and the stop, minus the end position itself
        start = end - nt if o.get("strand") != "-" else end + nt
    return None if start is None else (start, end)


def _within(hit_start: int, hit_end: int, strand: str | None, span: tuple[int, int]) -> bool:
    """Hit ``[hit_start, hit_end]`` (1-based closed; ``hit_end`` past the replicon end when it
    crosses the origin) lies in the ORF span: ``[start, end]`` on ``+``, ``[end, start]`` on
    ``-``; an ORF whose span wraps (``end < start`` on ``+``) covers ``[start, ...]`` and
    ``[..., end]``."""
    lo, hi = span if strand != "-" else (span[1], span[0])
    if lo <= hi:
        return lo <= hit_start and hit_end <= hi
    return hit_start >= lo or hit_end <= hi


def candidate_ids(
    hits: Sequence[Mapping[str, str | None]], orfs: Sequence[Mapping[str, str | None]]
) -> list[str | None]:
    """``novel_candidate_id`` per ``peptide_hits`` row (M7-06 Q6, back-filled in M7-15): the
    ``orf_id`` of the first ORF listing the row's peptide (I = L folded) on its contig and
    strand whose span holds the hit (``[start_longest, end]``, ``[end, start_longest]`` on
    ``-``; the ORF's own start when ``start_longest`` is NA), for rows of a novel locus class;
    else ``None``. Rows and ORFs are published-table cells (1-based closed). Pure.
    """
    owners: dict[tuple[str | None, str | None, str], list[tuple[str, tuple[int, int] | None]]]
    owners = {}
    for o in orfs:
        oid = o["orf_id"]
        if not oid:
            continue
        span = _orf_span(o)
        for pep in dict.fromkeys(_il(p) for p in (o["peptides"] or "").split(";") if p):
            owners.setdefault((o["contig"], o["strand"], pep), []).append((oid, span))

    def owner(h: Mapping[str, str | None]) -> str | None:
        if h["locus_class"] not in NOVEL_CLASSES:
            return None
        start, end = _cell_int(h.get("start")), _cell_int(h.get("end"))
        for oid, span in owners.get((h["contig"], h["strand"], _il(h["peptide"] or "")), []):
            if span is None or start is None or end is None:
                continue
            if _within(start, end, h["strand"], span):
                return oid
        return None

    return [owner(h) for h in hits]
