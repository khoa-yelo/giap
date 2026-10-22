"""Novel-ORF candidates, ORF reconstruction, IDs, second-round Promotech windows and OSTIR rows
(issue M7-06; GIAP_SPEC.md §7.6.1, §7.6.2, §7.6.4 "Neither", §7.1.4, §6.5).

Pure. Walks run in *transcription coordinates*: the plus string for ``+``, its reverse complement
for ``-`` (a minus codon ``[x, x + 3)`` maps to the genome codon ``[L - x - 3, L - x)``).
Circular replicons wrap; linear contig ends open the ORF. Reported codon coordinates are the
lowest 0-based genome coordinate of the codon, wrapped into ``[0, L)``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, cast, get_args

import pandas as pd
from Bio.Data import CodonTable

from giap.annotate.ostir import OstirRow
from giap.annotate.promoters import Window, dedupe_and_chunk, gene_windows
from giap.config.model import GiapConfig
from giap.io.protein_fasta import ProteinRecord
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel

NovelClass = Literal[
    "nterm_extension", "cterm_extension", "intergenic", "cds_outframe", "antisense",
    "pseudogene", "sorf",
]  # fmt: skip
NOVEL_CLASSES: tuple[str, ...] = get_args(NovelClass)
# Rule 5 / decision 3: the M3-14 priority order, highest first.
CLASS_PRIORITY: tuple[str, ...] = (
    "pseudogene", "sorf", "nterm_extension", "cterm_extension", "cds_outframe", "antisense",
    "intergenic",
)  # fmt: skip
EXCLUDING_CLASSES = frozenset({"cds_inframe", "rna_overlap"})
CALL_COLUMNS = ["orf_id", "n_peptides", "called", "class_fdr", "flags"]
_COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


@dataclass(frozen=True, slots=True)
class Orf:
    """One candidate novel ORF (key ``(contig, strand, stop)``)."""

    orf_id: str  # NOVEL_<contig>_<start1>_<p|m>
    cls: NovelClass
    contig: str
    strand: Literal["+", "-"]
    frame: int
    stop: int  # lowest 0-based coordinate of the stop codon (wrapped); last whole codon if open
    starts: tuple[int, ...]  # candidate start codons, 5'->3' (lowest coordinate of each)
    start_longest: int | None
    end_5p_open: bool  # upstream walk hit a linear contig end
    end_3p_open: bool  # downstream walk hit a linear contig end
    peptides: tuple[str, ...]  # distinct, I/L-folded, sorted
    best_q: float | None  # min input q of its peptides (None for de novo)
    overlapping_gene: str | None


def _il(s: str) -> str:
    return s.upper().replace("I", "L")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _stops(table: int) -> frozenset[str]:
    try:
        return frozenset(CodonTable.unambiguous_dna_by_id[table].stop_codons)
    except KeyError:
        raise ValueError(f"unknown translation table {table}") from None


# ---------------------------------------------------------------- rule 1: candidate rows


def candidate_rows(
    peptide_hits: pd.DataFrame,
    peptides: pd.DataFrame,
    contaminants: Sequence[ProteinRecord] | None,
) -> pd.DataFrame:
    """The retained C11 rows (all columns kept) plus ``q_value`` joined from C9 ``peptides``."""
    hits = peptide_hits
    pep = hits["peptide"].astype(str)
    banned = set(pep[hits["locus_class"].isin(EXCLUDING_CLASSES)])
    banned |= set(peptides.loc[peptides["is_contaminant"].astype(bool), "peptide"].astype(str))
    if contaminants:
        library = "#".join(_il(r.sequence) for r in contaminants)
        banned |= {p for p in set(pep) if _il(p) in library}
    keep = hits["locus_class"].isin(NOVEL_CLASSES) & ~pep.isin(banned)
    out = hits.loc[keep].copy()
    loci = out.drop_duplicates(["peptide", "contig", "start", "end", "strand"])
    n_loci = loci["peptide"].value_counts()
    out = out.loc[out["peptide"].map(n_loci).fillna(0).to_numpy() == 1]
    q = dict(zip(peptides["peptide"].astype(str), peptides["q_value"], strict=True))
    out["q_value"] = pd.to_numeric(out["peptide"].map(q), errors="coerce").astype("float64")
    return out.reset_index(drop=True)


# ---------------------------------------------------------------- rules 2-6: reconstruction


@dataclass(frozen=True, slots=True)
class _Strand:
    """A strand-oriented contig sequence for codon walks."""

    contig: str
    strand: str
    tseq: str
    n: int
    circular: bool

    def codon(self, p: int) -> str | None:
        """The codon at transcription coordinate ``p``; ``None`` off a linear contig."""
        n = self.n
        if not self.circular:
            return self.tseq[p : p + 3] if p >= 0 and p + 3 <= n else None
        q = p % n
        return (
            self.tseq[q : q + 3]
            if q + 3 <= n
            else "".join(self.tseq[(q + i) % n] for i in range(3))
        )

    def to_t(self, start: int, end: int) -> tuple[int, int]:
        return (start, end) if self.strand == "+" else (self.n - end, self.n - start)

    def codon_genome(self, t: int) -> int:
        """Lowest genome coordinate (wrapped) of the codon at transcription coordinate ``t``."""
        return t % self.n if self.strand == "+" else (self.n - t - 3) % self.n

    def no_stop(self) -> ValueError:
        return ValueError(f"no stop codon on {self.contig} strand {self.strand} in this frame")


def _strand_view(contig: str, strand: str, rep: Replicon, genome: Mapping[str, str]) -> _Strand:
    seq = genome[contig].upper()
    tseq = seq if strand == "+" else _revcomp(seq)
    return _Strand(contig, strand, tseq, rep.length, rep.topology == "circular")


def _walk_down(sv: _Strand, ts: int, te: int, stops: frozenset[str]) -> tuple[int, bool]:
    """``(end_t, open)``: exclusive end of the stop codon, or of the last whole codon if open."""
    p = ts + 3 * -(-(te - ts) // 3)
    while True:
        cod = sv.codon(p)
        if cod is None:
            return p, True
        if cod in stops:
            return p + 3, False
        p += 3
        if p - te > sv.n + 3:
            raise sv.no_stop()


def _walk_up(sv: _Strand, first: int, stops: frozenset[str]) -> tuple[int, bool]:
    """``(q, open)``: the previous in-frame stop codon (or the first invalid codon if open)."""
    q = first - 3
    while True:
        cod = sv.codon(q)
        if cod is None:
            return q, True
        if cod in stops:
            return q, False
        q -= 3
        if first - q > sv.n + 3:
            raise sv.no_stop()


def _overlapping_gene(cls: str, row: Mapping[str, Any], gm: GeneModel) -> str | None:
    if cls == "intergenic":
        return None
    gene = row.get("gene")
    if isinstance(gene, str) and gene:
        return gene
    if cls != "cds_outframe":
        return None
    iv = Interval(str(row["contig"]), int(row["start"]), int(row["end"]), row["strand"])
    for g in gm.overlapping(iv):
        if isinstance(g, Gene) and g.interval.strand == iv.strand:
            return g.locus_tag
    return None


def _best_q(qs: Sequence[float]) -> float | None:
    vals = [float(q) for q in qs if pd.notna(q)]
    return min(vals) if vals else None


def reconstruct(
    rows: pd.DataFrame, gm: GeneModel, genome: Mapping[str, str], cfg: GiapConfig
) -> list[Orf]:
    """Group candidate rows into ORFs by ``(contig, strand, stop)``; sorted by contig, ID."""
    if len(rows) == 0:
        return []
    stops = _stops(cfg.annotation.translation_table)
    start_codons = frozenset(c.upper() for c in cfg.annotation.start_codons)
    views: dict[tuple[str, str], _Strand] = {}
    groups: dict[tuple[str, str, int, bool], list[tuple[int, int, dict[str, Any]]]] = {}
    uniq = rows.drop_duplicates(["peptide", "contig", "start", "end", "strand"])
    for rec in uniq.to_dict("records"):
        contig, strand = str(rec["contig"]), str(rec["strand"])
        key = (contig, strand)
        if key not in views:
            views[key] = _strand_view(contig, strand, gm.replicons[contig], genome)
        sv = views[key]
        ts, te = sv.to_t(int(rec["start"]), int(rec["end"]))
        end_t, open3 = _walk_down(sv, ts, te, stops)
        stop_key = end_t if open3 else (end_t - 3) % sv.n
        groups.setdefault((contig, strand, stop_key, open3), []).append((end_t - ts, ts, rec))
    orfs = [
        _build(views[(c, s)], members, open3, gm, stops, start_codons)
        for (c, s, _, open3), members in groups.items()
    ]
    return sorted(orfs, key=lambda o: (o.contig, _id_pos(o), o.strand, o.stop))


def _id_pos(o: Orf) -> int:
    return int(o.orf_id.rsplit("_", 2)[1])


def _build(
    sv: _Strand,
    members: list[tuple[int, int, dict[str, Any]]],
    open3: bool,
    gm: GeneModel,
    stops: frozenset[str],
    start_codons: frozenset[str],
) -> Orf:
    members = sorted(members, key=lambda m: (-m[0], str(m[2]["peptide"])))
    dist, first, top = members[0]
    end_t = first + dist
    q, open5 = _walk_up(sv, first, stops)
    t_starts = [x for x in range(q + 3, first + 1, 3) if sv.codon(x) in start_codons]
    n = sv.n
    plus_low = t_starts[0] if t_starts else first
    lowest = plus_low % n if sv.strand == "+" else (n - end_t) % n
    classes = {str(m[2]["locus_class"]) for m in members}
    cls = next(c for c in CLASS_PRIORITY if c in classes)
    cls_row = next(m[2] for m in members if m[2]["locus_class"] == cls)
    starts = tuple(sv.codon_genome(x) for x in t_starts)
    return Orf(
        orf_id=f"NOVEL_{sv.contig}_{lowest + 1}_{'p' if sv.strand == '+' else 'm'}",
        cls=cast(NovelClass, cls),
        contig=sv.contig,
        strand=cast(Literal["+", "-"], sv.strand),
        frame=int(top["frame"]),
        stop=sv.codon_genome(end_t - 3),
        starts=starts,
        start_longest=starts[0] if starts else None,
        end_5p_open=open5,
        end_3p_open=open3,
        peptides=tuple(sorted({_il(str(m[2]["peptide"])) for m in members})),
        best_q=_best_q([m[2].get("q_value", float("nan")) for m in members]),
        overlapping_gene=_overlapping_gene(cls, cls_row, gm),
    )


# ---------------------------------------------------------------- rules 7-8: second round


def promoter_windows(
    orfs: Sequence[Orf], gm: GeneModel, genome: Mapping[str, str], cfg: GiapConfig
) -> list[Window]:
    """M2-13 ``gene_windows`` at every candidate start, unique by ``(contig, strand, end3p)``."""
    pcfg = cfg.annotation.promoters
    lo, hi = pcfg.window_3p_range
    seen: dict[tuple[str, str, int], Window] = {}
    for o in orfs:
        rep = gm.replicons[o.contig]
        for x in o.starts:
            anchor = Gene(o.orf_id, "cds", Interval(o.contig, x, x + 3, o.strand), "", None, None)
            for w in gene_windows(anchor, rep, genome[o.contig], (lo, hi), pcfg.window_nt):
                seen.setdefault((w.contig, w.strand, w.end3p), w)
    return list(seen.values())


def window_chunks(
    windows: Sequence[Window], chunk_size: int = 100_000
) -> tuple[list[list[tuple[str, str]]], pd.DataFrame]:
    """``dedupe_and_chunk`` of the second-round windows; the map gains ``orf_id`` (= owner)."""
    chunks, wmap = dedupe_and_chunk(windows, chunk_size)
    return chunks, wmap.assign(orf_id=wmap["owner"])


def ostir_rows(
    orfs: Sequence[Orf],
    gm: GeneModel,
    genome: Mapping[str, str],
    anti_sd: str,
    cfg: GiapConfig,
) -> list[OstirRow]:
    """One row per ORF with a start: flanks around all its starts, ``start``/``end`` 1-based."""
    up, down = cfg.annotation.ostir.upstream_nt, cfg.annotation.ostir.downstream_nt
    out: list[OstirRow] = []
    for o in orfs:
        if not o.starts:
            continue
        sv = _strand_view(o.contig, o.strand, gm.replicons[o.contig], genome)
        n = sv.n
        t0, t1 = (x if o.strand == "+" else n - x - 3 for x in (o.starts[0], o.starts[-1]))
        span = (t1 - t0) % n if sv.circular else t1 - t0
        lo, hi = t0 - up, t0 + span + down
        if sv.circular:
            seq = "".join(sv.tseq[p % n] for p in range(lo, hi))
            kept_up = up
        else:
            seq = sv.tseq[max(lo, 0) : min(hi, n)]
            kept_up = t0 - max(lo, 0)
        out.append(OstirRow(o.orf_id, seq, kept_up + 1, kept_up + 1 + span, anti_sd))
    return out


# ---------------------------------------------------------------- rule 10: pass-1 call


def minimal_call(orfs: Sequence[Orf], cfg: GiapConfig) -> pd.DataFrame:
    """§7.6.4 "Neither" for every ORF: ``called`` = distinct peptides >= the no-decoy minimum."""
    need = cfg.novel_orfs.min_peptides_without_decoys
    rows = []
    for o in orfs:
        flags = {"class_fdr_unestimated"}
        if not o.starts:
            flags.add("no_start_codon")
        if o.end_5p_open or o.end_3p_open:
            flags.add("partial")
        n = len(o.peptides)
        rows.append((o.orf_id, n, n >= need, "unestimated", ";".join(sorted(flags))))
    df = pd.DataFrame(rows, columns=CALL_COLUMNS)
    return df.astype({"n_peptides": "int64", "called": "bool"})
