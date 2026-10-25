"""Novel-ORF class evidence, RBS support, promoter/terminator evidence and ranking (issue M7-09;
GIAP_SPEC.md §7.6.3-§7.6.5, §9.4; ADR-0012, ADR-0013).

Pure. ``apply_evidence`` fills the M7-07 ``novel_orfs`` table: ``rbs_supported`` (the chosen
start's OSTIR score against the ``rbs_min_quantile`` quantile of detected genes, rule 2),
``promoter_supported``, the ADR-0012 terminator columns, the class labels, ``called`` (the
error-control rule 8 and the class rule) and ``rank`` (called ORFs only, 1 = best). Promoter and
terminator evidence never gate a call. Table positions are 0-based first bases in transcription
direction (M7-07); intervals are 0-based half-open.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from giap.annotate.terminators import NOT_PRIMARY, Terminator
from giap.config.model import GiapConfig
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel
from giap.novel.candidates import Orf, _revcomp
from giap.novel.outputs import is_na, novel_frame, order_flags

RBS_CLASSES = frozenset({"nterm_extension", "intergenic"})
RBS_MIN_GENES = 20  # decision 1: the M6 floor for a fitted distribution
TWO_PEPTIDE_CLASSES = frozenset({"cds_outframe", "antisense"})
TERMINATOR_WINDOW_NT = 300  # ADR-0012 d1
HOMOPOLYMER_MIN = 5
_HOMOPOLYMER_WINDOW = (-5, 8)  # decision 2: [junction - 5, junction + 8)
_STOPS = frozenset({"TAA", "TAG", "TGA"})
_CODON_RANK = {"ATG": 0, "GTG": 1, "TTG": 2}
_OTHER_CODON, _NO_CODON = 3, 4


# ---------------------------------------------------------------- rules 2-3: RBS support


def rbs_threshold(gene_table: pd.DataFrame, q: float, min_n: int = RBS_MIN_GENES) -> float | None:
    """The ``q`` quantile (linear) of ``tir_ostir`` over non-partial detected genes with a score;
    ``None`` when fewer than ``min_n`` such genes."""
    partial = gene_table["partial"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    detected = gene_table["detected"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    tir = pd.to_numeric(gene_table["tir_ostir"], errors="coerce").to_numpy(dtype=float)
    vals = tir[~partial & detected & ~np.isnan(tir)]
    if len(vals) < min_n:
        return None
    return float(np.quantile(vals, q, method="linear"))


def rbs_supported(tir: object, rbs_min: float | None) -> bool | None:
    """Rule 3: ``tir ≥ rbs_min``; ``None`` when either is NA."""
    if rbs_min is None or is_na(tir):
        return None
    return float(tir) >= rbs_min  # type: ignore[arg-type]


def denovo_only_unsupported(cls: str, has_denovo: bool, db_support: bool, rbs: bool | None) -> bool:
    """M7-17 rule 4 (ADR-0050): de novo peptides, no database PSM and no RBS support beyond the
    class rule (``rbs`` NA is no support; ``RBS_CLASSES`` are never exempt)."""
    return has_denovo and not db_support and not (rbs is True and cls not in RBS_CLASSES)


# ---------------------------------------------------------------- rule 4: class evidence


def class_evidence(
    row: Mapping[str, object], cfg: GiapConfig, rbs_min: float | None
) -> tuple[bool, list[str]]:
    """``(class rule of rule 4, labels)`` for a ``novel_orfs`` row plus ``start_upstream_nt``
    (nt from ``start_longest`` to the annotated start, positive = upstream) and
    ``spans_annotated_stop``. Rule 8 (the peptide minimum) is not part of it."""
    cls = str(row["class"])
    labels: set[str] = set()
    rbs = bool(rbs_supported(row.get("ostir_tir"), rbs_min))
    if rbs:
        labels.add("rbs_supported")
    if cls in RBS_CLASSES and rbs_min is None:
        labels.add("rbs_threshold_unestimated")
    has_start = not is_na(row.get("start_longest"))
    ok: bool
    if cls == "nterm_extension":
        up = row.get("start_upstream_nt")
        ok = has_start and not is_na(up) and int(up) > 0 and rbs  # type: ignore[call-overload]
    elif cls == "intergenic":
        n_aa = row.get("length_aa_best")
        long_enough = not is_na(n_aa) and int(n_aa) >= cfg.novel_orfs.min_codons  # type: ignore[call-overload]
        ok = has_start and long_enough and rbs
    elif cls in TWO_PEPTIDE_CLASSES:
        ok = int(row["n_peptides"]) >= 2  # type: ignore[call-overload]
    elif cls in ("pseudogene", "sorf"):
        ok = True
    elif cls == "cterm_extension":  # M7-16 rule 12: the ORF is the in-frame continuation
        ok = row.get("spans_annotated_stop") is True
        if ok:
            labels.add("readthrough_candidate")
    else:
        raise ValueError(f"unknown novel ORF class {cls!r}")
    flags = order_flags(labels)
    return ok, [] if flags is None else flags.split(";")


# ---------------------------------------------------------------- rule 7: terminators


def terminator_downstream(
    orf: Orf, terminators: Sequence[Terminator], rep: Replicon
) -> tuple[float | None, int | None]:
    """``(conf, distance)`` of the best same-strand call whose 5' end lies 0..300 nt after the
    stop codon's last base (wrapping on circular replicons); equal confidences → the nearer.
    BacTermFinder calls are skipped (M2-15)."""
    n, circular = rep.length, rep.topology == "circular"
    plus = orf.strand == "+"
    last = orf.stop + 2 if plus else orf.stop
    best: tuple[float, int] | None = None
    for t in terminators:
        iv = t.interval
        if t.source == NOT_PRIMARY or iv.contig != orf.contig or iv.strand != orf.strand:
            continue
        d = iv.start - last if plus else last - (iv.end - 1)
        if circular:
            d %= n
        if not 0 <= d <= TERMINATOR_WINDOW_NT:
            continue
        if best is None or t.confidence > best[0] or (t.confidence == best[0] and d < best[1]):
            best = (t.confidence, d)
    return (None, None) if best is None else best


# ---------------------------------------------------------------- rule 5: pseudogene junction


def _codon_at(genome: str, lo: int, n: int) -> str:
    return "".join(genome[(lo + i) % n] for i in range(3))


def pseudogene_junction(gene: Gene, genome: str, rep: Replicon) -> int | None:
    """Lowest coordinate (wrapped) of the first in-frame stop codon from the annotated start
    that is not the gene's last codon; ``None`` if there is none."""
    iv, n, g = gene.interval, rep.length, genome.upper()
    if iv.strand == "+":
        for lo in range(iv.start, iv.end - 5, 3):
            if _codon_at(g, lo, n) in _STOPS:
                return lo % n
    else:
        for lo in range(iv.end - 3, iv.start + 2, -3):
            if _revcomp(_codon_at(g, lo, n)) in _STOPS:
                return lo % n
    return None


def _contains_codon(iv: Interval, lo: int, rep: Replicon) -> bool:
    """``[lo, lo + 3)`` lies inside ``iv`` (along the origin on a circular replicon)."""
    if rep.topology == "circular":
        return (lo - iv.start) % rep.length + 3 <= iv.end - iv.start
    return iv.start <= lo and lo + 3 <= iv.end


def _spans_stop(iv: Interval, lo: int, rep: Replicon) -> bool:
    """``iv`` contains the stop codon ``[lo, lo + 3)`` and extends past it on the 3' side."""
    d = lo - iv.start
    if rep.topology == "circular":
        d %= rep.length
    if d < 0 or d + 3 > iv.end - iv.start:
        return False
    return d + 3 < iv.end - iv.start if iv.strand == "+" else d > 0


def _run_length(g: str, p: int, rep: Replicon) -> int:
    """Length of the maximal run of identical bases through position ``p`` (forward strand)."""
    n, circular = rep.length, rep.topology == "circular"

    def base(i: int) -> str | None:
        if circular:
            return g[i % n]
        return g[i] if 0 <= i < n else None

    b = base(p)
    if b is None:
        return 0
    lo = hi = p
    while hi - lo + 1 < n and base(lo - 1) == b:
        lo -= 1
    while hi - lo + 1 < n and base(hi + 1) == b:
        hi += 1
    return hi - lo + 1


def junction_labels(
    junction: int, peptides: Sequence[Interval], genome: str, rep: Replicon
) -> set[str]:
    """``junction_spanned`` (a peptide contains the junction codon) and ``homopolymer_junction``
    (a run of ≥ 5 identical bases overlaps [junction - 5, junction + 8))."""
    out: set[str] = set()
    if any(_contains_codon(iv, junction, rep) for iv in peptides):
        out.add("junction_spanned")
    g = genome.upper()
    a, b = _HOMOPOLYMER_WINDOW
    if any(_run_length(g, p, rep) >= HOMOPOLYMER_MIN for p in range(junction + a, junction + b)):
        out.add("homopolymer_junction")
    return out


# ---------------------------------------------------------------- rule 10: ranking


def _desc(v: object) -> tuple[int, float]:
    return (1, 0.0) if is_na(v) else (0, -float(v))  # type: ignore[arg-type]


def _asc(v: object) -> tuple[int, float]:
    return (1, 0.0) if is_na(v) else (0, float(v))  # type: ignore[arg-type]


def rank_key(row: Mapping[str, object]) -> tuple[object, ...]:
    """G13 + ADR-0012 d4: peptides ↓, q ↑, start codon ATG < GTG < TTG < other < none,
    TIR ↓, promoter ↓, terminator ↓, length ↓, ``orf_id`` ↑; NA last in every key."""
    c = row.get("start_codon_best")
    codon = _NO_CODON if is_na(c) else _CODON_RANK.get(str(c).upper(), _OTHER_CODON)
    return (
        _desc(row["n_peptides"]),
        _asc(row.get("best_peptide_q")),
        codon,
        _desc(row.get("ostir_tir")),
        _desc(row.get("promoter_score")),
        _desc(row.get("terminator_conf")),
        _desc(row.get("length_aa_best")),
        str(row["orf_id"]),
    )


# ---------------------------------------------------------------- apply_evidence


def _signed_offset(d: int, rep: Replicon) -> int:
    """``d`` on a circular replicon folded to (-n/2, n/2]."""
    if rep.topology != "circular":
        return d
    d %= rep.length
    return d - rep.length if d > rep.length // 2 else d


def _gene_geometry(
    r: Mapping[str, Any], gene: Gene | None, peptides: Sequence[Interval], rep: Replicon
) -> dict[str, Any]:
    """``start_upstream_nt`` and ``spans_annotated_stop`` of a row against its gene."""
    out: dict[str, Any] = {"start_upstream_nt": None, "spans_annotated_stop": False}
    if gene is None:
        return out
    g = gene.interval
    plus = g.strand == "+"
    if not is_na(r["start_longest"]):
        sl = int(r["start_longest"])
        out["start_upstream_nt"] = _signed_offset(g.start - sl if plus else sl - (g.end - 1), rep)
    stop_lo = (g.end - 3) % rep.length if plus else g.start
    out["spans_annotated_stop"] = any(_spans_stop(iv, stop_lo, rep) for iv in peptides)
    return out


def _gene_of(tag: object, gm: GeneModel) -> Gene | None:
    if is_na(tag):
        return None
    try:
        return gm.by_locus(str(tag))
    except KeyError:
        return None


def _row_evidence(
    r: Mapping[str, Any], o: Orf, gm: GeneModel, cfg: GiapConfig, rbs_min: float | None,
    intervals: Sequence[Interval],
) -> tuple[bool, list[str], Gene | None]:  # fmt: skip
    gene = _gene_of(r["overlapping_gene"], gm)
    rep = gm.replicons[o.contig]
    ok, labels = class_evidence({**r, **_gene_geometry(r, gene, intervals, rep)}, cfg, rbs_min)
    return ok, labels, gene


def class_rules(
    table: pd.DataFrame,
    orfs: Sequence[Orf],
    *,
    gene_table: pd.DataFrame,
    gm: GeneModel,
    genome: Mapping[str, str],
    cfg: GiapConfig,
    peptide_intervals: Mapping[str, Sequence[Interval]] | None = None,
) -> dict[str, bool]:
    """``orf_id`` -> the M7-09 rule-4 class rule of every ``orf_table`` row (no peptide minimum
    and no error control)."""
    by_id = {o.orf_id: o for o in orfs}
    rbs_min = rbs_threshold(gene_table, cfg.novel_orfs.rbs_min_quantile)
    peps = peptide_intervals or {}
    return {
        str(r["orf_id"]): _row_evidence(r, by_id[str(r["orf_id"])], gm, cfg, rbs_min,
                                        peps.get(str(r["orf_id"]), ()))[0]
        for r in table.to_dict("records")
    }  # fmt: skip


def _call(cell: object, q: float | None, rule: bool, n_peptides: int, need: int,
          alpha: float) -> bool:  # fmt: skip
    """Rule 6 (M7-11b): a numeric class FDR q calls at q <= alpha (no peptide minimum); an
    ``unestimated`` cell needs ``need`` peptides (§7.6.4 "Neither")."""
    if is_na(cell) or cell == "unestimated":
        return rule and n_peptides >= need
    return rule and (float(cell) if q is None else q) <= alpha  # type: ignore[arg-type]


def apply_evidence(
    table: pd.DataFrame,
    orfs: Sequence[Orf],
    *,
    gene_table: pd.DataFrame,
    terminators: Sequence[Terminator] | None,
    gm: GeneModel,
    genome: Mapping[str, str],
    cfg: GiapConfig,
    peptide_intervals: Mapping[str, Sequence[Interval]] | None = None,
    class_q: Mapping[str, float] | None = None,
    denovo_support: Mapping[str, tuple[bool, bool]] | None = None,
) -> pd.DataFrame:
    """The ``orf_table`` frame (same rows and order) with the M7-09 columns filled and its
    ``class_fdr`` cells kept. ``terminators`` is ``None`` when the layer is missing (its three
    columns NA); ``peptide_intervals`` maps an ``orf_id`` to its peptides' genome intervals;
    ``class_q`` gives the unformatted q of the numeric cells (the step passes it, M7-11b lead
    decision 6; without it the cell is read back); ``denovo_support`` maps an ``orf_id`` to
    ``db_supported``'s pair and sets ``denovo_only_unsupported`` (M7-17; no other effect)."""
    by_id = {o.orf_id: o for o in orfs}
    rbs_min = rbs_threshold(gene_table, cfg.novel_orfs.rbs_min_quantile)
    p_yes = cfg.annotation.promoters.promotech_yes
    t_yes = cfg.annotation.terminators.terminatornet_yes
    need = cfg.novel_orfs.min_peptides_without_decoys
    alpha = cfg.stats.alpha.novel_orfs
    peps = peptide_intervals or {}
    recs: list[dict[str, Any]] = [dict(r) for r in table.to_dict("records")]
    for r in recs:
        o = by_id[str(r["orf_id"])]
        rep = gm.replicons[o.contig]
        flags = set() if is_na(r["flags"]) else set(str(r["flags"]).split(";")) - {""}
        intervals = peps.get(o.orf_id, ())
        ok, labels, gene = _row_evidence(r, o, gm, cfg, rbs_min, intervals)
        flags |= set(labels)
        if o.cls == "pseudogene":
            j = None if gene is None else pseudogene_junction(gene, genome[o.contig], rep)
            if j is None:
                flags.add("junction_unknown")
            else:
                flags |= junction_labels(j, intervals, genome[o.contig], rep)
        r["rbs_supported"] = rbs_supported(r["ostir_tir"], rbs_min)
        dn, db = (denovo_support or {}).get(o.orf_id, (False, False))
        if denovo_only_unsupported(o.cls, dn, db, r["rbs_supported"]):
            flags.add("denovo_only_unsupported")
        ps = r["promoter_score"]
        r["promoter_supported"] = None if is_na(ps) else bool(float(ps) >= p_yes)
        if r["promoter_supported"]:
            flags.add("promoter_supported")
        if terminators is None:
            r["terminator_conf"] = r["terminator_distance_nt"] = r["terminator_downstream"] = None
        else:
            conf, dist = terminator_downstream(o, terminators, rep)
            r["terminator_conf"], r["terminator_distance_nt"] = conf, dist
            r["terminator_downstream"] = conf is not None and conf >= t_yes
            if r["terminator_downstream"]:
                flags.add("terminator_supported")
        q = None if class_q is None else class_q.get(o.orf_id)
        r["called"] = _call(r["class_fdr"], q, ok, int(r["n_peptides"]), need, alpha)
        r["flags"] = order_flags(flags)
        r["rank"] = None
    for k, r in enumerate(sorted((r for r in recs if r["called"]), key=rank_key), start=1):
        r["rank"] = k
    return novel_frame(recs)
