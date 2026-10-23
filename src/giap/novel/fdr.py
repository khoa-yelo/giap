"""Class FDR of novel ORFs (issue M7-11; GIAP_SPEC.md §7.6.4, §8.1; ADR-0009, ADR-0010).

Three regimes, chosen per ORF from its peptides' source (lead decision 1, database first):

* ``picked_decoy``: decoy rows whose stripped IDs are class-tagged ``makedb`` IDs. A target ORF is
  paired with the makedb six-frame entry of the same ``(contig, strand, stop)`` (``nterm_extension``
  entries by locus tag); scores are (best input q, then best score) and a picked target-decoy
  competition runs within each entry class, decoy wins rescaled by the class's paired share
  (ADR-0009) and the class unestimated below an 80 % share (decision 4).
* ``denovo_subgroup``: each decoy hit takes the locus class of its mirrored target interval; per
  class the ADR-0010 weighted list ``(W(e) + c) / T(e)`` over peptides, ORF q = the minimum of its
  peptides' q (decision 3).
* ``unestimated``: everything else (the pass-1 >= 2-peptide rule applies downstream).

Frame layouts: aggregated C9 peptides, ``db_id_map.tsv`` as a frame (1-based closed), mapping hit
frames with ``e_target`` / ``e_decoy`` (``tests/helpers/novel_fdr_fixture.py``). Pure: no I/O.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, Protocol

import numpy as np
import pandas as pd

from giap.config.model import InputConfig
from giap.genetable.detection import MIN_DECOY_PAIRING
from giap.makedb.build import ORF_CLASSES
from giap.mapping.placement_fdr import stratum_qvalues
from giap.model.coords import Interval
from giap.novel.candidates import Orf
from giap.stats.tda import picked_fdr_rescaled

__all__ = [
    "ClassStats",
    "PickedClassFdr",
    "Regime",
    "denovo_class_q",
    "orf_regimes",
    "parse_orf_id",
    "picked_class_fdr",
    "picked_class_q",
    "regime",
]

Regime = Literal["picked_decoy", "denovo_subgroup", "unestimated"]
TAGGED_CLASSES: tuple[str, ...] = ("nterm_extension", *sorted(ORF_CLASSES))
DEFAULT_DECOY_PREFIXES: tuple[str, ...] = tuple(InputConfig().decoy_prefixes)
_STOP_LEN = 3


class Classifier(Protocol):
    """The M3-14 ``LocusClassifier`` interface used here."""

    def classify(self, iv: Interval) -> tuple[str, str | None]: ...


@dataclass(frozen=True)
class ClassStats:
    """Per entry class: picked target and decoy wins, paired share of its decoys."""

    n_target_wins: int
    n_decoy_wins: int
    paired_share: float
    estimated: bool
    reason: str | None = None  # why not estimated: paired_share_low, decoys_unscored


@dataclass(frozen=True)
class PickedClassFdr:
    """``q`` per ORF (NaN unless estimated and won), ``reason`` per unestimated ORF
    (``decoys_not_makedb``, ``untagged_class``, ``no_makedb_pair``, ``paired_share_low``,
    ``decoys_unscored``, ``lost_pick``, ``no_target_score``), ORFs whose class differs from their
    entry's tag."""

    q: dict[str, float]
    reason: dict[str, str]
    class_mismatch: frozenset[str]
    classes: dict[str, ClassStats]


def _il(s: str) -> str:
    return s.upper().replace("I", "L")


def _ids(x: object) -> list[str]:
    """Protein IDs of one C9 cell: a ``;``-joined string or a list-like; NA = none."""
    if isinstance(x, str):
        return [s for s in x.split(";") if s]
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return []
    return [str(s) for s in x] if isinstance(x, Iterable) else []


def _strip(pid: str, prefixes: Sequence[str]) -> str:
    for p in sorted(prefixes, key=len, reverse=True):  # longest prefix first, as M4-12
        if pid.startswith(p):
            return pid[len(p) :]
    return pid


def parse_orf_id(db_id: str) -> tuple[str, str, int, int, str] | None:
    """``(class, contig, start0, end0, strand)`` of ``GIAPORF_{class}_{contig}_{s1}_{e1}_{p|m}``
    for a makedb-tagged class; None otherwise."""
    if not db_id.startswith("GIAPORF_"):
        return None
    rest = db_id[len("GIAPORF_") :]
    for cls in TAGGED_CLASSES:
        if rest.startswith(cls + "_"):
            parts = rest[len(cls) + 1 :].rsplit("_", 3)
            if len(parts) == 4 and parts[1].isdigit() and parts[2].isdigit() and parts[3] in "pm":
                return cls, parts[0], int(parts[1]) - 1, int(parts[2]), "+-"["pm".index(parts[3])]
    return None


def _decoy_tags(peptides: pd.DataFrame, prefixes: Sequence[str]) -> list[tuple[str, int]]:
    """``(stripped class-tagged ID, row)`` for every decoy row ID."""
    out = []
    dec = peptides["is_decoy"].astype(bool).to_numpy()
    for i in np.flatnonzero(dec):
        for pid in _ids(peptides["protein_ids"].iat[i]):
            sid = _strip(pid, prefixes)
            if parse_orf_id(sid) is not None:
                out.append((sid, int(i)))
    return out


def orf_regimes(
    peptides: pd.DataFrame,
    orfs: Sequence[Orf],
    decoy_prefixes: Sequence[str] = DEFAULT_DECOY_PREFIXES,
) -> dict[str, Regime]:
    """Regime per ORF (rule 1, decision 1): ``picked_decoy`` when the run has makedb-tagged
    decoys and the ORF has a database peptide, else ``denovo_subgroup`` when it has a de novo
    peptide, else ``unestimated``."""
    makedb = bool(_decoy_tags(peptides, decoy_prefixes))
    tgt = peptides[~peptides["is_decoy"].astype(bool)]
    src: dict[str, set[str]] = {}
    for p, s in zip(tgt["peptide"].astype(str), tgt["source"].astype(str), strict=True):
        src.setdefault(_il(p), set()).add(s)
    out: dict[str, Regime] = {}
    for o in orfs:
        s = set().union(*(src.get(p, set()) for p in o.peptides))
        if makedb and s - {"denovo"}:
            out[o.orf_id] = "picked_decoy"
        elif "denovo" in s:
            out[o.orf_id] = "denovo_subgroup"
        else:
            out[o.orf_id] = "unestimated"
    return out


def regime(
    peptides: pd.DataFrame,
    orfs: Sequence[Orf],
    decoy_prefixes: Sequence[str] = DEFAULT_DECOY_PREFIXES,
) -> Regime:
    """Regime of ``orfs``: the first of picked, de novo, unestimated held by any of them."""
    held = set(orf_regimes(peptides, orfs, decoy_prefixes).values())
    r: Regime
    for r in ("picked_decoy", "denovo_subgroup"):
        if r in held:
            return r
    return "unestimated"


# ---------------------------------------------------------------- picked (rule 2)


@dataclass(frozen=True)
class _Entry:
    cls: str
    contig: str
    start0: int
    end0: int
    strand: str
    locus_tag: str | None


def _entries(id_map: pd.DataFrame | None, tagged: Iterable[str]) -> dict[str, _Entry]:
    if id_map is None:
        parsed = {sid: parse_orf_id(sid) for sid in tagged}
        return {k: _Entry(*v, None) for k, v in parsed.items() if v is not None}
    m = id_map[(id_map["kind"].astype(str) == "orf") & id_map["class"].isin(TAGGED_CLASSES)]
    return {
        str(r.db_id): _Entry(
            str(r.cls), str(r.contig), int(r.start) - 1, int(r.end), str(r.strand),
            None if pd.isna(r.locus_tag) else str(r.locus_tag),
        )
        for r in m.rename(columns={"class": "cls"}).itertuples(index=False)
    }  # fmt: skip


def _stop(e: _Entry, lengths: Mapping[str, int] | None) -> int:
    """Lowest coordinate of the closing stop codon; wrapped modulo the replicon length if given
    (decision 8: without lengths an origin-crossing entry matches no ORF)."""
    s = e.end0 if e.strand == "+" else e.start0 - _STOP_LEN
    n = None if lengths is None else lengths.get(e.contig)
    return s % n if n else s


def _orf_key(o: Orf, lengths: Mapping[str, int] | None) -> tuple[str, str, int]:
    """Pairing key of ``o``; a 3'-open ORF's ``stop`` is its last whole codon, so the key moves
    one codon to the makedb segment's end (decisions 13-14), wrapped like the entry key."""
    s = o.stop + (_STOP_LEN if o.strand == "+" else -_STOP_LEN) * o.end_3p_open
    n = None if lengths is None else lengths.get(o.contig)
    return o.contig, o.strand, s % n if n else s


def _class_rank(
    keys: Sequence[str],
    rows_t: Mapping[str, set[int]],
    rows_d: Mapping[str, set[int]],
    qv: np.ndarray,
    sc: np.ndarray,
) -> dict[int, float] | None:
    """Dense rank (higher = better) of a class's target and decoy rows on a key every row has
    (decisions 11-12, as M4-24): q when every q is finite, broken by score only when every decoy
    row has one (an unscored target ranks last on it); else score alone when every score is
    finite; else None (``decoys_unscored``). Equal keys tie, and ties go to the decoy."""
    dec = {i for k in keys for i in rows_d.get(k, ())}
    rows = sorted(dec | {i for k in keys for i in rows_t.get(k, ())})
    if not rows:
        return {}
    q, s = qv[rows], sc[rows]
    if np.isfinite(q).all():
        if np.isfinite(sc[sorted(dec)]).all():
            key = np.column_stack([q, np.where(np.isnan(s), np.inf, -s)])
        else:
            key = q[:, None]
    elif np.isfinite(s).all():
        key = -s[:, None]
    else:
        return None
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    return dict(zip(rows, (len(uniq) - np.ravel(inv)).astype(np.float64).tolist(), strict=True))


def picked_class_fdr(
    orfs: Sequence[Orf],
    peptides: pd.DataFrame,
    id_map: pd.DataFrame | None,
    *,
    decoy_prefixes: Sequence[str] = DEFAULT_DECOY_PREFIXES,
    replicon_lengths: Mapping[str, int] | None = None,
) -> PickedClassFdr:
    """Picked class FDR per ORF with makedb decoys (rule 2, decisions 2, 4, 5, 7-9).

    A decoy is paired when its stripped ID is an ``id_map`` entry (``id_map=None``: every tagged
    ID is an entry); its class is the ID's tag. Ties go to the decoy.
    """
    q = {o.orf_id: math.nan for o in orfs}
    reason: dict[str, str] = {}
    tags = _decoy_tags(peptides, decoy_prefixes)
    if not tags:
        return PickedClassFdr(q, dict.fromkeys(q, "decoys_not_makedb"), frozenset(), {})
    qv = peptides["q_value"].astype("float64").to_numpy()
    sc = peptides["score"].astype("float64").to_numpy()
    all_tagged = {_strip(pid, decoy_prefixes) for ids in peptides["protein_ids"]
                  for pid in _ids(ids)}  # fmt: skip
    entries = _entries(id_map, all_tagged)
    by_key = {(e.contig, e.strand, _stop(e, replicon_lengths)): k
              for k, e in entries.items() if e.cls != "nterm_extension"}  # fmt: skip
    by_tag = {e.locus_tag: k for k, e in entries.items()
              if e.cls == "nterm_extension" and e.locus_tag}  # fmt: skip
    is_t = ~peptides["is_decoy"].astype(bool).to_numpy()
    is_t &= peptides["source"].astype(str).to_numpy() != "denovo"
    t_rows: dict[str, list[int]] = {}
    for j in np.flatnonzero(is_t).tolist():
        t_rows.setdefault(_il(str(peptides["peptide"].iat[j])), []).append(j)
    rows_t: dict[str, set[int]] = {}
    of_orf: dict[str, str] = {}
    mismatch = set()
    for o in orfs:
        if o.cls not in TAGGED_CLASSES:
            reason[o.orf_id] = "untagged_class"
            continue
        k = by_key.get(_orf_key(o, replicon_lengths))
        if k is None and o.cls == "nterm_extension" and o.overlapping_gene:
            k = by_tag.get(o.overlapping_gene)
        if k is None:
            reason[o.orf_id] = "no_makedb_pair"
            continue
        of_orf[o.orf_id] = k
        if entries[k].cls != o.cls:
            mismatch.add(o.orf_id)  # decision 9: the entry tag wins
        mine = [i for p in o.peptides for i in t_rows.get(p, [])]
        if mine:
            rows_t.setdefault(k, set()).update(mine)
    rows_d: dict[str, set[int]] = {}
    n_dec: dict[str, set[str]] = {}
    for sid, i in tags:
        n_dec.setdefault(parse_orf_id(sid)[0], set()).add(sid)  # type: ignore[index]
        if sid in entries:
            rows_d.setdefault(sid, set()).add(i)
    entry_q: dict[str, float] = {}
    classes: dict[str, ClassStats] = {}
    for cls in sorted({e.cls for e in entries.values()}):
        seen = rows_t.keys() | rows_d.keys()
        keys = sorted(k for k, e in entries.items() if e.cls == cls and k in seen)
        dec = n_dec.get(cls, set())
        f = len(dec & entries.keys()) / len(dec) if dec else 1.0
        rank = _class_rank(keys, rows_t, rows_d, qv, sc)
        why = "paired_share_low" if f < MIN_DECOY_PAIRING or not keys else None
        why = why or (None if rank is not None else "decoys_unscored")
        if why is not None or rank is None:
            classes[cls] = ClassStats(0, 0, f, False, why)
            continue
        t = [max((rank[i] for i in rows_t.get(k, ())), default=math.nan) for k in keys]
        d = [max((rank[i] for i in rows_d.get(k, ())), default=math.nan) for k in keys]
        res = picked_fdr_rescaled(t, d, f)
        classes[cls] = ClassStats(res.n_target_wins, res.n_decoy_wins, f, True)
        entry_q.update(zip(keys, res.q.tolist(), strict=True))
    for oid, k in of_orf.items():
        st = classes[entries[k].cls]
        if not st.estimated:
            reason[oid] = st.reason or "paired_share_low"
        elif k not in rows_t:
            reason[oid] = "no_target_score"
        elif math.isnan(entry_q[k]):
            reason[oid] = "lost_pick"
        else:
            q[oid] = entry_q[k]
    return PickedClassFdr(q, reason, frozenset(mismatch), classes)


def picked_class_q(
    orfs: Sequence[Orf], peptides: pd.DataFrame, id_map: pd.DataFrame | None
) -> dict[str, float]:
    """``orf_id`` -> picked class q (NaN when unestimated or the decoy wins)."""
    return picked_class_fdr(orfs, peptides, id_map).q


# ---------------------------------------------------------------- de novo (rule 3)


def denovo_class_q(
    orfs: Sequence[Orf],
    target_hits: pd.DataFrame,
    decoy_hits: pd.DataFrame,
    classifier: Classifier,
    pi: Mapping[str, float],
) -> dict[str, float]:
    """``orf_id`` -> de novo subgroup class q (rule 3, decision 3; NaN without a target hit).

    Per class (the ORF's for targets, the mirrored interval's for decoys at each peptide's best
    decoy cost): ``stratum_qvalues`` over peptides by E with decoy weights π̂ (missing = 1) and
    c = the largest π̂ over every searched peptide (ADR-0010).

    Raises:
        ValueError: if a target or decoy hit has no E (``e_target`` / ``e_decoy`` NaN).
    """
    for frame, col in ((target_hits, "e_target"), (decoy_hits, "e_decoy")):
        bad = frame["peptide"][frame[col].astype("float64").isna().to_numpy()]
        if len(bad):
            raise ValueError(f"{col} is missing for peptide {bad.iat[0]!r} (decision 11)")
    pif = {_il(k): float(v) for k, v in pi.items()}
    t_peps = [_il(p) for p in target_hits["peptide"].astype(str)]
    d_peps = [_il(p) for p in decoy_hits["peptide"].astype(str)]
    c = max((pif.get(p, 1.0) for p in (*t_peps, *d_peps)), default=1.0)
    e_t: dict[str, float] = {}
    for p, e in zip(t_peps, target_hits["e_target"].astype("float64"), strict=True):
        e_t[p] = min(e_t.get(p, math.inf), e)
    decoys: dict[tuple[str, str], float] = {}
    if len(decoy_hits):
        nd = decoy_hits["n_diff"].to_numpy()
        best = pd.Series(nd).groupby(d_peps).transform("min").to_numpy()
        cls_of: dict[Interval, str] = {}
        for i in np.flatnonzero(nd == best):
            r = decoy_hits.iloc[i]
            iv = Interval(str(r["contig"]), int(r["start"]), int(r["end"]), r["strand"])
            if iv not in cls_of:
                cls_of[iv] = classifier.classify(iv)[0]
            key = (cls_of[iv], d_peps[i])
            decoys[key] = min(decoys.get(key, math.inf), float(r["e_decoy"]))
    lists: dict[str, dict[str, float]] = {}
    for o in orfs:
        for p in o.peptides:
            if p in e_t:
                lists.setdefault(o.cls, {})[p] = e_t[p]
    pep_q: dict[tuple[str, str], float] = {}
    for cls, tl in lists.items():
        dk = [(p, e) for (k, p), e in decoys.items() if k == cls]
        qs = stratum_qvalues(list(tl.values()), [e for _, e in dk],
                             [pif.get(p, 1.0) for p, _ in dk], c)  # fmt: skip
        pep_q.update({(cls, p): float(v) for p, v in zip(tl, qs, strict=True)})
    return {
        o.orf_id: min((pep_q[(o.cls, p)] for p in o.peptides if (o.cls, p) in pep_q),
                      default=math.nan)
        for o in orfs
    }  # fmt: skip
