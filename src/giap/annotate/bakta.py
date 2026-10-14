"""Defensive parser for Bakta 1.12.2 JSON (GIAP_SPEC.md §7.1.1, §7.1.2, §6.5; M2-03).

``parse_bakta_json`` returns the shared M3-01 ``GeneModel`` plus the contig sequences, the Bakta
and DB versions and a count of skipped (unknown) feature types. Field-name variations are
tolerated: a feature's contig is ``sequence`` or ``contig``; sequence records are ``sequences``
or ``contigs`` with nucleotides in ``nt`` or ``sequence``. Bakta's 1-based closed ``start``/``stop``
become 0-based half-open intervals; ``start > stop`` (``edge``) is unwrapped to ``end > length`` on
circular replicons and is an error on linear ones. Every malformed input raises
``BaktaParseError`` (a ``ValueError``) naming the feature.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Literal, cast

from giap.io._open import read_json
from giap.model.coords import (
    Interval,
    Replicon,
    RepliconKind,
    Strand,
    from_one_based_closed,
    interval_sort_key,
)
from giap.model.genes import _NC_TYPES, Gene, GeneModel, GeneType, NcFeature

if TYPE_CHECKING:
    from _typeshed import StrPath

DbType = Literal["full", "light", "unknown"]

_HYPOTHETICAL = re.compile(r"^hypothetical protein")
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_GENE_TYPES = ("cds", "sorf")


class BaktaParseError(ValueError):
    """Malformed or inconsistent Bakta JSON."""


@dataclass(frozen=True, slots=True)
class BaktaAnnotation:
    model: GeneModel
    sequences: Mapping[str, str]
    bakta_version: str
    db_version: str
    db_type: DbType
    skipped_types: Mapping[str, int]


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _slice(seq: str, iv: Interval) -> str:
    """Strand-aware, origin-aware nucleotides of ``iv`` (``end`` may exceed ``len(seq)``)."""
    n = len(seq)
    s = seq[iv.start : iv.end] if iv.end <= n else seq[iv.start :] + seq[: iv.end - n]
    return _revcomp(s) if iv.strand == "-" else s


def _record_nt(rec: Mapping[str, Any], cid: str) -> str | None:
    nt = rec.get("nt")
    if nt is None:
        nt = rec.get("sequence")
    if nt is not None and not isinstance(nt, str):
        raise BaktaParseError(f"sequence {cid!r}: nt must be a string, got {type(nt).__name__}")
    return nt or None


def _replicons_from_json(
    records: list[Any],
    replicons: Mapping[str, Replicon] | None,
    genome: Mapping[str, str] | None,
) -> tuple[dict[str, Replicon], dict[str, str]]:
    reps: dict[str, Replicon] = {}
    seqs: dict[str, str] = {}
    for rec in records:
        if not isinstance(rec, dict) or not rec.get("id"):
            raise BaktaParseError(f"sequence record without id: {rec!r:.200}")
        cid = str(rec["id"])
        nt = _record_nt(rec, cid)
        if replicons is not None and cid not in replicons:
            raise BaktaParseError(f"sequence {cid!r} is not among the given replicons")
        if nt is None and genome is not None and cid in genome:
            nt = genome[cid]
        if nt is not None:
            seqs[cid] = nt
        elif replicons is None:
            raise BaktaParseError(
                f"sequence {cid!r} has no nucleotides in the JSON and no genome was given"
            )
        kind = str(rec.get("type", "")).lower()
        topology = str(rec.get("topology") or "linear").lower()
        if topology not in ("circular", "linear"):
            raise BaktaParseError(f"sequence {cid!r}: unknown topology {topology!r}")
        length = replicons[cid].length if replicons is not None else len(seqs[cid])
        declared = rec.get("length")
        if declared is not None and declared != length:
            raise BaktaParseError(
                f"sequence {cid!r}: JSON length {declared!r} != genome sequence length {length}"
            )
        if cid in seqs and len(seqs[cid]) != length:
            raise BaktaParseError(
                f"sequence {cid!r}: nucleotide length {len(seqs[cid])} != replicon length {length}"
            )
        reps[cid] = Replicon(
            cid,
            length,
            "circular" if topology == "circular" else "linear",
            cast(RepliconKind, kind if kind in ("chromosome", "plasmid") else "unknown"),
        )
    if replicons is None:
        return reps, seqs
    for cid in replicons:
        if cid not in seqs and genome is not None and cid in genome:
            seqs[cid] = genome[cid]
    return dict(replicons), seqs


def _int(f: Mapping[str, Any], key: str, name: str) -> int:
    v = f.get(key)
    if isinstance(v, bool) or not isinstance(v, int):
        raise BaktaParseError(f"feature {name!r}: {key} must be an int, got {v!r}")
    return v


def _opt_str(f: Mapping[str, Any], key: str, name: str) -> str | None:
    v = f.get(key)
    if v is not None and not isinstance(v, str):
        raise BaktaParseError(f"feature {name!r}: {key} must be a string, got {type(v).__name__}")
    return v


def _expert(f: Mapping[str, Any], name: str) -> tuple[str, ...]:
    raw = f.get("expert")
    if raw is None:
        return ()
    if not isinstance(raw, list) or not all(isinstance(e, dict) for e in raw):
        raise BaktaParseError(
            f"feature {name!r}: expert must be a list of objects, got {raw!r:.200}"
        )
    out: list[str] = []
    for e in raw:
        src = e.get("type") or e.get("source")
        if src and str(src) not in out:
            out.append(str(src))
    return tuple(out)


def _interval(
    f: Mapping[str, Any], name: str, replicons: Mapping[str, Replicon], *, coding: bool
) -> Interval:
    contig = f.get("sequence", f.get("contig"))
    if contig is None:
        raise BaktaParseError(f"feature {name!r}: no 'sequence' or 'contig' field")
    if contig not in replicons:
        raise BaktaParseError(f"feature {name!r}: unknown contig {contig!r}")
    raw = f.get("strand")
    if raw not in ("+", "-"):
        if coding or raw not in (None, "?", "."):
            raise BaktaParseError(f"feature {name!r}: invalid strand {raw!r}")
        raw = "."
    start, stop = _int(f, "start", name), _int(f, "stop", name)
    try:
        return from_one_based_closed(contig, start, stop, cast(Strand, raw), replicons[contig])
    except ValueError as exc:
        raise BaktaParseError(f"feature {name!r}: {exc}") from exc


def _gene(f: Mapping[str, Any], replicons: Mapping[str, Replicon], seqs: Mapping[str, str]) -> Gene:
    if not f.get("locus"):
        raise BaktaParseError(f"{f.get('type')} feature without locus: {f!r:.200}")
    locus = str(f["locus"])
    iv = _interval(f, locus, replicons, coding=True)
    ftype = cast(GeneType, str(f["type"]).lower())
    if ftype == "cds" and (f.get("pseudo") or "pseudogene" in f):
        ftype = "pseudogene"
    xrefs = f.get("db_xrefs") or []
    if not isinstance(xrefs, list | tuple):
        raise BaktaParseError(f"feature {locus!r}: db_xrefs must be a list, got {xrefs!r}")
    nt = _opt_str(f, "nt", locus) or (_slice(seqs[iv.contig], iv) if iv.contig in seqs else None)
    product = f.get("product") or ""
    hypothetical = (
        bool(f["hypothetical"]) if "hypothetical" in f else bool(_HYPOTHETICAL.match(product))
    )
    rep = replicons[iv.contig]
    at_end = rep.topology == "linear" and (iv.start == 0 or iv.end == rep.length)
    aa = _opt_str(f, "aa", locus) or ""
    start_type = f.get("start_type")
    # Pyrodigal partial start ("Edge") is not a start codon (M3-01 Q3).
    start_codon = None if start_type == "Edge" else start_type or (nt[:3] if nt else None)
    return Gene(
        locus_tag=locus,
        feature_type=ftype,
        interval=iv,
        aa=aa[:-1] if aa.endswith("*") else aa,
        nt=nt,
        start_codon=start_codon,
        gene=f.get("gene"),
        product=product,
        partial=bool(f.get("truncated")) or at_end,
        hypothetical=hypothetical,
        db_xrefs=tuple(str(x) for x in xrefs),
        expert=_expert(f, locus),
    )


def _nc(f: Mapping[str, Any], replicons: Mapping[str, Replicon]) -> NcFeature:
    fid = str(f.get("id") or f.get("locus") or "")
    anticodon = f.get("anti_codon")
    return NcFeature(
        id=fid,
        type=_NC_TYPES[str(f["type"]).lower()],
        interval=_interval(f, fid, replicons, coding=False),
        gene=f.get("gene"),
        product=f.get("product") or "",
        amino_acid=f.get("amino_acid"),
        anticodon=str(anticodon).upper().replace("U", "T") if anticodon else None,
        nt=_opt_str(f, "nt", fid),
    )


def _versions(data: Mapping[str, Any]) -> tuple[str, str, DbType]:
    version = data.get("version")
    if not isinstance(version, dict):
        return "", "", "unknown"
    raw_db = version.get("db")
    db: dict[str, Any] = raw_db if isinstance(raw_db, dict) else {}
    db_type = str(db.get("type", "")).lower()
    return (
        str(version.get("bakta") or ""),
        str(db.get("version") or ""),
        cast(DbType, db_type if db_type in ("full", "light") else "unknown"),
    )


def parse_bakta_json(
    path: StrPath,
    replicons: Mapping[str, Replicon] | None = None,
    genome: Mapping[str, str] | None = None,
) -> BaktaAnnotation:
    """Parse a (possibly gzipped) Bakta JSON; ``replicons`` (preflight's table) wins if given."""
    data = read_json(path)
    if not isinstance(data, dict):
        raise BaktaParseError(f"{path}: top-level JSON must be an object")
    records = data.get("sequences", data.get("contigs")) or []
    reps, seqs = _replicons_from_json(records, replicons, genome)
    genes: list[Gene] = []
    features: list[NcFeature] = []
    skipped: Counter[str] = Counter()
    seen: set[str] = set()
    for f in data.get("features") or []:
        if not isinstance(f, dict):
            raise BaktaParseError(f"feature must be an object: {f!r:.200}")
        t = str(f.get("type", "")).lower()
        if t in _GENE_TYPES:
            g = _gene(f, reps, seqs)
            if g.locus_tag in seen:
                raise BaktaParseError(f"duplicate locus tag {g.locus_tag}")
            seen.add(g.locus_tag)
            genes.append(g)
        elif t in _NC_TYPES:
            features.append(_nc(f, reps))
        else:
            skipped[str(f.get("type", ""))] += 1
    order = {c: k for k, c in enumerate(reps)}
    genes.sort(key=lambda g: interval_sort_key(g.interval, order, g.locus_tag))
    features.sort(key=lambda x: interval_sort_key(x.interval, order, x.id))
    bakta_version, db_version, db_type = _versions(data)
    return BaktaAnnotation(
        model=GeneModel(replicons=reps, genes=tuple(genes), features=tuple(features)),
        sequences=MappingProxyType(seqs),
        bakta_version=bakta_version,
        db_version=db_version,
        db_type=db_type,
        skipped_types=MappingProxyType(dict(skipped)),
    )
