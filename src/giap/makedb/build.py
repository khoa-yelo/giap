"""`giap makedb`: optional strain search database (issue M3-23; GIAP_SPEC.md §7.9.2, §5.1).

Writes ``strain_db.faa`` (Bakta CDS/sORF proteins, patched proteins, class-tagged six-frame ORFs,
contaminants; wrapped at 60, no ``*``, no decoys), ``db_id_map.tsv`` (one row per record, 1-based
closed coordinates, unwrapped across the origin) and ``makedb_manifest.json``.

Six-frame ORFs are the maximal stop-to-stop segments of each frame, split at ``X`` (lead decision
3), translated codon by codon from the genome (no I -> L, no start rule) with the M3-01
``translate11_text``. An ORF in frame with an annotated CDS/sORF is dropped; each such gene with
>= 7 stop-free upstream codons gets one ``nterm_extension`` entry = upstream part + the first 30
codons of the gene translated from the genome (GTG -> V; review decision 5), whose header/id_map
interval covers both (lead decision 2). The entry is undefined (not written) when the upstream
walk runs off a linear contig, meets an ``X`` or finds no stop within one circular turn (review
decision 4). The other ORFs take the M3-14 ``LocusClassifier`` class; ``rna_overlap``/
``pseudogene`` ORFs are omitted. A patch CDS with an internal ``*`` and an annotation without
sequences or features (``{}``) are rejected (ValueError; review decision 7).
"""

from __future__ import annotations

import datetime
import hashlib
import json
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import giap
from giap.annotate.bakta import interval_nt, parse_bakta_json
from giap.io._open import read_json
from giap.io.fasta import read_genome_fasta
from giap.io.gff3 import AttrValue, read_gff3
from giap.io.protein_fasta import read_protein_fasta, write_protein_fasta
from giap.io.tsv import write_tsv
from giap.makedb import ORF_ID_PREFIX
from giap.mapping.loci import LocusClassifier, in_frame
from giap.mapping.sixframe import SixFrameIndex, translate11_text
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel

__all__ = ["MakedbOptions", "MakedbResult", "OrfRecord", "build_strain_db", "six_frame_orfs"]

ORF_CLASSES = frozenset({"intergenic", "cds_outframe", "antisense"})
NTERM_PROTEIN_AA = 30
NTERM_MIN_UPSTREAM_AA = 7
ID_MAP_COLUMNS = ("db_id", "kind", "locus_tag", "class", "contig", "start", "end", "strand")
_CODING = frozenset({"cds", "sorf"})
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _shifts(rep: Replicon) -> tuple[int, ...]:
    return (0, rep.length, -rep.length) if rep.topology == "circular" else (0,)


@dataclass(frozen=True)
class MakedbOptions:
    annotation: Path
    genome: Path | None
    outdir: Path
    patch: Path | None = None
    six_frame_orfs: bool = False
    min_orf_aa: int = 20
    contaminants: Path | None = None  # library FASTA from `giap setup`; None: --no-contaminants


@dataclass(frozen=True)
class OrfRecord:
    """A six-frame ORF or ``nterm_extension`` entry; 0-based half-open, unwrapped."""

    cls: str
    contig: str
    start: int
    end: int
    strand: str
    aa: str
    locus_tag: str | None = None

    @property
    def db_id(self) -> str:
        sfx = "p" if self.strand == "+" else "m"
        return f"{ORF_ID_PREFIX}{self.cls}_{self.contig}_{self.start + 1}_{self.end}_{sfx}"


@dataclass(frozen=True)
class MakedbResult:
    fasta: Path
    id_map: Path
    manifest: Path
    counts: Mapping[str, int]


# ---------------------------------------------------------------- six-frame ORFs


def _strand_runs(t: str, circular: bool, min_aa: int) -> Iterator[tuple[int, int, str]]:
    """Stop/X-delimited segments >= ``min_aa`` of the three frames of ``t`` (strand coords)."""
    n = len(t)
    ext = t + t + t[:6] if circular else t
    for o in range(3):
        aa = translate11_text(ext, o)
        codes = np.frombuffer(aa.encode("ascii"), dtype=np.uint8)
        bnd = np.flatnonzero((codes == ord("*")) | (codes == ord("X")))
        if circular:
            starts, ends = bnd[:-1] + 1, bnd[1:]
            keep = (o + 3 * bnd[:-1] < n) & (3 * (ends - starts) <= n)
        else:
            full = np.concatenate(([-1], bnd, [len(aa)]))
            starts, ends = full[:-1] + 1, full[1:]
            keep = np.ones(len(starts), dtype=bool)
        keep &= ends - starts >= min_aa
        for a, b in zip(starts[keep].tolist(), ends[keep].tolist(), strict=True):
            s = o + 3 * a
            if s >= n:
                s -= n
            yield s, s + 3 * (b - a), aa[a:b]


def _raw_orfs(contig: str, seq: str, rep: Replicon, min_aa: int) -> Iterator[OrfRecord]:
    circ = rep.topology == "circular"
    n = len(seq)
    for s, e, aa in _strand_runs(seq, circ, min_aa):
        yield OrfRecord("", contig, s, e, "+", aa)
    for a, b, aa in _strand_runs(seq.translate(_COMP)[::-1], circ, min_aa):
        s, e = n - b, n - a
        if s < 0:
            s, e = s + n, e + n
        yield OrfRecord("", contig, s, e, "-", aa)


def _in_frame_coding(orf: OrfRecord, gm: GeneModel, rep: Replicon) -> bool:
    iv = Interval(orf.contig, orf.start, orf.end, orf.strand)  # type: ignore[arg-type]
    for g in gm.overlapping(iv):
        if not isinstance(g, Gene) or g.feature_type not in _CODING:
            continue
        gi = g.interval
        if gi.strand != orf.strand:
            continue
        for d in _shifts(rep):
            a, b = orf.start + d, orf.end + d
            if max(a, gi.start) < min(b, gi.end) and in_frame(a, b, gi):
                return True
    return False


def _upstream(seq: str, g: Gene, rep: Replicon) -> tuple[int, int] | None:
    """In-frame upstream part of ``g`` up to the first stop (unwrapped plus coordinates); None
    when the walk runs off a linear contig, meets an ``X`` or finds no stop within one turn."""
    n, circ, gi = rep.length, rep.topology == "circular", g.interval
    step = -3 if gi.strand == "+" else 3
    q0 = gi.start - 3 if gi.strand == "+" else gi.end
    q = q0
    while True:
        if not circ and not 0 <= q <= n - 3:
            return None
        if abs(q - q0) >= n:
            return None
        cod = interval_nt(seq, Interval(gi.contig, q % n, q % n + 3, gi.strand))
        aa = translate11_text(cod)
        if aa == "*":
            break
        if aa == "X":
            return None
        q += step
    a, b = (q + 3, gi.start) if gi.strand == "+" else (gi.end, q)
    if a < 0:
        a, b = a + n, b + n
    elif a >= n:
        a, b = a - n, b - n
    return (a, b) if b > a else None


def _nterm_entries(gm: GeneModel, genome: Mapping[str, str]) -> Iterator[OrfRecord]:
    for g in gm.genes:
        if g.feature_type not in _CODING or g.interval.strand not in ("+", "-"):
            continue
        rep, seq = gm.replicons[g.interval.contig], genome[g.interval.contig]
        up = _upstream(seq, g, rep)
        if up is None:
            continue
        a, b = up
        strand = g.interval.strand
        up_aa = translate11_text(interval_nt(seq, Interval(g.interval.contig, a, b, strand)))
        if len(up_aa) < NTERM_MIN_UPSTREAM_AA:
            continue
        gi = g.interval
        span = 3 * min(NTERM_PROTEIN_AA, (gi.end - gi.start) // 3 - 1)  # junction codons
        if strand == "+":
            s, e, js, je = a, b + span, gi.start, gi.start + span
        else:
            s, e, js, je = a - span, b, gi.end - span, gi.end
            if s < 0:
                s, e = s + rep.length, e + rep.length
        junction = translate11_text(_strand_nt(seq, js, je, strand, rep.length))
        yield OrfRecord("nterm_extension", gi.contig, s, e, strand, up_aa + junction, g.locus_tag)


def _strand_nt(seq: str, start: int, end: int, strand: str, n: int) -> str:
    """Nucleotides of unwrapped ``[start, end)`` (modulo ``n``) read on ``strand``."""
    if start < 0:
        start, end = start + n, end + n
    nt = (seq + seq)[start:end] if end > n else seq[start:end]
    return nt.translate(_COMP)[::-1] if strand == "-" else nt


def six_frame_orfs(
    index: SixFrameIndex, gm: GeneModel, genome: Mapping[str, str], min_aa: int
) -> list[OrfRecord]:
    """Class-tagged six-frame ORFs and ``nterm_extension`` entries, sorted (contig, start, end,
    strand). ``index`` (M3-01) supplies the replicons; sequences come from ``genome``."""
    reps = index.replicons
    classifier = LocusClassifier(gm, genome)
    out: list[OrfRecord] = list(_nterm_entries(gm, genome))
    for contig, seq in genome.items():
        rep = reps[contig]
        for orf in _raw_orfs(contig, seq, rep, min_aa):
            if _in_frame_coding(orf, gm, rep):
                continue
            iv = Interval(contig, orf.start, orf.end, orf.strand)  # type: ignore[arg-type]
            cls, _ = classifier.classify(iv)
            if cls in ORF_CLASSES:
                out.append(OrfRecord(cls, contig, orf.start, orf.end, orf.strand, orf.aa))
    order = {c: k for k, c in enumerate(reps)}
    out.sort(key=lambda o: (order[o.contig], o.start, o.end, o.strand == "-"))
    return out


# ---------------------------------------------------------------- database


def _protein(g: Gene) -> str:
    """Annotated protein (Bakta ``aa``), else the table-11 translation with an initiator M."""
    if g.aa:
        return g.aa
    return ("M" + translate11_text((g.nt or "")[3:])).rstrip("*")


def _attr(attrs: Mapping[str, AttrValue], key: str) -> str | None:
    v = attrs.get(key)
    if v is None or isinstance(v, str):
        return v
    return v[0] if v else None


def _gene_product(gm: GeneModel, tag: str | None) -> str:
    try:
        return gm.by_locus(tag).product if tag else ""
    except KeyError:
        return ""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _declared_lengths(path: Path) -> dict[str, int | None]:
    """Contig -> the length the Bakta JSON itself states; ``None`` when it states none.

    A record with neither ``length`` nor ``nt`` says nothing about its own size, so a ``--genome``
    FASTA cannot be checked against it (lead decision 13): ``build_strain_db`` refuses such an
    annotation instead of comparing the replicon length with the FASTA it was derived from.
    """
    data = read_json(path)
    out: dict[str, int | None] = {}
    if not isinstance(data, Mapping):
        return out
    for rec in data.get("sequences", data.get("contigs")) or []:
        if not isinstance(rec, Mapping) or not rec.get("id"):
            continue
        length, nt = rec.get("length"), rec.get("nt") or rec.get("sequence")
        if isinstance(length, int) and not isinstance(length, bool):
            out[str(rec["id"])] = length
        else:
            out[str(rec["id"])] = len(nt) if isinstance(nt, str) else None
    return out


Row = tuple[str, str, str | None, str | None, str | None, int | None, int | None, str | None]


def build_strain_db(opts: MakedbOptions) -> MakedbResult:
    """Build ``strain_db.faa``, ``db_id_map.tsv`` and ``makedb_manifest.json`` in ``opts.outdir``.

    Raises ValueError on an annotation without sequences or features (e.g. ``{}``), never writing
    an empty database."""
    if opts.min_orf_aa < 1:
        raise ValueError(f"min_orf_aa must be >= 1, got {opts.min_orf_aa}")
    fasta = (
        {r.id: r.sequence for r in read_genome_fasta(opts.genome)}
        if opts.genome is not None
        else None
    )
    # The genome also replaces nucleotides a recorded Bakta JSON does not carry (M3-24 decision 8).
    ann = parse_bakta_json(opts.annotation, genome=fasta)
    gm = ann.model
    if not gm.replicons:
        raise ValueError(f"{opts.annotation}: annotation has no sequences or features")
    genome = dict(ann.sequences)
    if fasta is not None:
        # The JSON's own lengths, never the replicon lengths: without `nt` those come from this
        # very FASTA, so comparing them would compare a value with itself (lead decision 13).
        declared = _declared_lengths(opts.annotation)
        for cid in gm.replicons:
            length = declared.get(cid)
            if length is None:
                raise ValueError(
                    f"{opts.annotation}: sequence {cid!r} declares no length and carries no "
                    f"nucleotides, so {opts.genome} cannot be checked against it"
                )
            if cid not in fasta or len(fasta[cid]) != length:
                raise ValueError(f"{opts.genome}: contig {cid!r} missing or of the wrong length")
            genome[cid] = fasta[cid]
    records: list[tuple[str, str]] = []
    rows: list[Row] = []

    def add(header: str, seq: str, row: Row) -> None:
        records.append((header, seq))
        rows.append(row)

    for g in gm.genes:
        if g.feature_type in _CODING:
            iv = g.interval
            row: Row = (g.locus_tag, g.feature_type, g.locus_tag, None, iv.contig, iv.start + 1,
                        iv.end, iv.strand)  # fmt: skip
            add(f"{g.locus_tag} {g.product or 'hypothetical protein'}", _protein(g), row)
    if opts.patch is not None:
        taken = {r[0] for r in rows}
        for f in read_gff3(opts.patch, gm.replicons):
            if f.type != "CDS":
                continue
            iv = f.interval
            tag = _attr(f.attributes, "locus_tag")
            pid = _attr(f.attributes, "ID") or tag or f"PATCH_{iv.contig}_{iv.start + 1}"
            rid = f"{pid}_patch" if pid in taken else pid
            taken.add(rid)
            nt = interval_nt(genome[iv.contig], iv)
            seq = ("M" + translate11_text(nt[3:])).removesuffix("*")
            if "*" in seq:
                raise ValueError(f"{opts.patch}: patch CDS {rid} has an internal stop codon")
            product = _attr(f.attributes, "product") or _gene_product(gm, tag)
            row = (rid, "patch", tag, None, iv.contig, iv.start + 1, iv.end, iv.strand)
            add(f"{rid} {product or 'hypothetical protein'}", seq, row)
    if opts.six_frame_orfs:
        from giap.mapping.sixframe import build_sixframe

        index = build_sixframe(genome, gm.replicons, wrap_nt=0)
        for o in six_frame_orfs(index, gm, genome, opts.min_orf_aa):
            header = o.db_id
            if o.cls == "nterm_extension":
                header += f" upstream of {o.locus_tag} + its first {NTERM_PROTEIN_AA} codons"
            row = (o.db_id, "orf", o.locus_tag, o.cls, o.contig, o.start + 1, o.end, o.strand)
            add(header, o.aa, row)
    if opts.contaminants is not None:
        for r in read_protein_fasta(opts.contaminants):
            row = (f"CONT_{r.id}", "contaminant", None, None, None, None, None, None)
            add(f"CONT_{r.id}", r.sequence.rstrip("*"), row)

    out = opts.outdir
    out.mkdir(parents=True, exist_ok=True)
    fasta_path, map_path, manifest_path = (
        out / "strain_db.faa",
        out / "db_id_map.tsv",
        out / "makedb_manifest.json",
    )
    write_protein_fasta(fasta_path, records)
    write_tsv(map_path, ID_MAP_COLUMNS, rows)
    counts = {k: sum(r[1] == k for r in rows) for k in ("cds", "sorf", "patch", "orf",
                                                        "contaminant")}  # fmt: skip
    inputs = {
        name: {"path": str(p), "sha256": _sha256(p)}
        for name, p in (
            ("annotation", opts.annotation),
            ("genome", opts.genome),
            ("patch", opts.patch),
            ("contaminants", opts.contaminants),
        )
        if p is not None
    }
    manifest = {
        "giap_version": giap.__version__,
        "date": datetime.date.today().isoformat(),
        "inputs": inputs,
        "options": {
            "six_frame_orfs": opts.six_frame_orfs,
            "min_orf_aa": opts.min_orf_aa,
            "contaminants": opts.contaminants is not None,
        },
        "counts": counts,
        "decoys": False,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return MakedbResult(fasta_path, map_path, manifest_path, counts)
