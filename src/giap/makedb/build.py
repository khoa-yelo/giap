"""`giap makedb`: optional strain search database (issue M3-23; GIAP_SPEC.md §7.9.2, §5.1).

Writes ``strain_db.faa`` (Bakta CDS/sORF proteins, patched proteins, class-tagged six-frame ORFs,
contaminants; wrapped at 60, no ``*``, no decoys), ``id_map.tsv`` (one row per record, 1-based
closed coordinates, unwrapped across the origin) and ``makedb_manifest.json``.

Six-frame ORFs are the maximal stop-to-stop segments of each frame, split at ``X`` (lead decision
3), translated codon by codon from the genome (no I -> L, no start rule) with the M3-01 kernel.
An ORF in frame with an annotated CDS/sORF is dropped; each such gene with >= 7 stop-free
upstream codons gets one ``nterm_extension`` entry = upstream part + the first 30 codons of the
annotated protein, whose header/id_map interval covers both (lead decision 2). The other ORFs
take the M3-14 ``classify_locus`` class; ``rna_overlap``/``pseudogene`` ORFs are omitted.
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
from giap.io.fasta import read_genome_fasta
from giap.io.gff3 import AttrValue, read_gff3
from giap.io.protein_fasta import read_protein_fasta, write_protein_fasta
from giap.io.tsv import write_tsv
from giap.mapping.loci import _build_index, _classify, _in_frame, _shifts
from giap.mapping.sixframe import _BASES, _NT_LUT, _TABLE11_TCAG, SixFrameIndex, _translate_kernel
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel

__all__ = ["MakedbOptions", "MakedbResult", "OrfRecord", "build_strain_db", "six_frame_orfs"]

ORF_CLASSES = frozenset({"intergenic", "cds_outframe", "antisense"})
NTERM_PROTEIN_AA = 30
NTERM_MIN_UPSTREAM_AA = 7
ID_MAP_COLUMNS = ("db_id", "kind", "locus_tag", "class", "contig", "start", "end", "strand")
_CODING = frozenset({"cds", "sorf"})
_AMBIG = 20  # code written by the M3-01 kernel for a codon with a non-ACGT base
_STOP = ord("*")
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _ascii_codons() -> np.ndarray:
    table = np.empty(64, dtype=np.uint8)
    for i, b0 in enumerate("TCAG"):
        for j, b1 in enumerate("TCAG"):
            for k, b2 in enumerate("TCAG"):
                idx = 16 * _BASES.index(b0) + 4 * _BASES.index(b1) + _BASES.index(b2)
                table[idx] = ord(_TABLE11_TCAG[16 * i + 4 * j + k])
    return table


_CODON_ASCII = _ascii_codons()


def _translate(nt: str, offset: int = 0) -> np.ndarray:
    """Table 11 codon by codon as ASCII bytes; a non-ACGT codon -> code ``_AMBIG``."""
    raw = np.frombuffer(nt.encode("ascii", errors="replace"), dtype=np.uint8)
    out: np.ndarray = _translate_kernel(raw, _NT_LUT, _CODON_ASCII, offset, False)
    return out


def _text(codes: np.ndarray) -> str:
    return np.where(codes == _AMBIG, ord("X"), codes).astype(np.uint8).tobytes().decode("ascii")


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
        return f"GIAPORF_{self.cls}_{self.contig}_{self.start + 1}_{self.end}_{sfx}"


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
        aa = _translate(ext, o)
        bnd = np.flatnonzero((aa == _STOP) | (aa == _AMBIG))
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
            yield s, s + 3 * (b - a), _text(aa[a:b])


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
            if max(a, gi.start) < min(b, gi.end) and _in_frame(a, b, gi):
                return True
    return False


def _upstream(seq: str, g: Gene, rep: Replicon) -> tuple[int, int] | None:
    """In-frame stop/X-free upstream part of ``g`` (unwrapped plus coordinates)."""
    n, circ, gi = rep.length, rep.topology == "circular", g.interval
    step = -3 if gi.strand == "+" else 3
    q = gi.start - 3 if gi.strand == "+" else gi.end
    while circ or 0 <= q <= n - 3:
        if abs(q - (gi.start - 3 if gi.strand == "+" else gi.end)) >= n:
            return None
        cod = interval_nt(seq, Interval(gi.contig, q % n, q % n + 3, gi.strand))
        if _translate(cod)[0] in (_STOP, _AMBIG):
            break
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
        up_aa = _text(_translate(interval_nt(seq, Interval(g.interval.contig, a, b, strand))))
        if len(up_aa) < NTERM_MIN_UPSTREAM_AA:
            continue
        prot = _protein(g)[:NTERM_PROTEIN_AA]
        span = 3 * len(prot)
        if strand == "+":
            s, e = a, b + span
        else:
            s, e = a - span, b
            if s < 0:
                s, e = s + rep.length, e + rep.length
        yield OrfRecord(
            "nterm_extension", g.interval.contig, s, e, strand, up_aa + prot, g.locus_tag
        )


def six_frame_orfs(
    index: SixFrameIndex, gm: GeneModel, genome: Mapping[str, str], min_aa: int
) -> list[OrfRecord]:
    """Class-tagged six-frame ORFs and ``nterm_extension`` entries, sorted (contig, start, end,
    strand). ``index`` (M3-01) supplies the replicons; sequences come from ``genome``."""
    reps = index.replicons
    loci_index = _build_index(gm)
    out: list[OrfRecord] = list(_nterm_entries(gm, genome))
    for contig, seq in genome.items():
        rep = reps[contig]
        for orf in _raw_orfs(contig, seq, rep, min_aa):
            if _in_frame_coding(orf, gm, rep):
                continue
            cls, _ = _classify(contig, orf.start, orf.end, orf.strand, gm, genome, loci_index)
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
    return ("M" + _text(_translate((g.nt or "")[3:]))).rstrip("*")


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


Row = tuple[str, str, str | None, str | None, str | None, int | None, int | None, str | None]


def build_strain_db(opts: MakedbOptions) -> MakedbResult:
    """Build ``strain_db.faa``, ``id_map.tsv`` and ``makedb_manifest.json`` in ``opts.outdir``."""
    if opts.min_orf_aa < 1:
        raise ValueError(f"min_orf_aa must be >= 1, got {opts.min_orf_aa}")
    ann = parse_bakta_json(opts.annotation)
    gm = ann.model
    genome = dict(ann.sequences)
    if opts.genome is not None:
        fasta = {r.id: r.sequence for r in read_genome_fasta(opts.genome)}
        for cid, rep in gm.replicons.items():
            if cid not in fasta or len(fasta[cid]) != rep.length:
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
            seq = ("M" + _text(_translate(nt[3:]))).rstrip("*")
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
        out / "id_map.tsv",
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
