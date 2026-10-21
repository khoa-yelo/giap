r"""K-12 pair and terminator truth from RegulonDB (GIAP_SPEC.md §15.4 steps 1-2; issue M6-11).

Pure parsing and labelling: everything here takes paths or already-parsed objects and returns
values, so the module is offline-testable (the network lives in ``validation/fetch_regulondb.py``).

§15.4 step 2, implemented literally: "A same-strand adjacent pair is positive when both genes lie
in one RegulonDB transcription unit. It is negative when they lie in different TUs and no TU
contains both. Everything else is ambiguous and excluded."

Adjacency comes from the recorded K-12 Bakta gene model (consecutive genes on a contig with equal
strand). RegulonDB names genes, so the gene set supplies coordinates and the names are mapped onto
Bakta locus tags by coordinates (within a tolerance) and then by gene name; the mapped share is
reported, and a share below 0.8 is escalated rather than worked around.

Usage::

    python validation/regulondb_truth.py \
        --regulondb-dir $GIAP_DATA_DIR/regulondb/14.5.0 \
        --bakta-json tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz \
        --genome tests/data/genomes/k12_mg1655.fna.gz \
        --out-dir $GIAP_DATA_DIR/regulondb/14.5.0/truth

``--genome`` is needed whenever the Bakta JSON omits ``sequences[].nt``, as the recorded one does.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from bisect import bisect_left, bisect_right
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from giap.annotate.bakta import parse_bakta_json
from giap.annotate.operons import adjacent_pairs
from giap.io.fasta import read_genome_fasta
from giap.model.coords import Interval, Replicon, Strand
from giap.model.genes import Gene, GeneModel

#: RegulonDB gene coordinates agree with Bakta exactly on K-12; a few bp of slack absorbs
#: annotation differences at the start codon without letting a neighbour win.
DEFAULT_TOLERANCE_BP = 30

#: Below this mapped share the result is escalated, not worked around (spec example 5).
MIN_MAPPED_SHARE = 0.8

#: Bakta non-coding features a RegulonDB RNA gene may map onto (lead decision 4).
RNA_FEATURE_TYPES = frozenset({"ncrna", "trna", "rrna"})

#: The three downloads ``build`` reads, each verified against ``provenance.json`` first.
DATASETS = ("TUSet", "TerminatorSet", "GeneProductSet")

PAIR_COLUMNS = ("gene_a", "gene_b", "label", "tu_ids")
TERMINATOR_COLUMNS = ("contig", "start", "end", "strand", "source_id")

_STRANDS: Mapping[str, Strand] = {"forward": "+", "reverse": "-", "+": "+", "-": "-"}


class RegulonDBFormatError(ValueError):
    """A RegulonDB download does not have the expected columns."""


class ProvenanceError(ValueError):
    """A download does not match its ``provenance.json`` record, or the record is missing."""


# ------------------------------------------------------------------ parsing


@dataclass(frozen=True, slots=True)
class RegulonDBGene:
    """A gene row of ``GeneProductSet`` with its interval in GIAP coordinates."""

    gene_id: str
    name: str
    bnumber: str
    interval: Interval


@dataclass(frozen=True, slots=True)
class TranscriptionUnit:
    """A row of ``TUSet``; ``genes`` holds RegulonDB gene names, in no guaranteed order."""

    tu_id: str
    name: str
    operon_id: str
    operon_name: str
    genes: tuple[str, ...]
    terminator_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RegulonDBTerminator:
    """A row of ``TerminatorSet`` with its interval in GIAP coordinates."""

    terminator_id: str
    interval: Interval
    terminator_type: str
    related_tus: tuple[str, ...]


def parse_rif(text: str) -> tuple[tuple[str, ...], tuple[dict[str, str], ...]]:
    """Parse a RegulonDB "rif-version 1" table into its column names and rows.

    Comment lines (``#``) and blank lines are dropped; the first remaining line is the header,
    whose fields look like ``3)operonId``. Short rows are padded with empty strings.
    """
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise RegulonDBFormatError("empty RegulonDB table")
    header = tuple(_column_name(f) for f in lines[0].split("\t"))
    rows = []
    for line in lines[1:]:
        fields = [f.strip() for f in line.split("\t")]
        if len(fields) > len(header):
            raise RegulonDBFormatError(
                f"row has {len(fields)} fields but the header has {len(header)}: {line!r:.120}"
            )
        fields += [""] * (len(header) - len(fields))
        rows.append(dict(zip(header, fields, strict=True)))
    return header, tuple(rows)


def _column_name(field: str) -> str:
    """``"3)operonId"`` -> ``"operonId"``; a field without the ``N)`` prefix is kept as written."""
    head, sep, tail = field.strip().partition(")")
    return tail.strip() if sep and head.strip().isdigit() else field.strip()


def _require(header: Sequence[str], needed: Iterable[str], what: str) -> None:
    missing = [c for c in needed if c not in header]
    if missing:
        raise RegulonDBFormatError(f"{what}: missing column(s) {', '.join(missing)}")


def _strand(value: str, what: str) -> Strand:
    try:
        return _STRANDS[value.strip().lower()]
    except KeyError:
        raise RegulonDBFormatError(f"{what}: unknown strand {value!r}") from None


def _interval(contig: str, left: str, right: str, strand: Strand, what: str) -> Interval:
    """RegulonDB is 1-based closed; GIAP is 0-based half-open (spec §6.5)."""
    try:
        start, end = int(left), int(right)
    except ValueError:
        raise RegulonDBFormatError(f"{what}: non-integer coordinates {left!r}, {right!r}") from None
    if start < 1 or end < start:
        raise RegulonDBFormatError(f"{what}: invalid 1-based interval {start}-{end}")
    return Interval(contig=contig, start=start - 1, end=end, strand=strand)


def _present(value: str) -> bool:
    """RegulonDB writes the literal ``None`` for an absent value as well as an empty field."""
    return bool(value.strip()) and value.strip() != "None"


def _split(value: str, sep: str) -> tuple[str, ...]:
    return tuple(p.strip() for p in value.split(sep) if p.strip())


def parse_gene_set(text: str, contig: str) -> tuple[RegulonDBGene, ...]:
    """Parse ``GeneProductSet``; rows without coordinates are dropped (they cannot be mapped)."""
    header, rows = parse_rif(text)
    _require(
        header, ("geneId", "geneName", "leftEndPos", "rightEndPos", "strand"), "GeneProductSet"
    )
    out = []
    for r in rows:
        if not all(_present(r[c]) for c in ("leftEndPos", "rightEndPos", "strand")):
            continue
        what = f"GeneProductSet {r['geneId']}"
        out.append(
            RegulonDBGene(
                gene_id=r["geneId"],
                name=r["geneName"],
                bnumber=r.get("bnumber", ""),
                interval=_interval(
                    contig, r["leftEndPos"], r["rightEndPos"], _strand(r["strand"], what), what
                ),
            )
        )
    return tuple(out)


def parse_tu_set(text: str) -> tuple[TranscriptionUnit, ...]:
    """Parse ``TUSet``; ``tuGenes`` is a ``;``-separated list of RegulonDB gene names."""
    header, rows = parse_rif(text)
    _require(header, ("id", "name", "tuGenes"), "TUSet")
    return tuple(
        TranscriptionUnit(
            tu_id=r["id"],
            name=r["name"],
            operon_id=r.get("operonId", ""),
            operon_name=r.get("operonName", ""),
            genes=_split(r["tuGenes"], ";"),
            terminator_ids=_split(r.get("terminatorIds", ""), ";"),
        )
        for r in rows
        if r["id"]
    )


def parse_terminator_set(text: str, contig: str) -> tuple[RegulonDBTerminator, ...]:
    """Parse ``TerminatorSet`` into intervals on ``contig``."""
    header, rows = parse_rif(text)
    _require(header, ("id", "leftEndPos", "rightEndPos", "strand"), "TerminatorSet")
    out = []
    for r in rows:
        if not all(_present(r[c]) for c in ("id", "leftEndPos", "rightEndPos", "strand")):
            continue
        what = f"TerminatorSet {r['id']}"
        out.append(
            RegulonDBTerminator(
                terminator_id=r["id"],
                interval=_interval(
                    contig, r["leftEndPos"], r["rightEndPos"], _strand(r["strand"], what), what
                ),
                terminator_type=r.get("type", ""),
                related_tus=_split(r.get("relatedTus", ""), "|"),
            )
        )
    return tuple(out)


# ------------------------------------------------------------------ mapping


@dataclass(frozen=True, slots=True)
class GeneMapping:
    """RegulonDB genes resolved onto Bakta locus tags (or RNA feature ids), with the shares.

    Keys are RegulonDB gene names, except that a name RegulonDB uses more than once is keyed by
    each copy's RegulonDB gene id, so every copy is kept (lead decision 5); ``names`` maps such a
    key back to its gene name (a key absent from ``names`` is the name itself).
    """

    locus_by_name: Mapping[str, str]
    how: Mapping[str, str]  # key -> "coordinates" | "name"
    unmapped: tuple[str, ...]  # keys of RegulonDB genes with no Bakta counterpart
    rna: frozenset[str] = frozenset()  # keys mapped onto a Bakta ncRNA / tRNA / rRNA feature
    names: Mapping[str, str] = field(default_factory=dict)

    @property
    def n_mapped(self) -> int:
        return len(self.locus_by_name)

    @property
    def mapped_share(self) -> float:
        """Share of all RegulonDB genes (coding and RNA) that found a Bakta counterpart."""
        total = self.n_mapped + len(self.unmapped)
        return self.n_mapped / total if total else 0.0

    @property
    def n_mapped_cds(self) -> int:
        return self.n_mapped - len(self.rna)

    @property
    def cds_mapped_share(self) -> float:
        """Share mapped onto a Bakta CDS among the RegulonDB genes not matched to an RNA feature."""
        total = self.n_mapped_cds + len(self.unmapped)
        return self.n_mapped_cds / total if total else 0.0

    def name_of(self, key: str) -> str:
        return self.names.get(key, key)


def chromosome_contig(model: GeneModel) -> str:
    """The contig RegulonDB coordinates refer to: the chromosome, or the only replicon."""
    chromosomes = [r.id for r in model.replicons.values() if r.kind == "chromosome"]
    if len(chromosomes) == 1:
        return chromosomes[0]
    if len(model.replicons) == 1:
        return next(iter(model.replicons))
    raise ValueError(
        f"cannot pick the RegulonDB contig among {sorted(model.replicons)}; pass --contig"
    )


@dataclass(frozen=True, slots=True)
class _Target:
    """A Bakta gene or RNA feature a RegulonDB gene may map onto."""

    ident: str
    interval: Interval
    name: str | None
    rna: bool


def _targets(model: GeneModel, contig: str) -> list[_Target]:
    genes = [
        _Target(g.locus_tag, g.interval, g.gene, False)
        for g in model.genes
        if g.interval.contig == contig
    ]
    rnas = [
        _Target(f.id, f.interval, f.gene, True)
        for f in model.features
        if f.interval.contig == contig and f.type in RNA_FEATURE_TYPES
    ]
    return sorted(genes + rnas, key=lambda t: (t.interval.start, t.interval.end, t.ident))


def _keys(rdb_genes: Sequence[RegulonDBGene]) -> list[str]:
    """Each gene's key: its name, or its RegulonDB gene id when the name is not unique."""
    counts = Counter(g.name for g in rdb_genes)
    return [g.name if g.name and counts[g.name] == 1 else g.gene_id for g in rdb_genes]


def _overlaps(a: Interval, b: Interval) -> bool:
    return a.start < b.end and b.start < a.end


def map_genes(
    rdb_genes: Sequence[RegulonDBGene],
    model: GeneModel,
    contig: str,
    tolerance_bp: int = DEFAULT_TOLERANCE_BP,
) -> GeneMapping:
    """Map RegulonDB genes onto Bakta CDS and RNA features, by coordinates first, then by name.

    A coordinate match needs the same strand and both ends within ``tolerance_bp`` (inclusive);
    the smallest total offset wins, globally and deterministically, and each Bakta gene or feature
    is claimed at most once. RNA genes match Bakta ncRNA / tRNA / rRNA features the same way
    (lead decision 4). The name fallback needs a unique name on both sides, the same strand and at
    least 1 bp of overlap (lead decision 5).
    """
    targets = _targets(model, contig)
    starts = [t.interval.start for t in targets]
    keys = _keys(rdb_genes)
    candidates: list[tuple[int, str, str]] = []
    for key, rg in zip(keys, rdb_genes, strict=True):
        lo, hi = rg.interval.start - tolerance_bp, rg.interval.start + tolerance_bp
        for t in targets[bisect_left(starts, lo) : bisect_right(starts, hi)]:
            if t.interval.strand != rg.interval.strand:
                continue
            d_end = abs(t.interval.end - rg.interval.end)
            if d_end > tolerance_bp:
                continue
            candidates.append((abs(t.interval.start - rg.interval.start) + d_end, key, t.ident))

    locus_by_name: dict[str, str] = {}
    how: dict[str, str] = {}
    claimed: set[str] = set()
    for _offset, key, ident in sorted(candidates):
        if key in locus_by_name or ident in claimed:
            continue
        locus_by_name[key] = ident
        how[key] = "coordinates"
        claimed.add(ident)

    by_name: dict[str, list[_Target]] = {}
    for t in targets:
        if t.name:
            by_name.setdefault(t.name.casefold(), []).append(t)
    for key, rg in sorted(zip(keys, rdb_genes, strict=True), key=lambda kg: kg[0]):
        if key in locus_by_name or key != rg.name:
            continue  # resolved already, or a duplicated name (no unambiguous name match)
        same_name = by_name.get(rg.name.casefold(), [])
        if len(same_name) != 1:
            continue
        (t,) = same_name
        if (
            t.ident in claimed
            or t.interval.strand != rg.interval.strand
            or not _overlaps(t.interval, rg.interval)
        ):
            continue
        locus_by_name[key] = t.ident
        how[key] = "name"
        claimed.add(t.ident)

    rna_ids = {t.ident for t in targets if t.rna}
    return GeneMapping(
        locus_by_name=locus_by_name,
        how=how,
        unmapped=tuple(sorted(k for k in set(keys) if k not in locus_by_name)),
        rna=frozenset(k for k, ident in locus_by_name.items() if ident in rna_ids),
        names={k: g.name for k, g in zip(keys, rdb_genes, strict=True) if k != g.name},
    )


# ----------------------------------------------------------------- the truth


@dataclass(frozen=True, slots=True)
class PairRow:
    """One same-strand adjacent pair and its §15.4 label."""

    gene_a: str
    gene_b: str
    label: str  # "positive" | "negative" | "ambiguous"
    tu_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class TerminatorRow:
    """One RegulonDB terminator in GIAP coordinates."""

    contig: str
    start: int
    end: int
    strand: Strand
    source_id: str


def adjacent_same_strand_pairs(model: GeneModel, contig: str) -> tuple[tuple[Gene, Gene], ...]:
    """The pipeline's adjacent pairs (M2-12, ``pairs_raw`` and UniOP) restricted to ``contig``.

    ``giap.annotate.operons.adjacent_pairs`` keeps same-strand pairs only and adds the wrap pair
    of a circular replicon (lead decision 2).
    """
    return tuple(
        (a, b)
        for a, b in adjacent_pairs(model.genes, model.replicons)
        if a.interval.contig == contig and a.interval.strand == b.interval.strand
    )


def tus_by_locus(
    tus: Sequence[TranscriptionUnit], mapping: GeneMapping
) -> dict[str, frozenset[str]]:
    """Locus tag -> the ids of the transcription units that contain it.

    A TU names its genes; a name RegulonDB uses more than once resolves to every mapped copy.
    """
    loci_by_name: dict[str, list[str]] = {}
    for key, locus in mapping.locus_by_name.items():
        loci_by_name.setdefault(mapping.name_of(key), []).append(locus)
    out: dict[str, set[str]] = {}
    for tu in tus:
        for name in tu.genes:
            for locus in loci_by_name.get(name, ()):
                out.setdefault(locus, set()).add(tu.tu_id)
    return {locus: frozenset(ids) for locus, ids in out.items()}


def _gap(a: Gene, b: Gene, replicon: Replicon) -> tuple[int, int] | None:
    """The intergenic gap from ``a`` to ``b`` on the unrolled line, or None when they touch.

    For the wrap pair of a circular replicon (``b`` starts before ``a``) the gap runs past the
    origin, so ``b``'s start is shifted by the replicon length.
    """
    lo, hi = a.interval.end, b.interval.start
    if replicon.topology == "circular" and b.interval.start < a.interval.start:
        hi += replicon.length
    return (lo, hi) if hi > lo else None


def _unmapped_between(
    a: Gene,
    b: Gene,
    unmapped_intervals: Sequence[tuple[int, int]],
    replicon: Replicon,
) -> bool:
    """True when a RegulonDB gene without a Bakta counterpart lies in the gap from ``a`` to ``b``.

    Strand-independent (lead decision 3) and circular-aware (lead decision 2): each interval is
    also tried one replicon length further on, which catches the part of a wrap gap past the
    origin.
    """
    gap = _gap(a, b, replicon)
    if gap is None:
        return False
    lo, hi = gap
    shifts = (0, replicon.length) if replicon.topology == "circular" else (0,)
    return any(start + k < hi and end + k > lo for start, end in unmapped_intervals for k in shifts)


def label_pair(
    a: Gene,
    b: Gene,
    tu_ids: Mapping[str, frozenset[str]],
    mapped_loci: frozenset[str],
    blocked: bool,
) -> PairRow:
    """Label one adjacent same-strand pair as GIAP_SPEC.md §15.4 step 2 says (lead decision 3).

    In order: **positive** when some RegulonDB TU holds both genes — even when a RegulonDB-only
    gene lies between them, and two overlapping TUs still count; otherwise **ambiguous** when
    either gene has no RegulonDB counterpart or a RegulonDB-only gene lies between them
    (``blocked``: the different-TU evidence is then unreliable); otherwise **negative** when each
    gene is in at least one TU; otherwise **ambiguous** (a gene in no TU at all).
    """
    a_tus = tu_ids.get(a.locus_tag, frozenset())
    b_tus = tu_ids.get(b.locus_tag, frozenset())
    both = a_tus & b_tus
    if both:
        return PairRow(a.locus_tag, b.locus_tag, "positive", tuple(sorted(both)))
    union = tuple(sorted(a_tus | b_tus))
    unmapped = a.locus_tag not in mapped_loci or b.locus_tag not in mapped_loci
    if unmapped or blocked or not (a_tus and b_tus):
        return PairRow(a.locus_tag, b.locus_tag, "ambiguous", union)
    return PairRow(a.locus_tag, b.locus_tag, "negative", union)


def pair_truth(
    model: GeneModel,
    tus: Sequence[TranscriptionUnit],
    rdb_genes: Sequence[RegulonDBGene],
    mapping: GeneMapping,
    contig: str,
) -> tuple[PairRow, ...]:
    """All same-strand adjacent pairs of ``contig`` (wrap pair included) with their labels."""
    by_locus = tus_by_locus(tus, mapping)
    mapped_loci = frozenset(mapping.locus_by_name.values())
    unmapped_keys = frozenset(mapping.unmapped)
    keys = _keys(rdb_genes)
    unmapped_intervals = sorted(
        (g.interval.start, g.interval.end)
        for key, g in zip(keys, rdb_genes, strict=True)
        if key in unmapped_keys
    )
    replicon = model.replicons[contig]
    return tuple(
        label_pair(
            a, b, by_locus, mapped_loci, _unmapped_between(a, b, unmapped_intervals, replicon)
        )
        for a, b in adjacent_same_strand_pairs(model, contig)
    )


def wrap_pairs(model: GeneModel, pairs: Sequence[PairRow]) -> list[dict[str, str]]:
    """The pairs that span the origin (the second gene starts before the first), for the report."""
    by_locus = {g.locus_tag: g for g in model.genes}
    return [
        {"gene_a": p.gene_a, "gene_b": p.gene_b, "label": p.label}
        for p in pairs
        if by_locus[p.gene_b].interval.start < by_locus[p.gene_a].interval.start
    ]


def terminator_truth(terminators: Sequence[RegulonDBTerminator]) -> tuple[TerminatorRow, ...]:
    """``TerminatorSet`` as truth rows, sorted by position (spec §15.4 step 4 consumes these)."""
    rows = [
        TerminatorRow(
            contig=t.interval.contig,
            start=t.interval.start,
            end=t.interval.end,
            strand=t.interval.strand,
            source_id=t.terminator_id,
        )
        for t in terminators
    ]
    return tuple(sorted(rows, key=lambda r: (r.contig, r.start, r.end, r.strand, r.source_id)))


def summarise(
    pairs: Sequence[PairRow], mapping: GeneMapping, n_terminators: int
) -> dict[str, object]:
    """Counts and shares for PROGRESS.md and for M6-12 (spec example 5)."""
    counts = {
        label: sum(1 for p in pairs if p.label == label)
        for label in ("positive", "negative", "ambiguous")
    }
    total = len(pairs)
    labelled = counts["positive"] + counts["negative"]
    return {
        "pairs_total": total,
        "pairs_positive": counts["positive"],
        "pairs_negative": counts["negative"],
        "pairs_ambiguous": counts["ambiguous"],
        "pairs_labelled": labelled,
        "excluded_share": counts["ambiguous"] / total if total else 0.0,
        "genes_mapped": mapping.n_mapped,
        "genes_mapped_cds": mapping.n_mapped_cds,
        "genes_mapped_rna": len(mapping.rna),
        "genes_unmapped": len(mapping.unmapped),
        "mapped_share": mapping.mapped_share,
        "cds_mapped_share": mapping.cds_mapped_share,
        "mapped_share_ok": mapping.mapped_share >= MIN_MAPPED_SHARE,
        "terminators": n_terminators,
    }


# -------------------------------------------------------------- the writers


def pair_tsv(pairs: Sequence[PairRow]) -> str:
    """``pair_truth.tsv``: the four §15.4 columns, TU ids joined with ``;``."""
    lines = ["\t".join(PAIR_COLUMNS)]
    lines += ["\t".join((p.gene_a, p.gene_b, p.label, ";".join(p.tu_ids))) for p in pairs]
    return "\n".join(lines) + "\n"


def terminator_tsv(rows: Sequence[TerminatorRow]) -> str:
    """``terminator_truth.tsv``; coordinates stay 0-based half-open (spec §6.5, internal file)."""
    lines = ["\t".join(TERMINATOR_COLUMNS)]
    lines += ["\t".join((r.contig, str(r.start), str(r.end), r.strand, r.source_id)) for r in rows]
    return "\n".join(lines) + "\n"


def verified_release(regulondb_dir: Path) -> tuple[dict[str, str], dict[str, object]]:
    """Read the three downloads after checking each against ``provenance.json`` (ADR-0017).

    Returns the decoded tables by dataset and the provenance record. Raises
    :class:`ProvenanceError` (a ``ValueError``) naming the dataset when its size or SHA-256 differs
    from the record or the record lacks it, and naming ``provenance.json`` when the sidecar is
    missing or unreadable — verify before use (lead decision 1).
    """
    sidecar = regulondb_dir / "provenance.json"
    try:
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise ProvenanceError(
            f"{sidecar}: no provenance record; re-run validation/fetch_regulondb.py"
        ) from None
    except (OSError, json.JSONDecodeError) as exc:
        raise ProvenanceError(f"{sidecar}: unreadable provenance record: {exc}") from exc
    files = meta.get("files") if isinstance(meta, dict) else None
    if not isinstance(files, list):
        raise ProvenanceError(f"{sidecar}: provenance record has no 'files' list")
    records = {str(f.get("dataset")): f for f in files if isinstance(f, dict)}
    tables: dict[str, str] = {}
    for dataset in DATASETS:
        rec = records.get(dataset)
        if rec is None:
            raise ProvenanceError(f"{dataset}: not recorded in {sidecar}")
        path = regulondb_dir / str(rec.get("file") or f"{dataset}.txt")
        data = path.read_bytes()
        if len(data) != rec.get("bytes"):
            raise ProvenanceError(
                f"{dataset}: {path} has {len(data)} bytes, provenance records {rec.get('bytes')}"
            )
        digest = hashlib.sha256(data).hexdigest()
        if digest != rec.get("sha256"):
            raise ProvenanceError(
                f"{dataset}: {path} sha256 {digest} != recorded {rec.get('sha256')}"
            )
        tables[dataset] = data.decode("utf-8")
    return tables, meta


def build(
    regulondb_dir: Path,
    bakta_json: Path,
    genome: Path | None = None,
    contig: str | None = None,
    tolerance_bp: int = DEFAULT_TOLERANCE_BP,
) -> tuple[tuple[PairRow, ...], tuple[TerminatorRow, ...], dict[str, object]]:
    """Verify and read the downloads and the recorded Bakta JSON; return pairs, terminators and
    the summary (which carries the release and the three verified SHA-256s).

    ``genome`` supplies the nucleotides when the recorded JSON omits them (M2 lead decision 4).
    """
    tables, meta = verified_release(regulondb_dir)
    seqs = {r.id: r.sequence for r in read_genome_fasta(genome)} if genome else None
    model = parse_bakta_json(bakta_json, genome=seqs).model
    name = contig or chromosome_contig(model)
    rdb_genes = parse_gene_set(tables["GeneProductSet"], name)
    tus = parse_tu_set(tables["TUSet"])
    terms = parse_terminator_set(tables["TerminatorSet"], name)
    mapping = map_genes(rdb_genes, model, name, tolerance_bp)
    pairs = pair_truth(model, tus, rdb_genes, mapping, name)
    rows = terminator_truth(terms)
    summary = summarise(pairs, mapping, len(rows))
    summary["contig"] = name
    summary["tolerance_bp"] = tolerance_bp
    summary["wrap_pairs"] = wrap_pairs(model, pairs)
    summary["regulondb_version"] = meta.get("regulondb_version")
    summary["release_date"] = meta.get("release_date")
    summary["sha256"] = {
        dataset: hashlib.sha256(tables[dataset].encode("utf-8")).hexdigest() for dataset in DATASETS
    }
    return pairs, rows, summary


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--regulondb-dir",
        type=Path,
        required=True,
        help="a release directory written by fetch_regulondb.py",
    )
    ap.add_argument("--bakta-json", type=Path, required=True, help="the recorded K-12 Bakta JSON")
    ap.add_argument(
        "--genome", type=Path, default=None, help="genome FASTA, when the Bakta JSON omits nt"
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--contig", default=None, help="override the chromosome contig name")
    ap.add_argument("--tolerance-bp", type=int, default=DEFAULT_TOLERANCE_BP)
    args = ap.parse_args(argv)

    pairs, rows, summary = build(
        args.regulondb_dir, args.bakta_json, args.genome, args.contig, args.tolerance_bp
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "pair_truth.tsv").write_text(pair_tsv(pairs), encoding="utf-8")
    (args.out_dir / "terminator_truth.tsv").write_text(terminator_tsv(rows), encoding="utf-8")
    (args.out_dir / "truth_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"regulondb_truth: {json.dumps(summary)}")
    if not summary["mapped_share_ok"]:
        print(
            f"regulondb_truth: mapped share {summary['mapped_share']:.3f} < {MIN_MAPPED_SHARE}; "
            "escalate (spec §15.4 example 5), do not work around it",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
