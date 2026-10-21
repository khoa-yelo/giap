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
import json
import sys
from bisect import bisect_left, bisect_right
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from giap.annotate.bakta import parse_bakta_json
from giap.io.fasta import read_genome_fasta
from giap.model.coords import Interval, Strand
from giap.model.genes import Gene, GeneModel

#: RegulonDB gene coordinates agree with Bakta exactly on K-12; a few bp of slack absorbs
#: annotation differences at the start codon without letting a neighbour win.
DEFAULT_TOLERANCE_BP = 30

#: Below this mapped share the result is escalated, not worked around (spec example 5).
MIN_MAPPED_SHARE = 0.8

PAIR_COLUMNS = ("gene_a", "gene_b", "label", "tu_ids")
TERMINATOR_COLUMNS = ("contig", "start", "end", "strand", "source_id")

_STRANDS: Mapping[str, Strand] = {"forward": "+", "reverse": "-", "+": "+", "-": "-"}


class RegulonDBFormatError(ValueError):
    """A RegulonDB download does not have the expected columns."""


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
    """RegulonDB gene names resolved onto Bakta locus tags, with the share that resolved."""

    locus_by_name: Mapping[str, str]
    how: Mapping[str, str]  # gene name -> "coordinates" | "name"
    unmapped: tuple[str, ...]  # RegulonDB gene names with no Bakta counterpart

    @property
    def n_mapped(self) -> int:
        return len(self.locus_by_name)

    @property
    def mapped_share(self) -> float:
        total = self.n_mapped + len(self.unmapped)
        return self.n_mapped / total if total else 0.0


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


def map_genes(
    rdb_genes: Sequence[RegulonDBGene],
    model: GeneModel,
    contig: str,
    tolerance_bp: int = DEFAULT_TOLERANCE_BP,
) -> GeneMapping:
    """Map RegulonDB gene names onto Bakta locus tags, by coordinates first and then by name.

    A coordinate match needs the same strand and both ends within ``tolerance_bp``; among the
    candidates the smallest total offset wins, globally and deterministically, and each Bakta gene
    is claimed at most once. Names are matched only when the name is unique on both sides.
    """
    bakta = [g for g in model.genes if g.interval.contig == contig]
    candidates: list[tuple[int, str, str]] = []
    by_start = sorted(bakta, key=lambda g: g.interval.start)
    starts = [g.interval.start for g in by_start]
    for rg in rdb_genes:
        lo, hi = rg.interval.start - tolerance_bp, rg.interval.start + tolerance_bp
        for g in _window(by_start, starts, lo, hi):
            if g.interval.strand != rg.interval.strand:
                continue
            if abs(g.interval.end - rg.interval.end) > tolerance_bp:
                continue
            offset = abs(g.interval.start - rg.interval.start) + abs(
                g.interval.end - rg.interval.end
            )
            candidates.append((offset, rg.name, g.locus_tag))

    locus_by_name: dict[str, str] = {}
    how: dict[str, str] = {}
    claimed: set[str] = set()
    for _offset, name, locus in sorted(candidates):
        if name in locus_by_name or locus in claimed:
            continue
        locus_by_name[name] = locus
        how[name] = "coordinates"
        claimed.add(locus)

    seen = Counter(g.name for g in rdb_genes if g.name)
    unresolved = sorted(n for n in seen if n not in locus_by_name)
    for name, locus in _unique_name_matches(unresolved, seen, bakta, claimed).items():
        locus_by_name[name] = locus
        how[name] = "name"
        claimed.add(locus)

    unmapped = tuple(n for n in unresolved if n not in locus_by_name)
    return GeneMapping(locus_by_name=locus_by_name, how=how, unmapped=unmapped)


def _window(genes: Sequence[Gene], starts: Sequence[int], lo: int, hi: int) -> Sequence[Gene]:
    """Genes whose start lies in ``[lo, hi]`` (``starts`` is ``genes``' sorted start list)."""
    return genes[bisect_left(starts, lo) : bisect_right(starts, hi)]


def _unique_name_matches(
    names: Sequence[str],
    occurrences: Mapping[str, int],
    bakta: Sequence[Gene],
    claimed: Iterable[str],
) -> dict[str, str]:
    """Name matches that are unambiguous on both sides and do not steal a claimed locus tag."""
    taken = set(claimed)
    by_name: dict[str, list[str]] = {}
    for g in bakta:
        if g.gene and g.locus_tag not in taken:
            by_name.setdefault(g.gene.casefold(), []).append(g.locus_tag)
    out: dict[str, str] = {}
    used: set[str] = set()
    for name in names:
        if occurrences.get(name, 0) != 1:
            continue
        hits = by_name.get(name.casefold(), ())
        if len(hits) == 1 and hits[0] not in used:
            out[name] = hits[0]
            used.add(hits[0])
    return out


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
    """Consecutive genes on ``contig`` (coordinate order) whose strands agree."""
    genes = sorted(
        (g for g in model.genes if g.interval.contig == contig),
        key=lambda g: (g.interval.start, g.interval.end, g.locus_tag),
    )
    return tuple((a, b) for a, b in pairwise(genes) if a.interval.strand == b.interval.strand)


def tus_by_locus(
    tus: Sequence[TranscriptionUnit], mapping: GeneMapping
) -> dict[str, frozenset[str]]:
    """Locus tag -> the ids of the transcription units that contain it."""
    out: dict[str, set[str]] = {}
    for tu in tus:
        for name in tu.genes:
            locus = mapping.locus_by_name.get(name)
            if locus is not None:
                out.setdefault(locus, set()).add(tu.tu_id)
    return {locus: frozenset(ids) for locus, ids in out.items()}


def _unmapped_between(a: Gene, b: Gene, unmapped_intervals: Sequence[tuple[int, int]]) -> bool:
    """True when a RegulonDB gene without a Bakta counterpart lies between ``a`` and ``b``."""
    lo, hi = a.interval.end, b.interval.start
    if hi <= lo:
        return False
    return any(start < hi and end > lo for start, end in unmapped_intervals)


def label_pair(
    a: Gene,
    b: Gene,
    tu_ids: Mapping[str, frozenset[str]],
    mapped_loci: frozenset[str],
    blocked: bool,
) -> PairRow:
    """Label one adjacent same-strand pair exactly as §15.4 step 2 and its degenerate table say.

    In order: a pair is **ambiguous** when either gene has no RegulonDB counterpart or when a
    RegulonDB-only gene lies between them (``blocked`` — the degenerate table states this without
    qualification, so it wins over the TU evidence); otherwise **positive** when some TU holds
    both (two overlapping TUs still count); otherwise **negative** when each gene is in at least
    one TU and none holds both; otherwise **ambiguous** (a gene in no TU at all).
    """
    a_tus = tu_ids.get(a.locus_tag, frozenset())
    b_tus = tu_ids.get(b.locus_tag, frozenset())
    union = tuple(sorted(a_tus | b_tus))
    unmapped = a.locus_tag not in mapped_loci or b.locus_tag not in mapped_loci
    if unmapped or blocked:
        return PairRow(a.locus_tag, b.locus_tag, "ambiguous", union)
    both = a_tus & b_tus
    if both:
        return PairRow(a.locus_tag, b.locus_tag, "positive", tuple(sorted(both)))
    if a_tus and b_tus:
        return PairRow(a.locus_tag, b.locus_tag, "negative", union)
    return PairRow(a.locus_tag, b.locus_tag, "ambiguous", union)


def pair_truth(
    model: GeneModel,
    tus: Sequence[TranscriptionUnit],
    rdb_genes: Sequence[RegulonDBGene],
    mapping: GeneMapping,
    contig: str,
) -> tuple[PairRow, ...]:
    """All same-strand adjacent pairs of ``contig`` with their §15.4 labels."""
    by_locus = tus_by_locus(tus, mapping)
    mapped_loci = frozenset(mapping.locus_by_name.values())
    unmapped_names = frozenset(mapping.unmapped)
    unmapped_intervals = sorted(
        (g.interval.start, g.interval.end) for g in rdb_genes if g.name in unmapped_names
    )
    return tuple(
        label_pair(a, b, by_locus, mapped_loci, _unmapped_between(a, b, unmapped_intervals))
        for a, b in adjacent_same_strand_pairs(model, contig)
    )


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
        "genes_unmapped": len(mapping.unmapped),
        "mapped_share": mapping.mapped_share,
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


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def build(
    regulondb_dir: Path,
    bakta_json: Path,
    genome: Path | None = None,
    contig: str | None = None,
    tolerance_bp: int = DEFAULT_TOLERANCE_BP,
) -> tuple[tuple[PairRow, ...], tuple[TerminatorRow, ...], dict[str, object]]:
    """Read the downloads and the recorded Bakta JSON; return pairs, terminators and the summary.

    ``genome`` supplies the nucleotides when the recorded JSON omits them (M2 lead decision 4).
    """
    seqs = {r.id: r.sequence for r in read_genome_fasta(genome)} if genome else None
    model = parse_bakta_json(bakta_json, genome=seqs).model
    name = contig or chromosome_contig(model)
    rdb_genes = parse_gene_set(_read(regulondb_dir / "GeneProductSet.txt"), name)
    tus = parse_tu_set(_read(regulondb_dir / "TUSet.txt"))
    terms = parse_terminator_set(_read(regulondb_dir / "TerminatorSet.txt"), name)
    mapping = map_genes(rdb_genes, model, name, tolerance_bp)
    pairs = pair_truth(model, tus, rdb_genes, mapping, name)
    rows = terminator_truth(terms)
    summary = summarise(pairs, mapping, len(rows))
    summary["contig"] = name
    summary["tolerance_bp"] = tolerance_bp
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
