"""Test-time K-12 excerpts and oracles for the MVP-04 giapsim tests.

Spec: docs/specs/MVP-04-ecoli-evidence-run.md.

Nothing built here is written into the repository (spec "Inputs", rule 3): every excerpt goes to a
pytest temporary directory. The excerpt is the first ``length`` bases of the committed
``tests/data/genomes/k12_mg1655.fna.gz`` (optionally rotated, so that a chosen position becomes
the new origin) with the recorded full-DB Bakta features that lie inside it, as a circular replicon.

The translation, ORF and TSV helpers are independent oracles (table 11 written out here), not
giapsim code.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import os
import shutil
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
K12_FASTA = REPO / "tests" / "data" / "genomes" / "k12_mg1655.fna.gz"
K12_DIR = REPO / "tests" / "data" / "recorded" / "bakta" / "1.12.2-full" / "k12"
K12_JSON = K12_DIR / "k12_mg1655.json.gz"
K12_CONTIG_MAP = K12_DIR / "contig_map.tsv"
CONTIG = "U00096.3"
K12_LENGTH = 4_641_652

#: The spec's input table (MVP-04 "Inputs").
SHA_FASTA_FILE = "4bfd7163f6f0de6394a2131cbfdbe75100daa55cbee227ec3e403a44640ac03c"
SHA_SEQUENCE_UPPER = "e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0"
SHA_JSON_FILE = "c08c6466edf99d16cbe79679a8971b6e06e081cfc6adb9cc6f1ea9d51d3f083d"

#: RegulonDB 14.5.0 as fetched by M6-11 (never committed, R11); tests skip when it is absent.
REGULONDB_DIR = (
    Path(os.environ.get("GIAP_DATA_DIR", REPO / "local" / "data")) / "regulondb" / "14.5.0"
)
TUSET = REGULONDB_DIR / "TUSet.txt"

STARTS = frozenset({"ATG", "GTG", "TTG"})
_BASES = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE11 = {
    a + b + c: _AA11[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
STOPS = frozenset(c for c, a in CODE11.items() if a == "*")
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
NOVEL_CLASSES = ("intergenic", "cds_outframe", "antisense")


# ---------------------------------------------------------------- sequence oracles


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def codons(nt: str) -> list[str]:
    return [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]


def translate(nt: str) -> str:
    """Plain table-11 translation (stops as ``*``, first codon not forced to M)."""
    return "".join(CODE11.get(c, "X") for c in codons(nt.upper()))


def il(s: str) -> str:
    return s.replace("I", "L")


def segment(seq: str, start1: int, end1: int, strand: str) -> str:
    """Bases of the 1-based closed ``[start1, end1]`` on ``strand`` (no origin wrap)."""
    nt = seq[start1 - 1 : end1]
    return revcomp(nt) if strand == "-" else nt


def read_fasta(path: Path) -> dict[str, str]:
    """``{first header word: upper-case sequence}`` of a plain or gzipped FASTA."""
    raw = path.read_bytes()
    text = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode()
    out: dict[str, list[str]] = {}
    cur = ""
    for line in text.splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif line.strip():
            out[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def header_of(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as fh:
        return fh.readline().rstrip("\n").split("\t")


def values(obj: Any) -> Iterator[Any]:
    """Every scalar inside a JSON document (for "the hash is recorded somewhere" checks)."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield k
            yield from values(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from values(v)
    else:
        yield obj


# ---------------------------------------------------------------- K-12 inputs


def k12_text() -> tuple[str, str]:
    """(header line without ``>``, upper-case sequence) of the committed K-12 FASTA."""
    text = gzip.decompress(K12_FASTA.read_bytes()).decode()
    lines = text.splitlines()
    return lines[0][1:], "".join(lines[1:]).upper()


def k12_json() -> dict[str, Any]:
    data: dict[str, Any] = json.loads(gzip.decompress(K12_JSON.read_bytes()))
    return data


def write_fasta_gz(path: Path, header: str, seq: str, width: int = 80) -> Path:
    body = "".join(seq[i : i + width] + "\n" for i in range(0, len(seq), width))
    path.write_bytes(gzip.compress(f">{header}\n{body}".encode(), mtime=0))
    return path


def write_json_gz(path: Path, data: dict[str, Any]) -> Path:
    path.write_bytes(gzip.compress(json.dumps(data).encode(), mtime=0))
    return path


@dataclass(frozen=True)
class Gene:
    """A recorded protein gene: 0-based half-open ``[start, end)``."""

    locus: str
    kind: str
    start: int
    end: int
    strand: str
    aa: str
    gene: str | None
    truncated: bool
    pseudo: bool

    @property
    def eligible(self) -> bool:
        """Counts in the eligible denominator (D8: truncated genes are ``partial_excluded``)."""
        return not self.truncated and not self.pseudo


def protein_genes(data: dict[str, Any]) -> list[Gene]:
    """CDS and sORF features of a Bakta JSON, in coordinate order."""
    out = []
    for f in data["features"]:
        if f.get("type") not in ("cds", "sorf"):
            continue
        out.append(
            Gene(
                str(f["locus"]),
                str(f["type"]),
                int(f["start"]) - 1,
                int(f["stop"]),
                str(f["strand"]),
                str(f.get("aa") or "").rstrip("*"),
                f.get("gene"),
                bool(f.get("truncated")),
                bool(f.get("pseudo") or f.get("pseudogene")),
            )
        )
    return sorted(out, key=lambda g: (g.start, g.end, g.strand, g.locus))


@dataclass(frozen=True)
class Excerpt:
    fasta: Path
    bakta_json: Path
    seq: str
    data: dict[str, Any]

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def genes(self) -> list[Gene]:
        return protein_genes(self.data)


def make_excerpt(outdir: Path, length: int = 200_000, rotate: int = 0) -> Excerpt:
    """The first ``length`` K-12 bases as a circular replicon (sequence rotated left by
    ``rotate``: base ``rotate`` becomes position 1) with the recorded features fully inside the
    window; features spanning ``rotate`` are dropped. Files keep the committed names' suffixes."""
    outdir.mkdir(parents=True, exist_ok=True)
    _, full = k12_text()
    win = full[:length]
    seq = win[rotate:] + win[:rotate]
    data = k12_json()
    feats = []
    for f in data["features"]:
        s1, e1 = int(f["start"]), int(f["stop"])
        if not (1 <= s1 <= e1 <= length):
            continue
        if rotate and s1 - 1 < rotate < e1:
            continue
        g = dict(f)
        g["start"] = (s1 - 1 - rotate) % length + 1
        g["stop"] = (e1 - 1 - rotate) % length + 1
        feats.append(g)
    feats.sort(key=lambda f: (int(f["start"]), int(f["stop"])))
    data["features"] = feats
    data["sequences"] = [dict(data["sequences"][0], length=length, topology="circular")]
    fasta = write_fasta_gz(
        outdir / "k12_excerpt.fna.gz", f"{CONTIG} K-12 excerpt (test time) [topology=circular]", seq
    )
    js = write_json_gz(outdir / "k12_excerpt.json.gz", data)
    shutil.copyfile(K12_CONTIG_MAP, outdir / "contig_map.tsv")
    return Excerpt(fasta, js, seq, data)


# ---------------------------------------------------------------- ORF oracle


def plus_orfs(seq: str) -> list[tuple[int, int, int]]:
    """Every plus-strand ORF of a linear ``seq``: (first start, last start, stop) 0-based codon
    starts, one per stop codon that has at least one in-frame start since the previous stop."""
    out = []
    for fr in range(3):
        first = last = None
        for j, c in enumerate(codons(seq[fr:])):
            p = fr + 3 * j
            if c in STOPS:
                if first is not None and last is not None:
                    out.append((first, last, p))
                first = last = None
            elif c in STARTS:
                if first is None:
                    first = p
                last = p
    return sorted(out)


def free_intergenic_orf(
    seq: str, data: dict[str, Any], lo_aa: int = 40, hi_aa: int = 140
) -> tuple[int, int, int]:
    """A plus-strand ORF (first start, last start, stop) whose longest form has ``lo_aa``-``hi_aa``
    residues and that overlaps no recorded feature of any type (with a 10-bp margin)."""
    spans = sorted((int(f["start"]) - 1, int(f["stop"])) for f in data["features"])
    for first, last, stop in plus_orfs(seq):
        n_aa = (stop - first) // 3
        if not lo_aa <= n_aa <= hi_aa or last + 3 >= stop:
            continue
        a, b = first - 10, stop + 13
        if any(s < b and a < e for s, e in spans):
            continue
        return first, last, stop
    raise AssertionError("no feature-free intergenic ORF in the excerpt")


# ---------------------------------------------------------------- RegulonDB-shaped fixtures

_TU_HEADER = (
    "1)id\t2)name\t3)operonId\t4)operonName\t5)tuGenes\t6)promoterId\t7)promoterName\t"
    "8)terminatorIds\t9)tuEvidence\t10)addEvidence\t11)confidenceLevel\t12)PMIDs "
)
_GENE_HEADER = (
    "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand\t7)productId\t"
    "8)productName\t9)geneEvidences\t10)genePMIDS\t11)confidenceLevel\t12)relatedBnumbers\t"
    "13)otherDbsIds "
)


def write_regulondb_fixture(
    outdir: Path,
    genes: list[tuple[str, int, int, str]],
    tus: list[tuple[str, list[str]]],
) -> Path:
    """Hand-made RegulonDB-shaped ``TUSet.txt`` + ``GeneProductSet.txt`` (never RegulonDB bytes).

    ``genes``: (name, leftEndPos, rightEndPos, ``+``/``-``), 1-based closed; ``tus``: (TU id,
    member gene names). Returns the ``TUSet.txt`` path (``GeneProductSet.txt`` sits beside it)."""
    outdir.mkdir(parents=True, exist_ok=True)
    gl = ["# hand-made RegulonDB-shaped fixture (MVP-04 tests)", _GENE_HEADER]
    for i, (name, left, right, strand) in enumerate(genes):
        sense = "forward" if strand == "+" else "reverse"
        gl.append(f"FIXGN{i:05d}\t{name}\tb9{i:03d}\t{left}\t{right}\t{sense}\tFIXPD{i:05d}"
                  "\tfixture product\t\t\tNone\t\t ")  # fmt: skip
    (outdir / "GeneProductSet.txt").write_text("\n".join(gl) + "\n")
    tl = ["# hand-made RegulonDB-shaped fixture (MVP-04 tests)", _TU_HEADER]
    for tid, members in tus:
        tl.append(f"{tid}\t{'-'.join(members)}\tOP{tid}\t{'-'.join(members)}\t{';'.join(members)}"
                  "\t\t\t\tCOMP-AINF:W\t\tW\t ")  # fmt: skip
    path = outdir / "TUSet.txt"
    path.write_text("\n".join(tl) + "\n")
    return path


# ---------------------------------------------------------------- logging capture


@contextmanager
def capture_warnings() -> Iterator[list[logging.LogRecord]]:
    """Every WARNING-or-higher record emitted anywhere while the block runs."""
    records: list[logging.LogRecord] = []

    class _H(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(record)

    h = _H(level=logging.WARNING)
    root = logging.getLogger()
    old = root.level
    root.addHandler(h)
    root.setLevel(min(old, logging.WARNING) if old else logging.WARNING)
    try:
        yield records
    finally:
        root.removeHandler(h)
        root.setLevel(old)


# ---------------------------------------------------------------- output readers


def tu_of(outdir: Path) -> dict[str, str]:
    return {r["locus_tag"]: r["tu_id"] for r in read_tsv(outdir / "truth" / "genes.tsv")}


def abundance_of(outdir: Path) -> dict[str, float]:
    return {r["locus_tag"]: float(r["log_abundance"]) for r in read_tsv(outdir / "abundance.tsv")}


def tu_number(tu_id: str) -> int:
    digits = "".join(ch for ch in tu_id.rsplit("_", 1)[-1] if ch.isdigit())
    return int(digits)


def island_genes(outdir: Path) -> dict[str, set[str]]:
    """``island_id -> loci`` from ``truth/islands.tsv`` (``first_tu_id``..``last_tu_id``, the
    ``TU_<n>`` ids of ``truth/genes.tsv``)."""
    tus = tu_of(outdir)
    out = {}
    for r in read_tsv(outdir / "truth" / "islands.tsv"):
        lo, hi = tu_number(r["first_tu_id"]), tu_number(r["last_tu_id"])
        out[r["island_id"]] = {g for g, t in tus.items() if lo <= tu_number(t) <= hi}
    return out


def true_gene_loci(outdir: Path) -> set[str]:
    """Loci with at least one correct (not wrong, not decoy) peptide in ``truth/peptides.tsv``."""
    return {
        r["locus_tag"]
        for r in read_tsv(outdir / "truth" / "peptides.tsv")
        if r["is_wrong"] == "false" and r["origin"] != "decoy"
    }


# ---------------------------------------------------------------- novel-ORF rule oracles


def _overlaps(a0: int, a1: int, b0: int, b1: int) -> bool:
    return a0 < b1 and b0 < a1


def _in_frame(s0: int, e0: int, strand: str, g: Gene) -> bool:
    """Same strand and the same reading frame as ``g`` (would share its peptides)."""
    if g.strand != strand:
        return False
    return (s0 - g.start) % 3 == 0 if strand == "+" else (g.end - e0) % 3 == 0


def novel_rule_violations(
    rows: list[dict[str, str]],
    k12: str,
    sample: str,
    genes: list[Gene],
    len_aa: tuple[int, int] = (30, 150),
) -> list[str]:
    """Violations of spec pipeline step 4 by the ``intergenic`` / ``cds_outframe`` /
    ``antisense`` rows of ``truth/novel_orfs.tsv``: each is an existing ORF of the unchanged K-12
    sequence (ATG/GTG/TTG ... first in-frame stop), ``len_aa`` residues, no origin crossing (Q7),
    not in frame with an overlapping same-strand protein gene, protein not part of an annotated
    protein (I = L), and its class matches its position."""
    hay = "|".join(il(g.aa) for g in genes)
    length = len(k12)
    bad = []
    for r in rows:
        cls = r["class"]
        if cls not in NOVEL_CLASSES:
            continue
        s1, e1, strand, aa = int(r["start"]), int(r["end"]), r["strand"], r["aa"]
        tag = f"{r['orf_id']} {cls} {s1}-{e1}{strand}"
        if not (1 <= s1 < e1 <= length):
            bad.append(f"{tag}: crosses the origin or lies outside the replicon")
            continue
        nt = segment(k12, s1, e1, strand)
        cods = codons(nt)
        if len(nt) % 3 or cods[0] not in STARTS or cods[-1] not in STOPS:
            bad.append(f"{tag}: not start..stop on K-12")
            continue
        if STOPS.intersection(cods[1:-1]):
            bad.append(f"{tag}: internal stop (not the first in-frame stop)")
        n_aa = len(cods) - 1
        if not len_aa[0] <= n_aa <= len_aa[1]:
            bad.append(f"{tag}: {n_aa} aa outside {len_aa}")
        expect = {"M" + translate(nt[3:-3]), "M" + translate(segment(sample, s1, e1, strand)[3:-3])}
        if aa not in expect:
            bad.append(f"{tag}: aa is not the translation of its span")
        if len(aa) > 1 and il(aa[1:]) in hay:
            bad.append(f"{tag}: protein is part of an annotated protein (I = L)")
        s0, e0 = s1 - 1, e1
        near = [g for g in genes if _overlaps(s0, e0, g.start, g.end)]
        if any(_in_frame(s0, e0, strand, g) for g in near):
            bad.append(f"{tag}: in frame with an overlapping annotated protein gene")
        if cls == "intergenic" and near:
            bad.append(f"{tag}: intergenic ORF overlaps {near[0].locus}")
        if cls == "cds_outframe" and not any(g.strand == strand for g in near):
            bad.append(f"{tag}: cds_outframe ORF overlaps no same-strand CDS")
        if cls == "antisense" and not any(g.strand != strand for g in near):
            bad.append(f"{tag}: antisense ORF overlaps no opposite-strand CDS")
    return bad


def nterm_violations(rows: list[dict[str, str]], k12: str, genes: list[Gene]) -> list[str]:
    """Violations of the ``nterm_extension`` rule (spec step 4, example 7): the row runs from an
    in-frame start codon >= 5 codons upstream of an annotated CDS start to that CDS's stop, with no
    in-frame stop in the extension; on the minus strand upstream means higher coordinates."""
    plus_end = {(g.end, "+"): g for g in genes}
    minus_start = {(g.start, "-"): g for g in genes}
    bad = []
    for r in rows:
        if r["class"] != "nterm_extension":
            continue
        s1, e1, strand = int(r["start"]), int(r["end"]), r["strand"]
        tag = f"{r['orf_id']} {s1}-{e1}{strand}"
        g = plus_end.get((e1, "+")) if strand == "+" else minus_start.get((s1 - 1, "-"))
        if g is None or g.kind != "cds":
            bad.append(f"{tag}: shares its stop with no annotated CDS")
            continue
        ext = g.start - (s1 - 1) if strand == "+" else e1 - g.end
        if ext < 15 or ext % 3:
            bad.append(f"{tag}: extension of {ext} nt (need >= 15, a multiple of 3)")
            continue
        nt = segment(k12, s1, e1, strand)
        if nt[:3] not in STARTS:
            bad.append(f"{tag}: true start {nt[:3]} is not ATG/GTG/TTG")
        if STOPS.intersection(codons(nt[:ext])):
            bad.append(f"{tag}: in-frame stop between the true and the annotated start")
    return bad
