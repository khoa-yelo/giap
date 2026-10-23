"""giapsim on a real genome: the MVP-04 K-12 evidence inputs (MVP-04-ecoli-evidence-run.md).

``simulate_on_genome`` reads a real genome FASTA and its recorded Bakta JSON (never modified) and
writes the simulated inputs of a GIAP evidence run plus their truth:

1. **Read** (rule 1): the input hashes are checked first (``verify_inputs``); the CDS and sORF
   features become ``Feat`` (truncated CDSs and pseudogenes are not expressed), contig IDs go
   through ``contig_map.tsv`` next to the JSON, the topology comes from the FASTA header tag.
2. **Sample strain** (step 2, Q4): ``strain.derive_reference`` is applied to the genome and its
   output is the *sample*; GIAP gets the official genome. SNVs only, so the coordinates agree.
   ``truth/alt_sites.tsv`` has the M7-12 columns with ``ref_aa`` the genome's residue and
   ``alt_aa`` the sample's (swapped roles); SNVs in a gene's initiator codon are left out (the
   initiator reads M on both sides in GIAP).
3. **TUs and abundance** (step 3, Q3): connected components of same-strand genes sharing a
   RegulonDB TU (match: same strand and same stop, the right end on ``+``, the left end on ``-``;
   gene coordinates from ``GeneProductSet.txt`` next to ``TUSet.txt``), else the gap rule
   (adjacent in genome order, same strand, gap <= 40 bp) with a WARNING.
   ``layout.layout_from_simulation`` numbers them ``TU_<n>`` (the replicon is handed to it as
   linear, so a RegulonDB TU with an interleaved gene keeps its min..max bounds; no TU crosses
   the origin) and ``expression.expression_for_simulation`` gives the TU-level AR(1) plus the
   planted islands; the natural-log abundance of the expressed genes is centred to
   ``LOG_ABUNDANCE[0]``.
4. **Novel ORFs** (step 4, Q5, Q7): ``intergenic`` / ``cds_outframe`` / ``antisense`` ORFs are
   *found* on the unchanged genome (ATG/GTG/TTG to the first in-frame stop, ``novel_len_aa``,
   one per stop, no origin crossing, in frame with no overlapping same-strand protein gene, no
   overlap with an RNA gene or a same-strand pseudogene, protein not part of an annotated one
   with I = L, start..stop intact on the sample too); ``nterm_extension`` CDSs get an in-frame
   upstream start >= 5 codons before the annotated one with no stop between (``Feat.true_start``;
   upstream = higher coordinates on ``-``). A class short of supply takes what exists, with a
   WARNING. Novel ORFs (not extensions) take abundances from the upper tercile of the genes'.
5. **Peptides** (step 5): ``proteomics.simulate_proteomics`` on the sample with the abundance
   above (``log_abundance`` hook). The scored null decoys (``decoys = "scored_null"``, step 6)
   are part 2; this module runs the M3-19 decoys.
6. **Outputs** (step 8): ``peptides.tsv``, ``abundance.tsv``, ``config.yaml`` (``seed`` and
   ``sample_name`` only), ``truth/{peptides,genes,islands,alt_sites,novel_orfs}.tsv`` (+
   ``truth/novel_orfs.json``, ``truth/proteomics.json``), ``truth/sample_genome.fna`` (debugging
   only, never a GIAP input) and ``provenance.json``. Everything is staged in a temporary
   directory and moved into ``outdir`` only when it is complete; a ``peptides.tsv`` without
   target rows is refused (D9).

Draws use ``numpy.random.default_rng([seed, 40, k])`` (k = 0 TU expression, 1 novel ORFs, 2 novel
abundances) next to the existing giapsim streams. Nothing is imported from ``giap``.

CLI: ``python -m giapsim.real_genome --genome G --annotation J --seed S --outdir D
[--tu-set TUSet.txt]`` runs ``scenarios.k12_evidence()`` with the hash check on.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import math
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Literal

import numpy as np

from giapsim.codons import CODE_11, STOPS, codons_of, revcomp
from giapsim.expression import ISLAND_COLUMNS, ExpressionConfig, expression_for_simulation
from giapsim.features import STARTS, Feat
from giapsim.layout import layout_from_simulation
from giapsim.proteomics import (
    LOG_ABUNDANCE,
    ProteomicsConfig,
    _codon_mask,
    _cover,
    simulate_proteomics,
)
from giapsim.scenario_truth import (
    ALT_COLUMNS,
    NOVEL_COLUMNS,
    alt_site_rows,
    novel_report,
    novel_rows,
)
from giapsim.strain import StrainConfig, derive_reference

__all__ = [
    "EXPECTED_HASHES",
    "GAP_RULE_BP",
    "NOVEL_FOUND_CLASSES",
    "RealGenomeConfig",
    "RealGenomeResult",
    "main",
    "simulate_on_genome",
    "verify_inputs",
]

_LOG = logging.getLogger(__name__)

#: The MVP-04 input table: (file sha256, upper-case sequence sha256) of the committed K-12 FASTA
#: and the file sha256 of the recorded full-DB Bakta JSON.
SHA_FASTA_FILE = "4bfd7163f6f0de6394a2131cbfdbe75100daa55cbee227ec3e403a44640ac03c"
SHA_SEQUENCE_UPPER = "e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0"
SHA_JSON_FILE = "c08c6466edf99d16cbe79679a8971b6e06e081cfc6adb9cc6f1ea9d51d3f083d"
EXPECTED_HASHES: Mapping[str, str] = {
    "genome_file": SHA_FASTA_FILE,
    "genome_sequence": SHA_SEQUENCE_UPPER,
    "annotation_file": SHA_JSON_FILE,
}
GAP_RULE_BP = 40
NOVEL_FOUND_CLASSES = ("intergenic", "cds_outframe", "antisense")
NTERM_MIN_CODONS = 5
SAMPLE_NAME = "k12_mvp04"
_STREAM = 40  # numpy SeedSequence key of this module (sub-keys 0, 1, 2)
_RNA_TYPES = frozenset({"trna", "rrna", "tmrna", "ncrna"})


@dataclass(frozen=True)
class RealGenomeConfig:
    """Knobs of the real-genome simulation (spec "Interface"); ``scenarios.k12_evidence``."""

    proteomics: ProteomicsConfig = field(default_factory=ProteomicsConfig)
    strain: StrainConfig = field(
        default_factory=lambda: StrainConfig(snv_rate_per_kb=0.5, indel_rate_per_kb=0.0)
    )
    expression: ExpressionConfig = field(default_factory=ExpressionConfig)
    tu_source: Literal["regulondb", "gap_rule"] = "regulondb"
    n_novel_per_class: int = 300
    novel_len_aa: tuple[int, int] = (30, 150)
    sample_name: str = SAMPLE_NAME

    def __post_init__(self) -> None:
        object.__setattr__(self, "novel_len_aa", tuple(int(x) for x in self.novel_len_aa))
        if self.strain.indel_rate_per_kb != 0.0:
            raise ValueError("real-genome sample strains are SNV-only (indel_rate_per_kb = 0)")
        if self.tu_source not in ("regulondb", "gap_rule"):
            raise ValueError(f"unknown tu_source {self.tu_source!r}")
        lo, hi = self.novel_len_aa
        if not 1 <= lo <= hi:
            raise ValueError("novel_len_aa must satisfy 1 <= lo <= hi")
        if self.n_novel_per_class < 0:
            raise ValueError("n_novel_per_class must be >= 0")


@dataclass(frozen=True)
class RealGenomeResult:
    """Where the outputs are and what ``provenance.json`` records."""

    outdir: Path
    paths: Mapping[str, Path]
    provenance: Mapping[str, Any]


# ---------------------------------------------------------------- inputs


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _maybe_gunzip(raw: bytes) -> bytes:
    return gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw


def _read_fasta(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """(``{first header word: upper-case sequence}``, ``{id: topology}``) of a plain or gzipped
    FASTA; a header carrying ``[topology=circular]`` is circular, any other linear."""
    text = _maybe_gunzip(Path(path).read_bytes()).decode()
    seqs: dict[str, list[str]] = {}
    topo: dict[str, str] = {}
    cur = ""
    for line in text.splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
            topo[cur] = "circular" if "[topology=circular]" in line else "linear"
        elif line.strip():
            seqs[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}, topo


def _sequence_sha(seqs: Mapping[str, str]) -> str:
    return _sha256_bytes("".join(seqs.values()).encode())


def verify_inputs(genome_fasta: Path, bakta_json: Path) -> dict[str, str]:
    """Check the genome file, its upper-case sequence and the annotation JSON against
    ``EXPECTED_HASHES`` (spec rule 1); ``ValueError("input_hash_mismatch:<file name>")`` on the
    first difference. Returns the hashes."""
    genome_fasta, bakta_json = Path(genome_fasta), Path(bakta_json)
    got = {"genome_file": _sha256(genome_fasta)}
    if got["genome_file"] != EXPECTED_HASHES["genome_file"]:
        raise ValueError(f"input_hash_mismatch:{genome_fasta.name}")
    got["genome_sequence"] = _sequence_sha(_read_fasta(genome_fasta)[0])
    if got["genome_sequence"] != EXPECTED_HASHES["genome_sequence"]:
        raise ValueError(f"input_hash_mismatch:{genome_fasta.name}")
    got["annotation_file"] = _sha256(bakta_json)
    if got["annotation_file"] != EXPECTED_HASHES["annotation_file"]:
        raise ValueError(f"input_hash_mismatch:{bakta_json.name}")
    return got


def _contig_map(bakta_json: Path) -> dict[str, str]:
    path = Path(bakta_json).parent / "contig_map.tsv"
    if not path.is_file():
        return {}
    with path.open(newline="", encoding="utf-8") as fh:
        return {r["original_id"]: r["contig"] for r in csv.DictReader(fh, delimiter="\t")}


def _feats(data: Mapping[str, Any], seqs: Mapping[str, str]) -> tuple[list[Feat], list[Feat]]:
    """(protein genes: CDS and sORF, RNA genes) of a Bakta JSON, 0-based half-open, in genome
    order; ``nt`` is read from the genome. Truncated CDSs and pseudogenes get a
    ``pseudo_reason`` (never expressed; D8)."""
    prot, rna = [], []
    for f in data["features"]:
        typ = str(f.get("type", ""))
        cid = str(f.get("contig", f.get("sequence")))
        s0, e0 = int(f["start"]) - 1, int(f["stop"])
        if e0 <= s0:
            e0 += len(seqs[cid])
        if typ.lower() in _RNA_TYPES:
            rna.append(Feat(typ, cid, s0, e0, str(f.get("strand", "+")), locus=str(f.get("locus"))))
            continue
        if typ not in ("cds", "sorf"):
            continue
        strand = str(f["strand"])
        seq = seqs[cid]
        nt = seq[s0:e0] if e0 <= len(seq) else seq[s0:] + seq[: e0 - len(seq)]
        nt = nt if strand == "+" else revcomp(nt)
        pseudo = bool(f.get("pseudo") or f.get("pseudogene"))
        reason = "truncated" if f.get("truncated") else "pseudo" if pseudo else ""
        extra: dict[str, Any] = {"gene": f.get("gene"), "product": f.get("product"),
                                 "start_type": f.get("start_type") or nt[:3]}  # fmt: skip
        if pseudo:
            extra["pseudo"] = True
        prot.append(Feat(typ, cid, s0, e0, strand, nt=nt, aa=str(f.get("aa") or "").rstrip("*"),
                         extra=extra, pseudo_reason=reason, locus=str(f["locus"])))  # fmt: skip
    order = {c: i for i, c in enumerate(seqs)}
    prot.sort(key=lambda g: (order[g.contig], g.start, g.end, g.strand, g.locus))
    return prot, rna


# ---------------------------------------------------------------- TUs


def gap_rule_tus(genes: Sequence[Feat], gap_bp: int = GAP_RULE_BP) -> dict[str, str]:
    """``locus -> group``: genes adjacent in genome order, on one contig and strand, with a gap
    of at most ``gap_bp`` (overlaps included) share a group (no joining across the origin)."""
    out: dict[str, str] = {}
    k = 0
    for i, g in enumerate(genes):
        if i:
            a = genes[i - 1]
            if not (a.contig == g.contig and a.strand == g.strand and g.start - a.end <= gap_bp):
                k += 1
        out[g.locus] = f"G{k}"
    return out


def _rdb_rows(path: Path) -> list[list[str]]:
    rows = []
    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("1)") or not line.strip():
                continue
            rows.append([x.strip() for x in line.rstrip("\n").split("\t")])
    return rows


def regulondb_tus(genes: Sequence[Feat], tu_set: Path) -> dict[str, str]:
    """``locus -> group``: connected components of same-strand genes sharing a RegulonDB TU. A
    gene matches a RegulonDB gene (``GeneProductSet.txt`` beside ``tu_set``) on the same strand
    with the same stop: the right end on ``+``, the left end on ``-`` (1-based). Unmatched genes
    and genes in no TU are singletons."""
    keys: dict[str, list[tuple[str, int]]] = {}
    for r in _rdb_rows(Path(tu_set).parent / "GeneProductSet.txt"):
        if len(r) < 6 or not r[3].isdigit() or not r[4].isdigit():
            continue
        strand = {"forward": "+", "reverse": "-"}.get(r[5].lower())
        if strand is None:
            continue
        keys.setdefault(r[1], []).append((strand, int(r[4]) if strand == "+" else int(r[3])))
    by_key: dict[tuple[str, int], list[str]] = {}
    for g in genes:
        by_key.setdefault((g.strand, g.end if g.strand == "+" else g.start + 1), []).append(g.locus)
    parent = {g.locus: g.locus for g in genes}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for r in _rdb_rows(tu_set):
        if len(r) < 5:
            continue
        by_strand: dict[str, list[str]] = {}
        for name in (n.strip() for n in r[4].split(";")):
            for key in keys.get(name, ()):
                by_strand.setdefault(key[0], []).extend(by_key.get(key, ()))
        for loci in by_strand.values():
            for x in loci[1:]:
                a, b = find(loci[0]), find(x)
                if a != b:
                    parent[max(a, b)] = min(a, b)
    return {g.locus: find(g.locus) for g in genes}


def _layout(
    genes: Sequence[Feat], groups: Mapping[str, str], seqs: Mapping[str, str], tmp: Path
) -> Any:
    """``layout_from_simulation`` on a minimal Bakta-like JSON of ``genes`` (replicons linear,
    see the module doc) and a ``truth/genes.tsv`` of ``groups``."""
    feats = [{"type": g.kind, "sequence": g.contig, "start": g.start + 1, "stop": g.end,
              "strand": g.strand, "locus": g.locus, "gene": g.extra.get("gene"),
              "product": g.extra.get("product"), "nt": g.nt, "aa": g.aa,
              "start_type": g.extra.get("start_type"), "pseudo": bool(g.extra.get("pseudo"))}
             for g in genes]  # fmt: skip
    data = {"sequences": [{"id": c, "length": len(s), "topology": "linear", "type": "chromosome"}
                          for c, s in seqs.items()], "features": feats}  # fmt: skip
    js, tsv = tmp / "layout.json", tmp / "layout_genes.tsv"
    js.write_text(json.dumps(data), encoding="utf-8")
    _write_tsv(tsv, ["locus_tag", "tu_id"], [[g.locus, groups[g.locus]] for g in genes])
    return layout_from_simulation(SimpleNamespace(bakta_json=js, genes_truth=tsv))


def _island_rows(tus: Any, cfg: ExpressionConfig) -> list[list[Any]]:
    """``truth/islands.tsv`` rows (``expression.ISLAND_COLUMNS``) of the planted islands, with the
    layout's ``TU_<n>`` ids (TU order within each replicon, as ``expression``)."""
    num = [int(str(t).rsplit("_", 1)[-1]) for t in tus["tu_id"]]
    ordered = tus.iloc[np.argsort(num, kind="stable")].reset_index(drop=True)
    by_rep: dict[str, list[str]] = {}
    for c, t in zip(ordered["contig"].astype(str), ordered["tu_id"].astype(str), strict=True):
        by_rep.setdefault(c, []).append(t)
    stationary = cfg.tu_sd / math.sqrt(1.0 - cfg.ar1_rho**2)
    rows = []
    for k, isl in enumerate(cfg.islands):
        ids = by_rep.get(isl.replicon, [])[isl.first_tu : isl.first_tu + isl.n_tus]
        rows.append([f"island_{k + 1}", isl.replicon, isl.first_tu, isl.n_tus, float(isl.effect),
                     float(isl.effect * stationary), ids[0], ids[-1]])  # fmt: skip
    return rows


# ---------------------------------------------------------------- novel ORFs

_NT = np.full(256, 4, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    _NT[_b] = _i


def _code(c: str) -> int:
    return 16 * "ACGT".index(c[0]) + 4 * "ACGT".index(c[1]) + "ACGT".index(c[2])


_STOP_CODES = np.array(sorted(_code(c) for c in STOPS))
_START_CODES = np.array(sorted(_code(c) for c in STARTS))


def _codon_codes(t: str, off: int) -> np.ndarray:
    a = _NT[np.frombuffer(t.encode(), dtype=np.uint8)]
    m = (len(t) - off) // 3
    c = a[off : off + 3 * m].reshape(m, 3)
    code = 16 * c[:, 0] + 4 * c[:, 1] + c[:, 2]
    out: np.ndarray = np.where((c == 4).any(axis=1), -1, code)
    return out


def _hit(cov: np.ndarray, s0: np.ndarray, e0: np.ndarray) -> np.ndarray:
    out: np.ndarray = cov[e0] > cov[s0]
    return out


def _translate(nt: str) -> str:
    return "".join(CODE_11.get(c, "X") for c in codons_of(nt))


def _il(s: str) -> str:
    return s.replace("I", "L")


def found_orf_candidates(
    ref: str, sample: str, genes: Sequence[Feat], rna: Sequence[Feat], len_aa: tuple[int, int]
) -> dict[str, list[tuple[int, int, str]]]:
    """Candidate novel ORFs (0-based ``s0``, ``e0``, strand) of one linear-scanned contig by
    class, in genome order: one per stop codon (its longest start..stop form in ``len_aa``), on
    ``ref`` and intact on ``sample`` (same coordinates), in frame with no overlapping same-strand
    protein gene, overlapping no RNA gene and no same-strand pseudogene. ``cds_outframe`` overlaps
    a same-strand protein gene, else ``antisense`` an opposite-strand one, else ``intergenic``."""
    n = len(ref)
    lo, hi = len_aa
    spans = [(g.start, g.end, g.strand) for g in genes]
    mask = _codon_mask(n, spans)
    cov = _cover(n, spans)
    pcov = _cover(n, [(g.start, g.end, g.strand) for g in genes if g.extra.get("pseudo")])
    rcov = _cover(n, [(f.start, f.end, "+") for f in rna])["+"]
    out: dict[str, list[tuple[int, int, str]]] = {c: [] for c in NOVEL_FOUND_CLASSES}
    for strand in "+-":
        other = "-" if strand == "+" else "+"
        tr, ts = (ref, sample) if strand == "+" else (revcomp(ref), revcomp(sample))
        for off in range(3):
            ck, cs = _codon_codes(tr, off), _codon_codes(ts, off)
            stops = np.flatnonzero(np.isin(ck, _STOP_CODES))
            starts = np.flatnonzero(np.isin(ck, _START_CODES))
            if not len(stops) or not len(starts):
                continue
            prev = np.concatenate([[-1], stops[:-1]])
            k = np.searchsorted(starts, np.maximum(prev + 1, stops - hi))
            ok = k < len(starts)
            s = starts[np.minimum(k, len(starts) - 1)]
            ok &= (s < stops) & (stops - s >= lo) & (stops - s <= hi)
            s, j = s[ok], stops[ok]
            s_stop = np.concatenate([[0], np.cumsum(np.isin(cs, _STOP_CODES))])
            ok = np.isin(cs[s], _START_CODES) & np.isin(cs[j], _STOP_CODES)
            ok &= s_stop[j] == s_stop[s]
            cm = np.concatenate([[0], np.cumsum(mask[strand][off::3][: len(ck)])])
            ok &= cm[j + 1] == cm[s]  # in frame with no protein gene
            s, j = s[ok], j[ok]
            t0, t1 = off + 3 * s, off + 3 * j + 3
            s0, e0 = (t0, t1) if strand == "+" else (n - t1, n - t0)
            keep = ~_hit(rcov, s0, e0) & ~_hit(pcov[strand], s0, e0)
            same, anti = _hit(cov[strand], s0, e0), _hit(cov[other], s0, e0)
            for cls, sel in (("cds_outframe", same), ("antisense", ~same & anti),
                             ("intergenic", ~same & ~anti)):  # fmt: skip
                sel = sel & keep
                pairs = zip(s0[sel].tolist(), e0[sel].tolist(), strict=True)
                out[cls] += [(a, b, strand) for a, b in pairs]
    for c in out:
        out[c].sort()
    return out


def _segment(seq: str, s0: int, e0: int, strand: str) -> str:
    nt = seq[s0:e0]
    return nt if strand == "+" else revcomp(nt)


def draw_found_orfs(
    rng: np.random.Generator,
    cands: Mapping[str, Sequence[tuple[int, int, str]]],
    contig: str,
    sample: str,
    hay: str,
    n: int,
) -> list[list[Any]]:
    """``n`` ORFs per class in random order whose protein (I = L, initiator dropped) is part of no
    annotated protein (``hay``); rows (class, contig, start1, end, strand, sample protein)."""
    rows: list[list[Any]] = []
    for cls in NOVEL_FOUND_CLASSES:
        pool = cands.get(cls, ())
        got = []
        for i in rng.permutation(len(pool)).tolist() if n else []:
            s0, e0, strand = pool[i]
            aa = "M" + _translate(_segment(sample, s0, e0, strand)[3:-3])
            if _il(aa[1:]) in hay:
                continue
            got.append([cls, contig, s0 + 1, e0, strand, aa])
            if len(got) == n:
                break
        if len(got) < n:
            _LOG.warning("novel ORFs: only %d %s ORFs on %s for %d requested; taking %d",
                         len(got), cls, contig, n, len(got))  # fmt: skip
        rows += sorted(got, key=lambda r: (r[2], r[3], r[4]))
    return rows


def nterm_options(g: Feat, ref: str, sample: str) -> list[int]:
    """1-based true starts (first base of the start codon in transcription direction) of ``g``:
    in-frame start codons >= ``NTERM_MIN_CODONS`` codons upstream of the annotated start with no
    stop in between, on ``ref`` and ``sample`` (linear: the walk stops at the contig end)."""
    n = len(ref)
    tr, ts = (ref, sample) if g.strand == "+" else (revcomp(ref), revcomp(sample))
    a = g.start if g.strand == "+" else n - g.end
    out = []
    p, k = a - 3, 1
    while p >= 0:
        ck, cs = tr[p : p + 3], ts[p : p + 3]
        if ck in STOPS or cs in STOPS:
            break
        if k >= NTERM_MIN_CODONS and ck in STARTS and cs in STARTS:
            out.append(p + 1 if g.strand == "+" else n - p)
        p, k = p - 3, k + 1
    return out


def draw_nterm(
    rng: np.random.Generator,
    genes: Sequence[Feat],
    sample_genes: Mapping[str, Feat],
    ref: Mapping[str, str],
    sample: Mapping[str, str],
    n: int,
) -> dict[str, int]:
    """``locus -> 1-based true start`` for up to ``n`` expressed, complete CDSs with an upstream
    start (one drawn uniformly among ``nterm_options``)."""
    pool = [g for g in genes if g.kind == "cds" and not g.pseudo_reason
            and not sample_genes[g.locus].extra.get("pseudo")
            and g.end <= len(ref[g.contig])]  # fmt: skip
    out: dict[str, int] = {}
    for i in rng.permutation(len(pool)).tolist() if n else []:
        g = pool[i]
        opts = nterm_options(g, ref[g.contig], sample[g.contig])
        if opts:
            out[g.locus] = opts[int(rng.integers(len(opts)))]
            if len(out) == n:
                break
    if len(out) < n:
        _LOG.warning("novel ORFs: only %d nterm_extension CDSs for %d requested", len(out), n)
    return out


# ---------------------------------------------------------------- outputs


def _write_tsv(path: Path, header: Sequence[str], rows: Sequence[Sequence[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _expressed(f: Feat) -> bool:
    """``proteomics.expressed_proteins``' choice of annotated proteins."""
    return f.kind in ("cds", "sorf") and not f.pseudo_reason and not f.extra.get("pseudo")


def _git_sha() -> str:
    try:
        proc = subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=30, check=False,
        )  # fmt: skip
    except OSError:
        return "unknown"
    return proc.stdout.strip() if proc.returncode == 0 else "unknown"


def _publish(stage: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for src in sorted(stage.rglob("*")):
        if src.is_dir():
            continue
        dst = outdir / src.relative_to(stage)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))


def simulate_on_genome(
    cfg: RealGenomeConfig,
    seed: int,
    outdir: Path,
    *,
    genome_fasta: Path,
    bakta_json: Path,
    tu_set: Path | None,
    check_hashes: bool = True,
) -> RealGenomeResult:
    """Simulated evidence-run inputs on a real genome (module doc); ``check_hashes=False`` lets a
    test-time excerpt through rule 1. Nothing is written when a step fails."""
    genome_fasta, bakta_json, outdir = Path(genome_fasta), Path(bakta_json), Path(outdir)
    hashes = verify_inputs(genome_fasta, bakta_json) if check_hashes else {}
    raw_seqs, topology_raw = _read_fasta(genome_fasta)
    cmap = _contig_map(bakta_json)
    ref = {cmap.get(c, c): s for c, s in raw_seqs.items()}
    topology = {cmap.get(c, c): t for c, t in topology_raw.items()}
    hashes.setdefault("genome_file", _sha256(genome_fasta))
    hashes.setdefault("genome_sequence", _sequence_sha(raw_seqs))
    hashes.setdefault("annotation_file", _sha256(bakta_json))
    data = json.loads(_maybe_gunzip(bakta_json.read_bytes()))
    genes, rna = _feats(data, ref)

    # step 2: the sample is derive_reference(genome); SNVs only, so coordinates agree
    sample, sample_feats, strain_rows = derive_reference(cfg.strain, seed, ref, genes)
    if any(len(sample[c]) != len(ref[c]) for c in ref):
        raise ValueError("the sample strain changed a contig length (SNVs only)")
    by_locus = {f.locus: f for f in sample_feats}

    # step 3: TUs and abundance
    tu_hashes: dict[str, str] = {}
    if cfg.tu_source == "regulondb" and tu_set is not None and Path(tu_set).is_file():
        groups, tu_source = regulondb_tus(genes, Path(tu_set)), "regulondb"
        gps = Path(tu_set).parent / "GeneProductSet.txt"
        tu_hashes = {"tu_set": _sha256(Path(tu_set)), "gene_product_set": _sha256(gps)}
    else:
        groups, tu_source = gap_rule_tus(genes), "gap_rule"
        _LOG.warning("tu_source = gap_rule: no RegulonDB TUSet.txt (adjacent same-strand genes "
                     "with a gap <= %d bp share a TU)", GAP_RULE_BP)  # fmt: skip
    outdir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".giapsim-real-", dir=outdir.parent) as tmp:
        stage = Path(tmp) / "out"
        (stage / "truth").mkdir(parents=True)
        lay = _layout(genes, groups, ref, Path(tmp))
        erng = np.random.default_rng([seed, _STREAM, 0])
        ln = expression_for_simulation(lay.genes, lay.tus, cfg.expression, erng)
        tu_of = {str(g): str(t) for t, gs in zip(lay.tus["tu_id"], lay.tus["genes"], strict=True)
                 for g in gs}  # fmt: skip
        value = dict(zip(lay.genes["locus_tag"].astype(str), ln.tolist(), strict=True))

        # step 4: novel ORFs, found on the unchanged genome
        nrng = np.random.default_rng([seed, _STREAM, 1])
        hay = "|".join(_il(g.aa) for g in genes)
        found: list[list[Any]] = []
        for cid in ref:
            on = [g for g in genes if g.contig == cid]
            rna_on = [f for f in rna if f.contig == cid]
            cands = found_orf_candidates(ref[cid], sample[cid], on, rna_on, cfg.novel_len_aa)
            found += draw_found_orfs(nrng, cands, cid, sample[cid], hay, cfg.n_novel_per_class)
        novel = [[f"NOVEL_{i + 1:04d}", *r] for i, r in enumerate(found)]
        nterm = draw_nterm(nrng, genes, by_locus, ref, sample, cfg.n_novel_per_class)
        feats = [replace(f, true_start=nterm[f.locus], extra=dict(f.extra))
                 if f.locus in nterm else f for f in sample_feats]  # fmt: skip

        # abundance: expressed genes centred to LOG_ABUNDANCE[0], novel ORFs upper tercile
        expressed = [f for f in feats if _expressed(f)]
        gv = np.array([value[f.locus] for f in expressed], dtype=float)
        gv = np.round(gv - (gv.mean() if len(gv) else 0.0) + LOG_ABUNDANCE[0], 6)
        cut = float(np.quantile(gv, 2 / 3)) if len(gv) else 0.0
        upper = np.sort(gv[gv >= cut])
        arng = np.random.default_rng([seed, _STREAM, 2])
        nv = upper[arng.integers(len(upper), size=len(novel))] if len(upper) else np.zeros(0)
        la = np.concatenate([gv, nv])

        # step 5: peptides
        paths = simulate_proteomics(cfg.proteomics, seed, stage, feats, novel, sample, topology,
                                    log_abundance=la)  # fmt: skip
        with paths["peptides"].open(encoding="utf-8") as fh:
            n_target = sum(1 for r in csv.DictReader(fh, delimiter="\t")
                           if r["is_decoy"] == "false")  # fmt: skip
        if n_target == 0:
            raise ValueError("no_target_peptides: refusing a peptides.tsv without target rows (D9)")
        paths["proteins"].unlink()  # the sample's proteome is not a GIAP input

        # truth tables
        truth = stage / "truth"
        peps = _truth_peptides(paths["truth"])
        alt = [r for r in alt_site_rows(strain_rows, feats, sample, ref, peps)
               if not _initiator(by_locus[str(r[3])], int(r[1]))]  # fmt: skip
        _write_tsv(truth / "alt_sites.tsv", ALT_COLUMNS, alt)
        nrows = novel_rows(novel, feats, sample, set(), peps)
        _write_tsv(truth / "novel_orfs.tsv", NOVEL_COLUMNS, nrows)
        requested = {c: cfg.n_novel_per_class for c in (*NOVEL_FOUND_CLASSES, "nterm_extension")}
        report = novel_report(requested, nrows)
        (truth / "novel_orfs.json").write_text(json.dumps(report, indent=2) + "\n")
        _write_tsv(truth / "genes.tsv", ["locus_tag", "tu_id", "contig", "start", "end", "strand"],
                   [[g.locus, tu_of[g.locus], g.contig, g.start + 1, g.end, g.strand]
                    for g in genes])  # fmt: skip
        islands = _island_rows(lay.tus, cfg.expression)
        _write_tsv(truth / "islands.tsv", ISLAND_COLUMNS, islands)
        (truth / "sample_genome.fna").write_text("".join(
            f">{c} sample strain (giapsim MVP-04; debugging only, never a GIAP input)\n"
            + "".join(s[i : i + 80] + "\n" for i in range(0, len(s), 80))
            for c, s in sample.items()))  # fmt: skip
        (stage / "config.yaml").write_text(f"seed: {int(seed)}\nsample_name: {cfg.sample_name}\n")
        n_snv = sum(1 for r in strain_rows if r[5] == "snv")
        sizes: dict[str, int] = {}
        for t in tu_of.values():
            sizes[t] = sizes.get(t, 0) + 1
        n_multi = sum(v > 1 for v in sizes.values())
        prov: dict[str, Any] = {
            "generator": "giapsim.real_genome (MVP-04 part 1)",
            "giapsim_git_sha": _git_sha(),
            "seed": int(seed),
            "b0": float(cfg.proteomics.detect_logit[0]),
            "detect_logit": [float(x) for x in cfg.proteomics.detect_logit],
            "tu_source": tu_source,
            "decoys": "unscored",
            "decoy_collisions": 0,
            "inputs": {
                "genome": {
                    "name": genome_fasta.name,
                    "sha256_file": hashes["genome_file"],
                    "sha256_sequence_upper": hashes["genome_sequence"],
                    "hash_checked": check_hashes,
                },
                "annotation": {"name": bakta_json.name, "sha256_file": hashes["annotation_file"]},
                "tu_set": tu_hashes or None,
            },
            "counts": {
                "n_snv": n_snv,
                "n_alt_sites": len(alt),
                "n_tus": len(sizes),
                "n_multi_gene_tus": n_multi,
                "n_islands": len(islands),
                "n_target_rows": n_target,
                "novel_orfs": report,
            },
        }
        (stage / "provenance.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")
        _publish(stage, outdir)
    rel = {k: outdir / p.relative_to(stage) for k, p in paths.items() if k != "proteins"}
    for name in ("alt_sites", "novel_orfs", "genes", "islands"):
        rel[name] = outdir / "truth" / f"{name}.tsv"
    rel["provenance"], rel["config"] = outdir / "provenance.json", outdir / "config.yaml"
    rel["sample_genome"] = outdir / "truth" / "sample_genome.fna"
    return RealGenomeResult(outdir, rel, prov)


def _initiator(g: Feat, codon_start1: int) -> bool:
    """Whether the 1-based lowest base ``codon_start1`` is ``g``'s annotated initiator codon."""
    return codon_start1 - 1 == (g.start if g.strand == "+" else g.end - 3)


def _truth_peptides(path: Path) -> dict[str, list[tuple[str, int, int, str]]]:
    """Correct database peptides of ``truth/peptides.tsv`` by locus (``scenario_truth`` shape)."""
    by: dict[str, list[tuple[str, int, int, str]]] = {}
    with path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["origin"] == "decoy" or r["is_wrong"] == "true":
                continue
            by.setdefault(r["locus_tag"], []).append(
                (r["contig"], int(r["start"]), int(r["end"]), _il(r["peptide"]))
            )
    return by


# ---------------------------------------------------------------- CLI


def main(argv: Sequence[str] | None = None) -> int:
    """The input step: ``scenarios.k12_evidence()`` at ``--seed`` with the hash check on."""
    ap = argparse.ArgumentParser(prog="python -m giapsim.real_genome", description=__doc__.split(
        "\n", 1)[0])  # fmt: skip
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--annotation", type=Path, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--tu-set", type=Path, default=None, help="RegulonDB TUSet.txt (optional)")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    from giapsim.scenarios import k12_evidence

    try:
        res = simulate_on_genome(k12_evidence(), args.seed, args.outdir, genome_fasta=args.genome,
                                 bakta_json=args.annotation, tu_set=args.tu_set)  # fmt: skip
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(res.provenance["counts"], sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
