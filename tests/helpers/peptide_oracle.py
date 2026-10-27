"""Independent oracle for the gippasim proteomics tests (issue M3-19; test-author code).

Not part of ``gippasim``: it re-derives everything the simulated peptide tables claim from the
genome itself (FASTA + ``truth/features.tsv`` + ``truth/novel_orfs.tsv``), so that the tests check
the simulator against the genome rather than against its own internals.

* ``trypsin_p`` is a reference Trypsin/P digest (cleave after every K/R, also before P).
* ``true_proteins`` translates every expressed true protein from the genome: non-pseudo CDS and
  sORFs from their *true* start (N-terminal extensions included) and the planted/unannotated
  novel ORFs; the initiator codon reads as M, the stop codon is dropped.
* ``maps_to`` checks that a peptide is the translation of a 1-based closed genome span.
* ``sixframe_count`` counts the distinct genomic loci whose six-frame translation is a peptide.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any

BASES = "TCAG"
_AA_11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE: dict[str, str] = {
    a + b + c: _AA_11[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}
STARTS = frozenset({"ATG", "GTG", "TTG"})
STOPS = frozenset({"TAA", "TAG", "TGA"})
EMPTY = frozenset({"", "NA", "na", "NaN", "nan", "."})
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def translate(nt: str) -> str:
    return "".join(CODE.get(nt[i : i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def il(s: str) -> str:
    return s.replace("I", "L")


def empty(v: Any) -> bool:
    return v is None or (isinstance(v, float) and math.isnan(v)) or str(v).strip() in EMPTY


def truthy(v: Any) -> bool:
    return str(v).strip().lower() in ("true", "1")


# ---------------------------------------------------------------- digest


def trypsin_p(
    protein: str, length: tuple[int, int] = (7, 30), missed: int = 1
) -> set[tuple[int, int]]:
    """Spans ``[a, b)`` of Trypsin/P products: cleave after K/R (also before P)."""
    cuts = [0] + [i + 1 for i, c in enumerate(protein[:-1]) if c in "KR"] + [len(protein)]
    cuts = sorted(set(cuts))
    out = set()
    for i in range(len(cuts) - 1):
        for m in range(missed + 1):
            j = i + 1 + m
            if j >= len(cuts):
                break
            a, b = cuts[i], cuts[j]
            if length[0] <= b - a <= length[1]:
                out.add((a, b))
    return out


def trypsin_p_peptides(
    protein: str, length: tuple[int, int] = (7, 30), missed: int = 1
) -> set[str]:
    return {protein[a:b] for a, b in trypsin_p(protein, length, missed)}


# ---------------------------------------------------------------- files


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif cur is not None and line.strip():
            out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


@dataclass(frozen=True)
class Protein:
    locus: str
    kind: str  # cds | sorf | novel:<class>
    contig: str
    start: int  # 0-based, unwrapped (end may exceed the contig length across the origin)
    end: int
    strand: str
    aa: str
    annotated_start: int  # 0-based annotated span (differs from start/end for N-ext genes)
    annotated_end: int


@dataclass
class Run:
    out: Path
    meta: dict[str, Any] = field(default_factory=dict)

    @cached_property
    def genome(self) -> dict[str, str]:
        return read_fasta(self.out / "genome.fasta")

    @cached_property
    def bakta(self) -> dict[str, Any]:
        data: dict[str, Any] = json.loads((self.out / "bakta_like.json").read_text())
        return data

    @cached_property
    def topology(self) -> dict[str, str]:
        return {s["id"]: s.get("topology", "circular") for s in self.bakta["sequences"]}

    @cached_property
    def features(self) -> list[dict[str, str]]:
        return rows(self.out / "truth" / "features.tsv")

    @cached_property
    def novel(self) -> list[dict[str, str]]:
        p = self.out / "truth" / "novel_orfs.tsv"
        return rows(p) if p.exists() else []

    @cached_property
    def peptides(self) -> list[dict[str, str]]:
        return rows(self.out / "peptides.tsv")

    @cached_property
    def truth(self) -> list[dict[str, str]]:
        return rows(self.out / "truth" / "peptides.tsv")

    @cached_property
    def truth_by_peptide(self) -> dict[str, list[dict[str, str]]]:
        d: dict[str, list[dict[str, str]]] = {}
        for r in self.truth:
            d.setdefault(r["peptide"], []).append(r)
        return d

    @cached_property
    def json_aa(self) -> dict[str, str]:
        return {
            f["locus"]: f["aa"]
            for f in self.bakta["features"]
            if f["type"] in ("cds", "sorf") and f.get("aa")
        }

    @cached_property
    def pseudo_loci(self) -> set[str]:
        return {f["feature_id"] for f in self.features if not empty(f.get("pseudo_reason"))}

    @cached_property
    def proteins(self) -> dict[str, Protein]:
        return true_proteins(self)

    @cached_property
    def proteome_il(self) -> str:
        """Every true protein (I/L-equivalent), ``|``-separated, for substring tests."""
        return "|".join(il(p.aa) for p in self.proteins.values())

    @cached_property
    def proteome_and_annotation_il(self) -> str:
        """True proteins plus every annotated protein (I/L-equivalent)."""
        return self.proteome_il + "|" + "|".join(il(a) for a in self.json_aa.values())

    def decoy_rows(self) -> list[dict[str, str]]:
        return [r for r in self.peptides if truthy(r.get("is_decoy"))]

    def target_rows(self) -> list[dict[str, str]]:
        return [r for r in self.peptides if not truthy(r.get("is_decoy"))]

    def is_wrong(self, peptide: str) -> bool:
        rs = self.truth_by_peptide[peptide]
        return any(truthy(r["is_wrong"]) for r in rs if r["origin"] != "decoy")


def genome_slice(seq: str, s0: int, e0: int, strand: str) -> str:
    n = len(seq)
    if s0 < 0:
        s0, e0 = s0 + n, e0 + n
    if s0 >= n:
        s0, e0 = s0 - n, e0 - n
    part = seq[s0:e0] if e0 <= n else seq[s0:] + seq[: e0 - n]
    return revcomp(part) if strand == "-" else part


def _protein_from_genome(seq: str, s0: int, e0: int, strand: str) -> str:
    nt = genome_slice(seq, s0, e0, strand)
    aa = translate(nt)
    assert nt[:3] in STARTS, f"no start codon at {s0}-{e0}{strand}: {nt[:3]}"
    assert nt[-3:] in STOPS, f"no stop codon at {s0}-{e0}{strand}: {nt[-3:]}"
    return "M" + aa[1:-1]


def true_proteins(run: Run) -> dict[str, Protein]:
    out: dict[str, Protein] = {}
    for f in run.features:
        if f["type"] not in ("cds", "sorf") or not empty(f.get("pseudo_reason")):
            continue
        s0, e0, strand = int(f["start"]) - 1, int(f["end"]), f["strand"]
        ts, te = s0, e0
        if not empty(f.get("true_start")):
            t = int(f["true_start"])
            if strand == "+":
                ts = t - 1
            else:
                te = t
        seq = run.genome[f["contig"]]
        aa = _protein_from_genome(seq, ts, te, strand)
        out[f["feature_id"]] = Protein(
            f["feature_id"], f["type"], f["contig"], ts, te, strand, aa, s0, e0
        )
    for r in run.novel:
        s0, e0 = int(r["start"]) - 1, int(r["end"])
        aa = _protein_from_genome(run.genome[r["contig"]], s0, e0, r["strand"])
        assert aa == r["aa"], (r["orf_id"], aa, r["aa"])
        out[r["orf_id"]] = Protein(
            r["orf_id"], "novel:" + r["class"], r["contig"], s0, e0, r["strand"], aa, s0, e0
        )
    return out


def maps_to(peptide: str, genome: dict[str, str], contig: str, start1: int, end1: int,
            strand: str) -> bool:  # fmt: skip
    """``peptide`` is the translation of ``contig:start1-end1`` (1-based closed) on ``strand``;
    an initiator M may be read from a GTG/TTG start codon."""
    if (end1 - start1 + 1) != 3 * len(peptide):
        return False
    nt = genome_slice(genome[contig], start1 - 1, end1, strand)
    aa = translate(nt)
    if aa[1:] != peptide[1:]:
        return False
    return aa[0] == peptide[0] or (peptide[0] == "M" and nt[:3] in STARTS)


_WRAP_AA = 40  # longest peptide the six-frame search handles across the origin
_FRAME_CACHE: dict[tuple[int, str], list[tuple[str, str, int, str, str, int, int]]] = {}


def _frames(
    contig: str, seq: str, circular: bool
) -> list[tuple[str, str, int, str, str, int, int]]:
    """(contig, strand, offset, aa, aa with I->L, length, extended length), cached per sequence."""
    key = (hash(seq), contig)
    if key not in _FRAME_CACHE:
        ext = seq + seq[: 3 * _WRAP_AA] if circular else seq
        rc = revcomp(ext)
        out = []
        for off in range(3):
            for strand, src in (("+", ext), ("-", rc)):
                aa = translate(src[off:])
                out.append((contig, strand, off, aa, il(aa), len(seq), len(ext)))
        if len(_FRAME_CACHE) > 64:
            _FRAME_CACHE.clear()
        _FRAME_CACHE[key] = out
    return _FRAME_CACHE[key]


def sixframe_loci(peptide: str, genome: dict[str, str], topology: dict[str, str],
                  il_equal: bool) -> set[tuple[str, int, str]]:  # fmt: skip
    """Distinct loci (contig, 0-based plus-strand start, strand) translating to ``peptide``."""
    assert len(peptide) <= _WRAP_AA
    q = il(peptide) if il_equal else peptide
    loci = set()
    for contig, seq in genome.items():
        circular = topology.get(contig, "circular") == "circular"
        for _, strand, off, aa, aa_il, n, ext_len in _frames(contig, seq, circular):
            t = aa_il if il_equal else aa
            i = t.find(q)
            while i != -1:
                s = off + 3 * i if strand == "+" else ext_len - off - 3 * (i + len(q))
                if circular:
                    s %= n
                loci.add((contig, s, strand))
                i = t.find(q, i + 1)
    return loci


# ---------------------------------------------------------------- statistics


def phi(x: float) -> float:
    return math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)


def sf(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2))


def mixture_pep(s: float, pi: float, mix: tuple[float, float, float, float]) -> float:
    mu1, sd1, mu0, sd0 = mix
    f0 = phi((s - mu0) / sd0) / sd0
    f1 = phi((s - mu1) / sd1) / sd1
    return pi * f0 / (pi * f0 + (1 - pi) * f1)


def mixture_fdr(s: float, pi: float, mix: tuple[float, float, float, float]) -> float:
    mu1, sd1, mu0, sd0 = mix
    s0 = sf((s - mu0) / sd0)
    s1 = sf((s - mu1) / sd1)
    return pi * s0 / (pi * s0 + (1 - pi) * s1)


def fnum(v: Any) -> float | None:
    return None if empty(v) else float(v)


# ---------------------------------------------------------------- GIPPA's own exact matcher


def gippa_exact_counts(run: Run, peptides: list[str]) -> dict[str, int]:
    """Distinct target loci per peptide as GIPPA's merged exact matcher (M3-02 six-frame index +
    ``gippa.mapping.exact``, I = L) finds them; origin copies are counted once."""
    import numpy as np

    from gippa.mapping.decoys import build_decoys
    from gippa.mapping.exact import match_exact
    from gippa.mapping.sixframe import build_sixframe
    from gippa.model.coords import Replicon

    reps = {}
    for s in run.bakta["sequences"]:
        kind = s.get("type") if s.get("type") in ("chromosome", "plasmid") else "unknown"
        reps[s["id"]] = Replicon(s["id"], int(s["length"]), s["topology"], kind)
    idx = build_sixframe(run.genome, reps)
    dec = {"reverse": build_decoys(idx, "reverse", 0, np.random.default_rng(0))}
    peps = sorted(set(peptides))
    hits = match_exact(peps, idx, dec, dict.fromkeys(peps, "reverse"))
    loci: dict[str, set[tuple[str, int, int, str]]] = {p: set() for p in peps}
    for h in hits[hits["search"] == "target"].itertuples():
        n = reps[h.contig].length
        s = int(h.start) % n if reps[h.contig].topology == "circular" else int(h.start)
        loci[h.peptide].add((h.contig, s, int(h.end) - int(h.start), h.strand))
    return {p: len(v) for p, v in loci.items()}


def find_truth_meta(out: Path, key: str = "pi_wrong_realised") -> dict[str, Any]:
    """The truth metadata mapping (JSON or YAML under ``out/truth`` or ``out``) holding ``key``."""
    import yaml

    def walk(d: Any) -> dict[str, Any] | None:
        if isinstance(d, dict):
            if key in d:
                return d
            for v in d.values():
                got = walk(v)
                if got is not None:
                    return got
        return None

    for p in sorted([*(out / "truth").glob("*"), *out.glob("*")]):
        if p.suffix not in (".json", ".yaml", ".yml") or not p.is_file():
            continue
        try:
            data = (
                json.loads(p.read_text()) if p.suffix == ".json" else yaml.safe_load(p.read_text())
            )
        except Exception:
            continue
        got = walk(data)
        if got is not None:
            return got
    raise AssertionError(f"no truth metadata with {key!r} under {out}")
