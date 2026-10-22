"""Hand-built genomes, C9/C11 frames, mirrors and a brute-force ORF oracle for M7-06 tests.

Independent of ``giap.novel``: the oracle walks codon by codon on the strand-oriented sequence
("transcription coordinates": the plus string for ``+``, its reverse complement for ``-``) and maps
back to 0-based genome coordinates (docs/specs/M7-06-novel-candidates-and-orfs.md rules 2-6).

Conventions pinned here (and in the tests that import them):

* a codon is reported by its lowest 0-based genome coordinate, wrapped into ``[0, L)``;
  on ``-`` the codon ``[x, x + 3)`` of the reverse complement maps to ``[L - x - 3, L - x)``;
* ``starts`` are listed 5'->3' along the ORF's strand (descending coordinate on ``-``);
* the ORF ID is ``NOVEL_<contig>_<lowest + 1>_<p|m>`` with ``lowest`` the lowest genome coordinate
  of the longest interval (from ``start_longest``, else the first peptide, to the stop codon);
* C11 ``frame`` is the M3 label (``+1..+3`` from ``start % 3``, ``-1..-3`` from ``(L - end) % 3``).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel

COMP = str.maketrans("ACGT", "TGCA")
# Stops in all three phases on both strands every 22 nt (as tests/unit/test_loci.py).
WALL = "TAAATAAATAATTATTTATTTA"
FILLER = "GCA"  # no stop and no ATG/GTG/TTG in its own phase on +; no stop in any phase on -
STOPS = frozenset({"TAA", "TAG", "TGA"})
START_CODONS: tuple[str, ...] = ("ATG", "GTG", "TTG")
CANDIDATE_CLASSES: tuple[str, ...] = (
    "nterm_extension", "cterm_extension", "intergenic", "cds_outframe", "antisense",
    "pseudogene", "sorf",
)  # fmt: skip

_BASES = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODON = {
    a + b + c: _AA11[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def il(s: str) -> str:
    return s.upper().replace("I", "L")


def translate(nt: str) -> str:
    return "".join(_CODON.get(nt[i : i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def wrapped(seq: str, lo: int, hi: int) -> str:
    """``seq[lo, hi)`` read circularly (``lo`` may be negative, ``hi`` may exceed ``len``)."""
    n = len(seq)
    return "".join(seq[p % n] for p in range(lo, hi))


def interval_nt(seq: str, start: int, end: int, strand: str) -> str:
    """The 5'->3' nucleotides of ``[start, end)`` on ``strand`` (circular read)."""
    s = wrapped(seq, start, end)
    return s if strand == "+" else revcomp(s)


class Seq:
    """A mutable genome: WALL background, stop-free filler and planted codons (wrapping)."""

    def __init__(self, n: int, background: str | None = None) -> None:
        self.n = n
        self.b = list(background) if background is not None else list(
            (WALL * (n // len(WALL) + 1))[:n]
        )  # fmt: skip

    def put(self, start: int, s: str) -> Seq:
        for k, c in enumerate(s):
            self.b[(start + k) % self.n] = c
        return self

    def fill(self, start: int, end: int) -> Seq:
        """Filler on ``[start, end)``; its codon phase starts at ``start``."""
        for k, i in enumerate(range(start, end)):
            self.b[i % self.n] = FILLER[k % 3]
        return self

    def text(self) -> str:
        return "".join(self.b)


def frame_label(start: int, end: int, strand: str, length: int) -> int:
    if strand == "+":
        return 1 + start % 3
    return -(1 + ((length - end) % length) % 3)


# ---------------------------------------------------------------- C9 / C11 rows


def hit_row(
    peptide: str,
    contig: str,
    start: int,
    end: int,
    strand: str,
    locus_class: str,
    length: int,
    *,
    gene: str | None = None,
    uniqueness: str = "unique",
    route: str = "exact",
    flags: str = "",
) -> dict[str, Any]:
    """One C11 ``peptide_hits`` row (contract columns of giap.genetable.contract)."""
    return {
        "peptide": peptide,
        "peptide_len": len(peptide),
        "source": "database",
        "best_q_input": 0.001,
        "pep_input": 0.01,
        "score_input": 50.0,
        "contig": contig,
        "start": start,
        "end": end,
        "strand": strand,
        "frame": frame_label(start, end, strand, length),
        "crosses_origin": end > length,
        "route": route,
        "n_diff": 0,
        "diff_types": [],
        "placement_q": 0.001,
        "length_class": None,
        "locus_class": locus_class,
        "uniqueness": uniqueness,
        "gene": gene,
        "novel_candidate_id": None,
        "flags": flags,
    }


HIT_COLUMNS: tuple[str, ...] = tuple(hit_row("A", "c", 0, 3, "+", "intergenic", 10))


def hits_frame(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(list(rows), columns=list(HIT_COLUMNS))
    return df.astype({"start": "int64", "end": "int64", "peptide_len": "int64"})


def peptides_frame(
    peptides: Iterable[str],
    q: Mapping[str, float | None] | None = None,
    contaminant: Iterable[str] = (),
) -> pd.DataFrame:
    """C9 ``peptides`` rows: one per distinct peptide string (input order kept)."""
    q = q or {}
    con = set(contaminant)
    rows = []
    for p in dict.fromkeys(peptides):
        qv = q.get(p, 0.001)
        rows.append(
            {
                "peptide": p,
                "q_value": float("nan") if qv is None else qv,
                "pep": 0.01,
                "score": 50.0,
                "is_decoy": False,
                "is_contaminant": p in con,
                "protein_ids": [],
                "modified_peptides": [p],
                "n_rows": 1,
                "n_runs": 1,
                "intensity": 1.0e6,
                "source": "database" if qv is not None else "denovo",
                "search_db": None,
            }
        )
    cols = ["peptide", "q_value", "pep", "score", "is_decoy", "is_contaminant", "protein_ids",
            "modified_peptides", "n_rows", "n_runs", "intensity", "source",
            "search_db"]  # fmt: skip
    return pd.DataFrame(rows, columns=cols)


def cds(tag: str, contig: str, start: int, end: int, strand: str, kind: str = "cds") -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type=kind,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" + "A" * max((end - start) // 3 - 2, 0),
        nt=None,
        start_codon="ATG",
    )


@dataclass
class Case:
    """A genome, its model and the C11/C9 inputs of one test scenario."""

    genome: dict[str, str]
    replicons: dict[str, Replicon]
    hits: list[dict[str, Any]]
    genes: list[Gene] = field(default_factory=list)
    q: dict[str, float | None] = field(default_factory=dict)
    contaminant: tuple[str, ...] = ()

    @property
    def gm(self) -> GeneModel:
        return GeneModel(self.replicons, tuple(self.genes), ())

    def hits_df(self) -> pd.DataFrame:
        return hits_frame(self.hits)

    def peptides_df(self) -> pd.DataFrame:
        return peptides_frame([h["peptide"] for h in self.hits], self.q, self.contaminant)


# ---------------------------------------------------------------- mirrors


def mirror_interval(start: int, end: int, length: int) -> tuple[int, int]:
    """``[a, b)`` -> ``[L - b, L - a)``, shifted by ``+L`` when it would start below 0."""
    a, b = length - end, length - start
    if a < 0:
        a, b = a + length, b + length
    return a, b


def mirror_codon(x: int, length: int) -> int:
    """Lowest coordinate of the mirror of the codon ``[x, x + 3)``."""
    return (length - x - 3) % length


def mirror_case(case: Case) -> Case:
    """Reverse-complement every contig; mirror hits and genes (strands swapped)."""
    genome = {c: revcomp(s) for c, s in case.genome.items()}
    flip = {"+": "-", "-": "+"}
    hits = []
    for h in case.hits:
        n = case.replicons[h["contig"]].length
        a, b = mirror_interval(h["start"], h["end"], n)
        hits.append(
            hit_row(h["peptide"], h["contig"], a, b, flip[h["strand"]], h["locus_class"], n,
                    gene=h["gene"], uniqueness=h["uniqueness"], route=h["route"],
                    flags=h["flags"])
        )  # fmt: skip
    genes = []
    for g in case.genes:
        iv = g.interval
        n = case.replicons[iv.contig].length
        a, b = mirror_interval(iv.start, iv.end, n)
        genes.append(cds(g.locus_tag, iv.contig, a, b, flip[iv.strand], g.feature_type))
    return Case(genome, dict(case.replicons), hits, genes, dict(case.q), case.contaminant)


# ---------------------------------------------------------------- oracle


@dataclass(frozen=True)
class OracleOrf:
    orf_id: str
    contig: str
    strand: str
    stop: int | None  # None when the 3' end is open
    starts: tuple[int, ...]
    start_longest: int | None
    end_5p_open: bool
    end_3p_open: bool
    peptides: tuple[str, ...]
    # transcription coordinates (for window/OSTIR checks)
    t_starts: tuple[int, ...]
    t_first: int


def _codon(seq: str, p: int, circular: bool) -> str | None:
    n = len(seq)
    if not circular and (p < 0 or p + 3 > n):
        return None
    return wrapped(seq, p, p + 3)


def oracle_orfs(
    hits: Sequence[Mapping[str, Any]],
    genome: Mapping[str, str],
    replicons: Mapping[str, Replicon],
    start_codons: Sequence[str] = START_CODONS,
    max_walk: int | None = None,
) -> list[OracleOrf]:
    """Rules 1-6 on intergenic-only rows: drop multi-locus peptides, walk, group, ID.

    Raises ``ValueError`` when a circular frame has no stop for a whole turn, or when a walk
    exceeds ``max_walk`` nt (property tests use a margin below one turn, so the exact turn
    boundary, which the spec leaves to the implementation, is never probed).
    """
    loci: dict[str, set[tuple[str, int, int, str]]] = {}
    for h in hits:
        loci.setdefault(h["peptide"], set()).add((h["contig"], h["start"], h["end"], h["strand"]))
    groups: dict[tuple[Any, ...], list[tuple[int, int, str]]] = {}
    tseqs: dict[tuple[str, str], str] = {}
    for pep, ls in loci.items():
        if len(ls) != 1:
            continue
        ((contig, s, e, strand),) = ls
        rep = replicons[contig]
        n = rep.length
        circ = rep.topology == "circular"
        tseq = tseqs.setdefault(
            (contig, strand), genome[contig] if strand == "+" else revcomp(genome[contig])
        )
        ts, te = (s, e) if strand == "+" else mirror_interval(s, e, n)
        p = te
        stop_t: int | None = None
        while True:
            cod = _codon(tseq, p, circ)
            if cod is None:
                end_t = p
                break
            if cod in STOPS:
                stop_t, end_t = p, p + 3
                break
            p += 3
            if p - te > (n + 3 if max_walk is None else max_walk):
                raise ValueError(f"no stop codon on {contig} {strand} frame")
        key = (contig, strand, ("stop", stop_t % n) if stop_t is not None else ("open", end_t))
        groups.setdefault(key, []).append((ts, end_t, pep))
    out = []
    for (contig, strand, (kind, _)), members in groups.items():
        rep = replicons[contig]
        n = rep.length
        circ = rep.topology == "circular"
        tseq = tseqs[(contig, strand)]
        first, end_t, _ = max(members, key=lambda m: m[1] - m[0])
        q = first - 3
        open5 = False
        while True:
            cod = _codon(tseq, q, circ)
            if cod is None:
                open5 = True
                break
            if cod in STOPS:
                break
            q -= 3
            if first - q > (n + 3 if max_walk is None else max_walk):
                raise ValueError(f"no stop codon on {contig} {strand} frame")
        t_starts = tuple(
            x for x in range(q + 3, first + 1, 3) if _codon(tseq, x, circ) in start_codons
        )
        lo_t = t_starts[0] if t_starts else first
        stop_t = end_t - 3 if kind == "stop" else None
        if strand == "+":
            stop = None if stop_t is None else stop_t % n
            starts = tuple(x % n for x in t_starts)
            lowest = lo_t % n
        else:
            stop = None if stop_t is None else mirror_codon(stop_t, n)
            starts = tuple(mirror_codon(x, n) for x in t_starts)
            lowest = (n - end_t) % n
        peps = tuple(sorted({il(m[2]) for m in members}))
        out.append(
            OracleOrf(
                orf_id=f"NOVEL_{contig}_{lowest + 1}_{'p' if strand == '+' else 'm'}",
                contig=contig,
                strand=strand,
                stop=stop,
                starts=starts,
                start_longest=starts[0] if starts else None,
                end_5p_open=open5,
                end_3p_open=kind == "open",
                peptides=peps,
                t_starts=t_starts,
                t_first=first,
            )
        )
    return out


def orf_view(o: Any) -> tuple[Any, ...]:
    """The fields an ``Orf`` and an ``OracleOrf`` share (stop omitted when the 3' end is open)."""
    return (
        o.orf_id, o.contig, o.strand, None if o.end_3p_open else o.stop, tuple(o.starts),
        o.start_longest, bool(o.end_5p_open), bool(o.end_3p_open), tuple(o.peptides),
    )  # fmt: skip
