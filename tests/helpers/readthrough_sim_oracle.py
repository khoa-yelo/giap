"""An independent oracle for the giapsim read-through knob (M7-16 rule 14; PR 16b-3).

Written from the spec, not from ``giapsim``: it reads only the files a build writes
(``genome.fasta``, ``truth/features.tsv``, ``bakta_like.json``, ``truth/peptides.tsv``,
``truth/readthrough.tsv``) and the table-11 code. Interpretations pinned here (reported to the
lead as spec ambiguities):

* **Eligible genes.** Non-pseudo ``cds`` and ``sorf`` features (rule 1's feature set; a pseudogene
  is not expressed). The terminal codon is the gene's last codon (``codon_start`` as in the
  Definitions: lowest 0-based coordinate, wrapped). Scanning on the gene strand, codon ``k``
  after the annotated stop starts ``3k`` nt further in the gene direction. The gene is eligible
  when the first stop of that scan is at ``5 <= k <= 60``; on a linear replicon every scanned
  codon must lie inside ``[0, n)``, and a scan that leaves the replicon first is not eligible; on
  a circular replicon the scan wraps.
* **ext_stop** is the new stop codon's lowest 0-based coordinate, wrapped (the ``codon_start``
  convention): ``codon_start + 3k`` on ``+``, ``codon_start - 3k`` on ``-``, both mod ``n``.
* **The expressed protein** is the annotated protein, then the residue read at the stop (``W`` at
  TGA, ``Q`` at TAA/TAG), then the ``k - 1`` sense codons up to the new stop.
* **n_true_peptides_spanning** counts the distinct correct (``is_wrong`` false) peptides of
  ``truth/peptides.tsv`` with a row at the gene whose interval holds the annotated stop codon.
  They are the detected peptides (the only ones a run can see); since the residue read at the
  stop is never K/R and the extension holds >= 4 more residues, a peptide over the stop always
  extends past it.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from giapsim.codons import CODE_11, STOPS

READTHROUGH_COLUMNS = [
    "locus_tag",
    "contig",
    "strand",
    "codon_start",
    "codon",
    "residue",
    "ext_stop",
    "n_true_peptides_spanning",
]
RESIDUE = {"TGA": "W", "TAA": "Q", "TAG": "Q"}
K_MIN, K_MAX = 5, 60


@dataclass(frozen=True)
class Genome:
    seqs: Mapping[str, str]
    circular: Mapping[str, bool]

    def n(self, contig: str) -> int:
        return len(self.seqs[contig])


@dataclass(frozen=True)
class Eligible:
    locus_tag: str
    kind: str  # cds | sorf
    contig: str
    strand: str
    codon_start: int
    codon: str
    k: int  # codons from the annotated stop to the new stop
    ext_stop: int


def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_genome(out: Path) -> Genome:
    seqs: dict[str, list[str]] = {}
    circ: dict[str, bool] = {}
    name = ""
    for ln in (out / "genome.fasta").read_text().splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            seqs[name] = []
            circ[name] = "topology=linear" not in ln
        elif ln:
            seqs[name].append(ln.strip())
    return Genome({k: "".join(v) for k, v in seqs.items()}, circ)


def tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def codon_at(g: Genome, contig: str, a: int, strand: str) -> str | None:
    """Gene-strand codon whose lowest coordinate is ``a`` (wrapped on circular replicons); None
    off a linear replicon."""
    seq, n = g.seqs[contig], g.n(contig)
    if g.circular[contig]:
        a %= n
    elif a < 0 or a + 3 > n:
        return None
    c = "".join(seq[(a + i) % n] for i in range(3))
    return c if strand == "+" else rc(c)


def terminal_codon_start(x: object, n: int) -> int:
    """0-based lowest coordinate of a ``truth/features.tsv`` row's last codon (wrapped)."""
    end, start1 = int(x.end), int(x.start)  # type: ignore[attr-defined]
    return (end - 3) % n if x.strand == "+" else (start1 - 1) % n  # type: ignore[attr-defined]


def step(strand: str, k: int) -> int:
    return 3 * k if strand == "+" else -3 * k


def eligible(out: Path) -> dict[str, Eligible]:
    g = read_genome(out)
    feats = tsv(out / "truth" / "features.tsv")
    found: dict[str, Eligible] = {}
    for x in feats.itertuples(index=False):
        if x.type not in ("cds", "sorf") or x.pseudo_reason:
            continue
        n = g.n(x.contig)
        cs = terminal_codon_start(x, n)
        stop = codon_at(g, x.contig, cs, x.strand)
        if stop not in STOPS:
            continue
        for k in range(1, K_MAX + 1):
            c = codon_at(g, x.contig, cs + step(x.strand, k), x.strand)
            if c is None:  # left a linear replicon before the next stop
                break
            if c in STOPS:
                if k >= K_MIN:
                    ext = (cs + step(x.strand, k)) % n
                    found[x.feature_id] = Eligible(
                        x.feature_id, x.type, x.contig, x.strand, cs, stop, k, ext
                    )
                break
    return found


def annotated_proteins(out: Path) -> dict[str, str]:
    doc = json.loads((out / "bakta_like.json").read_text())
    return {f["locus"]: f.get("aa") or "" for f in doc["features"] if f.get("locus")}


def extension(g: Genome, e: Eligible) -> str:
    """The ``k - 1`` residues between the annotated stop and the new stop."""
    out = []
    for j in range(1, e.k):
        c = codon_at(g, e.contig, e.codon_start + step(e.strand, j), e.strand)
        assert c is not None and c not in STOPS
        out.append(CODE_11[c])
    return "".join(out)


def expressed(out: Path, e: Eligible, g: Genome | None = None) -> str:
    g = g or read_genome(out)
    return annotated_proteins(out)[e.locus_tag] + RESIDUE[e.codon] + extension(g, e)


def readthrough_table(out: Path) -> pd.DataFrame:
    df = tsv(out / "truth" / "readthrough.tsv")
    for col in ("codon_start", "ext_stop", "n_true_peptides_spanning"):
        df[col] = df[col].astype(int)
    return df


def holds(s0: int, e0: int, codon_start: int, n: int, circular: bool) -> bool:
    """The 0-based interval ``[s0, e0)`` (``e0`` may exceed ``n``) holds the codon in frame."""
    xs = (codon_start, codon_start + n) if circular else (codon_start,)
    return any(s0 <= x and x + 3 <= e0 and (x - s0) % 3 == 0 for x in xs)


def translate_interval(
    g: Genome, contig: str, s0: int, e0: int, strand: str, stop_at: int, residue: str
) -> str:
    """Translation of ``[s0, e0)`` on ``strand`` with the codon at ``stop_at`` (lowest
    coordinate, wrapped) read as ``residue``."""
    n = g.n(contig)
    lows = list(range(s0, e0, 3)) if strand == "+" else list(range(e0 - 3, s0 - 3, -3))
    aa = []
    for a in lows:
        c = codon_at(g, contig, a, strand)
        assert c is not None
        aa.append(residue if a % n == stop_at % n else CODE_11[c])
    return "".join(aa)


def spanning_peptides(out: Path, locus: str, codon_start: int, contig: str) -> set[str]:
    g = read_genome(out)
    t = tsv(out / "truth" / "peptides.tsv")
    t = t[(t["is_wrong"] == "false") & (t["locus_tag"] == locus)]
    n, circ = g.n(contig), g.circular[contig]
    return {
        r.peptide
        for r in t.itertuples(index=False)
        if holds(int(r.start) - 1, int(r.end), codon_start, n, circ)
    }
