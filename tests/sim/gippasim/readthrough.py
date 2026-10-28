"""Planted stop-codon read-through (M7-16 rule 14; ADR-0058).

``SimConfig.n_readthrough = n > 0`` picks ``n`` non-pseudo CDS / sORFs whose next in-frame stop
lies ``k`` in ``[5, 60]`` codons past the annotated stop (inside a linear replicon; wrapping on a
circular one). The genome and the annotation are untouched: the expressed protein reads the
annotated stop as ``W`` (TGA) or ``Q`` (TAA, TAG) and continues for the ``k - 1`` sense codons
before the new stop. The draw uses its own stream, ``default_rng([seed, 10])``, so ``n = 0``
builds are byte-identical to the builds before the knob.

``truth/readthrough.tsv`` (written only when ``n > 0``): ``codon_start`` and ``ext_stop`` are the
annotated and the new stop's lowest 0-based coordinates (wrapped); ``n_true_peptides_spanning``
counts the distinct detected correct peptides of the gene whose interval holds the annotated
stop codon in frame.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from gippasim.codons import CODE_11, STOPS, revcomp
from gippasim.features import Feat

__all__ = [
    "K_MAX",
    "K_MIN",
    "READTHROUGH_COLUMNS",
    "RESIDUE",
    "ReadThrough",
    "eligible_readthrough",
    "pick_readthrough",
    "readthrough_rows",
]

READTHROUGH_COLUMNS = (
    "locus_tag", "contig", "strand", "codon_start", "codon", "residue", "ext_stop",
    "n_true_peptides_spanning",
)  # fmt: skip
RESIDUE: Mapping[str, str] = {"TGA": "W", "TAA": "Q", "TAG": "Q"}
K_MIN, K_MAX = 5, 60
_STREAM = 10  # numpy SeedSequence key of the read-through draw


@dataclass(frozen=True)
class ReadThrough:
    """A planted read-through gene: its terminal stop and the ``k - 1`` extension residues."""

    locus: str
    contig: str
    strand: str
    codon_start: int  # 0-based lowest coordinate of the annotated stop, wrapped
    codon: str
    k: int  # codons from the annotated stop to the new stop
    ext_stop: int  # 0-based lowest coordinate of the new stop, wrapped
    extension: str  # residues after the read stop residue, before the new stop (k - 1)

    @property
    def residue(self) -> str:
        return RESIDUE[self.codon]


def _codon(seq: str, a: int, strand: str, circular: bool) -> str | None:
    n = len(seq)
    if circular:
        a %= n
    elif a < 0 or a + 3 > n:
        return None
    c = "".join(seq[(a + i) % n] for i in range(3))
    return c if strand == "+" else revcomp(c)


def eligible_readthrough(
    feats: Sequence[Feat], seqs: Mapping[str, str], circular: Mapping[str, bool]
) -> list[ReadThrough]:
    """Every eligible gene, in feature order."""
    out = []
    for f in feats:
        if f.kind not in ("cds", "sorf") or f.pseudo_reason or f.extra.get("pseudo"):
            continue
        seq, circ = seqs[f.contig], circular[f.contig]
        n = len(seq)
        cs = (f.end - 3) % n if f.strand == "+" else f.start % n
        stop = _codon(seq, cs, f.strand, circ)
        if stop not in STOPS:
            continue
        sign = 3 if f.strand == "+" else -3
        ext: list[str] = []
        for k in range(1, K_MAX + 1):
            c = _codon(seq, cs + sign * k, f.strand, circ)
            if c is None:
                break
            if c in STOPS:
                if k >= K_MIN:
                    out.append(ReadThrough(f.locus, f.contig, f.strand, cs, stop, k,
                                           (cs + sign * k) % n, "".join(ext)))  # fmt: skip
                break
            ext.append(CODE_11[c])
    return out


def pick_readthrough(
    n: int,
    seed: int,
    feats: Sequence[Feat],
    seqs: Mapping[str, str],
    circular: Mapping[str, bool],
) -> list[ReadThrough]:
    """``n`` eligible genes drawn without replacement, in feature order; ``ValueError`` when
    fewer than ``n`` are eligible (no silent shortfall)."""
    if n <= 0:
        return []
    elig = eligible_readthrough(feats, seqs, circular)
    if len(elig) < n:
        raise ValueError(f"n_readthrough = {n} but only {len(elig)} eligible genes")
    rng = np.random.default_rng([seed, _STREAM])
    idx = sorted(rng.choice(len(elig), size=n, replace=False).tolist())
    return [elig[i] for i in idx]


def _holds(s0: int, e0: int, cs: int, n: int, circular: bool) -> bool:
    xs = (cs, cs + n) if circular else (cs,)
    return any(s0 <= x and x + 3 <= e0 and (x - s0) % 3 == 0 for x in xs)


def readthrough_rows(
    planted: Sequence[ReadThrough],
    truth: Sequence[Mapping[str, str]],
    lengths: Mapping[str, int],
    circular: Mapping[str, bool],
) -> list[list[object]]:
    """``truth/readthrough.tsv`` rows from the ``truth/peptides.tsv`` records."""
    rows: list[list[object]] = []
    for p in planted:
        n, circ = lengths[p.contig], circular[p.contig]
        peps = {
            r["peptide"]
            for r in truth
            if r["is_wrong"] == "false" and r["locus_tag"] == p.locus
            and _holds(int(r["start"]) - 1, int(r["end"]), p.codon_start, n, circ)
        }  # fmt: skip
        rows.append([p.locus, p.contig, p.strand, p.codon_start, p.codon, p.residue, p.ext_stop,
                     len(peps)])  # fmt: skip
    return rows
