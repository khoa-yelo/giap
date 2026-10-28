"""Test helpers for M11-14w (docs/specs/M11-14w-wrong-genetic-code-warning.md).

- ``synthetic_genome``: the spec's fixture generator (test plan): ``n_genes`` genes of
  ``n_codons`` codons (ATG, sense codons, TAA) drawn from base frequencies at GC ``gc`` with the
  three stops rejected; for code 4 (TGA = Trp) and code 25 (TGA = Gly; M11-16 item 8) each
  sense codon becomes TGA with probability ``tga_rate``.
  Genes are separated by ``spacer`` bp random spacers (same base frequencies) and half of them
  are reverse-complemented. ``numpy.random.default_rng(seed)``.
- ``reference_runs`` / ``reference_cover``: a pure-Python, codon-by-codon brute-force reference
  of the spec's run definition (six frames, linear record, upper-cased, non-ACGT codon breaks a
  run, the stop is not covered, open runs at the record end count, minus-strand ``[a, b)`` maps
  to ``[n-b, n-a)``).
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

S11 = frozenset({"TAA", "TAG", "TGA"})
S4 = frozenset({"TAA", "TAG"})
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_BASES = np.array(list("ACGT"))
_STOP_SET = frozenset({"TAA", "TAG", "TGA"})


def revcomp(seq: str) -> str:
    """Reverse complement; any non-ACGT character stays itself (it breaks runs either way)."""
    return seq.translate(_COMP)[::-1]


def _strand_runs(seq: str, stops: Iterable[str], min_codons: int) -> list[tuple[int, int, int]]:
    """``(frame, start, end)`` runs of >= ``min_codons`` whole codons on ``seq`` itself."""
    stop = frozenset(stops)
    out: list[tuple[int, int, int]] = []
    n = len(seq)
    for frame in range(3):
        start: int | None = None
        i = frame
        while i + 3 <= n:
            codon = seq[i : i + 3]
            broken = codon in stop or any(c not in "ACGT" for c in codon)
            if broken:
                if start is not None and (i - start) // 3 >= min_codons:
                    out.append((frame, start, i))
                start = None
            elif start is None:
                start = i
            i += 3
        if start is not None and (i - start) // 3 >= min_codons:
            out.append((frame, start, i))  # open run at the record end
    return out


def reference_runs(
    seq: str, stops: Iterable[str], min_codons: int
) -> list[tuple[str, int, int, int]]:
    """``(strand, frame, start, end)`` in record coordinates for both strands."""
    s = seq.upper()
    n = len(s)
    plus = [("+", f, a, b) for f, a, b in _strand_runs(s, stops, min_codons)]
    minus = [("-", f, n - b, n - a) for f, a, b in _strand_runs(revcomp(s), stops, min_codons)]
    return plus + minus


def reference_cover(seq: str, stops: Iterable[str], min_codons: int) -> tuple[int, int]:
    """``(positions covered by at least one run, number of runs)`` for one record."""
    runs = reference_runs(seq, stops, min_codons)
    covered = np.zeros(len(seq), dtype=bool)
    for _, _, a, b in runs:
        covered[a:b] = True
    return int(covered.sum()), len(runs)


def _probs(gc: float) -> np.ndarray:
    return np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])


def _sense_codons(rng: np.random.Generator, m: int, p: np.ndarray) -> list[str]:
    codons = ["".join(r) for r in _BASES[rng.choice(4, size=(m, 3), p=p)]]
    bad = [i for i, c in enumerate(codons) if c in _STOP_SET]
    while bad:
        redraw = ["".join(r) for r in _BASES[rng.choice(4, size=(len(bad), 3), p=p)]]
        for i, c in zip(bad, redraw, strict=True):
            codons[i] = c
        bad = [i for i in bad if codons[i] in _STOP_SET]
    return codons


def synthetic_genome(
    seed: int,
    gc: float,
    code: int,
    *,
    tga_rate: float = 0.01,
    n_genes: int = 300,
    n_codons: int = 400,
    spacer: int = 60,
) -> str:
    """One record: the spec's code-4 (``code=4``), code-25 (``code=25``, M11-16) or code-11
    (``code=11``) fixture; codes 4 and 25 give the same sequence at the same arguments."""
    rng = np.random.default_rng(seed)
    p = _probs(gc)
    minus = set(rng.permutation(n_genes)[: n_genes // 2].tolist())
    parts: list[str] = []
    for i in range(n_genes):
        parts.append("".join(_BASES[rng.choice(4, size=spacer, p=p)]))
        sense = _sense_codons(rng, n_codons - 2, p)
        if code in (4, 25):
            hit = rng.random(len(sense)) < tga_rate
            sense = ["TGA" if h else c for c, h in zip(sense, hit, strict=True)]
        gene = "ATG" + "".join(sense) + "TAA"
        parts.append(revcomp(gene) if i in minus else gene)
    parts.append("".join(_BASES[rng.choice(4, size=spacer, p=p)]))
    return "".join(parts)


def random_dna(seed: int, n: int, gc: float = 0.5) -> str:
    rng = np.random.default_rng(seed)
    return "".join(_BASES[rng.choice(4, size=n, p=_probs(gc))])
