"""Planted selenocysteine codons (M7-16c rule 23; ``SimConfig.n_selenoproteins``).

``n`` distinct annotated, non-pseudo, non-origin-crossing CDSs with at least 21 residues each get
one codon ``i`` (``10 <= i <= len(aa) - 11``) rewritten to TGA on the gene strand; residue ``i``
of the annotated protein becomes ``U`` and the expressed protein carries ``U`` there
(``Feat.sec_residues``). The codon overlaps no other annotated feature (its span, N-terminal
extensions included) and no other simulator edit (planted novel ORFs and the RBS windows before
them and before N-terminal true starts). Every draw comes from ``default_rng([seed,
_SEC_STREAM])``; with ``n = 0`` nothing is drawn and nothing changes. Sec only, no Pyl (Q23).
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from giapsim.codons import revcomp
from giapsim.features import Feat

__all__ = ["SEC_COLUMNS", "plant_selenoproteins"]

SEC_COLUMNS = ["locus_tag", "contig", "strand", "codon_start", "residue"]
_SEC_STREAM = 12  # numpy SeedSequence key of the Sec draws (1-9 and 40 are taken)
_EDGE_BP = 20  # RBS window kept clear before planted starts (spacer <= 10 + AGGAGG + slack)
_MIN_AA, _MARGIN = 21, 10


def _cover(n: int, spans: Sequence[tuple[int, int]]) -> np.ndarray:
    """Per-base count of the (unwrapped, possibly origin-crossing) ``spans`` on a contig of
    ``n`` bases."""
    d = np.zeros(n + 1, dtype=np.int64)
    for s, e in spans:
        s, e = max(s, -n), min(e, 2 * n)
        for a, b in ((s, e), (s - n, e - n), (s + n, e + n)):
            a, b = max(a, 0), min(b, n)
            if a < b:
                d[a] += 1
                d[b] -= 1
    return np.cumsum(d[:-1])


def _codon_start(f: Feat, i: int) -> int:
    return f.start + 3 * i if f.strand == "+" else f.end - 3 * i - 3


def plant_selenoproteins(
    n: int,
    seed: int,
    seqs: dict[str, bytearray],
    feats: Sequence[Feat],
    novel: Sequence[Sequence[Any]],
) -> list[list[Any]]:
    """Plant ``n`` Sec codons in ``seqs`` and ``feats`` (in place); the truth rows sorted by
    (contig, codon_start).

    Raises:
        ValueError: if fewer than ``n`` genes are eligible (before any edit).
    """
    if n <= 0:
        return []
    spans: dict[str, list[tuple[int, int]]] = {c: [] for c in seqs}
    for f in feats:
        spans[f.contig].append(f.span())
        if f.kind == "cds" and f.true_start and not f.pseudo_reason:
            s, e = f.span()
            spans[f.contig].append((s - _EDGE_BP, s) if f.strand == "+" else (e, e + _EDGE_BP))
    for r in novel:
        spans[str(r[2])].append((int(r[3]) - 1 - _EDGE_BP, int(r[4]) + _EDGE_BP))
    cover = {c: _cover(len(s), spans[c]) for c, s in seqs.items()}

    eligible: list[tuple[Feat, list[int]]] = []
    for f in feats:
        if (
            f.kind != "cds" or f.pseudo_reason or f.extra.get("pseudo") or f.extra.get("edge")
            or f.end > len(seqs[f.contig]) or len(f.aa) < _MIN_AA
        ):  # fmt: skip
            continue
        cov = cover[f.contig]
        res = [
            i
            for i in range(_MARGIN, len(f.aa) - _MARGIN)
            if (cov[_codon_start(f, i) : _codon_start(f, i) + 3] == 1).all()
        ]
        if res:
            eligible.append((f, res))
    if len(eligible) < n:
        raise ValueError(f"n_selenoproteins={n} but only {len(eligible)} eligible genes")

    rng = np.random.default_rng([seed, _SEC_STREAM])
    rows: list[list[Any]] = []
    for k in sorted(rng.choice(len(eligible), size=n, replace=False).tolist()):
        f, res = eligible[k]
        i = res[int(rng.integers(len(res)))]
        cs = _codon_start(f, i)
        seqs[f.contig][cs : cs + 3] = ("TGA" if f.strand == "+" else revcomp("TGA")).encode()
        f.aa = f.aa[:i] + "U" + f.aa[i + 1 :]
        if len(f.nt) == f.end - f.start:
            f.nt = f.nt[: 3 * i] + "TGA" + f.nt[3 * i + 3 :]
        f.sec_residues = (*f.sec_residues, i)
        rows.append([f.locus, f.contig, f.strand, cs, i])
    return sorted(rows, key=lambda r: (r[1], r[3]))
