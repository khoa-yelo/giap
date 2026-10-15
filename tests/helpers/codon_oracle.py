"""Independent codon-usage oracle for the M2-17 tests (docs/specs/M2-17-codon-usage.md).

Pure Python (no codon-bias, no numpy), written from the definitions:

* translation table 11 (sense codons identical to table 1; stops TAA TAG TGA);
* scored sequence: codons 2 .. n-1 of a CDS (start and stop dropped, lead decision Q1); ``None``
  for pseudogenes, partial genes, ``len % 3 != 0`` or fewer than 10 codons in the full CDS;
* CAI (Sharp & Li 1987) with codon-bias' conventions: reference counts + 1 pseudocount per sense
  codon, weight = count / max count in the synonymous family (Met, Trp weight 1), score = geometric
  mean of the weights over the sense codons of the sequence;
* Fop (Ikemura 1981) with the same weights: optimal when weight >= 0.8; share of sense codons;
* tAI (dos Reis 2004, prokaryote s-values, codon-bias' conventions): W(codon) = sum over tRNAs
  pairing at positions 1-2 of (1 - s) * copies; w = W / max W; zero weights get the geometric mean
  of the non-zero ones; codons no tRNA can read have no weight and are skipped; score = geometric
  mean over the weighted sense codons;
* GC3, ``at_2_10``, synonymous shares, rare codons (< 10 % strict) as in the spec.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Mapping

_B = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE11: dict[str, str] = {
    a + b + c: _AA11[16 * i + 4 * j + k]
    for i, a in enumerate(_B)
    for j, b in enumerate(_B)
    for k, c in enumerate(_B)
}
SENSE: tuple[str, ...] = tuple(sorted(c for c, aa in CODE11.items() if aa != "*"))
FAMILY: dict[str, tuple[str, ...]] = {}
for _c in SENSE:
    FAMILY.setdefault(CODE11[_c], ())
    FAMILY[CODE11[_c]] += (_c,)
_COMP = str.maketrans("ACGT", "TGCA")

# dos Reis et al. 2004 s-values: (anticodon first base, codon third base) -> (s, min degeneracy)
_SVAL: dict[tuple[str, str], tuple[float, int]] = {
    ("A", "T"): (0.0, 1), ("G", "C"): (0.0, 1), ("T", "A"): (0.0, 1), ("C", "G"): (0.0, 1),
    ("G", "T"): (0.41, 2), ("A", "C"): (0.28, 3), ("A", "A"): (0.9999, 3), ("T", "G"): (0.68, 2),
    ("C", "A"): (0.89, 2),  # prokaryote only (lysidine-modified CAT reads ATA)
}  # fmt: skip


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def codons(seq: str) -> list[str]:
    """In-frame codons (a trailing partial codon is dropped), upper case."""
    s = seq.upper()
    return [s[i : i + 3] for i in range(0, len(s) - len(s) % 3, 3)]


def sense_codons(seq: str) -> list[str]:
    return [c for c in codons(seq) if c in CODE11 and CODE11[c] != "*"]


def scored(nt: str | None, feature_type: str = "cds", partial: bool = False) -> str | None:
    if nt is None or feature_type == "pseudogene" or partial:
        return None
    if len(nt) % 3 or len(nt) < 30:
        return None
    return nt[3:-3]


def at_2_10(nt: str) -> float | None:
    if len(nt) < 30:
        return None
    w = nt[3:30].upper()
    acgt = sum(w.count(b) for b in "ACGT")
    return None if acgt == 0 else (w.count("A") + w.count("T")) / acgt


def gc3(seq: str) -> float:
    cs = codons(seq)
    thirds = [c[2] for c in cs if c[2] in "ACGT"]
    return sum(b in "GC" for b in thirds) / len(thirds)


def shares(seqs: Iterable[str]) -> dict[str, float]:
    n: Counter[str] = Counter()
    for s in seqs:
        n.update(sense_codons(s))
    out: dict[str, float] = {}
    for fam in FAMILY.values():
        tot = sum(n[c] for c in fam)
        if tot:
            out.update({c: n[c] / tot for c in fam})
    return out


def rare(sh: Mapping[str, float], cutoff: float = 0.10) -> frozenset[str]:
    return frozenset(c for c, v in sh.items() if v < cutoff)


def rare_2_10(nt: str, rare_set: frozenset[str]) -> int | None:
    if len(nt) < 30:
        return None
    return sum(c in rare_set for c in codons(nt[3:30]))


def _family_weights(n: Mapping[str, float]) -> dict[str, float]:
    """Frequencies within the family first, then relative to the family maximum (codon-bias'
    order of operations, which matters for Fop's ``>= 0.8`` at the boundary, e.g. 4/9 vs 5/9)."""
    w: dict[str, float] = {}
    for fam in FAMILY.values():
        tot = sum(n[c] for c in fam)
        freq = {c: n[c] / tot for c in fam}
        top = max(freq.values())
        w.update({c: freq[c] / top for c in fam})
    return w


def cai_weights(ref: Iterable[str], pseudocount: int = 1) -> dict[str, float]:
    n: Counter[str] = Counter()
    for s in ref:
        n.update(sense_codons(s))
    return _family_weights({c: n[c] + pseudocount for c in SENSE})


def geomean(weights: Mapping[str, float], seq: str) -> float | None:
    """Geometric mean over the weighted sense codons; ``None`` when there is none (NaN in
    codon-bias, reported as missing)."""
    cs = [c for c in sense_codons(seq) if c in weights]
    if not cs:
        return None
    return math.exp(sum(math.log(weights[c]) for c in cs) / len(cs))


def cai(seq: str, ref: Iterable[str]) -> float | None:
    return geomean(cai_weights(ref), seq)


def fop(seq: str, ref: Iterable[str], thresh: float = 0.8) -> float:
    w = cai_weights(ref)
    cs = sense_codons(seq)
    return sum(w[c] >= thresh for c in cs) / len(cs)


def tai_weights(trna_counts: Mapping[str, int]) -> dict[str, float]:
    raw: dict[str, float] = {}
    for c in SENSE:
        deg = len(FAMILY[CODE11[c]])
        pairs = []
        for anti, gcn in trna_counts.items():
            a = anti.upper().replace("U", "T")
            sv = _SVAL.get((a[0], c[2]))
            if revcomp(a)[:2] == c[:2] and sv is not None and deg >= sv[1]:
                pairs.append((1 - sv[0]) * gcn)
        if pairs:  # codons no tRNA can read have no weight and are left out of the score
            raw[c] = sum(pairs)
    top = max(raw.values())
    w = {c: v / top for c, v in raw.items()}
    nz = [v for v in w.values() if v > 0]
    g = math.exp(sum(math.log(v) for v in nz) / len(nz))
    return {c: (v if v > 0 else g) for c, v in w.items()}


def tai(seq: str, trna_counts: Mapping[str, int]) -> float | None:
    return geomean(tai_weights(trna_counts), seq)
