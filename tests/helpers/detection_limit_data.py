"""Test data and oracles for the M6-02 detection-limit model
(docs/specs/M6-02-detection-limit-model.md).

Everything is built by hand so the oracles do not reuse the code under test:

- peptides are tryptic by construction: interior residues never K or R, the last residue K, so a
  ``trypsin`` / ``trypsin_p`` digest of the concatenation returns exactly the planted peptides;
  I and L are never used (no I/L-equivalence question in matching);
- GRAVY is the Kyte-Doolittle mean computed here from a literal table (the scale M4-10 uses);
- the Poisson-binomial oracle enumerates all 2^n outcomes (n <= 12);
- C11 ``peptide_hits`` rows carry every contract column (``giap.genetable.contract.TABLES``).
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# Kyte & Doolittle (1982) hydropathy, the scale of Biopython ProtParam ``gravy`` (M4-10).
KD: Mapping[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}  # fmt: skip
HYDROPHOBIC = "AVFMC"
HYDROPHILIC = "DENQSTGHYW"

#: The M6-02 toy model: logit p_j = -12 + 0.6 a - 0.02 L_j (GRAVY coefficient 0).
TOY_COEF: Mapping[str, float] = {
    "intercept": -12.0,
    "abundance_log2": 0.6,
    "length": -0.02,
    "gravy": 0.0,
}
#: Gene u of the examples: three observable peptides of length 10, 20, 30.
TOY_PEPS: tuple[tuple[int, float], ...] = ((10, 0.0), (20, 0.0), (30, 0.0))

STEP = "operons.detection_limit"


def gravy(pep: str) -> float:
    return sum(KD[c] for c in pep) / len(pep)


def expit(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def peptide_p(coef: Mapping[str, float], a: float, length: int, g: float) -> float:
    """Per-peptide probability of the linear peptide model (hand oracle)."""
    eta = (
        coef["intercept"]
        + coef["abundance_log2"] * a
        + coef["length"] * length
        + coef.get("gravy", 0.0) * g
    )
    return expit(eta)


def pb_sf_enum(p: Sequence[float], k: int) -> float:
    """P(X >= k) by enumerating every outcome (closed form; n <= 12)."""
    n = len(p)
    if n > 12:
        raise ValueError("enumeration oracle is for n <= 12")
    total = 0.0
    for s in itertools.product((0, 1), repeat=n):
        if sum(s) >= k:
            total += math.prod(q if b else 1.0 - q for q, b in zip(p, s, strict=True))
    return total


def p_detect_oracle(
    coef: Mapping[str, float], a: float, peps: Sequence[tuple[int, float]], k_min: int
) -> float:
    if not peps:
        return math.nan
    return pb_sf_enum([peptide_p(coef, a, n, g) for n, g in peps], k_min)


def p_miss_oracle(
    coef: Mapping[str, float],
    a_v: float,
    ratios: Sequence[float],
    peps: Sequence[tuple[int, float]],
    k_min: int,
) -> float:
    if len(ratios) == 0 or not peps:
        return math.nan
    return float(np.mean([1.0 - p_detect_oracle(coef, a_v + r, peps, k_min) for r in ratios]))


# ------------------------------------------------------------------------------ sequences


def make_peptide(rng: np.random.Generator, length: int, frac_hydrophobic: float) -> str:
    """A tryptic peptide of ``length`` (>= 2): ``length - 1`` residues without K/R/I/L, then K."""
    m = length - 1
    hyd = rng.random(m) < frac_hydrophobic
    a = rng.integers(0, len(HYDROPHOBIC), m)
    b = rng.integers(0, len(HYDROPHILIC), m)
    body = [HYDROPHOBIC[i] if h else HYDROPHILIC[j] for h, i, j in zip(hyd, a, b, strict=True)]
    return "".join(body) + "K"


def make_peptides(
    rng: np.random.Generator, n: int, lengths: tuple[int, int] = (7, 30)
) -> list[str]:
    """``n`` distinct tryptic peptides with lengths uniform in the inclusive range."""
    out: list[str] = []
    seen: set[str] = set()
    while len(out) < n:
        pep = make_peptide(rng, int(rng.integers(lengths[0], lengths[1] + 1)), float(rng.random()))
        if pep not in seen:
            seen.add(pep)
            out.append(pep)
    return out


def protein(peptides: Iterable[str], spacers: bool = False) -> str:
    """Concatenation of tryptic peptides; ``spacers`` interleaves unobservable ``GK``."""
    parts: list[str] = []
    for pep in peptides:
        parts.append(pep)
        if spacers:
            parts.append("GK")  # length 2: never observable
    return "".join(parts)


# ------------------------------------------------------------------------------ C11 rows


def hit_row(
    peptide: str, gene: str | None, uniqueness: str = "unique", route: str = "exact"
) -> dict[str, object]:
    """One C11 ``peptide_hits`` row (internal representation) with every contract column."""
    assigning = gene is not None
    return {
        "peptide": peptide,
        "peptide_len": len(peptide),
        "source": "database",
        "best_q_input": 0.001,
        "pep_input": 0.01,
        "score_input": 10.0,
        "contig": "chr1",
        "start": 0,
        "end": 3 * len(peptide),
        "strand": "+",
        "frame": 0,
        "crosses_origin": False,
        "route": route,
        "n_diff": 0 if route == "exact" else 1,
        "diff_types": [] if route == "exact" else ["substitution"],
        "placement_q": 0.001,
        "length_class": None,
        "locus_class": "cds_inframe" if assigning else "intergenic",
        "uniqueness": uniqueness,
        "gene": gene,
        "novel_candidate_id": None,
        "flags": [],
    }


HIT_COLUMNS: tuple[str, ...] = tuple(hit_row("AAAAAAK", "g").keys())


def hits_frame(rows: Sequence[Mapping[str, object]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame({c: pd.Series([], dtype=object) for c in HIT_COLUMNS})
    return pd.DataFrame(list(rows), columns=list(HIT_COLUMNS))


def runs_frame(rows: Sequence[tuple[str, str]]) -> pd.DataFrame:
    """C16 ``peptide_runs`` (M1-06 columns) from (peptide, run) presence pairs."""
    return pd.DataFrame(
        {
            "peptide": pd.Series([p for p, _ in rows], dtype=object),
            "run": pd.Series([r for _, r in rows], dtype=object),
            "intensity": pd.Series([1000.0] * len(rows), dtype="float64"),
            "n_rows": pd.Series([1] * len(rows), dtype="int64"),
        }
    )


# ------------------------------------------------------------------------------ genes


def genes_frame(
    tags: Sequence[str],
    tus: Sequence[str],
    aas: Sequence[str],
    detected: Sequence[bool],
    abundance: Sequence[float],
    tier: str = "mayu",
    extra: Mapping[str, Sequence[object]] | None = None,
) -> pd.DataFrame:
    """The gene frame M6-02 reads: ``locus_tag``, ``tu_id``, ``aa``, ``detected``,
    ``abundance_log2``, ``detection_tier`` (plus optional covariate columns)."""
    df = pd.DataFrame(
        {
            "locus_tag": pd.Series(list(tags), dtype=object),
            "tu_id": pd.Series(list(tus), dtype=object),
            "aa": pd.Series(list(aas), dtype=object),
            "detected": pd.Series(list(detected), dtype=bool),
            "abundance_log2": pd.Series(list(abundance), dtype="float64"),
            "detection_tier": pd.Series([tier] * len(tags), dtype=object),
        }
    )
    for k, v in (extra or {}).items():
        df[k] = list(v)
    return df


@dataclass
class PeptideSim:
    """A simulated genome for the peptide model with known per-peptide truth."""

    genes: pd.DataFrame
    hits: pd.DataFrame
    peptides: dict[str, list[str]]  # locus_tag -> its observable peptides
    observed: dict[str, np.ndarray]  # locus_tag -> bool per peptide
    beta: Mapping[str, float]
    truth_p: dict[str, np.ndarray] = field(default_factory=dict)

    def peps(self, tag: str) -> list[tuple[int, float]]:
        return [(len(p), gravy(p)) for p in self.peptides[tag]]


def simulate_peptides(
    rng: np.random.Generator,
    n_genes: int,
    beta: Mapping[str, float],
    a_mu: float = 15.0,
    a_sd: float = 2.5,
    n_peps: tuple[int, int] = (5, 20),
    tu_size: int = 2,
    k_min: int = 1,
    tier: str = "mayu",
    min_total: int = 0,
) -> PeptideSim:
    """Genes with tryptic peptides, logistic peptide truth in (a, L, GRAVY), detected = >= k_min
    observed peptides. ``min_total`` keeps adding genes until the detected genes carry at least
    that many observable peptides (example 7: 5,000 training rows)."""
    tags: list[str] = []
    tus: list[str] = []
    aas: list[str] = []
    det: list[bool] = []
    ab: list[float] = []
    peptides: dict[str, list[str]] = {}
    observed: dict[str, np.ndarray] = {}
    truth: dict[str, np.ndarray] = {}
    rows: list[dict[str, object]] = []
    n_train = 0
    g = 0
    while g < n_genes or n_train < min_total:
        tag = f"G{g:05d}"
        a = float(rng.normal(a_mu, a_sd))
        peps = make_peptides(rng, int(rng.integers(n_peps[0], n_peps[1] + 1)))
        p = np.array([peptide_p(beta, a, len(s), gravy(s)) for s in peps])
        o = rng.random(len(peps)) < p
        d = int(o.sum()) >= k_min
        tags.append(tag)
        tus.append(f"TU{g // tu_size:05d}")
        aas.append(protein(peps))
        det.append(d)
        ab.append(a)
        peptides[tag] = peps
        observed[tag] = o
        truth[tag] = p
        rows.extend(hit_row(s, tag) for s, seen in zip(peps, o, strict=True) if seen)
        if d:
            n_train += len(peps)
        g += 1
    genes = genes_frame(tags, tus, aas, det, ab, tier=tier)
    return PeptideSim(genes, hits_frame(rows), peptides, observed, dict(beta), truth)


def auc(y: Sequence[bool], s: Sequence[float]) -> float:
    """Mann-Whitney AUC (ties count 1/2)."""
    yy = np.asarray(y, dtype=bool)
    ss = np.asarray(s, dtype=float)
    pos, neg = ss[yy], ss[~yy]
    if pos.size == 0 or neg.size == 0:
        return math.nan
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (pos.size * neg.size))
