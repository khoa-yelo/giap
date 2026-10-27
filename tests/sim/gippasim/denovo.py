"""Simulated de novo reads, their errors and reference peptides (issue M3-20).

GIPPA_SPEC.md §11.3 "De novo errors", §5.3 (reference peptides) and §8.4. Every non-decoy row of
M3-19's ``peptides.tsv`` gets one de novo read of the same spectrum. A spectrum's true peptide is
its database peptide when that identification is correct, else a random correct one. Scores are
standard logistic draws; their pooled deciles (edges at the k/10 quantiles, right-closed) set the
probability ``precision_by_decile[b]`` that the read is the true peptide (I/L re-drawn), else it
carries one error drawn with ``error_weights``. All draws come from
``numpy.random.default_rng([seed, 4])``, so M3-19 outputs are unchanged.

The block lists are gippasim's own copy of §7.2.4 (independent of the code under test).
"""

from __future__ import annotations

import csv
import itertools
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType

import numpy as np

__all__ = [
    "BLOCKS",
    "DEFAULT_ERROR_WEIGHTS",
    "ERROR_TYPES",
    "NEAR_BLOCKS",
    "DenovoConfig",
    "apply_error",
    "exchange_block",
    "misread_modification",
    "redraw_il",
    "reverse_stretch",
    "simulate_denovo",
    "swap_adjacent",
]

AA20 = "ACDEFGHIKLMNPQRSTVWY"
_STREAM = 4  # numpy SeedSequence key (1: features, 2: signals, 3: proteomics)
_LEN = (7, 40)
_RANDOM_LEN = (7, 20)
_NA = "NA"

#: GIPPA_SPEC.md §7.2.4: mass-equivalent compositions (either direction) and near-equal ones.
BLOCKS: tuple[tuple[str, str], ...] = (
    ("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"),
)  # fmt: skip
NEAR_BLOCKS: tuple[tuple[str, str], ...] = (("W", "AD"), ("W", "EG"), ("R", "GV"))
#: Modification misreads: A + Oxidation = S, N + Deamidated = D, G + Methyl = A.
MOD_READS: Mapping[str, str] = MappingProxyType({"A": "S", "N": "D", "G": "A"})

ERROR_TYPES: tuple[str, ...] = (
    "substitution", "swap", "block", "near_block", "reversed_stretch", "modification", "random",
)  # fmt: skip
DEFAULT_ERROR_WEIGHTS: Mapping[str, float] = MappingProxyType({
    "substitution": 0.3, "block": 0.2, "swap": 0.15, "near_block": 0.05,
    "reversed_stretch": 0.1, "modification": 0.1, "random": 0.1,
})  # fmt: skip
DEFAULT_PRECISION: tuple[float, ...] = tuple(
    float(x) for x in np.round(np.linspace(0.2, 0.95, 10), 6)
)


@dataclass(frozen=True)
class DenovoConfig:
    """De novo knobs (docs/specs/M3-20-sim-denovo-errors.md)."""

    precision_by_decile: tuple[float, ...] = DEFAULT_PRECISION  # bin 0 = lowest scores
    error_weights: Mapping[str, float] = field(default_factory=lambda: DEFAULT_ERROR_WEIGHTS)

    def __post_init__(self) -> None:
        p = tuple(float(x) for x in self.precision_by_decile)
        if len(p) != 10 or not all(0.0 <= x <= 1.0 for x in p):
            raise ValueError("precision_by_decile must be 10 values in [0, 1]")
        w = {k: float(v) for k, v in dict(self.error_weights).items()}
        if set(w) - set(ERROR_TYPES):
            raise ValueError(f"error_weights keys must be among {ERROR_TYPES}")
        if any(v < 0 for v in w.values()) or sum(w.values()) <= 0:
            raise ValueError("error_weights must be >= 0 with a positive sum")
        object.__setattr__(self, "precision_by_decile", p)
        object.__setattr__(self, "error_weights", MappingProxyType(w))


# ---------------------------------------------------------------- pure error operators


def _il(s: str) -> str:
    return s.replace("I", "L")


def swap_adjacent(pep: str, i: int) -> str:
    """``pep`` with residues ``i`` and ``i + 1`` exchanged."""
    if not 0 <= i < len(pep) - 1:
        raise ValueError(f"no adjacent pair at {i}")
    return pep[:i] + pep[i + 1] + pep[i] + pep[i + 2 :]


def reverse_stretch(pep: str, start: int, length: int) -> str:
    """``pep`` with ``pep[start:start + length]`` reversed (3-6 residues)."""
    if not 3 <= length <= 6 or not 0 <= start <= len(pep) - length:
        raise ValueError(f"bad stretch {start}+{length}")
    return pep[:start] + pep[start : start + length][::-1] + pep[start + length :]


def exchange_block(pep: str, start: int, old: str, new: str) -> str:
    """``pep`` with the block ``old`` at ``start`` replaced by ``new``."""
    if pep[start : start + len(old)] != old:
        raise ValueError(f"no {old!r} at {start}")
    return pep[:start] + new + pep[start + len(old) :]


def misread_modification(pep: str, i: int) -> str:
    """``pep`` with the modified residue ``i`` read as its mass-equivalent residue."""
    if pep[i] not in MOD_READS:
        raise ValueError(f"no modelled misread of {pep[i]!r}")
    return pep[:i] + MOD_READS[pep[i]] + pep[i + 1 :]


def redraw_il(pep: str, rng: np.random.Generator) -> str:
    """``pep`` with every I/L re-drawn as I or L at random."""
    pos = [i for i, c in enumerate(pep) if c in "IL"]
    if not pos:
        return pep
    out = list(pep)
    for i, u in zip(pos, rng.random(len(pos)).tolist(), strict=True):
        out[i] = "I" if u < 0.5 else "L"
    return "".join(out)


def _orders(s: str) -> list[str]:
    """Every residue order of ``s`` with each L also written as I."""
    out = {"".join(p) for p in itertools.permutations(s)}
    return sorted({x for o in out for x in _il_spellings(o)})


def _il_spellings(s: str) -> list[str]:
    opts = [("I", "L") if c == "L" else (c,) for c in s]
    return ["".join(p) for p in itertools.product(*opts)]


def _block_sites(pep: str, table: tuple[tuple[str, str], ...]) -> list[tuple[int, str, str]]:
    sites = []
    for a, b in table:
        for old, new in ((a, b), (b, a)):
            if not _LEN[0] <= len(pep) - len(old) + len(new) <= _LEN[1]:
                continue
            for o in _orders(old):
                sites += [(i, o, new) for i in range(len(pep) - len(o) + 1) if pep.startswith(o, i)]
    return sites


def apply_error(pep: str, error_type: str, rng: np.random.Generator) -> str:
    """``pep`` after one random error of ``error_type``; ``ValueError`` if it has no site."""
    n = len(pep)
    if error_type == "substitution":
        i = int(rng.integers(n))
        alt = [c for c in AA20 if _il(c) != _il(pep[i])]
        return pep[:i] + alt[int(rng.integers(len(alt)))] + pep[i + 1 :]
    if error_type == "random":
        while True:
            k = int(rng.integers(_RANDOM_LEN[0], _RANDOM_LEN[1] + 1))
            out = "".join(AA20[j] for j in rng.integers(20, size=k).tolist())
            if _il(out) != _il(pep):
                return out
    if error_type in ("block", "near_block"):
        sites = _block_sites(pep, BLOCKS if error_type == "block" else NEAR_BLOCKS)
        if not sites:
            raise ValueError(f"no {error_type} site")
        i, old, new = sites[int(rng.integers(len(sites)))]
        perm = sorted({"".join(p) for p in itertools.permutations(new)})
        return exchange_block(pep, i, old, perm[int(rng.integers(len(perm)))])
    if error_type == "swap":
        idx = [i for i in range(n - 1) if _il(pep[i]) != _il(pep[i + 1])]
        if not idx:
            raise ValueError("no swap site")
        return swap_adjacent(pep, idx[int(rng.integers(len(idx)))])
    if error_type == "reversed_stretch":
        st = [
            (i, k)
            for k in range(3, 7)
            for i in range(n - k + 1)
            if _il(pep[i : i + k]) != _il(pep[i : i + k][::-1])
        ]
        if not st:
            raise ValueError("no reversible stretch")
        return reverse_stretch(pep, *st[int(rng.integers(len(st)))])
    if error_type == "modification":
        idx = [i for i, c in enumerate(pep) if c in MOD_READS]
        if not idx:
            raise ValueError("no modifiable residue")
        return misread_modification(pep, idx[int(rng.integers(len(idx)))])
    raise ValueError(f"unknown error type {error_type!r}")


def _draw_error(pep: str, w: np.ndarray, rng: np.random.Generator) -> tuple[str, str]:
    """One error by weight; types without a site in ``pep`` are dropped and the draw repeated."""
    w = w.copy()
    while w.sum() > 0:
        et = ERROR_TYPES[int(rng.choice(len(w), p=w / w.sum()))]
        try:
            return et, apply_error(pep, et, rng)
        except ValueError:
            w[ERROR_TYPES.index(et)] = 0.0
    return "substitution", apply_error(pep, "substitution", rng)


# ---------------------------------------------------------------- simulation


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.write_text("\n".join("\t".join(r) for r in [header, *rows]) + "\n")


def decile_bins(scores: np.ndarray) -> np.ndarray:
    """Pooled score deciles: edges at the k/10 quantiles, right-closed (as M3-11)."""
    if scores.size == 0:
        return np.zeros(0, dtype=np.int64)
    edges = np.unique(np.quantile(scores, np.arange(1, 10) / 10))
    out: np.ndarray = np.searchsorted(edges, scores, side="left")
    return out


def simulate_denovo(cfg: DenovoConfig, seed: int, out: Path) -> dict[str, Path]:
    """Write ``peptides_denovo.tsv``, ``reference_peptides.tsv``, ``truth/denovo.tsv`` and
    ``truth/denovo_precision.tsv`` from M3-19's ``peptides.tsv`` and its truth in ``out``."""
    rng = np.random.default_rng([seed, _STREAM])
    db = _read(out / "peptides.tsv")
    header = list(db[0]) if db else []
    targets = [r for r in db if r["is_decoy"] != "true"]
    wrong: dict[str, bool] = {}
    for t in _read(out / "truth" / "peptides.tsv"):
        if t["origin"] != "decoy":
            wrong[t["peptide"]] = wrong.get(t["peptide"], False) or t["is_wrong"] == "true"
    correct_peps = sorted(p for p, w in wrong.items() if not w)
    n = len(targets)
    scores = np.round(rng.logistic(size=n), 6)
    bins = decile_bins(scores)
    p_ok = np.array(cfg.precision_by_decile)[bins] if n else np.zeros(0)
    ok = rng.random(n) < p_ok
    w = np.array([cfg.error_weights.get(t, 0.0) for t in ERROR_TYPES])
    dn_rows, truth = [], []
    for r, s, good in zip(targets, scores.tolist(), ok.tolist(), strict=True):
        true = r["peptide"]
        if wrong.get(true, False) and correct_peps:
            true = correct_peps[int(rng.integers(len(correct_peps)))]
        et, pep = ("none", redraw_il(true, rng)) if good else _draw_error(true, w, rng)
        row = {**r, "peptide": pep, "protein_ids": _NA, "q_value": _NA, "pep": _NA}
        row.update(score=f"{s:.10g}", source="denovo", search_db=_NA)
        dn_rows.append([row[h] for h in header])
        truth.append([r["run"], r["spectrum_id"], pep, true, et, "true" if good else "false"])
    ref = [[{**r, "source": "database", "search_db": "same"}[h] for h in header] for r in targets]
    counts = np.bincount(bins, minlength=10) if n else np.zeros(10, dtype=np.int64)
    prec = [
        [str(b), f"{cfg.precision_by_decile[b]:.6g}", str(int(c))]
        for b, c in enumerate(counts.tolist())
    ]
    paths = {
        "denovo_peptides": out / "peptides_denovo.tsv",
        "reference_peptides": out / "reference_peptides.tsv",
        "denovo": out / "truth" / "denovo.tsv",
        "denovo_precision": out / "truth" / "denovo_precision.tsv",
    }
    _write(paths["denovo_peptides"], header, dn_rows)
    _write(paths["reference_peptides"], header, ref)
    cols = ["run", "spectrum_id", "peptide", "true_peptide", "error_type", "correct"]
    _write(paths["denovo"], cols, truth)
    _write(paths["denovo_precision"], ["score_bin", "precision_true", "n"], prec)
    return paths
