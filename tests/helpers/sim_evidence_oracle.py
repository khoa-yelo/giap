"""Independent oracle for gippasim evidence (docs/specs/M4-09-sim-evidence.md; test-author code).

Written from the spec text, not from ``gippasim.evidence``:

* Trypsin/P digestion (cut after every K/R, also before P), 7-30 aa, at most one missed cleavage
  for true and decoy source peptides; ``n_observable`` is the M4-10 count (0 missed, distinct);
* decoys: the reversed peptide with its C-terminal K/R kept in place;
* PEP = pi0 f0 / (pi0 f0 + (1 - pi0) f1), f0 = N(mu0, 1), f1 = N(mu1, 1), pi0 = ``pi_wrong``;
  q(s) = mean PEP over the target rows with score >= s;
* the exact mean and variance of PEP under each mixture component (numerical integration), for
  the calibration z-tests of the PEP sum;
* genomic coordinates of an amino-acid span of a layout gene (gene ``nt`` is 5'->3' on the gene
  strand; codon ``i`` codes residue ``i``): + strand ``[start + 3a, start + 3b)``, - strand
  ``[end - 3b, end - 3a)`` (0-based half-open);
* locating the M3 tables (``<work>/<name>`` or ``<work>/mapping/<name>``) for the M4-06 loader.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from covariates_oracle import n_observable
from peptide_oracle import trypsin_p_peptides

M3_NAMES = {"peptide_hits": "peptide_hits.tsv", "protein_hits": "protein_hits.tsv",
            "unmapped": "unmapped.tsv"}  # fmt: skip
LOCUS_CLASS = {"cds": "cds_inframe", "sorf": "sorf", "pseudogene": "pseudogene"}


# ----------------------------------------------------------------------------- sequences


def is_na(v: object) -> bool:
    if v is None:
        return True
    if isinstance(v, list | tuple | np.ndarray):
        return False
    try:
        return bool(pd.isna(v))
    except (TypeError, ValueError):
        return False


def as_list(v: object) -> list[Any]:
    if v is None or (not isinstance(v, list | tuple | np.ndarray) and is_na(v)):
        return []
    if isinstance(v, str):
        return [v]
    return list(v)  # type: ignore[call-overload]


def n_obs(aa: object) -> int:
    """M4-10 observable peptides: distinct Trypsin/P, 0 missed cleavages, 7-30 aa."""
    return n_observable(aa, "trypsin_p", (7, 30)) if isinstance(aa, str) and aa else 0


def peptides_of(aa: object, missed: int = 1) -> set[str]:
    """Trypsin/P products of 7-30 aa with at most ``missed`` missed cleavages."""
    return trypsin_p_peptides(aa, (7, 30), missed) if isinstance(aa, str) and aa else set()


def digestible(aa: object) -> bool:
    return bool(peptides_of(aa))


def reverse_keep_cterm(p: str) -> str:
    """Decoy of a tryptic peptide: reversed, with a C-terminal K/R kept at the C terminus."""
    return p[:-1][::-1] + p[-1] if p and p[-1] in "KR" else p[::-1]


def decoys_of(aa: object) -> set[str]:
    return {reverse_keep_cterm(p) for p in peptides_of(aa)}


def offsets(aa: str, peptide: str, max_diff: int = 0) -> list[int]:
    """Start offsets where ``peptide`` matches ``aa`` with at most ``max_diff`` substitutions."""
    n, k = len(aa), len(peptide)
    out = []
    for a in range(0, n - k + 1):
        if max_diff == 0:
            if aa.startswith(peptide, a):
                out.append(a)
        elif sum(x != y for x, y in zip(aa[a : a + k], peptide, strict=True)) <= max_diff:
            out.append(a)
    return out


def span_coords(start: int, end: int, strand: str, a: int, b: int) -> tuple[int, int]:
    """0-based half-open genomic span of residues ``[a, b)`` of a gene."""
    if strand == "+":
        return start + 3 * a, start + 3 * b
    return end - 3 * b, end - 3 * a


# ----------------------------------------------------------------------------- scores


def _log_phi(x: np.ndarray) -> np.ndarray:
    out: np.ndarray = -0.5 * x * x - 0.5 * math.log(2 * math.pi)
    return out


def pep(scores: Any, pi: float, mu0: float, mu1: float) -> np.ndarray:
    s = np.asarray(scores, dtype=float)
    if pi <= 0:
        return np.zeros_like(s)
    if pi >= 1:
        return np.ones_like(s)
    l0 = math.log(pi) + _log_phi(s - mu0)
    l1 = math.log1p(-pi) + _log_phi(s - mu1)
    out: np.ndarray = np.exp(l0 - np.logaddexp(l0, l1))
    return out


def q_from_targets(scores: Any, target_scores: Any, target_peps: Any) -> np.ndarray:
    """q(s) = mean PEP of the targets with score >= s (NaN where no target reaches s)."""
    s = np.asarray(scores, dtype=float)
    ts = np.asarray(target_scores, dtype=float)
    tp = np.asarray(target_peps, dtype=float)
    cum = np.cumsum(tp[np.argsort(-ts, kind="stable")])
    # number of targets with score >= s: searchsorted on the ascending copy
    asc = np.sort(ts)
    k = len(ts) - np.searchsorted(asc, s, side="left")
    out = np.full(s.shape, np.nan)
    ok = k > 0
    out[ok] = cum[k[ok] - 1] / k[ok]
    return out


def pep_moments(pi: float, mu0: float, mu1: float) -> tuple[float, float, float, float]:
    """(E0, V0, E1, V1): mean and variance of PEP(S) for S ~ N(mu0, 1) and S ~ N(mu1, 1)."""
    from scipy import integrate, stats

    def mom(mu: float) -> tuple[float, float]:
        f = stats.norm(mu, 1).pdf
        lo, hi = mu - 12, mu + 12
        m1 = integrate.quad(lambda s: float(pep([s], pi, mu0, mu1)[0]) * f(s), lo, hi,
                            limit=200)[0]  # fmt: skip
        m2 = integrate.quad(lambda s: float(pep([s], pi, mu0, mu1)[0]) ** 2 * f(s), lo, hi,
                            limit=200)[0]  # fmt: skip
        return m1, max(m2 - m1 * m1, 0.0)

    e0, v0 = mom(mu0)
    e1, v1 = mom(mu1)
    return e0, v0, e1, v1


# ----------------------------------------------------------------------------- files


def m3_path(work: Path, key: str) -> Path:
    name = M3_NAMES[key]
    nested = work / "mapping" / name
    return nested if nested.exists() else work / name


def work_paths(work: Path) -> dict[str, Path]:
    from gippa.io.work_tables import default_paths

    paths = default_paths(work)
    for key in M3_NAMES:
        paths[key] = m3_path(work, key)
    return paths


def load(work: Path) -> tuple[Any, list[Any]]:
    from gippa.io.work_tables import load_m4_inputs

    return load_m4_inputs(work_paths(work))


def errors(issues: Iterable[Any]) -> list[Any]:
    return [i for i in issues if i.severity == "error"]


def tables_in(work: Path) -> dict[str, bytes]:
    """Every file under ``work`` (relative path -> bytes)."""
    return {p.relative_to(work).as_posix(): p.read_bytes() for p in sorted(work.rglob("*"))
            if p.is_file()}  # fmt: skip


# ----------------------------------------------------------------------------- frames


def gene_rows(layout: Any) -> dict[str, Mapping[str, Any]]:
    return {str(r["locus_tag"]): r for r in layout.genes.to_dict("records")}


def targets(peptides: pd.DataFrame) -> pd.DataFrame:
    return peptides[~peptides["is_decoy"].astype(bool)]


def decoys(peptides: pd.DataFrame) -> pd.DataFrame:
    return peptides[peptides["is_decoy"].astype(bool)]


def false_set(truth: pd.DataFrame) -> set[str]:
    t = truth[truth["is_false"].map(lambda v: bool(v) if not is_na(v) else False)]
    return set(t["peptide"].astype(str))


def n_obs_share(layout: Any, tags: Sequence[str]) -> float:
    """Share of the genome-wide observable-peptide count held by ``tags``."""
    w = {str(t): n_obs(a) for t, a in zip(layout.genes["locus_tag"], layout.genes["aa"],
                                          strict=True)}  # fmt: skip
    tot = sum(w.values())
    return sum(w[t] for t in tags) / tot
