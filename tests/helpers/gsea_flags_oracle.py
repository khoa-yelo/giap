"""Independent references for the M8-05b GSEA row flags (M8-05 spec, decision 9).

* ``delta_bruteforce``: the ``genomic_cluster`` statistic with P_w summed over *every* pair of
  ranked TUs (no closed form), so it checks the implementation's closed form as well as C.
* ``z_counts``: the hypergeometric z of the detected count from the four gene counts.
* ``flags_frame``: an M8-05 gene frame (``helpers.gsea_data.frame``) with the ``contig`` column
  of the gene table, one contig per TU.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping, Sequence

import pandas as pd

from helpers.gsea_data import Gene, frame

PHI0 = 0.5
DMAX = 10
TAU = 0.15
Z0 = 1.0
#: the decision-9 order of the ``flags`` cell
FLAG_ORDER = (
    "ci_unestimated",
    "unestimated",
    "size_filtered",
    "genomic_cluster",
    "detection_confounded",
)


def _dist(a: int, b: int, n: int, circular: bool) -> int:
    d = abs(a - b)
    return min(d, n - d) if circular else d


def _w(d: int) -> float:
    return PHI0**d if 1 <= d <= DMAX else 0.0


def delta_bruteforce(
    positions: Mapping[str, Sequence[int]],
    n_ranked: Mapping[str, int],
    circular: Mapping[str, bool],
) -> float:
    """delta = (C - E[C]) / (1 - k/n) by enumeration of all pairs (decision 9a)."""
    k = sum(len(v) for v in positions.values())
    n = sum(n_ranked.values())
    if k < 2 or k >= n:
        return 0.0
    s = 0.0
    for r, pos in positions.items():
        for a, b in itertools.combinations(sorted(pos), 2):
            s += _w(_dist(a, b, n_ranked[r], circular.get(r, False)))
    pw = 0.0
    for r, nr in n_ranked.items():
        for i, j in itertools.combinations(range(nr), 2):
            pw += _w(_dist(i, j, nr, circular.get(r, False)))
    c = 2.0 * s / k
    ec = 2.0 * (k - 1) * pw / (n * (n - 1))
    return (c - ec) / (1.0 - k / n)


def z_counts(M: int, N: int, n: int, x: int) -> float:
    """Hypergeometric z of the detected count; 0 when the variance is 0 (decision 9b)."""
    if M <= 1:
        return 0.0
    f = N / M
    v = n * f * (1 - f) * (M - n) / (M - 1)
    return 0.0 if v <= 0 else (x - n * f) / math.sqrt(v)


def detection_flag(z: float, es: float) -> bool:
    if math.isnan(es) or es == 0:
        return abs(z) >= Z0
    return z * es < 0 and abs(z) >= Z0


def flags_frame(
    tus: Mapping[str, Sequence[Gene]], contig_of: Mapping[str, str] | str = "chr"
) -> pd.DataFrame:
    """``helpers.gsea_data.frame`` plus ``contig`` (a TU -> contig map, or one contig for all)."""
    df = frame(tus)
    if isinstance(contig_of, str):
        df["contig"] = contig_of
    else:
        df["contig"] = df["tu_id"].map(dict(contig_of)).astype(object)
    return df


def ring_frame(n: int, values: Sequence[float] | None = None, contig: str = "chr") -> pd.DataFrame:
    """``n`` singleton TUs TU_1 .. TU_n on one contig, gene g<i> scoring with a distinct value."""
    vals = [20.0 + 0.01 * i for i in range(n)] if values is None else list(values)
    return flags_frame(
        {f"TU_{i + 1}": [Gene(f"g{i + 1}", float(vals[i]))] for i in range(n)}, contig
    )


def replicons_frame(topology: Mapping[str, str]) -> pd.DataFrame:
    """A C1 replicons frame (``contig``, ``length``, ``topology``, ``kind``)."""
    return pd.DataFrame(
        {
            "contig": list(topology),
            "length": [1_000_000] * len(topology),
            "topology": list(topology.values()),
            "kind": ["chromosome"] * len(topology),
        }
    )
