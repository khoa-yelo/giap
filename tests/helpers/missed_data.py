"""Hand-built pair frames for M6-06 missed operons (docs/specs/M6-06-missed-operons.md).

The frames have the M6-04 ``pairs.tsv`` columns (``pairs_data.PAIR_COLS``) without going through
``pair_frame``: each pair gets its own two genes, ``abundance_log2_a = 10`` and
``abundance_log2_b = 10 + delta`` so ``delta_log2 = delta`` and the score is ``s = -|delta|``.
Distance bins follow ``operons.distance_bins_bp`` = [-inf, 0, 50, 150, 300, inf]: ``BP[b]`` is an
``intergenic_bp`` inside bin ``b``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from giap.config import GiapConfig
from helpers.pairs_data import PAIR_COLS, TEST_COLS

#: one ``intergenic_bp`` inside each of the five default distance bins
BP: tuple[int, ...] = (-10, 20, 100, 200, 500)
BASE_LOG2 = 10.0


def block(
    cls: str,
    deltas: Iterable[float],
    bin_: int = 1,
    *,
    detected: tuple[bool, bool] = (True, True),
    flags: Sequence[str] = (),
) -> pd.DataFrame:
    """Pairs of class ``cls`` with the given ``delta_log2``, all in distance bin ``bin_``."""
    d = np.asarray(list(deltas), dtype=np.float64)
    n = d.size
    out = pd.DataFrame({c: [np.nan] * n for c in PAIR_COLS})
    out["class"] = cls
    out["strand"] = "." if cls in ("CV", "DV") else "+"
    out["intergenic_bp"] = np.full(n, BP[bin_], dtype=np.int64)
    out["uniop_p"] = 0.1 if cls == "NO" else (0.9 if cls == "OP" else np.nan)
    out["detected_a"], out["detected_b"] = detected
    out["p_detect_a"] = out["p_detect_b"] = 0.9
    out["abundance_log2_a"] = BASE_LOG2
    out["abundance_log2_b"] = BASE_LOG2 + d
    out["delta_log2"] = d
    out["test"] = "none"
    for c in TEST_COLS:
        out[c] = None if c in ("method", "null_model") else np.nan
    out["n_perm"] = pd.array([pd.NA] * n, dtype="Int64")
    out["flags"] = [list(flags) for _ in range(n)]
    return out


def stack(*blocks: pd.DataFrame) -> pd.DataFrame:
    """Concatenate blocks and number pairs ``PAIR_1..n`` with genes ``G<n>a``/``G<n>b``."""
    df = pd.concat(blocks, ignore_index=True)
    k = np.arange(1, len(df) + 1)
    df["pair_id"] = [f"PAIR_{i}" for i in k]
    df["gene_a"] = [f"G{i}a" for i in k]
    df["gene_b"] = [f"G{i}b" for i in k]
    return df[list(PAIR_COLS)]


def levels(spec: Mapping[int, Sequence[tuple[float, int]]], cls: str) -> list[pd.DataFrame]:
    """Blocks from ``{bin: [(abs_delta, count), ...]}``; signs alternate so only |delta| counts."""
    out = []
    for bin_, items in spec.items():
        for ad, count in items:
            sign = np.where(np.arange(count) % 2 == 0, 1.0, -1.0)
            out.append(block(cls, sign * ad, bin_))
    return out


# M6-01 examples 4-6 as pairs: |delta| levels 0.1 < 0.5 < 1.0 stand for the scores 1 > 0 > -1.
# At |delta| <= 0.1: 40 targets, 5 + 2 CV decoys, 1 + 1 SW decoys.
TARGETS = {0: [(0.1, 30), (0.5, 40), (1.0, 130)], 1: [(0.1, 10), (0.5, 20), (1.0, 70)]}
CV = {0: [(0.1, 5), (0.5, 15), (1.0, 80)], 1: [(0.1, 2), (0.5, 13), (1.0, 85)]}
SW = {0: [(0.1, 1), (1.0, 49)], 1: [(0.1, 1), (1.0, 19)]}
EX1 = 0.3712871287128713  # (200*6/101 + 100*3/101)/40, M6-01 example 4
EX2 = 0.4341736694677871  # (200*2/51 + 100*2/21)/40, M6-01 example 6


def example_frame(sandwich: bool = True) -> pd.DataFrame:
    parts = levels(TARGETS, "NO") + levels(CV, "CV")
    if sandwich:
        parts += levels(SW, "SW")
    return stack(*parts)


def cfg(**operons: Any) -> GiapConfig:
    """Default config with ``operons`` overrides (``decoys``, ``pi0``, ``distance_bins_bp``)."""
    return GiapConfig.model_validate({"seed": 1, "operons": operons} if operons else {"seed": 1})


def genes_for(pairs: pd.DataFrame, **cols: Mapping[str, Any]) -> pd.DataFrame:
    """A §C-G-shaped gene table for every gene of ``pairs``; ``cols`` override per locus tag."""
    tags = list(dict.fromkeys([*pairs["gene_a"], *pairs["gene_b"]]))
    default = {
        "partial": False, "abundance_imputed": False, "abundance_type": "ibaq",
        "abundance_scale": "molar",
    }  # fmt: skip
    data: dict[str, list[Any]] = {"locus_tag": tags}
    for col, dflt in default.items():
        over = cols.get(col, {})
        data[col] = [over.get(t, dflt) for t in tags]
    return pd.DataFrame(data)


def null_frame(rng: np.random.Generator, n: int = 2000, planted: int = 0) -> pd.DataFrame:
    """``n`` NO pairs and ``n`` CV and ``n`` SW decoys, Delta ~ N(0, 1), bins 1 and 2 at random.

    The first ``planted`` NO pairs are co-abundant with |Delta| <= 1e-3 SD (example 6)."""
    d = rng.normal(size=n)
    d[:planted] = rng.uniform(-1e-3, 1e-3, planted)
    parts = []
    for cls, delta in (("NO", d), ("CV", rng.normal(size=n)), ("SW", rng.normal(size=n))):
        b = rng.choice((1, 2), size=n)
        parts += [block(cls, delta[b == k], k) for k in (1, 2)]
    df = pd.concat(parts, ignore_index=True)
    is_no = df["class"].to_numpy() == "NO"
    flag = is_no & df["delta_log2"].isin(d[:planted]).to_numpy()
    out = stack(df)
    out.attrs["planted"] = set(out.loc[flag, "pair_id"])
    return out
