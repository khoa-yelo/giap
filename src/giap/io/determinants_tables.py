"""Writers for ``tables/determinants_{features,model,attribution}.tsv`` (GIAP_SPEC.md §7.8.3,
§9.3, §9.5; docs/specs/M8-10-determinants-step-and-run.md).

Columns ``giap.determinants.features.FEATURE_COLUMNS``, ``giap.determinants.models.
TABLE_COLUMNS`` and ``giap.determinants.attribution.ATTRIBUTION_COLUMNS``; through M1-04
``write_tsv`` (numbers to 6 significant digits, booleans ``true``/``false``, missing ``NA``). An
empty ``flags`` cell is written ``NA``. Schemas: ``schemas/tables/determinants_*.schema.json``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

from giap.determinants.attribution import ATTRIBUTION_COLUMNS, SHAP_COLUMNS
from giap.determinants.features import FEATURE_COLUMNS, KINDS, Kind
from giap.determinants.models import TABLE_COLUMNS
from giap.io.tsv import write_tsv

__all__ = [
    "write_determinants_attribution", "write_determinants_features", "write_determinants_model",
]  # fmt: skip

_MODEL_KINDS: Mapping[str, Kind] = {
    "model": "str", "feature_or_block": "str", "metric": "str", "value": "float",
    "ci_low": "float", "ci_high": "float", "n": "int", "n_perm": "int", "flags": "str",
}  # fmt: skip
_ATTRIBUTION_KINDS: Mapping[str, Kind] = {
    "locus_tag": "str", "abundance_log2": "float", "predicted_log2": "float", "residual": "float",
    **dict.fromkeys(SHAP_COLUMNS, "float"), "dominant_block": "str", "top_decile": "bool",
    "flags": "str",
}  # fmt: skip


def _cell(v: object, kind: Kind) -> object:
    if v is None or pd.isna(v) or (kind == "str" and v == ""):
        return None
    if kind == "int":
        return int(v)  # type: ignore[call-overload]
    if kind == "float":
        return float(v)  # type: ignore[arg-type]
    if kind == "bool":
        return bool(v)
    return str(v)


def _write(path: Path, df: pd.DataFrame, columns: Sequence[str], kinds: Mapping[str, Kind]) -> None:
    rows = [
        [_cell(v, kinds[c]) for c, v in zip(columns, r, strict=True)]
        for r in df[list(columns)].itertuples(index=False, name=None)
    ]
    write_tsv(path, columns, rows)


def write_determinants_features(path: Path, df: pd.DataFrame) -> None:
    """``feature_frame`` output (M8-07) to ``path``."""
    _write(path, df, FEATURE_COLUMNS, KINDS)


def write_determinants_model(path: Path, df: pd.DataFrame) -> None:
    """``determinants_model.tsv`` rows (§9.3 plus ``flags``), in the given order, to ``path``."""
    _write(path, df, TABLE_COLUMNS, _MODEL_KINDS)


def write_determinants_attribution(path: Path, df: pd.DataFrame) -> None:
    """``attribution_table`` output (M8-09) to ``path``; no rows -> header only."""
    _write(path, df, ATTRIBUTION_COLUMNS, _ATTRIBUTION_KINDS)
