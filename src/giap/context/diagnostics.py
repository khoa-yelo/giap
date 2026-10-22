"""Per-replicon island diagnostics: l-hat, l' and the lag-1 Moran's I of the detrended TU scores.

ADR-0021 **Amendment 1** (decision 2a) and M5-06 decision 48 dropped the ``uncalibrated`` flag
and put these three numbers per replicon in its place. ``IslandScan`` dies with the step, so
:mod:`giap.steps.context_stats` persists this frame as ``work/context/island_diagnostics.parquet``
and the §9.4 run report will consume it rather than recompute (M5-09 lead decision 1).

``IslandScan`` reports ``l_hat = 0`` both for a replicon that was not scanned and for one whose
``context.islands.block_len`` is fixed, and ``block_len = 0`` when it was not scanned. Printed as
a number, l-hat = 0 reads "no autocorrelation" — the false reassurance Amendment 1 refuses — so
here those sentinels are **NA, never 0**, and :data:`BLOCK_RULES` keeps the two reasons apart:
``auto`` (l-hat and l' >= 1, finite Moran's I, calibration applies), ``fixed`` (l-hat NA, l' the
configured value, finite Moran's I) and ``not_scanned`` (all missing, flag ``few_values``).
Only ``auto`` is inside the calibrated regime. See docs/methods/context_tables.md.

Pure functions; no randomness and no file I/O.
"""

from __future__ import annotations

__all__ = [
    "BLOCK_RULES",
    "DIAGNOSTICS_COLUMNS",
    "block_rule",
    "calibration_applies",
    "island_diagnostics",
    "render_l_hat",
]

import math
from typing import Any

import pandas as pd

from giap.config.model import GiapConfig
from giap.context.islands import IslandScan

#: Columns of the persisted frame, in order.
DIAGNOSTICS_COLUMNS: tuple[str, ...] = (
    "contig", "scanned", "block_rule", "l_hat", "block_len", "morans_i_lag1", "flags",
)  # fmt: skip
#: How the bootstrap block length of a replicon was obtained.
BLOCK_RULES: tuple[str, ...] = ("auto", "fixed", "not_scanned")
#: What :func:`render_l_hat` prints when l-hat was never estimated (never the character ``0``).
_NOT_ESTIMATED = {
    "fixed": "not estimated (block length fixed by configuration)",
    "not_scanned": "not estimated (replicon not scanned)",
}


def block_rule(*, scanned: bool, configured: int | str) -> str:
    """:data:`BLOCK_RULES` entry of one replicon: ``not_scanned`` when it never entered the
    scan, else ``fixed`` when ``context.islands.block_len`` is an integer, else ``auto``."""
    if not scanned:
        return "not_scanned"
    return "fixed" if isinstance(configured, int) else "auto"


def calibration_applies(rule: str) -> bool:
    """Whether the ADR-0021 FWER calibration covers a replicon with this block rule: True only
    for ``auto``, since Amendment 1's evidence is conditional on l' = min(20, 2 l-hat + 1, T).
    An unknown rule raises ``ValueError``."""
    if rule not in BLOCK_RULES:
        raise ValueError(f"unknown block rule {rule!r}; expected one of {list(BLOCK_RULES)}")
    return rule == "auto"


def render_l_hat(l_hat: Any, rule: str) -> str:
    """l-hat as the §9.4 run report should print it: the integer when it exists, else a
    sentence saying it was not estimated and why. ``0``, None and NaN all mean "never
    estimated", and the text never contains the character ``0`` (ADR-0021 Amendment 1)."""
    if rule not in BLOCK_RULES:
        raise ValueError(f"unknown block rule {rule!r}; expected one of {list(BLOCK_RULES)}")
    missing = l_hat is None or l_hat is pd.NA
    if not missing:
        value = float(l_hat)
        missing = math.isnan(value) or value < 1
    if missing:
        return _NOT_ESTIMATED.get(rule, _NOT_ESTIMATED["not_scanned"])
    return str(int(l_hat))


def _positive(value: int) -> int | None:
    """An ``IslandScan`` length: ``0`` is the "not estimated" sentinel, so it becomes None."""
    return int(value) if int(value) >= 1 else None


def island_diagnostics(scan: IslandScan, cfg: GiapConfig) -> pd.DataFrame:
    """:data:`DIAGNOSTICS_COLUMNS` for a :func:`giap.context.islands.scan_replicons` result,
    one row per replicon in C1 order: ``l_hat`` and ``block_len`` as nullable ``Int64`` (NA,
    never 0, when not estimated), ``morans_i_lag1`` as ``float64`` (NaN when not scanned) and
    the replicon's run ``flags``. ``cfg.context.islands.block_len`` tells ``fixed`` from
    ``auto``."""
    configured = cfg.context.islands.block_len
    rules = [block_rule(scanned=bool(ok), configured=configured) for ok in scan.scanned]
    estimated = [
        _positive(v) if rule == "auto" else None for rule, v in zip(rules, scan.l_hat, strict=True)
    ]
    return pd.DataFrame(
        {
            "contig": pd.Series([str(ts.contig) for ts in scan.series], dtype=object),
            "scanned": pd.Series([bool(ok) for ok in scan.scanned], dtype=bool),
            "block_rule": pd.Series(rules, dtype=object),
            "l_hat": pd.Series(estimated, dtype="Int64"),
            "block_len": pd.Series([_positive(v) for v in scan.block_len], dtype="Int64"),
            "morans_i_lag1": pd.Series([float(v) for v in scan.morans_i_lag1], dtype="float64"),
            "flags": pd.Series([tuple(f) for f in scan.flags], dtype=object),
        },
        columns=list(DIAGNOSTICS_COLUMNS),
    )
