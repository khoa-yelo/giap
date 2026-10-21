"""Per-replicon island diagnostics: l-hat, l' and the lag-1 Moran's I of the detrended TU scores.

ADR-0021 **Amendment 1** (decision 2a) and M5-06 decision 48 dropped the ``uncalibrated`` flag:
"no per-replicon statistic certifies the domain, so a flag would give false reassurance ...
Instead: docs/methods/islands.md states the domain, and the §9.4 run report carries l-hat, l'
(``block_len``) and the lag-1 Moran's I per replicon so a reader can see where their data sit."
``IslandScan`` dies with the step, so :mod:`giap.steps.context_stats` persists this frame as
``work/context/island_diagnostics.parquet``; when the §9.4 run report is built it **consumes
that parquet** rather than recomputing anything (M5-09 lead decision 1).

**0 is not a value here.** ``IslandScan.l_hat`` is ``0`` both when the replicon was not scanned
and when ``context.islands.block_len`` fixes the block length, and ``IslandScan.block_len`` is
``0`` when the replicon was not scanned. Printed as a number, l-hat = 0 reads "no autocorrelation",
which is exactly the false reassurance Amendment 1 refuses. This module turns those sentinels
into missing values and keeps the two reasons apart in :data:`BLOCK_RULES`:

===============  ======  ===========  =========  ===============  ===========
``block_rule``   l-hat       l'           the lag-1 Moran's I         calibration      flags
===============  ======  ===========  =========  ===============  ===========
``auto``         ≥ 1     ≥ 1          finite     applies          none
``fixed``        NA      configured   finite     does **not**     none
``not_scanned``  NA      NA           NaN        does **not**     ``few_values``
===============  ======  ===========  =========  ===============  ===========

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
    """Which :data:`BLOCK_RULES` entry describes one replicon.

    Parameters
    ----------
    scanned : bool
        Whether the replicon entered the scan (``IslandScan.scanned``).
    configured : int or str
        ``context.islands.block_len``: an integer fixes l', ``"auto"`` estimates l-hat.

    Returns
    -------
    str
        ``"not_scanned"``, ``"fixed"`` or ``"auto"``.
    """
    if not scanned:
        return "not_scanned"
    return "fixed" if isinstance(configured, int) else "auto"


def calibration_applies(rule: str) -> bool:
    """Whether the ADR-0021 FWER calibration covers a replicon with this block rule.

    Parameters
    ----------
    rule : str
        One of :data:`BLOCK_RULES`.

    Returns
    -------
    bool
        True only for ``"auto"``: Amendment 1's FWER evidence is conditional on ADR-0021
        decision 1's rule l' = min(20, 2 l-hat + 1, T), so an overridden ``block_len`` — or a
        replicon that was never scanned — is outside it.

    Raises
    ------
    ValueError
        If ``rule`` is not a :data:`BLOCK_RULES` member.
    """
    if rule not in BLOCK_RULES:
        raise ValueError(f"unknown block rule {rule!r}; expected one of {list(BLOCK_RULES)}")
    return rule == "auto"


def render_l_hat(l_hat: Any, rule: str) -> str:
    """l-hat as the §9.4 run report should print it.

    Parameters
    ----------
    l_hat : int, float, None or pandas.NA
        The estimate; ``0``, ``None`` and NaN all mean "never estimated".
    rule : str
        One of :data:`BLOCK_RULES`.

    Returns
    -------
    str
        The integer as text when l-hat exists, otherwise a sentence saying it was not estimated
        and why. Never the character ``0``: printed as a number, l-hat = 0 would read "no
        autocorrelation" (ADR-0021 Amendment 1).

    Raises
    ------
    ValueError
        If ``rule`` is not a :data:`BLOCK_RULES` member.
    """
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
    """l-hat, l' and the lag-1 Moran's I per replicon, with 0 replaced by a missing value.

    Parameters
    ----------
    scan : giap.context.islands.IslandScan
        The scan of :func:`giap.context.islands.scan_replicons`; its replicons are in C1 order.
    cfg : giap.config.model.GiapConfig
        ``context.islands.block_len`` tells ``fixed`` from ``auto``.

    Returns
    -------
    pandas.DataFrame
        :data:`DIAGNOSTICS_COLUMNS`, one row per replicon in C1 order: ``scanned`` (bool),
        ``block_rule`` (:data:`BLOCK_RULES`), ``l_hat`` and ``block_len`` (nullable ``Int64``,
        NA and never 0 when not estimated), ``morans_i_lag1`` (``float64``, NaN when the
        replicon was not scanned) and ``flags`` (the replicon's run flags, e.g.
        ``few_values``).
    """
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
