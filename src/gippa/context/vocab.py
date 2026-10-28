"""Controlled flag vocabulary of the region-and-context tables (GIPPA_SPEC.md §9.4).

§9.4 is closed ("extend only through the schema"): the two M5 table schemas constrain their
``flags`` column to exactly :data:`CONTEXT_FLAGS`, and the writers reject anything else. The set
is the §9.4 "Region and context" row **minus** ``uncalibrated`` (ADR-0021 Amendment 1 and M5-06
decision 48 drop it as false reassurance; :mod:`gippa.context.diagnostics` replaces it), plus the
M4 additions, the M5 additions of lead decision 36 and the ADR-0061 flags (M11-14d). The legacy
``regions_stats.schema.json`` pattern still lists ``uncalibrated`` and is left alone.
"""

from __future__ import annotations

__all__ = ["CONTEXT_FLAGS", "flags_pattern", "require_known_flags"]

from collections.abc import Iterable

#: §9.4 "Region and context" flags, minus the dropped ``uncalibrated``.
_SPEC_FLAGS: tuple[str, ...] = (
    "small_region",
    "dosage_plus_regulation",
    "descriptive",
    "no_oric",
)
#: M4 additions (``regions_stats`` and ``detection_summary``; M4-16 / M4-17 / M4-18).
_M4_FLAGS: tuple[str, ...] = (
    "few_shifts",
    "split_contigs",
    "hl_variance_nonpositive",
    "no_reference_replicon",
    "depth_adjusted",
    "ci_few_units",
)
#: M5 additions (lead decision 36): no statistic, too few values, no CI.
_M5_FLAGS: tuple[str, ...] = ("not_estimable", "few_values", "ci_unestimated")

#: ADR-0061 (M11-14d): TU-block Moran's I rows assume exact TU boundaries (every run) and are
#: anti-conservative with distance-model TUs; only ``autocorrelation`` rows carry them.
_ADR0061_FLAGS: tuple[str, ...] = ("tu_block_tu_error", "tu_source_distance_model")

#: The closed vocabulary, in §9.4 order then milestone order.
CONTEXT_FLAGS: tuple[str, ...] = (*_SPEC_FLAGS, *_M4_FLAGS, *_M5_FLAGS, *_ADR0061_FLAGS)

_KNOWN = frozenset(CONTEXT_FLAGS)


def flags_pattern() -> str:
    """Frictionless ``pattern`` constraint for a ``;``-joined ``flags`` cell.

    Returns
    -------
    str
        A regular expression matching :data:`CONTEXT_FLAGS` members joined by ``;``.
    """
    alt = "|".join(CONTEXT_FLAGS)
    return f"(?:{alt})(?:;(?:{alt}))*"


def require_known_flags(flags: Iterable[str], *, where: str) -> tuple[str, ...]:
    """Check ``flags`` against the closed vocabulary and return them as a tuple of ``str``.

    Parameters
    ----------
    flags : iterable of str
        Flags of one table row.
    where : str
        Name used in the error message (the table being written).

    Returns
    -------
    tuple of str
        ``flags`` as strings, in the given order.

    Raises
    ------
    ValueError
        If a flag is not in :data:`CONTEXT_FLAGS`; the message names the offending flag.
    """
    out = tuple(str(f) for f in flags)
    unknown = [f for f in out if f not in _KNOWN]
    if unknown:
        raise ValueError(
            f"{where}: flag(s) outside the §9.4 vocabulary: {', '.join(sorted(unknown))}"
        )
    return out
