"""Variant-module QC counts for ``detection_summary.tsv`` (GIAP_SPEC.md §7.5.1 row 1, §9.3;
docs/specs/M7-14-te2-proxy-and-id-mixups.md rule 3).

``id_mixups``: a peptide searched against this strain's own Bakta proteins (``database``,
``search_db = same``) cannot contradict the genome, except at selenocysteine sites and non-ATG
starts, which carry the M7-05 flags ``stop_codon`` / ``start_codon``. Every other ALT site such a
peptide supports is an identification mix-up. Confidence does not matter (lead decision 1).
Pure functions.
"""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

__all__ = ["EXEMPT_FLAGS", "SUPPORT_COLUMNS", "id_mixups", "site_support"]

#: the M7-05 flags of the expected same-db ALT sites (Sec sites, non-ATG starts)
EXEMPT_FLAGS = frozenset({"start_codon", "stop_codon"})
#: one row per (ALT site, supporting peptide)
SUPPORT_COLUMNS = ("variant_id", "peptide", "source", "search_db")


def _flag_set(v: object) -> set[str]:
    if isinstance(v, str):
        return {x for x in v.split(";") if x}
    if isinstance(v, Iterable):
        return {str(x) for x in v}
    return set()


def id_mixups(sites: pd.DataFrame, support: pd.DataFrame) -> tuple[int, int]:
    """``(count, n)``: ``n`` = ALT sites of ``sites`` (``variants.tsv`` rows) supported by at
    least one ``database`` peptide with ``search_db = same`` in ``support`` (columns
    :data:`SUPPORT_COLUMNS`); ``count`` = those flagged neither ``start_codon`` nor
    ``stop_codon``. Support rows of unknown sites are ignored. Pure."""
    if not len(sites) or not len(support):
        return 0, 0
    same = support["source"].astype(str).eq("database") & support["search_db"].astype(str).eq(
        "same"
    )
    same_ids = set(support.loc[same, "variant_id"].astype(str))
    n = count = 0
    for vid, flags in zip(sites["variant_id"].astype(str), sites["flags"], strict=True):
        if vid in same_ids:
            n += 1
            count += not (EXEMPT_FLAGS & _flag_set(flags))
    return int(count), int(n)


def site_support(
    variant_ids: Iterable[str],
    alt_peptides: Iterable[Iterable[str]],
    sources: dict[str, tuple[str, str]],
) -> pd.DataFrame:
    """The long ``support`` frame of :func:`id_mixups`: per site (``variant_ids`` aligned with
    ``alt_peptides``, the ALT peptides of each site), one row per distinct peptide with its
    ``(source, search_db)`` from ``sources``. Rows sorted by site order, then peptide. Pure."""
    rows = [
        (str(vid), p, *sources[p])
        for vid, peps in zip(variant_ids, alt_peptides, strict=True)
        for p in sorted(set(peps))
    ]
    df = pd.DataFrame(rows, columns=list(SUPPORT_COLUMNS))
    return df.astype(dict.fromkeys(SUPPORT_COLUMNS, object))
