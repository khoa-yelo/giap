"""Class-specific false-ALT rates from paired de novo and database results (GIAP_SPEC.md §7.5.4,
§8.1 row "ALT call is real", §7.5.5 ``CLASS_ERR``; docs/specs/M7-10-class-error-rates.md).

A de novo row pairs with a confident (``q_value <= q_max``) target reference row of the same
``(run, spectrum_id)``; it is an observation of every class ``(alt_type, mass_equivalent)`` of
the ALT sites its peptide supports, and an error when the two peptides differ after I/L folding.
The rate of a class is ``errors / n``; a class without observations is ``unestimated``. Pure.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

import pandas as pd

__all__ = [
    "UNESTIMATED",
    "ClassRate",
    "alt_support",
    "annotate_sites",
    "class_rates",
    "paired_observations",
    "site_class",
]

UNESTIMATED = "unestimated"
_KEY = ["run", "spectrum_id"]
_CLASS = ["alt_type", "mass_equivalent"]


@dataclass(frozen=True, slots=True)
class ClassRate:
    alt_type: str
    mass_equivalent: bool
    n: int  # paired ALT-supporting de novo observations
    errors: int
    rate: float | None  # errors / n; None when n == 0


def site_class(alt_type: object, flags: object) -> tuple[str, bool]:
    """``(alt_type, "mass_equivalent" in flags)``; frameshifts are ``(frameshift, False)``."""
    t = str(alt_type)
    fl = tuple(flags) if isinstance(flags, Iterable) and not isinstance(flags, str) else ()
    return t, t != "frameshift" and "mass_equivalent" in fl


def alt_support(sites: pd.DataFrame, alt_peptides: Sequence[Iterable[str]]) -> pd.DataFrame:
    """One row per (supporting peptide, site): ``variant_id``, ``peptide``, ``alt_type``,
    ``mass_equivalent``; ``alt_peptides`` is aligned with the rows of ``sites``."""
    flags = sites["flags"] if "flags" in sites.columns else [()] * len(sites)
    cols = zip(sites["variant_id"], sites["alt_type"], flags, alt_peptides, strict=True)
    rows = [(vid, str(p), *site_class(t, f)) for vid, t, f, peps in cols for p in peps]
    return pd.DataFrame(rows, columns=["variant_id", "peptide", *_CLASS]).astype(
        {"mass_equivalent": bool}
    )


def _keyed(df: pd.DataFrame) -> pd.DataFrame:
    """``run`` (NA as ``""``), ``spectrum_id`` and ``peptide`` as text; NA spectrum IDs dropped."""
    has = df["spectrum_id"].notna().to_numpy(dtype=bool)
    out = df.loc[has, ["run", "spectrum_id", "peptide"]].copy()
    out["run"] = out["run"].astype(object).where(out["run"].notna(), "").astype(str)
    for c in ("spectrum_id", "peptide"):
        out[c] = out[c].astype(str)
    return out.reset_index(drop=True)


def _fold(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace("I", "L", regex=False)


def paired_observations(
    denovo: pd.DataFrame,
    reference: pd.DataFrame,
    alt_support: pd.DataFrame,
    q_max: float = 0.01,
) -> pd.DataFrame:
    """One row per (distinct de novo ``(run, spectrum_id, peptide)``, class) of rules 1-3:
    ``run``, ``spectrum_id``, ``peptide``, ``reference_peptide``, ``alt_type``,
    ``mass_equivalent`` and ``error``; rows repeated over charges count once (lead decision 8).

    ``reference`` rows count when not decoys and ``q_value <= q_max``; a table without the
    ``q_value`` / ``is_decoy`` columns is taken as already filtered (C17, lead decision 5).
    """
    keep = pd.Series(True, index=reference.index)
    if "q_value" in reference.columns:
        keep &= pd.to_numeric(reference["q_value"], errors="coerce").le(q_max)
    if "is_decoy" in reference.columns:
        keep &= ~reference["is_decoy"].astype("boolean").fillna(False).astype(bool)
    ref = _keyed(reference[keep.to_numpy(dtype=bool)])
    ref = ref.drop_duplicates(_KEY, keep="last").rename(columns={"peptide": "reference_peptide"})
    classes = alt_support[["peptide", *_CLASS]].drop_duplicates().astype({"peptide": str})
    dn = _keyed(denovo).drop_duplicates()  # one observation per (run, spectrum_id, peptide)
    obs = dn.merge(ref, on=_KEY).merge(classes, on="peptide")
    obs["mass_equivalent"] = obs["mass_equivalent"].astype(bool)
    obs["error"] = (_fold(obs["peptide"]) != _fold(obs["reference_peptide"])).astype(bool)
    return obs.reset_index(drop=True)


def class_rates(obs: pd.DataFrame) -> list[ClassRate]:
    """Per class of ``obs``: ``n``, ``errors`` and ``rate = errors / n`` (rule 4)."""
    if len(obs) == 0:
        return []
    g = obs.astype({"error": bool}).groupby(_CLASS, sort=True)["error"].agg(["size", "sum"])
    return [
        ClassRate(str(t), bool(me), int(n), int(e), int(e) / int(n))
        for (t, me), n, e in zip(g.index, g["size"], g["sum"], strict=True)
    ]


def annotate_sites(sites: pd.DataFrame, rates: Sequence[ClassRate]) -> pd.DataFrame:
    """A copy of ``sites`` with ``class_error_rate`` (float or ``unestimated``) and
    ``class_error_n`` (Int64, NA when unestimated) of each row's class."""
    known = {(r.alt_type, bool(r.mass_equivalent)): r for r in rates if r.n > 0}
    flags = sites["flags"] if "flags" in sites.columns else [()] * len(sites)
    hits = [known.get(site_class(t, f)) for t, f in zip(sites["alt_type"], flags, strict=True)]
    out = sites.copy()
    rate: list[Any] = [UNESTIMATED if h is None else h.rate for h in hits]
    out["class_error_rate"] = pd.Series(rate, index=out.index, dtype=object)
    n = pd.array([pd.NA if h is None else h.n for h in hits], dtype="Int64")
    out["class_error_n"] = pd.Series(n, index=out.index)
    return out
