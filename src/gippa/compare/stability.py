"""Jaccard overlap of five result families between runs (M9-22; GIPPA_SPEC.md §7.9.3).

Tables are read through ``RunView`` (the I/O layer); the set extraction and the Jaccard index are
pure. A family whose table is absent in a run gives NA with ``missing_in:<label>``; two empty sets
give NA with ``both_empty``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Hashable, Mapping, Sequence
from itertools import combinations

import pandas as pd

from gippa.config import GippaConfig
from gippa.io.run_view import RunView

__all__ = ["COLUMNS", "FAMILIES", "family_set", "jaccard_row", "stability"]

COLUMNS = ("family", "run_a", "run_b", "n_a", "n_b", "n_shared", "jaccard", "flags")
FAMILIES = ("islands", "operon_flags", "enriched_categories", "novel_orfs", "alt_sites")
KeySet = frozenset[Hashable]


def _le(frame: pd.DataFrame, col: str, alpha: float) -> pd.DataFrame:
    return frame[(frame[col] <= alpha).fillna(False).astype(bool)]


def _islands(t: Mapping[str, pd.DataFrame], cfg: GippaConfig) -> KeySet:
    sel = _le(t["islands"], "p_value", cfg.stats.alpha.context)
    return frozenset(g for cell in sel["genes"].dropna() for g in str(cell).split(";") if g)


def _operon_flags(t: Mapping[str, pd.DataFrame], cfg: GippaConfig) -> KeySet:
    sel = _le(t["operon_flags"], "q_value", cfg.stats.alpha.operon_flags)
    return frozenset(zip(sel["gene_a"], sel["gene_b"], strict=True))


def _enriched(t: Mapping[str, pd.DataFrame], cfg: GippaConfig) -> KeySet:
    alpha = cfg.stats.alpha.enrichment
    out: set[Hashable] = set()
    if (ora := t.get("enrichment_detection")) is not None:
        sel = _le(ora, "q_value", alpha)
        out |= {("enrichment_detection", c, term, d)
                for c, term, d in zip(sel["collection"], sel["term"], sel["direction"],
                                      strict=True)}  # fmt: skip
    if (gsea := t.get("enrichment_gsea")) is not None:
        sel = _le(gsea, "q_value", alpha)
        out |= {("enrichment_gsea", c, term, "up" if nes > 0 else "down")
                for c, term, nes in zip(sel["collection"], sel["term"], sel["nes"],
                                        strict=True)}  # fmt: skip
    return frozenset(out)


def _novel(t: Mapping[str, pd.DataFrame], cfg: GippaConfig) -> KeySet:
    sel = t["novel_orfs"]
    return frozenset(sel.loc[sel["called"].fillna(False).astype(bool), "orf_id"])


def _alt(t: Mapping[str, pd.DataFrame], cfg: GippaConfig) -> KeySet:
    sel = t["variants"][t["variants"]["call"] == "ALT"]
    keys = zip(sel["contig"], sel["codon_start"], sel["strand"], sel["alt_aa"], strict=True)
    return frozenset((c, int(s), st, aa) for c, s, st, aa in keys)


Extract = Callable[[Mapping[str, pd.DataFrame], GippaConfig], KeySet]
#: family -> (tables read, extractor, any one table suffices)
_SPEC: dict[str, tuple[tuple[str, ...], Extract, bool]] = {
    "islands": (("islands",), _islands, False),
    "operon_flags": (("operon_flags",), _operon_flags, False),
    "enriched_categories": (("enrichment_detection", "enrichment_gsea"), _enriched, True),
    "novel_orfs": (("novel_orfs",), _novel, False),
    "alt_sites": (("variants",), _alt, False),
}


def family_set(view: RunView, family: str, cfg: GippaConfig) -> KeySet | None:
    """The family's key set of one run, or ``None`` when its table is absent (or the
    ``function`` module is off, for ``enriched_categories``)."""
    keys, extract, any_one = _SPEC[family]
    if family == "enriched_categories" and "function" in view.modules_off():
        return None
    tables = {k: t for k in keys if (t := view.table(k)) is not None}
    if not tables or (not any_one and len(tables) < len(keys)):
        return None
    return extract(tables, cfg)


def jaccard_row(
    family: str, la: str, lb: str, a: KeySet | None, b: KeySet | None
) -> dict[str, object]:
    """|A & B| / |A | B|; NA with ``missing_in`` (absent table) or ``both_empty``."""
    row: dict[str, object] = {"family": family, "run_a": la, "run_b": lb, "jaccard": math.nan,
                              "n_a": None if a is None else len(a),
                              "n_b": None if b is None else len(b), "n_shared": None,
                              "flags": None}  # fmt: skip
    missing = [f"missing_in:{lab}" for lab, s in ((la, a), (lb, b)) if s is None]
    if missing:
        row["flags"] = ";".join(missing)
    elif a is not None and b is not None:
        row["n_shared"] = len(a & b)
        union = len(a | b)
        if union == 0:
            row["flags"] = "both_empty"
        else:
            row["jaccard"] = len(a & b) / union
    return row


def stability(runs: Mapping[str, RunView], labels: Sequence[str], cfg: GippaConfig) -> pd.DataFrame:
    """One row per family x unordered pair (families in ``FAMILIES`` order, pairs in label
    order) with the ``compare_stability`` columns."""
    sets = {(lab, fam): family_set(runs[lab], fam, cfg) for lab in labels for fam in FAMILIES}
    rows = [
        jaccard_row(fam, la, lb, sets[la, fam], sets[lb, fam])
        for fam in FAMILIES for la, lb in combinations(labels, 2)
    ]  # fmt: skip
    out = pd.DataFrame(rows, columns=list(COLUMNS))
    return out.astype({"n_a": "Int64", "n_b": "Int64", "n_shared": "Int64"})
