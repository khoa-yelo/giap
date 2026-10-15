"""Accepted evidence and per-gene evidence counts (docs/specs/M4-11-accepted-evidence.md).

Spec §7.3.1 (five conditions, distinct peptides), §7.2.6 (uniqueness, locus classes), §7.2.7
(genome relation), §9.3 (``genes.tsv`` evidence columns), §9.4 (gene flags), §6.5 (origin).

Lead decisions: ``pseudogene`` and ``cterm_extension`` placements assign to the gene (detection);
NA ``placement_q`` without ``same_db`` is rejected; C9 decoy rows are ignored; the C9 q wins over
C11 ``best_q_input``; ``nterm_extension`` and ``cterm_extension`` placements are not
``abundance_eligible``; contaminant peptides passing the five conditions are accepted and carried
to the gene row as the boolean ``contaminant_evidence`` (and ``is_contaminant`` on accepted rows).
When a (gene, peptide) pair has several accepted placements the representative row is the exact
one first, then the leftmost; every accepted placement of the pair is kept in its ``placements``
column and ``coverage_frac`` is the union of all accepted placements of the gene (review decision
5). Everything here is pure and vectorised.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from giap.genetable import contract

if TYPE_CHECKING:
    from giap.config import GiapConfig

# The M4-06 assigning classes plus ``pseudogene``, which the contract lets carry a gene and which
# assigns for detection (lead decision 1); derived from the contract, not duplicated.
ASSIGNING_CLASSES: tuple[str, ...] = (*contract.ASSIGNING_CLASSES, "pseudogene")
EXTENSION_CLASSES = ("nterm_extension", "cterm_extension")
ACCEPTED_COLUMNS = (
    "gene",
    "peptide",
    "q_value",
    "pep",
    "score",
    "intensity",
    "uniqueness",
    "route",
    "locus_class",
    "abundance_eligible",
    "start",
    "end",
    "contig",
    "crosses_origin",
    "is_contaminant",
)
# Every accepted placement of the pair as sorted tuples of ``_PLACEMENT_FIELDS`` (decision 5).
PLACEMENTS_COLUMN = "placements"
_PLACEMENT_FIELDS: dict[str, Any] = {
    "contig": str,
    "start": "int64",
    "end": "int64",
    "route": str,
    "uniqueness": str,
    "locus_class": str,
}
COUNT_COLUMNS = (
    "n_peptides_unique",
    "n_peptides_multi_identical",
    "n_peptides_shared",
    "n_peptides_tolerant",
)
EVIDENCE_COLUMNS = (
    "locus_tag",
    *COUNT_COLUMNS,
    "coverage_frac",
    "best_peptide_q",
    "detected_2pep",
    "n_accepted",
    "evidence_flags",
    "contaminant_evidence",
)
_DTYPES = {
    "q_value": "float64",
    "pep": "float64",
    "score": "float64",
    "intensity": "float64",
    "abundance_eligible": bool,
    "start": "int64",
    "end": "int64",
    "crosses_origin": bool,
    "is_contaminant": bool,
}


def _has_flag(cells: pd.Series, flag: str) -> np.ndarray:
    def one(cell: Any) -> bool:
        if isinstance(cell, (list, tuple, np.ndarray)):
            return flag in [str(f) for f in cell]
        return isinstance(cell, str) and flag in cell.split(";")

    return np.fromiter((one(c) for c in cells), dtype=bool, count=len(cells))


def _c9_values(hits: pd.DataFrame, peptides: pd.DataFrame | None) -> pd.DataFrame:
    """Per C11 row: q, pep, score, intensity, is_contaminant (C9 wins; ValueError if absent)."""
    if peptides is None:
        return pd.DataFrame(
            {
                "q_value": hits["best_q_input"].astype("float64").to_numpy(),
                "pep": hits["pep_input"].astype("float64").to_numpy(),
                "score": hits["score_input"].astype("float64").to_numpy(),
                "intensity": np.full(len(hits), np.nan),
                "is_contaminant": np.zeros(len(hits), dtype=bool),
            },
            index=hits.index,
        )
    target = peptides[~peptides["is_decoy"].astype(bool)] if len(peptides) else peptides
    c9 = target.drop_duplicates("peptide").set_index("peptide")
    keys = hits["peptide"].astype(str)
    missing = sorted(set(keys) - set(c9.index.astype(str)))
    if missing:
        raise ValueError(f"C11 peptides without a C9 row: {', '.join(missing[:10])}")
    out = pd.DataFrame(index=hits.index)
    for col in ("q_value", "pep", "score", "intensity"):
        out[col] = keys.map(c9[col].astype("float64")).astype("float64")
    out["is_contaminant"] = keys.map(c9["is_contaminant"].astype(bool)).astype(bool)
    return out


def _empty_accepted() -> pd.DataFrame:
    cols = (*ACCEPTED_COLUMNS, PLACEMENTS_COLUMN)
    df = pd.DataFrame({c: pd.Series([], dtype=object) for c in cols})
    return df.astype(_DTYPES)


def accepted_peptides(
    peptide_hits: pd.DataFrame, peptides: pd.DataFrame | None, cfg: GiapConfig
) -> pd.DataFrame:
    """One row per accepted (gene, peptide) pair (§7.3.1), with C9 values joined."""
    hits = peptide_hits
    if hits.empty:
        return _empty_accepted()
    c9 = _c9_values(hits, peptides)
    q = c9["q_value"].to_numpy()
    pq = hits["placement_q"].astype("float64").to_numpy()
    related = cfg.genome_relation == "related"
    route = hits["route"].astype(str).to_numpy()
    ok = (
        (np.isnan(q) | (q <= cfg.input.q_value_max))
        & ((pq <= cfg.stats.alpha.placement) | _has_flag(hits["flags"], "same_db"))
        & hits["gene"].notna().to_numpy()
        & hits["locus_class"].isin(ASSIGNING_CLASSES).to_numpy()
        & hits["uniqueness"].isin(("unique", "multi_identical")).to_numpy()
        & ((route == "exact") | related)
    )
    if not ok.any():
        return _empty_accepted()
    sel = hits.loc[ok, ["gene", "peptide", "uniqueness", "route", "locus_class"]].copy()
    sel = sel.astype({"gene": str, "peptide": str})
    sel[["start", "end", "contig", "crosses_origin"]] = hits.loc[
        ok, ["start", "end", "contig", "crosses_origin"]
    ]
    sel = pd.concat([sel, c9.loc[ok]], axis=1)
    sel["abundance_eligible"] = (
        (sel["uniqueness"] == "unique")
        & ~sel["locus_class"].isin(EXTENSION_CLASSES)
        & ((sel["route"] == "exact") | ((sel["route"] == "tolerant") & related))
    )
    sel["_tol"] = sel["route"] != "exact"
    sel = sel.sort_values(["gene", "peptide", "_tol", "start", "end"], kind="mergesort")
    fields = zip(*(sel[c].astype(t).tolist() for c, t in _PLACEMENT_FIELDS.items()), strict=True)
    grouped: dict[tuple[str, str], set[tuple[Any, ...]]] = {}
    for key, pl in zip(zip(sel["gene"], sel["peptide"], strict=True), fields, strict=True):
        grouped.setdefault(key, set()).add(pl)
    sel = sel.drop_duplicates(["gene", "peptide"]).reset_index(drop=True)
    sel[PLACEMENTS_COLUMN] = [
        tuple(sorted(grouped[k])) for k in zip(sel["gene"], sel["peptide"], strict=True)
    ]
    return sel[[*ACCEPTED_COLUMNS, PLACEMENTS_COLUMN]].astype(_DTYPES)


def _placement_rows(accepted: pd.DataFrame) -> pd.DataFrame:
    """One row per accepted placement (all of a pair's placements); else the accepted rows."""
    if PLACEMENTS_COLUMN not in accepted.columns or accepted.empty:
        return accepted
    ex = accepted[PLACEMENTS_COLUMN].explode().dropna()
    if ex.empty:
        return accepted
    fields = list(_PLACEMENT_FIELDS)
    rest = [c for c in accepted.columns if c not in (*fields, PLACEMENTS_COLUMN)]
    out = accepted.loc[ex.index, rest].reset_index(drop=True)
    parts = pd.DataFrame(ex.tolist(), columns=fields).astype(_PLACEMENT_FIELDS)
    return pd.concat([out, parts], axis=1)


def _wrap(
    owner: np.ndarray, start: np.ndarray, end: np.ndarray, length: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split intervals at the origin into parts inside ``[0, length)`` (vectorised §6.5)."""
    s = np.mod(start, length)
    e = s + (end - start)
    over = e > length
    own = np.concatenate([owner, owner[over]])
    ps = np.concatenate([s, np.zeros(int(over.sum()), dtype=np.int64)])
    pe = np.concatenate([np.minimum(e, length), (e - length)[over]])
    return own, ps, pe


def _coverage_bp(genes: pd.DataFrame, accepted: pd.DataFrame, length: pd.Series) -> np.ndarray:
    """Covered bp per gene: |union of accepted placements ∩ gene| on circular coordinates."""
    n = len(genes)
    covered = np.zeros(n, dtype=np.int64)
    if accepted.empty or n == 0:
        return covered
    pos = pd.Series(np.arange(n), index=genes["locus_tag"].astype(str))
    gi = accepted["gene"].map(pos)
    keep = gi.notna().to_numpy()
    gi_arr = gi[keep].to_numpy(dtype=np.int64)
    gcontig = genes["contig"].astype(str).to_numpy()
    same = accepted["contig"].astype(str).to_numpy()[keep] == gcontig[gi_arr]
    gi_arr = gi_arr[same]
    acc = accepted[keep][same]
    glen = genes["contig"].astype(str).map(length).to_numpy(dtype=np.int64)
    a_own, a_s, a_e = _wrap(
        gi_arr,
        acc["start"].to_numpy(np.int64),
        acc["end"].to_numpy(np.int64),
        glen[gi_arr],
    )
    g_own, g_s, g_e = _wrap(
        np.arange(n),
        genes["start"].to_numpy(np.int64),
        genes["end"].to_numpy(np.int64),
        glen,
    )
    parts = pd.DataFrame({"g": a_own, "s": a_s, "e": a_e}).merge(
        pd.DataFrame({"g": g_own, "gs": g_s, "ge": g_e}), on="g"
    )
    parts["s"] = np.maximum(parts["s"], parts["gs"])
    parts["e"] = np.minimum(parts["e"], parts["ge"])
    parts = parts[parts["e"] > parts["s"]].sort_values(["g", "s", "e"], kind="mergesort")
    if parts.empty:
        return covered
    reach = parts.groupby("g")["e"].cummax()
    prev = reach.groupby(parts["g"]).shift(1).fillna(-1).to_numpy(np.int64)
    s, e = parts["s"].to_numpy(np.int64), parts["e"].to_numpy(np.int64)
    contrib = np.maximum(0, e - np.maximum(s, prev))
    np.add.at(covered, parts["g"].to_numpy(np.int64), contrib)
    return covered


def _distinct(acc: pd.DataFrame, mask: pd.Series, pos: pd.Series, n: int) -> np.ndarray:
    """Distinct peptides per gene among rows of ``acc`` selected by ``mask``."""
    sub = acc.loc[mask, ["gene", "peptide"]].drop_duplicates()
    counts = sub["gene"].map(pos).dropna().astype(np.int64)
    return np.bincount(counts.to_numpy(), minlength=n)[:n]


def gene_evidence(
    genes: pd.DataFrame,
    peptide_hits: pd.DataFrame,
    accepted: pd.DataFrame,
    replicons: pd.DataFrame,
    min_peptides: int = 2,
) -> pd.DataFrame:
    """One row per C2 gene in C2 order: §9.3 evidence columns, ``n_accepted``, flags (§9.4)."""
    n = len(genes)
    tags = genes["locus_tag"].astype(str)
    pos = pd.Series(np.arange(n), index=tags.to_numpy())
    acc = accepted.astype({"gene": str, "peptide": str}) if len(accepted) else accepted
    acc = _placement_rows(acc)
    true = pd.Series(True, index=acc.index)
    n_acc = _distinct(acc, true, pos, n)
    n_uni = _distinct(acc, acc["uniqueness"] == "unique", pos, n)
    n_multi = _distinct(acc, acc["uniqueness"] == "multi_identical", pos, n)
    n_tol = _distinct(acc, acc["route"] == "tolerant", pos, n)
    hits = peptide_hits[peptide_hits["gene"].notna()] if len(peptide_hits) else peptide_hits
    hits = hits.astype({"gene": str, "peptide": str}) if len(hits) else hits
    n_sh = _distinct(hits, hits["uniqueness"] == "shared", pos, n)
    nterm = _distinct(acc, acc["locus_class"] == "nterm_extension", pos, n) > 0
    cont = _distinct(acc, acc["is_contaminant"].astype(bool), pos, n) > 0
    best = np.full(n, np.nan)
    if len(acc):
        mins = acc.groupby("gene")["q_value"].min()
        best = mins.reindex(tags).to_numpy(dtype="float64")
    length = replicons.set_index(replicons["contig"].astype(str))["length"].astype("int64")
    glen = (genes["end"].astype("int64") - genes["start"].astype("int64")).to_numpy()
    cov = _coverage_bp(genes, acc, length) / np.where(glen > 0, glen, 1)
    flags: list[list[str]] = [[] for _ in range(n)]
    for mask, names in (
        (n_multi > 0, ("multi_identical",)),
        (n_sh > 0, ("shared_evidence",)),
        (n_tol > 0, ("tolerant_evidence", "strain_difference")),
        (nterm, ("nterm_extension_evidence",)),
    ):
        for i in np.flatnonzero(mask):
            flags[i].extend(names)
    partial = genes["partial"].astype(bool).to_numpy() if n else np.zeros(0, dtype=bool)
    return pd.DataFrame(
        {
            "locus_tag": tags.to_numpy(dtype=object),
            "n_peptides_unique": n_uni,
            "n_peptides_multi_identical": n_multi,
            "n_peptides_shared": n_sh,
            "n_peptides_tolerant": n_tol,
            "coverage_frac": cov.astype("float64"),
            "best_peptide_q": best,
            "detected_2pep": ~partial & (n_acc >= min_peptides),
            "n_accepted": n_acc,
            "evidence_flags": pd.Series(flags, dtype=object),
            "contaminant_evidence": cont,
        },
        columns=list(EVIDENCE_COLUMNS),
    )
