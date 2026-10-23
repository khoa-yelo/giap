"""Detection tiers and gene-level FDR (docs/specs/M4-12-detection-tiers.md).

Spec §2.1 decision 4, §7.3 (eligible = non-partial genes), §7.3.2 (tier table, picked steps 1-4,
MAYU), §9.3 (``detection_status``; ``detection_summary.tsv`` section ``tier``). ADR-0009: picked
FDR(t) = (D(t)/f + 1) / T(t) with f the paired share of the decoys passing q <= ``q_value_max``.
Lead decisions: the pairing cut-off is max(``mapping.proteins.min_identity``, ``min_coverage``);
contaminant evidence counts toward detection in every tier; in the ``mayu`` tier a warning
``unspecific_search`` fires when more than half of the accepted non-C-terminal peptides do not end
in K/R. Review decisions 5-8: PEP-only input gets peptide q from the running mean of the sorted
PEPs; all-contaminant decoys are neither paired nor unpaired (decision 12); the longest decoy
prefix is stripped first; ``Detection.protein_gene_map`` is the map ``detect`` used. Everything
here is pure.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.stats.mayu import mayu_gene_fdr
from giap.stats.tda import picked_fdr_rescaled

if TYPE_CHECKING:
    from giap.config import GiapConfig

__all__ = [
    "CURVE_THRESHOLDS",
    "MIN_DECOY_PAIRING",
    "Detection",
    "Tier",
    "choose_tier",
    "detect",
    "pair_decoys",
    "protein_gene_map",
]

Tier = Literal["decoy_picked", "mayu", "two_peptide", "protein_list"]
MIN_DECOY_PAIRING = 0.80
CURVE_THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.1)
UNSPECIFIC_SHARE = 0.5
TABLE_COLUMNS = ("locus_tag", "detected", "detection_status", "gene_q_value", "detection_tier")
SUMMARY_COLUMNS = ("section", "group_type", "group", "metric", "value", "ci_low", "ci_high", "n",
                   "flags")  # fmt: skip


@dataclass(frozen=True, slots=True)
class Detection:
    """Tier, per-gene detection table (C2 order), ``tier`` summary rows, warnings and the
    protein -> gene map (``protein_id``, ``gene``, ``method``) ``detect`` used."""

    tier: Tier
    table: pd.DataFrame
    summary_rows: pd.DataFrame
    warnings: tuple[str, ...]
    protein_gene_map: pd.DataFrame


def _ids(cell: Any) -> list[str]:
    """Members of a C9 ``protein_ids`` cell; protein groups ``A;B`` are split."""
    if cell is None or isinstance(cell, float):
        return []
    items = [cell] if isinstance(cell, str) else list(cell)
    return [m for x in items for m in str(x).split(";") if m]


def _is_na(x: Any) -> bool:
    return x is None or x is pd.NA or (isinstance(x, float) and np.isnan(x))


def _contaminants(prefixes: Sequence[str] | None) -> tuple[str, ...]:
    """Configured contaminant prefixes; ``None`` -> ``input.contaminant_prefixes`` default."""
    if prefixes is None:
        from giap.config.model import InputConfig

        prefixes = InputConfig().contaminant_prefixes
    return tuple(prefixes)


def protein_gene_map(
    genes: pd.DataFrame,
    peptides: pd.DataFrame | None,
    peptide_hits: pd.DataFrame,
    protein_hits: pd.DataFrame,
    *,
    accepted: pd.DataFrame | None = None,
    contaminant_prefixes: Sequence[str] | None = None,
) -> pd.DataFrame:
    """protein_id -> gene by (1) C2 locus tag, (2) C12 rank-1 hit gene, (3) peptide vote.

    The vote counts distinct unique target peptides per gene among the peptides whose C9
    ``protein_ids`` contain the protein (``accepted`` rows when given, else C11 unique placements
    with a gene); a tie leaves the protein unmapped. Contaminant-prefixed IDs are never mapped.
    Only mapped proteins are returned (columns ``protein_id``, ``gene``, ``method``);
    ``contaminant_prefixes`` defaults to ``input.contaminant_prefixes`` of the default config.
    """
    cont = _contaminants(contaminant_prefixes)
    tags = genes["locus_tag"].astype(str).tolist()
    out: dict[str, tuple[str, str]] = {t: (t, "locus_tag") for t in tags if not t.startswith(cont)}
    if len(protein_hits):
        ph = protein_hits[(protein_hits["rank"].astype(int) == 1) & protein_hits["gene"].notna()]
        by_id: dict[str, set[str]] = {}
        pids, hit_genes = ph["protein_id"].astype(str).tolist(), ph["gene"].astype(str).tolist()
        for pid, g in zip(pids, hit_genes, strict=True):
            by_id.setdefault(pid, set()).add(g)
        for pid in sorted(by_id):
            if pid not in out and not pid.startswith(cont) and len(by_id[pid]) == 1:
                out[pid] = (next(iter(by_id[pid])), "protein_hit")
    if peptides is not None and len(peptides):
        if accepted is not None:
            src = accepted[accepted["uniqueness"] == "unique"]
        elif len(peptide_hits):
            src = peptide_hits[
                peptide_hits["gene"].notna() & (peptide_hits["uniqueness"] == "unique")
            ]
        else:
            src = pd.DataFrame({"gene": [], "peptide": []})
        pairs = src[["gene", "peptide"]].astype(str).drop_duplicates()
        tgt = peptides[~peptides["is_decoy"].astype(bool)]
        tpep = tgt["peptide"].astype(str).tolist()
        pid_of = dict(zip(tpep, tgt["protein_ids"].tolist(), strict=True))
        votes: dict[str, dict[str, int]] = {}
        for g, p in zip(pairs["gene"].tolist(), pairs["peptide"].tolist(), strict=True):
            for pid in set(_ids(pid_of.get(p))):
                if pid not in out and not pid.startswith(cont):
                    v = votes.setdefault(pid, {})
                    v[g] = v.get(g, 0) + 1
        for pid in sorted(votes):
            ranked = sorted(votes[pid].items(), key=lambda kv: -kv[1])
            if len(ranked) == 1 or ranked[0][1] > ranked[1][1]:
                out[pid] = (ranked[0][0], "peptide_vote")
    return pd.DataFrame(
        [(k, g, m) for k, (g, m) in out.items()], columns=["protein_id", "gene", "method"]
    ).astype(object)


def _strip(pid: str, prefixes: Sequence[str]) -> str:
    """Remove the longest matching decoy prefix (``REV__X`` -> ``X`` whatever the order)."""
    for pre in sorted(prefixes, key=len, reverse=True):
        if pid.startswith(pre):
            return pid[len(pre) :]
    return pid


def pair_decoys(
    peptides: pd.DataFrame | None,
    pmap: pd.DataFrame,
    decoy_prefixes: Sequence[str],
    *,
    q_value_max: float | None = None,
    contaminant_prefixes: Sequence[str] | None = None,
) -> pd.DataFrame:
    """One row per decoy peptide with ``protein_ids`` (and q <= ``q_value_max`` when both exist):
    ``peptide, q_value, pep, score, gene``; ``gene`` NA when the stripped IDs do not all map to
    one and the same gene (unpaired). Contaminant decoys (every stripped ID starts with a
    contaminant prefix, default ``input.contaminant_prefixes``) get no row: neither paired nor
    unpaired; a mixed decoy pairs on its non-contaminant IDs only (decision 12)."""
    cols = ["peptide", "q_value", "pep", "score", "gene"]
    if peptides is None or not len(peptides):
        return pd.DataFrame({c: pd.Series([], dtype=object) for c in cols})
    d = peptides[peptides["is_decoy"].astype(bool)]
    d = d[np.array([len(_ids(c)) > 0 for c in d["protein_ids"].tolist()], dtype=bool)]
    if q_value_max is not None:
        q = d["q_value"].astype("float64")
        d = d[q.isna() | (q <= q_value_max)]
    cont = _contaminants(contaminant_prefixes)
    # Decision 12: only all-contaminant decoys leave the counts; a mixed decoy pairs on its
    # non-contaminant IDs.
    stripped = [
        [m for m in (_strip(x, decoy_prefixes) for x in _ids(c)) if not m.startswith(cont)]
        for c in d["protein_ids"].tolist()
    ]
    keep = np.array([len(ms) > 0 for ms in stripped], dtype=bool)
    d = d[keep]
    lookup = dict(zip(pmap["protein_id"].astype(str), pmap["gene"], strict=True))
    genes: list[Any] = []
    for ms in (s for s, k in zip(stripped, keep, strict=True) if k):
        hit = [lookup.get(m) for m in ms]
        found = {str(g) for g in hit if not _is_na(g)}
        ok = len(found) == 1 and all(not _is_na(g) for g in hit)
        genes.append(found.pop() if ok else None)
    out = pd.DataFrame({"peptide": np.asarray(d["peptide"].astype(str).tolist(), dtype=object)})
    for c in ("q_value", "pep", "score"):
        out[c] = d[c].astype("float64").to_numpy()
    out["gene"] = pd.Series(genes, dtype=object)
    return out


def choose_tier(
    peptides: pd.DataFrame | None,
    peptide_hits: pd.DataFrame,
    protein_hits: pd.DataFrame,
    pairing: pd.DataFrame,
    *,
    accepted: pd.DataFrame | None = None,
    min_pairing: float = MIN_DECOY_PAIRING,
) -> tuple[Tier, list[str]]:
    """Tier of the §7.3.2 table and warnings (``decoy_pairing_low``).

    q / PEP presence is read from ``accepted`` when given, else from C9 targets (or C11 without C9).
    """
    warns: list[str] = []
    if len(pairing):
        share = float(pairing["gene"].notna().mean())
        if share >= min_pairing:
            return "decoy_picked", warns
        return "mayu", ["decoy_pairing_low"]
    if accepted is not None:
        qp = accepted[["q_value", "pep"]]
    elif peptides is not None:
        qp = peptides.loc[~peptides["is_decoy"].astype(bool), ["q_value", "pep"]]
    else:
        qp = peptide_hits[["best_q_input", "pep_input"]] if len(peptide_hits) else pd.DataFrame()
    if len(qp) and bool(qp.astype("float64").notna().any(axis=None)):
        return "mayu", warns
    if len(peptide_hits):
        return "two_peptide", warns
    return "protein_list", warns


def _key_rank(q: Iterable[float], pep: Iterable[float], score: Iterable[float]) -> npt.NDArray[Any]:
    """Dense rank of (q asc, PEP asc, score desc), NA last; higher = better."""
    a = np.column_stack(
        [np.asarray(list(q), dtype=np.float64), np.asarray(list(pep), dtype=np.float64),
         -np.asarray(list(score), dtype=np.float64)]
    )  # fmt: skip
    a = np.where(np.isnan(a), np.inf, a) + 0.0
    if not len(a):
        return np.zeros(0, dtype=np.float64)
    uniq, inv = np.unique(a, axis=0, return_inverse=True)
    return (len(uniq) - np.ravel(inv)).astype(np.float64)


def _picked_q(
    tags: Sequence[str], accepted: pd.DataFrame, pairing: pd.DataFrame, f: float
) -> npt.NDArray[np.float64]:
    """Picked gene q (ADR-0009: each decoy win counts 1/f); NaN for decoy wins / no target."""
    paired = pairing[pairing["gene"].notna()]
    n_t = len(accepted)
    cols = ("q_value", "pep", "score")
    rank = _key_rank(*(np.r_[accepted[c].astype(float), paired[c].astype(float)] for c in cols))
    best = np.full((2, len(tags)), np.nan)
    for k, (g, r) in enumerate(((accepted["gene"], rank[:n_t]), (paired["gene"], rank[n_t:]))):
        s = pd.Series(r).groupby(np.asarray(g.astype(str).tolist(), dtype=object)).max()
        best[k] = s.reindex(list(tags)).to_numpy(dtype=np.float64)
    return picked_fdr_rescaled(best[0], best[1], f).q


def _pep_q(pep: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Peptide q from PEP (decision 5): mean of the PEPs <= PEP_i over the distinct peptides
    (the running mean the MAYU-with-PEP path sums); NaN where PEP is NaN."""
    out = np.full(pep.size, np.nan)
    ok = ~np.isnan(pep)
    if not ok.any():
        return out
    srt = np.sort(pep[ok])
    run = np.cumsum(srt) / np.arange(1, srt.size + 1)
    out[ok] = run[np.searchsorted(srt, pep[ok], side="right") - 1]
    return out


def _unspecific(genes: pd.DataFrame, accepted: pd.DataFrame) -> bool:
    """More than half of accepted peptides (C-terminal placements excepted) not ending in K/R."""
    if not len(accepted):
        return False
    g = genes.set_index(genes["locus_tag"].astype(str))
    gene = accepted["gene"].astype(str)
    plus = gene.map(g["strand"].astype(str)).to_numpy() == "+"
    g_end, g_start = gene.map(g["end"]).to_numpy(), gene.map(g["start"]).to_numpy()
    cterm = np.where(plus, accepted["end"].to_numpy() >= g_end - 3,
                     accepted["start"].to_numpy() <= g_start + 3)  # fmt: skip
    peps = accepted.loc[~cterm, "peptide"].astype(str).drop_duplicates()
    if not len(peps):
        return False
    return float((~peps.str[-1].isin(["K", "R"])).mean()) > UNSPECIFIC_SHARE


def _row(group_type: str, group: str, metric: str, value: float) -> dict[str, Any]:
    return {"section": "tier", "group_type": group_type, "group": group, "metric": metric,
            "value": float(value), "ci_low": np.nan, "ci_high": np.nan, "n": pd.NA,
            "flags": ""}  # fmt: skip


def _validate(accepted: pd.DataFrame, pairing: pd.DataFrame) -> None:
    pep = np.r_[accepted["pep"].astype(float), pairing["pep"].astype(float)]
    if ((pep < 0) | (pep > 1)).any():
        raise ValueError("pep outside [0, 1] among accepted peptides or decoys")


def detect(
    genes: pd.DataFrame,
    evidence: pd.DataFrame,
    accepted: pd.DataFrame,
    peptides: pd.DataFrame | None,
    peptide_hits: pd.DataFrame,
    protein_hits: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: GiapConfig,
) -> Detection:
    """Detection per C2 gene in the run's tier (§7.3.2) with the ``tier`` summary rows."""
    alpha = cfg.stats.alpha.identification
    tags = genes["locus_tag"].astype(str).tolist()
    n = len(tags)
    eligible = ~genes["partial"].astype(bool).to_numpy() if n else np.zeros(0, dtype=bool)
    ev = evidence.set_index(evidence["locus_tag"].astype(str)).reindex(tags)
    pmap = protein_gene_map(genes, peptides, peptide_hits, protein_hits, accepted=accepted,
                            contaminant_prefixes=cfg.input.contaminant_prefixes)  # fmt: skip
    pairing = pair_decoys(peptides, pmap, cfg.input.decoy_prefixes,
                          q_value_max=cfg.input.q_value_max,
                          contaminant_prefixes=cfg.input.contaminant_prefixes)  # fmt: skip
    min_pair = max(cfg.mapping.proteins.min_identity, cfg.mapping.proteins.min_coverage)
    tier, warns = choose_tier(peptides, peptide_hits, protein_hits, pairing, accepted=accepted,
                              min_pairing=min_pair)  # fmt: skip
    _validate(accepted, pairing)
    if not eligible.any():
        warns.append("no_eligible_genes")
    n_paired = int(pairing["gene"].notna().sum())
    n_unpaired = len(pairing) - n_paired
    f = n_paired / len(pairing) if len(pairing) else np.nan
    n_acc = ev["n_accepted"].fillna(0).astype(int).to_numpy()
    q = np.full(n, np.nan)
    status = np.full(n, "not_detected", dtype=object)
    curve: list[dict[str, Any]] = []
    elig_tags = [t for t, e in zip(tags, eligible, strict=True) if e]
    q_floor = np.nan
    if tier == "decoy_picked" and elig_tags:
        q[eligible] = _picked_q(elig_tags, accepted, pairing, f)
        n_win = int(np.isfinite(q).sum())
        q_floor = 1.0 / n_win if n_win else np.nan  # smallest picked q: no decoy win, all T wins
        if q_floor > alpha:  # Khoa 2026-10-22: flag a too-small genome; q is reported unchanged
            warns.append("q_floor_above_alpha")
    elif tier == "mayu" and elig_tags:
        sizes = covariates.set_index(covariates["locus_tag"].astype(str))["n_observable_peptides"]
        a = sizes.reindex(elig_tags).astype(float).to_numpy()
        if not np.isfinite(a).all():
            raise ValueError("n_observable_peptides must be finite for every eligible gene")
        uniq = accepted.drop_duplicates("peptide")
        pep = uniq["pep"].astype(float).to_numpy()
        if np.isnan(pep).any() and not np.isnan(pep).all():
            warns.append("pep_partial")
        use_pep = pep if pep.size and not np.isnan(pep).any() else None
        pq = uniq["q_value"].astype(float).to_numpy()
        best = ev["best_peptide_q"].astype(float).to_numpy()[eligible]
        if pq.size and np.isnan(pq).all() and not np.isnan(pep).all():
            pq = _pep_q(pep)  # decision 5: PEP-only input
            by_pep = dict(zip(uniq["peptide"].astype(str).tolist(), pq.tolist(), strict=True))
            gq = accepted["peptide"].astype(str).map(by_pep).astype(float)
            gbest = gq.groupby(accepted["gene"].astype(str).to_numpy()).min()
            best = gbest.reindex(elig_tags).to_numpy(dtype=np.float64)
        res = mayu_gene_fdr(best, pq, use_pep, a, cfg.detection.mayu_size_bins)
        q[eligible] = res.gene_q
        for c in CURVE_THRESHOLDS:
            row = res.curve[res.curve["t"] <= c + 1e-12].tail(1)
            for m in ("fdr", "fp", "n_genes", "fdr_2pep"):
                curve.append(_row("q_threshold", repr(c), m,
                                  float(row[m].iloc[0]) if len(row) else np.nan))  # fmt: skip
        if _unspecific(genes, accepted):
            warns.append("unspecific_search")
    if tier in ("decoy_picked", "mayu"):
        status[q <= alpha] = "detected"
        for c in CURVE_THRESHOLDS if elig_tags else ():
            hit = q[q <= c + 1e-12]
            curve.append(_row("q_threshold", repr(c), "n_detected", hit.size))
            if tier == "decoy_picked":
                curve.append(_row("q_threshold", repr(c), "fdr",
                                  hit.max() if hit.size else np.nan))  # fmt: skip
    elif tier == "two_peptide":
        status[n_acc == 1] = "single_peptide_unestimated"
        status[ev["detected_2pep"].fillna(False).astype(bool).to_numpy()] = "detected"
    else:
        ph = protein_hits
        top = ph[(ph["rank"].astype(int) == 1) & ph["gene"].notna()] if len(ph) else ph
        hit_genes = set(top["gene"].astype(str)) if len(top) else set()
        status[np.array([t in hit_genes for t in tags], dtype=bool)] = "detected"
    status[~eligible] = "partial_excluded"
    q[~eligible] = np.nan
    detected = status == "detected"
    table = pd.DataFrame({"locus_tag": np.array(tags, dtype=object), "detected": detected,
                          "detection_status": status, "gene_q_value": q,
                          "detection_tier": np.full(n, tier, dtype=object)},
                         columns=list(TABLE_COLUMNS))  # fmt: skip
    n_2pep = int((ev["detected_2pep"].fillna(False).astype(bool).to_numpy() & eligible).sum())
    tier_rows = [
        _row("tier", tier, "n_detected", int(detected.sum())),
        _row("tier", tier, "n_detected_2pep", n_2pep),
        _row("tier", tier, "n_decoys_paired", n_paired),
        _row("tier", tier, "n_decoys_unpaired", n_unpaired),
        _row("tier", tier, "decoy_paired_share", f),
        _row("tier", tier, "decoy_win_weight", 1.0 / f if f and f > 0 else np.nan),
        _row("tier", tier, "q_floor", q_floor),
    ]
    summary = pd.DataFrame([*tier_rows, *curve], columns=list(SUMMARY_COLUMNS))
    return Detection(
        tier=tier, table=table, summary_rows=summary, warnings=tuple(warns), protein_gene_map=pmap
    )
