"""Simulated replicates for the M4-12 ``gene_fdr`` calibration (docs/specs/M4-12-...md).

One replicate: gippasim layout (M4-07, ``n_genes`` chromosome genes) -> expression (M4-08 defaults)
-> evidence (M4-09, ``pi_wrong``, false targets on random genes in proportion to
``n_observable``, matched decoys ``rev_<locus_tag>``, every value emitted) -> M4-10
``covariate_table`` -> M4-11 ``accepted_peptides`` / ``gene_evidence`` -> M4-12 ``detect``.
The input variants only mask what the evidence simulator emits (``gippasim.evidence`` ``emit``
masks values and drops decoy rows the same way):

- ``picked``: everything (decoy rows with ``protein_ids``) -> tier ``decoy_picked``;
- ``mayu_pep``: decoy rows dropped -> ``mayu`` with F(t) = sum PEP;
- ``mayu_nopep``: decoy rows dropped, C9 ``pep`` and C11 ``pep_input`` NA -> ``mayu`` with
  F(t) = t |P(t)|;
- ``two_peptide``: decoys dropped, q / PEP / score NA in C9 and C11 -> ``two_peptide``
  (``unestimated``; recorded, not gated);
- ``loss<pct>``: ``picked`` with ``floor(pct/100 * n)`` of the n decoys that pass
  q <= ``q_value_max`` re-pointed to ``rev_LOST_<k>`` (unpairable; ADR-0009 random decoy loss).

FDP = (detected genes with ``present = false``) / max(1, detected), ``present`` from the M4-09
``gene_truth`` (>= 1 true peptide emitted).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

LEVELS = (0.01, 0.05, 0.10)
N_GENES_FAST, N_GENES_FULL = 2000, 4500
SCORED = frozenset({"q_value", "pep", "score", "is_decoy", "intensity", "protein_ids"})
TIER_OF = {"picked": "decoy_picked", "mayu_pep": "mayu", "mayu_nopep": "mayu",
           "two_peptide": "two_peptide"}  # fmt: skip


def simulate(rng: np.random.Generator, pi: float, n_genes: int) -> tuple[Any, Any]:
    import gippasim.evidence as ev
    import gippasim.expression as ex
    import gippasim.layout as la

    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=n_genes), rng)
    expr = ex.simulate_expression(lay, ex.ExpressionConfig(), rng)
    e = ev.simulate_evidence(lay, expr, ev.EvidenceConfig(pi_wrong=pi, emit=SCORED), rng)
    return lay, e


def masked(e: Any, variant: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(C9, C11) of one input variant (copies; the replicate is not mutated)."""
    p = e.peptides.copy()
    h = e.peptide_hits.copy()
    if variant == "picked" or variant.startswith("loss"):
        return p, h
    p = p[~p["is_decoy"].astype(bool)].reset_index(drop=True)
    if variant in ("mayu_nopep", "two_peptide"):
        p["pep"] = np.nan
        h["pep_input"] = np.nan
    if variant == "two_peptide":
        p["q_value"] = np.nan
        p["score"] = np.nan
        h["best_q_input"] = np.nan
        h["score_input"] = np.nan
    return p, h


def lose_decoys(peptides: pd.DataFrame, share: float, q_max: float,
                rng: np.random.Generator) -> tuple[pd.DataFrame, int]:  # fmt: skip
    """Re-point floor(share n) of the n decoys with q <= q_max to unpairable IDs."""
    p = peptides.copy()
    cand = np.flatnonzero((p["is_decoy"].astype(bool) & (p["q_value"] <= q_max)).to_numpy())
    k = math.floor(share * len(cand))
    lost = rng.choice(cand, size=k, replace=False) if k else np.zeros(0, dtype=int)
    ids = p["protein_ids"].tolist()
    for j, i in enumerate(lost.tolist()):
        ids[i] = [f"rev_LOST_{j}"]
    p["protein_ids"] = pd.Series(ids, index=p.index, dtype=object)
    return p, k


def protein_hits(e: Any) -> pd.DataFrame:
    """C12 of the replicate: an empty contract frame when the simulator wrote none (as
    ``load_m4_inputs`` returns for the header-only ``protein_hits.tsv``)."""
    import abundance_oracle as ao

    return e.protein_hits if e.protein_hits is not None else ao.protein_hits_frame([])


def fdp(table: pd.DataFrame, present: pd.Series) -> float:
    det = table["detected"].astype(bool).to_numpy()
    pres = present.reindex(table["locus_tag"].astype(str)).astype(bool).to_numpy()
    return float((det & ~pres).sum()) / max(1, int(det.sum()))


def run_variant(
    lay: Any,
    e: Any,
    peptides: pd.DataFrame,
    hits: pd.DataFrame,
    levels: Sequence[float] = LEVELS,
    **cfg_kw: Any,
) -> tuple[dict[float, float], set[str]]:
    """Detection rerun at each level (``stats.alpha.identification`` = q); returns the FDP per
    level and the tiers chosen."""
    from gippa.config import GippaConfig
    from gippa.genetable.covariates import covariate_table
    from gippa.genetable.detection import detect
    from gippa.genetable.evidence import accepted_peptides, gene_evidence

    present = e.gene_truth.set_index(e.gene_truth["locus_tag"].astype(str))["present"]
    base = GippaConfig(**cfg_kw)
    acc = accepted_peptides(hits, peptides, base)
    ev = gene_evidence(lay.genes, hits, acc, lay.replicons,
                       min_peptides=base.detection.min_peptides_sequences_only)  # fmt: skip
    cov = covariate_table(lay.genes, base)
    out: dict[float, float] = {}
    tiers: set[str] = set()
    for q in levels:
        stats = {**cfg_kw.get("stats", {}), "alpha": {"identification": q}}
        c = GippaConfig(**{**cfg_kw, "stats": stats})
        d = detect(lay.genes, ev, acc, peptides, hits, protein_hits(e), cov, c)
        out[q] = fdp(d.table, present)
        tiers.add(d.tier)
    return out, tiers


def rep_index(stream: str, rng: np.random.Generator, R: int) -> int:
    """Replicate number of a harness rng (``replicate_rng(stream, i)``, i < 2R: both blocks)."""
    import _criteria as C

    entropy = int(rng.bit_generator.seed_seq.entropy)  # type: ignore[attr-defined]
    for i in range(2 * R):
        if C.suite_seed(stream, i) == entropy:
            return i
    raise LookupError(f"rng is not a replicate of {stream!r} (R = {R})")
