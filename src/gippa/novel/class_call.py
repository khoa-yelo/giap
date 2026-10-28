"""Novel-ORF class FDR per ORF, its class summaries and the TSV cell (issue M7-11b;
GIPPA_SPEC.md §7.6.4, §8.1, §9.3; docs/specs/M7-11b-class-fdr-wiring.md rules 1-5 and 8).

Pure: wires ``gippa.novel.fdr`` to the step's frames. The picked regime pairs ORFs with the IDs
the search reported (``makedb_id_map``: ``nterm_extension`` locus tags from the run's own
``nterm_entries``); the de novo regime reads E and π̂ from ``placement_costs`` (de novo peptides
only, Q11). No cascade between regimes (Q6).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import pandas as pd

from gippa.config.model import GippaConfig
from gippa.makedb.build import ID_MAP_COLUMNS, nterm_entries
from gippa.model.genes import GeneModel
from gippa.novel import fdr as F
from gippa.novel.candidates import Orf
from gippa.novel.fdr import Classifier, Regime

Reason = Literal["no_decoys", "decoys_not_makedb", "untagged_class", "no_makedb_pair",
                 "paired_share_low", "decoys_unscored", "no_target_score", "lost_pick",
                 "no_target_hit", "denovo_inputs_missing"]  # fmt: skip
_OFF_LIST = frozenset({"untagged_class", "no_makedb_pair", "no_target_score"})


@dataclass(frozen=True)
class OrfClassFdr:
    orf_id: str
    regime: Regime
    fdr_class: str  # the list it competed in: the entry tag (picked), else the ORF class
    q: float  # NaN unless estimated and won (picked) / with a target hit (de novo)
    reason: str | None  # a ``Reason``; None iff q is finite
    class_mismatch: bool  # M7-11 decision 9


@dataclass(frozen=True)
class ClassSummary:
    regime: Literal["picked_decoy", "denovo_subgroup"]
    cls: str
    estimated: bool
    reason: str | None
    n_targets: int  # picked: ORFs entering the list; de novo: target peptides
    n_target_wins: int | None  # picked only
    n_decoys: int  # picked: decoy wins; de novo: decoy peptides of the class
    paired_share: float | None  # picked only
    c: float | None  # de novo only (ADR-0010)
    q_floor: float | None  # picked 1 / n_target_wins, de novo c / n_targets (rule 8)
    floor_above_alpha: bool


@dataclass(frozen=True)
class DenovoInputs:
    target_hits: pd.DataFrame  # peptide, e_target (de novo peptides with a finite e_target)
    decoy_hits: pd.DataFrame  # decoy HIT_COLUMNS rows of de novo peptides + e_decoy
    pi: Mapping[str, float]  # de novo peptide -> pi_denovo (NaN dropped: missing = 1)
    classifier: Classifier


@dataclass(frozen=True)
class ClassFdrResult:
    orfs: dict[str, OrfClassFdr]
    classes: tuple[ClassSummary, ...]  # sorted by (regime, cls)
    n_decoy_rows: int
    decoys_not_makedb: bool


def makedb_id_map(
    peptides: pd.DataFrame, gm: GeneModel, genome: Mapping[str, str],
    decoy_prefixes: Sequence[str],
) -> pd.DataFrame:  # fmt: skip
    """``ID_MAP_COLUMNS`` rows (1-based closed, ``kind = orf``) for every distinct class-tagged
    ID in ``protein_ids`` (decoys stripped); ``nterm_extension`` rows get the locus tag of the
    ``nterm_entries`` record with that ID (NA without one, decision 7)."""
    ids = {F._strip(p, decoy_prefixes) for c in peptides["protein_ids"] for p in F._ids(c)}
    parsed = {i: v for i in sorted(ids) if (v := F.parse_orf_id(i)) is not None}
    tags: dict[str, str | None] = {}
    if any(v[0] == "nterm_extension" for v in parsed.values()):
        tags = {e.db_id: e.locus_tag for e in nterm_entries(gm, genome)}
    rows = [(i, "orf", tags.get(i) if c == "nterm_extension" else None, c, ctg, s + 1, e, st)
            for i, (c, ctg, s, e, st) in parsed.items()]  # fmt: skip
    return pd.DataFrame(rows, columns=list(ID_MAP_COLUMNS)).astype({"start": "int64",
                                                                    "end": "int64"})  # fmt: skip


def denovo_inputs(costs: pd.DataFrame, hits: pd.DataFrame, classifier: Classifier) -> DenovoInputs:
    """The de novo lists' inputs from ``placement_costs`` and the hit tables, restricted to
    ``source = denovo`` peptides (Q11); decoy rows get their peptide's ``e_decoy``."""
    dn = costs[costs["source"].astype(str).to_numpy() == "denovo"]
    et = dn["e_target"].astype("float64")
    th = dn.loc[et.notna().to_numpy(), ["peptide", "e_target"]].reset_index(drop=True)
    is_d = (hits["search"].astype(str) == "decoy") & hits["peptide"].isin(set(dn["peptide"]))
    dh = hits[is_d.to_numpy()].drop(columns=["e_decoy", "e_target"], errors="ignore")
    dh = dh.merge(dn[["peptide", "e_decoy"]], on="peptide", how="left")
    pi = {str(p): float(v) for p, v in zip(dn["peptide"], dn["pi_denovo"], strict=True)
          if pd.notna(v)}  # fmt: skip
    return DenovoInputs(th, dh, pi, classifier)


def _floor(q_floor: float | None, alpha: float) -> tuple[float | None, bool]:
    return q_floor, q_floor is not None and q_floor > alpha


def class_fdr(
    orfs: Sequence[Orf], peptides: pd.DataFrame, *, gm: GeneModel, genome: Mapping[str, str],
    cfg: GippaConfig, denovo: DenovoInputs | None,
) -> ClassFdrResult:  # fmt: skip
    """Class FDR of every ORF in its regime (rules 1-4) and one summary per class list (rule 8);
    ``denovo=None`` makes the de novo ORFs ``denovo_inputs_missing``."""
    prefixes = tuple(cfg.input.decoy_prefixes)
    alpha = cfg.stats.alpha.novel_orfs
    n_dec = int(peptides["is_decoy"].astype(bool).sum())
    makedb = bool(F._decoy_tags(peptides, prefixes))
    regimes = F.orf_regimes(peptides, orfs, prefixes)
    out: dict[str, OrfClassFdr] = {}
    summaries: list[ClassSummary] = []
    why: str = "decoys_not_makedb" if n_dec else "no_decoys"
    for o in orfs:
        if regimes[o.orf_id] == "unestimated":
            out[o.orf_id] = OrfClassFdr(o.orf_id, "unestimated", o.cls, math.nan, why, False)
    picked = [o for o in orfs if regimes[o.orf_id] == "picked_decoy"]
    if picked:
        id_map = makedb_id_map(peptides, gm, genome, prefixes)
        lengths = {k: r.length for k, r in gm.replicons.items()}
        pr = F.picked_class_fdr(picked, peptides, id_map, decoy_prefixes=prefixes,
                                replicon_lengths=lengths)  # fmt: skip
        entries = F._entries(id_map, ())
        by_key = {(e.contig, e.strand, F._stop(e, lengths)): e.cls for e in entries.values()
                  if e.cls != "nterm_extension"}  # fmt: skip
        for o in picked:
            cls = by_key.get(F._orf_key(o, lengths), o.cls) if o.orf_id in pr.class_mismatch \
                else o.cls  # fmt: skip
            q = pr.q[o.orf_id]
            out[o.orf_id] = OrfClassFdr(o.orf_id, "picked_decoy", cls, q,
                                        None if math.isfinite(q) else pr.reason[o.orf_id],
                                        o.orf_id in pr.class_mismatch)  # fmt: skip
        for cls, st in sorted(pr.classes.items()):
            n_t = sum(1 for r in out.values() if r.regime == "picked_decoy"
                      and r.fdr_class == cls and r.reason not in _OFF_LIST)  # fmt: skip
            wins = st.n_target_wins if st.estimated else 0
            qf, above = _floor(1 / wins if wins else None, alpha)
            summaries.append(ClassSummary("picked_decoy", cls, st.estimated, st.reason, n_t,
                                          st.n_target_wins, st.n_decoy_wins, st.paired_share,
                                          None, qf, above))  # fmt: skip
    dn = [o for o in orfs if regimes[o.orf_id] == "denovo_subgroup"]
    if dn and denovo is None:
        for o in dn:
            out[o.orf_id] = OrfClassFdr(o.orf_id, "denovo_subgroup", o.cls, math.nan,
                                        "denovo_inputs_missing", False)  # fmt: skip
    elif dn and denovo is not None:
        d = denovo
        res = F.denovo_class_fdr(dn, d.target_hits, d.decoy_hits, d.classifier, d.pi)
        for o in dn:
            q = res.q[o.orf_id]
            r_dn = None if math.isfinite(q) else "no_target_hit"
            out[o.orf_id] = OrfClassFdr(o.orf_id, "denovo_subgroup", o.cls, q, r_dn, False)
        for cls, ds in sorted(res.classes.items()):
            qf, above = _floor(ds.c / ds.n_targets if ds.n_targets else None, alpha)
            summaries.append(ClassSummary("denovo_subgroup", cls, True, None, ds.n_targets,
                                          None, ds.n_decoys, None, ds.c, qf, above))  # fmt: skip
    summaries.sort(key=lambda s: (s.regime, s.cls))
    return ClassFdrResult(out, tuple(summaries), n_dec, n_dec > 0 and not makedb)


def class_fdr_cell(r: OrfClassFdr) -> float | str:
    """The ``novel_orfs.tsv`` ``class_fdr`` cell (rule 5): q, 1.0 for ``lost_pick`` (Q1), else
    ``"unestimated"``."""
    if math.isfinite(r.q):
        return float(r.q)
    return 1.0 if r.reason == "lost_pick" else "unestimated"


def db_supported(
    orfs: Sequence[Orf], peptides: pd.DataFrame, sources: pd.Series
) -> dict[str, tuple[bool, bool]]:
    """``orf_id`` -> ``(has_denovo, db_support)`` from the ORF's own peptides' resolved C9
    ``sources`` (index-aligned with ``peptides``); decoy and contaminant rows are ignored and a
    folded key takes the union of its rows' sources (M7-17 rules 1-3)."""
    keep = ~(peptides["is_decoy"].astype(bool) | peptides["is_contaminant"].astype(bool))
    src: dict[str, set[str]] = {}
    for p, s in zip(peptides.loc[keep, "peptide"].astype(str), sources[keep], strict=True):
        src.setdefault(F._il(p), set()).add(str(s))
    out: dict[str, tuple[bool, bool]] = {}
    for o in orfs:
        got: set[str] = set().union(*(src.get(p, set()) for p in o.peptides))
        out[o.orf_id] = ("denovo" in got, "database" in got)
    return out
