"""Inputs and an independent oracle for the M4-16 detection-summary tests.

Gene frames carry only ``genes.tsv`` columns (GIAP_SPEC.md §9.3) that the merged producers write:
``locus_tag``, ``contig``, ``feature_type``, ``hypothetical``, ``partial`` (C2, M4-06/M4-07),
``tu_id`` (C6, giapsim ``layers_frame``), ``detected`` / ``detection_status`` (M4-12 two-peptide
rule applied to M4-11 ``gene_evidence``), ``p_detect`` and ``cog_category``.

Stand-ins for producers that are not on main yet (documented, interface-faithful):
- ``p_detect`` (M4-13, PR #83): the simulator's ``p_detect_true`` (M4-08); NA for partial genes.
- ``cog_category`` (M4-15 ``function_columns``): C3 ``COG`` xrefs mapped through a fixed test
  category table (some IDs map to two letters); NA without a COG xref.
- detectability diagnostics: rows shaped exactly like ``calibration_diagnostics`` on
  ``origin/m4/detectability`` (7 columns, no ``n`` / ``flags``; Wald CIs on slope/intercept).
- tier rows: the M4-12 ``summary_rows`` shape (section ``tier``; ``tier`` and ``q_threshold``
  groups), same 7 columns.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

SUMMARY_COLUMNS = ("section", "group_type", "group", "metric", "value", "ci_low", "ci_high", "n",
                   "flags")  # fmt: skip
PART_COLUMNS = ("section", "group_type", "group", "metric", "value", "ci_low", "ci_high")
SECTIONS = ("tier", "basic", "detectability", "qc")
BASIC_METRICS = ("n_genes", "n_detected", "pct_detected", "expected", "o_minus_e", "o_over_e")
FEATURE_CLASSES = ("cds", "sorf", "pseudogene", "hypothetical")
FEW_UNITS_FLAG = "ci_few_units"  # docs/methods/bootstrap.md: CIs from < 10 TUs
FEW_UNITS_MIN = 10
GENE_COLUMNS = ("locus_tag", "contig", "feature_type", "hypothetical", "partial", "tu_id",
                "detected", "detection_status", "p_detect", "cog_category")  # fmt: skip


# --------------------------------------------------------------------------- hand-built genes


def gene(
    locus_tag: str,
    tu_id: str,
    detected: bool,
    p_detect: float | None,
    *,
    contig: str = "chr1",
    feature_type: str = "cds",
    hypothetical: bool = False,
    cog_category: str | None = None,
    status: str | None = None,
) -> dict[str, Any]:
    """One genes.tsv row (summary columns only); ``status`` defaults from ``detected``."""
    st = status or ("detected" if detected else "not_detected")
    return {
        "locus_tag": locus_tag,
        "contig": contig,
        "feature_type": feature_type,
        "hypothetical": hypothetical,
        "partial": st == "partial_excluded",
        "tu_id": tu_id,
        "detected": bool(detected),
        "detection_status": st,
        "p_detect": math.nan if p_detect is None else float(p_detect),
        "cog_category": cog_category,
    }


def genes_frame(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame([dict(r) for r in rows], columns=list(GENE_COLUMNS))
    return df.astype(
        {"hypothetical": bool, "partial": bool, "detected": bool, "p_detect": "float64"}
    )


# --------------------------------------------------------------------------- oracle


def is_na(v: Any) -> bool:
    if v is None:
        return True
    try:
        return bool(pd.isna(v))
    except (TypeError, ValueError):
        return False


def cog_letters(v: Any) -> list[str]:
    if is_na(v) or str(v) == "":
        return ["none"]
    return sorted(set(str(v)))


def group_members(genes: pd.DataFrame) -> dict[tuple[str, str], np.ndarray]:
    """(group_type, group) -> boolean mask over ``genes`` rows (eligible genes only)."""
    elig = (genes["detection_status"].astype(str) != "partial_excluded").to_numpy()
    out: dict[tuple[str, str], np.ndarray] = {}

    def add(key: tuple[str, str], mask: np.ndarray) -> None:
        m = mask & elig
        if m.any():
            out[key] = m

    contig = genes["contig"].astype(str).to_numpy()
    for c in sorted(set(contig[elig])):
        add(("replicon", c), contig == c)
    ft = genes["feature_type"].astype(str).to_numpy()
    for f in ("cds", "sorf", "pseudogene"):
        add(("feature_class", f), ft == f)
    add(("feature_class", "hypothetical"), genes["hypothetical"].astype(bool).to_numpy())
    letters = [cog_letters(v) for v in genes["cog_category"]]
    for letter in sorted({x for ls in letters for x in ls}):
        add(("cog_category", letter), np.array([letter in ls for ls in letters], dtype=bool))
    return out


def point_estimates(genes: pd.DataFrame, mask: np.ndarray) -> dict[str, float]:
    d = genes["detected"].astype(bool).to_numpy()[mask]
    p = genes["p_detect"].astype(float).to_numpy()[mask]
    n, o, e = int(mask.sum()), float(d.sum()), float(p.sum())
    return {
        "n_genes": float(n),
        "n_detected": o,
        "pct_detected": 100.0 * o / n,
        "expected": e,
        "o_minus_e": o - e,
        "o_over_e": o / e if e > 0 else math.nan,
        "n_tu": float(len(set(genes["tu_id"].astype(str).to_numpy()[mask]))),
    }


def tu_sums(genes: pd.DataFrame, mask: np.ndarray) -> dict[str, float]:
    """Per TU: sum of d - p over the group's genes in that TU."""
    sub = genes[mask]
    r = sub["detected"].astype(float) - sub["p_detect"].astype(float)
    return {str(k): float(v) for k, v in r.groupby(sub["tu_id"].astype(str)).sum().items()}


def flag_set(v: Any) -> set[str]:
    if isinstance(v, (list, tuple, set, np.ndarray)):
        return {str(x) for x in v if str(x)}
    if is_na(v):
        return set()
    return {x for x in str(v).split(";") if x}


def rows_of(df: pd.DataFrame, group_type: str, group: str) -> pd.DataFrame:
    sel = (df["group_type"].astype(str) == group_type) & (df["group"].astype(str) == group)
    return df[sel]


def metric(df: pd.DataFrame, group_type: str, group: str, name: str) -> pd.Series:
    r = rows_of(df, group_type, group)
    r = r[r["metric"].astype(str) == name]
    assert len(r) == 1, (group_type, group, name, len(r))
    return r.iloc[0]


# --------------------------------------------------------------------------- QC inputs (C9/C11/C13)


def qc_peptide(peptide: str, *, search_db: str = "same", is_decoy: bool = False) -> dict[str, Any]:
    import evidence_oracle as eo

    row = eo.pep(peptide, 0.001, is_decoy=is_decoy)
    row["search_db"] = search_db
    return row


def qc_hit(peptide: str, *, route: str = "exact", flags: Sequence[str] = ()) -> dict[str, Any]:
    import evidence_oracle as eo

    return eo.hit(peptide, "GSIM_00005", 1000, route=route, flags=list(flags),
                  n_diff=1 if route == "tolerant" else None)  # fmt: skip


def peptides_frame(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    import evidence_oracle as eo

    return eo.peptides_frame(rows)  # type: ignore[no-any-return]


def hits_frame(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    import evidence_oracle as eo

    return eo.hits_frame(rows)  # type: ignore[no-any-return]


UNMAPPED_COLUMNS = ("sequence_id", "sequence", "kind", "label", "reason", "contaminant_match")


def unmapped_frame(rows: Iterable[tuple[str, str, str]]) -> pd.DataFrame:
    """C13 rows from (sequence_id, kind, label)."""
    recs = []
    for sid, kind, label in rows:
        reason = "contaminant_id" if label == "contaminant" else "no_hit"
        match = "CONT_P1" if label == "contaminant" else None
        recs.append((sid, "PEPTIDEK" if kind == "peptide" else "MPROTEIN", kind, label, reason,
                     match))  # fmt: skip
    return pd.DataFrame(recs, columns=list(UNMAPPED_COLUMNS))


# --------------------------------------------------------------------------- other summary parts


def part_row(section: str, group_type: str, group: str, metric_: str, value: float,
             ci: tuple[float, float] = (math.nan, math.nan)) -> dict[str, Any]:  # fmt: skip
    return {"section": section, "group_type": group_type, "group": str(group), "metric": metric_,
            "value": float(value), "ci_low": float(ci[0]), "ci_high": float(ci[1])}  # fmt: skip


def detectability_rows(d: Sequence[bool], p: Sequence[float], n_bins: int = 10) -> pd.DataFrame:
    """Rows shaped like M4-13 ``calibration_diagnostics`` (origin/m4/detectability)."""
    dd = np.asarray(d, dtype=float)
    pp = np.asarray(p, dtype=float)
    rows = []
    order = np.argsort(pp, kind="stable")
    for i, b in enumerate(np.array_split(order, n_bins), start=1):
        rows.append(part_row("detectability", "bin", str(i), "mean_p",
                             float(pp[b].mean()) if b.size else math.nan))  # fmt: skip
        rows.append(part_row("detectability", "bin", str(i), "observed_rate",
                             float(dd[b].mean()) if b.size else math.nan))  # fmt: skip
        rows.append(part_row("detectability", "bin", str(i), "n", float(b.size)))
    rows.append(part_row("detectability", "overall", "all", "brier",
                         float(np.mean((pp - dd) ** 2))))  # fmt: skip
    rows.append(part_row("detectability", "overall", "all", "calibration_intercept", 0.05,
                         (-0.2, 0.3)))  # fmt: skip
    rows.append(part_row("detectability", "overall", "all", "calibration_slope", 0.97,
                         (0.81, 1.13)))  # fmt: skip
    rows.append(part_row("detectability", "overall", "all", "n_detected", float(dd.sum())))
    rows.append(part_row("detectability", "overall", "all", "n_eligible", float(dd.size)))
    return pd.DataFrame(rows, columns=list(PART_COLUMNS))


def tier_rows(n_detected: int, n_detected_2pep: int) -> pd.DataFrame:
    """Rows shaped like M4-12 ``summary_rows`` (section ``tier``)."""
    rows = [
        part_row("tier", "tier", "two_peptide", "n_detected", n_detected),
        part_row("tier", "tier", "two_peptide", "n_detected_2pep", n_detected_2pep),
    ]
    for q in ("0.001", "0.005", "0.01", "0.02", "0.05", "0.1"):
        rows.append(part_row("tier", "q_threshold", q, "n_detected", n_detected))
        rows.append(part_row("tier", "q_threshold", q, "fdr", math.nan))
    return pd.DataFrame(rows, columns=list(PART_COLUMNS))


# --------------------------------------------------------------------------- simulated fixture

# Fixed test COG category table (stands in for the COG 2020 table of M4-15 / lead decision 2).
_LETTERS = "JKLDVTMNUOCGEFHIPQS"


def cog_category_of(cog_id: str) -> str:
    k = int(cog_id[3:])
    first = _LETTERS[k % len(_LETTERS)]
    if k % 7 == 0:  # some COGs carry two categories (e.g. "KT")
        second = _LETTERS[(k // 7) % len(_LETTERS)]
        return "".join(sorted({first, second}))
    return first


def cog_categories(genes: pd.DataFrame, xrefs: pd.DataFrame) -> pd.Series:
    cog = xrefs[xrefs["db"].astype(str) == "COG"].drop_duplicates("locus_tag")
    m = dict(zip(cog["locus_tag"].astype(str), cog["id"].astype(str), strict=True))
    vals = [cog_category_of(m[t]) if t in m else None for t in genes["locus_tag"].astype(str)]
    return pd.Series(vals, index=genes.index, dtype=object)


def two_peptide_status(evidence: pd.DataFrame, partial: pd.Series) -> tuple[pd.Series, pd.Series]:
    """M4-12 two_peptide tier on M4-11 ``gene_evidence``: (detected, detection_status)."""
    det2 = evidence["detected_2pep"].astype(bool).to_numpy()
    n_acc = evidence["n_accepted"].astype(int).to_numpy()
    part = partial.astype(bool).to_numpy()
    single = np.where(n_acc == 1, "single_peptide_unestimated", "not_detected")
    status = np.where(det2, "detected", single)
    status = np.where(part, "partial_excluded", status)
    detected = det2 & ~part
    return (pd.Series(detected, index=evidence.index, dtype=bool),
            pd.Series(status, index=evidence.index, dtype=object))  # fmt: skip


@dataclass(frozen=True)
class SimFixture:
    inputs: Any  # giap.io.work_tables.M4Inputs (loaded through load_m4_inputs)
    genes: pd.DataFrame  # genes.tsv summary columns
    detectability: pd.DataFrame
    tier: pd.DataFrame


def _cfg() -> Any:
    from giap.config import GiapConfig

    return GiapConfig(genome_relation="same", input={"q_value_max": 0.01},
                      stats={"alpha": {"placement": 0.01}})  # fmt: skip


def sim_fixture(work: Path, seed: int = 7, n_genes: int = 300) -> SimFixture:
    """M4-07 layout + M4-08 expression + perfect C9/C11 written to ``work`` and read back with
    M4-06 ``load_m4_inputs``; M4-11 evidence gives the two-peptide detection columns."""
    import evidence_oracle as eo

    from giap.genetable.evidence import accepted_peptides, gene_evidence
    from giap.io.work_tables import default_paths, load_m4_inputs
    from giapsim.contract_tables import layers_frame, write_contract_tables
    from giapsim.expression import ExpressionConfig, simulate_expression
    from giapsim.layout import LayoutConfig, simulate_layout

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "unit" / "genetable"))
    import m4_contract_fixture as mf

    layout = simulate_layout(
        LayoutConfig(n_genes_chromosome=n_genes, plasmids=(25,), linear_contigs=(15,)),
        np.random.default_rng(seed),
    )
    expr = simulate_expression(layout, ExpressionConfig(), np.random.default_rng(seed + 1))
    detected_true = dict(zip(expr.genes["locus_tag"].astype(str),
                             expr.genes["detected_true"].astype(bool), strict=True))  # fmt: skip
    peps, hits = eo.perfect_evidence(layout.genes, layout.replicons, detected_true,
                                     np.random.default_rng(seed + 2))  # fmt: skip
    peps = peps.copy()
    peps["search_db"] = "same"
    hits = hits.copy()
    # two engine/protein mismatches and one same-db peptide without any placement
    first = hits["peptide"].drop_duplicates().iloc[:2]
    hits["flags"] = [list(f) + (["engine_protein_mismatch"] if p in set(first) else [])
                     for f, p in zip(hits["flags"], hits["peptide"], strict=True)]  # fmt: skip
    peps = pd.concat([peps, eo.peptides_frame([qc_peptide("NOWHEREPEPK")])], ignore_index=True)
    unmapped = unmapped_frame([("NOWHEREPEPK", "peptide", "unmapped"),
                               ("cont_1", "peptide", "contaminant"),
                               ("cont_2", "protein", "contaminant")])  # fmt: skip
    protein_hits = pd.DataFrame({c: pd.Series([], dtype=object)
                                 for c in mf.PROTEIN_HITS_COLUMNS})  # fmt: skip

    write_contract_tables(layout, work)
    paths = default_paths(work)
    frames = {"peptides": peps, "peptide_hits": hits, "protein_hits": protein_hits,
              "unmapped": unmapped}  # fmt: skip
    mf.write_work(frames, {k: v for k, v in paths.items() if k in frames})
    inputs, issues = load_m4_inputs(paths)
    errors = [i for i in issues if getattr(i, "severity", "error") == "error"]
    assert not errors, errors[:5]

    g = inputs.genes
    acc = accepted_peptides(inputs.peptide_hits, inputs.peptides, _cfg())
    ev = gene_evidence(g, inputs.peptide_hits, acc, inputs.replicons)
    detected, status = two_peptide_status(ev.reset_index(drop=True), g["partial"])
    tu = layers_frame(layout).set_index("locus_tag")["tu_id"]
    p_true = expr.genes.set_index("locus_tag")["p_detect_true"]
    tags = g["locus_tag"].astype(str)
    p = tags.map(p_true).astype(float).to_numpy()
    p = np.where(status.to_numpy() == "partial_excluded", np.nan, p)
    genes = pd.DataFrame(
        {
            "locus_tag": tags.to_numpy(),
            "contig": g["contig"].astype(str).to_numpy(),
            "feature_type": g["feature_type"].astype(str).to_numpy(),
            "hypothetical": g["hypothetical"].astype(bool).to_numpy(),
            "partial": g["partial"].astype(bool).to_numpy(),
            "tu_id": tags.map(tu).astype(str).to_numpy(),
            "detected": detected.to_numpy(),
            "detection_status": status.to_numpy(),
            "p_detect": p,
            "cog_category": cog_categories(g, inputs.xrefs).to_numpy(),
        }
    )
    elig = genes["detection_status"] != "partial_excluded"
    diag = detectability_rows(genes.loc[elig, "detected"].to_numpy(),
                              genes.loc[elig, "p_detect"].to_numpy())  # fmt: skip
    tier = tier_rows(int(genes["detected"].sum()), int(ev["detected_2pep"].astype(bool).sum()))
    return SimFixture(inputs, genes, diag, tier)
