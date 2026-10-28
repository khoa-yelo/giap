"""Preflight rules part B T-PRE-11..22: inputs, proteins, peptides, abundance, evidence tier.

GIPPA_SPEC.md §1.1, §7.0, §7.3.2, §2.1 decision 3; ADR-0005 items 3-4; ADR-0006. Pure functions
over parsed records, validated frames and reader ``Issue`` lists; no file I/O. Thresholds come
only from the config: the decoy share is a value in ``peptides.content`` and never warns
(ADR-0005 item 4).
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from typing import Literal

import numpy as np
import pandas as pd

from gippa.config import GippaConfig
from gippa.config.model import AbundanceScale, InputConfig
from gippa.io.abundance import (
    AbundanceTable,
    Scale,
    abundance_stats,
    resolve_scale,
    validate_abundance,
)
from gippa.io.issues import Issue
from gippa.io.protein_fasta import ProteinRecord, route_by_length
from gippa.preflight.checks import LEVEL_ORDER, Check, InputName, JsonValue, Level, capped

__all__ = [
    "Tier",
    "abundance_checks",
    "candidate_tier",
    "inputs_checks",
    "peptide_checks",
    "protein_checks",
    "tier_check",
]

Tier = Literal["decoy_picked", "mayu", "two_peptide", "protein_list"]

DEFERRED = "deferred: no proteins.faa; IDs are matched to locus tags after annotation"
_NO_MATCH = "deferred: no abundance ID matches a protein; assumed keyed by locus tags"
_SCALES: tuple[AbundanceScale, ...] = ("auto", "log2", "linear")
# Peptide content issues whose rows are dropped from the evidence frame.
_DROPPING = frozenset({"non_letters", "short_peptide", "long_peptide", "duplicate_spectrum_row"})
# Columns that rank PSMs for target-decoy picking (§7.3.2: q-value, then PEP, then score).
_RANKING = ("q_value", "pep", "score")
# Abundance reader/validator issues judged by ``abundance.values``; all others are shape issues.
_VALUE_CODES = frozenset(
    {"non_numeric_values", "non_finite", "negative_on_linear", "too_many_missing"}
)


def _check(
    input_name: InputName,
    rule_id: str,
    level: Level,
    rule: str,
    value: JsonValue,
    message: str = "",
) -> Check:
    return Check(
        rule_id=rule_id, input=input_name, level=level, rule=rule, value=value, message=message
    )


def _issue_json(issues: Iterable[Issue]) -> list[JsonValue]:
    return [
        {"code": i.code, "severity": i.severity, "n_rows": i.n_rows, "message": i.message}
        for i in issues
    ]


def _issue_level(issues: Iterable[Issue]) -> Level:
    """The most severe issue severity (``"ok"`` when there are none)."""
    return max((i.severity for i in issues), key=LEVEL_ORDER.index, default="ok")


def _messages(issues: Iterable[Issue]) -> str:
    return "; ".join(i.message for i in issues)


# ---------------------------------------------------------------- inputs


def inputs_checks(
    has_peptides: bool, has_proteins: bool, *, n_usable_peptides: int | None = None
) -> list[Check]:
    """§1.1: ``proteins.faa`` is required unless peptides are given; peptides whose every row
    was dropped by validation (``n_usable_peptides == 0``) do not count as evidence."""
    value: dict[str, JsonValue] = {"peptides": has_peptides, "proteins": has_proteins}
    if n_usable_peptides is not None:
        value["n_usable_peptides"] = n_usable_peptides
    usable = has_peptides and n_usable_peptides != 0
    msg = ""
    if not (has_peptides or has_proteins):
        msg = "no evidence: give --peptides and/or --proteins"
    elif not (usable or has_proteins):
        msg = "no usable evidence: every peptide row was dropped and no --proteins was given"
    return [
        _check(
            "inputs",
            "inputs.evidence_present",
            "error" if msg else "ok",
            "--proteins required unless --peptides is given with usable rows → else error",
            value,
            msg,
        )
    ]


# ---------------------------------------------------------------- proteins


def protein_checks(
    records: Sequence[ProteinRecord], parse_error: str | None, peptide_len: tuple[int, int]
) -> list[Check]:
    """Evaluate every protein rule; only ``proteins.ids`` when the file gave no records."""
    ids_rule = "non-empty, unique IDs; UniProt sp|ACC|NAME and A;B groups parsed → else error"
    if parse_error is not None:
        return [_check("proteins", "proteins.ids", "error", ids_rule, parse_error, parse_error)]
    if not records:
        return [_check("proteins", "proteins.ids", "error", ids_rule, 0, "no FASTA records")]

    counts = Counter(r.id for r in records)
    dups: list[JsonValue] = [i for i, n in counts.items() if n > 1]
    dups_v, msg = capped(dups)
    if dups:
        ids_check = _check(
            "proteins", "proteins.ids", "error", ids_rule, dups_v, f"duplicate protein IDs {msg}"
        )
    else:
        ids_check = _shared_members_check(records, ids_rule)
    out = [ids_check]

    nonstandard: Counter[str] = Counter()
    other: Counter[str] = Counter()
    for r in records:
        nonstandard.update(dict(r.nonstandard))
        other.update(dict(r.other))
    stops = sum(r.internal_stops for r in records)
    alphabet: dict[str, JsonValue] = {
        "nonstandard": dict(sorted(nonstandard.items())),
        "internal_stops": stops,
    }
    if other:
        alphabet["other"] = dict(sorted(other.items()))
    out.append(
        _check(
            "proteins",
            "proteins.alphabet",
            "warn" if nonstandard or other or stops else "ok",
            "20 aa + U; X, B, Z, J and other letters counted; internal '*' counted → warn",
            alphabet,
        )
    )

    lo, hi = peptide_len
    routes = Counter[str]()
    dropped: list[JsonValue] = []
    for r in records:
        route = route_by_length(len(r.sequence), peptide_len)
        routes[route] += 1
        if route == "dropped":
            dropped.append(r.id)
    classes: dict[str, JsonValue] = {
        "protein": routes["protein"],
        "peptide": routes["peptide"],
        "dropped": routes["dropped"],
    }
    msg = ""
    if dropped:
        classes["dropped_ids"], msg = capped(dropped)
    out.append(
        _check(
            "proteins",
            "proteins.length_classes",
            "warn" if dropped else "info",
            f"≤ {hi} aa → peptide path, > {hi} aa → proteins; < {lo} aa dropped → warn",
            classes,
            msg,
        )
    )

    by_seq: dict[str, list[str]] = {}
    for r in records:
        by_seq.setdefault(r.sequence, []).append(r.id)
    groups: list[JsonValue] = [list(ids) for ids in by_seq.values() if len(ids) > 1]
    groups_v, msg = capped(groups)
    out.append(
        _check(
            "proteins",
            "proteins.duplicate_sequences",
            "info" if groups else "ok",
            "same sequence under several IDs → info",
            groups_v,
            msg,
        )
    )
    return out


def _shared_members_check(records: Sequence[ProteinRecord], rule: str) -> Check:
    """``proteins.ids`` for unique IDs: records whose groups share a member are info."""
    owners: dict[str, list[str]] = {}
    for r in records:
        for m in dict.fromkeys(r.members):
            owners.setdefault(m, []).append(r.id)
    shared = {m: ids for m, ids in owners.items() if len(ids) > 1}
    if not shared:
        return _check("proteins", "proteins.ids", "ok", rule, len(records))
    sharing: list[JsonValue] = list(dict.fromkeys(i for ids in shared.values() for i in ids))
    members: list[JsonValue] = [m for m in sorted(shared)]
    members_v, m_msg = capped(members)
    records_v, r_msg = capped(sharing)
    msgs = [f"{len(sharing)} protein records share a group member"]
    msgs += [f"shared members: {m_msg}"] if m_msg else []
    msgs += [f"records: {r_msg}"] if r_msg else []
    return _check(
        "proteins",
        "proteins.ids",
        "info",
        rule,
        {"n_records": len(records), "shared_members": members_v, "records": records_v},
        "; ".join(msgs),
    )


# ---------------------------------------------------------------- peptides


def _column_values(df: pd.DataFrame, col: str) -> np.ndarray:
    if col not in df.columns:
        return np.empty(0, dtype=np.float64)
    v = np.asarray(df[col].to_numpy(dtype="float64", na_value=np.nan), dtype=np.float64)
    return v[~np.isnan(v)]


def _decoy_with_ids(df: pd.DataFrame) -> bool:
    """Whether any decoy row names its protein IDs (the column test for ``decoy_picked``)."""
    if "is_decoy" not in df.columns or "protein_ids" not in df.columns:
        return False
    decoy = df["is_decoy"].to_numpy(dtype=bool, na_value=False)
    ids = df["protein_ids"].notna().to_numpy(dtype=bool)
    return bool((decoy & ids).any())


def _has_ranking(df: pd.DataFrame) -> bool:
    return any(_column_values(df, col).size for col in _RANKING)


def peptide_checks(
    issues: Sequence[Issue], df: pd.DataFrame | None, cfg: InputConfig
) -> list[Check]:
    """Schema and content rules over ``validate_peptides`` output; none without peptides.
    Error issues are schema findings; all other issues are content (§7.0)."""
    if df is None:
        return []
    schema = [i for i in issues if i.severity == "error"]
    content = [i for i in issues if i.severity != "error"]
    schema_level = _issue_level(schema)
    out = [
        _check(
            "peptides",
            "peptides.schema",
            schema_level,
            "required 'peptide' column; types and ranges of optional columns (§5.3) → else error",
            _issue_json(schema),
            _messages(schema),
        )
    ]
    if schema_level == "error" and df.empty:
        return out

    n_rows = len(df)
    decoy = (
        df["is_decoy"].to_numpy(dtype=bool, na_value=False)
        if "is_decoy" in df.columns
        else np.zeros(n_rows, dtype=bool)
    )
    q = _column_values(df, "q_value")
    target_q = _column_values(df.loc[~decoy], "q_value")
    q_fail = bool(target_q.size) and not bool((target_q <= cfg.q_value_max).any())
    no_ranking = _decoy_with_ids(df) and not _has_ranking(df)
    level = _issue_level(content)
    if (q_fail or no_ranking) and level != "error":
        level = "warn"
    msgs = [_messages(content)] if content else []
    if q_fail:
        msgs.append(f"no target row has q_value ≤ input.q_value_max ({cfg.q_value_max:g})")
    if no_ranking:
        msgs.append(
            "decoy rows but no ranking column (q_value, pep or score): decoy_picked is"
            " unavailable, the tier falls back"
        )
    out.append(
        _check(
            "peptides",
            "peptides.content",
            level,
            "non-letters after stripping or length outside input.peptide_len dropped; no target"
            " row with q_value ≤ input.q_value_max → warn; decoys without a ranking column → warn;"
            " decoy share reported",
            {
                "n_rows": n_rows,
                "n_dropped": sum(i.n_rows for i in content if i.code in _DROPPING),
                "decoy_share": float(decoy.mean()) if n_rows else None,
                "q_min": float(q.min()) if q.size else None,
                "q_max": float(q.max()) if q.size else None,
                "issues": _issue_json(content),
            },
            "; ".join(msgs),
        )
    )
    return out


def candidate_tier(df: pd.DataFrame | None, has_proteins: bool) -> Tier:
    """§7.3.2 tier the columns allow; the ≥ 80 % decoy-pairing test is applied in M4.
    ``decoy_picked`` needs decoy rows with protein IDs and a ranking column."""
    if df is None or df.empty:
        return "protein_list"
    if _decoy_with_ids(df) and _has_ranking(df):
        return "decoy_picked"
    if _column_values(df, "q_value").size or _column_values(df, "pep").size:
        return "mayu"
    return "two_peptide"


def tier_check(tier: Tier) -> Check:
    """Info check naming the candidate detection tier (§7.3.2)."""
    value: dict[str, JsonValue] = {"tier": tier}
    if tier == "decoy_picked":
        value["provisional"] = True
    return _check(
        "peptides",
        "peptides.evidence_tier",
        "info",
        "which detection tier the columns allow (§7.3.2)",
        value,
        "decoy pairing (≥ 80 % of decoys pair with a gene) is checked after annotation"
        if tier == "decoy_picked"
        else "",
    )


# ---------------------------------------------------------------- abundance


def _join_check(table: AbundanceTable, proteins: Sequence[ProteinRecord] | None) -> Check:
    rule = "abundance IDs missing from proteins, and proteins missing from abundance → warn"
    if proteins is None:
        return _check("abundance", "abundance.join", "info", rule, DEFERRED)
    prot_members = {m for p in proteins for m in p.members}
    abund_members = {m for i in table.ids for m in i.split(";") if m}
    not_in_prot: list[JsonValue] = [
        i for i in table.ids if i and not prot_members.intersection(i.split(";"))
    ]
    if len(not_in_prot) == sum(1 for i in table.ids if i):
        return _check("abundance", "abundance.join", "info", rule, _NO_MATCH)
    not_in_abund: list[JsonValue] = [
        p.id for p in proteins if not abund_members.intersection(p.members)
    ]
    a_v, a_msg = capped(not_in_prot)
    p_v, p_msg = capped(not_in_abund)
    return _check(
        "abundance",
        "abundance.join",
        "warn" if not_in_prot or not_in_abund else "ok",
        rule,
        {"abundance_not_in_proteins": a_v, "proteins_not_in_abundance": p_v},
        "; ".join(m for m in (a_msg, p_msg) if m),
    )


def abundance_checks(
    table: AbundanceTable | None,
    issues: Sequence[Issue],
    proteins: Sequence[ProteinRecord] | None,
    scale: Scale,
    requested_scale: str,
    peptides_have_intensity: bool,
    cfg: GippaConfig,
) -> list[Check]:
    """Evaluate every abundance rule over ``read_abundance`` output; only ``abundance.shape``
    when the table is unusable. The scale is resolved once here from ``requested_scale``;
    ``scale`` is accepted for call compatibility and ignored."""
    if requested_scale not in _SCALES:
        raise ValueError(f"requested_scale must be one of {_SCALES}, got {requested_scale!r}")
    requested: AbundanceScale = requested_scale
    shape_issues = [i for i in issues if i.code not in _VALUE_CODES]
    shape_rule = "first column = ID, one numeric column, ≥ 1 row → else error; extra columns warn"
    if table is None:
        return [
            _check(
                "abundance",
                "abundance.shape",
                "error",
                shape_rule,
                _issue_json(shape_issues),
                _messages(shape_issues) or "abundance table is unusable",
            )
        ]
    out = [
        _check(
            "abundance",
            "abundance.shape",
            _issue_level(shape_issues),
            shape_rule,
            {
                "id_column": table.id_column,
                "value_column": table.value_column,
                "ignored_columns": list(table.ignored_columns),
                "issues": _issue_json(shape_issues),
            },
            _messages(shape_issues),
        ),
        _join_check(table, proteins),
    ]

    frac = cfg.preflight.abundance_missing_warn_frac
    scale_issues = validate_abundance(table, requested, frac)
    value_issues = [i for i in issues if i.code in _VALUE_CODES]
    value_issues += [i for i in scale_issues if i.code in _VALUE_CODES]
    stats = abundance_stats(table)
    resolved = resolve_scale(table.values, requested)
    out.append(
        _check(
            "abundance",
            "abundance.values",
            _issue_level(value_issues),
            f"non-numeric → error; negatives on a linear scale → error; > {frac:.0%} missing →"
            " warn; zeros and NA counted",
            {
                "n_rows": stats.n_rows,
                "n_na": stats.n_na,
                "n_zero": stats.n_zero,
                "n_negative": stats.n_negative,
                "n_bad": sum(table.bad_rows),
                "missing_frac": stats.missing_frac,
                "issues": _issue_json(value_issues),
            },
            _messages(value_issues),
        )
    )

    undetermined = [i for i in scale_issues if i.code not in _VALUE_CODES]
    msgs = [_messages(undetermined)] if undetermined else []
    if peptides_have_intensity:
        msgs.append("peptide intensities given: the abundance table is not used (§2.1 decision 3)")
    out.append(
        _check(
            "abundance",
            "abundance.scale",
            "info",
            "log2 or linear auto-detected; --abundance-scale overrides",
            {
                "scale": resolved,
                "source": "auto" if requested == "auto" else "cli",
                "used": not peptides_have_intensity,
            },
            "; ".join(msgs),
        )
    )
    return out
