"""Promoter / terminator evidence "yes", abundance patterns and the seven operon flags
(GIAP_SPEC.md §7.4.5, §7.4.7, §9.3; docs/specs/M6-07-flags-and-operons-table.md, pass 1).

:func:`merge_tests` writes the M6-05 (``Discordance.rows``) and M6-06 (``MissedOperons.rows``)
test columns into the M6-04 pair frame; :func:`flag_pairs` turns the merged frame into
``operon_flags.tsv``. The flag is a pure function of (class, promoter yes, terminator yes,
pattern), :func:`flag_for`, applied first match wins.

Promotech and Promoter Calculator see sigma-70 promoters only, and TerminatorNet and TransTermHP see
intrinsic terminators only. A missing call therefore never proves a missing element: the
``evidence_note`` never asserts absence, and a layer the run did not have (named by the caller,
``layer_missing:<layer>``) is said to be missing rather than read as "no".

Pure and deterministic: no I/O, no RNG.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.config import GiapConfig
from giap.operons.missed import IMPUTED, UNESTIMATED, _bool, _flags, _per_gene
from giap.operons.pairs import PAIR_COLUMNS, TEST_COLUMNS

__all__ = [
    "FLAGS",
    "FLAG_COLUMNS",
    "PATTERNS",
    "TU_ERROR_FLAG",
    "TU_SOURCE_FLAG",
    "abundance_pattern",
    "add_tu_source_flag",
    "evidence_yes",
    "flag_for",
    "flag_pairs",
    "merge_tests",
]

FLAGS: tuple[str, ...] = (
    "internal_promoter_supported", "internal_terminator_supported", "unexplained_discordance",
    "conflict_terminator_in_operon", "novel_operon_candidate", "coabundant_despite_terminator",
    "other",
)  # fmt: skip
PATTERNS: tuple[str, ...] = (
    "downstream_higher", "upstream_higher", "discordant", "concordant", "coabundant", "none",
)  # fmt: skip
#: §9.3 ``operon_flags.tsv`` columns, in table order
FLAG_COLUMNS: tuple[str, ...] = (
    "pair_id", "gene_a", "gene_b", "uniop_class", "promoter_between_yes",
    "terminator_between_yes", "abundance_pattern", "flag", "q_value", "lfdr", "evidence_note",
)  # fmt: skip
LAYER_MISSING = "layer_missing:"
#: ADR-0061 decision 1: rows whose calibration does not hold with distance-model TUs
TU_SOURCE_FLAG = "tu_source_distance_model"
#: ADR-0061 decision 4: TU-block Moran's I rows assume exact TU boundaries (every run)
TU_ERROR_FLAG = "tu_block_tu_error"
_TU_TABLES = ("autocorrelation", "pairs")
_DISCORDANT = frozenset(PATTERNS[:3])
#: the §7.4.5 table: class, promoter yes (None = any), terminator yes, patterns, flag
_TABLE: tuple[tuple[str, bool | None, bool, frozenset[str], str], ...] = (
    ("OP", True, False, frozenset({"downstream_higher"}), FLAGS[0]),
    ("OP", False, True, frozenset({"upstream_higher"}), FLAGS[1]),
    ("OP", False, False, _DISCORDANT, FLAGS[2]),
    ("OP", None, True, frozenset({"concordant"}), FLAGS[3]),
    ("NO", False, False, frozenset({"coabundant"}), FLAGS[4]),
    ("NO", None, True, frozenset({"coabundant"}), FLAGS[5]),
)


def flag_for(cls: str, promoter_yes: bool, terminator_yes: bool, pattern: str) -> str | None:
    """The §7.4.5 flag of one pair, first match wins; ``None`` = no flag row."""
    for c, p, t, pats, flag in _TABLE:
        if c == cls and p in (None, promoter_yes) and t == terminator_yes and pattern in pats:
            return flag
    if cls in ("OP", "NO") and (pattern in _DISCORDANT or pattern == "coabundant"):
        return "other"
    return None


def _has(pairs: pd.DataFrame, token: str) -> npt.NDArray[np.bool_]:
    return np.array([token in _flags(v) for v in pairs["flags"]], dtype=np.bool_)


def evidence_yes(pairs: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame:
    """``promoter_between_yes`` / ``terminator_between_yes``: score at or above the config
    cut-off; NA, a missing layer or an ``unestimated`` row is never "yes"."""
    ann = cfg.annotation
    bad = _has(pairs, UNESTIMATED)
    out = {}
    for name, col, layer, cut in (
        ("promoter_between_yes", "promoter_between", "promoters", ann.promoters.promotech_yes),
        ("terminator_between_yes", "terminator_between", "terminators",
         ann.terminators.terminatornet_yes),
    ):  # fmt: skip
        v = pd.to_numeric(pairs[col], errors="coerce").to_numpy(dtype=float)
        ok = np.nan_to_num(v, nan=-math.inf) >= cut
        out[name] = ok & ~bad & ~_has(pairs, LAYER_MISSING + layer)
    return pd.DataFrame(out, index=pairs.index, columns=list(out) or None).astype(bool)


def abundance_pattern(pairs: pd.DataFrame, alpha: float) -> pd.Series:
    """One of :data:`PATTERNS` per pair (spec step 2); an ``unestimated`` row is ``none``."""
    q = pd.to_numeric(pairs["q_value"], errors="coerce").to_numpy(dtype=float)
    d = pd.to_numeric(pairs["delta_log2"], errors="coerce").to_numpy(dtype=float)
    det_a, det_b = _bool(pairs["detected_a"]), _bool(pairs["detected_b"])
    cls = pairs["class"].to_numpy(dtype=object)
    sig, tested = q <= alpha, ~np.isnan(q) & ~_has(pairs, UNESTIMATED)
    op, both = (cls == "OP") & tested, det_a & det_b
    lo = op & sig & ((both & (d < 0)) | (det_a & ~det_b))
    hi = op & sig & ((both & (d > 0)) | (~det_a & det_b))
    conds = [
        hi,
        lo,
        op & sig & (det_a | det_b),
        op & ~sig & both,
        (cls == "NO") & tested & sig,
    ]
    pat = np.select(conds, list(PATTERNS[:5]), default="none")
    return pd.Series(pat, index=pairs.index, dtype=object)


def _num(x: object, fmt: str = ".3g") -> str:
    return format(float(x), fmt)  # type: ignore[arg-type]


def _note(r: pd.Series, pattern: str) -> str:
    """A short sentence naming the evidence used; never asserts absence (§7.4.5 caveat)."""
    fl = _flags(r["flags"])
    parts = []
    for layer, col in (("promoter", "promoter_between"), ("terminator", "terminator_between")):
        v = float(r[col]) if pd.notna(r[col]) else math.nan
        if LAYER_MISSING + layer + "s" in fl:
            parts.append(f"{layer} layer missing")
        elif math.isnan(v):
            parts.append(f"no {layer} call between")
        else:
            parts.append(f"{layer} {_num(v, '.2f')} between")
    parts.append(f"q={_num(r['q_value'])}")
    d = float(r["delta_log2"]) if pd.notna(r["delta_log2"]) else math.nan
    side = pattern.split("_")[0]
    if pattern.endswith("_higher") and math.isnan(d):
        parts.append(f"only the {side} gene detected")
    elif pattern.endswith("_higher"):
        parts.append(f"{side} higher by {_num(abs(d))} log2")
    elif not math.isnan(d):
        parts.append(f"{pattern}, |delta| {_num(abs(d))} log2")
    if IMPUTED in fl:
        parts.append("abundance imputed")
    return "; ".join(parts)


def flag_pairs(pairs: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame:
    """``operon_flags.tsv``: the flagged rows of the merged pair frame, by ``pair_id`` number."""
    ev = evidence_yes(pairs, cfg)
    pat = abundance_pattern(pairs, cfg.stats.alpha.operon_flags)
    flag = [
        flag_for(c, bool(p), bool(t), s)
        for c, p, t, s in zip(
            pairs["class"], ev["promoter_between_yes"], ev["terminator_between_yes"], pat,
            strict=True,
        )
    ]  # fmt: skip
    keep = np.array([f is not None for f in flag], dtype=bool)
    sel = pairs.loc[keep]
    out = pd.DataFrame(
        {
            "pair_id": sel["pair_id"].astype(object),
            "gene_a": sel["gene_a"].astype(object),
            "gene_b": sel["gene_b"].astype(object),
            "uniop_class": sel["class"].astype(object),
            "promoter_between_yes": ev.loc[keep, "promoter_between_yes"].astype(bool),
            "terminator_between_yes": ev.loc[keep, "terminator_between_yes"].astype(bool),
            "abundance_pattern": pat[keep].astype(object),
            "flag": pd.Series([f for f in flag if f is not None], index=sel.index, dtype=object),
            "q_value": pd.to_numeric(sel["q_value"], errors="coerce").astype(float),
            "lfdr": pd.to_numeric(sel["lfdr"], errors="coerce").astype(float),
            "evidence_note": pd.Series(
                [_note(r, p) for (_, r), p in zip(sel.iterrows(), pat[keep], strict=True)],
                index=sel.index,
                dtype=object,
            ),
        },
        columns=list(FLAG_COLUMNS),
    )
    order = np.argsort([int(str(p).rsplit("_", 1)[-1]) for p in out["pair_id"]], kind="stable")
    return out.iloc[order].reset_index(drop=True)


def merge_tests(
    pairs: pd.DataFrame,
    discordance_rows: pd.DataFrame,
    missed_rows: pd.DataFrame,
    genes: pd.DataFrame,
    *,
    missing_layers: Sequence[str] = (),
) -> pd.DataFrame:
    """The M6-04 pair frame with the M6-05/06 test columns written in (``which`` dropped, lead
    decision 4), flags united, ``abundance_imputed`` (decision 6) and, on OP and NO rows,
    ``layer_missing:<layer>`` for each layer the caller names as missing (Q1)."""
    out = pairs.copy()
    parts = [r.drop(columns=["which"], errors="ignore") for r in (discordance_rows, missed_rows)]
    rows = pd.concat([r for r in parts if len(r)] or parts[:1], ignore_index=True)
    rows = rows.drop_duplicates("pair_id", keep="first").set_index("pair_id")
    hit = out["pair_id"].isin(rows.index).to_numpy()
    ids = out.loc[hit, "pair_id"]
    for col in ("test", *TEST_COLUMNS) if hit.any() else ():
        if col in ("test", "method", "null_model"):
            out[col] = out[col].astype(object)
        vals = rows.loc[ids, col].to_numpy()
        out.loc[hit, col] = pd.array(vals, dtype="Int64") if col == "n_perm" else vals
    extra = {str(k): _flags(v) for k, v in rows["flags"].items()}
    imputed = _per_gene(genes, IMPUTED, out["gene_a"]) | _per_gene(genes, IMPUTED, out["gene_b"])
    op_no = out["class"].isin(["OP", "NO"]).to_numpy()
    merged = []
    for pid, f, imp, on in zip(out["pair_id"], out["flags"], imputed, op_no, strict=True):
        fl = [*_flags(f), *extra.get(str(pid), [])]
        fl += [IMPUTED] if imp else []
        fl += [LAYER_MISSING + m for m in missing_layers] if on else []
        merged.append(list(dict.fromkeys(fl)))
    out["flags"] = pd.Series(merged, index=out.index, dtype=object)
    return out[list(PAIR_COLUMNS)]


def _tu_targets(frame: pd.DataFrame, table: str) -> npt.NDArray[np.bool_]:
    if table == "autocorrelation":
        hit = (frame["statistic"] == "morans_i") & (frame["null_model"] == "tu_block")
    else:
        hit = frame["test"] == "both_detected_t"
    return np.asarray(hit.fillna(False), dtype=np.bool_)


def _with_tokens(cell: object, add: Sequence[str]) -> object:
    """``cell`` with the ``add`` tokens appended (existing tokens kept, none duplicated), in the
    cell's own form: a ``;``-joined string, a tuple or a list."""
    old = _flags(cell)
    new = [t for t in add if t not in old]
    if not new:
        return cell
    tokens = [*old, *new]
    if isinstance(cell, str):
        return ";".join(tokens)
    return tuple(tokens) if isinstance(cell, tuple) else tokens


def add_tu_source_flag(frame: pd.DataFrame, table: str, operon_source: str | None) -> pd.DataFrame:
    """ADR-0061 row flags of ``autocorrelation`` or ``pairs`` rows (a new frame; pure).

    Every ``morans_i`` row with ``null_model = tu_block`` gets :data:`TU_ERROR_FLAG` whatever
    the source (decision 4); when ``operon_source`` is ``distance_model``, those rows and every
    ``both_detected_t`` pair row also get :data:`TU_SOURCE_FLAG` (decision 1). Only the
    ``flags`` cells of the target rows change. ``ValueError`` for any other ``table``.
    """
    if table not in _TU_TABLES:
        raise ValueError(f"add_tu_source_flag: unknown table {table!r}; use one of {_TU_TABLES}")
    add = [TU_ERROR_FLAG] if table == "autocorrelation" else []
    if operon_source == "distance_model":
        add.append(TU_SOURCE_FLAG)
    out = frame.copy()
    if not add or not len(out):
        return out
    hit = _tu_targets(out, table)
    cells = [_with_tokens(c, add) if h else c for c, h in zip(out["flags"], hit, strict=True)]
    out["flags"] = pd.Series(cells, index=out.index, dtype=object)
    return out
