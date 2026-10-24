"""Matching helpers of the M7 variants gate (docs/specs/M7-13-calibration-gate.md rules 3, 4, 9;
interface pinned in ``test_m7_gate_helpers.py``). Pure functions: written ``variants.tsv`` and
giapsim truth frames in, integer counts out. Coordinates are the 1-based lowest base of the codon,
as both ``variants.tsv`` and ``truth/alt_sites.tsv`` write them (spec "Degenerate cases": minus
strand).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from types import MappingProxyType

import pandas as pd

#: replicates of the 1M novel-ORF gate (lead decision 1; ADR-0048 pending Khoa)
GATE_R: Mapping[str, int] = MappingProxyType({"fast": 10, "full": 300})
#: rule 4: a codon within this many nt of an indel junction is masked
INDEL_MASK_NT = 3
#: giapsim de novo error types that are mass-equivalent (``near_block`` and ``swap`` are not)
ISOBARIC_ERRORS: frozenset[str] = frozenset({"block", "modification"})

PRED_COLUMNS = ("contig", "codon_start", "strand", "alt_aa", "mass_equivalent")
SiteKey = tuple[str, int, str, str]  # (contig, codon_start, strand, folded ALT)


def fold_il(s: str) -> str:
    """``s`` with I, L and J folded to one letter (L); other residues unchanged."""
    return str(s).replace("I", "L").replace("J", "L")


def _flags(value: object) -> frozenset[str]:
    return frozenset(f for f in str(value).split(";") if f and f != "NA")


def _site(contig: object, codon_start: object, strand: object, alt_aa: object) -> SiteKey:
    return (str(contig), int(str(codon_start)), str(strand), fold_il(str(alt_aa)))


def g1_predictions(variants: pd.DataFrame) -> pd.DataFrame:
    """G1 prediction set: ``high`` substitution rows, one per (site, folded ALT residue).

    ``mass_equivalent`` is true when any row of the group carries the flag.
    """
    sel = variants[(variants["confidence"] == "high") & (variants["alt_type"] == "substitution")]
    rows: dict[SiteKey, list[object]] = {}
    for c, pos, s, alt, flags in zip(sel["contig"], sel["codon_start"], sel["strand"],
                                     sel["alt_aa"], sel["flags"], strict=True):  # fmt: skip
        key = _site(c, pos, s, alt)
        me = "mass_equivalent" in _flags(flags)
        if key in rows:
            rows[key][4] = bool(rows[key][4]) or me
        else:
            rows[key] = [key[0], key[1], key[2], str(alt), me]
    out = pd.DataFrame(list(rows.values()), columns=list(PRED_COLUMNS))
    return out.astype({"codon_start": "int64", "mass_equivalent": bool})


def _masked(contig: str, codon_start: int, junctions: Mapping[str, list[int]]) -> bool:
    lo, hi = codon_start - INDEL_MASK_NT, codon_start + 2 + INDEL_MASK_NT
    return any(lo <= j <= hi for j in junctions.get(contig, ()))


def _recall_eligible(ref_aa: object, alt_aa: object, n_true: object) -> bool:
    """Rule 3: >= 2 true peptides and MS-observable (not an I <-> L exchange)."""
    return int(str(n_true)) >= 2 and not {str(ref_aa), str(alt_aa)} <= {"I", "L"}


def match_variants(
    pred: pd.DataFrame, alt_sites: pd.DataFrame, indel_sites: pd.DataFrame
) -> dict[str, int]:
    """Counts ``{n_pred, tp_p, n_truth, tp_r, masked}`` of one replicate (rules 3-4).

    Masked predictions (codon ``[c, c+2]`` within ``INDEL_MASK_NT`` of a junction on the same
    contig) are dropped before matching and counted in ``masked``.
    """
    junctions: dict[str, list[int]] = {}
    for c, j in zip(indel_sites["contig"], indel_sites["junction"], strict=True):
        junctions.setdefault(str(c), []).append(int(j))
    kept: list[SiteKey] = []
    masked = 0
    for c, pos, s, alt in zip(pred["contig"], pred["codon_start"], pred["strand"],
                              pred["alt_aa"], strict=True):  # fmt: skip
        key = _site(c, pos, s, alt)
        if _masked(key[0], key[1], junctions):
            masked += 1
        else:
            kept.append(key)
    truth_all: set[SiteKey] = set()
    truth_recall: list[SiteKey] = []
    for c, pos, s, ref, alt, n in zip(alt_sites["contig"], alt_sites["codon_start"],
                                      alt_sites["strand"], alt_sites["ref_aa"],
                                      alt_sites["alt_aa"], alt_sites["n_true_peptides"],
                                      strict=True):  # fmt: skip
        key = _site(c, pos, s, alt)
        truth_all.add(key)
        if _recall_eligible(ref, alt, n):
            truth_recall.append(key)
    predicted = set(kept)
    return {
        "n_pred": len(kept),
        "tp_p": sum(k in truth_all for k in kept),
        "n_truth": len(truth_recall),
        "tp_r": sum(k in predicted for k in truth_recall),
        "masked": masked,
    }


def _high(variants: pd.DataFrame) -> pd.DataFrame:
    return variants[variants["confidence"] == "high"]


def g2_violations(variants: pd.DataFrame) -> int:
    """``high`` rows flagged neither ``start_codon`` nor ``stop_codon`` (G2)."""
    high = _high(variants)
    return sum(not ({"start_codon", "stop_codon"} & _flags(f)) for f in high["flags"])


def start_codon_high(variants: pd.DataFrame) -> int:
    """``high`` rows flagged ``start_codon`` (recorded, rule 10)."""
    return sum("start_codon" in _flags(f) for f in _high(variants)["flags"])


def a6a_unflagged(
    variants: pd.DataFrame,
    site_peptides: Mapping[tuple[str, int, str, str, str], Iterable[str]],
    denovo_truth: pd.DataFrame,
) -> dict[str, int]:
    """A6a: substitution/block rows supported only by isobaric de novo errors, and how many of
    them lack ``mass_equivalent`` (``{n_rows, n_unflagged}``)."""
    errors: dict[str, list[str]] = {}
    for p, e in zip(denovo_truth["peptide"], denovo_truth["error_type"], strict=True):
        errors.setdefault(fold_il(str(p)), []).append(str(e))
    n_rows = n_unflagged = 0
    for c, pos, s, t, alt, flags in zip(variants["contig"], variants["codon_start"],
                                        variants["strand"], variants["alt_type"],
                                        variants["alt_aa"], variants["flags"],
                                        strict=True):  # fmt: skip
        if t not in ("substitution", "block"):
            continue
        peps = site_peptides.get((str(c), int(str(pos)), str(s), str(t), str(alt)), ())
        reads = [e for p in peps for e in errors.get(fold_il(p), ())]
        if reads and all(e in ISOBARIC_ERRORS for e in reads):
            n_rows += 1
            n_unflagged += "mass_equivalent" not in _flags(flags)
    return {"n_rows": n_rows, "n_unflagged": n_unflagged}


def a6_high_swaps(variants: pd.DataFrame, min_peptides: int = 2) -> dict[str, int]:
    """A6b/A6c (rule 9): swap tokens folded per site; ``folded`` counts the tokens with a ``high``
    row or ``>= min_peptides`` ALT peptides in total, ``unfolded`` the ``high`` swap rows."""
    swaps = variants[variants["alt_type"] == "swap"]
    groups: dict[SiteKey, tuple[bool, int]] = {}
    for c, pos, s, alt, conf, n_alt in zip(swaps["contig"], swaps["codon_start"],
                                           swaps["strand"], swaps["alt_aa"],
                                           swaps["confidence"], swaps["n_alt"],
                                           strict=True):  # fmt: skip
        key = _site(c, pos, s, alt)
        high, n = groups.get(key, (False, 0))
        groups[key] = (high or conf == "high", n + int(str(n_alt)))
    return {
        "folded": sum(h or n >= min_peptides for h, n in groups.values()),
        "unfolded": int((swaps["confidence"] == "high").sum()),
    }


# --------------------------------------------------------------------------- novel ORFs (M7-13b)

#: the four gated novel-ORF classes (rule 5)
NOVEL_CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
NOVEL_COUNT_KEYS = ("n_pred", "tp_p", "n_truth", "tp_r", "start_ok", "mismatch", "rbs_planted",
                    "rbs_pass")  # fmt: skip
NOVEL_PRED_COLUMNS = ("orf_id", "class", "contig", "strand", "stop", "start_longest",
                      "start_best_ostir", "rbs_supported")  # fmt: skip
NovelKey = tuple[str, str, int]  # (contig, strand, stop: 1-based lowest base of the stop codon)


def _wrap(pos: int, length: int | None) -> int:
    return pos if length is None else (pos - 1) % length + 1


def novel_stop(strand: str, end: int, length: int | None = None) -> int:
    """1-based lowest base of the stop codon from ``novel_orfs.tsv`` ``end`` (the last base of the
    stop in transcription direction): ``end - 2`` on ``+``, ``end`` on ``-``; wrapped into
    ``1..length`` when ``length`` is given (origin-crossing keys, rule 5)."""
    return _wrap(int(end) - 2 if str(strand) == "+" else int(end), length)


def _true(v: object) -> bool:
    if isinstance(v, str):
        return v.strip().lower() == "true"
    return False if pd.isna(v) else bool(v)  # type: ignore[call-overload]


def _int_or_na(v: object) -> object:
    if isinstance(v, str):
        return pd.NA if v in ("", "NA") else int(v)
    return pd.NA if pd.isna(v) else int(v)  # type: ignore[call-overload]


def _pred_frame(
    rows: Iterable[Mapping[str, object]], lengths: Mapping[str, int] | None
) -> pd.DataFrame:
    lens = lengths or {}
    data = [
        [str(r["orf_id"]), str(r["class"]), str(r["contig"]), str(r["strand"]),
         novel_stop(str(r["strand"]), int(str(r["end"])), lens.get(str(r["contig"]))),
         _int_or_na(r["start_longest"]), _int_or_na(r["start_best_ostir"]),
         _true(r["rbs_supported"])]
        for r in rows
    ]  # fmt: skip
    out = pd.DataFrame(data, columns=list(NOVEL_PRED_COLUMNS))
    return out.astype({"stop": "int64", "start_longest": "Int64", "start_best_ostir": "Int64",
                       "rbs_supported": bool})  # fmt: skip


def g3_predictions(
    novel_orfs: pd.DataFrame, lengths: Mapping[str, int] | None = None
) -> pd.DataFrame:
    """G3 prediction set: ``called = true`` rows of the four classes (``NOVEL_PRED_COLUMNS``)."""
    called = pd.Series([_true(c) for c in novel_orfs["called"]], index=novel_orfs.index, dtype=bool)
    sel = novel_orfs[called & novel_orfs["class"].isin(NOVEL_CLASSES)]
    return _pred_frame(sel.to_dict("records"), lengths)


def g4_predictions(
    novel_orfs: pd.DataFrame,
    side: pd.DataFrame,
    level: float | None,
    *,
    regime: str | None = None,
    lengths: Mapping[str, int] | None = None,
) -> pd.DataFrame:
    """G4 prediction set: side rows with ``class_rule_met`` and ``q <= level`` (NaN q never enters;
    ``level=None``: no q filter), of ``regime`` when given, of the four classes by side ``class``;
    coordinates from ``novel_orfs`` (ValueError when a side row is missing there)."""
    s = side[side["class_rule_met"].astype(bool) & side["class"].isin(NOVEL_CLASSES)]
    if regime is not None:
        s = s[s["regime"] == regime]
    if level is not None:
        s = s[s["q"].astype(float) <= float(level)]  # NaN compares false
    by_id = {str(r["orf_id"]): r for r in novel_orfs.to_dict("records")}
    rows = []
    for oid, cls in zip(s["orf_id"], s["class"], strict=True):
        if str(oid) not in by_id:
            raise ValueError(f"g4_predictions: side row {oid!r} is not in novel_orfs.tsv")
        rows.append({**by_id[str(oid)], "class": str(cls)})
    return _pred_frame(rows, lengths)


def match_novel(
    pred: pd.DataFrame,
    truth: pd.DataFrame,
    lengths: Mapping[str, int] | None = None,
    min_peptides: int = 2,
) -> dict[str, dict[str, int]]:
    """Counts ``NOVEL_COUNT_KEYS`` per class of ``NOVEL_CLASSES`` (rule 5; interface pinned in
    ``test_m7_gate_novel_helpers.py``): match by ``(contig, strand, stop)`` and class."""
    lens = lengths or {}
    truth_rows: dict[str, dict[NovelKey, dict[str, object]]] = {c: {} for c in NOVEL_CLASSES}
    recall: dict[str, list[NovelKey]] = {c: [] for c in NOVEL_CLASSES}
    for r in truth.to_dict("records"):
        cls = str(r["class"])
        if cls not in truth_rows:
            continue
        key = (
            str(r["contig"]),
            str(r["strand"]),
            _wrap(int(r["stop"]), lens.get(str(r["contig"]))),
        )
        truth_rows[cls].setdefault(key, r)
        if int(r["n_true_peptides"]) >= min_peptides:
            recall[cls].append(key)
    first: dict[str, dict[NovelKey, dict[str, object]]] = {c: {} for c in NOVEL_CLASSES}
    out = {c: dict.fromkeys(NOVEL_COUNT_KEYS, 0) for c in NOVEL_CLASSES}
    for p in pred.to_dict("records"):
        cls = str(p["class"])
        if cls not in out:
            continue
        key = (str(p["contig"]), str(p["strand"]), int(p["stop"]))
        out[cls]["n_pred"] += 1
        if key in truth_rows[cls]:
            first[cls].setdefault(key, p)
        elif any(key in truth_rows[c] for c in NOVEL_CLASSES if c != cls):
            out[cls]["mismatch"] += 1
    for cls in NOVEL_CLASSES:
        o = out[cls]
        o["tp_p"] = len(first[cls])
        o["n_truth"] = len(recall[cls])
        o["tp_r"] = sum(k in first[cls] for k in recall[cls])
        for key, p in first[cls].items():
            t = truth_rows[cls][key]
            starts = (p["start_longest"], p["start_best_ostir"])
            want = int(str(t["true_start"]))
            o["start_ok"] += any(
                _int_or_na(s) is not pd.NA and _int_or_na(s) == want for s in starts
            )
            if _true(t["rbs_planted"]):
                o["rbs_planted"] += 1
                o["rbs_pass"] += bool(p["rbs_supported"])
    return out


def _per_class(values: Iterable[object]) -> dict[str, int]:
    out = dict.fromkeys(NOVEL_CLASSES, 0)
    for c in values:
        if c in out:
            out[str(c)] += 1
    return out


def list_targets(side: pd.DataFrame, regime: str) -> dict[str, int]:
    """Side rows of ``regime`` per ``fdr_class`` that entered the class list (finite ``q`` or
    ``reason = "lost_pick"``): the rule 6 vacuity input."""
    s = side[side["regime"] == regime]
    q = s["q"].astype(float).to_numpy()
    entered = [math.isfinite(x) or r == "lost_pick" for x, r in zip(q, s["reason"], strict=True)]
    return _per_class(s.loc[entered, "fdr_class"])


def reason_counts(side: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Non-null ``reason`` values per side ``class`` (the four classes)."""
    out: dict[str, dict[str, int]] = {c: {} for c in NOVEL_CLASSES}
    for c, r in zip(side["class"], side["reason"], strict=True):
        if c in out and r is not None and not (isinstance(r, float) and math.isnan(r)):
            out[str(c)][str(r)] = out[str(c)].get(str(r), 0) + 1
    return out


def fallback_calls(side: pd.DataFrame) -> dict[str, int]:
    """``called`` rows with a non-finite ``q`` per side ``class`` (the >= 2-peptide fallback; not
    error-controlled, recorded and not gated)."""
    q = side["q"].astype(float).to_numpy()
    fb = [bool(c) and not math.isfinite(x) for c, x in zip(side["called"], q, strict=True)]
    return _per_class(side.loc[fb, "class"])
