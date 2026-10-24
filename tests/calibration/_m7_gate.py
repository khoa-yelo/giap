"""Matching helpers of the M7 variants gate (docs/specs/M7-13-calibration-gate.md rules 3, 4, 9;
interface pinned in ``test_m7_gate_helpers.py``). Pure functions: written ``variants.tsv`` and
giapsim truth frames in, integer counts out. Coordinates are the 1-based lowest base of the codon,
as both ``variants.tsv`` and ``truth/alt_sites.tsv`` write them (spec "Degenerate cases": minus
strand).
"""

from __future__ import annotations

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
