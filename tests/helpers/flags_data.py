"""Hand-built merged pair frames, TU tables and gene tables for M6-07 (flags and ``operons.tsv``;
docs/specs/M6-07-flags-and-operons-table.md; test author).

``row`` builds one ``pairs.tsv`` row (M6-04 columns, ``helpers.pairs_data.PAIR_COLS``) in the
state M6-07 receives it: the M6-05 / M6-06 test columns already merged in. ``frame`` numbers the
rows ``PAIR_1..n`` with genes ``G<n>a`` (upstream) / ``G<n>b`` (downstream) unless the row names
its genes. A tested OP row carries ``ci_unestimated`` (M6-07 lead decision 5), as M6-05 writes it.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from gippa.config import GippaConfig
from helpers.pairs_data import PAIR_COLS, TEST_COLS

NAN = math.nan
BASE_LOG2 = 10.0
TEST_OF = {"OP": "both_detected_t", "NO": "missed_operon_decoy"}


def row(
    cls: str,
    *,
    prom: float = NAN,
    term: float = NAN,
    delta: float = NAN,
    q: float = NAN,
    lfdr: float = NAN,
    detected: tuple[bool, bool] = (True, True),
    flags: Sequence[str] | None = None,
    gene_a: str | None = None,
    gene_b: str | None = None,
    test: str | None = None,
) -> dict[str, Any]:
    """One merged pair row; ``delta`` sets ``abundance_log2_b = 10 + delta`` when both are
    detected; an undetected gene has NaN abundance (and then ``delta_log2`` is NaN)."""
    det_a, det_b = detected
    a = BASE_LOG2 if det_a else NAN
    if det_b:
        b = BASE_LOG2 + (0.0 if math.isnan(delta) else delta) if det_a else BASE_LOG2
    else:
        b = NAN
    d = b - a if det_a and det_b else NAN
    tested = not math.isnan(q)
    if flags is None:
        flags = ["ci_unestimated"] if (tested and cls == "OP") else []
    out: dict[str, Any] = {c: NAN for c in PAIR_COLS}
    out.update(
        {
            "class": cls,
            "gene_a": gene_a,
            "gene_b": gene_b,
            "strand": "." if cls in ("CV", "DV") else "+",
            "intergenic_bp": 40,
            "uniop_p": 0.9 if cls == "OP" else (0.1 if cls == "NO" else NAN),
            "promoter_between": prom,
            "tx_rate_between": NAN,
            "terminator_between": term,
            "detected_a": det_a,
            "detected_b": det_b,
            "p_detect_a": 0.9,
            "p_detect_b": 0.9,
            "abundance_log2_a": a,
            "abundance_log2_b": b,
            "delta_log2": d,
            "test": test if test is not None else (TEST_OF.get(cls, "none") if tested else "none"),
            "method": None,
            "null_model": None,
            "q_value": q,
            "lfdr": lfdr,
            "flags": list(flags),
        }
    )
    if tested:
        out["p_value"] = q if cls == "OP" else NAN
    return out


def frame(*rows: Mapping[str, Any], start: int = 1) -> pd.DataFrame:
    """Rows -> a merged pair frame in ``PAIR_COLS`` order, ``pair_id`` ``PAIR_<start + k>``."""
    recs = [dict(r) for r in rows]
    for k, r in enumerate(recs, start=start):
        r["pair_id"] = f"PAIR_{k}"
        r["gene_a"] = r.get("gene_a") or f"G{k}a"
        r["gene_b"] = r.get("gene_b") or f"G{k}b"
    if not recs:
        return pd.DataFrame({c: pd.Series(dtype=object) for c in PAIR_COLS})
    df = pd.DataFrame(recs)
    df["intergenic_bp"] = df["intergenic_bp"].astype("int64")
    df["n_perm"] = pd.array([pd.NA] * len(df), dtype="Int64")
    for c in ("method", "null_model"):
        df[c] = df[c].astype(object)
    return df[list(PAIR_COLS)]


def cfg(
    promotech_yes: float | None = None,
    terminatornet_yes: float | None = None,
    alpha: float | None = None,
) -> GippaConfig:
    """Default config with the three M6-07 knobs optionally overridden. ``terminatornet_yes``
    defaults to 0.5, the threshold the M6-07 examples were written for, not the calibrated default
    (ADR-0022, M6-12 part C)."""
    if terminatornet_yes is None:
        terminatornet_yes = 0.5
    data: dict[str, Any] = {"seed": 1}
    ann: dict[str, Any] = {}
    if promotech_yes is not None:
        ann["promoters"] = {"promotech_yes": promotech_yes}
    if terminatornet_yes is not None:
        ann["terminators"] = {"terminatornet_yes": terminatornet_yes}
    if ann:
        data["annotation"] = ann
    if alpha is not None:
        data["stats"] = {"alpha": {"operon_flags": alpha}}
    return GippaConfig.model_validate(data)


# ------------------------------------------------------------------- M6-05 / M6-06 row frames

DISC_COLS: tuple[str, ...] = ("pair_id", "test", *TEST_COLS, "flags")
MISSED_COLS: tuple[str, ...] = ("pair_id", "test", *TEST_COLS, "which", "flags")


def rows_frame(records: Sequence[Mapping[str, Any]], *, missed: bool = False) -> pd.DataFrame:
    """A ``Discordance.rows`` (or, with ``missed``, ``MissedOperons.rows``) shaped frame.

    Each record names ``pair_id`` and any of the test columns, ``which`` and ``flags``."""
    cols = MISSED_COLS if missed else DISC_COLS
    out = []
    for r in records:
        base: dict[str, Any] = {c: NAN for c in cols}
        base.update({"test": "none", "method": None, "null_model": None, "flags": []})
        if missed:
            base["which"] = ""
        base.update(r)
        out.append(base)
    if not out:
        return pd.DataFrame({c: pd.Series(dtype=object) for c in cols})
    df = pd.DataFrame(out)
    df["n_perm"] = pd.array([pd.NA] * len(df), dtype="Int64")
    return df[list(cols)]


# ------------------------------------------------------------------------- TUs and gene table


def tu(
    tu_id: str, genes: Sequence[str], start: int, end: int, strand: str = "+", contig: str = "c1"
) -> dict[str, Any]:
    """One C5 ``tus.parquet`` row (0-based half-open; ``genes`` in transcription order)."""
    return {
        "tu_id": tu_id,
        "contig": contig,
        "start": start,
        "end": end,
        "strand": strand,
        "genes": list(genes),
        "upstream_start": max(start - 100, 0),
        "upstream_end": start,
    }


def tus(*rows: Mapping[str, Any]) -> pd.DataFrame:
    cols = ["tu_id", "contig", "start", "end", "strand", "genes", "upstream_start", "upstream_end"]
    if not rows:
        return pd.DataFrame({c: pd.Series(dtype=object) for c in cols})
    df = pd.DataFrame([dict(r) for r in rows])
    for c in ("start", "end", "upstream_start", "upstream_end"):
        df[c] = df[c].astype("int64")
    return df[cols]


def gene_table(
    abundance: Mapping[str, float],
    tags: Sequence[str],
    *,
    partial: Sequence[str] = (),
    imputed: Sequence[str] = (),
    detected: Mapping[str, bool] | None = None,
) -> pd.DataFrame:
    """Gene-table rows for ``tags``: detected iff the tag has an abundance (unless overridden)."""
    det = detected or {}
    return pd.DataFrame(
        {
            "locus_tag": list(tags),
            "partial": [t in partial for t in tags],
            "detected": [bool(det.get(t, t in abundance)) for t in tags],
            "abundance_log2": [float(abundance.get(t, NAN)) for t in tags],
            "abundance_imputed": [t in imputed for t in tags],
            "abundance_type": ["ibaq"] * len(tags),
            "abundance_scale": ["molar"] * len(tags),
        }
    )


def chain_pairs(
    tu_genes: Sequence[Sequence[str]],
    abundance: Mapping[str, float],
    *,
    q: Mapping[tuple[str, str], float] | None = None,
    prom: Mapping[tuple[str, str], float] | None = None,
    term: Mapping[tuple[str, str], float] | None = None,
) -> pd.DataFrame:
    """Merged pair frame of consecutive + strand genes across ``tu_genes``: OP inside a TU, NO
    between TUs; ``delta_log2`` = abundance(b) - abundance(a) (NaN unless both detected)."""
    q, prom, term = q or {}, prom or {}, term or {}
    flat = [(g, k) for k, gs in enumerate(tu_genes) for g in gs]
    recs = []
    for (a, ka), (b, kb) in itertools.pairwise(flat):
        cls = "OP" if ka == kb else "NO"
        det = (a in abundance, b in abundance)
        r = row(
            cls,
            prom=prom.get((a, b), 0.1),
            term=term.get((a, b), 0.1),
            q=q.get((a, b), NAN),
            detected=det,
            gene_a=a,
            gene_b=b,
        )
        if all(det):
            r["abundance_log2_a"], r["abundance_log2_b"] = abundance[a], abundance[b]
            r["delta_log2"] = abundance[b] - abundance[a]
        elif det[0]:
            r["abundance_log2_a"] = abundance[a]
        elif det[1]:
            r["abundance_log2_b"] = abundance[b]
        recs.append(r)
    return frame(*recs)


def as_lists(df: pd.DataFrame) -> list[list[str]]:
    """The ``flags`` column as plain lists of str (whatever sequence type the frame holds)."""
    return [
        [str(f) for f in (v if isinstance(v, list | tuple | np.ndarray) else [])]
        for v in df["flags"]
    ]
