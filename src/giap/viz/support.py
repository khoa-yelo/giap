"""Plot support tables from the run's inputs (M9-03; ADR-0043; spec §10.2 V1, V2).

``replicon_landmarks`` places the oriC and the terminus of every replicon, and ``gc_skew`` tiles
every replicon with GC-skew windows (M9-03); ``variant_windows`` cuts the reference protein
around each V12 top ALT site (M9-18). All compute in 0-based half-open coordinates (protein
positions 1-based) and return the table **as written** (1-based positions, 1-based closed
windows), converting once where the frame is built (lead decision 3);
``giap.io.plot_support.write_support`` writes them unchanged.

- **oriC.** ``giap.context.oriter.oric_position`` on the contig (the ori-ter analysis's own rule,
  so the V1/V2 marks and the replichores agree). A non-integral midpoint (an odd-length call,
  or the circular mean of several, with its float error) goes to the nearest base, ties up:
  ``floor(o + 1/2) mod L`` (0-based), written ``+ 1``. The terminus is ``(o0 + L // 2) mod L``.
- **GC skew.** (G - C) / (G + C) per window, case-insensitive, other letters not counted; NA when
  G + C = 0. The cumulative skew is the running sum along the contig, an NA window adding 0
  (decision 5). Windows of ``window_bp`` from position 1, the last shorter, never wrapping.
- **Variant windows.** The gene's protein (``genes.parquet`` ``aa``, transcription order on both
  strands) from ``aa_pos - VARIANT_FLANK_AA`` to ``aa_pos + VARIANT_FLANK_AA``, clipped to the
  protein, for the at most ``MAX_TOP_SITES`` ``call = ALT`` sites of ``top_alt_sites`` (M9-18
  lead decisions 1, 2). A gene absent from ``genes`` or with an empty ``aa`` keeps its row with
  NA ``window_start_aa`` and ``ref_window`` (lead decision 6).

Pure: no file I/O.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np
import pandas as pd

from giap.context.oriter import oric_position

__all__ = [
    "LANDMARK_COLUMNS", "MAX_TOP_SITES", "SKEW_COLUMNS", "VARIANT_FLANK_AA",
    "VARIANT_WINDOW_COLUMNS", "WINDOW_BP", "class_error_value", "gc_skew", "replicon_landmarks",
    "top_alt_sites", "variant_windows",
]  # fmt: skip

#: GC-skew window (lead decision 1: the V1 peptide-density bin; a constant, not a config key)
WINDOW_BP = 10_000
LANDMARK_COLUMNS: tuple[str, ...] = (
    "contig", "length_bp", "topology", "oric_position", "ter_position", "status",
)  # fmt: skip
SKEW_COLUMNS: tuple[str, ...] = ("contig", "start", "end", "gc_skew", "cumulative_gc_skew")
#: residues shown each side of a V12 site (M9-18 lead decision 2: a display choice)
VARIANT_FLANK_AA = 10
#: V12 panel (b) sites (§10.2 V12: "up to 30 top ALT sites")
MAX_TOP_SITES = 30
VARIANT_WINDOW_COLUMNS: tuple[str, ...] = (
    "variant_id", "gene", "aa_pos", "window_start_aa", "ref_window",
)  # fmt: skip


def _oric0(features: pd.DataFrame | None, contig: str, length: int) -> int | None:
    """0-based oriC base of a circular replicon, or None (``no_oric``)."""
    o = oric_position(features, contig, length)
    return None if o is None else math.floor(o + 0.5) % length


def replicon_landmarks(replicons: pd.DataFrame, features: pd.DataFrame | None) -> pd.DataFrame:
    """One row per replicon (C1 order): length, topology, 1-based oriC and terminus, status.

    Parameters
    ----------
    replicons : pandas.DataFrame
        C1 ``replicons`` (``contig``, ``length``, ``topology``).
    features : pandas.DataFrame or None
        C15 ``features`` (0-based half-open ``oric`` rows); None behaves like no oriC call.

    Returns
    -------
    pandas.DataFrame
        ``LANDMARK_COLUMNS``; positions ``Int64`` with NA on linear and ``no_oric`` replicons.
    """
    rows: list[tuple[str, int, str, int | None, int | None, str]] = []
    for contig, length, topology in zip(
        replicons["contig"], replicons["length"], replicons["topology"], strict=True
    ):
        c, n, topo = str(contig), int(length), str(topology)
        o0 = _oric0(features, c, n) if topo == "circular" else None
        if o0 is None:
            rows.append((c, n, topo, None, None, "no_oric" if topo == "circular" else "linear"))
        else:
            rows.append((c, n, topo, o0 + 1, (o0 + n // 2) % n + 1, "called"))
    out = pd.DataFrame(rows, columns=list(LANDMARK_COLUMNS))
    return out.astype({"length_bp": "int64", "oric_position": "Int64", "ter_position": "Int64"})


def _windows(seq: str, window_bp: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """0-based half-open window starts and ends, and the window skew and cumulative skew."""
    n = len(seq)
    starts = np.arange(0, n, window_bp, dtype=np.int64)
    ends = np.minimum(starts + window_bp, n)
    codes = np.frombuffer(seq.upper().encode("ascii", errors="replace"), dtype=np.uint8)
    g = np.add.reduceat((codes == ord("G")).astype(np.int64), starts)
    c = np.add.reduceat((codes == ord("C")).astype(np.int64), starts)
    total = g + c
    with np.errstate(invalid="ignore", divide="ignore"):
        skew = np.where(total > 0, (g - c) / np.where(total > 0, total, 1), np.nan)
    return starts, ends, skew, np.cumsum(np.nan_to_num(skew, nan=0.0))


def gc_skew(genome: Mapping[str, str], replicons: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    """GC-skew windows of every replicon (C1 order, then ``start``), 1-based closed.

    Parameters
    ----------
    genome : Mapping[str, str]
        Contig -> sequence; every replicon's contig must be present.
    replicons : pandas.DataFrame
        C1 ``replicons`` (``contig``), giving the row order; linear replicons included.
    window_bp : int
        Window width (``WINDOW_BP`` in a run).

    Returns
    -------
    pandas.DataFrame
        ``SKEW_COLUMNS``; ``gc_skew`` NaN when the window holds no G or C.
    """
    if window_bp < 1:
        raise ValueError(f"window_bp must be positive, got {window_bp}")
    parts: list[pd.DataFrame] = []
    for contig in replicons["contig"].astype(str):
        if contig not in genome:
            raise ValueError(f"replicon {contig!r} is not in the genome")
        seq = genome[contig]
        if not seq:
            continue
        starts, ends, skew, cum = _windows(seq, window_bp)
        cols = {"contig": contig, "start": starts + 1, "end": ends,
                "gc_skew": skew, "cumulative_gc_skew": cum}  # fmt: skip
        parts.append(pd.DataFrame(cols))
    if not parts:
        empty = {"contig": "object", "start": "int64", "end": "int64",
                 "gc_skew": "float64", "cumulative_gc_skew": "float64"}  # fmt: skip
        return pd.DataFrame({k: pd.Series([], dtype=v) for k, v in empty.items()})
    return pd.concat(parts, ignore_index=True)


def class_error_value(value: object) -> float | None:
    """A ``variants.tsv`` ``class_error_rate`` cell as a number; ``None`` for NA, ``unestimated``
    or any other non-numeric text."""
    if value is None or value is pd.NA:
        return None
    try:
        rate = float(str(value))
    except ValueError:
        return None
    return rate if math.isfinite(rate) else None


def top_alt_sites(variants: pd.DataFrame, limit: int = MAX_TOP_SITES) -> pd.DataFrame:
    """The ``call = ALT`` rows in the V12 top-site order, at most ``limit`` (lead decision 1):
    ``confidence`` high first, ``n_alt`` descending, ``class_error_rate`` ascending (NA last),
    then ``variant_id``. A copy with a fresh index; ``variants`` is not modified."""
    alt = variants[variants["call"].astype(str) == "ALT"]
    rates = [class_error_value(v) for v in alt["class_error_rate"]]
    keys = [
        (0 if str(conf) == "high" else 1, -int(n), r is None, 0.0 if r is None else r, str(vid))
        for conf, n, r, vid in zip(
            alt["confidence"], alt["n_alt"], rates, alt["variant_id"], strict=True
        )
    ]
    order = sorted(range(len(alt)), key=keys.__getitem__)[:limit]
    return alt.iloc[order].reset_index(drop=True)


def _proteins(genes: pd.DataFrame) -> dict[str, str]:
    """``locus_tag`` -> its first non-empty ``aa``."""
    out: dict[str, str] = {}
    for tag, aa in zip(genes["locus_tag"], genes["aa"], strict=True):
        if isinstance(aa, str) and aa and str(tag) not in out:
            out[str(tag)] = aa
    return out


def variant_windows(variants: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """The reference protein window of every V12 top ALT site, as written (1-based).

    Parameters
    ----------
    variants : pandas.DataFrame
        ``variants.tsv`` as ``RunView.table("variants")`` reads it.
    genes : pandas.DataFrame
        C2 ``genes.parquet`` (``locus_tag``, ``aa``).

    Returns
    -------
    pandas.DataFrame
        ``VARIANT_WINDOW_COLUMNS`` in top-site order; ``aa_pos`` / ``window_start_aa`` ``Int64``,
        ``window_start_aa`` and ``ref_window`` NA when the protein is unknown (lead decision 6).
    """
    proteins = _proteins(genes)
    rows: list[tuple[str, str, int, int | None, str | None]] = []
    for _, site in top_alt_sites(variants).iterrows():
        gene, pos = str(site["gene"]), int(site["aa_pos"])
        aa = proteins.get(gene, "")
        lo, hi = max(1, pos - VARIANT_FLANK_AA), min(len(aa), pos + VARIANT_FLANK_AA)
        if lo > hi:  # protein unknown (or shorter than aa_pos): kept, unaligned
            rows.append((str(site["variant_id"]), gene, pos, None, None))
        else:
            rows.append((str(site["variant_id"]), gene, pos, lo, aa[lo - 1 : hi]))
    out = pd.DataFrame(rows, columns=list(VARIANT_WINDOW_COLUMNS))
    return out.astype({"variant_id": object, "gene": object, "aa_pos": "Int64",
                       "window_start_aa": "Int64", "ref_window": object})  # fmt: skip
