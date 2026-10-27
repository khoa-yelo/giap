"""Promoter Calculator regions, CSV parser, ``tx_rates`` layer and TU transcription features
(GIPPA_SPEC.md §7.1.1, §7.1.4, §9.3; M2-14).

The tool joins every FASTA line into one sequence, so each call gets exactly one region on one
sequence line, 5'→3' in TU orientation. Its ``TSS`` column is the 0-based +1 base for ``+`` rows and
the 1-based forward position for ``-`` rows; ``PcHit.tss`` is 1-based forward (region coordinates).
Only ``+`` hits (the TU's strand) with the TSS inside the region are kept. Region ``[a, b)`` maps
a ``+`` TU's TSS to ``a + tss - 1`` and a ``-`` TU's to ``b - tss``, reduced modulo the replicon
length on circular replicons.
"""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from gippa.annotate.operons import Tu, region_sequence
from gippa.model.coords import Replicon
from gippa.model.genes import Gene

TX_COLUMNS = ("contig", "strand", "tss", "tx_rate", "tu_id")
FEATURE_COLUMNS = ("tu_id", "tx_rate_max", "tx_rate_sum", "tss_distance_bp")
_NEEDED = ("TSS", "strand", "Tx_rate")


@dataclass(frozen=True, slots=True)
class PcHit:
    """One Promoter Calculator prediction; ``tss`` 1-based in forward region coordinates."""

    tss: int
    strand: str
    tx_rate: float


def region_fasta_text(tu: Tu, genome: Mapping[str, str]) -> str:
    """``>{tu_id}\\n{seq}\\n``: one record, one sequence line, in TU orientation."""
    seq = region_sequence(tu.upstream, genome[tu.contig])
    text = f">{tu.tu_id}\n{seq}\n"
    assert text.count(">") == 1 and text.count("\n") == 2, tu.tu_id
    return text


def parse_pc_csv(path: Path) -> list[PcHit]:
    """Tool CSV (CRLF or LF; header-only or empty = no hits) → hits in file order."""
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            return []
        for col in _NEEDED:
            if col not in reader.fieldnames:
                raise ValueError(f"{path}: Promoter Calculator CSV lacks column {col!r}")
        hits = []
        for r in reader:
            strand = r["strand"].strip()
            tss = int(r["TSS"]) + (1 if strand == "+" else 0)
            hits.append(PcHit(tss=tss, strand=strand, tx_rate=float(r["Tx_rate"])))
    return hits


def _genome_tss(tu: Tu, tss: int, replicons: Mapping[str, Replicon]) -> int:
    g = tu.upstream.start + tss - 1 if tu.strand == "+" else tu.upstream.end - tss
    rep = replicons.get(tu.contig)
    return g % rep.length if rep is not None and rep.topology == "circular" else g


def tx_layer(
    tus: Sequence[Tu], hits: Mapping[str, Sequence[PcHit]], replicons: Mapping[str, Replicon]
) -> pd.DataFrame:
    """Kept hits as genome TSSs: ``contig, strand, tss`` (0-based), ``tx_rate``, ``tu_id``."""
    rows = []
    for tu in tus:
        n = tu.upstream.end - tu.upstream.start
        for h in hits.get(tu.tu_id, ()):
            if h.strand == "+" and 1 <= h.tss <= n:
                g = _genome_tss(tu, h.tss, replicons)
                rows.append((tu.contig, tu.strand, g, float(h.tx_rate), tu.tu_id))
    df = pd.DataFrame(rows, columns=list(TX_COLUMNS))
    return df.astype({"contig": object, "strand": object, "tss": "int64", "tx_rate": "float64"})


def _signed(d: int, length: int | None) -> int:
    """``d`` reduced to (-length/2, length/2] on a ring; unchanged without a length."""
    if length is None:
        return d
    d %= length
    return d - length if d > length // 2 else d


def tu_tx_features(
    tus: Sequence[Tu], tx: pd.DataFrame, replicons: Mapping[str, Replicon]
) -> pd.DataFrame:
    """One row per TU (TU order): max / sum of kept Tx rates (0 without a hit) and the distance
    from the max-rate TSS to the first start codon (positive = upstream; NA without a hit)."""
    groups = {str(k): g for k, g in tx.groupby("tu_id", sort=False)} if len(tx) else {}
    rows: list[tuple[str, float, float, int | None]] = []
    for tu in tus:
        sub = groups.get(tu.tu_id)
        if sub is None or sub.empty:
            rows.append((tu.tu_id, 0.0, 0.0, None))
            continue
        rates = sub["tx_rate"].to_numpy(dtype=float)
        tss = int(sub["tss"].iloc[int(rates.argmax())])
        rep = replicons.get(tu.contig)
        length = rep.length if rep is not None and rep.topology == "circular" else None
        d = tu.interval.start - tss if tu.strand == "+" else tss - (tu.interval.end - 1)
        rows.append((tu.tu_id, float(rates.max()), float(rates.sum()), _signed(d, length)))
    df = pd.DataFrame(rows, columns=list(FEATURE_COLUMNS))
    dtypes = {"tx_rate_max": "float64", "tx_rate_sum": "float64", "tss_distance_bp": "Int64"}
    return df.astype(dtypes)


def _tx_gap(up: Gene, down: Gene, L: int) -> tuple[int, int]:
    """Origin ``lo`` and length ``span`` of the gap in transcription order (across the origin)."""
    if up.interval.strand == "+":
        return up.interval.end, _signed(down.interval.start - up.interval.end, L)
    return up.interval.start - 1, _signed(up.interval.start - down.interval.end, L)


def tx_rate_between(up: Gene, down: Gene, tx: pd.DataFrame, replicon: Replicon) -> float | None:
    """Max ``tx_rate`` of TSSs on the pair's strand strictly between ``up``'s 3' end and ``down``'s
    start codon (transcription order, across the origin); ``None`` when there is none."""
    if _tx_gap(up, down, replicon.length)[1] <= 0:
        return None
    sel = tx[(tx["contig"] == replicon.id) & (tx["strand"] == up.interval.strand)]
    return tx_rate_between_arrays(
        up, down, sel["tss"].to_numpy(), sel["tx_rate"].to_numpy(), replicon
    )


def tx_rate_between_arrays(
    up: Gene, down: Gene, tss: np.ndarray, tx_rate: np.ndarray, replicon: Replicon
) -> float | None:
    """``tx_rate_between`` on the ``tss`` / ``tx_rate`` arrays of the pair's contig and strand."""
    L = replicon.length
    lo, span = _tx_gap(up, down, L)
    if span <= 0:
        return None
    off = (tss - lo) % L if up.interval.strand == "+" else (lo - tss) % L
    rates = tx_rate[off < span]
    return float(rates.max()) if len(rates) else None
