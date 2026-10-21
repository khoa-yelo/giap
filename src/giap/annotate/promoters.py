"""Targeted Promotech windows, promoter layer, per-gene promoter features and ``promoters.bed``
(GIAP_SPEC.md §7.1.1, §7.1.4, §6.5; M2-13).

S = first base of the start codon (0-based), E = gene interval end. ``+``: window = genome
``[S + k - 40, S + k)``, ``end3p = S + k``; ``-``: window = revcomp of ``[E - k, E - k + 40)``,
``end3p = E - k + 1`` (``end3p`` is the 1-based coordinate of the window's 3'-most base, in
``[1, length]``). Circular replicons wrap; linear contigs drop windows leaving ``[0, length)``.

In *transcription coordinates* ``tx(p) = p`` on ``+`` and ``-p`` on ``-``, and a gene's anchor is
``A = S`` (``+``) or ``-(E + 1)`` (``-``), so a window of the gene has ``k = tx(end3p) - A``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from giap.annotate.operons import Tu
from giap.io.bed import BedRecord
from giap.io.fasta import write_fasta
from giap.model.coords import Interval, Replicon, interval_sort_key
from giap.model.genes import Gene

WINDOW_NT = 40
K_RANGE = (-200, 20)
MAP_COLUMNS = ["header", "contig", "strand", "end3p", "owner", "k"]
LAYER_COLUMNS = ["contig", "strand", "end3p", "score"]
FEATURE_COLUMNS = ["locus_tag", "promoter_score", "promoter_score_tu", "n_promoter_windows_tu"]
_COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")

_Track = tuple[np.ndarray, np.ndarray]  # sorted end3p, scores (one contig and strand)


@dataclass(frozen=True, slots=True)
class Window:
    """One 40-nt Promotech window, 5'→3' on the gene strand."""

    contig: str
    strand: str
    end3p: int  # 1-based genome coordinate of the window's 3'-most base
    seq: str
    owner: str  # locus tag
    k: int


def check_promoter_mode(mode: str) -> None:
    """Genome-wide Promotech is a valid config but not implemented (M2-13 Q2, revised)."""
    if mode == "genome":
        raise ValueError("annotation.promoters.mode: genome is not implemented in this version")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _slice(seq: str, lo: int, width: int, n: int, circular: bool) -> str | None:
    if not circular:
        return seq[lo : lo + width] if lo >= 0 and lo + width <= n else None
    lo %= n
    if lo + width <= n:
        return seq[lo : lo + width]
    return "".join(seq[(lo + i) % n] for i in range(width))


def gene_windows(
    gene: Gene,
    replicon: Replicon,
    seq: str,
    k_range: tuple[int, int] = K_RANGE,
    width: int = WINDOW_NT,
) -> list[Window]:
    """The windows whose 3' end lies at ``S + k`` for ``k`` in ``k_range`` (k ascending)."""
    iv, n = gene.interval, replicon.length
    if iv.strand not in "+-":
        raise ValueError(f"gene {gene.locus_tag}: strand must be + or -, got {iv.strand!r}")
    circular = replicon.topology == "circular"
    out: list[Window] = []
    for k in range(k_range[0], k_range[1] + 1):
        lo = iv.start + k - width if iv.strand == "+" else iv.end - k
        w = _slice(seq, lo, width, n, circular)
        if w is None:
            continue
        p0 = lo + width - 1 if iv.strand == "+" else lo  # 3'-most base, 0-based
        w = w if iv.strand == "+" else _revcomp(w)
        out.append(Window(iv.contig, iv.strand, p0 % n + 1, w, gene.locus_tag, k))
    return out


def dedupe_and_chunk(
    windows: Iterable[Window], chunk_size: int = 100_000
) -> tuple[list[list[tuple[str, str]]], pd.DataFrame]:
    """Unique sequences as ``w<index>`` (first-seen order) in chunks, plus the full map frame."""
    if chunk_size < 1:
        raise ValueError(f"chunk_size must be >= 1, got {chunk_size}")
    header_of: dict[str, str] = {}
    rows: list[tuple[str, str, str, int, str, int]] = []
    for w in windows:
        h = header_of.setdefault(w.seq, f"w{len(header_of)}")
        rows.append((h, w.contig, w.strand, w.end3p, w.owner, w.k))
    recs = [(h, s) for s, h in header_of.items()]
    chunks = [recs[i : i + chunk_size] for i in range(0, len(recs), chunk_size)]
    wmap = pd.DataFrame(rows, columns=MAP_COLUMNS).astype({"end3p": "int64", "k": "int64"})
    return chunks, wmap


def write_chunks(outdir: Path, chunks: Sequence[Sequence[tuple[str, str]]]) -> list[Path]:
    """Write ``windows_000.fa``, ``windows_001.fa``, …; every sequence must be exactly 40 nt."""
    for chunk in chunks:
        for h, s in chunk:
            if len(s) != WINDOW_NT:
                raise ValueError(f"window {h}: Promotech needs {WINDOW_NT} nt, got {len(s)}")
    outdir.mkdir(parents=True, exist_ok=True)
    paths = [outdir / f"windows_{i:03d}.fa" for i in range(len(chunks))]
    for path, chunk in zip(paths, chunks, strict=True):
        write_fasta(path, chunk, width=WINDOW_NT)
    return paths


def promoter_layer(window_map: pd.DataFrame, scores: Mapping[str, float]) -> pd.DataFrame:
    """One row per unique (contig, strand, end3p) with the score of its window."""
    missing = sorted(set(window_map["header"]) - set(scores))
    if missing:
        raise ValueError(f"no Promotech score for {len(missing)} window(s), e.g. {missing[:3]}")
    uniq = window_map.drop_duplicates(["contig", "strand", "end3p"])
    return pd.DataFrame(
        {
            "contig": uniq["contig"].to_numpy(),
            "strand": uniq["strand"].to_numpy(),
            "end3p": uniq["end3p"].to_numpy(dtype=np.int64),
            "score": np.array([float(scores[h]) for h in uniq["header"]], dtype=float),
        }
    )


# ---------------------------------------------------------------- per-gene features


def _tracks(window_scores: pd.DataFrame) -> dict[tuple[str, str], _Track]:
    out: dict[tuple[str, str], _Track] = {}
    for (contig, strand), sub in window_scores.groupby(["contig", "strand"], sort=True):
        p = sub["end3p"].to_numpy(dtype=np.int64)
        order = np.argsort(p, kind="stable")
        out[(str(contig), str(strand))] = (p[order], sub["score"].to_numpy(dtype=float)[order])
    return out


def _anchor(g: Gene) -> int:
    return g.interval.start if g.interval.strand == "+" else -(g.interval.end + 1)


def _period(contig: str, replicons: Mapping[str, Replicon]) -> int | None:
    rep = replicons.get(contig)
    return rep.length if rep is not None and rep.topology == "circular" else None


def _select(
    tracks: Mapping[tuple[str, str], _Track],
    contig: str,
    strand: str,
    t0: int,
    span: tuple[int, int],
    n: int | None,
) -> np.ndarray:
    """Scores of windows with ``tx(end3p)`` in ``[t0 + a, t0 + b]`` (wrapped when ``n``)."""
    a, b = span
    track = tracks.get((contig, strand))
    if track is None or a > b:
        return np.empty(0)
    lo, hi = (t0 + a, t0 + b) if strand == "+" else (-(t0 + b), -(t0 + a))
    ranges = [(lo, hi)]
    if n is not None:
        lo1 = (lo - 1) % n + 1
        hi1 = lo1 + (hi - lo)
        ranges = [(1, n)] if hi - lo + 1 >= n else [(lo1, min(hi1, n)), (1, hi1 - n)]
    p, s = track
    parts = [s[np.searchsorted(p, x, "left") : np.searchsorted(p, y, "right")] for x, y in ranges]
    return np.concatenate(parts)


def _max(x: np.ndarray) -> float:
    return float(x.max()) if x.size else float("nan")


def gene_promoter_features(
    genes: Sequence[Gene],
    tus: Sequence[Tu],
    window_scores: pd.DataFrame,
    yes: float = 0.5,
    k_range: tuple[int, int] = K_RANGE,
    *,
    replicons: Mapping[str, Replicon],
) -> pd.DataFrame:
    """``promoter_score``, ``promoter_score_tu`` and ``n_promoter_windows_tu`` per gene.

    ``replicons`` is required: circular contigs wrap, linear ones (or unknown contigs) do not.
    """
    tracks = _tracks(window_scores)
    k0, k1 = k_range
    own_span = (0, k1 - k0)

    def own(g: Gene) -> float:
        iv = g.interval
        n = _period(iv.contig, replicons)
        return _max(_select(tracks, iv.contig, iv.strand, _anchor(g) + k0, own_span, n))

    by_tag = {g.locus_tag: g for g in genes}
    first_of = {tag: tu.genes[0] for tu in tus for tag in tu.genes}
    rows: list[tuple[str, float, float, int]] = []
    for g in genes:
        f = by_tag.get(first_of.get(g.locus_tag, g.locus_tag), g)
        iv = g.interval
        n = _period(iv.contig, replicons)
        gap = _anchor(g) - _anchor(f)
        gap = gap % n if n is not None else gap
        tu = _select(tracks, iv.contig, iv.strand, _anchor(f) + k0, (0, gap + k1 - k0), n)
        rows.append((g.locus_tag, own(g), own(f), int((tu >= yes).sum())))
    return pd.DataFrame(rows, columns=FEATURE_COLUMNS).astype({"n_promoter_windows_tu": "int64"})


def promoter_tracks(window_scores: pd.DataFrame) -> dict[tuple[str, str], _Track]:
    """Window scores indexed per (contig, strand), for repeated ``promoter_between_tracks``."""
    return _tracks(window_scores)


def promoter_between(
    up: Gene,
    down: Gene,
    window_scores: pd.DataFrame,
    replicon: Replicon,
    k_range: tuple[int, int] = K_RANGE,
) -> float | None:
    """Max score of ``down``'s windows whose 3' end lies after ``up``'s end; ``None`` if none."""
    return promoter_between_tracks(up, down, _tracks(window_scores), replicon, k_range)


def promoter_between_tracks(
    up: Gene,
    down: Gene,
    tracks: Mapping[tuple[str, str], _Track],
    replicon: Replicon,
    k_range: tuple[int, int] = K_RANGE,
) -> float | None:
    """``promoter_between`` on ``promoter_tracks(window_scores)`` (built once per run)."""
    iv = down.interval
    n = replicon.length if replicon.topology == "circular" else None
    t0 = _anchor(down) + k_range[0]
    bound = up.interval.end + 1 if iv.strand == "+" else -up.interval.start
    a = bound - t0
    if n is not None:
        a = (a + n // 2) % n - n // 2
    got = _select(tracks, iv.contig, iv.strand, t0, (max(a, 0), k_range[1] - k_range[0]), n)
    return float(got.max()) if got.size else None


# ---------------------------------------------------------------- promoters.bed


def promoters_bed(
    window_scores: pd.DataFrame, replicons: Mapping[str, Replicon], yes: float = 0.5
) -> list[BedRecord]:
    """Windows ``>= yes`` merged per contig and strand when overlapping (touching stays apart)."""
    return [b for b, _ in promoter_records(window_scores, replicons, yes)]


def promoter_records(
    window_scores: pd.DataFrame, replicons: Mapping[str, Replicon], yes: float = 0.5
) -> list[tuple[BedRecord, float]]:
    """``promoters_bed`` records paired with the tool's own score (the maximum window score)."""
    keep = window_scores[window_scores["score"] >= yes]
    merged: list[tuple[Interval, float]] = []
    for (contig, strand), sub in keep.groupby(["contig", "strand"], sort=True):
        rep = replicons[str(contig)]
        n, circ = rep.length, rep.topology == "circular"
        starts = [p - WINDOW_NT if strand == "+" else p - 1 for p in sub["end3p"].tolist()]
        wins = sorted(
            zip([s % n if circ else s for s in starts], sub["score"].tolist(), strict=True)
        )
        runs: list[tuple[int, int, float]] = []
        for lo, x in wins:
            if runs and lo < runs[-1][1]:
                r0, r1, rx = runs[-1]
                runs[-1] = (r0, max(r1, lo + WINDOW_NT), max(rx, x))
            else:
                runs.append((lo, lo + WINDOW_NT, x))
        if circ and len(runs) > 1 and runs[-1][1] - n > runs[0][0]:  # overlap across the origin
            _, f1, fx = runs.pop(0)
            r0, r1, rx = runs[-1]
            runs[-1] = (r0, max(r1, f1 + n), max(rx, fx))
        merged += [(Interval(str(contig), lo, hi, strand), x) for lo, hi, x in runs]
    order = {cid: i for i, cid in enumerate(replicons)}
    merged.sort(key=lambda m: interval_sort_key(m[0], order))
    return [
        (BedRecord(iv, f"promoter_{i}", round(1000 * float(x))), float(x))
        for i, (iv, x) in enumerate(merged, start=1)
    ]
