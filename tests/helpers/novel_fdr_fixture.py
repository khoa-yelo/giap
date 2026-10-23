"""Builders for the M7-11 class-FDR tests (docs/specs/M7-11-novel-class-fdr.md).

Conventions pinned by the test author (the spec leaves the frame layouts open):

* ``peptides`` is the aggregated C9 frame (``giap.io.peptides.AGG_COLUMNS``): ``protein_ids`` is a
  list of IDs, decoy rows have ``is_decoy`` true and prefixed IDs (``rev_GIAPORF_...``).
* ``id_map`` is ``db_id_map.tsv`` as a frame (``giap.makedb.build.ID_MAP_COLUMNS``): 1-based
  closed ``start``/``end``, ``kind = orf``, ``locus_tag`` set for ``nterm_extension`` entries only.
* A plus-strand makedb ORF ``[s, e)`` (0-based, stop-to-stop, stops excluded) is closed by the stop
  codon at ``e``; a minus-strand one by the stop codon at ``s - 3``. ``Orf.stop`` is that lowest
  coordinate, so the ORF of ``GIAPORF_cds_outframe_c1_101_400_p`` has ``stop = 400``.
* Picked scores are encoded by a rank (1 = best): input ``q_value = rank * 1e-5`` (all <= 0.01, so
  no ``q_value_max`` filter can drop a row) and ``score = 1000 - rank`` (consistent with q).
* De novo hit frames carry ``giap.mapping.hits.HIT_COLUMNS`` plus the M3-10 ranking score E:
  ``e_target`` on target rows (``search = target``) and ``e_decoy`` on decoy rows
  (``search = decoy``, the mirrored target interval in ``contig/start/end/strand``), one decoy row
  per peptide (its best decoy cost), ``n_diff = 0``.
* ``classify`` stubs: ``ContigClassifier`` maps a contig name to a locus class (duck-typed
  ``LocusClassifier``) for the property and calibration tests; the unit examples use the real one.
* An unestimated ORF is absent from the returned dict or maps to NaN (``q_of``).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from giap.mapping.hits import HIT_COLUMNS
from giap.model.coords import Interval
from giap.novel.candidates import Orf

ALPHA = 0.01  # stats.alpha.novel_orfs
LETTERS = "ACDEFGHKMNPQRSTVWY"  # no I/L: I/L folding never changes a name
ID_MAP_COLUMNS = ("db_id", "kind", "locus_tag", "class", "contig", "start", "end", "strand")


def q_of(qs: Mapping[str, float], orf_id: str) -> float:
    """Class q of ``orf_id``; NaN when unestimated (absent or NaN)."""
    v = qs.get(orf_id)
    return math.nan if v is None else float(v)


def callable_at(qs: Mapping[str, float], orf_id: str, alpha: float = ALPHA) -> bool:
    q = q_of(qs, orf_id)
    return not math.isnan(q) and q <= alpha


def pep_name(k: int, head: str = "PEPT") -> str:
    """Distinct peptide for counter ``k`` (no I or L): ``head`` + 8 coded residues."""
    digits = []
    for _ in range(8):
        digits.append(LETTERS[k % len(LETTERS)])
        k //= len(LETTERS)
    return head + "".join(digits)


def make_orf(
    orf_id: str,
    cls: str,
    contig: str,
    strand: str,
    stop: int,
    peptides: Sequence[str],
    *,
    overlapping_gene: str | None = None,
    best_q: float | None = None,
) -> Orf:
    start = stop - 300 if strand == "+" else stop + 300
    return Orf(
        orf_id=orf_id,
        cls=cls,  # type: ignore[arg-type]
        contig=contig,
        strand=strand,  # type: ignore[arg-type]
        frame=0,
        stop=stop,
        starts=(start,),
        start_longest=start,
        end_5p_open=False,
        end_3p_open=False,
        peptides=tuple(sorted({p.upper().replace("I", "L") for p in peptides})),
        best_q=best_q,
        overlapping_gene=overlapping_gene,
        first_peptide=start,
    )


def db_id(cls: str, contig: str, start0: int, end0: int, strand: str) -> str:
    return f"GIAPORF_{cls}_{contig}_{start0 + 1}_{end0}_{'p' if strand == '+' else 'm'}"


def pep_row(
    peptide: str,
    ids: Sequence[str],
    rank: float,
    *,
    is_decoy: bool,
    source: str = "database",
    q: float | None = None,
    score: float | None = None,
) -> dict[str, Any]:
    return {
        "peptide": peptide,
        "q_value": rank * 1e-5 if q is None else q,
        "pep": math.nan,
        "score": 1000.0 - rank if score is None else score,
        "is_decoy": is_decoy,
        "is_contaminant": False,
        "protein_ids": list(ids),
        "modified_peptides": [],
        "n_rows": 1,
        "n_runs": 1,
        "intensity": 1.0,
        "source": source,
        "search_db": "same",
    }


def peptides_frame(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    cols = ["peptide", "q_value", "pep", "score", "is_decoy", "is_contaminant", "protein_ids",
            "modified_peptides", "n_rows", "n_runs", "intensity", "source",
            "search_db"]  # fmt: skip
    df = pd.DataFrame(list(rows), columns=cols)
    df["is_decoy"] = df["is_decoy"].astype(bool)
    df["is_contaminant"] = df["is_contaminant"].astype(bool)
    for c in ("q_value", "pep", "score", "intensity"):
        df[c] = df[c].astype("float64")
    return df


@dataclass
class PickedWorld:
    """A makedb search: entries in ``id_map``, target ORFs and decoy rows in C9."""

    orfs: list[Orf] = field(default_factory=list)
    id_rows: list[tuple[Any, ...]] = field(default_factory=list)
    pep_rows: list[dict[str, Any]] = field(default_factory=list)
    _k: int = 0
    _slot: dict[str, int] = field(default_factory=dict)

    def _pep(self) -> str:
        self._k += 1
        return pep_name(self._k)

    def entry(
        self,
        cls: str,
        contig: str,
        strand: str = "+",
        *,
        locus_tag: str | None = None,
        in_map: bool = True,
        start0: int | None = None,
        length: int = 300,
    ) -> tuple[str, int]:
        """A makedb ORF entry; returns ``(db_id, stop)`` with ``stop`` the closing stop codon's
        lowest coordinate (``Orf.stop`` of the matching ORF)."""
        if start0 is None:
            i = self._slot.get(contig, 0)
            self._slot[contig] = i + 1
            start0 = 1000 + 600 * i
        end0 = start0 + length
        did = db_id(cls, contig, start0, end0, strand)
        if in_map:
            self.id_rows.append((did, "orf", locus_tag, cls, contig, start0 + 1, end0, strand))
        return did, (end0 if strand == "+" else start0 - 3)

    def target(
        self,
        cls: str,
        contig: str,
        rank: float | Sequence[float],
        *,
        strand: str = "+",
        entry: tuple[str, int] | None = None,
        locus_tag: str | None = None,
        stop: int | None = None,
        q: float | None = None,
        score: float | None = None,
    ) -> tuple[Orf, str]:
        """A target ORF (one peptide per rank) matched to ``entry`` (a new one by default);
        returns ``(orf, db_id)``."""
        did, est = entry or self.entry(cls, contig, strand, locus_tag=locus_tag)
        ranks = [rank] if isinstance(rank, int | float) else list(rank)
        peps = []
        for r in ranks:
            p = self._pep()
            peps.append(p)
            self.pep_rows.append(pep_row(p, [did], r, is_decoy=False, q=q, score=score))
        st = est if stop is None else stop
        orf = make_orf(f"NOVEL_{contig}_{len(self.orfs) + 1}_{strand}", cls, contig, strand, st,
                       peps, overlapping_gene=locus_tag)  # fmt: skip
        self.orfs.append(orf)
        return orf, did

    def decoy(
        self,
        did: str,
        rank: float,
        *,
        prefix: str = "rev_",
        q: float | None = None,
        score: float | None = None,
    ) -> str:
        self._k += 1
        p = pep_name(self._k, "DECY")
        self.pep_rows.append(pep_row(p, [prefix + did], rank, is_decoy=True, q=q, score=score))
        return p

    @property
    def peptides(self) -> pd.DataFrame:
        return peptides_frame(self.pep_rows)

    @property
    def id_map(self) -> pd.DataFrame:
        df = pd.DataFrame(self.id_rows, columns=list(ID_MAP_COLUMNS))
        return df.astype({"start": "int64", "end": "int64"})


# ---------------------------------------------------------------- de novo


class ContigClassifier:
    """Duck-typed ``LocusClassifier``: the class of an interval is that of its contig."""

    def __init__(self, classes: Mapping[str, str]) -> None:
        self.classes = dict(classes)
        self.calls = 0

    def classify(self, iv: Interval) -> tuple[str, str | None]:
        self.calls += 1
        return self.classes[iv.contig], None


def hit_rows(
    peptides: Sequence[str],
    intervals: Sequence[Interval],
    e: Sequence[float],
    *,
    search: str,
) -> pd.DataFrame:
    """Hit rows (``HIT_COLUMNS`` + ``e_target`` / ``e_decoy``) for one search."""
    rows = []
    for p, iv in zip(peptides, intervals, strict=True):
        rows.append({
            "peptide": p, "query": p.replace("I", "L"), "search": search,
            "decoy_type": "shuffle" if search == "decoy" else None, "contig": iv.contig,
            "frame_id": 0, "aa_start": 0, "aa_end": len(p), "start": iv.start, "end": iv.end,
            "strand": iv.strand, "frame": 0, "route": "exact", "n_diff": 0, "diff_types": "",
            "blocks": None,
        })  # fmt: skip
    df = pd.DataFrame(rows, columns=list(HIT_COLUMNS))
    for c in ("frame_id", "aa_start", "aa_end", "start", "end", "frame", "n_diff"):
        df[c] = df[c].astype("int64")
    df["e_target" if search == "target" else "e_decoy"] = np.asarray(e, dtype=np.float64)
    return df


def iv_on(contig: str, k: int, strand: str = "-") -> Interval:
    """A 30-nt interval on ``contig`` (slot ``k``)."""
    return Interval(contig, 100 + 60 * k, 130 + 60 * k, strand)  # type: ignore[arg-type]
