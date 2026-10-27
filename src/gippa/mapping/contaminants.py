"""Unmapped bucket and contaminants (GIPPA_SPEC.md §7.2.8, §9.3 ``mapping/unmapped.tsv``; M3-15).

Sequences that land nowhere are matched against the universal contaminant library: peptides by
exact substring with I = L (Aho-Corasick), proteins by ``blastp`` (>= 90 % identity over >= 50 %
of the query, per HSP, both inclusive). A library match labels a row ``contaminant`` and keeps the
mapping reason; otherwise a row whose engine protein IDs all carry a contaminant prefix is a
``contaminant`` with reason ``contaminant_id``; everything else is ``unmapped``.

Input frames (M3-15 lead decisions):

- ``peptides``: aggregated query peptides (``peptide``, ``protein_ids`` = list of strings);
- ``placements``: M3-10 per-peptide frame; presence = a target hit, ``accepted`` = decision;
- ``protein_state``: ``protein_id``, ``sequence``, ``assigned``, ``has_alignment``.

Engine decoys, peptides with ``q_value > q_value_max``, accepted placements and assigned proteins
never enter the bucket. Rows are sorted by ``kind`` then ``sequence_id``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import ahocorasick
import pandas as pd

from gippa.config.model import InputConfig
from gippa.io.protein_fasta import ProteinRecord, resolve_protein_id
from gippa.mapping.tblastn import OUTFMT

Reason = Literal[
    "no_hit", "rejected_placement_fdr", "below_identity", "too_short", "contaminant_id"
]
Kind = Literal["peptide", "protein"]
Label = Literal["contaminant", "unmapped"]

_N_COLS = len(OUTFMT.split()) - 1
_COL = {name: i for i, name in enumerate(OUTFMT.split()[1:])}
_DEFAULT_Q_VALUE_MAX: float = InputConfig().q_value_max


@dataclass(frozen=True, slots=True)
class UnmappedRow:
    """One row of ``mapping/unmapped.tsv`` (§9.3)."""

    sequence_id: str
    sequence: str
    kind: Kind
    label: Label
    reason: Reason
    contaminant_match: str | None


def peptide_reason(
    has_target_hit: bool, accepted: bool, length: int, tolerant_min_len: int
) -> Reason | None:
    """Mapping reason of a peptide; ``None`` when placed and accepted (never in the bucket).

    ``rejected_placement_fdr``: a target hit that failed placement FDR; ``too_short``: no target
    hit and shorter than the tolerant-search minimum; ``no_hit``: searched, nothing found.
    """
    if has_target_hit:
        return None if accepted else "rejected_placement_fdr"
    return "too_short" if length < tolerant_min_len else "no_hit"


def _il(seq: str) -> str:
    return seq.upper().replace("I", "L")


def match_contaminant_peptides(
    peptides: Sequence[str], library: Sequence[ProteinRecord]
) -> dict[str, str]:
    """Peptide -> lowest library ID whose sequence contains it (I = L); unmatched absent."""
    keys: dict[str, list[str]] = {}
    for pep in peptides:
        if pep:
            keys.setdefault(_il(pep), []).append(pep)
    if not keys or not library:
        return {}
    automaton = ahocorasick.Automaton(ahocorasick.STORE_ANY)
    for key in keys:
        automaton.add_word(key, key)
    automaton.make_automaton()
    best: dict[str, str] = {}
    for rec in library:
        for _, key in automaton.iter(_il(rec.sequence)):
            cur = best.get(key)
            if cur is None or rec.id < cur:
                best[key] = rec.id
    return {pep: best[key] for key, peps in keys.items() if key in best for pep in peps}


def blastp_cmd(query: Path, library: Path, threads: int) -> list[str]:
    """``blastp -query Q -subject LIB -num_threads N -seg no -outfmt OUTFMT``.

    The library FASTA is searched directly (``-subject``: a few hundred entries, no database
    build); the tabular columns are tblastn's 15 (``gippa.mapping.tblastn.OUTFMT``).
    """
    return [
        "blastp",
        "-query",
        str(query),
        "-subject",
        str(library),
        "-num_threads",
        str(threads),
        "-seg",
        "no",
        "-outfmt",
        OUTFMT,
    ]


@dataclass(frozen=True, slots=True)
class _BlastpHit:
    qseqid: str
    sseqid: str
    identity: float
    coverage: float
    bitscore: float

    def key(self) -> tuple[float, float, float, str]:
        """Lead decision 2: highest bitscore, identity, coverage, then lowest subject ID."""
        return (-self.bitscore, -self.identity, -self.coverage, self.sseqid)


def _parse_blastp_row(cols: list[str], where: str) -> _BlastpHit:
    if len(cols) != _N_COLS:
        raise ValueError(f"{where}: {len(cols)} columns, expected {_N_COLS} (outfmt {OUTFMT!r})")
    try:
        qstart, qend = int(cols[_COL["qstart"]]), int(cols[_COL["qend"]])
        qlen = int(cols[_COL["qlen"]])
        if qlen <= 0:
            raise ValueError(f"qlen {qlen} <= 0")
        return _BlastpHit(
            qseqid=resolve_protein_id(cols[_COL["qseqid"]])[0],
            sseqid=resolve_protein_id(cols[_COL["sseqid"]])[0],
            identity=float(cols[_COL["pident"]]) / 100.0,
            coverage=(abs(qend - qstart) + 1) / qlen,
            bitscore=float(cols[_COL["bitscore"]]),
        )
    except ValueError as exc:
        raise ValueError(f"{where}: {exc}") from exc


def _read_blastp(path: Path) -> list[_BlastpHit]:
    out: list[_BlastpHit] = []
    with path.open(encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n\r")
            if not line.strip() or line.startswith("#"):
                continue
            out.append(_parse_blastp_row(line.split("\t"), f"{path.name} line {lineno}"))
    return out


def _best_hits(
    hits: Iterable[_BlastpHit], min_identity: float, min_coverage: float
) -> dict[str, str]:
    """Query -> subject of its best passing HSP (thresholds inclusive, per HSP)."""
    best: dict[str, _BlastpHit] = {}
    for h in hits:
        if h.identity < min_identity or h.coverage < min_coverage:
            continue
        cur = best.get(h.qseqid)
        if cur is None or h.key() < cur.key():
            best[h.qseqid] = h
    return {q: h.sseqid for q, h in sorted(best.items())}


def match_contaminant_proteins(
    blastp_tsv: Path, min_identity: float = 0.90, min_coverage: float = 0.50
) -> dict[str, str]:
    """Protein ID -> library ID from a tabular blastp file (``OUTFMT``); unmatched absent.

    Identity = ``pident / 100`` and coverage = ``(qend - qstart + 1) / qlen`` of one HSP.
    Raises ``ValueError`` naming the line on a malformed row.
    """
    return _best_hits(_read_blastp(blastp_tsv), min_identity, min_coverage)


def _prefixed_id(ids: Sequence[str], prefixes: tuple[str, ...]) -> str | None:
    """The first ID when every ID carries a contaminant prefix, else ``None``."""
    if not ids or not prefixes:
        return None
    return ids[0] if all(i.startswith(prefixes) for i in ids) else None


def _as_ids(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and value != value):
        return []
    if isinstance(value, str):
        return [x for x in value.split(";") if x]
    if isinstance(value, Iterable):
        return [str(x) for x in value if str(x)]
    return [str(value)]


def _label(
    seq_id: str,
    sequence: str,
    kind: Kind,
    reason: Reason,
    match: str | None,
    prefixed: str | None,
) -> UnmappedRow:
    if match is not None:
        return UnmappedRow(seq_id, sequence, kind, "contaminant", reason, match)
    if prefixed is not None:
        return UnmappedRow(seq_id, sequence, kind, "contaminant", "contaminant_id", prefixed)
    return UnmappedRow(seq_id, sequence, kind, "unmapped", reason, None)


def _placement_state(placements: pd.DataFrame) -> dict[str, bool]:
    """Peptide -> accepted (any accepted row) for every peptide with a target hit."""
    state: dict[str, bool] = {}
    for pep, acc in zip(placements["peptide"], placements["accepted"], strict=True):
        state[str(pep)] = state.get(str(pep), False) or bool(acc)
    return state


def _searched(peptides: pd.DataFrame, q_value_max: float) -> pd.DataFrame:
    """Rows that were searched: not an engine decoy and ``q_value <= q_value_max`` or no q.

    M3-15 review decision 1; a missing q-value (NaN) is not above the cut-off.
    """
    decoy = peptides["is_decoy"].astype(bool)
    q = pd.to_numeric(peptides["q_value"], errors="coerce")
    return peptides[~decoy & ~(q > q_value_max)]


def _peptide_rows(
    peptides: pd.DataFrame,
    placements: pd.DataFrame,
    pep_matches: Mapping[str, str],
    prefixes: tuple[str, ...],
    tolerant_min_len: int,
    q_value_max: float,
) -> list[UnmappedRow]:
    state = _placement_state(placements)
    searched = _searched(peptides, q_value_max)
    ids_of: dict[str, list[str]] = {}
    for pep, ids in zip(searched["peptide"], searched["protein_ids"], strict=True):
        cur = ids_of.setdefault(str(pep), [])
        cur.extend(i for i in _as_ids(ids) if i not in cur)
    rows: list[UnmappedRow] = []
    for pep, ids in ids_of.items():
        reason = peptide_reason(pep in state, state.get(pep, False), len(pep), tolerant_min_len)
        if reason is None:
            continue
        rows.append(
            _label(pep, pep, "peptide", reason, pep_matches.get(pep), _prefixed_id(ids, prefixes))
        )
    return rows


def _protein_rows(
    protein_state: pd.DataFrame, prot_matches: Mapping[str, str], prefixes: tuple[str, ...]
) -> list[UnmappedRow]:
    rows: list[UnmappedRow] = []
    seen: set[str] = set()
    cols = ("protein_id", "sequence", "assigned", "has_alignment")
    for pid, seq, assigned, aligned in zip(*(protein_state[c] for c in cols), strict=True):
        pid = str(pid)
        if bool(assigned) or pid in seen:
            continue
        seen.add(pid)
        reason: Reason = "below_identity" if bool(aligned) else "no_hit"
        rows.append(
            _label(
                pid,
                str(seq),
                "protein",
                reason,
                prot_matches.get(pid),
                _prefixed_id([pid], prefixes),
            )
        )
    return rows


def unmapped_table(
    peptides: pd.DataFrame,
    placements: pd.DataFrame,
    protein_state: pd.DataFrame,
    pep_matches: Mapping[str, str],
    prot_matches: Mapping[str, str],
    contaminant_prefixes: Sequence[str],
    *,
    tolerant_min_len: int = 10,
    q_value_max: float = _DEFAULT_Q_VALUE_MAX,
) -> list[UnmappedRow]:
    """Rows of ``mapping/unmapped.tsv`` for every peptide and protein that landed nowhere.

    ``pep_matches``/``prot_matches`` come from :func:`match_contaminant_peptides` and
    :func:`match_contaminant_proteins` (empty when the library is absent); ``tolerant_min_len``
    is ``mapping.tolerant.min_len``. Only searched peptides enter the bucket: engine decoys
    (``is_decoy``) and peptides with ``q_value > q_value_max`` (``input.q_value_max``) are
    excluded; ``q_value == q_value_max`` and a missing q-value are kept.
    """
    prefixes = tuple(p for p in contaminant_prefixes if p)
    rows = _peptide_rows(peptides, placements, pep_matches, prefixes, tolerant_min_len, q_value_max)
    rows += _protein_rows(protein_state, prot_matches, prefixes)
    return sorted(rows, key=lambda r: (r.kind, r.sequence_id))
