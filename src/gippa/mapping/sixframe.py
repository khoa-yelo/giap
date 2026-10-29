"""Six-frame translation and frame <-> genome coordinate maps (GIPPA_SPEC.md §7.2.1, §6.5).

Every contig is translated in all six frames with the run's genetic code, 4, 11 or 25 (ADR-0067;
``*`` at the table's stops; a codon with any non-ACGT base -> ``X``). I is folded into L in
targets and queries. Circular replicons get their first ``wrap_nt`` nucleotides appended so that
peptides can span the origin; hits are mapped back modulo the replicon length (origin-crossing
intervals stay unwrapped, §6.5) and deduplicated with :func:`dedupe_loci`.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
from numba import njit

from gippa.model.coords import Interval, Replicon, wrap_position
from gippa.model.genetic_code import codon_table

if TYPE_CHECKING:
    from gippa.mapping.decoys import DecoyType
    from gippa.mapping.initiator import StartCells
    from gippa.mapping.readthrough import StopCells

_log = logging.getLogger(__name__)
NO_START_CELLS_EVENT = "mapping.index_without_start_cells"

ALPHABET: str = "ACDEFGHKLMNPQRSTVWY*X"  # codes 0..18 residues (L includes I), 19 '*', 20 'X'
STOP: int = 19
AMBIG: int = 20
NOMATCH: int = 21  # query letters outside the residue alphabet (B J O U X Z *) never match
READTHROUGH: int = 22  # search-text code of a consumable (annotated) stop; never in frames
MAX_PEPTIDE_AA: int = 40  # longest mapped peptide (mapping.exact / tolerant max_diff classes)
MAX_DIFF: int = 2  # largest tolerant difference budget k (the +-2 band of M3-07)
# Circular overlap: a window of max_len + k residues crossing the origin (M3-07 review decision 2).
WRAP_NT: int = 3 * (MAX_PEPTIDE_AA + MAX_DIFF)

_BASES = "ACGT"  # nucleotide codes 0..3; 4 = any other byte
_NT_OTHER = 4


def _letter_codes(*, query: bool) -> np.ndarray:
    """256-entry byte -> code table; I -> L; unknown letters -> NOMATCH (query) or AMBIG."""
    table = np.full(256, NOMATCH if query else AMBIG, dtype=np.uint8)
    letters = ALPHABET[:STOP] if query else ALPHABET
    for code, letter in enumerate(letters):
        table[ord(letter)] = table[ord(letter.lower())] = code
    for src in ("I", "i"):
        table[ord(src)] = ALPHABET.index("L")
    return table


def _nt_codes() -> np.ndarray:
    table = np.full(256, _NT_OTHER, dtype=np.uint8)
    for code, base in enumerate(_BASES):
        table[ord(base)] = table[ord(base.lower())] = code
    return table


@cache
def _codon_luts(table: int) -> tuple[np.ndarray, np.ndarray]:
    """64-entry tables indexed by ``16*b0 + 4*b1 + b2`` (ACGT codes) -> amino-acid code, and ->
    ASCII residue letter (no I -> L), for ``table`` (raises for an unsupported table)."""
    codes = np.empty(64, dtype=np.uint8)
    letters = np.empty(64, dtype=np.uint8)
    aa_code = _letter_codes(query=False)
    for codon, aa in codon_table(table).items():
        idx = 16 * _BASES.index(codon[0]) + 4 * _BASES.index(codon[1]) + _BASES.index(codon[2])
        codes[idx], letters[idx] = aa_code[ord(aa)], ord(aa)
    codes.flags.writeable = letters.flags.writeable = False
    return codes, letters


_QUERY_LUT = _letter_codes(query=True)
_TARGET_LUT = _letter_codes(query=False)
_NT_LUT = _nt_codes()
_DECODE_LUT = np.full(256, ord("#"), dtype=np.uint8)
_DECODE_LUT[: len(ALPHABET)] = np.frombuffer(ALPHABET.encode("ascii"), dtype=np.uint8)


@njit(cache=True, nogil=True)
def _translate_kernel(
    nt: np.ndarray, nt_lut: np.ndarray, codon_lut: np.ndarray, offset: int, minus: bool
) -> np.ndarray:  # pragma: no cover - compiled by numba; exercised through translate
    """Translate raw ASCII ``nt`` from ``offset``; ``minus`` reads the reverse complement."""
    n = nt.shape[0]
    m = (n - offset) // 3 if n > offset else 0
    out = np.empty(m, dtype=np.uint8)
    for i in range(m):
        if minus:
            p = n - 1 - offset - 3 * i
            c0 = nt_lut[nt[p]]
            c1 = nt_lut[nt[p - 1]]
            c2 = nt_lut[nt[p - 2]]
            if c0 > 3 or c1 > 3 or c2 > 3:
                out[i] = 20
            else:
                out[i] = codon_lut[16 * (3 - c0) + 4 * (3 - c1) + (3 - c2)]
        else:
            p = offset + 3 * i
            c0 = nt_lut[nt[p]]
            c1 = nt_lut[nt[p + 1]]
            c2 = nt_lut[nt[p + 2]]
            if c0 > 3 or c1 > 3 or c2 > 3:
                out[i] = 20
            else:
                out[i] = codon_lut[16 * c0 + 4 * c1 + c2]
    return out


def _translate(nt: np.ndarray, offset: int, minus: bool, table: int) -> np.ndarray:
    out: np.ndarray = _translate_kernel(nt, _NT_LUT, _codon_luts(table)[0], offset, minus)
    return out


def _encode(text: str, lut: np.ndarray) -> np.ndarray:
    raw = np.frombuffer(text.encode("ascii", errors="replace"), dtype=np.uint8)
    out: np.ndarray = lut[raw]
    return out


def encode_query(peptide: str) -> np.ndarray:
    """Encode a query peptide (I -> L); letters outside the 19 residues -> ``NOMATCH``."""
    return _encode(peptide, _QUERY_LUT)


def encode_target(aa: str) -> np.ndarray:
    """Encode a target amino-acid string over ``ALPHABET`` (I -> L; unknown -> ``AMBIG``)."""
    return _encode(aa, _TARGET_LUT)


def decode(codes: np.ndarray) -> str:
    """Decode codes back to letters; ``AMBIG`` -> ``X``, ``NOMATCH`` and every code >= 21 -> ``#``.

    ``#`` lies outside the protein alphabet, so a decoded query never seems to match a target X.
    """
    out: np.ndarray = _DECODE_LUT[np.asarray(codes, dtype=np.uint8)]
    return out.tobytes().decode("ascii")


def translate(nt: bytes, table: int = 11) -> np.ndarray:
    """Translate ``nt`` in frame 0 with ``table`` to codes; a trailing partial codon is dropped."""
    return _translate(np.frombuffer(nt, dtype=np.uint8), 0, False, table)


def translate_text(nt: str, table: int = 11, offset: int = 0) -> str:
    """``table`` codon by codon from ``offset`` as text: no I -> L, no start-codon rule, ``*`` at
    stops, a codon with a non-ACGT base -> ``X``; a trailing partial codon is dropped."""
    raw = np.frombuffer(nt.encode("ascii", errors="replace"), dtype=np.uint8)
    codes: np.ndarray = _translate_kernel(raw, _NT_LUT, _codon_luts(table)[1], offset, False)
    letters = np.where(codes == AMBIG, ord("X"), codes).astype(np.uint8)
    return letters.tobytes().decode("ascii")


def translate11(nt: bytes) -> np.ndarray:
    """:func:`translate` with table 11 (alias kept, M11-16 decision D13)."""
    return translate(nt, 11)


def translate11_text(nt: str, offset: int = 0) -> str:
    """:func:`translate_text` with table 11 (alias kept, M11-16 decision D13)."""
    return translate_text(nt, 11, offset)


@dataclass(frozen=True, slots=True)
class Frame:
    """One translated frame of the extended (circularly wrapped) sequence of ``ext_nt`` nucleotides.

    ``offset`` (0..2) counts from the 5' end of the strand of the *extended* sequence, so on the
    minus strand of a circular replicon it is not the replicon phase. Callers that need the TSV
    ``frame`` use :meth:`SixFrameIndex.frame_label` on the mapped interval, never ``offset``.
    """

    contig: str
    strand: Literal["+", "-"]
    offset: int
    aa: np.ndarray
    ext_nt: int


@dataclass(frozen=True)
class SixFrameIndex:
    """Six frames per contig in FASTA order (``+0, +1, +2, -0, -1, -2``), translated with the
    genetic code ``translation_table``."""

    replicons: Mapping[str, Replicon]
    frames: tuple[Frame, ...]
    translation_table: int = 11

    def to_interval(self, frame_id: int, aa_start: int, aa_end: int) -> Interval:
        """Map amino acids ``[aa_start, aa_end)`` of a frame to a genomic interval.

        Raises ValueError when the window lies outside the frame or spans more than the replicon.
        """
        frame = self.frames[frame_id]
        if not 0 <= aa_start < aa_end <= len(frame.aa):
            raise ValueError(
                f"aa window [{aa_start}, {aa_end}) outside frame {frame_id} "
                f"of length {len(frame.aa)}"
            )
        length = self.replicons[frame.contig].length
        if 3 * (aa_end - aa_start) > length:
            raise ValueError(
                f"aa window [{aa_start}, {aa_end}) spans {3 * (aa_end - aa_start)} nt, "
                f"longer than replicon {frame.contig!r} ({length} nt)"
            )
        if frame.strand == "+":
            start = frame.offset + 3 * aa_start
        else:
            start = frame.ext_nt - frame.offset - 3 * aa_end
        start = wrap_position(start, self.replicons[frame.contig])
        return Interval(frame.contig, start, start + 3 * (aa_end - aa_start), frame.strand)

    def frame_label(self, iv: Interval) -> int:
        """TSV ``frame``: phase (1..3) from the strand's 5' end of the replicon, signed."""
        if iv.strand == "+":
            return 1 + iv.start % 3
        length = self.replicons[iv.contig].length
        return -(1 + ((length - iv.end) % length) % 3)


def build_sixframe(
    genome: Mapping[str, str],
    replicons: Mapping[str, Replicon],
    wrap_nt: int = WRAP_NT,
    table: int = 11,
) -> SixFrameIndex:
    """Translate every contig of ``genome`` (FASTA order) in six frames with ``table``."""
    _codon_luts(table)  # an unsupported table raises, even for an empty genome
    frames: list[Frame] = []
    for contig, seq in genome.items():
        rep = replicons[contig]
        if rep.length != len(seq):
            raise ValueError(
                f"contig {contig!r}: sequence length {len(seq)} != replicon length {rep.length}"
            )
        raw = seq.encode("ascii", errors="replace")
        if rep.topology == "circular":
            raw += raw[: min(wrap_nt, rep.length)]
        nt = np.frombuffer(raw, dtype=np.uint8)
        for strand in ("+", "-"):
            for offset in range(3):
                aa = _translate(nt, offset, strand == "-", table)
                frames.append(Frame(contig, strand, offset, aa, len(raw)))
    return SixFrameIndex(replicons=dict(replicons), frames=tuple(frames), translation_table=table)


def dedupe_loci(ivs: Iterable[Interval], replicons: Mapping[str, Replicon]) -> list[Interval]:
    """Reduce starts modulo the replicon length and drop duplicates (first occurrence kept)."""
    seen: dict[Interval, None] = {}
    for iv in ivs:
        start = wrap_position(iv.start, replicons[iv.contig])
        seen.setdefault(Interval(iv.contig, start, start + iv.length_bp, iv.strand), None)
    return list(seen)


# 2: + consumable stop cells (M7-16, ADR-0058); 3: + translation_table (ADR-0067);
# 4: + initiator start sites and cells (M11-20, ADR-0068)
_INDEX_FORMAT = 4


def _frame_meta(frames: tuple[Frame, ...]) -> list[list[object]]:
    return [[f.contig, f.strand, f.offset, f.ext_nt] for f in frames]


def _pack(frames: tuple[Frame, ...]) -> tuple[np.ndarray, np.ndarray]:
    bounds = np.zeros(len(frames) + 1, dtype=np.int64)
    np.cumsum([len(f.aa) for f in frames], out=bounds[1:])
    aa = np.concatenate([f.aa for f in frames]) if frames else np.zeros(0, dtype=np.uint8)
    return aa.astype(np.uint8, copy=False), bounds


def _unpack(meta: list[list[object]], aa: np.ndarray, bounds: np.ndarray) -> tuple[Frame, ...]:
    out: list[Frame] = []
    for i, (contig, strand, offset, ext_nt) in enumerate(meta):
        lo, hi = int(bounds[i]), int(bounds[i + 1])
        out.append(Frame(str(contig), strand, int(offset), aa[lo:hi].copy(), int(ext_nt)))  # type: ignore[arg-type,call-overload]
    return tuple(out)


def _pack_cells(cells: tuple[np.ndarray, ...], n: int) -> tuple[np.ndarray, np.ndarray]:
    if len(cells) != n:
        raise ValueError(f"{len(cells)} cell arrays for {n} frames")
    bounds = np.zeros(n + 1, dtype=np.int64)
    np.cumsum([len(c) for c in cells], out=bounds[1:])
    flat = np.concatenate([np.asarray(c, dtype=np.int64) for c in cells]) if n else bounds[:0]
    return flat.astype(np.int64, copy=False), bounds


def save_index(
    path: Path,
    index: SixFrameIndex,
    decoys: Mapping[DecoyType, tuple[Frame, ...]],
    cells: StopCells | None = None,
    starts: StartCells | None = None,
) -> None:
    """Write the ``sixframe_index`` work file (``.npz``, format 4): its translation table, target
    frames, decoy frames, replicons, the consumable stop cells and the initiator start cells of
    the target and each decoy type (none without ``cells`` / ``starts``)."""
    frame_meta = _frame_meta(index.frames)
    n = len(index.frames)
    none = tuple(np.zeros(0, dtype=np.int64) for _ in range(n))
    meta = {
        "format": _INDEX_FORMAT,
        "translation_table": index.translation_table,
        "replicons": [asdict(r) for r in index.replicons.values()],
        "frames": frame_meta,
        "decoys": list(decoys),
        "stop_sites": [] if cells is None else [asdict(s) for s in cells.sites],
        "start_sites": [] if starts is None else [asdict(s) for s in starts.sites],
    }
    arrays: dict[str, np.ndarray] = {"meta": np.array(json.dumps(meta))}
    arrays["target_aa"], arrays["target_bounds"] = _pack(index.frames)
    arrays["target_rt"], arrays["target_rt_bounds"] = _pack_cells(
        none if cells is None else cells.target, n
    )
    arrays["target_init"], arrays["target_init_bounds"] = _pack_cells(
        none if starts is None else starts.target, n
    )
    for kind, frames in decoys.items():
        if _frame_meta(frames) != frame_meta:
            raise ValueError(f"{kind} decoy frames do not mirror the index frames")
        arrays[f"{kind}_aa"], arrays[f"{kind}_bounds"] = _pack(frames)
        if cells is not None and kind not in cells.decoys:
            raise ValueError(f"stop cells have no {kind} decoy cells")
        dc = none if cells is None else cells.decoys[kind]
        arrays[f"{kind}_rt"], arrays[f"{kind}_rt_bounds"] = _pack_cells(dc, n)
        if starts is not None and kind not in starts.decoys:
            raise ValueError(f"start cells have no {kind} decoy cells")
        di = none if starts is None else starts.decoys[kind]
        arrays[f"{kind}_init"], arrays[f"{kind}_init_bounds"] = _pack_cells(di, n)
    with Path(path).open("wb") as fh:
        np.savez_compressed(fh, **arrays)  # type: ignore[arg-type]


def _check_format(path: Path, meta: Mapping[str, object], table: int) -> int:
    """The index's translation table; raises unless it is ``table`` (format 2 is table 11, D3).
    Formats 2 and 3 load with no initiator start cells (:func:`load_start_cells`)."""
    fmt = meta.get("format")
    if fmt == 1:
        raise ValueError(
            f"{path}: six-frame index format 1 has no stop cells; rebuild the six-frame index"
        )
    built = 11 if fmt == 2 else meta.get("translation_table")
    if fmt not in (2, 3, _INDEX_FORMAT) or not isinstance(built, int) or isinstance(built, bool):
        raise ValueError(f"{path}: unsupported sixframe_index format {fmt!r}")
    if built != table:
        raise ValueError(
            f"{path}: six-frame index built for translation table {built}; this run uses "
            f"{table}; rebuild the six-frame index"
        )
    return built


def load_index(
    path: Path, *, table: int
) -> tuple[SixFrameIndex, dict[DecoyType, tuple[Frame, ...]]]:
    """Read a work file written by :func:`save_index` for the run's ``table``.

    Raises:
        ValueError: on an index built for another table, a format-2 file unless ``table`` is 11,
            a format-1 file ("rebuild the six-frame index") or an unknown format.
    """
    with np.load(Path(path), allow_pickle=False) as npz:
        meta = json.loads(str(npz["meta"]))
        built = _check_format(path, meta, table)
        replicons = {r["id"]: Replicon(**r) for r in meta["replicons"]}
        frames = _unpack(meta["frames"], npz["target_aa"], npz["target_bounds"])
        decoys: dict[DecoyType, tuple[Frame, ...]] = {
            kind: _unpack(meta["frames"], npz[f"{kind}_aa"], npz[f"{kind}_bounds"])
            for kind in meta["decoys"]
        }
    return SixFrameIndex(replicons=replicons, frames=frames, translation_table=built), decoys


def load_stop_cells(path: Path, *, table: int) -> StopCells:
    """The stop sites and cells of a work file built for the run's ``table``.

    Raises:
        ValueError: as :func:`load_index`.
    """
    from gippa.mapping.readthrough import StopCells, StopSite

    def split(flat: np.ndarray, bounds: np.ndarray) -> tuple[np.ndarray, ...]:
        return tuple(flat[bounds[i] : bounds[i + 1]].copy() for i in range(len(bounds) - 1))

    with np.load(Path(path), allow_pickle=False) as npz:
        meta = json.loads(str(npz["meta"]))
        _check_format(path, meta, table)
        sites = tuple(
            StopSite(**{**s, "locus_tags": tuple(s["locus_tags"])}) for s in meta["stop_sites"]
        )
        target = split(npz["target_rt"], npz["target_rt_bounds"])
        decoys = {k: split(npz[f"{k}_rt"], npz[f"{k}_rt_bounds"]) for k in meta["decoys"]}
    return StopCells(sites=sites, target=target, decoys=decoys)


def _split(flat: np.ndarray, bounds: np.ndarray) -> tuple[np.ndarray, ...]:
    return tuple(flat[bounds[i] : bounds[i + 1]].copy() for i in range(len(bounds) - 1))


def load_start_cells(path: Path, *, table: int) -> StartCells:
    """The initiator start sites and cells of a work file built for the run's ``table``. A
    format-2 or format-3 file (built before M11-20) has none: one WARNING event
    ``mapping.index_without_start_cells``.

    Raises:
        ValueError: as :func:`load_index`.
    """
    from gippa.mapping.initiator import StartCells, StartSite

    with np.load(Path(path), allow_pickle=False) as npz:
        meta = json.loads(str(npz["meta"]))
        _check_format(path, meta, table)
        n = len(meta["frames"])
        if meta["format"] != _INDEX_FORMAT:
            none = tuple(np.zeros(0, dtype=np.int64) for _ in range(n))
            _log.warning(
                "%s: six-frame index format %s was built before M11-20; rebuild the six-frame "
                "index for initiator-M credit", path, meta["format"],
                extra={"event": NO_START_CELLS_EVENT,
                       "data": {"path": str(path), "format": meta["format"]}},
            )  # fmt: skip
            return StartCells(sites=(), target=none, decoys={k: none for k in meta["decoys"]})
        sites = tuple(
            StartSite(**{**s, "locus_tags": tuple(s["locus_tags"])}) for s in meta["start_sites"]
        )
        target = _split(npz["target_init"], npz["target_init_bounds"])
        decoys = {k: _split(npz[f"{k}_init"], npz[f"{k}_init_bounds"]) for k in meta["decoys"]}
    return StartCells(sites=sites, target=target, decoys=decoys)
