"""Wrong-genetic-code warning (M11-14w; docs/specs/M11-14w-wrong-genetic-code-warning.md).

GIAP 1.0 supports genetic code 11 only (Q-K30). A genome left at the default table 11 that
probably uses code 4 (TGA = Trp) is detected from the genome alone: long open reading frames
(``>= MIN_CODONS`` codons, six frames) cover far less of the genome when TGA is read as a stop
than when it is read as a sense codon. The rule compares the two coverages on exact integers
(``20 * bp11 < 11 * bp4``, i.e. ratio < 0.55). Pure (numpy only, no I/O, no randomness).

Coordinates are 0-based half-open (GIAP_SPEC.md §6.5); every record is treated as linear.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

from giap.preflight.checks import Check, JsonValue, Level

__all__ = [
    "MIN_CODONS",
    "MIN_LONG_ORFS",
    "RATIO_WARN",
    "RULE_ID",
    "genetic_code_check",
    "long_orf_cover",
]

MIN_CODONS = 300
#: shown in ``value`` and the rule text; the decision is ``20 * bp11 < 11 * bp4``
RATIO_WARN = 0.55
_WARN_NUM, _WARN_DEN = 11, 20
MIN_LONG_ORFS = 20

RULE_ID = "genome.genetic_code"
RULE_TEXT = (
    f"long ORFs (>= {MIN_CODONS} codons, six frames): bp with TGA as stop / bp with TGA as "
    f"sense < {RATIO_WARN} -> warn (needs >= {MIN_LONG_ORFS} such ORFs)"
)
MSG_TABLE = "not evaluated: annotation.translation_table is not 11"
MSG_FEW = (
    f"not evaluated: fewer than {MIN_LONG_ORFS} open reading frames of at least {MIN_CODONS} codons"
)

S11: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})
S4: frozenset[str] = frozenset({"TAA", "TAG"})

_INVALID = 4  # any non-ACGT base
_BAD_CODON = 64  # a codon with a non-ACGT base
_CODE = np.full(256, _INVALID, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _CODE[ord(_b)] = _CODE[ord(_b.lower())] = _i
_COMP = np.array([3, 2, 1, 0, _INVALID], dtype=np.uint8)

_U8 = npt.NDArray[np.uint8]
_I64 = npt.NDArray[np.int64]


def _encode(seq: str) -> _U8:
    """Bases as 0-3 (A, C, G, T, either case) and 4 for anything else, one per character."""
    raw = np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)
    return _CODE[raw]


def _codon_index(codon: str) -> int | None:
    idx = 0
    for base in codon.upper():
        pos = "ACGT".find(base)
        if pos < 0:
            return None
        idx = 4 * idx + pos
    return idx if len(codon) == 3 else None


def _breaker_table(stops: frozenset[str]) -> npt.NDArray[np.bool_]:
    """Codon index (0-63, 64 for a non-ACGT codon) -> whether the codon ends a run."""
    table = np.zeros(_BAD_CODON + 1, dtype=np.bool_)
    table[_BAD_CODON] = True
    for codon in stops:
        idx = _codon_index(codon)
        if idx is not None:
            table[idx] = True
    return table


def _strand_runs(codes: _U8, breaker: npt.NDArray[np.bool_], min_codons: int) -> tuple[_I64, _I64]:
    """Start and end of every run of ``>= min_codons`` whole codons, all three frames, on
    ``codes`` itself (half-open, the breaking codon excluded; open runs at the end count)."""
    n = len(codes)
    empty = np.zeros(0, dtype=np.int64)
    if n < 3:
        return empty, empty
    c = codes.astype(np.int16)
    idx = 16 * c[:-2] + 4 * c[1:-1] + c[2:]
    bad = (codes[:-2] == _INVALID) | (codes[1:-1] == _INVALID) | (codes[2:] == _INVALID)
    idx[bad] = _BAD_CODON
    brk = breaker[idx]  # codon starting at position i breaks a run
    starts: list[_I64] = []
    ends: list[_I64] = []
    for frame in range(3):
        fb = brk[frame::3]  # whole codons of this frame
        m = len(fb)
        stops_at = np.flatnonzero(fb)
        bounds = np.concatenate(([-1], stops_at, [m])).astype(np.int64)
        first = bounds[:-1] + 1  # codon index of the run's first codon
        last = bounds[1:]  # codon index of the breaking codon (or m when open)
        keep = (last - first) >= min_codons
        starts.append(frame + 3 * first[keep])
        ends.append(frame + 3 * last[keep])
    return np.concatenate(starts), np.concatenate(ends)


def long_orf_cover(seq: str, stops: frozenset[str], min_codons: int) -> tuple[int, int]:
    """``(positions covered by at least one run, number of runs)`` for one record over both
    strands and all three frames per strand. A run is a maximal stretch of whole codons with no
    codon from ``stops`` (upper-case codons) and no non-ACGT base, of ``>= min_codons`` codons;
    a minus-strand run ``[a, b)`` on the reverse complement maps to ``[n-b, n-a)``."""
    n = len(seq)
    if n < 3:
        return 0, 0
    codes = _encode(seq)
    breaker = _breaker_table(stops)
    a_plus, b_plus = _strand_runs(codes, breaker, min_codons)
    a_rc, b_rc = _strand_runs(_COMP[codes[::-1]], breaker, min_codons)
    starts = np.concatenate((a_plus, n - b_rc))
    ends = np.concatenate((b_plus, n - a_rc))
    if len(starts) == 0:
        return 0, 0
    delta = np.bincount(starts, minlength=n + 1) - np.bincount(ends, minlength=n + 1)
    covered = int(np.count_nonzero(np.cumsum(delta[:n]) > 0))
    return covered, len(starts)


def _decide(bp11: int, bp4: int, n4: int) -> Level:
    """The rule on exact integers: ``info`` below ``MIN_LONG_ORFS`` long ORFs, else ``warn``
    when ``bp11 / bp4 < 0.55`` (as ``20 * bp11 < 11 * bp4``), else ``ok``."""
    if n4 < MIN_LONG_ORFS:
        return "info"
    return "warn" if _WARN_DEN * bp11 < _WARN_NUM * bp4 else "ok"


def _percent(part: int, whole: int) -> int:
    return (200 * part + whole) // (2 * whole) if whole else 0


def _message(bp11: int, bp4: int, genome_bp: int) -> str:
    q11, q4 = _percent(bp11, genome_bp), _percent(bp4, genome_bp)
    ratio = f"0.{(100 * bp11) // bp4:02d}"  # truncated, never rounded up to the cut
    return (
        "This genome probably uses genetic code 4 (TGA = Trp; Mycoplasma, Spiroplasma, "
        "Mesoplasma, Ureaplasma) or another code in which TGA is not a stop. Long open reading "
        f"frames cover {q4} % of the genome when TGA is read as a codon but only {q11} % when it "
        f"is a stop (ratio {ratio} < {RATIO_WARN}). GIAP 1.0 supports genetic code 11 only: "
        "genes will be cut at TGA, proteins and peptides that contain W will not map, and "
        "strain_difference and nterm_extension_evidence flags will be artefacts. This warning "
        "does not change any result."
    )


def genetic_code_check(seqs: Sequence[str], table: int) -> Check:
    """The ``genome.genetic_code`` check for the genome records ``seqs`` at translation table
    ``table``; the level is ``ok``, ``info`` or ``warn``, never ``error``."""
    if table != 11:
        return Check(RULE_ID, "genome", "info", RULE_TEXT, None, MSG_TABLE)
    bp11 = bp4 = n4 = 0
    for seq in seqs:
        bp11 += long_orf_cover(seq, S11, MIN_CODONS)[0]
        b, k = long_orf_cover(seq, S4, MIN_CODONS)
        bp4, n4 = bp4 + b, n4 + k
    value: dict[str, JsonValue] = {
        "min_codons": MIN_CODONS,
        "bp_tga_stop": bp11,
        "bp_tga_sense": bp4,
        "n_orfs_tga_sense": n4,
        "ratio": bp11 / bp4 if bp4 else None,
        "threshold": RATIO_WARN,
    }
    level = _decide(bp11, bp4, n4)
    message = {"info": MSG_FEW, "ok": ""}.get(level)
    if message is None:
        message = _message(bp11, bp4, sum(len(s) for s in seqs))
    return Check(RULE_ID, "genome", level, RULE_TEXT, value, message)
