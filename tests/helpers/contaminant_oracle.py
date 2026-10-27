"""Independent oracles for M3-15 (contaminants and the unmapped bucket).

Brute force on purpose: no automaton, no pandas tricks; used by the unit, property and TE5-proxy
tests to check ``gippa.mapping.contaminants`` without sharing its code.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

# Column order of the tabular blastp output (same 15 columns as tblastn, M3-13).
OUTFMT_COLUMNS: tuple[str, ...] = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore", "qlen", "nident", "positive",
)  # fmt: skip

_BASES = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODON = {
    a + b + c: _AA11[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def il(s: str) -> str:
    """I = L equivalence: every I becomes L."""
    return s.replace("I", "L")


def oracle_peptide_matches(
    peptides: Iterable[str], library: Sequence[tuple[str, str]]
) -> dict[str, str]:
    """Peptide -> lowest library ID containing it (I = L); unmatched peptides absent."""
    out: dict[str, str] = {}
    for pep in peptides:
        q = il(pep)
        ids = sorted(lid for lid, seq in library if q and q in il(seq))
        if ids:
            out[pep] = ids[0]
    return out


def passes(row: Mapping[str, str], min_identity: float, min_coverage: float) -> bool:
    """One blastp HSP passes when pident/100 >= min_identity and (qend-qstart+1)/qlen >= cov."""
    ident = float(row["pident"]) / 100.0
    lo, hi = sorted((int(row["qstart"]), int(row["qend"])))
    cov = (hi - lo + 1) / int(row["qlen"])
    return ident >= min_identity and cov >= min_coverage


def oracle_passing(
    rows: Iterable[Mapping[str, str]], min_identity: float = 0.90, min_coverage: float = 0.50
) -> dict[str, set[str]]:
    """Query ID -> set of library IDs with at least one passing HSP."""
    out: dict[str, set[str]] = {}
    for r in rows:
        if passes(r, min_identity, min_coverage):
            out.setdefault(r["qseqid"], set()).add(r["sseqid"])
    return out


def oracle_best(
    rows: Iterable[Mapping[str, str]], min_identity: float = 0.90, min_coverage: float = 0.50
) -> dict[str, str]:
    """Query ID -> library ID of the best passing HSP: highest bitscore, then identity, then
    coverage, then lowest subject ID (M3-15 lead decision 2)."""
    best: dict[str, tuple[float, float, float, str]] = {}
    for r in rows:
        if not passes(r, min_identity, min_coverage):
            continue
        lo, hi = sorted((int(r["qstart"]), int(r["qend"])))
        key = (
            -float(r["bitscore"]),
            -float(r["pident"]),
            -(hi - lo + 1) / int(r["qlen"]),
            r["sseqid"],
        )
        q = r["qseqid"]
        if q not in best or key < best[q]:
            best[q] = key
    return {q: k[3] for q, k in best.items()}


def blast_line(**kw: object) -> str:
    """One tab-separated blastp row in ``OUTFMT_COLUMNS`` order (sensible defaults)."""
    qlen = int(kw.get("qlen", 100))  # type: ignore[call-overload]
    qstart = int(kw.get("qstart", 1))  # type: ignore[call-overload]
    qend = int(kw.get("qend", qlen))  # type: ignore[call-overload]
    length = qend - qstart + 1
    pident = float(kw.get("pident", 100.0))  # type: ignore[arg-type]
    nident = round(length * pident / 100.0)
    vals: dict[str, object] = {
        "qseqid": "q1", "sseqid": "s1", "pident": f"{pident:.3f}", "length": length,
        "mismatch": length - nident, "gapopen": 0, "qstart": qstart, "qend": qend,
        "sstart": 1, "send": length, "evalue": "1e-50", "bitscore": "200",
        "qlen": qlen, "nident": nident, "positive": nident,
    }  # fmt: skip
    vals.update({k: v for k, v in kw.items() if k in vals and k not in ("qlen", "qstart", "qend")})
    vals["pident"] = f"{pident:.3f}"
    return "\t".join(str(vals[c]) for c in OUTFMT_COLUMNS)


def tryptic_digest(seq: str, min_len: int = 7, max_len: int = 40) -> list[str]:
    """Fully tryptic peptides (after K/R, not before P; no missed cleavages), sorted."""
    peps: list[str] = []
    start = 0
    for i, aa in enumerate(seq):
        nxt = seq[i + 1] if i + 1 < len(seq) else ""
        if aa in "KR" and nxt != "P":
            peps.append(seq[start : i + 1])
            start = i + 1
    peps.append(seq[start:])
    return sorted({p for p in peps if min_len <= len(p) <= max_len and "*" not in p})


def six_frames(nt: str, circular: bool = True) -> list[str]:
    """Six-frame translations (code 11) of ``nt``; circular contigs get a 120-nt wrap."""
    nt = nt.upper()
    if circular:
        nt = nt + nt[:120]
    rc = nt.translate(_COMP)[::-1]
    frames = []
    for s in (nt, rc):
        for off in range(3):
            frames.append("".join(_CODON.get(s[i : i + 3], "X") for i in range(off, len(s) - 2, 3)))
    return frames


def hits_genome(peptide: str, frames_il: Sequence[str]) -> bool:
    """Exact (I = L) substring of any six-frame translation (already passed through ``il``)."""
    q = il(peptide)
    return any(q in f for f in frames_il)
