"""Recall of planted peptides with 1-2 differences, tolerant route (M3-07; GIAP_SPEC.md §7.2.4).

Claim (§7.2.4 "Seed and verify", "Completeness requirement"): the q-gram threshold
``t = (L - q + 1) - k (q + 1)`` is lossless, so a peptide within ``k`` differences of a genome
window is never missed. giapsim genomes (circular, 50 kb and 100 kb; seeds
``seeds.calibration_base`` + i); every annotated protein is digested with trypsin (no missed
cleavages) and peptides of 10-30 aa get 1..k planted differences (k = 1 at 10-11 aa, 2 at >= 12),
each drawn from substitution, adjacent swap, mass-equivalent block and near-equal block (both
directions; substitution when the window has no block partner). ``I`` is never used as a
replacement (I = L is free).

Criteria, per planted peptide with ``n`` planted operations and origin interval ``o``
(lead decision 6 after code review, PR #48). The *origin cost* ``c_o`` is the oracle cost
(``oracle_tolerant.align_window``, unbanded, independent of the implementation) of the best
window of the origin's frame whose ends are within ``k`` residues of the planted ends; the
*brute-force minimum* ``C`` is the unbanded DP minimum over every window of every target frame
(``brute_tolerant.semi_global_costs``):

- planting sanity (exact): ``c_o <= n``;
- lossless (exact, no tolerance): a target row exists and the reported minimal cost equals ``C``;
  whenever ``c_o == C`` a reported locus overlaps ``o`` in the same codon phase (free window ends
  may shift the minimal window by a residue or two). A miss of the origin is allowed only when
  the oracle shows a strictly cheaper window elsewhere (``C < c_o``), not because of the planted
  operation count;
- recall: the share of planted peptides whose origin is reported is ``>= mapping_recall_min``;
- every difference type was planted (the criterion is not vacuous for any type).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "property"))
from brute_tolerant import INF, semi_global_costs
from oracle_tolerant import align_window

from giap.config.model import TolerantConfig
from giap.io.fasta import read_genome_fasta
from giap.mapping.decoys import build_decoys
from giap.mapping.sixframe import build_sixframe, decode, encode_query, translate11
from giap.mapping.tolerant import match_tolerant
from giap.model import Interval, Replicon
from giap.model.genes import gene_model_from_bakta_json
from giap.rng import for_step
from giapsim.simulate import SimConfig, simulate

pytestmark = pytest.mark.calibration

RESIDUES = "ACDEFGHKLMNPQRSTVWY"  # no I: I -> L would be free
_BLOCK_RAW = (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
_NEAR_RAW = (("W", "AD"), ("W", "EG"), ("R", "GV"))


def _pairs(raw: tuple[tuple[str, str], ...]) -> list[tuple[str, str]]:
    """(genome segment, query segment), every residue order, both directions."""
    out = set()
    for a, b in raw:
        for x in {a, a[::-1]}:
            for y in {b, b[::-1]}:
                out |= {(x, y), (y, x)}
    return sorted(out)


SEGS = {"block": _pairs(_BLOCK_RAW), "near_block": _pairs(_NEAR_RAW)}


def _tryptic(aa: str) -> list[tuple[int, str]]:
    out, a = [], 0
    for m in re.finditer(r"[KR](?!P)", aa):
        out.append((a, aa[a : m.end()]))
        a = m.end()
    if a < len(aa):
        out.append((a, aa[a:]))
    return out


def _plant(window: str, n_ops: int, rng: np.random.Generator) -> tuple[str, list[str]]:
    """``window`` with ``n_ops`` non-overlapping operations; returns (query, planted types)."""
    used = np.zeros(len(window), dtype=bool)  # genome residues touched (with a 1-residue margin)
    edits: list[tuple[int, int, str, str]] = []  # (start, end, query segment, type)
    for _ in range(n_ops):
        kind = str(rng.choice(["substitution", "swap", "block", "near_block"]))
        options: list[tuple[int, int, str]] = []
        if kind in SEGS:
            for g, qseg in SEGS[kind]:
                for m in re.finditer(f"(?={g})", window):
                    options.append((m.start(), m.start() + len(g), qseg))
        elif kind == "swap":
            for j in range(len(window) - 1):
                if window[j] != window[j + 1] and "I" not in window[j : j + 2]:
                    options.append((j, j + 2, window[j + 1] + window[j]))
        options = [
            o for o in options if not used[max(o[0] - 1, 0) : min(o[1] + 1, len(window))].any()
        ]
        if not options:  # no partner in this window: a substitution instead
            kind = "substitution"
            free = [j for j in range(len(window)) if not used[max(j - 1, 0) : j + 2].any()]
            j = int(rng.choice(free))
            new = str(rng.choice([r for r in RESIDUES if r != window[j].replace("I", "L")]))
            options = [(j, j + 1, new)]
        s, e, qseg = options[int(rng.integers(len(options)))]
        used[s:e] = True
        edits.append((s, e, qseg, kind))
    query = window
    for s, e, qseg, _ in sorted(edits, reverse=True):
        query = query[:s] + qseg + query[e:]
    return query, [kind for *_, kind in edits]


def _overlaps(o: Interval, f: tuple[str, int, int, str], n: int) -> bool:
    """Same contig and strand, same codon phase, overlapping modulo the replicon length ``n``."""
    contig, start, length, strand = f
    if (contig, strand) != (o.contig, o.strand):
        return False
    for shift in (-n, 0, n):
        s = start + shift
        if s < o.end and o.start < s + length and (s - o.start) % 3 == 0:
            return True
    return False


_RC = str.maketrans("ACGTN", "TGCAN")


def _origin_cost(genome: str, o: Interval, query: str, k: int) -> int:
    """Oracle cost of the best window whose ends are within ``k`` residues of the origin's."""
    n = len(genome)
    a, b = o.start - 3 * k, o.end + 3 * k
    nt = "".join(genome[x % n] for x in range(a, b))
    if o.strand == "-":
        nt = nt.translate(_RC)[::-1]
    text = decode(translate11(nt.encode("ascii")))
    q = decode(encode_query(query))
    m = (o.end - o.start) // 3
    best = INF
    for s in range(0, 2 * k + 1):
        for e in range(k + m - k, min(k + m + k, len(text)) + 1):
            if e - s < len(q) - k or e - s > len(q) + 2 * k:
                continue
            got = align_window(q, text, s, e, None, k)
            if got is not None:
                best = min(best, got[0])
    return best


def _brute_min(query: str, frames: tuple, k: int) -> int:  # type: ignore[type-arg]
    qc = encode_query(query)
    return min(int(semi_global_costs(qc, f.aa, k).min()) for f in frames)


@pytest.mark.parametrize(("length", "i"), [(50_000, 0), (100_000, 1)])
def test_planted_tolerant_peptide_recall(
    length: int, i: int, tmp_path: Path, tolerances: dict[str, object]
) -> None:
    seed = int(tolerances["seeds"]["calibration_base"]) + 700 + i  # type: ignore[index]
    res = simulate(SimConfig(genome_length=length), seed=seed, outdir=tmp_path)
    (rec,) = read_genome_fasta(res.genome_fasta)
    rep = Replicon(rec.id, len(rec.sequence), "circular")
    model = gene_model_from_bakta_json(res.bakta_json, {rec.id: rep})
    idx = build_sixframe({rec.id: rec.sequence}, {rec.id: rep})
    rng = for_step("calibration_tolerant_recall", seed)
    planted: dict[str, tuple[int, Interval]] = {}
    clash: set[str] = set()
    types_seen: set[str] = set()
    for g in model.genes:
        iv = g.interval
        for off, pep in _tryptic(g.aa):
            if not 10 <= len(pep) <= 30 or off == 0 or "*" in pep:
                continue
            k = 1 if len(pep) <= 11 else 2
            n_ops = int(rng.integers(1, k + 1))
            query, kinds = _plant(pep, n_ops, rng)
            if not 10 <= len(query) <= 40 or (len(query) <= 11 and n_ops > 1):
                continue  # a block changed the length class: k no longer covers n_ops
            s = iv.start + 3 * off if iv.strand == "+" else iv.end - 3 * (off + len(pep))
            s %= rep.length
            origin = Interval(rec.id, s, s + 3 * len(pep), iv.strand)
            if query in planted:
                clash.add(query)
            planted[query] = (n_ops, origin)
            types_seen |= set(kinds)
    for q in clash:
        del planted[q]
    assert len(planted) > 200
    assert types_seen == {"substitution", "swap", "block", "near_block"}
    queries = sorted(planted)
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", seed))}
    df = match_tolerant(
        queries, idx, decoys, dict.fromkeys(queries, "reverse"), set(queries), TolerantConfig()
    )
    tgt = df[df["search"] == "target"]
    found: dict[str, set[tuple[str, int, int, str]]] = {}  # (contig, start mod L, length, strand)
    cost: dict[str, int] = {}
    for r in tgt.itertuples(index=False):
        found.setdefault(r.peptide, set()).add(
            (r.contig, int(r.start) % rep.length, int(r.end) - int(r.start), r.strand)
        )
        cost[r.peptide] = int(r.n_diff)
    missed, wrong_min, bad_plant, hit = [], [], [], 0
    for q, (n_ops, o) in planted.items():
        k_q = 1 if len(q) <= 11 else 2  # k of the query's own length class
        c_o = _origin_cost(rec.sequence, o, q, k_q)
        if c_o > n_ops:
            bad_plant.append((q, n_ops, c_o, o))
        c_min = _brute_min(q, idx.frames, k_q)
        if q not in cost or cost[q] != c_min:
            wrong_min.append((q, cost.get(q), c_min))
        overlapped = any(_overlaps(o, f, rep.length) for f in found.get(q, set()))
        if overlapped:
            hit += 1
        elif c_o == c_min:
            missed.append((q, n_ops, o, c_o, cost.get(q)))
    assert not bad_plant, bad_plant[:10]  # the planted window is within n_ops of the query
    assert not wrong_min, wrong_min[:10]  # reported minimum == brute-force minimum
    assert not missed, missed[:10]  # lossless: the origin is reported whenever it is minimal
    recall = hit / len(planted)
    assert recall >= float(tolerances["mapping_recall_min"]), (recall, hit, len(planted))  # type: ignore[arg-type]
