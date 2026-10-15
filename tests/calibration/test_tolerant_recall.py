"""Recall of planted peptides with 1-2 differences, tolerant route (M3-07; GIAP_SPEC.md §7.2.4).

Claim (§7.2.4 "Seed and verify", "Completeness requirement"): the q-gram threshold
``t = (L - q + 1) - k (q + 1)`` is lossless, so a peptide within ``k`` differences of a genome
window is never missed. giapsim genomes (circular, 50 kb and 100 kb; seeds
``seeds.calibration_base`` + i); every annotated protein is digested with trypsin (no missed
cleavages) and peptides of 10-30 aa get 1..k planted differences (k = 1 at 10-11 aa, 2 at >= 12),
each drawn from substitution, adjacent swap, mass-equivalent block and near-equal block (both
directions; substitution when the window has no block partner). ``I`` is never used as a
replacement (I = L is free).

Criteria, per planted peptide with ``n`` planted operations and origin interval ``o``:

- lossless (exact, no tolerance): a target row exists with minimal cost ``<= n``, and either a
  reported locus overlaps ``o`` in the same codon phase (free window ends may shift the minimal
  window by a residue or two at the planted end) or the reported minimal cost is ``< n`` (a
  strictly cheaper locus elsewhere);
- recall: the share of planted peptides whose origin is reported is ``>= mapping_recall_min``;
- every difference type was planted (the criterion is not vacuous for any type).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from giap.config.model import TolerantConfig
from giap.io.fasta import read_genome_fasta
from giap.mapping.decoys import build_decoys
from giap.mapping.sixframe import build_sixframe
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
    missed, hit = [], 0
    for q, (n_ops, o) in planted.items():
        if any(_overlaps(o, f, rep.length) for f in found.get(q, set())):
            hit += 1
        elif q not in cost or cost[q] >= n_ops:
            missed.append((q, n_ops, o, cost.get(q)))
    assert not missed, missed[:10]  # lossless: zero false negatives
    recall = hit / len(planted)
    assert recall >= float(tolerances["mapping_recall_min"]), (recall, hit, len(planted))  # type: ignore[arg-type]
