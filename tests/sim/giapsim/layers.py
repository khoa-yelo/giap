"""Precomputed annotation layers and their truth (issue M2-05; GIAP_SPEC.md §11.3, M2-02 formats).

Each layer lists the true elements planted by the simulator plus false positives: for a layer
with ``n`` true rows, ``Poisson(rate * n)`` false rows drawn with the simulator's RNG. The truth
table of a layer has the layer's columns plus ``true_element`` and exactly the layer's rows.
Coordinates are 1-based closed (M2-02); operon-pair probabilities are a logistic function of
the intergenic distance (clipped to +-1 kb) plus noise (UniOP-like: high within TUs, low between).

Values follow what the parse steps write (M2-12, M2-15, M2-18, M2-19): terminators from
``TerminatorNet``; regions ``<type>_<n>`` numbered per type in contig order, then start and end;
prophages ``geNOMAD``, BGCs ``antiSMASH`` (score ``NA``), CRISPR arrays ``Bakta`` (score ``NA``,
no ``cas`` genes), plasmids ``Bakta+geNOMAD`` (the whole replicon). Every region with an interval,
false positives included, lists the genes with >= 50 % of their CDS inside (M2-19 membership).
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np

from giap.annotate.precomputed import PRECOMPUTED_COLUMNS
from giap.io.tsv import write_tsv

LAYERS = ("operon_pairs", "promoters", "tx_rates", "terminators", "rbs", "regions")
TRUTH_LAYERS = ("promoters", "terminators", "rbs", "regions")
_TERM_LEN = 26


@dataclass
class Truth:
    """True elements, 1-based closed coordinates."""

    tus: list[list[Any]] = field(default_factory=list)
    promoters: list[tuple[str, str, int]] = field(default_factory=list)  # contig, strand, tss
    terminators: list[tuple[str, int, int, str]] = field(default_factory=list)
    rbs: list[tuple[str, int, str]] = field(default_factory=list)  # locus, spacer, start codon
    regions: list[tuple[str, str, int, int]] = field(default_factory=list)  # kind, contig, a, b


def _site(rng: np.random.Generator, lengths: Mapping[str, int], width: int) -> tuple[str, int]:
    ids = list(lengths)
    w = np.array([lengths[c] for c in ids], dtype=float)
    cid = ids[int(rng.choice(len(ids), p=w / w.sum()))]
    return cid, int(rng.integers(1, lengths[cid] - width + 2))


def _strand(rng: np.random.Generator) -> str:
    return "+" if rng.random() < 0.5 else "-"


def _r(x: float) -> float:
    return round(float(x), 4)


def build(
    rng: np.random.Generator,
    rate: float,
    truth: Truth,
    lengths: Mapping[str, int],
    genes: Sequence[Any],
) -> dict[str, list[list[Any]]]:
    """Layer rows (each ending with ``true_element``); ``genes`` are the coding Feats."""
    order = {c: k for k, c in enumerate(lengths)}
    out: dict[str, list[list[Any]]] = {}

    prom: list[list[Any]] = [
        [c, s, t, _r(rng.uniform(0.6, 1.0)), True] for c, s, t in truth.promoters
    ]
    for _ in range(int(rng.poisson(rate * len(prom)))):
        c, p = _site(rng, lengths, 1)
        prom.append([c, _strand(rng), p, _r(rng.uniform(0.5, 0.9)), False])
    prom.sort(key=lambda r: (order[r[0]], r[2], r[1]))
    out["promoters"] = prom
    out["tx_rates"] = [[c, s, p, _r(rng.lognormal(0.0, 1.0))] for c, s, p, _, _ in prom]

    term: list[list[Any]] = [[c, a, b, s, _r(rng.uniform(0.6, 1.0)), "TerminatorNet", True]
            for c, a, b, s in truth.terminators]  # fmt: skip
    for _ in range(int(rng.poisson(rate * len(term)))):
        c, p = _site(rng, lengths, _TERM_LEN)
        term.append([c, p, p + _TERM_LEN - 1, _strand(rng), _r(rng.uniform(0.5, 0.9)),
                     "TerminatorNet", False])  # fmt: skip
    term.sort(key=lambda r: (order[r[0]], r[1], r[3]))
    out["terminators"] = term

    def rbs_row(locus: str, sp: int, codon: str, true: bool) -> list[Any]:
        sd = rng.uniform(-11.0, -8.0) if true else rng.uniform(-6.0, -2.0)
        tir = rng.lognormal(8.0 if true else 6.0, 0.8)
        dg = [_r(rng.normal(m, 1.0)) for m in (-5.0, 0.5, -1.0)]
        return [locus, _r(tir), _r(sd), *dg, sp, codon, true]

    rbs = [rbs_row(locus, sp, codon, True) for locus, sp, codon in truth.rbs]
    planted = {r[0] for r in rbs}
    free = [g for g in genes if g.locus not in planted]
    pool = free or list(genes)
    n_fp = int(rng.poisson(rate * len(rbs)))
    idx: Any = rng.choice(len(pool), size=n_fp, replace=n_fp > len(pool)) if pool and n_fp else []
    for i in idx:
        g = pool[int(i)]
        rbs.append(rbs_row(g.locus, int(rng.integers(3, 15)), g.extra["start_type"], False))
    rank = {g.locus: k for k, g in enumerate(genes)}
    rbs.sort(key=lambda r: (rank[r[0]], not r[-1]))
    out["rbs"] = rbs

    regions: list[list[Any]] = [
        [kind, c, a, b, _TOOLS[kind], _score(rng, kind, 0.7, 1.0), None, None, True]
        for kind, c, a, b in truth.regions
    ]
    for _ in range(int(rng.poisson(rate * len(regions)))):
        kind = ("prophage", "bgc")[int(rng.integers(2))]
        width = int(rng.integers(5_000, 15_001))
        c, p = _site(rng, lengths, min(width, min(lengths.values())))
        b = min(lengths[c], p + width - 1)
        regions.append([kind, c, p, b, _TOOLS[kind], _score(rng, kind, 0.3, 0.7), None, None,
                        False])  # fmt: skip
    regions.sort(key=lambda r: (order[r[1]], r[2], r[3], r[0]))
    numbered: Counter[str] = Counter()
    out["regions"] = []
    for r in regions:
        numbered[r[0]] += 1
        if r[0] != "crispr":
            r[6] = members(genes, r[1], r[2] - 1, r[3], lengths[r[1]]) or None
        out["regions"].append([f"{r[0]}_{numbered[r[0]]}", *r])

    out["operon_pairs"] = _operon_pairs(rng, genes, order)
    return out


_TOOLS = {
    "prophage": "geNOMAD",
    "plasmid": "Bakta+geNOMAD",
    "bgc": "antiSMASH",
    "crispr": "Bakta",
}
_NO_SCORE = frozenset({"bgc", "crispr"})
_MAX_DIST = 1_000.0  # operon-pair distances are clipped before the exponential


def _score(rng: np.random.Generator, kind: str, lo: float, hi: float) -> float | None:
    """Region score in ``[lo, hi)``; ``None`` (NA) for types whose tool reports none. The draw
    is made either way, so the RNG stream does not depend on the type."""
    x = _r(rng.uniform(lo, hi))
    return None if kind in _NO_SCORE else x


def members(genes: Sequence[Any], contig: str, r0: int, r1: int, n: int) -> list[str]:
    """Loci of ``genes`` with >= 50 % of their CDS in ``[r0, r1)`` (0-based) on ``contig`` of
    length ``n``; gene ends may exceed ``n`` (origin-crossing genes)."""
    out = []
    for g in genes:
        if g.contig != contig:
            continue
        ov = sum(max(0, min(g.end, r1 + k) - max(g.start, r0 + k)) for k in (-n, 0, n))
        if 2 * ov >= g.end - g.start:
            out.append(g.locus)
    return out


def _operon_pairs(
    rng: np.random.Generator, genes: Sequence[Any], order: Mapping[str, int]
) -> list[list[Any]]:
    ordered = sorted(genes, key=lambda g: (order[g.contig], g.start, g.end))
    rows = []
    for a, b in pairwise(ordered):
        if a.contig != b.contig or a.strand != b.strand:
            continue
        d = b.start - a.end
        x = float(np.clip(d, -_MAX_DIST, _MAX_DIST))
        p = 1.0 / (1.0 + np.exp((x - 60.0) / 15.0)) + rng.normal(0.0, 0.08)
        rows.append([a.locus, b.locus, _r(min(0.999, max(0.001, p))), "uniop"])
    return rows


def write(
    out: Path, rows: Mapping[str, list[list[Any]]]
) -> tuple[dict[str, Path], dict[str, Path]]:
    """Write ``precomputed/<layer>.tsv`` and ``truth/<layer>.tsv``; returns both path maps."""
    (out / "precomputed").mkdir(parents=True, exist_ok=True)
    pre, tru = {}, {}
    for key in LAYERS:
        cols = list(PRECOMPUTED_COLUMNS[key])
        data = rows[key]
        pre[key] = out / "precomputed" / f"{key}.tsv"
        if key in TRUTH_LAYERS:
            write_tsv(pre[key], cols, [r[:-1] for r in data])
            tru[key] = out / "truth" / f"{key}.tsv"
            write_tsv(tru[key], [*cols, "true_element"], data)
        else:
            write_tsv(pre[key], cols, data)
    return pre, tru


TUS_COLUMNS = [
    "tu_id",
    "contig",
    "start",
    "end",
    "strand",
    "genes",
    "promoter",
    "tss",
    "terminator",
]


def write_config(path: Path, bakta: Path, layers: Mapping[str, Path]) -> None:
    """``config.yaml`` with every ``annotation.precomputed`` key set (absolute paths)."""
    import json

    lines = ["annotation:", "  precomputed:", f"    bakta_json: {json.dumps(str(bakta.resolve()))}"]
    lines += [f"    {k}: {json.dumps(str(p.resolve()))}" for k, p in layers.items()]
    path.write_text("\n".join(lines) + "\n")


_START_CODONS = ("ATG", "GTG", "TTG")
_COMP = str.maketrans("ACGT", "TGCA")


def rbs_novel_rows(
    rng: np.random.Generator, text: Mapping[str, str], truth: set[tuple[str, str, int]]
) -> list[list[Any]]:
    """The ``rbs_novel`` layer (M7-08 rule 6, ADR-0024): one row per start codon on both strands
    of every replicon (origin-crossing codons left out), ``start`` = the 1-based first base of
    the codon in transcription direction; ``tir_ostir`` lognormal(8, 0.8) at the 0-based
    ``truth`` keys ``(contig, strand, first base)``, lognormal(6, 0.8) elsewhere."""
    keys: list[tuple[str, str, int, str]] = []
    for c, seq in text.items():
        mine: list[tuple[str, str, int, str]] = []
        s = seq.upper()
        for i in range(len(s) - 2):
            cod = s[i : i + 3]
            if cod in _START_CODONS:
                mine.append((c, "+", i, cod))
            rc = cod[::-1].translate(_COMP)
            if rc in _START_CODONS:
                mine.append((c, "-", i + 2, rc))
        keys += sorted(mine, key=lambda k: (k[2], k[1]))
    n = len(keys)
    true = np.array([(c, st, p) in truth for c, st, p, _ in keys], dtype=bool)
    tir = np.where(true, rng.lognormal(8.0, 0.8, n), rng.lognormal(6.0, 0.8, n))
    sd = np.where(true, rng.uniform(-11.0, -8.0, n), rng.uniform(-6.0, -2.0, n))
    dg = [rng.normal(m, 1.0, n) for m in (-5.0, 0.5, -1.0)]
    spacer = rng.integers(3, 15, n)
    return [
        [c, st, p + 1, _r(tir[k]), _r(sd[k]), _r(dg[0][k]), _r(dg[1][k]), _r(dg[2][k]),
         int(spacer[k]), cod]
        for k, (c, st, p, cod) in enumerate(keys)
    ]  # fmt: skip


def write_rbs_novel(out: Path, rows: Sequence[Sequence[Any]]) -> Path:
    """``precomputed/rbs_novel.tsv`` (columns of ``PRECOMPUTED_COLUMNS['rbs_novel']``)."""
    path = out / "precomputed" / "rbs_novel.tsv"
    write_tsv(path, list(PRECOMPUTED_COLUMNS["rbs_novel"]), rows)
    return path
