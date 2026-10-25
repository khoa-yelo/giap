"""Toy genomes and gene models of docs/specs/M7-16-stop-codon-readthrough.md (examples 1-8, edge
cases). Test tooling only: no ``giap.mapping.readthrough`` import, so the helpers build inputs
before the implementation exists.

- :func:`toy_c1`: contig ``c1`` (linear, 126 nt, frame ``+0`` reads
  ``MSDKPLTEAGRWQLNEFK*VPAGHTLDMRK*AAGLKEPSTR*``; TGA at p = 18, TAA at 30, TAG at 41) with
  ``G1`` = CDS ``+`` ``[0, 57)``.
- :func:`toy_c2`: contig ``c2`` (linear, 300 nt) with ``G2`` = CDS ``-`` ``[120, 210)``, stop
  ``[120, 123)`` (gene-strand TAA).
- :func:`gene`, :func:`model`: Gene / GeneModel builders; :func:`write_annotation`: the
  ``parse_bakta`` work tables (genes, features, replicons) of a gene model.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from giap.model import Replicon
from giap.model.coords import Interval
from giap.model.genes import Gene, GeneModel

#: one fixed codon per residue (table 11); ``*`` is chosen per stop
CODON: Mapping[str, str] = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG",
    "R": "CGT", "S": "TCT", "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip
_COMP = str.maketrans("ACGTN", "TGCAN")

C1_AA = "MSDKPLTEAGRWQLNEFK*VPAGHTLDMRK*AAGLKEPSTR*"
C1_STOPS = {18: "TGA", 30: "TAA", 41: "TAG"}


def rc(nt: str) -> str:
    return nt.translate(_COMP)[::-1]


def back_translate(aa: str, stops: Mapping[int, str]) -> str:
    return "".join(stops[i] if a == "*" else CODON[a] for i, a in enumerate(aa))


def gene(
    tag: str, contig: str, start: int, end: int, strand: str, feature_type: str = "cds",
    partial: bool = False,
) -> Gene:  # fmt: skip
    return Gene(
        tag, feature_type, Interval(contig, start, end, strand), "", None, None,  # type: ignore[arg-type]
        partial=partial,
    )  # fmt: skip


def model(replicons: Mapping[str, Replicon], genes: Sequence[Gene]) -> GeneModel:
    order = {c: i for i, c in enumerate(replicons)}
    gs = sorted(genes, key=lambda g: (order[g.interval.contig], g.interval.start,
                                      g.interval.end, g.interval.strand, g.locus_tag))  # fmt: skip
    return GeneModel(replicons=dict(replicons), genes=tuple(gs), features=())


def toy_c1() -> tuple[dict[str, str], dict[str, Replicon], GeneModel]:
    nt = back_translate(C1_AA, C1_STOPS)
    reps = {"c1": Replicon("c1", len(nt), "linear")}
    return {"c1": nt}, reps, model(reps, [gene("G1", "c1", 0, 57, "+")])


def toy_c2() -> tuple[dict[str, str], dict[str, Replicon], GeneModel]:
    """300 nt; ``G2`` ``-`` ``[120, 210)``: gene-strand ATG ... TAA, the stop at ``[120, 123)``."""
    rng = np.random.default_rng(8)
    body = "".join(CODON[a] for a in rng.choice(list("ACDEFGHKLMNPQRSTVWY"), size=28))
    gene_strand = "ATG" + body + "TAA"  # 90 nt on the gene strand
    left = "".join(CODON[a] for a in rng.choice(list("ACDEFGHKLMNPQRSTVWY"), size=40))
    right = "".join(CODON[a] for a in rng.choice(list("ACDEFGHKLMNPQRSTVWY"), size=30))
    nt = left + rc(gene_strand) + right
    assert len(nt) == 300 and nt[120:123] == rc("TAA")
    reps = {"c2": Replicon("c2", 300, "linear")}
    return {"c2": nt}, reps, model(reps, [gene("G2", "c2", 120, 210, "-")])


# ---------------------------------------------------------------- parse_bakta work tables


def _write(path: Path, cols: Sequence[tuple[str, Any]], rows: Sequence[tuple[Any, ...]]) -> None:
    arrays = [pa.array([r[i] for r in rows], type=t) for i, (_, t) in enumerate(cols)]
    pq.write_table(pa.Table.from_arrays(arrays, names=[c for c, _ in cols]), path)


_S, _I, _B = pa.string(), pa.int64(), pa.bool_()
GENE_COLS = [
    ("locus_tag", _S), ("gene", _S), ("product", _S), ("feature_type", _S), ("contig", _S),
    ("start", _I), ("end", _I), ("strand", _S), ("partial", _B), ("crosses_origin", _B),
    ("hypothetical", _B), ("aa", _S), ("nt", _S), ("start_codon", _S),
    ("ribosomal_protein", _B), ("amr", _B),
]  # fmt: skip
FEATURE_COLS = [
    ("feature_id", _S), ("type", _S), ("contig", _S), ("start", _I), ("end", _I),
    ("strand", _S), ("gene", _S), ("product", _S), ("amino_acid", _S), ("anticodon", _S),
]  # fmt: skip
REPLICON_COLS = [("contig", _S), ("length", _I), ("topology", _S), ("kind", _S)]


def write_annotation(gm: GeneModel, out: Path) -> dict[str, Path]:
    """``genes``, ``features`` and ``replicons`` Parquet files (the ``parse_bakta`` columns)."""
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for g in gm.genes:
        iv = g.interval
        rep = gm.replicons[iv.contig]
        rows.append((g.locus_tag, g.gene, g.product, g.feature_type, iv.contig, iv.start, iv.end,
                     iv.strand, g.partial, iv.end > rep.length, g.hypothetical, g.aa, g.nt,
                     g.start_codon, False, False))  # fmt: skip
    paths = {k: out / f"{k}.parquet" for k in ("genes", "features", "replicons")}
    _write(paths["genes"], GENE_COLS, rows)
    _write(paths["features"], FEATURE_COLS, [])
    reps = [(r.id, r.length, r.topology, r.kind) for r in gm.replicons.values()]
    _write(paths["replicons"], REPLICON_COLS, reps)
    return paths


def write_genome(genome: Mapping[str, str], path: Path) -> Path:
    with path.open("w") as fh:
        for cid, seq in genome.items():
            fh.write(f">{cid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    return path
