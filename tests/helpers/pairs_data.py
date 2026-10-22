"""Hand-built inputs for the M6-04 pair table (docs/specs/M6-04-pair-table.md; test author).

``BASE`` is the spec's example layout on one circular 10 kb contig ``c1``: G1 + [100, 400),
G2 + [450, 900), G3 - [950, 1400), G4 - [1500, 1900), G5 + [2000, 2400). ``pairs_raw`` rows are
keyed in genome order (left, right), as M2-12 ``build_pairs`` writes them, wrap pair included.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd

from giap.config import GiapConfig
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

L = 10_000
BASE = (
    ("G1", 100, 400, "+"), ("G2", 450, 900, "+"), ("G3", 950, 1400, "-"),
    ("G4", 1500, 1900, "-"), ("G5", 2000, 2400, "+"),
)  # fmt: skip
TEST_COLS = (
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "null_model",
    "n_perm", "lfdr",
)  # fmt: skip
PAIR_COLS = (
    "pair_id", "class", "gene_a", "gene_b", "strand", "intergenic_bp", "uniop_p",
    "promoter_between", "tx_rate_between", "terminator_between", "detected_a", "detected_b",
    "p_detect_a", "p_detect_b", "abundance_log2_a", "abundance_log2_b", "delta_log2", "test",
    *TEST_COLS, "flags",
)  # fmt: skip


def gene(tag: str, start: int, end: int, strand: str, contig: str = "c1", **kw: Any) -> Gene:
    iv = Interval(contig, start, end, strand)  # type: ignore[arg-type]
    return Gene(tag, "cds", iv, "M", None, "ATG", **kw)


def genes_of(rows: Sequence[tuple[str, int, int, str]], contig: str = "c1") -> list[Gene]:
    return [gene(t, s, e, st, contig) for t, s, e, st in rows]


def reps(*specs: tuple[str, int, str]) -> dict[str, Replicon]:
    return {c: Replicon(c, n, topo) for c, n, topo in specs}  # type: ignore[arg-type]


def ring(length: int = L) -> dict[str, Replicon]:
    return reps(("c1", length, "circular"))


def pairs_raw(probs: Mapping[tuple[str, str], float]) -> pd.DataFrame:
    rows = [(a, b, float(p), "uniop") for (a, b), p in probs.items()]
    return pd.DataFrame(rows, columns=["gene_a", "gene_b", "uniop_p", "source"])


def gene_tbl(
    genes: Sequence[Gene], abundance: Mapping[str, float] | None = None, **cols: Mapping[str, Any]
) -> pd.DataFrame:
    ab = abundance or {}
    tags = [g.locus_tag for g in genes]
    det = cols.get("detected", {})
    pdet = cols.get("p_detect", {})
    return pd.DataFrame(
        {
            "locus_tag": tags,
            "detected": [det.get(t, t in ab) for t in tags],
            "p_detect": [pdet.get(t, 0.9) for t in tags],
            "abundance_log2": [ab.get(t, float("nan")) for t in tags],
        }
    )


def cfg(threshold: float = 0.5) -> GiapConfig:
    return GiapConfig.model_validate(
        {"seed": 1, "annotation": {"uniop": {"pair_threshold": threshold}}}
    )


def rna_frame(
    rows: Sequence[tuple[str, str, int, int, str]] = (), contig: str = "c1"
) -> pd.DataFrame:
    return pd.DataFrame(
        [(i, t, contig, s, e, st) for i, t, s, e, st in rows],
        columns=["feature_id", "type", "contig", "start", "end", "strand"],
    )


def evidence(
    windows: pd.DataFrame | None = None,
    terminators: Sequence[Any] | None = None,
    tx: pd.DataFrame | None = None,
    rna: pd.DataFrame | None = None,
) -> Any:
    from giap.operons.pairs import Evidence

    return Evidence(
        windows=windows,
        terminators=None if terminators is None else tuple(terminators),
        tx=tx,
        rna=rna_frame() if rna is None else rna,
    )


def frame(
    genes: Sequence[Gene],
    replicons: Mapping[str, Replicon],
    probs: Mapping[tuple[str, str], float],
    ev: Any = None,
    tbl: pd.DataFrame | None = None,
    threshold: float = 0.5,
) -> pd.DataFrame:
    from giap.operons.pairs import pair_frame

    return pair_frame(
        gene_tbl(genes) if tbl is None else tbl, list(genes), replicons, pairs_raw(probs),
        evidence() if ev is None else ev, cfg(threshold),
    )  # fmt: skip


def row(df: pd.DataFrame, a: str, b: str, cls: str | None = None) -> pd.Series:
    sel = (df["gene_a"] == a) & (df["gene_b"] == b)
    if cls is not None:
        sel &= df["class"] == cls
    hit = df[sel]
    assert len(hit) == 1, f"expected one {cls or ''} row {a}->{b}, got {len(hit)}:\n{df}"
    return hit.iloc[0]
