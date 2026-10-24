"""K-12 plot survey, batch 1 item 1: the V13 ORF-peptide overlay matches on an I/L-folded key.

``locus._draw_orf`` matched ``novel_orfs.peptides`` to ``peptide_hits.peptide`` by exact string,
so a peptide that differs only in I/L was not overlaid (59 % missed on K-12). The house
convention (``giap.mapping``: ``seq.upper().replace("I", "L")``) folds I into L; the overlay
must use the same key, in both directions. Other residue differences still do not match.

Assumed (as the V13 tests do): the ORF's peptides are highlighted on the ``"pileup"`` axes with
artists coloured ``#CC79A7``; a highlighted run spans the placement's ``[start, end)``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import numpy as np
import pytest
from helpers import viz_locus as vl
from matplotlib.collections import Collection
from matplotlib.figure import Figure

from giap.config import GiapConfig, resolve_config
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Skip

ORF = "#CC79A7"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


def _orf(peptides: str) -> dict[str, object]:
    """One called ``+`` novel ORF on chr1 [5001, 5600] (1-based closed)."""
    return {
        "orf_id": "NOVEL_chr1_5001_p", "class": "intergenic", "contig": "chr1", "strand": "+",
        "frame": 1, "start_longest": 5_001, "start_best_ostir": 5_001, "end": 5_600,
        "length_aa_longest": 199, "length_aa_best": 199, "start_codon_best": "ATG",
        "n_peptides": len(peptides.split(";")), "peptides": peptides, "best_peptide_q": 0.0001,
        "class_fdr": "0.004", "called": True, "rank": 1,
    }  # fmt: skip


def _figure(tmp: Path, cfg: GiapConfig, peptides: str, hits: list[dict[str, object]]) -> Figure:
    import giap.viz.novel_orfs  # noqa: F401  (registers V13)

    out = vl.build_out(tmp, novel_orfs=[_orf(peptides)], hits=hits)
    result = base.registered()["V13"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    (view,) = [i for i in result if i.suffix and i.suffix.startswith("locus_")]
    return view.figure


def _hex(c: Any) -> str | None:
    try:
        return mcolors.to_hex(c).upper()
    except (TypeError, ValueError):
        return None


def _orf_spans(fig: Figure) -> list[tuple[float, float]]:
    """x extents of the ``#CC79A7`` artists on the pileup axes."""
    (pile,) = [ax for ax in fig.axes if ax.get_label() == "pileup"]
    spans = []
    for c in pile.collections:
        if not isinstance(c, Collection):
            continue
        cols = {_hex(x) for x in [*c.get_facecolors(), *c.get_edgecolors()]}
        if ORF not in cols:
            continue
        for p in c.get_paths():
            v = np.asarray(p.vertices)
            if len(v):
                spans.append((float(v[:, 0].min()), float(v[:, 0].max())))
    for p in pile.patches:
        if ORF in {_hex(p.get_facecolor()), _hex(p.get_edgecolor())}:
            e = p.get_extents()
            x = pile.transData.inverted().transform([[e.x0, e.y0], [e.x1, e.y1]])
            spans.append((float(x[0, 0]), float(x[1, 0])))
    return spans


def _covered(spans: list[tuple[float, float]], s0: int, e: int) -> bool:
    return any(a <= s0 + 0.5 and b >= e - 0.5 for a, b in spans)


def test_orf_with_i_overlays_a_placement_with_l(tmp_path: Path, cfg: GiapConfig) -> None:
    hits = [vl.hit("PEPTLDEK", "chr1", 5_101, 5_124)]
    spans = _orf_spans(_figure(tmp_path, cfg, "PEPTIDEK", hits))
    assert _covered(spans, 5_100, 5_124), spans


def test_orf_with_l_overlays_a_placement_with_i(tmp_path: Path, cfg: GiapConfig) -> None:
    hits = [vl.hit("SAMIQERK", "chr1", 5_201, 5_224)]
    spans = _orf_spans(_figure(tmp_path, cfg, "SAMLQERK", hits))
    assert _covered(spans, 5_200, 5_224), spans


def test_every_il_variant_of_an_orf_peptide_is_overlaid(tmp_path: Path, cfg: GiapConfig) -> None:
    """Two placements of two I/L spellings and one exact match: all three runs are purple."""
    hits = [
        vl.hit("LLKPEPR", "chr1", 5_101, 5_121),
        vl.hit("ILKPEPR", "chr1", 5_301, 5_321),
        vl.hit("GGAWNK", "chr1", 5_401, 5_418),
    ]
    spans = _orf_spans(_figure(tmp_path, cfg, "IIKPEPR;GGAWNK", hits))
    for s0, e in ((5_100, 5_121), (5_300, 5_321), (5_400, 5_418)):
        assert _covered(spans, s0, e), (s0, e, spans)


def test_other_residue_differences_do_not_match(tmp_path: Path, cfg: GiapConfig) -> None:
    """Guard: only I/L is folded (e.g. K/Q, the other near-isobaric pair, is not)."""
    hits = [vl.hit("PEPTIDEQ", "chr1", 5_101, 5_124)]
    spans = _orf_spans(_figure(tmp_path, cfg, "PEPTIDEK", hits))
    assert not _covered(spans, 5_100, 5_124), spans
