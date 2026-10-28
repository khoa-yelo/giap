"""M7-12 lead decision 23: makedb ORF-entry decoys leave gene-level pairing.

A decoy whose every stripped ID is a makedb six-frame entry (``GIPPAORF_…``) belongs to the M7-11
class FDR, not to a gene: like a contaminant decoy (M4-12 decision 12) it is neither paired nor
unpaired, so a makedb search keeps the ``decoy_picked`` tier and the gene q values. A decoy that
mixes a gene ID with an ORF entry pairs on the gene ID.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_fixtures as fx
from detection_tiers_fixtures import Builder

ORF_IDS = [f"rev_GIPPAORF_intergenic_chr1_{100 + 30 * i}_{190 + 30 * i}_p" for i in range(3)]


def _eight_gene_decoys(extra: list[list[str]]) -> Builder:
    b = Builder(10)
    t = b.tags or []
    for i, g in enumerate(t):
        b.target(g, 0.001, 0.01, 100.0 + i)
    for i in range(8):
        b.decoy([f"rev_{t[i]}"], 0.001, 0.01, 1.0 + i)
    for i, ids in enumerate(extra):
        b.decoy(ids, 0.001, 0.01, 20.0 + i)
    return b


def test_orf_entry_decoys_are_neither_paired_nor_unpaired() -> None:
    r = _eight_gene_decoys([[x] for x in ORF_IDS]).run()
    assert r.det.tier == "decoy_picked"
    assert "decoy_pairing_low" not in r.det.warnings
    assert fx.tier_metric(r.det, "n_decoys_paired") == 8
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0


def test_orf_entry_decoys_do_not_change_gene_q() -> None:
    base = _eight_gene_decoys([]).run()
    with_orf = _eight_gene_decoys([[x] for x in ORF_IDS]).run()
    pd.testing.assert_frame_equal(with_orf.det.table, base.det.table)


def test_a_decoy_mixing_a_gene_and_an_orf_entry_pairs_on_the_gene() -> None:
    b = Builder(10)
    t = b.tags or []
    mixed = _eight_gene_decoys([[f"rev_{t[9]}", ORF_IDS[0]]]).run()
    assert fx.tier_metric(mixed.det, "n_decoys_paired") == 9
    assert fx.tier_metric(mixed.det, "n_decoys_unpaired") == 0
