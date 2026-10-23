"""M7-12 lead decision 28 (amends decision 23; code review #225 B1).

An all-``GIAPORF_`` decoy leaves gene-level pairing only when none of its stripped IDs maps to a
gene in the protein-gene map. An ``nterm_extension`` entry is voted to its gene by its target
peptides, so its decoy pairs with that gene, as on main; otherwise the gene's target would compete
without its decoy (anti-conservative). An all-ORF decoy that maps to no gene is still excluded.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd
from hypothesis import given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_fixtures as fx
import evidence_oracle as eo
from detection_tiers_fixtures import Builder

from giap.genetable.detection import pair_decoys

EXT = "GIAPORF_nterm_extension_chr1_50_190_p"
UNMAPPED = [f"GIAPORF_intergenic_chr1_{5000 + 30 * i}_{5090 + 30 * i}_p" for i in range(3)]


def _nterm_builder(ext_decoy: str, extra: list[list[str]] | None = None) -> Builder:
    """The reviewer's repro: ten genes, an extension-only target placed in gene 0 (its
    ``protein_ids`` is the ORF entry only), eight gene decoys and one best-scoring decoy."""
    b = Builder(10)
    t = b.tags or []
    for i, g in enumerate(t):
        b.target(g, 0.001, 0.01, 100.0 + i)
    b.target(t[0], 0.001, 0.01, 5.0, protein_ids=[EXT])
    for i in range(1, 9):
        b.decoy([f"rev_{t[i]}"], 0.001, 0.01, 1.0 + i)
    b.decoy([ext_decoy], 0.001, 0.01, 500.0)
    for i, ids in enumerate(extra or []):
        b.decoy(ids, 0.001, 0.01, 20.0 + i)
    return b


def test_nterm_extension_decoy_pairs_with_its_gene() -> None:
    r = _nterm_builder(f"rev_{EXT}").run()
    assert r.det.tier == "decoy_picked"
    assert fx.tier_metric(r.det, "n_decoys_paired") == 9
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0


def test_nterm_extension_decoy_wins_its_gene_as_on_main() -> None:
    """On main the decoy pairs and outscores gene 0, so gene 0 gets no q (NA), not 0.1."""
    r = _nterm_builder(f"rev_{EXT}").run()
    t = Builder(10).tags or []
    assert math.isnan(r.q(t[0])), r.q(t[0])


def test_nterm_extension_decoy_acts_like_a_decoy_of_its_gene() -> None:
    """Same gene table as when the best decoy names gene 0's own protein."""
    t = Builder(10).tags or []
    ext = _nterm_builder(f"rev_{EXT}").run()
    gene = _nterm_builder(f"rev_{t[0]}").run()
    pd.testing.assert_frame_equal(ext.det.table, gene.det.table)


def test_unmapped_orf_decoys_stay_excluded_next_to_a_mapped_one() -> None:
    r = _nterm_builder(f"rev_{EXT}", [[f"rev_{x}"] for x in UNMAPPED]).run()
    assert r.det.tier == "decoy_picked"
    assert "decoy_pairing_low" not in r.det.warnings
    assert fx.tier_metric(r.det, "n_decoys_paired") == 9
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0
    base = _nterm_builder(f"rev_{EXT}").run()
    pd.testing.assert_frame_equal(r.det.table, base.det.table)


# ------------------------------------------------------------------- pair_decoys directly


def _pmap(rows: list[tuple[str, str]]) -> pd.DataFrame:
    return pd.DataFrame({"protein_id": [p for p, _ in rows], "gene": [g for _, g in rows]})


def _decoys(id_lists: list[list[str]]) -> pd.DataFrame:
    rows = [
        eo.pep(fx.peptide_name(k + 1), 0.001, score=1.0 + k, is_decoy=True, protein_ids=ids)
        for k, ids in enumerate(id_lists)
    ]
    return eo.peptides_frame(rows)


def test_pair_decoys_maps_an_orf_entry_voted_to_a_gene() -> None:
    pmap = _pmap([("G1", "g1"), (EXT, "g1")])
    out = pair_decoys(_decoys([[f"rev_{EXT}"], [f"rev_{UNMAPPED[0]}"]]), pmap, ["rev_"])
    assert len(out) == 1
    assert out["gene"].tolist() == ["g1"]


ORF_POOL = [f"GIAPORF_nterm_extension_chr1_{10 * i}_{10 * i + 90}_p" for i in range(6)]
GENES = ["g0", "g1", "g2"]


@settings(max_examples=60, deadline=None)
@given(
    voted=st.dictionaries(st.sampled_from(ORF_POOL), st.sampled_from(GENES), max_size=6),
    decoys=st.lists(st.sampled_from(ORF_POOL), min_size=1, max_size=12),
)
def test_single_orf_decoy_pairs_iff_its_entry_is_voted(
    voted: dict[str, str], decoys: list[str]
) -> None:
    """Property: a single-ID ORF decoy gets a row iff its entry is in ``pmap``, and that row
    pairs with the voted gene; unvoted ones get no row (neither paired nor unpaired)."""
    pmap = _pmap([("G0", "g0"), *voted.items()])
    out = pair_decoys(_decoys([[f"rev_{x}"] for x in decoys]), pmap, ["rev_"])
    want = [voted[x] for x in decoys if x in voted]
    assert out["gene"].tolist() == want
