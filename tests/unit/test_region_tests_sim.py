"""giapsim layout extension owed to M4-17: per-contig region placement (``RegionSpec.contig``).

docs/specs/M4-17-region-tests.md "Note from M4-07 review" and M4-07 lead decision 6: a region can
be placed on a named replicon (e.g. a 60-gene plasmid with its own AMR set); ``contig=None`` (the
default) keeps the chromosome. A plasmid AMR gene is then in two regions (``amr_<n>`` and
``plasmid_<n>``, M2-20 example 4). Region IDs stay ``<type>_<n>`` numbered by contig order then
start (§6.5).
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

import giapsim.layout as la


def _cfg(**kw: object) -> object:
    return la.LayoutConfig(**kw)


def test_region_spec_has_contig_defaulting_to_chromosome() -> None:
    fields = {f.name: f for f in dataclasses.fields(la.RegionSpec)}
    assert "contig" in fields
    assert la.RegionSpec("bgc", 5).contig is None


def test_scattered_set_on_plasmid() -> None:
    spec = la.RegionSpec("amr", 6, scattered=True, contig="p1")
    cfg = _cfg(n_genes_chromosome=300, plasmids=(60,), regions=(*la.DEFAULT_REGIONS, spec))
    lay = la.simulate_layout(cfg, np.random.default_rng(1))
    reg = lay.regions.set_index("region_id")
    amr = reg[reg["region_type"] == "amr"]
    assert sorted(amr["contig"]) == ["chr1", "p1"]
    # numbered by contig order: the chromosome set first
    assert amr.loc["amr_1", "contig"] == "chr1" and amr.loc["amr_2", "contig"] == "p1"
    contig_of = dict(zip(lay.genes["locus_tag"], lay.genes["contig"], strict=True))
    members = lay.region_genes.loc[lay.region_genes["region_id"] == "amr_2", "locus_tag"]
    assert len(members) == 6 and {contig_of[g] for g in members} == {"p1"}
    chr_members = lay.region_genes.loc[lay.region_genes["region_id"] == "amr_1", "locus_tag"]
    assert {contig_of[g] for g in chr_members} == {"chr1"}
    # every plasmid AMR gene is in amr_2 and plasmid_1 (region_ids sorted, M4-07 decision 3)
    truth = lay.truth.set_index("locus_tag")
    for g in members:
        assert list(truth.loc[g, "region_ids"]) == ["amr_2", "plasmid_1"]
    assert lay.genes.set_index("locus_tag").loc[list(members), "amr"].all()


def test_block_region_on_plasmid_is_consecutive() -> None:
    spec = la.RegionSpec("bgc", 12, contig="p1")
    cfg = _cfg(n_genes_chromosome=200, plasmids=(40,), regions=(spec,))
    lay = la.simulate_layout(cfg, np.random.default_rng(2))
    members = set(lay.region_genes.loc[lay.region_genes["region_id"] == "bgc_1", "locus_tag"])
    p1 = lay.genes[lay.genes["contig"] == "p1"].reset_index(drop=True)
    pos = sorted(int(i) for i in np.flatnonzero(p1["locus_tag"].isin(members)))
    assert len(pos) == 12
    n = len(p1)
    # one run of consecutive genes (allowing a wrap on the circular plasmid)
    gaps = sum(1 for a, b in zip(pos, [*pos[1:], pos[0] + n], strict=True) if b - a != 1)
    assert gaps == 1
    assert lay.regions.set_index("region_id").loc["bgc_1", "contig"] == "p1"


def test_unknown_contig_rejected() -> None:
    spec = la.RegionSpec("amr", 3, scattered=True, contig="p9")
    with pytest.raises(ValueError, match="p9"):
        la.simulate_layout(_cfg(plasmids=(30,), regions=(spec,)), np.random.default_rng(0))


def test_capacity_checked_per_contig() -> None:
    spec = la.RegionSpec("bgc", 31, contig="p1")
    with pytest.raises(ValueError):
        la.simulate_layout(_cfg(plasmids=(30,), regions=(spec,)), np.random.default_rng(0))


def test_default_regions_stay_on_chromosome() -> None:
    lay = la.simulate_layout(_cfg(plasmids=(30,)), np.random.default_rng(3))
    reg = lay.regions[lay.regions["region_type"] != "plasmid"]
    assert set(reg["contig"]) == {"chr1"}
