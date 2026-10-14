"""Gene model interface shared by M3 and M2-03 (issue M3-01; GIAP_SPEC.md §6.5, §7.1.2, d16).

``giap.model.genes`` defines ``Gene``, ``NcFeature`` and ``GeneModel`` and one minimal
constructor, ``gene_model_from_bakta_json``, that reads the Bakta-like JSON of giapsim. Bakta JSON
is 1-based closed; everything inside GIAP is 0-based half-open, origin-crossing features on
circular replicons are stored unwrapped (``end > length``) and linear contigs never wrap.

Ambiguities resolved here (spec silent):
- ``by_locus`` on an unknown locus tag raises ``KeyError``.
- unknown contigs, duplicate locus tags and wrapping on linear contigs raise ``ValueError``.
- ``overlapping`` order is not pinned (compared as sets), only its content and determinism.
- ``crosses_origin`` is not a ``Gene`` field; it is ``gene.interval.crosses_origin(replicon)``.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import random
from collections.abc import Mapping
from pathlib import Path
from typing import Any, get_args

import pytest

from giap.model import Interval, Replicon, interval_sort_key
from giap.model.genes import (
    Gene,
    GeneModel,
    GeneType,
    NcFeature,
    NcType,
    gene_model_from_bakta_json,
)

CIRC = {"chr1": Replicon("chr1", 50_000, "circular", "chromosome")}
LIN = {"chr1": Replicon("chr1", 50_000, "linear", "unknown")}
TWO = {
    "chr1": Replicon("chr1", 50_000, "circular", "chromosome"),
    "p1": Replicon("p1", 5_000, "circular", "plasmid"),
}


def _write(path: Path, features: list[dict[str, Any]], reps: Mapping[str, Replicon]) -> Path:
    data = {
        "sequences": [
            {"id": r.id, "length": r.length, "topology": r.topology, "type": r.kind}
            for r in reps.values()
        ],
        "features": features,
    }
    path.write_text(json.dumps(data, indent=1) + "\n")
    return path


def _cds(**kw: Any) -> dict[str, Any]:
    f: dict[str, Any] = {
        "type": "cds",
        "contig": "chr1",
        "start": 101,
        "stop": 400,
        "strand": "+",
        "locus": "GSIM_00005",
        "aa": "M" + "A" * 98,
        "nt": "ATG" + "GCT" * 98 + "TAA",
        "start_type": "ATG",
    }
    f.update(kw)
    return f


def _load(
    tmp_path: Path, features: list[dict[str, Any]], reps: Mapping[str, Replicon] = CIRC
) -> GeneModel:
    return gene_model_from_bakta_json(_write(tmp_path / "bakta.json", features, reps), reps)


def _ids(items: list[Gene | NcFeature]) -> set[str]:
    return {x.locus_tag if isinstance(x, Gene) else x.id for x in items}


# --- types and dataclass shape ---------------------------------------------------------------


def test_literal_vocabularies() -> None:
    assert set(get_args(GeneType)) == {"cds", "sorf", "pseudogene"}
    assert set(get_args(NcType)) == {
        "rrna",
        "trna",
        "tmrna",
        "ncrna",
        "ncrna_region",
        "crispr",
        "oric",
        "orit",
        "oriv",
        "gap",
    }


def test_gene_fields_and_defaults() -> None:
    iv = Interval("chr1", 100, 400, "+")
    g = Gene("GSIM_00005", "cds", iv, "MA", "ATGGCTTAA", "ATG")
    assert [f.name for f in dataclasses.fields(Gene)] == [
        "locus_tag",
        "feature_type",
        "interval",
        "aa",
        "nt",
        "start_codon",
        "gene",
        "product",
        "partial",
        "hypothetical",
        "db_xrefs",
    ]
    assert g.gene is None
    assert g.product == ""
    assert g.partial is False
    assert g.hypothetical is False
    assert g.db_xrefs == ()
    assert g.interval == iv


def test_gene_and_ncfeature_frozen_slots() -> None:
    g = Gene("L1", "cds", Interval("chr1", 0, 9, "+"), "MA", None, None)
    nc = NcFeature("r1", "rrna", Interval("chr1", 10, 20, "+"))
    with pytest.raises(dataclasses.FrozenInstanceError):
        g.aa = "X"  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        nc.id = "x"  # type: ignore[misc]
    assert not hasattr(g, "__dict__")
    assert not hasattr(nc, "__dict__")
    assert [f.name for f in dataclasses.fields(NcFeature)] == ["id", "type", "interval"]
    assert hash(g) == hash(dataclasses.replace(g))


def test_gene_model_frozen(tmp_path: Path) -> None:
    m = _load(tmp_path, [_cds()])
    assert dataclasses.is_dataclass(m)
    with pytest.raises(dataclasses.FrozenInstanceError):
        m.genes = ()  # type: ignore[misc]
    assert isinstance(m.genes, tuple)
    assert isinstance(m.features, tuple)


def test_gene_model_direct_construction() -> None:
    """M2-03's ``parse_bakta_json`` builds the same type directly from its fields."""
    g = Gene("L1", "cds", Interval("chr1", 100, 400, "+"), "MA", None, "ATG")
    nc = NcFeature("t1", "trna", Interval("chr1", 500, 576, "-"))
    m = GeneModel(replicons=CIRC, genes=(g,), features=(nc,))
    assert m.by_locus("L1") == g
    assert set(m.overlapping(Interval("chr1", 0, 1000, "."))) == {g, nc}
    assert m.overlapping(Interval("chr1", 600, 1000, ".")) == []


# --- spec examples 1-9 -----------------------------------------------------------------------


def test_example1_cds_plus(tmp_path: Path) -> None:
    m = _load(tmp_path, [_cds()])
    (g,) = m.genes
    assert g.locus_tag == "GSIM_00005"
    assert g.feature_type == "cds"
    assert g.interval == Interval("chr1", 100, 400, "+")
    assert g.aa == "M" + "A" * 98
    assert g.nt == "ATG" + "GCT" * 98 + "TAA"
    assert g.start_codon == "ATG"
    assert g.partial is False
    assert g.interval.crosses_origin(m.replicons["chr1"]) is False
    assert m.features == ()


def test_example2_sequence_key_equals_contig_key(tmp_path: Path) -> None:
    f = _cds()
    del f["contig"]
    f["sequence"] = "chr1"
    a = _load(tmp_path, [_cds()])
    b = gene_model_from_bakta_json(_write(tmp_path / "b.json", [f], CIRC), CIRC)
    assert a.genes == b.genes


def test_example3_pseudo_flag(tmp_path: Path) -> None:
    m = _load(tmp_path, [_cds(pseudo=True, aa="")])
    (g,) = m.genes
    assert g.feature_type == "pseudogene"
    assert g.aa == ""


def test_pseudogene_key_marks_pseudogene(tmp_path: Path) -> None:
    m = _load(tmp_path, [_cds(pseudogene={"inference": "x"}, aa="MA*K")])
    (g,) = m.genes
    assert g.feature_type == "pseudogene"
    assert g.aa == "MA*K"  # internal stop kept


def test_example4_sorf_minus(tmp_path: Path) -> None:
    m = _load(
        tmp_path,
        [_cds(type="sorf", start=2001, stop=2075, strand="-", locus="GSIM_00010", aa="M" * 24)],
    )
    (g,) = m.genes
    assert g.feature_type == "sorf"
    assert g.interval == Interval("chr1", 2000, 2075, "-")


def test_example5_origin_crossing_circular(tmp_path: Path) -> None:
    m = _load(tmp_path, [_cds(start=49901, stop=200)])
    (g,) = m.genes
    assert g.interval == Interval("chr1", 49900, 50200, "+")
    assert g.interval.crosses_origin(m.replicons["chr1"]) is True


def test_example6_linear_never_wraps(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="GSIM_00005"):
        _load(tmp_path, [_cds(start=49901, stop=200)], LIN)


def test_example7_rrna(tmp_path: Path) -> None:
    m = _load(tmp_path, [{"type": "rRNA", "contig": "chr1", "start": 5001, "stop": 6542,
                          "strand": "+", "id": "r1"}])  # fmt: skip
    assert m.genes == ()
    assert m.features == (NcFeature("r1", "rrna", Interval("chr1", 5000, 6542, "+")),)


def test_example8_giapsim_bakta_like(tmp_path: Path) -> None:
    import giapsim

    res = giapsim.simulate(giapsim.SimConfig(genome_length=50_000), seed=1, outdir=tmp_path)
    reps = {"chr1": Replicon("chr1", 50_000, "circular", "chromosome")}
    m = gene_model_from_bakta_json(res.bakta_json, reps)
    assert len(m.genes) == res.n_genes
    with res.genes_truth.open() as fh:
        truth = {r["locus_tag"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert set(truth) == {g.locus_tag for g in m.genes}
    for g in m.genes:
        t = truth[g.locus_tag]
        assert len(g.aa) == int(t["length_aa"])
        assert g.interval == Interval("chr1", int(t["start"]) - 1, int(t["end"]), t["strand"])
        assert g.start_codon == t["start_codon"]
        assert g.feature_type == "cds"
        assert not g.aa.endswith("*")


def test_example9_overlap_across_origin(tmp_path: Path) -> None:
    m = _load(tmp_path, [_cds(start=49901, stop=200)])
    hits = m.overlapping(Interval("chr1", 49990, 50010, "+"))
    assert hits == [m.genes[0]]


# --- reader details --------------------------------------------------------------------------


def test_nc_type_mapping_case_insensitive(tmp_path: Path) -> None:
    pairs = {
        "tRNA": "trna",
        "tmRNA": "tmrna",
        "rRNA": "rrna",
        "ncRNA": "ncrna",
        "ncRNA-region": "ncrna_region",
        "crispr": "crispr",
        "oriC": "oric",
        "oriT": "orit",
        "oriV": "oriv",
        "gap": "gap",
    }
    feats = [
        {"type": t, "contig": "chr1", "start": 1001 + 100 * i, "stop": 1050 + 100 * i,
         "strand": "+", "id": f"f{i}"}
        for i, t in enumerate(pairs)
    ]  # fmt: skip
    feats.append(_cds(type="CDS"))
    m = _load(tmp_path, feats)
    got = {f.id: f.type for f in m.features}
    assert got == {f"f{i}": v for i, v in enumerate(pairs.values())}
    assert [g.feature_type for g in m.genes] == ["cds"]


def test_unknown_types_ignored(tmp_path: Path) -> None:
    m = _load(tmp_path, [_cds(), {"type": "mystery", "contig": "chr1", "start": 1, "stop": 9,
                                  "strand": "+", "id": "x"}])  # fmt: skip
    assert [g.locus_tag for g in m.genes] == ["GSIM_00005"]
    assert m.features == ()


def test_unknown_contig_is_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="chrX"):
        _load(tmp_path, [_cds(contig="chrX")])


def test_duplicate_locus_tags_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="GSIM_00005"):
        _load(tmp_path, [_cds(), _cds(start=1001, stop=1300)])


def test_terminal_stop_stripped(tmp_path: Path) -> None:
    m = _load(tmp_path, [_cds(aa="MAK*")])
    assert m.genes[0].aa == "MAK"


def test_optional_fields(tmp_path: Path) -> None:
    f = _cds(
        gene="rpsL",
        product="30S ribosomal protein S12",
        db_xrefs=["UniRef:UniRef90_P0A7S3", "SO:0001217"],
    )
    del f["start_type"]
    del f["nt"]
    (g,) = _load(tmp_path, [f]).genes
    assert g.gene == "rpsL"
    assert g.product == "30S ribosomal protein S12"
    assert g.db_xrefs == ("UniRef:UniRef90_P0A7S3", "SO:0001217")
    assert isinstance(g.db_xrefs, tuple)
    assert g.start_codon is None
    assert g.nt is None
    assert g.hypothetical is False


def test_truncated_means_partial(tmp_path: Path) -> None:
    (g,) = _load(tmp_path, [_cds(truncated="5'")]).genes
    assert g.partial is True


@pytest.mark.parametrize(
    ("extra", "expected"),
    [
        ({"product": "hypothetical protein"}, True),
        ({"product": "hypothetical protein YbcX"}, True),
        ({"product": "conserved hypothetical protein"}, False),
        ({"product": "DNA gyrase"}, False),
        ({"product": "DNA gyrase", "hypothetical": True}, True),
        ({}, False),
    ],
)
def test_hypothetical(tmp_path: Path, extra: dict[str, Any], expected: bool) -> None:
    (g,) = _load(tmp_path, [_cds(**extra)]).genes
    assert g.hypothetical is expected


def test_feature_at_contig_end_does_not_cross(tmp_path: Path) -> None:
    (g,) = _load(tmp_path, [_cds(start=49701, stop=50000)]).genes
    assert g.interval == Interval("chr1", 49700, 50000, "+")
    assert g.interval.crosses_origin(CIRC["chr1"]) is False


def test_minus_strand_origin_crossing_unwrapped(tmp_path: Path) -> None:
    (g,) = _load(tmp_path, [_cds(start=49951, stop=150, strand="-")]).genes
    assert g.interval == Interval("chr1", 49950, 50150, "-")
    assert g.interval.crosses_origin(CIRC["chr1"]) is True


def test_nc_feature_origin_crossing(tmp_path: Path) -> None:
    m = _load(tmp_path, [{"type": "crispr", "sequence": "chr1", "start": 49990, "stop": 30,
                          "strand": "+", "id": "c1"}])  # fmt: skip
    assert m.features[0].interval == Interval("chr1", 49989, 50030, "+")


def test_empty_features(tmp_path: Path) -> None:
    m = _load(tmp_path, [])
    assert m.genes == ()
    assert m.features == ()
    assert m.overlapping(Interval("chr1", 0, 50_000, ".")) == []


def test_replicons_kept_in_input_order(tmp_path: Path) -> None:
    reps = {"p1": TWO["p1"], "chr1": TWO["chr1"]}
    m = _load(tmp_path, [_cds()], reps)
    assert list(m.replicons) == ["p1", "chr1"]
    assert m.replicons["chr1"] == TWO["chr1"]


# --- sort order and determinism --------------------------------------------------------------


def _mixed() -> list[dict[str, Any]]:
    return [
        _cds(locus="B", start=1001, stop=1300),
        _cds(locus="A", start=1001, stop=1300, strand="-"),
        _cds(locus="C", start=1001, stop=1300),
        _cds(locus="D", start=1001, stop=1600),
        _cds(locus="E", start=101, stop=400),
        _cds(locus="F", contig="p1", start=11, stop=310),
        _cds(locus="G", sequence="chr1", start=49901, stop=200),
    ]


def test_genes_sorted_contig_order_start_end_strand_locus(tmp_path: Path) -> None:
    reps = {"p1": TWO["p1"], "chr1": TWO["chr1"]}  # p1 first in input FASTA
    for i in range(5):
        feats = _mixed()
        random.Random(i).shuffle(feats)
        m = gene_model_from_bakta_json(_write(tmp_path / f"{i}.json", feats, reps), reps)
        assert [g.locus_tag for g in m.genes] == ["F", "E", "B", "C", "A", "D", "G"]
        order = {c: k for k, c in enumerate(reps)}
        keys = [interval_sort_key(g.interval, order, g.locus_tag) for g in m.genes]
        assert keys == sorted(keys)


def test_deterministic_regardless_of_json_order(tmp_path: Path) -> None:
    feats = [
        *_mixed(),
        {"type": "tRNA", "contig": "chr1", "start": 3001, "stop": 3076, "strand": "-", "id": "t1"},
        {"type": "rRNA", "contig": "chr1", "start": 5001, "stop": 6542, "strand": "+", "id": "r1"},
    ]
    a = gene_model_from_bakta_json(_write(tmp_path / "a.json", feats, TWO), TWO)
    feats.reverse()
    b = gene_model_from_bakta_json(_write(tmp_path / "b.json", feats, TWO), TWO)
    assert a.genes == b.genes
    assert a.features == b.features
    q = Interval("chr1", 0, 50_000, ".")
    assert a.overlapping(q) == b.overlapping(q) == a.overlapping(q)


# --- lookups ---------------------------------------------------------------------------------


def test_by_locus(tmp_path: Path) -> None:
    m = gene_model_from_bakta_json(_write(tmp_path / "a.json", _mixed(), TWO), TWO)
    assert m.by_locus("D").interval == Interval("chr1", 1000, 1600, "+")
    assert m.by_locus("F").interval.contig == "p1"
    with pytest.raises(KeyError):
        m.by_locus("NOPE")


def test_overlapping_half_open_both_strands(tmp_path: Path) -> None:
    m = _load(tmp_path, [
        _cds(locus="P", start=101, stop=400),                 # [100, 400) +
        _cds(locus="M", start=1001, stop=1300, strand="-"),   # [1000, 1300) -
        {"type": "tRNA", "contig": "chr1", "start": 391, "stop": 466, "strand": "-",
         "id": "t1"},                                         # [390, 466) -
    ])  # fmt: skip
    assert _ids(m.overlapping(Interval("chr1", 400, 410, "+"))) == {"t1"}
    assert _ids(m.overlapping(Interval("chr1", 90, 100, "+"))) == set()
    assert _ids(m.overlapping(Interval("chr1", 99, 101, "-"))) == {"P"}
    assert _ids(m.overlapping(Interval("chr1", 399, 400, "."))) == {"P", "t1"}
    assert _ids(m.overlapping(Interval("chr1", 1299, 1300, "+"))) == {"M"}
    assert _ids(m.overlapping(Interval("chr1", 1300, 2000, "-"))) == set()
    assert _ids(m.overlapping(Interval("chr1", 0, 50_000, "."))) == {"P", "M", "t1"}


def test_overlapping_other_contig_not_returned(tmp_path: Path) -> None:
    m = gene_model_from_bakta_json(_write(tmp_path / "a.json", _mixed(), TWO), TWO)
    assert _ids(m.overlapping(Interval("p1", 0, 5_000, "."))) == {"F"}


def test_overlapping_origin_aware(tmp_path: Path) -> None:
    m = _load(tmp_path, [
        _cds(locus="X", start=49901, stop=200),   # [49900, 50200) = [49900,50000) + [0,200)
        _cds(locus="Y", start=301, stop=600),     # [300, 600)
        _cds(locus="Z", start=49001, stop=49300), # [49000, 49300)
    ])  # fmt: skip
    assert _ids(m.overlapping(Interval("chr1", 100, 150, "-"))) == {"X"}  # wrapped part
    assert _ids(m.overlapping(Interval("chr1", 49950, 49960, "+"))) == {"X"}
    assert _ids(m.overlapping(Interval("chr1", 200, 300, "+"))) == set()  # half-open gap
    # origin-crossing query [49250, 50350) = [49250, 50000) + [0, 350)
    assert _ids(m.overlapping(Interval("chr1", 49250, 50350, "+"))) == {"X", "Y", "Z"}
    # origin-crossing query touching only the wrapped part
    assert _ids(m.overlapping(Interval("chr1", 49950, 50001, "+"))) == {"X"}
