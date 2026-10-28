"""M2-20 ``annotation/genes.gff3``: the merged model with source-tagged layers (spec table
"`annotation/genes.gff3` (M1-05 writer), one source per layer in column 2"; example 7;
GIPPA_SPEC.md §7.1.8, §9.5 "a GFF3 parser round trip").

The file written by ``context_layers`` on the SIM100K producer chain (``helpers.layers_chain``)
is read back with the M1-05 ``read_gff3``.

Settled by the test author (lead question 1 proposes the SO types used below):

* genes: source ``Bakta``; type ``CDS`` for ``cds`` and ``sorf`` (phase 0), ``pseudogene`` for
  pseudogenes (no phase); attributes ``ID`` = ``locus_tag`` = the locus tag, ``feature_type``,
  ``gene`` / ``product`` when set. Coordinates equal ``genes.parquet``.
* IDs are unique per file (lead decision 8, §9.5): Bakta's gene IDs as given (the locus tag);
  every ID GIPPA writes carries its feature kind as a prefix — ``operon:<tu_id>``,
  ``promoter:<promoters.bed name>``, ``terminator:<terminators.bed name>``,
  ``region:<region_id>`` — so a region ``plasmid_1`` never collides with the contig
  ``plasmid_1``'s ``Is_circular`` region line (``ID=plasmid_1``, M1-05 writer).
* TUs: type ``operon``, source ``UniOP`` or ``distance_model``, ``ID=operon:TU_<n>``,
  ``genes`` = the TU's genes in transcription order, ``upstream=<start1>-<end1>`` (1-based closed
  of ``upstream_start`` / ``upstream_end``, unwrapped).
* promoters: source ``Promotech``, type ``promoter``, one feature per merged ``promoters.bed``
  record with the same interval; ``ID=promoter:promoter_<n>``; column 6 holds the tool's own
  score — the maximum Promotech probability of the merged windows, not the BED score / 1000
  (lead decision 12).
* terminators: type ``terminator``, column 2 the terminator's ``source``,
  ``ID=terminator:<source>_<n>`` (the ``terminators.bed`` name), ``confidence`` attribute.
* regions: ``prophage`` / ``plasmid`` / ``biosynthetic_gene_cluster`` / ``CRISPR`` for
  prophage / plasmid / bgc / crispr, ``ID=region:<region_id>``, ``genes`` when the region has
  members; AMR regions are not written.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from helpers import layers_oracle as orc
from helpers.layers_chain import Chain, build_chain, fresh_steps

from gippa.config import GippaConfig
from gippa.io.gff3 import GffFeature, read_gff3
from gippa.io.tsv import format_value

REGION_TYPES = {
    "prophage": "prophage", "plasmid": "plasmid", "bgc": "biosynthetic_gene_cluster",
    "crispr": "CRISPR",
}  # fmt: skip
TERM_SOURCES = {"TerminatorNet", "TransTermHP", "BacTermFinder"}


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(tmp_path_factory.mktemp("layers_gff_chain"))


@pytest.fixture(scope="module")
def gff(chain: Chain, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("layers_gff")
    outs = {"annotation": out / "annotation", "layers": out / "layers.parquet"}
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        reg.run_step("context_layers", GippaConfig(seed=1), dict(chain.inputs), outs)
    return outs["annotation"] / "genes.gff3"


@pytest.fixture(scope="module")
def feats(chain: Chain, gff: Path) -> list[GffFeature]:
    return read_gff3(gff, orc.replicons(chain.frame("replicons")))


def _attr(f: GffFeature, key: str) -> list[str]:
    v = f.attributes.get(key)
    if v is None:
        return []
    return [v] if isinstance(v, str) else list(v)


def _one(f: GffFeature, key: str) -> str:
    v = _attr(f, key)
    assert len(v) == 1, (key, v)
    return v[0]


def test_sequence_region_pragmas_and_circular_regions(chain: Chain, gff: Path) -> None:
    lines = gff.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "##gff-version 3"
    for r in chain.frame("replicons").itertuples(index=False):
        assert f"##sequence-region {r.contig} 1 {r.length}" in lines
        circ = [x for x in lines if x.startswith(f"{r.contig}\t") and "Is_circular=true" in x]
        assert (len(circ) == 1) == (r.topology == "circular"), r.contig


def test_example_7_gene_features_round_trip(chain: Chain, feats: list[GffFeature]) -> None:
    genes = chain.frame("genes")
    bakta = [f for f in feats if f.source == "Bakta" and f.type in ("CDS", "pseudogene")]
    assert len(bakta) == len(genes)
    by_id = {_one(f, "ID"): f for f in bakta}
    for g in genes.itertuples(index=False):
        f = by_id[g.locus_tag]
        assert (f.interval.contig, f.interval.start, f.interval.end, f.interval.strand) == (
            g.contig, int(g.start), int(g.end), g.strand,
        )  # fmt: skip
        want_type = "pseudogene" if g.feature_type == "pseudogene" else "CDS"
        assert f.type == want_type, g.locus_tag
        assert f.phase == (0 if want_type == "CDS" else None), g.locus_tag
        assert _one(f, "locus_tag") == g.locus_tag
        assert _one(f, "feature_type") == g.feature_type
        if isinstance(g.gene, str) and g.gene:
            assert _one(f, "gene") == g.gene
        if isinstance(g.product, str) and g.product:
            assert ",".join(_attr(f, "product")) == g.product


def test_tu_features(chain: Chain, feats: list[GffFeature]) -> None:
    tus = chain.frame("tus")
    ops = {_one(f, "ID"): f for f in feats if f.type == "operon"}
    assert set(ops) == {f"operon:{t}" for t in tus["tu_id"]}
    for t in tus.itertuples(index=False):
        f = ops[f"operon:{t.tu_id}"]
        assert f.source in {"UniOP", "distance_model"}
        assert (f.interval.contig, f.interval.start, f.interval.end, f.interval.strand) == (
            t.contig, int(t.start), int(t.end), t.strand,
        )  # fmt: skip
        assert _attr(f, "genes") == list(t.genes)
        assert _one(f, "upstream") == f"{int(t.upstream_start) + 1}-{int(t.upstream_end)}"


def test_promoter_features_follow_promoters_bed(chain: Chain, feats: list[GffFeature]) -> None:
    from gippa.annotate.promoters import promoters_bed

    reps = orc.replicons(chain.frame("replicons"))
    yes = GippaConfig().annotation.promoters.promotech_yes
    recs = promoters_bed(chain.frame("promoter_windows"), reps, yes)
    prom = [f for f in feats if f.type == "promoter"]
    assert prom and len(prom) == len(recs)
    assert {f.source for f in prom} == {"Promotech"}
    got = {_one(f, "ID"): f for f in prom}
    tool = orc.promoter_record_scores(chain.frame("promoter_windows"), recs, reps, yes)
    assert any(round(1000 * x) / 1000 != x for x in tool.values()), "full-precision scores"
    for r in recs:
        f = got[f"promoter:{r.name}"]
        assert f.interval == r.interval
        assert f.score == float(format_value(tool[r.name])), r.name  # not r.score / 1000


def test_terminator_features(chain: Chain, feats: list[GffFeature]) -> None:
    from gippa.annotate.terminators import terminators_bed

    reps = orc.replicons(chain.frame("replicons"))
    recs = terminators_bed(orc.terminators(chain.frame("terminators")), reps)
    term = [f for f in feats if f.type == "terminator"]
    assert len(term) == len(recs)
    assert {f.source for f in term} <= TERM_SOURCES
    got = {_one(f, "ID"): f for f in term}
    for r in recs:
        f = got[f"terminator:{r.name}"]
        assert f.interval == r.interval
        assert f.source == r.name.rsplit("_", 1)[0]
        assert abs(float(_one(f, "confidence")) * 1000 - r.score) <= 0.5 + 1e-9


def test_region_features(chain: Chain, feats: list[GffFeature]) -> None:
    regions = chain.frame("regions")
    by_id = {_one(f, "ID"): f for f in feats if f.type in REGION_TYPES.values()}
    for r in regions.itertuples(index=False):
        if r.region_type == "amr":
            assert f"region:{r.region_id}" not in by_id
            continue
        f = by_id[f"region:{r.region_id}"]
        assert f.type == REGION_TYPES[r.region_type]
        assert (f.interval.contig, f.interval.start, f.interval.end) == (
            r.contig, int(r.start), int(r.end),
        )  # fmt: skip
        assert _attr(f, "genes") == list(r.genes if r.genes is not None else [])
    assert {"region:prophage_1", "region:plasmid_1", "region:bgc_1", "region:crispr_1"} <= set(
        by_id
    )


def test_no_other_feature_types(feats: list[GffFeature]) -> None:
    allowed = {"CDS", "pseudogene", "operon", "promoter", "terminator", *REGION_TYPES.values()}
    assert {f.type for f in feats} <= allowed


def test_gff3_layers_without_optional_inputs(chain: Chain, tmp_path: Path) -> None:
    """Only genes and TUs when every optional layer is missing (`genome` is required, dec. 7)."""
    outs = {"annotation": tmp_path / "annotation", "layers": tmp_path / "layers.parquet"}
    keep = {k: chain.inputs[k] for k in ("genes", "replicons", "pairs_raw", "tus", "genome")}
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        reg.run_step("context_layers", GippaConfig(seed=1), keep, outs)
    got = read_gff3(outs["annotation"] / "genes.gff3", orc.replicons(chain.frame("replicons")))
    assert {f.type for f in got} == {"CDS", "pseudogene", "operon"}
    n = len([f for f in got if f.source == "Bakta" and f.type in ("CDS", "pseudogene")])
    assert n == len(pd.read_parquet(chain.inputs["genes"]))


# ----------------------------------------------------------------------------- unique IDs (dec. 8)


def _raw_attrs(gff: Path) -> list[tuple[str, dict[str, str]]]:
    """(type, column-9 attributes) of every feature line, as written (no model)."""
    out = []
    for line in gff.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        attrs = dict(p.split("=", 1) for p in cols[8].split(";") if p)
        out.append((cols[2], attrs))
    return out


def test_ids_unique_per_file_with_a_contig_named_like_a_region(chain: Chain, gff: Path) -> None:
    """Lead decision 8 (§9.5): sim100k has a circular contig ``plasmid_1`` (region line
    ``ID=plasmid_1``) and a plasmid region ``plasmid_1``; no ID repeats and every Parent=
    reference resolves."""
    assert "plasmid_1" in set(chain.frame("replicons")["contig"])
    assert "plasmid_1" in set(chain.frame("regions")["region_id"])
    rows = _raw_attrs(gff)
    ids = [a["ID"] for _, a in rows if "ID" in a]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    assert not dup, f"duplicate GFF3 IDs: {dup}"
    known = set(ids)
    for _, a in rows:
        for parent in a.get("Parent", "").split(",") if "Parent" in a else []:
            assert parent in known, f"unresolved Parent={parent}"


def test_ids_carry_the_feature_kind(chain: Chain, gff: Path) -> None:
    """GIPPA-written IDs: ``<kind>:<id>``; Bakta gene IDs are the locus tags as given."""
    prefix = {"operon": "operon:", "promoter": "promoter:", "terminator": "terminator:",
              **{t: "region:" for t in REGION_TYPES.values()}}  # fmt: skip
    tags = set(chain.frame("genes")["locus_tag"])
    contigs = set(chain.frame("replicons")["contig"])
    seen = set()
    for ftype, a in _raw_attrs(gff):
        fid = a["ID"]
        if ftype in ("CDS", "pseudogene"):
            assert fid in tags, fid
        elif ftype == "region":
            assert fid in contigs, fid  # the writer's Is_circular line
        else:
            assert fid.startswith(prefix[ftype]), (ftype, fid)
            seen.add(ftype)
    assert {"operon", "promoter", "terminator", "prophage", "plasmid"} <= seen
