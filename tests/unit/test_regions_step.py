"""Step ``region_layer`` (M2-19; M2-index step contract row ``region_layer``; M2-18 decision 3).

The step is registered by M2-19 and must match annotate.nf's ``stepIn('region_layer', core +
optional)`` keys: required ``genes``, ``features``, ``replicons`` (the ``core`` list), optional
``genomad`` and ``antismash`` (directories), outputs ``regions`` → ``regions.parquet`` and
``region_genes`` → ``region_genes.parquet``; ``annotation.precomputed.regions`` replaces it.

Inputs are real producer outputs: ``parse_bakta`` (M2-04) run on giapsim's SIM100K (the genome
of the recorded sim100k runs), the recorded sim100k geNOMAD directory (``genome_summary/`` with
prefix ``genome``, header-only) and antiSMASH JSON (no areas); and, for the non-empty case, a
geNOMAD directory and antiSMASH JSON written in exactly the recorded formats (headers copied from
the recordings, area/feature keys as in ``tests/data/recorded/antismash``).

Pinned where the spec is silent:

* ``regions.parquet`` columns ⊇ the ``Region`` fields (``start``/``end`` 0-based half-open, null
  for AMR; ``score`` null for NA; ``genes`` a list or a ``;``-joined string);
  ``region_genes.parquet`` columns ⊇ the ``RegionGene`` fields.
* The antiSMASH input is the directory ``antismash/``; the step reads its single top-level
  ``*.json`` (real runs: ``bakta.json`` next to ``*.gbk``/``*.zip``; the NF stub:
  ``<meta_id>.json`` with ``{"records": []}``).
* A missing optional layer is logged as a warning containing ``layer_missing:genomad`` /
  ``layer_missing:antismash`` (the run flags themselves are M2-20's).
"""

from __future__ import annotations

import importlib
import json
import logging
import shutil
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow.parquet as pq
import pytest
from helpers.regions_oracle import ANTISMASH_REC, GENOMAD_REC, brute_members

import giapsim
from giap.annotate.bakta import parse_bakta_json
from giap.config import GiapConfig
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

PARSE_OUT = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}
REGION_COLS = {"region_id", "region_type", "contig", "start", "end", "source_tool", "score",
               "genes", "notes"}  # fmt: skip
GENE_COLS = {"region_id", "locus_tag", "overlap_frac", "subset", "core_biosynthetic"}
TAX = "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;"
PROVIRUS = (58_001, 67_503)  # 1-based closed, as geNOMAD writes it
AREA = (20_000, 30_000)  # antiSMASH JSON: 0-based half-open


def _fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"giap.steps.{n}"))
    return reg


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield _fresh_steps(mp)


class Sim:
    def __init__(self, work: Path, genes: list[Gene], crispr: Interval) -> None:
        self.work = work
        self.genes = genes
        self.crispr = crispr
        self.core = {k: work / PARSE_OUT[k] for k in ("genes", "features", "replicons")}


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Sim:
    out = tmp_path_factory.mktemp("sim100k")
    res = giapsim.simulate(giapsim.SIM100K, giapsim.SIM100K_SEED, out / "sim")
    work = out / "work"
    with pytest.MonkeyPatch.context() as mp:
        reg = _fresh_steps(mp)
        reg.run_step(
            "parse_bakta",
            GiapConfig(seed=1),
            {"bakta_json": res.bakta_json, "genome": res.genome_fasta},
            {k: work / v for k, v in PARSE_OUT.items()},
        )
    model = parse_bakta_json(res.bakta_json).model
    (arr,) = [f.interval for f in model.features if f.type == "crispr"]
    return Sim(work, list(model.genes), arr)


def _header(name: str) -> str:
    path = GENOMAD_REC / "sim100k" / "genome_summary" / f"genome_{name}.tsv"
    return path.read_text().splitlines()[0]


def genomad_dir(root: Path, sim: Sim) -> Path:
    """geNOMAD output in the recorded layout: one provirus on chr1, plasmid_1 as a plasmid."""
    d = root / "genomad" / "genome_summary"
    d.mkdir(parents=True)
    a, b = PROVIRUS
    name = f"chr1|provirus_{a}_{b}"
    (d / "genome_virus_summary.tsv").write_text(
        _header("virus_summary") + "\n"
        f"{name}\t{b - a + 1}\tProvirus\t{a}-{b}\t9\t11\t0.9512\tNA\t3\t12.5\t{TAX}\n"
    )
    (d / "genome_plasmid_summary.tsv").write_text(
        _header("plasmid_summary") + "\n"
        "plasmid_1\t8000\tNo terminal repeats\t2\t11\t0.9871\tNA\t1\t3.2\tNA\tNA\n"
    )
    g = next(g for g in sim.genes if g.interval.contig == "chr1" and a <= g.interval.start
             and g.interval.end <= b)  # fmt: skip
    s1, e1 = g.interval.start + 1, g.interval.end
    strand = "1" if g.interval.strand == "+" else "-1"
    (d / "genome_virus_genes.tsv").write_text(
        _header("virus_genes") + "\n"
        f"{name}_1\t{s1}\t{e1}\t{e1 - s1 + 1}\t{strand}\t0.5\t11\tNone\tGENOMAD.1.VV\t1e-20\t90\t0"
        "\t0\t1\t1\tCaudoviricetes\tNA\tNA\tNA\tPhage major capsid protein\n"
    )
    return root / "genomad"


def antismash_dir(root: Path, sim: Sim, name: str = "bakta.json") -> tuple[Path, list[str]]:
    """antiSMASH output dir: ``<name>.json`` with one area on chr1, plus files the step ignores."""
    a, b = AREA
    inside = [g for g in sim.genes if g.interval.contig == "chr1" and a <= g.interval.start
              and g.interval.end <= b]  # fmt: skip
    assert inside, "SIM100K has genes in the test area"
    core = inside[0]
    feats: list[dict[str, Any]] = [
        {"type": "region", "location": f"[{a}:{b}]", "qualifiers": {"product": ["NRPS"]}}
    ]
    for g in inside:
        q: dict[str, Any] = {"locus_tag": [g.locus_tag]}
        if g is core:
            q["gene_kind"] = ["biosynthetic"]
        loc = f"[{g.interval.start}:{g.interval.end}]({g.interval.strand})"
        feats.append({"location": loc, "type": "CDS", "qualifiers": q})
    data = json.loads((ANTISMASH_REC / "sim100k" / "bakta.json").read_text())
    rec = next(r for r in data["records"] if r["id"] == "chr1")
    rec["areas"] = [{"start": a, "end": b, "products": ["NRPS", "T1PKS"]}]
    rec["features"] = feats
    d = root / "antismash"
    d.mkdir(parents=True)
    (d / name).write_text(json.dumps(data))
    (d / "chr1.region001.gbk").write_text("LOCUS chr1\n//\n")
    (d / "bakta.zip").write_bytes(b"PK\x05\x06" + b"\x00" * 18)
    return d, [core.locus_tag]


def _run(reg: ModuleType, inputs: dict[str, Path], out: Path) -> tuple[list[dict], list[dict]]:
    outs = {"regions": out / "regions.parquet", "region_genes": out / "region_genes.parquet"}
    reg.run_step("region_layer", GiapConfig(seed=1), inputs, outs)
    regs = pq.read_table(outs["regions"]).to_pylist()
    rgs = pq.read_table(outs["region_genes"]).to_pylist()
    return regs, rgs


def _genes(cell: Any) -> list[str]:
    if cell is None:
        return []
    return cell.split(";") if isinstance(cell, str) else list(cell)


# ----------------------------------------------------------------------------- registry contract


def test_region_layer_is_registered_with_the_annotate_nf_keys(steps: ModuleType) -> None:
    spec = steps.get("region_layer")
    assert spec.kind == "py"
    assert sorted(spec.inputs) == ["features", "genes", "replicons"]
    assert sorted(spec.optional_inputs) == ["antismash", "genomad"]
    assert sorted(spec.outputs) == ["region_genes", "regions"]
    assert spec.precomputed == "regions"
    assert spec.module == "core"


def test_genomad_and_antismash_have_no_python_step(steps: ModuleType) -> None:
    """M2-18 lead decision 3: parsing happens inside region_layer."""
    py = {s.name for s in steps.all_steps() if s.kind == "py"}
    assert not py & {"genomad", "antismash", "parse_genomad", "parse_antismash"}
    assert "region_layer" in py


# ----------------------------------------------------------------------------- runs


def test_core_only_run_has_bakta_regions_and_flags_missing_layers(
    steps: ModuleType, sim: Sim, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.WARNING):
        regs, rgs = _run(steps, dict(sim.core), tmp_path)
    assert regs and set(regs[0]) >= REGION_COLS
    types = {r["region_type"] for r in regs}
    assert types == {"plasmid", "crispr"}  # SIM100K: plasmid_1 (Bakta kind plasmid), one array
    (p,) = [r for r in regs if r["region_type"] == "plasmid"]
    assert (p["region_id"], p["contig"], p["start"], p["end"]) == (
        "plasmid_1",
        "plasmid_1",
        0,
        8000,
    )
    assert (p["source_tool"], p["score"]) == ("Bakta", None)
    on_plasmid = {g.locus_tag for g in sim.genes if g.interval.contig == "plasmid_1"}
    assert set(_genes(p["genes"])) == on_plasmid
    (c,) = [r for r in regs if r["region_type"] == "crispr"]
    assert (c["region_id"], c["start"], c["end"]) == ("crispr_1", sim.crispr.start, sim.crispr.end)
    assert {x["locus_tag"] for x in rgs if x["region_id"] == "plasmid_1"} == on_plasmid
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "layer_missing:genomad" in msgs and "layer_missing:antismash" in msgs


def test_recorded_sim100k_genomad_and_antismash_outputs_add_nothing(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    """Real sim100k recordings: geNOMAD summaries header-only, antiSMASH JSON without areas."""
    gdir = tmp_path / "in" / "genomad"
    shutil.copytree(GENOMAD_REC / "sim100k", gdir)
    adir = tmp_path / "in" / "antismash"
    adir.mkdir(parents=True)
    shutil.copy(ANTISMASH_REC / "sim100k" / "bakta.json", adir / "bakta.json")
    regs, _ = _run(steps, {**sim.core, "genomad": gdir, "antismash": adir}, tmp_path / "out")
    assert {r["region_type"] for r in regs} == {"plasmid", "crispr"}


def test_nf_stub_antismash_output_is_accepted(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    adir = tmp_path / "antismash"
    adir.mkdir()
    (adir / "sim100k.json").write_text('{"records": []}\n')  # ANTISMASH stub block
    regs, _ = _run(steps, {**sim.core, "antismash": adir}, tmp_path / "out")
    assert not any(r["region_type"] == "bgc" for r in regs)


def test_full_run_builds_every_type(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    gdir = genomad_dir(tmp_path / "in", sim)
    adir, core = antismash_dir(tmp_path / "in", sim)
    regs, rgs = _run(steps, {**sim.core, "genomad": gdir, "antismash": adir}, tmp_path / "out")
    assert set(regs[0]) >= REGION_COLS
    assert set(rgs[0]) >= GENE_COLS
    ids = {r["region_id"]: r for r in regs}
    assert set(ids) == {"prophage_1", "plasmid_1", "bgc_1", "crispr_1"}

    pro = ids["prophage_1"]
    a, b = PROVIRUS
    assert (pro["contig"], pro["start"], pro["end"]) == ("chr1", a - 1, b)
    assert (pro["source_tool"], pro["score"], pro["notes"]) == ("geNOMAD", 0.9512, TAX)
    rep = Replicon("chr1", 100_000, "circular", "chromosome")
    want = brute_members(sim.genes, Interval("chr1", a - 1, b, "."), rep)
    assert set(_genes(pro["genes"])) == set(want)
    pro_rows = {x["locus_tag"]: x for x in rgs if x["region_id"] == "prophage_1"}
    assert set(pro_rows) == set(want)
    for tag, frac in want.items():
        assert pro_rows[tag]["overlap_frac"] == pytest.approx(frac)
        assert pro_rows[tag]["subset"] in {"integrase", "repressor", "structural", "lysis", "cargo"}
    marked = next(g for g in sim.genes if g.interval.contig == "chr1" and a <= g.interval.start
                  and g.interval.end <= b)  # fmt: skip
    assert pro_rows[marked.locus_tag]["subset"] == "structural"  # geNOMAD marker description

    pl = ids["plasmid_1"]
    assert (pl["source_tool"], pl["score"]) == ("Bakta+geNOMAD", 0.9871)

    bg = ids["bgc_1"]
    assert (bg["start"], bg["end"], bg["source_tool"], bg["score"]) == (*AREA, "antiSMASH", None)
    assert bg["notes"] == "products=NRPS,T1PKS"
    bg_rows = {x["locus_tag"]: x for x in rgs if x["region_id"] == "bgc_1"}
    assert {t for t, x in bg_rows.items() if x["core_biosynthetic"]} == set(core)
    assert all(x["subset"] is None for x in bg_rows.values())


def test_rerun_is_byte_identical(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    gdir = genomad_dir(tmp_path / "in", sim)
    adir, _ = antismash_dir(tmp_path / "in", sim)
    inputs = {**sim.core, "genomad": gdir, "antismash": adir}
    _run(steps, inputs, tmp_path / "a")
    _run(steps, inputs, tmp_path / "b")
    for name in ("regions.parquet", "region_genes.parquet"):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes(), name


def test_region_types_config_restricts_the_step(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    from giap.config import resolve_config

    cfg = resolve_config(overrides=["seed=1", "context.region_types=[bgc]"])
    adir, _ = antismash_dir(tmp_path / "in", sim)
    outs = {"regions": tmp_path / "regions.parquet", "region_genes": tmp_path / "rg.parquet"}
    steps.run_step("region_layer", cfg, {**sim.core, "antismash": adir}, outs)
    assert {r["region_type"] for r in pq.read_table(outs["regions"]).to_pylist()} == {"bgc"}
