"""Steps ``uniop_input``, ``uniop`` (ext) and ``operon_pairs`` (M2-12; M2-index step contract;
GIAP_SPEC.md §6.2, §6.4, §6.6, §7.1.3, §9.4).

Runs the registered steps through ``giap.steps.registry.run_step`` on giapsim's SIM100K genome
after M2-04's ``parse_bakta``, the local executor with the simulator's ``precomputed/`` layers
(M2-05), and checks TU assembly against ``truth/tus.tsv``.

Conventions pinned here where the contract table is silent (see the PR description):

* Step keys follow the work-directory file names: ``uniop_input`` → ``uniop_cds`` (GFF) and
  ``uniop_map`` (TSV); ``uniop`` (ext) ``uniop_cds`` → ``uniop_pred``; ``operon_pairs`` →
  ``pairs_raw``, ``tus``, ``tu_upstream`` (directory of ``<tu_id>.fa``). ``uniop_pred`` is an
  optional input; ``operon_pairs`` may also require ``replicons`` and ``genome`` (topology and
  sequence are needed for wrap pairs and upstream FASTA).
* ``uniop_input`` and ``uniop`` are replaced by ``annotation.precomputed.operon_pairs``;
  ``operon_pairs`` itself is not (``precomputed=None``) so TUs are still built from the layer.
* ``pairs_raw.parquet`` has at least ``gene_a, gene_b, uniop_p, source``; ``tus.parquet`` at least
  ``tu_id, contig, start, end, strand, genes`` with 0-based half-open coordinates (work-directory
  tables stay internal, §6.5) and ``genes`` as a native Parquet list column of locus tags in TU
  order (lead decision; ``;`` only in TSV outputs, ADR-0006).
* The distance-model parameters come from ``calibration.yaml`` (d₀ = 50, s = 20).
"""

from __future__ import annotations

import importlib
import itertools
import math
import shutil
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

import giapsim
from giap.annotate.bakta import parse_bakta_json
from giap.annotate.operons import (
    adjacent_pairs,
    build_tus,
    intergenic_distance,
    write_tu_regions,
)
from giap.annotate.uniop import PairProb
from giap.config import GiapConfig, resolve_config
from giap.io.fasta import read_genome_fasta
from giap.io.tsv import read_tsv
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
UNIOP_IN_OUT = {"uniop_cds": "uniop_cds.gff", "uniop_map": "uniop_map.tsv"}
PAIRS_OUT = {"pairs_raw": "pairs_raw.parquet", "tus": "tus.parquet", "tu_upstream": "tu_upstream"}
D0, S = 50.0, 20.0


def logistic(d: float) -> float:
    return 1.0 / (1.0 + math.exp((d - D0) / S))


def _fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    """A fresh registry holding every real step (robust to earlier temporary registries)."""
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":  # reloading the registry itself would drop the patched dict
            importlib.reload(importlib.import_module(f"giap.steps.{n}"))
    return reg


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield _fresh_steps(mp)


class Sim:
    def __init__(self, res: Any, work: Path) -> None:
        self.res = res
        self.work = work
        model = parse_bakta_json(res.bakta_json).model
        self.genes: list[Gene] = list(model.genes)
        self.reps = dict(model.replicons)
        order = [r.id for r in read_genome_fasta(res.genome_fasta)]
        self.reps = {c: self.reps[c] for c in order}
        self.seq = {r.id: r.sequence for r in read_genome_fasta(res.genome_fasta)}
        self.avail: dict[str, Path] = {k: work / v for k, v in PARSE_OUT.items()}
        self.avail["genome"] = res.genome_fasta


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
    return Sim(res, work)


def _inputs(reg: ModuleType, name: str, avail: dict[str, Path]) -> dict[str, Path]:
    spec = reg.get(name)
    missing = [k for k in spec.inputs if k not in avail]
    assert not missing, f"{name}: required inputs {missing} not producible here"
    ins = {k: avail[k] for k in spec.inputs}
    ins |= {k: avail[k] for k in spec.optional_inputs if k in avail}
    return ins


def _run_uniop_input(reg: ModuleType, sim: Sim, out: Path, seed: int = 1) -> dict[str, Path]:
    outs = {k: out / v for k, v in UNIOP_IN_OUT.items()}
    reg.run_step("uniop_input", GiapConfig(seed=seed), _inputs(reg, "uniop_input", sim.avail), outs)
    return outs


def _run_pairs(
    reg: ModuleType,
    sim: Sim,
    out: Path,
    extra: dict[str, Path],
    cfg: GiapConfig | None = None,
) -> dict[str, Path]:
    outs = {k: out / v for k, v in PAIRS_OUT.items()}
    ins = _inputs(reg, "operon_pairs", {**sim.avail, **extra})
    reg.run_step("operon_pairs", cfg or GiapConfig(seed=1), ins, outs)
    return outs


def _gff_ids(gff: Path) -> list[tuple[str, str, str]]:
    """(prodigal ID, contig, strand) per CDS line in file order."""
    out = []
    for ln in gff.read_text().splitlines():
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        out.append((f[8].split(";")[0].removeprefix("ID="), f[0], f[6]))
    return out


def _fake_pred(gff: Path, path: Path, p: float) -> Path:
    ids = _gff_ids(gff)
    rows = [
        f"{a[0]}\t{b[0]}\t{p}" for a, b in itertools.pairwise(ids) if a[1] == b[1] and a[2] == b[2]
    ]
    path.write_text("\n".join(["gene_a\tgene_b\tp", *rows]) + "\n")
    return path


def _pairs(path: Path) -> dict[tuple[str, str], tuple[float, str]]:
    t = pq.read_table(path).to_pylist()
    assert t == [] or {"gene_a", "gene_b", "uniop_p", "source"} <= set(t[0])
    return {(r["gene_a"], r["gene_b"]): (r["uniop_p"], r["source"]) for r in t}


def _genes_cell(v: Any) -> tuple[str, ...]:
    assert isinstance(v, list), f"tus.parquet genes must be a native list column, got {v!r}"
    return tuple(v)


def _tus(path: Path) -> list[tuple[Any, ...]]:
    schema = pq.read_schema(path)
    assert {"tu_id", "contig", "start", "end", "strand", "genes"} <= set(schema.names)
    gt = schema.field("genes").type
    assert pa.types.is_list(gt) or pa.types.is_large_list(gt), f"genes: {gt}"
    assert pa.types.is_string(gt.value_type) or pa.types.is_large_string(gt.value_type)
    return [
        (r["tu_id"], r["contig"], r["start"], r["end"], r["strand"], _genes_cell(r["genes"]))
        for r in pq.read_table(path).to_pylist()
    ]


def _expected_tus(sim: Sim, pairs: list[PairProb], threshold: float = 0.5) -> list[Any]:
    return build_tus(sim.genes, pairs, threshold, sim.reps)


def _tu_rows(tus: list[Any]) -> list[tuple[Any, ...]]:
    return [(t.tu_id, t.contig, t.interval.start, t.interval.end, t.strand, t.genes) for t in tus]


# ---------------------------------------------------------------- declarations (contract table)


def test_step_declarations(steps: ModuleType) -> None:
    ui = steps.get("uniop_input")
    assert (ui.kind, ui.module, ui.precomputed) == ("py", "core", "operon_pairs")
    assert "genes" in ui.inputs and set(ui.inputs) <= {"genes", "replicons"}
    assert ui.outputs == ("uniop_cds", "uniop_map")

    ext = steps.get("uniop")
    assert (ext.kind, ext.module, ext.precomputed, ext.func) == (
        "ext",
        "core",
        "operon_pairs",
        None,
    )
    assert ext.inputs == ("uniop_cds",)
    assert ext.outputs == ("uniop_pred",)

    op = steps.get("operon_pairs")
    assert (op.kind, op.module) == ("py", "core")
    assert op.precomputed is None  # "pairs only; TUs still built"
    assert "uniop_pred" in op.optional_inputs and "uniop_pred" not in op.inputs
    assert {"genes", "uniop_map"} <= set(op.inputs) <= {"genes", "uniop_map", "replicons", "genome"}
    assert op.outputs == ("pairs_raw", "tus", "tu_upstream")


def test_topological_order(steps: ModuleType) -> None:
    order = steps.topo_order(["operon_pairs", "uniop", "uniop_input", "parse_bakta"])
    assert order == ["parse_bakta", "uniop_input", "uniop", "operon_pairs"]


# ---------------------------------------------------------------- uniop_input


def test_uniop_input_writes_gff_and_map_for_every_gene(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    outs = _run_uniop_input(steps, sim, tmp_path)
    header, rows = read_tsv(outs["uniop_map"])
    assert header == ["prodigal_id", "locus_tag"]
    assert sorted(r[1] for r in rows) == sorted(g.locus_tag for g in sim.genes)  # type: ignore[type-var]
    ids = _gff_ids(outs["uniop_cds"])
    assert [i for i, _, _ in ids] == [r[0] for r in rows]
    # k restarts at 1 on each contig, in FASTA order
    contigs = [c for _, c, _ in ids]
    assert contigs == sorted(contigs, key=list(sim.reps).index)
    for pid, contig, _ in ids:
        assert pid.startswith(f"{contig}_")
    for contig in sim.reps:
        ks = [int(i.rsplit("_", 1)[1]) for i, c, _ in ids if c == contig]
        assert ks == list(range(1, len(ks) + 1))


# ---------------------------------------------------------------- operon_pairs with UniOP output


def test_operon_pairs_uses_uniop_and_distance_model_for_wrap(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    ui = _run_uniop_input(steps, sim, tmp_path / "ui")
    pred = _fake_pred(ui["uniop_cds"], tmp_path / "uniop.pred", 0.9)
    outs = _run_pairs(
        steps, sim, tmp_path / "op", {"uniop_map": ui["uniop_map"], "uniop_pred": pred}
    )
    got = _pairs(outs["pairs_raw"])
    adj = adjacent_pairs(sim.genes, sim.reps)
    assert set(got) == {(a.locus_tag, b.locus_tag) for a, b in adj}
    wraps = _wrap_pairs(sim)
    plasmid = next(r.id for r in sim.reps.values() if r.kind == "plasmid")
    assert any(a.interval.contig == plasmid for a, _ in wraps.values())  # two same-strand genes
    for (a, b), (p, src) in got.items():
        if (a, b) in wraps:
            ga, gb = wraps[a, b]
            assert src == "distance_model"
            d = intergenic_distance(ga, gb, sim.reps[ga.interval.contig])
            assert p == pytest.approx(logistic(d))
        else:
            assert (p, src) == (pytest.approx(0.9), "uniop")
    expect_pairs = [PairProb(a, b, p, src) for (a, b), (p, src) in got.items()]  # type: ignore[arg-type]
    assert _tus(outs["tus"]) == _tu_rows(_expected_tus(sim, expect_pairs))


def test_tu_upstream_fastas_match_library(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    ui = _run_uniop_input(steps, sim, tmp_path / "ui")
    pred = _fake_pred(ui["uniop_cds"], tmp_path / "uniop.pred", 0.9)
    outs = _run_pairs(
        steps, sim, tmp_path / "op", {"uniop_map": ui["uniop_map"], "uniop_pred": pred}
    )
    got_pairs = _pairs(outs["pairs_raw"])
    pairs = [PairProb(a, b, p, src) for (a, b), (p, src) in got_pairs.items()]  # type: ignore[arg-type]
    tus = _expected_tus(sim, pairs)
    ref = tmp_path / "ref"
    ref.mkdir()
    write_tu_regions(ref, tus, sim.seq)
    files = sorted(p.name for p in outs["tu_upstream"].iterdir())
    assert files == sorted(f"{t.tu_id}.fa" for t in tus)
    for name in files:
        assert (outs["tu_upstream"] / name).read_bytes() == (ref / name).read_bytes(), name


# ---------------------------------------------------------------- fallback (example 9)


def _wrap_pairs(sim: Sim) -> dict[tuple[str, str], tuple[Gene, Gene]]:
    """(last, first) of every circular contig whose end genes share a strand."""
    out = {}
    for cid, rep in sim.reps.items():
        on = sorted(
            (g for g in sim.genes if g.interval.contig == cid),
            key=lambda g: (g.interval.start, g.interval.end),
        )
        circular = rep.topology == "circular" and len(on) >= 2
        if circular and on[-1].interval.strand == on[0].interval.strand:
            out[on[-1].locus_tag, on[0].locus_tag] = (on[-1], on[0])
    return out


def _check_all_distance_model(sim: Sim, pairs_raw: Path) -> None:
    got = _pairs(pairs_raw)
    adj = adjacent_pairs(sim.genes, sim.reps)
    assert set(got) == {(a.locus_tag, b.locus_tag) for a, b in adj}
    for a, b in adj:
        p, src = got[a.locus_tag, b.locus_tag]
        assert src == "distance_model"
        d = intergenic_distance(a, b, sim.reps[a.interval.contig])
        assert p == pytest.approx(logistic(d)), (a.locus_tag, b.locus_tag, d)


def test_operon_pairs_without_uniop_pred_falls_back(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    ui = _run_uniop_input(steps, sim, tmp_path / "ui")
    outs = _run_pairs(steps, sim, tmp_path / "op", {"uniop_map": ui["uniop_map"]})
    _check_all_distance_model(sim, outs["pairs_raw"])
    assert len(_tus(outs["tus"])) >= 1


@pytest.mark.parametrize("content", ["", "gene_a\tgene_b\tp\n"], ids=["empty", "header_only"])
def test_operon_pairs_with_empty_uniop_pred_falls_back(
    steps: ModuleType, sim: Sim, tmp_path: Path, content: str
) -> None:
    ui = _run_uniop_input(steps, sim, tmp_path / "ui")
    pred = tmp_path / "uniop.pred"
    pred.write_text(content)
    outs = _run_pairs(
        steps, sim, tmp_path / "op", {"uniop_map": ui["uniop_map"], "uniop_pred": pred}
    )
    _check_all_distance_model(sim, outs["pairs_raw"])


def test_operon_pairs_threshold_from_config(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    ui = _run_uniop_input(steps, sim, tmp_path / "ui")
    pred = _fake_pred(ui["uniop_cds"], tmp_path / "uniop.pred", 0.9)
    cfg = resolve_config(None, ["annotation.uniop.pair_threshold=0.95", "seed=1"])
    outs = _run_pairs(
        steps,
        sim,
        tmp_path / "op",
        {"uniop_map": ui["uniop_map"], "uniop_pred": pred},
        cfg,
    )
    got = _pairs(outs["pairs_raw"])
    pairs = [PairProb(a, b, p, src) for (a, b), (p, src) in got.items()]  # type: ignore[arg-type]
    assert _tus(outs["tus"]) == _tu_rows(_expected_tus(sim, pairs, 0.95))
    linked = [k for k, (p, _) in got.items() if p >= 0.95]
    tus = _tus(outs["tus"])
    assert len(tus) == len(sim.genes) - len(linked)


# ---------------------------------------------------------------- determinism (§6.6, TE7)


def test_steps_are_byte_reproducible_and_seed_free(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    def run(tag: str, seed: int) -> Path:
        d = tmp_path / tag
        ui = _run_uniop_input(steps, sim, d / "ui", seed)
        pred = _fake_pred(ui["uniop_cds"], d / "uniop.pred", 0.8)
        _run_pairs(
            steps,
            sim,
            d / "op",
            {"uniop_map": ui["uniop_map"], "uniop_pred": pred},
            GiapConfig(seed=seed),
        )
        return d

    a, b = run("a", 1), run("b", 99)
    fa = sorted(p.relative_to(a) for p in a.rglob("*") if p.is_file())
    fb = sorted(p.relative_to(b) for p in b.rglob("*") if p.is_file())
    assert fa == fb and len(fa) > 5
    for rel in fa:
        assert (a / rel).read_bytes() == (b / rel).read_bytes(), str(rel)


# ---------------------------------------------------------------- precomputed layer (§6.4, M2-02)


def test_precomputed_operon_pairs_layer_still_builds_tus(sim: Sim, tmp_path: Path) -> None:
    from giap.runner.local import run_local

    cfg = resolve_config(sim.res.config_yaml)
    work = tmp_path / "work"
    work.mkdir()
    shutil.copy(sim.res.genome_fasta, work / "genome")
    # until the ext ``bakta`` step is declared, parse_bakta reads ``workdir/bakta_json``
    shutil.copy(sim.res.bakta_json, work / "bakta_json")
    with pytest.MonkeyPatch.context() as mp:
        _fresh_steps(mp)
        run_local(cfg, work, ["operon_pairs"])
    assert not (work / "uniop_cds").exists()  # uniop_input replaced by the layer
    got = _pairs(work / "pairs_raw")
    header, rows = read_tsv(sim.res.precomputed["operon_pairs"])
    assert header == ["gene_a", "gene_b", "uniop_p", "source"]
    for a, b, p, src in rows:
        assert got[a, b] == (pytest.approx(float(p)), src)  # type: ignore[index, arg-type]
    adj = {(a.locus_tag, b.locus_tag) for a, b in adjacent_pairs(sim.genes, sim.reps)}
    assert set(got) == adj
    extra = adj - {(r[0], r[1]) for r in rows}
    assert extra == set(_wrap_pairs(sim))  # the simulator's layer never lists wrap pairs
    assert all(got[k][1] == "distance_model" for k in extra)
    tus = _tus(work / "tus")
    assert sorted(g for t in tus for g in t[5]) == sorted(g.locus_tag for g in sim.genes)
    assert (work / "tu_upstream").is_dir()
    assert len(list((work / "tu_upstream").glob("*.fa"))) == len(tus)


# ---------------------------------------------------------------- simulator truth (M2-05)


def test_oracle_pairs_reproduce_simulator_tus(sim: Sim) -> None:
    """With p = 1 inside a true TU and 0 elsewhere, ``build_tus`` recovers ``truth/tus.tsv``."""
    header, rows = read_tsv(sim.res.truth["tus"])
    truth = [dict(zip(header, r, strict=True)) for r in rows]
    tu_of = {g: i for i, t in enumerate(truth) for g in str(t["genes"]).split(";")}
    pairs = [
        PairProb(
            a.locus_tag,
            b.locus_tag,
            1.0 if tu_of.get(a.locus_tag, -1) == tu_of.get(b.locus_tag, -2) else 0.0,
            "uniop",
        )
        for a, b in adjacent_pairs(sim.genes, sim.reps)
    ]
    tus = build_tus(sim.genes, pairs, 0.5, sim.reps)
    built = {t.genes: t for t in tus}
    assert len(truth) > 10
    for t in truth:
        genes = tuple(str(t["genes"]).split(";"))
        assert genes in built, f"truth TU {t['tu_id']} {genes} not built"
        b = built[genes]
        assert (b.contig, b.strand) == (t["contig"], t["strand"])
        assert (b.interval.start + 1, b.interval.end) == (int(t["start"]), int(t["end"]))  # type: ignore[arg-type]
    # genes outside every true TU (blocks, sORFs, pseudogenes) are singletons
    in_truth = set(tu_of)
    for t in tus:
        if not set(t.genes) & in_truth:
            assert len(t.genes) == 1


def test_simulator_layer_pairs_are_giap_adjacent_pairs(sim: Sim) -> None:
    _, rows = read_tsv(sim.res.precomputed["operon_pairs"])
    adj = adjacent_pairs(sim.genes, sim.reps)
    assert {(r[0], r[1]) for r in rows} <= {(a.locus_tag, b.locus_tag) for a, b in adj}
