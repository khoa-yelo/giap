"""Steps ``ostir_input``, ``ostir`` (ext) and ``rbs_layer`` (M2-16; M2-index step contract;
GIAP_SPEC.md §6.2, §6.5, §7.1.6, §9.3).

Runs the registered steps through ``run_step`` on giapsim's SIM100K after M2-04's ``parse_bakta``.

Contract pinned here (test author; the implementer keeps ``workflow/tests/data/annotate/
step_registry.json`` and ``annotate.nf`` consistent with it):

* ``ostir_input``: inputs ``genes``, ``replicons``, ``genome``, ``annotation_meta`` (the M2-index
  row lists only ``genes.parquet`` and ``annotation_meta.json``, but the -50/+50 windows need the
  genome sequence and the replicon topology, like ``operon_pairs``); output ``ostir_input``
  (``ostir_input.csv``); ``precomputed="rbs"``. ``annotate.nf`` must pass ``reps`` and
  ``kv('genome', genome)`` to ``stepIn('ostir_input', ...)``.
* ``ostir`` (ext): ``ostir_input`` → ``ostir``; ``precomputed="rbs"``.
* ``rbs_layer``: inputs ``ostir``, ``ostir_input``; output ``rbs`` (``rbs.parquet``);
  ``precomputed="rbs"``. Columns exactly the precomputed ``rbs`` layer's (M2-index decision 1):
  ``locus_tag, tir_ostir, sd_dg, mrna_dg, spacing_dg, standby_dg, spacer_bp, start_codon``
  (float64 x 5, int64, string). One row per gene that has an OSTIR row at its requested start,
  in ``ostir_input.csv`` order; genes without one are absent (→ ``NA`` downstream, example 8).
  ``start_codon`` is the gene's (``seq[start-1:start+2]`` of its input row), not OSTIR's column.
* ``ostir_input`` writes one row per gene of ``genes.parquet`` (CDS, sORF, pseudogene) in that
  file's §6.5 order, anti-SD from ``annotation_meta.json`` ``anti_sd``, windows from
  ``annotation.ostir.upstream_nt`` / ``downstream_nt``; byte-identical on rerun.
"""

from __future__ import annotations

import csv
import importlib
import json
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

import giapsim
from giap.annotate.precomputed import PRECOMPUTED_COLUMNS
from giap.config import GiapConfig
from giap.io.fasta import read_genome_fasta
from giap.model.coords import Interval, Replicon

PARSE_OUT = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}
OSTIR_HEADER = (
    "name,start_codon,start_position,expression,RBS_distance_bp,dG_total,dG_rRNA:mRNA,dG_mRNA,"
    "dG_spacing,dG_standby,dG_start_codon"
)


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
    def __init__(self, res: Any, work: Path) -> None:
        self.res = res
        self.work = work
        self.avail: dict[str, Path] = {k: work / v for k, v in PARSE_OUT.items()}
        self.avail["genome"] = res.genome_fasta
        self.genes = pq.read_table(work / "genes.parquet").to_pylist()
        self.reps = {
            r["contig"]: Replicon(r["contig"], r["length"], r["topology"], r["kind"])
            for r in pq.read_table(work / "replicons.parquet").to_pylist()
        }
        self.seq = {r.id: r.sequence for r in read_genome_fasta(res.genome_fasta)}
        self.meta = json.loads((work / "annotation_meta.json").read_text())


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


def _run_input(
    reg: ModuleType, sim: Sim, out: Path, cfg: GiapConfig | None = None, **over: Path
) -> Path:
    spec = reg.get("ostir_input")
    avail = {**sim.avail, **over}
    ins = {k: avail[k] for k in spec.inputs}
    target = out / "ostir_input.csv"
    reg.run_step("ostir_input", cfg or GiapConfig(seed=1), ins, {"ostir_input": target})
    return target


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- declarations


def test_step_declarations(steps: ModuleType) -> None:
    oi = steps.get("ostir_input")
    assert (oi.kind, oi.module, oi.precomputed) == ("py", "core", "rbs")
    assert set(oi.inputs) == {"genes", "replicons", "genome", "annotation_meta"}
    assert oi.optional_inputs == ()
    assert oi.outputs == ("ostir_input",)

    ext = steps.get("ostir")
    assert (ext.kind, ext.module, ext.precomputed, ext.func) == ("ext", "core", "rbs", None)
    assert ext.inputs == ("ostir_input",)
    assert ext.outputs == ("ostir",)

    rl = steps.get("rbs_layer")
    assert (rl.kind, rl.module, rl.precomputed) == ("py", "core", "rbs")
    assert set(rl.inputs) == {"ostir", "ostir_input"}
    assert rl.optional_inputs == ()
    assert rl.outputs == ("rbs",)


def test_topological_order(steps: ModuleType) -> None:
    order = steps.topo_order(["rbs_layer", "ostir", "ostir_input", "parse_bakta"])
    assert order == ["parse_bakta", "ostir_input", "ostir", "rbs_layer"]


# ---------------------------------------------------------------- ostir_input


def test_ostir_input_one_row_per_gene_in_gene_order(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    from giap.annotate.ostir import ostir_row, write_ostir_input
    from giap.model.genes import Gene

    path = _run_input(steps, sim, tmp_path)
    rows = _rows(path)
    assert path.read_text().splitlines()[0] == "id,seq,start,end,anti-Shine-Dalgarno"
    assert [r["id"] for r in rows] == [g["locus_tag"] for g in sim.genes]  # §6.5 order, all types
    assert {g["feature_type"] for g in sim.genes} >= {"cds"}

    expect = []
    for g in sim.genes:
        gene = Gene(
            locus_tag=g["locus_tag"],
            feature_type=g["feature_type"],
            interval=Interval(g["contig"], g["start"], g["end"], g["strand"]),
            aa=g["aa"] or "",
            nt=g["nt"],
            start_codon=g["start_codon"],
        )
        rep = sim.reps[g["contig"]]
        expect.append(ostir_row(gene, rep, sim.seq[g["contig"]], sim.meta["anti_sd"]))
    ref = tmp_path / "ref.csv"
    write_ostir_input(ref, expect)
    assert path.read_bytes() == ref.read_bytes()


def test_ostir_input_start_codon_identity(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    rows = _rows(_run_input(steps, sim, tmp_path))
    by_id = {g["locus_tag"]: g for g in sim.genes}
    circular = {c for c, r in sim.reps.items() if r.topology == "circular"}
    for r in rows:
        g = by_id[r["id"]]
        s = int(r["start"])
        assert r["start"] == r["end"]
        if g["start_codon"]:
            assert r["seq"][s - 1 : s + 2] == g["start_codon"], r["id"]
        if g["contig"] in circular:
            assert (len(r["seq"]), s) == (100, 51), r["id"]
        assert r["anti-Shine-Dalgarno"] == sim.meta["anti_sd"]


def test_ostir_input_covers_plasmid_sorf_and_pseudogenes(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    """SIM100K has a plasmid, sORFs and pseudogenes (M3-18); every one of them gets a row."""
    rows = {r["id"]: r for r in _rows(_run_input(steps, sim, tmp_path))}
    plasmid = [g for g in sim.genes if sim.reps[g["contig"]].kind == "plasmid"]
    other = [g for g in sim.genes if g["feature_type"] in ("sorf", "pseudogene")]
    assert plasmid and {g["feature_type"] for g in other} == {"sorf", "pseudogene"}
    for g in plasmid + other:
        assert g["locus_tag"] in rows


def test_ostir_input_byte_identical_on_rerun(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    a = _run_input(steps, sim, tmp_path / "a")
    b = _run_input(steps, sim, tmp_path / "b", GiapConfig(seed=99))
    assert a.read_bytes() == b.read_bytes()


def test_ostir_input_uses_the_anti_sd_from_annotation_meta(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    meta = dict(sim.meta, anti_sd="TCCTCCACT")
    mp = tmp_path / "annotation_meta.json"
    mp.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    rows = _rows(_run_input(steps, sim, tmp_path, annotation_meta=mp))
    assert rows and {r["anti-Shine-Dalgarno"] for r in rows} == {"TCCTCCACT"}


def test_ostir_input_honours_window_config(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    cfg = GiapConfig.model_validate(
        {"seed": 1, "annotation": {"ostir": {"upstream_nt": 30, "downstream_nt": 40}}}
    )
    rows = _rows(_run_input(steps, sim, tmp_path, cfg))
    by_id = {g["locus_tag"]: g for g in sim.genes}
    for r in rows:
        if sim.reps[by_id[r["id"]]["contig"]].topology == "circular":
            assert (len(r["seq"]), r["start"], r["end"]) == (70, "31", "31")


# ---------------------------------------------------------------- rbs_layer


def _fake_ostir(inp: Path, out: Path, drop: str, wrong_codon: str) -> list[dict[str, str]]:
    """What OSTIR 1.1.3 writes for ``inp`` (header and column order from the recorded runs):
    one row per id at its requested start, except ``drop``; an extra candidate row at another
    position for the first id; OSTIR's ``start_codon`` column wrong for ``wrong_codon``."""
    rows = _rows(inp)
    lines = [OSTIR_HEADER]
    for k, r in enumerate(rows):
        if r["id"] == drop:
            continue
        s = int(r["start"])
        codon = "CTG" if r["id"] == wrong_codon else r["seq"][s - 1 : s + 2]
        if k == 0:
            lines.append(f"{r['id']},GTG,{s + 12},1.5,4,9.0,-1.0,-2.0,0.5,0.0,-0.07")
        lines.append(
            f"{r['id']},{codon},{s},{100.0 + k},{k % 13},{-1.0 * k},{-0.5 * k},{-0.25 * k},"
            f"{0.125 * k},{-0.0625 * (k % 3)},-1.194"
        )
    out.write_text("\n".join(lines) + "\n")
    return rows


def _run_rbs(reg: ModuleType, inp: Path, ostir: Path, out: Path) -> Path:
    target = out / "rbs.parquet"
    reg.run_step(
        "rbs_layer", GiapConfig(seed=1), {"ostir": ostir, "ostir_input": inp}, {"rbs": target}
    )
    return target


def test_rbs_layer_features_and_schema(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    inp = _run_input(steps, sim, tmp_path)
    rows = _rows(inp)
    drop, wrong = rows[3]["id"], rows[5]["id"]
    _fake_ostir(inp, tmp_path / "ostir.csv", drop, wrong)
    rbs = _run_rbs(steps, inp, tmp_path / "ostir.csv", tmp_path)

    table = pq.read_table(rbs)
    assert tuple(table.column_names) == PRECOMPUTED_COLUMNS["rbs"]
    types = dict(zip(table.column_names, table.schema.types, strict=True))
    assert types["locus_tag"] == pa.string() and types["start_codon"] == pa.string()
    for c in ("tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg"):
        assert types[c] == pa.float64(), c
    assert types["spacer_bp"] == pa.int64()

    got = table.to_pylist()
    kept = [(k, r) for k, r in enumerate(rows) if r["id"] != drop]
    assert [g["locus_tag"] for g in got] == [r["id"] for _, r in kept]  # input order; drop absent
    for (k, r), g in zip(kept, got, strict=True):
        s = int(r["start"])
        assert g == {
            "locus_tag": r["id"],
            "tir_ostir": 100.0 + k,
            "sd_dg": -0.5 * k,
            "mrna_dg": -0.25 * k,
            "spacing_dg": 0.125 * k,
            "standby_dg": -0.0625 * (k % 3),
            "spacer_bp": k % 13,
            "start_codon": r["seq"][s - 1 : s + 2],  # the gene's codon, also for `wrong`
        }
    assert next(g for g in got if g["locus_tag"] == wrong)["start_codon"] != "CTG"


def test_rbs_layer_byte_identical_on_rerun(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    inp = _run_input(steps, sim, tmp_path)
    rows = _rows(inp)
    _fake_ostir(inp, tmp_path / "ostir.csv", rows[0]["id"], rows[1]["id"])
    a = _run_rbs(steps, inp, tmp_path / "ostir.csv", tmp_path / "a")
    b = _run_rbs(steps, inp, tmp_path / "ostir.csv", tmp_path / "b")
    assert a.read_bytes() == b.read_bytes()


def test_rbs_layer_header_only_ostir_output(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    """The OSTIR stub (``-stub-run``) writes only the header: an empty layer, not a crash."""
    inp = _run_input(steps, sim, tmp_path)
    ostir = tmp_path / "ostir.csv"
    ostir.write_text(OSTIR_HEADER + "\n")
    table = pq.read_table(_run_rbs(steps, inp, ostir, tmp_path))
    assert table.num_rows == 0
    assert tuple(table.column_names) == PRECOMPUTED_COLUMNS["rbs"]
