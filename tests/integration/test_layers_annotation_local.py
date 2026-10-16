"""M2-20 on SIM100K with every precomputed layer, local executor (docs/specs/M2-20-context-layers.md
examples 1, 5, 7, 8, 9 and the truth-check acceptance criterion; GIAP_SPEC.md §6.4, §7.1.8, §9.1).

``giap.runner.local.run_local(cfg, work, ["context_layers"], roots)`` runs "up to
context_layers" on the simulator's ``config.yaml`` (``bakta_json`` and every M2-02 layer
precomputed), then ``publish`` copies the published keys to ``OUT`` (§9.1 ``OUT/annotation/``).
No tools run (seconds).

Settled by the test author:

* Example 1 "all §7.1.8 files exist": every file except ``terminators_crosscheck.tsv`` (written
  "when TransTermHP ran"; with precomputed terminators it did not) and ``id_map.tsv`` (only with a
  GFF3/GBFF ``--annotation``). ``annotation/tools/`` is absent or empty (no tool ran).
* Truth check: a truth TU (``truth/tus.tsv``) is reproduced exactly when ``annotation/tus.tsv``
  has a TU with the same contig, strand and gene set; every planted region
  (``truth/regions.tsv`` ``true_element = true``) appears in ``annotation/regions.tsv`` with the
  same id, type, contig, start and end.
* Example 5 with precomputed regions: ``--modules core`` drops prophage and BGC ids from
  ``region_ids`` (those modules are off).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from giap.config import GiapConfig, resolve_config
from giap.io.gff3 import read_gff3
from giap.io.tsv import read_tsv
from giap.model.coords import Replicon

FILES = (
    "genes.gff3", "regions.tsv", "regions.bed", "promoters.bed", "terminators.bed",
    "operon_pairs.tsv", "tus.tsv",
)  # fmt: skip
SCHEMA_DIR = Path(__file__).resolve().parents[2] / "src" / "giap" / "schemas" / "annotation"


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim

    return giapsim.simulate(
        giapsim.SIM100K, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim"
    )


def _run(sim: Any, out: Path, *overrides: str) -> tuple[Path, GiapConfig]:
    """Run up to context_layers into ``out/work`` and publish to ``out``."""
    from giap.runner import local
    from giap.steps import load_step_modules

    load_step_modules()
    cfg = resolve_config(sim.config_yaml, overrides)
    roots = {"bakta_json": Path(cfg.annotation.precomputed.bakta_json), "genome": sim.genome_fasta}
    local.run_local(cfg, out / "work", ["context_layers"], roots)
    local.publish(out / "work", out)
    return out / "annotation", cfg


@pytest.fixture(scope="module")
def run_a(sim: Any, tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, GiapConfig]:
    return _run(sim, tmp_path_factory.mktemp("A"))


def _replicons(work: Path) -> dict[str, Replicon]:
    df = pd.read_parquet(work / "replicons")
    return {r.contig: Replicon(r.contig, int(r.length), r.topology, r.kind)
            for r in df.itertuples(index=False)}  # fmt: skip


def _tree(root: Path) -> dict[str, bytes]:
    return {
        str(p.relative_to(root)): p.read_bytes()
        for p in sorted(root.rglob("*"))
        if p.is_file() and "tools" not in p.relative_to(root).parts
    }


def test_context_layers_is_planned_without_tools(sim: Any) -> None:
    from giap.runner import local
    from giap.steps import load_step_modules

    load_step_modules()
    cfg = resolve_config(sim.config_yaml)
    order = local.plan_local(cfg, ["context_layers"], {"bakta_json", "genome"})
    assert order[-1] == "context_layers"
    assert "parse_bakta" in order and "operon_pairs" in order


def test_example_1_every_file_published(run_a: tuple[Path, GiapConfig]) -> None:
    ann, _ = run_a
    for name in FILES:
        assert (ann / name).is_file(), f"OUT/annotation/{name}"
    assert not (ann / "terminators_crosscheck.tsv").exists()
    assert not (ann / "id_map.tsv").exists()
    tools = ann / "tools"
    assert not tools.exists() or not any(tools.iterdir()), "no tool ran"
    assert (ann.parent / "work" / "layers").is_file()


def test_example_7_gff3_round_trip(run_a: tuple[Path, GiapConfig]) -> None:
    ann, _ = run_a
    work = ann.parent / "work"
    reps = _replicons(work)
    feats = read_gff3(ann / "genes.gff3", reps)
    genes = pd.read_parquet(work / "genes")
    genes_f = [f for f in feats if f.source == "Bakta" and f.type in ("CDS", "pseudogene")]
    got = {(f.attributes["ID"], f.interval.contig, f.interval.start, f.interval.end,
            f.interval.strand) for f in genes_f}  # fmt: skip
    want = {(g.locus_tag, g.contig, int(g.start), int(g.end), g.strand)
            for g in genes.itertuples(index=False)}  # fmt: skip
    assert got == want
    text = (ann / "genes.gff3").read_text(encoding="utf-8").splitlines()
    for r in reps.values():
        assert f"##sequence-region {r.id} 1 {r.length}" in text
        assert any(x.startswith(f"{r.id}\t") and "Is_circular=true" in x for x in text) == (
            r.topology == "circular"
        )


def test_example_8_rerun_byte_identical(
    sim: Any, run_a: tuple[Path, GiapConfig], tmp_path: Path
) -> None:
    ann_b, _ = _run(sim, tmp_path / "B")
    ta, tb = _tree(run_a[0]), _tree(ann_b)
    assert set(ta) == set(tb)
    for name in ta:
        assert ta[name] == tb[name], name
    lay_a = (run_a[0].parent / "work" / "layers").read_bytes()
    assert lay_a == (tmp_path / "B" / "work" / "layers").read_bytes()


@pytest.mark.parametrize("table", ["regions", "operon_pairs", "tus"])
def test_example_9_tsvs_validate(run_a: tuple[Path, GiapConfig], table: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    found = [p for p in (SCHEMA_DIR / f"{table}.json", SCHEMA_DIR / f"{table}.schema.json")
             if p.exists()]  # fmt: skip
    assert len(found) == 1, found
    d = {k: v for k, v in json.loads(found[0].read_text()).items() if k != "coordinates"}
    res = frictionless.Resource(
        path=f"{table}.tsv", basepath=str(run_a[0]), schema=frictionless.Schema.from_descriptor(d)
    )
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_truth_tus_reproduced_from_precomputed_pairs(
    sim: Any, run_a: tuple[Path, GiapConfig]
) -> None:
    head, rows = read_tsv(sim.truth["tus"])
    ix = {c: head.index(c) for c in ("contig", "strand", "genes")}
    truth = [(r[ix["contig"]], r[ix["strand"]], frozenset(str(r[ix["genes"]]).split(";")))
             for r in rows if r[ix["genes"]]]  # fmt: skip
    h2, got_rows = read_tsv(run_a[0] / "tus.tsv")
    j = {c: h2.index(c) for c in ("contig", "strand", "genes")}
    got = {(r[j["contig"]], r[j["strand"]], frozenset(str(r[j["genes"]]).split(";")))
           for r in got_rows}  # fmt: skip
    hit = sum(t in got for t in truth)
    assert truth and hit / len(truth) >= 0.95, f"{hit}/{len(truth)} truth TUs reproduced"


def test_every_planted_region_in_regions_tsv(sim: Any, run_a: tuple[Path, GiapConfig]) -> None:
    head, rows = read_tsv(sim.truth["regions"])
    cols = ("region_id", "region_type", "contig", "start", "end")
    t = head.index("true_element")
    planted = {tuple(r[head.index(c)] for c in cols) for r in rows if r[t] == "true"}
    assert planted
    h2, got = read_tsv(run_a[0] / "regions.tsv")
    have = {tuple(r[h2.index(c)] for c in cols) for r in got}
    assert planted <= have, sorted(planted - have)


def test_layers_cover_every_gene_with_precomputed_values(run_a: tuple[Path, GiapConfig]) -> None:
    work = run_a[0].parent / "work"
    lay = pd.read_parquet(work / "layers")
    genes = pd.read_parquet(work / "genes")
    assert lay["locus_tag"].tolist() == genes["locus_tag"].tolist()
    for c in ("promoter_score", "tx_rate_max_tu", "terminator_conf_downstream", "tir_ostir", "gc3",
              "uniop_p_next"):  # fmt: skip
        assert lay[c].notna().any(), c
    ids = {x for cell in lay["region_ids"] for x in cell}
    assert any(x.startswith("prophage_") for x in ids) and any(x.startswith("bgc_") for x in ids)


def test_example_5_core_only_drops_prophage_and_bgc_ids(sim: Any, tmp_path: Path) -> None:
    ann, _ = _run(sim, tmp_path / "CORE", "modules=[core]")
    lay = pd.read_parquet(ann.parent / "work" / "layers")
    ids = {x for cell in lay["region_ids"] for x in cell}
    assert ids, "plasmid / CRISPR regions remain"
    assert not {x for x in ids if x.startswith(("prophage_", "bgc_"))}
