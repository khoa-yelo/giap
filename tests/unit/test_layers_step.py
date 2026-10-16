"""M2-20 step ``context_layers``: registry contract, ``annotation/*`` files, ``layers.parquet``,
determinism, ``layer_missing`` logging, native copies, schemas (docs/specs/M2-20-context-layers.md;
GIAP_SPEC.md §6.2, §6.3, §6.6/R7, §7.1.8, §9.1-§9.5).

Inputs: the merged producers run on giapsim SIM100K (``helpers.layers_chain``) and the
simulator's precomputed layers (M2-02 TSVs).

Pinned / settled by the test author (see the spec questions and annotate.nf):

* Registry keys follow ``annotate.nf``'s ``CONTEXT_LAYERS`` call: outputs ``layers``
  (``layers.parquet``) and ``annotation`` (the ``annotation/`` directory) — the spec's
  ``annotation_dir`` is called ``annotation`` there. Required inputs ``genes``, ``replicons``,
  ``pairs_raw``, ``tus`` (spec); optional the spec's list plus the two extra keys annotate.nf
  stages: ``features`` (the ``core`` list) and ``promoters`` (``pre('promoters')`` when
  ``annotation.precomputed.promoters`` is set). ``precomputed`` None, module ``core``.
* With a precomputed layer the key holds the M2-02 TSV, not a parquet file: under the local
  executor every output key of a replaced step points to the layer file (``promoter_windows``,
  ``tx_rates``, ``terminators`` and ``terminators_crosscheck``, ``rbs``, ``regions`` and
  ``region_genes``); under Nextflow the promoter TSV arrives as ``promoters``. The step reads it
  when ``annotation.precomputed.<layer>`` is set.
* ``annotation/terminators_crosscheck.tsv`` is a byte copy of the ``terminator_layer`` table and is
  written only when that table was produced (TransTermHP ran), never from a precomputed
  terminators layer. ``annotation/id_map.tsv`` is a byte copy of the M2-11 ``id_map`` input and
  only exists with it.
* The §9.3 TSVs and the three BEDs equal, byte for byte, the merged writers (M2-12
  ``write_operon_pairs_tsv`` / ``write_tus_tsv``, M2-19 ``write_regions_tsv`` / ``regions_bed``,
  M2-13 ``promoters_bed`` with ``promotech_yes``, M2-15 ``terminators_bed``) applied to the
  producer files.
* ``layer_missing:<layer>`` run flags are logged as warnings by the step (the M2-19
  ``region_layer`` convention; ``giap.steps`` events reach ``logs/events.jsonl``).
* Schemas: ``src/giap/schemas/annotation/<table>.json`` (``<table>.schema.json`` accepted).
* Lead decision 4: ``gc_local_dev`` needs the sequence, so ``genome`` (the genome FASTA, the key
  annotate.nf already uses for ``parse_bakta`` / ``ostir_input``) is a required input and
  annotate.nf stages ``kv('genome', genome)`` for ``context_layers``; layers.parquet holds a float
  ``gc_local_dev`` column.
* Lead decision 5: annotate.nf passes ``terminators_crosscheck`` (the ``terminator_layer``
  output ``terminators_crosscheck.tsv``) to ``CONTEXT_LAYERS``.
"""

from __future__ import annotations

import json
import logging
import re
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from helpers import layers_oracle as orc
from helpers.layers_chain import (
    ANNOTATION_FILES,
    REC,
    Chain,
    build_chain,
    fresh_steps,
    genome_seqs,
)

from giap.config import GiapConfig

REPO = Path(__file__).resolve().parents[2]
SUBWF = REPO / "workflow" / "subworkflows" / "local" / "annotate.nf"
FIXTURE = REPO / "workflow" / "tests" / "data" / "annotate" / "step_registry.json"
SCHEMA_DIR = REPO / "src" / "giap" / "schemas" / "annotation"
REQUIRED = {"genes", "replicons", "pairs_raw", "tus", "genome"}
OPTIONAL = {
    "features", "promoter_windows", "promoters", "tx_rates", "terminators",
    "terminators_crosscheck", "rbs", "codon", "regions", "region_genes", "id_map", "tool_dirs",
}  # fmt: skip
SCHEMA_TABLES = {
    "regions": ("region_id", "region_type", "contig", "start", "end", "source_tool", "score",
                "genes", "notes"),
    "operon_pairs": ("gene_a", "gene_b", "uniop_p", "source"),
    "tus": ("tu_id", "contig", "start", "end", "strand", "genes", "upstream_start",
            "upstream_end"),
    "id_map": ("source_id", "locus_tag", "contig", "start", "end", "strand", "match"),
    "terminators_crosscheck": ("contig", "n_terminatornet", "n_transtermhp", "n_matched",
                               "share_matched"),
}  # fmt: skip
TIMESTAMP = re.compile(rb"20\d\d-[01]\d-[0-3]\d[T ][0-2]\d:[0-5]\d")


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(tmp_path_factory.mktemp("layers_step_chain"))


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield fresh_steps(mp)


def run(out: Path, inputs: dict[str, Any], cfg: GiapConfig | None = None) -> tuple[Path, Path]:
    """Run ``context_layers`` into ``out``; returns (annotation dir, layers.parquet)."""
    outs = {"annotation": out / "annotation", "layers": out / "layers.parquet"}
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        reg.run_step("context_layers", cfg or GiapConfig(seed=1), inputs, outs)
    return outs["annotation"], outs["layers"]


@pytest.fixture(scope="module")
def full(chain: Chain, tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    return run(tmp_path_factory.mktemp("layers_full"), dict(chain.inputs))


def _schema(table: str) -> dict[str, Any]:
    found = [p for p in (SCHEMA_DIR / f"{table}.json", SCHEMA_DIR / f"{table}.schema.json")
             if p.exists()]  # fmt: skip
    assert len(found) == 1, f"expected one {table} schema in {SCHEMA_DIR}, got {found}"
    return json.loads(found[0].read_text(encoding="utf-8"))


def _frictionless_valid(path: Path, table: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in _schema(table).items() if k != "coordinates"}
    res = frictionless.Resource(
        path=path.name, basepath=str(path.parent), schema=frictionless.Schema.from_descriptor(d)
    )
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def _tree(root: Path) -> dict[str, bytes]:
    return {
        str(p.relative_to(root)): p.read_bytes()
        for p in sorted(root.rglob("*"))
        if p.is_file() and "tools" not in p.relative_to(root).parts
    }


# ----------------------------------------------------------------------------- registry


def test_registered_with_the_annotate_nf_keys(steps: ModuleType) -> None:
    spec = steps.get("context_layers")
    assert spec.kind == "py"
    assert set(spec.inputs) == REQUIRED
    assert set(spec.optional_inputs) == OPTIONAL
    assert set(spec.outputs) == {"annotation", "layers"}
    assert spec.precomputed is None
    assert spec.module == "core"


def test_annotate_nf_out_keys_and_staged_keys_are_declared(steps: ModuleType) -> None:
    text = re.sub(r"(?m)//.*$", "", SUBWF.read_text(encoding="utf-8"))
    m = re.search(r"CONTEXT_LAYERS\(\s*cl\.files\s*,\s*cl\.keys\s*,\s*\[([^\]]*)\]", text)
    assert m, "annotate.nf calls CONTEXT_LAYERS(cl.files, cl.keys, [...])"
    assert set(re.findall(r"'(\w+)=", m.group(1))) == set(steps.get("context_layers").outputs)
    staged = set(re.findall(r"kv\(\s*'(\w+)'", text)) | set(re.findall(r"pre\(\s*'(\w+)'", text))
    core = {"genes", "features", "replicons"}
    layer_keys = staged & {"pairs_raw", "tus", "promoter_windows", "promoters", "tx_rates",
                           "terminators", "rbs", "codon", "regions", "region_genes"}  # fmt: skip
    spec = steps.get("context_layers")
    assert core | layer_keys <= set(spec.inputs) | set(spec.optional_inputs)


def test_step_registry_fixture_has_context_layers() -> None:
    reg = json.loads(FIXTURE.read_text(encoding="utf-8"))
    assert reg["context_layers"] == {
        "inputs": sorted(REQUIRED), "optional_inputs": sorted(OPTIONAL),
        "outputs": ["annotation", "layers"],
    }  # fmt: skip


# ----------------------------------------------------------------------------- gc_local_dev


def test_layers_parquet_has_float_gc_local_dev(chain: Chain, full: tuple[Path, Path]) -> None:
    """Lead decision 4: the column is in layers.parquet, float, and equals the window oracle on
    the ``genome`` input."""
    df = pd.read_parquet(full[1])
    assert "gc_local_dev" in df.columns
    assert pd.api.types.is_float_dtype(df["gc_local_dev"]), df["gc_local_dev"].dtype
    want = orc.gc_local_dev(chain.frame("genes"), chain.frame("replicons"), genome_seqs(chain))
    for tag, g, w in zip(df["locus_tag"], df["gc_local_dev"], want, strict=True):
        assert orc.same(float(g), w) or abs(float(g) - w) < 1e-12, (tag, g, w)
    assert df["gc_local_dev"].notna().all()


# ----------------------------------------------------------------------------- annotate.nf wiring


def _context_layers_staged_keys() -> set[str]:
    """Keys staged into ``stepIn('context_layers', <expr>)``: literal ``kv('k', ...)`` /
    ``pre('k')`` in ``<expr>``, plus those added to the ``layers`` / ``core`` lists when
    ``<expr>`` uses them."""
    text = re.sub(r"(?m)//.*$", "", SUBWF.read_text(encoding="utf-8"))
    m = re.search(r"def\s+(\w+)\s*=\s*stepIn\(\s*'context_layers'\s*,", text)
    assert m, "annotate.nf builds the CONTEXT_LAYERS inputs with stepIn('context_layers', ...)"
    depth, i = 1, m.end()
    while depth:
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    expr = text[m.end() : i - 1]

    var_key = dict(re.findall(r"def\s+(\w+)\s*=\s*kv\(\s*'(\w+)'", text))

    def keys(s: str) -> set[str]:
        lit = set(re.findall(r"\bkv\(\s*'(\w+)'", s)) | set(re.findall(r"\bpre\(\s*'(\w+)'", s))
        return lit | {var_key[w] for w in re.findall(r"\b(\w+)\b", s) if w in var_key}

    out = keys(expr)
    if re.search(r"\blayers\b", expr):
        for blk in re.findall(r"\blayers\s*\+=\s*\[(.*?)\]\s*$", text[: m.start()], re.M | re.S):
            out |= keys(blk)
    if re.search(r"\bcore\b", expr):
        core = re.search(r"def\s+core\s*=\s*\[(.*?)\]\s*$", text, re.M)
        assert core, "annotate.nf defines core = [...]"
        out |= keys(core.group(1))
    return out


def test_annotate_nf_passes_terminators_crosscheck_to_context_layers() -> None:
    """Lead decision 5: the terminator_layer cross-check table reaches CONTEXT_LAYERS."""
    assert "terminators_crosscheck" in _context_layers_staged_keys()
    text = re.sub(r"(?m)//.*$", "", SUBWF.read_text(encoding="utf-8"))
    assert re.search(
        r"kv\(\s*'terminators_crosscheck'\s*,\s*pick\(\s*TERMINATOR_LAYER\.out\.outputs\s*,"
        r"\s*'terminators_crosscheck\.tsv'\s*\)\s*\)",
        text,
    ), "terminators_crosscheck is picked from TERMINATOR_LAYER's terminators_crosscheck.tsv"


def test_annotate_nf_passes_genome_to_context_layers() -> None:
    """Lead decision 4: gc_local_dev needs the genome FASTA; every required key is staged."""
    staged = _context_layers_staged_keys()
    assert "genome" in staged
    assert staged >= REQUIRED, sorted(REQUIRED - staged)


# ----------------------------------------------------------------------------- outputs


def test_every_annotation_file_and_layers_parquet_written(full: tuple[Path, Path]) -> None:
    ann, lay = full
    for name in ANNOTATION_FILES:
        assert (ann / name).is_file(), name
    assert not (ann / "id_map.tsv").exists(), "id_map.tsv only with a GFF3/GBFF --annotation"
    assert lay.is_file()
    from giap.annotate.layers import LAYER_COLUMNS

    assert tuple(pd.read_parquet(lay).columns) == tuple(LAYER_COLUMNS)


def test_layers_parquet_equals_gene_layers(chain: Chain, full: tuple[Path, Path]) -> None:
    from helpers.layers_chain import call_gene_layers, read_layers_frames

    got = pd.read_parquet(full[1])
    want = call_gene_layers(read_layers_frames(chain), genome_seqs(chain))
    assert got["locus_tag"].tolist() == want["locus_tag"].tolist()
    for c in want.columns:
        for a, b in zip(got[c].tolist(), want[c].tolist(), strict=True):
            if c in ("region_ids", "flags"):
                assert list(a) == list(b), c
            else:
                assert orc.same(a, b), (c, a, b)


@pytest.mark.parametrize(
    "name", ["operon_pairs.tsv", "tus.tsv", "regions.tsv", "regions.bed", "promoters.bed",
             "terminators.bed"],
)  # fmt: skip
def test_tables_equal_the_merged_writers(
    chain: Chain, full: tuple[Path, Path], tmp_path: Path, name: str
) -> None:
    want = orc.expected_annotation_files(
        chain.inputs, tmp_path, GiapConfig().annotation.promoters.promotech_yes
    )
    assert (full[0] / name).read_bytes() == want[name].read_bytes(), name


@pytest.mark.parametrize("table", ["regions", "operon_pairs", "tus"])
def test_section_9_3_headers(full: tuple[Path, Path], table: str) -> None:
    head = (full[0] / f"{table}.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert tuple(head.split("\t")) == SCHEMA_TABLES[table]


def test_terminators_crosscheck_is_the_producer_table(
    chain: Chain, full: tuple[Path, Path]
) -> None:
    got = (full[0] / "terminators_crosscheck.tsv").read_bytes()
    assert got == chain.inputs["terminators_crosscheck"].read_bytes()


def test_no_crosscheck_input_no_crosscheck_file(chain: Chain, tmp_path: Path) -> None:
    ann, _ = run(tmp_path, chain.subset("terminators_crosscheck"))
    assert not (ann / "terminators_crosscheck.tsv").exists()


def test_id_map_copied_when_given(chain: Chain, tmp_path: Path) -> None:
    """M2-11 producer: simulator GFF3 -> read_user_cds -> build_id_map -> write_id_map."""
    from giap.annotate.user_annotation import build_id_map, read_user_cds, write_id_map

    reps = orc.replicons(chain.frame("replicons"))
    gff = chain.sim.annotation_gff3 or chain.sim.gff3
    assert gff is not None, "giapsim writes an annotation GFF3"
    cds = read_user_cds(Path(gff), "gff3", reps, {c: c for c in reps})
    rows = build_id_map(cds, orc.genes(chain.frame("genes")), list(reps))
    idm = tmp_path / "id_map.tsv"
    write_id_map(idm, rows)
    ann, _ = run(tmp_path / "out", {**chain.inputs, "id_map": idm})
    assert (ann / "id_map.tsv").read_bytes() == idm.read_bytes()
    _frictionless_valid(ann / "id_map.tsv", "id_map")


# ----------------------------------------------------------------------------- determinism


def test_example_8_rerun_is_byte_identical(chain: Chain, tmp_path: Path) -> None:
    """TE7 / R7: byte-identical annotation/* (tools/ excluded) and layers.parquet, also when
    written under different absolute paths."""
    a_ann, a_lay = run(tmp_path / "first_run", dict(chain.inputs))
    b_ann, b_lay = run(tmp_path / "second" / "deeper_run", dict(chain.inputs))
    ta, tb = _tree(a_ann), _tree(b_ann)
    assert set(ta) == set(tb)
    for name in ta:
        assert ta[name] == tb[name], name
    assert a_lay.read_bytes() == b_lay.read_bytes()


def test_no_timestamps_or_absolute_paths(chain: Chain, full: tuple[Path, Path]) -> None:
    ann, lay = full
    needles = {str(chain.root).encode(), str(ann.parent).encode(), str(Path.home()).encode()}
    for name, data in {**_tree(ann), "layers.parquet": lay.read_bytes()}.items():
        assert not TIMESTAMP.search(data), f"{name} holds a timestamp"
        for n in needles:
            assert n not in data, f"{name} holds the absolute path {n!r}"


# ----------------------------------------------------------------------------- missing layers


@pytest.mark.parametrize(
    ("cfg_data", "drop", "flag"),
    [
        ({"annotation": {"promoter_calculator": {"enabled": False}}}, ("tx_rates",),
         "layer_missing:tx_rates"),
        ({"modules": ["core"]}, (), "layer_missing:genomad"),
        ({"modules": ["core"]}, (), "layer_missing:antismash"),
    ],
)  # fmt: skip
def test_missing_optional_layer_is_logged(
    chain: Chain, tmp_path: Path, caplog: pytest.LogCaptureFixture,
    cfg_data: dict[str, Any], drop: tuple[str, ...], flag: str,
) -> None:  # fmt: skip
    cfg = GiapConfig.model_validate({"seed": 1, **cfg_data})
    with caplog.at_level(logging.WARNING, logger="giap.steps"):
        _, lay = run(tmp_path, chain.subset(*drop), cfg)
    assert any(flag in r.getMessage() for r in caplog.records), [
        r.getMessage() for r in caplog.records
    ]
    if "tx_rates" in drop:
        df = pd.read_parquet(lay)
        assert df[["tx_rate_max_tu", "tx_rate_sum_tu", "tss_distance_bp"]].isna().all().all()


def test_only_required_inputs(chain: Chain, tmp_path: Path) -> None:
    """Every optional layer missing: the step still writes layers.parquet (layer columns NA) and
    the files it can (genes.gff3, operon_pairs.tsv, tus.tsv)."""
    keep = {k: chain.inputs[k] for k in REQUIRED}
    ann, lay = run(tmp_path, keep)
    for name in ("genes.gff3", "operon_pairs.tsv", "tus.tsv"):
        assert (ann / name).is_file(), name
    df = pd.read_parquet(lay)
    assert len(df) == len(chain.frame("genes"))
    for c in ("promoter_score", "tx_rate_max_tu", "terminator_conf_downstream", "tir_ostir", "cai"):
        assert df[c].isna().all(), c
    assert all(len(list(x)) == 0 for x in df["region_ids"])


# ----------------------------------------------------------------------------- precomputed layers


@pytest.mark.parametrize("key", ["promoters", "promoter_windows"])
def test_precomputed_promoters_tsv(chain: Chain, tmp_path: Path, key: str) -> None:
    """annotate.nf stages pre('promoters') as ``promoters``; the local executor gives the same
    TSV as ``promoter_windows``; promoters.bed then comes from the TSV."""
    from giap.annotate.precomputed import load_precomputed
    from giap.annotate.promoters import promoters_bed
    from giap.io.bed import write_bed

    tsv = chain.sim.precomputed["promoters"]
    cfg = GiapConfig.model_validate(
        {"seed": 1, "annotation": {"precomputed": {"promoters": str(tsv)}}}
    )
    ann, lay = run(tmp_path / "out", {**chain.subset("promoter_windows"), key: tsv}, cfg)
    reps = orc.replicons(chain.frame("replicons"))
    frame = load_precomputed("promoters", tsv, reps)[["contig", "strand", "end3p", "score"]]
    want = tmp_path / "promoters.bed"
    write_bed(want, promoters_bed(frame, reps, cfg.annotation.promoters.promotech_yes), reps)
    assert (ann / "promoters.bed").read_bytes() == want.read_bytes()
    assert pd.read_parquet(lay)["promoter_score"].notna().any()


def test_local_executor_style_precomputed_layers(chain: Chain, tmp_path: Path) -> None:
    """Every precomputed TSV under the keys the local executor uses (all outputs of a replaced
    step point to the layer file)."""
    pre = chain.sim.precomputed
    layers = ("promoters", "tx_rates", "terminators", "rbs", "regions")
    cfg = GiapConfig.model_validate(
        {"seed": 1, "annotation": {"precomputed": {k: str(pre[k]) for k in layers}}}
    )
    inputs = {
        **{k: chain.inputs[k] for k in ("genes", "features", "replicons", "pairs_raw", "tus",
                                        "codon", "genome")},
        "promoter_windows": pre["promoters"], "tx_rates": pre["tx_rates"],
        "terminators": pre["terminators"], "terminators_crosscheck": pre["terminators"],
        "rbs": pre["rbs"], "regions": pre["regions"], "region_genes": pre["regions"],
    }  # fmt: skip
    ann, lay = run(tmp_path, inputs, cfg)
    for name in ANNOTATION_FILES:
        if name != "terminators_crosscheck.tsv":
            assert (ann / name).is_file(), name
    assert not (ann / "terminators_crosscheck.tsv").exists(), "TransTermHP did not run"
    df = pd.read_parquet(lay)
    ids = {x for cell in df["region_ids"] for x in cell}
    from giap.io.tsv import read_tsv

    head, rows = read_tsv(pre["regions"])
    planted = {str(r[0]) for r in rows if r[head.index("genes")]}
    assert ids == planted, "region_ids from the precomputed regions' gene lists"
    for c in ("tir_ostir", "promoter_score", "terminator_conf_downstream"):
        assert df[c].notna().any(), c
    head, _ = read_tsv(ann / "regions.tsv")
    assert tuple(head) == SCHEMA_TABLES["regions"]


# ----------------------------------------------------------------------------- native copies


def test_copy_native_byte_copies(tmp_path: Path) -> None:
    """annotation/tools/<tool>/: byte copies, no gzip, no rename (§9.1 "untouched")."""
    from giap.annotate.layers import copy_native

    srcs = {
        "genomad": REC / "genomad" / "1.12.0" / "sim100k",
        "terminatornet": REC / "terminatornet" / "721803095459" / "sim100k",
        "promotech": REC / "promotech" / "56251ad9" / "annotate-real-1847900",
    }
    dest = tmp_path / "annotation" / "tools"
    copy_native(srcs, dest)
    assert sorted(p.name for p in dest.iterdir()) == sorted(srcs)
    for tool, src in srcs.items():
        want = {str(p.relative_to(src)): p.read_bytes() for p in src.rglob("*") if p.is_file()}
        got = {str(p.relative_to(dest / tool)): p.read_bytes()
               for p in (dest / tool).rglob("*") if p.is_file()}  # fmt: skip
        assert got == want, tool


def test_no_tool_dirs_no_native_copies(full: tuple[Path, Path]) -> None:
    tools = full[0] / "tools"
    assert not tools.exists() or not any(tools.iterdir())


# ----------------------------------------------------------------------------- schemas


@pytest.mark.parametrize("table", sorted(SCHEMA_TABLES))
def test_schema_fields_match_the_columns(table: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    d = _schema(table)
    assert tuple(f["name"] for f in d["fields"]) == SCHEMA_TABLES[table]
    assert "NA" in d.get("missingValues", [])
    frictionless.Schema.from_descriptor({k: v for k, v in d.items() if k != "coordinates"})
    if table in ("regions", "tus", "id_map"):
        assert d.get("coordinates") == "1-based-closed", table


@pytest.mark.parametrize("table", ["regions", "operon_pairs", "tus", "terminators_crosscheck"])
def test_example_9_written_tsvs_validate(full: tuple[Path, Path], table: str) -> None:
    _frictionless_valid(full[0] / f"{table}.tsv", table)


def test_regions_schema_enums(full: tuple[Path, Path]) -> None:
    fields = {f["name"]: f for f in _schema("regions")["fields"]}
    assert set(fields["region_type"]["constraints"]["enum"]) == {
        "prophage", "plasmid", "bgc", "crispr", "amr",
    }  # fmt: skip
    assert set(_schema("operon_pairs")["fields"][3]["constraints"]["enum"]) >= {
        "uniop", "distance_model",
    }  # fmt: skip


def test_copies_of_inputs_are_not_moved(chain: Chain, tmp_path: Path) -> None:
    """The step reads its inputs; it never moves or rewrites them."""
    before = {k: p.read_bytes() for k, p in chain.inputs.items() if p.is_file()}
    run(tmp_path, dict(chain.inputs))
    assert {k: p.read_bytes() for k, p in chain.inputs.items() if p.is_file()} == before
