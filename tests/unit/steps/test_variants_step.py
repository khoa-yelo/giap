"""M7-04 the ``variants`` step on a giapsim producer chain
(docs/specs/M7-04-variants-step-and-run.md: rules 1-4 and 7, examples 1, 3, 4, 5 and 7 at step
level, every degenerate row; lead decisions 1-3; GIAP_SPEC.md §6.2 step row, §6.5, §6.7, §7.5.5,
§9.1).

Interface pinned here (test-author reading of the spec):

* ``@step("variants", inputs={mapping, genome, genes, replicons, peptides},
  optional_inputs={hits_tolerant, reference_peptides}, outputs={pileup, variants,
  variants_vcf})``: kind py, module core, no precomputed layer, no optional outputs.
* Rule 2's order on the step's own inputs. ``peptide_hits.tsv`` of the ``mapping`` directory
  (1-based closed, §6.5) becomes 0-based half-open; M7-01 ``align_hits`` (with ``hits_tolerant``)
  -> M7-02 ``eligible`` -> ``pileup`` (the frame with its ``attrs``, lead decision 2) ->
  ``sites`` -> M7-03 ``vcf_records`` -> ``write_pileup``, ``write_variants`` and ``write_vcf``.
  The oracle below calls exactly that chain. Per peptide it builds the M7-02 ``Support`` from the
  ``peptide_hits`` row (``source``, ``placement_q``, ``uniqueness``) and ``search_db`` from C9.
  Under the default ``genome_relation: related`` the oracle's ``search_db`` choice cannot change
  the result. The step's three files must equal the oracle's byte for byte.
* VCF header: ``giap_version = giap.__version__``; ``genome_sha256`` is the SHA-256 of the
  ``genome`` input file's bytes (as ``loci`` computes it for ``mappings.probed``); the contigs
  follow the ``replicons`` table order; the sample is ``cfg.sample_name``. The step never
  resolves ``auto``: ``giap run`` does that (lead decision 1).
* ``pileup.tsv.gz`` is gzip whatever the output path is called. The local executor's output path
  is ``work/pileup`` with no ``.gz``, and the publisher copies its bytes to
  ``tables/pileup.tsv.gz``.
* Rule 1, no ``hits_tolerant``: frameshift rows are left out (the chain runs on the other rows).
  A WARNING names the run flag ``layer_missing:hits_tolerant``.
* Rule 7: one §6.7 event ``variants.summary`` on the step logger (``extra={"event": ...,
  "data": ...}``) whose plain-JSON ``data`` holds ``aligned_hits`` (number of alignments),
  ``covered_codons`` (pileup rows), ``sites_by_call`` and ``sites_by_confidence`` (value ->
  number of ``variants.tsv`` rows; values with no rows may be absent or 0).

The fixture is a 50-kb giapsim genome with a strain layer (SNVs and 1-nt indels, M3-18). As in
tests/integration/test_align_sim100k.py, the peptides come from the sample strain and are mapped
to the reference genome, so the tolerant route gives substitution and frameshift rows. Preflight
and the mapping chain run through ``run_step`` (``helpers.mapping_steps.drive``) with a test
``$GIAP_DB_DIR``, never the real one.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import json
import logging
import shutil
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq
import pytest
from helpers.mapping_steps import drive, fresh_steps, preflight, roots_of

import giap
from giap.config import GiapConfig, resolve_config
from giap.variants.pileup import PILEUP_COLUMNS
from giap.variants.vcf_records import VARIANTS_COLUMNS

SEED = 11
SCHEMA_DIR = Path(__file__).resolve().parents[3] / "src" / "giap" / "schemas" / "tables"
REQUIRED = {"mapping", "genome", "genes", "replicons", "peptides"}
OPTIONAL = {"hits_tolerant", "reference_peptides", "peptide_spectra"}  # M7-10 C17
OUTPUTS = {
    "pileup": "tables/pileup.tsv.gz",
    "variants": "tables/variants.tsv",
    "variants_vcf": "tables/variants.vcf",
    "variants_qc": "work/variants_qc.parquet",  # M7-14: unpublished work intermediate
}
STEP_KEYS = ("mapping", "genome", "genes", "replicons", "peptides", "hits_tolerant")
EVENT = "variants.summary"


def _cfg(**over: Any) -> GiapConfig:
    return GiapConfig.model_validate({"seed": SEED, "sample_name": "S1", **over})


# ----------------------------------------------------------------------------- fixtures


def _chain(root: Path, *, strain: bool) -> dict[str, Path]:
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SimConfig, simulate
    from giapsim.strain import StrainConfig

    st = StrainConfig(snv_rate_per_kb=5.0, indel_rate_per_kb=0.5, indel_len=(1, 1))
    cfg_sim = SimConfig(genome_length=50_000, proteomics=ProteomicsConfig(),
                        strain=st if strain else None, write_precomputed=False)  # fmt: skip
    res = simulate(cfg_sim, SEED, root / "sim")
    if strain:
        assert res.reference_fasta is not None and res.reference_json is not None
        res = dataclasses.replace(
            res, genome_fasta=res.reference_fasta, bakta_json=res.reference_json
        )
    run = preflight(res, root / "pf")
    cfg = resolve_config(run.config_yaml, [f"seed={SEED}"])
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, root / "work", roots_of(run))
    out = {k: avail[k] for k in STEP_KEYS}
    out["strain_truth"] = root / "sim" / "truth" / "strain_differences.tsv"
    return out


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """Strain layer: sample-strain peptides on the reference genome."""
    return _chain(Path(tmp_path_factory.mktemp("variants_chain")), strain=True)


@pytest.fixture(scope="module")
def same_chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """Example 3: database peptides of the genome's own proteins (``search_db: same``)."""
    return _chain(Path(tmp_path_factory.mktemp("variants_same")), strain=False)


def step_inputs(chain: Mapping[str, Path], *, drop: tuple[str, ...] = ()) -> dict[str, Path]:
    return {k: chain[k] for k in STEP_KEYS if k not in drop}


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


class Run:
    def __init__(self, paths: dict[str, Path], records: list[logging.LogRecord],
                 inputs: dict[str, Path]) -> None:  # fmt: skip
        self.paths = paths
        self.records = records
        self.inputs = inputs

    def messages(self) -> list[str]:
        return [r.getMessage() for r in self.records]

    def frame(self, key: str) -> pd.DataFrame:
        return pd.read_csv(self.paths[key], sep="\t", dtype=str, keep_default_na=False)


def run_variants_step(
    inputs: Mapping[str, Path],
    out: Path,
    cfg: GiapConfig,
    outputs: Mapping[str, Path] | None = None,
) -> Run:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    paths = dict(outputs) if outputs is not None else {k: out / v for k, v in OUTPUTS.items()}
    rec = _Recorder()
    log = logging.getLogger("giap.steps.variants")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        reg.run_step("variants", cfg, dict(inputs), paths)
    finally:
        log.removeHandler(rec)
        log.setLevel(old)
    return Run(paths, rec.records, dict(inputs))


# ----------------------------------------------------------------------------- oracle


def _hits(mapping: Path) -> pd.DataFrame:
    ph = pd.read_csv(mapping / "peptide_hits.tsv", sep="\t", dtype=str, keep_default_na=False,
                     na_values=["NA"])  # fmt: skip
    for col in ("start", "end", "n_diff"):
        ph[col] = ph[col].astype("int64")
    ph["start"] -= 1  # §6.5: 1-based closed -> 0-based half-open
    ph["placement_q"] = ph["placement_q"].astype(float)
    return ph


def _is_fs(v: object) -> bool:
    return isinstance(v, str) and "frameshift" in v.split(";")


def _replicons(path: Path) -> dict[str, Any]:
    from giap.model import Replicon

    rows = pq.read_table(path).to_pylist()
    return {r["contig"]: Replicon(r["contig"], r["length"], r["topology"], r["kind"]) for r in rows}


def _gene_model(genes: Path, reps: Mapping[str, Any]) -> Any:
    from giap.model.coords import Interval
    from giap.model.genes import Gene, GeneModel

    gs = tuple(
        Gene(r["locus_tag"], r["feature_type"],
             Interval(r["contig"], r["start"], r["end"], r["strand"]), r["aa"] or "", r["nt"],
             r["start_codon"], r["gene"], r["product"] or "", bool(r["partial"]),
             bool(r["hypothetical"]))
        for r in pq.read_table(genes).to_pylist()
    )  # fmt: skip
    return GeneModel(replicons=dict(reps), genes=gs, features=())


def _genome(path: Path) -> dict[str, str]:
    from giap.io.fasta import read_genome_fasta

    return {r.id: r.sequence for r in read_genome_fasta(path)}


def oracle(inputs: Mapping[str, Path], cfg: GiapConfig, out: Path) -> dict[str, Any]:
    """Rule 2's chain on the step's own inputs, written with the M7-03 writers."""
    from giap.io.variant_tables import write_pileup, write_variants
    from giap.io.vcf import write_vcf
    from giap.variants.align import align_hits
    from giap.variants.pileup import Support, eligible, pileup, sites
    from giap.variants.vcf_records import FORMAT_FIELDS, INFO_FIELDS, vcf_records

    ph = _hits(inputs["mapping"])
    ht = pd.read_parquet(inputs["hits_tolerant"]) if "hits_tolerant" in inputs else None
    if ht is None:
        ph = ph[~ph["diff_types"].map(_is_fs)].reset_index(drop=True)
    genome = _genome(inputs["genome"])
    reps = _replicons(inputs["replicons"])
    gm = _gene_model(inputs["genes"], reps)
    pep = pd.read_parquet(inputs["peptides"])
    sdb = {str(p): str(s) for p, s in zip(pep["peptide"], pep["search_db"], strict=True)}
    support = {}
    for r in ph.to_dict("records"):
        s = sdb.get(r["peptide"], "same")
        support.setdefault(r["peptide"], Support(r["peptide"], r["source"],
                           s if s in ("same", "other") else "same", float(r["placement_q"]),
                           r["uniqueness"]))  # fmt: skip
    aligns = align_hits(ph, ht, genome, reps)
    el = eligible(aligns, support, cfg)
    pile = pileup(el, gm, genome, cfg, support=support)
    st = sites(pile, cfg)
    recs, table = vcf_records(st, pile, genome, reps)
    out.mkdir(parents=True, exist_ok=True)
    paths = {k: out / Path(v).name for k, v in OUTPUTS.items()}
    write_pileup(paths["pileup"], pile)
    write_variants(paths["variants"], table)
    sha = hashlib.sha256(Path(inputs["genome"]).read_bytes()).hexdigest()
    write_vcf(paths["variants_vcf"], recs, list(reps.values()), giap_version=giap.__version__,
              genome_sha256=sha, sample=cfg.sample_name, info_fields=INFO_FIELDS,
              format_fields=FORMAT_FIELDS, sequences=genome)  # fmt: skip
    return {"paths": paths, "aligns": aligns, "pile": pile, "table": table, "records": recs}


# ----------------------------------------------------------------------------- helpers


def _gunzip_text(path: Path) -> str:
    raw = path.read_bytes()
    assert raw[:2] == b"\x1f\x8b", f"{path.name} is not gzip"
    return gzip.decompress(raw).decode("utf-8")


def _vcf_lines(path: Path) -> tuple[list[str], list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return [x for x in lines if x.startswith("#")], [x for x in lines if not x.startswith("#")]


def _chrom_line(path: Path) -> list[str]:
    head, _ = _vcf_lines(path)
    chrom = [x for x in head if x.startswith("#CHROM")]
    assert len(chrom) == 1, head
    return chrom[0].split("\t")


def _validate(path: Path, stem: str) -> Any:
    frictionless = pytest.importorskip("frictionless")
    d = json.loads((SCHEMA_DIR / f"{stem}.schema.json").read_text(encoding="utf-8"))
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    return res.validate()


def _event(run: Run) -> dict[str, Any]:
    found = [r for r in run.records if getattr(r, "event", None) == EVENT]
    assert len(found) == 1, f"expected one {EVENT} event, got {len(found)}: {run.messages()}"
    data = getattr(found[0], "data", None)
    assert isinstance(data, dict), data
    return data


@pytest.fixture(scope="module")
def full(chain: dict[str, Path], tmp_path_factory: pytest.TempPathFactory) -> Run:
    """Example 1 at step level: every input, sample name ``S1``."""
    out = Path(tmp_path_factory.mktemp("variants_full"))
    return run_variants_step(step_inputs(chain), out, _cfg())


@pytest.fixture(scope="module")
def full_oracle(full: Run, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    return oracle(full.inputs, _cfg(), Path(tmp_path_factory.mktemp("variants_oracle")))


# ----------------------------------------------------------------------------- registry


def test_variants_step_is_registered_with_the_spec_keys() -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get("variants")
    assert spec.kind == "py"
    assert spec.module == "core"  # §1.3: Module 5 is core, no --modules switch
    assert spec.precomputed is None
    assert set(spec.inputs) == REQUIRED
    assert set(spec.optional_inputs) == OPTIONAL
    assert set(spec.outputs) == set(OUTPUTS)
    assert spec.optional_outputs == ()


def test_variants_runs_after_its_producers() -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    order = reg.topo_order([s.name for s in reg.all_steps()])
    for producer in ("parse_bakta", "map_tolerant", "loci"):
        assert order.index(producer) < order.index("variants"), producer


def test_local_plan_runs_variants_after_loci_with_peptides() -> None:
    from giap.runner.local import plan_local
    from giap.steps import load_step_modules

    load_step_modules()
    order = plan_local(_cfg(), ["variants"], {"genome", "peptides", "bakta_json"})
    assert "variants" in order, order
    assert order.index("loci") < order.index("variants")
    assert order.index("map_tolerant") < order.index("variants")


def test_local_plan_leaves_variants_out_without_peptides() -> None:
    """Degenerate row "no --peptides": no Module 5 step, and no error (as MAPPING today)."""
    from giap.runner.local import plan_local
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    reg.get("variants")  # registered, so its absence below is the skip, not a missing step
    order = plan_local(_cfg(), ["variants"], {"genome", "bakta_json"})
    assert "variants" not in order, order


# ----------------------------------------------------------------------------- example 1


def test_example_1_three_files_written_with_their_headers(full: Run) -> None:
    for key, path in full.paths.items():
        assert path.is_file(), f"{key} not written"
    pile = _gunzip_text(full.paths["pileup"]).splitlines()
    assert pile[0].split("\t") == list(PILEUP_COLUMNS)
    assert len(pile) > 1, "the fixture covers codons"
    var = full.paths["variants"].read_text(encoding="utf-8").splitlines()
    assert var[0].split("\t") == list(VARIANTS_COLUMNS)
    head, _ = _vcf_lines(full.paths["variants_vcf"])
    assert head[0] == "##fileformat=VCFv4.3"
    assert f"##source=giap-{giap.__version__}" in head


def _tables(paths: dict[str, Path]) -> dict[str, Path]:
    """The published tables; the oracle (``run_variants``) writes no ``variants_qc`` (M7-14)."""
    return {k: v for k, v in paths.items() if k != "variants_qc"}


def _without_flags(key: str, path: Path) -> object:
    """The file with M7-05's flags removed: the ``flags`` column of ``variants.tsv`` and the
    ``FLAGS`` INFO entry of the VCF. M7-05 fills both; test_variants_step_flags.py checks them."""
    if key == "variants":
        rows = [ln.split("\t") for ln in path.read_text(encoding="utf-8").splitlines()]
        j = rows[0].index("flags")
        return [r[:j] + r[j + 1 :] for r in rows]
    if key == "variants_vcf":
        out = []
        for ln in path.read_text(encoding="utf-8").splitlines():
            if not ln.startswith("#"):
                f = ln.split("\t")
                f[7] = ";".join(x for x in f[7].split(";") if not x.startswith("FLAGS=")) or "."
                ln = "\t".join(f)
            out.append(ln)
        return out
    return path.read_bytes()


def test_example_1_tables_equal_the_rule_2_chain(full: Run, full_oracle: dict[str, Any]) -> None:
    for key, path in _tables(full.paths).items():
        want = full_oracle["paths"][key]
        assert _without_flags(key, path) == _without_flags(key, want), key


def test_example_1_fixture_exercises_alt_and_frameshift_sites(full: Run) -> None:
    df = full.frame("variants")
    assert (df["alt_type"] == "substitution").any(), "fixture meant to give substitution sites"
    assert (df["alt_type"] == "frameshift").any(), "frameshift sites need hits_tolerant"
    _, body = _vcf_lines(full.paths["variants_vcf"])
    assert body, "substitution sites give VCF records"


def test_example_1_substitutions_sit_on_planted_snvs(chain: dict[str, Path], full: Run) -> None:
    """The 1-based → 0-based conversion of rule 1: most substitution sites land on the codon of a
    planted strain SNV (``pos_reference`` is 1-based on the reference genome). An off-by-one
    start would move every site."""
    truth = pd.read_csv(chain["strain_truth"], sep="\t")
    snv = {(r.contig, int(r.pos_reference) - 1) for r in truth[truth["type"] == "snv"].itertuples()}
    df = full.frame("variants")
    sub = df[df["alt_type"] != "frameshift"]
    assert len(sub) >= 10

    def on_snv(r: Any) -> bool:
        width = len(r.ref_codon)
        start0 = int(r.codon_start) - 1  # the TSV is 1-based (§6.5)
        lo = start0 if r.strand == "+" or width == 3 else start0 - 3
        return any((r.contig, p) in snv for p in range(lo, lo + width))

    hits = sum(on_snv(r) for r in sub.itertuples())
    assert hits >= 0.8 * len(sub), f"{hits} of {len(sub)} substitution sites on planted SNVs"


def test_example_1_vcf_reads_in_pysam(full: Run, chain: dict[str, Path]) -> None:
    pysam = pytest.importorskip("pysam")
    reps = _replicons(chain["replicons"])
    with pysam.VariantFile(str(full.paths["variants_vcf"])) as vf:
        assert list(vf.header.samples) == ["S1"]
        assert {c: vf.header.contigs[c].length for c in vf.header.contigs} == {
            c: r.length for c, r in reps.items()
        }
        n = sum(1 for _ in vf)
    _, body = _vcf_lines(full.paths["variants_vcf"])
    assert n == len(body)


def test_vcf_reference_is_the_genome_sha256(full: Run, chain: dict[str, Path]) -> None:
    head, _ = _vcf_lines(full.paths["variants_vcf"])
    sha = hashlib.sha256(chain["genome"].read_bytes()).hexdigest()
    assert f"##reference=sha256:{sha}" in head
    contigs = [x for x in head if x.startswith("##contig=")]
    reps = _replicons(chain["replicons"])
    assert contigs == [f"##contig=<ID={r.id},length={r.length}>" for r in reps.values()]


def test_pileup_is_gzip_even_without_a_gz_output_path(
    chain: dict[str, Path], full: Run, tmp_path: Path
) -> None:
    """``--executor local`` stores outputs at ``work/<key>``. The publisher copies the bytes of
    ``work/pileup`` to ``tables/pileup.tsv.gz``, so the step must gzip whatever the path is
    called."""
    outs = {k: tmp_path / "work" / k for k in OUTPUTS}
    run = run_variants_step(step_inputs(chain), tmp_path, _cfg(), outputs=outs)
    assert _gunzip_text(run.paths["pileup"]) == _gunzip_text(full.paths["pileup"])
    assert run.paths["variants"].read_bytes() == full.paths["variants"].read_bytes()
    assert run.paths["variants_vcf"].read_bytes() == full.paths["variants_vcf"].read_bytes()


# ----------------------------------------------------------------------------- example 7


def test_example_7_both_tsvs_validate(full: Run) -> None:
    for key, stem in (("pileup", "pileup"), ("variants", "variants")):
        report = _validate(full.paths[key], stem)
        assert report.valid, (key, report.flatten(["rowNumber", "fieldName", "type", "note"]))


# ----------------------------------------------------------------------------- examples 4, 5


def test_example_4_vcf_sample_column_is_the_config_sample_name(
    chain: dict[str, Path], tmp_path: Path
) -> None:
    run = run_variants_step(step_inputs(chain), tmp_path, _cfg(sample_name="run_A"))
    cols = _chrom_line(run.paths["variants_vcf"])
    assert len(cols) == 10 and cols[-1] == "run_A", cols


def test_example_5_resolved_sample_name_with_spaces_reaches_the_vcf(
    chain: dict[str, Path], tmp_path: Path
) -> None:
    from giap.config import resolve_sample_name

    cfg = _cfg(sample_name="my sample")
    name = resolve_sample_name(cfg, None, None)
    assert name == "my_sample"
    run = run_variants_step(
        step_inputs(chain), tmp_path, cfg.model_copy(update={"sample_name": name})
    )
    assert _chrom_line(run.paths["variants_vcf"])[-1] == "my_sample"


def test_unresolved_auto_never_reaches_the_vcf(chain: dict[str, Path], tmp_path: Path) -> None:
    """Lead decision 1: ``giap run`` resolves ``auto`` before any step runs; a step that is handed
    ``auto`` anyway fails instead of writing a VCF sample called ``auto``."""
    from giap.steps.registry import StepError

    with pytest.raises(StepError):
        run_variants_step(step_inputs(chain), tmp_path, _cfg(sample_name="auto"))
    assert not (tmp_path / "tables" / "variants.vcf").exists()


# ----------------------------------------------------------------------------- example 3


def test_example_3_same_strain_database_peptides_give_no_alt_except_non_atg_starts(
    same_chain: dict[str, Path], tmp_path: Path
) -> None:
    """Database peptides against the strain's own proteins: a ``variants.tsv`` row only where the
    peptide really contradicts the genome, i.e. the M of an annotated non-ATG start."""
    run = run_variants_step(step_inputs(same_chain), tmp_path, _cfg())
    pile = _gunzip_text(run.paths["pileup"]).splitlines()
    assert len(pile) > 100, "the same-strain fixture covers codons"
    genes = pq.read_table(same_chain["genes"]).to_pandas()
    starts = {}
    for g in genes.itertuples():
        first = g.start if g.strand == "+" else g.end - 3  # lowest coordinate of codon 1
        starts[(g.contig, first, g.strand)] = g.start_codon
    df = run.frame("variants")
    for r in df.itertuples():
        key = (r.contig, int(r.codon_start) - 1, r.strand)
        assert key in starts, f"ALT site off any gene start: {r}"
        assert starts[key] != "ATG" and r.alt_aa == "M", r
    want = oracle(run.inputs, _cfg(), tmp_path / "oracle")
    for key, path in _tables(run.paths).items():
        assert _without_flags(key, path) == _without_flags(key, want["paths"][key]), key


# ----------------------------------------------------------------------------- rule 7


def test_rule_7_summary_event_counts(full: Run, full_oracle: dict[str, Any]) -> None:
    data = _event(full)
    json.dumps(data)  # §6.7 events are plain JSON (no numpy scalars)
    assert data["aligned_hits"] == len(full_oracle["aligns"])
    assert data["covered_codons"] == len(full_oracle["pile"])
    df = full.frame("variants")
    for col, key in (("call", "sites_by_call"), ("confidence", "sites_by_confidence")):
        want = df[col].value_counts().to_dict()
        got = {k: v for k, v in data[key].items() if v}
        assert got == want, (key, got, want)
    assert any(r.levelno == logging.INFO and str(len(full_oracle["pile"])) in r.getMessage()
               for r in full.records), full.messages()  # fmt: skip


# ----------------------------------------------------------------------------- degenerate rows


def test_no_accepted_placement_gives_header_only_outputs(
    chain: dict[str, Path], tmp_path: Path
) -> None:
    mapping = tmp_path / "mapping"
    shutil.copytree(chain["mapping"], mapping)
    hits = mapping / "peptide_hits.tsv"
    hits.write_text(hits.read_text(encoding="utf-8").splitlines()[0] + "\n", encoding="utf-8")
    ins = step_inputs(chain) | {"mapping": mapping}
    run = run_variants_step(ins, tmp_path / "out", _cfg())
    assert _gunzip_text(run.paths["pileup"]) == "\t".join(PILEUP_COLUMNS) + "\n"
    assert run.paths["variants"].read_text(encoding="utf-8") == "\t".join(VARIANTS_COLUMNS) + "\n"
    head, body = _vcf_lines(run.paths["variants_vcf"])
    assert body == []
    assert head[0] == "##fileformat=VCFv4.3"
    assert _chrom_line(run.paths["variants_vcf"])[-1] == "S1"
    data = _event(run)
    assert data["aligned_hits"] == 0 and data["covered_codons"] == 0


def test_without_hits_tolerant_frameshift_rows_are_skipped_with_a_run_flag(
    chain: dict[str, Path], tmp_path: Path
) -> None:
    ph = _hits(chain["mapping"])
    assert ph["diff_types"].map(_is_fs).any(), "fixture meant to have frameshift rows"
    ins = step_inputs(chain, drop=("hits_tolerant",))
    run = run_variants_step(ins, tmp_path / "out", _cfg())
    warned = [r for r in run.records if "layer_missing:hits_tolerant" in r.getMessage()]
    assert warned and all(r.levelno >= logging.WARNING for r in warned), run.messages()
    df = run.frame("variants")
    assert not (df["alt_type"] == "frameshift").any()
    assert (df["alt_type"] == "substitution").any(), "the other rows still count"
    want = oracle(ins, _cfg(), tmp_path / "oracle")
    for key, path in _tables(run.paths).items():
        assert _without_flags(key, path) == _without_flags(key, want["paths"][key]), key


def test_hits_that_disagree_with_the_genome_fail_the_step(
    chain: dict[str, Path], tmp_path: Path
) -> None:
    """M7-01 raises; the step fails (``giap run`` then exits 4 naming ``variants``)."""
    from giap.steps.registry import StepError

    comp = str.maketrans("ACGT", "TGCA")
    fa = tmp_path / "genome.fa"
    fa.write_text(
        "".join(f">{c}\n{s.translate(comp)}\n" for c, s in _genome(chain["genome"]).items()),
        encoding="utf-8",
    )
    ins = step_inputs(chain) | {"genome": fa}
    with pytest.raises(StepError) as err:
        run_variants_step(ins, tmp_path / "out", _cfg())
    assert err.value.name == "variants"
    log = err.value.write_log(tmp_path / "logs")
    assert log == tmp_path / "logs" / "variants.log"
    for key in OUTPUTS:
        assert not (tmp_path / "out" / OUTPUTS[key]).exists(), key


def test_te7_same_inputs_give_byte_identical_files(
    chain: dict[str, Path], full: Run, tmp_path: Path
) -> None:
    again = run_variants_step(step_inputs(chain), tmp_path, _cfg(seed=SEED + 1))
    for key in OUTPUTS:
        assert again.paths[key].read_bytes() == full.paths[key].read_bytes(), key


def test_reference_peptides_is_accepted_and_unused_in_pass_1(
    chain: dict[str, Path], full: Run, tmp_path: Path
) -> None:
    """Declared now, used only by M7-10: staging it changes nothing."""
    ref = tmp_path / "reference_peptides.tsv"
    ref.write_text("peptide\tq_value\nPEPTIDEK\t0.001\n", encoding="utf-8")
    run = run_variants_step(step_inputs(chain) | {"reference_peptides": ref}, tmp_path, _cfg())
    for key in OUTPUTS:
        assert run.paths[key].read_bytes() == full.paths[key].read_bytes(), key
