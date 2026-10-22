"""M7-04 through a real ``giap run`` on sim100k (docs/specs/M7-04-variants-step-and-run.md
examples 1-4 and 7, degenerate rows "two runs, same seed" and "the step fails", lead decisions 1
and 2; GIAP_SPEC.md §6.4 "byte-identical tables", §7.5.5, §9.1, §9.5).

Inputs: the giapsim SIM100K genome with M3-19 proteomics, M3-20 de novo reads and every
precomputed layer (as in tests/integration/test_m2_21_local_vs_nextflow.py), ``--modules core``
and one seed. The test DB comes from tests/helpers/mapping_steps.py, never the real one.

* Example 1 + 4: ``--executor local`` with the de novo table gzipped as ``run_A.tsv.gz`` and
  ``sample_name`` left at ``auto``. It exits 0 and writes ``OUT/tables/pileup.tsv.gz`` (gzip),
  ``tables/variants.tsv`` (rows: the de novo errors give ALT sites) and ``tables/variants.vcf``.
  The VCF reads in pysam with the sample ``run_A``. ``OUT/logs/events.jsonl`` holds one
  ``variants.summary`` event.
* Example 7: both TSVs validate against their schemas. Lead decision 2: the run writes
  ``OUT/schemas/datapackage.json`` (§9.1 layout) with the ``pileup`` and ``variants`` resources,
  whose ``tables/...`` paths resolve against ``OUT``.
* TE7: a second run with the same seed gives byte-identical files.
* Example 3: the database table (``search_db: same``, the genome's own proteins) gives
  ``variants.tsv`` rows only at annotated non-ATG starts (M read at a GTG/TTG codon).
* A failing ``variants`` step: ``giap run`` exits 4 and names ``variants`` and
  ``logs/variants.log``.
* Example 2: ``--executor nextflow --profile farmshare_local`` gives byte-identical files
  (needs Nextflow).
"""

from __future__ import annotations

import dataclasses
import gzip
import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner

pytestmark = [pytest.mark.slow]

SEED = 7
TABLES = ("tables/pileup.tsv.gz", "tables/variants.tsv", "tables/variants.vcf")
SCHEMA_DIR = Path(__file__).resolve().parents[2] / "src" / "giap" / "schemas" / "tables"


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(giapsim.SIM100K, proteomics=ProteomicsConfig(),
                              denovo=DenovoConfig(), write_precomputed=True)  # fmt: skip
    return giapsim.simulate(cfg, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")


@pytest.fixture(scope="module")
def run_a(sim: Any, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The de novo table, gzipped, under the name ``run_A.tsv.gz`` (example 4)."""
    src = Path(sim.outdir) / "peptides_denovo.tsv"
    assert src.is_file(), "giapsim writes peptides_denovo.tsv with SimConfig.denovo"
    dst = tmp_path_factory.mktemp("in") / "run_A.tsv.gz"
    dst.write_bytes(gzip.compress(src.read_bytes(), mtime=0))
    return dst


@pytest.fixture(scope="module")
def db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.mapping_steps import write_test_db

    return Path(write_test_db(tmp_path_factory.mktemp("test_db_m7_04")))


def _invoke(sim: Any, peptides: Path, out: Path, db: Path, *extra: str) -> Any:
    from giap.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(peptides)]
    args += ["--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", "core", "--seed", str(SEED), *extra]
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIAP_DB_DIR", str(db))
        mp.setenv("NXF_ANSI_LOG", "false")
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _run(sim: Any, peptides: Path, out: Path, db: Path, *extra: str) -> Path:
    res = _invoke(sim, peptides, out, db, *extra)
    assert res.exit_code == 0, res.output
    return out


@pytest.fixture(scope="module")
def local(sim: Any, run_a: Path, db: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _run(sim, run_a, tmp_path_factory.mktemp("LOCAL") / "OUT", db, "--executor", "local")


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _validate(path: Path, stem: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    d = json.loads((SCHEMA_DIR / f"{stem}.schema.json").read_text(encoding="utf-8"))
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, (stem, report.flatten(["rowNumber", "fieldName", "type", "note"]))


def _digests(out: Path) -> dict[str, bytes]:
    return {rel: (out / rel).read_bytes() for rel in TABLES}


# ----------------------------------------------------------------------------- examples 1, 4


def test_example_1_three_tables_published(local: Path) -> None:
    for rel in TABLES:
        assert (local / rel).is_file(), f"{rel} not published"
    raw = (local / "tables" / "pileup.tsv.gz").read_bytes()
    assert raw[:2] == b"\x1f\x8b", "tables/pileup.tsv.gz is not gzip"
    from giap.variants.pileup import PILEUP_COLUMNS
    from giap.variants.vcf_records import VARIANTS_COLUMNS

    pile = gzip.decompress(raw).decode("utf-8").splitlines()
    assert pile[0].split("\t") == list(PILEUP_COLUMNS)
    assert len(pile) > 1
    var = _read(local / "tables" / "variants.tsv")
    assert var.columns.tolist() == list(VARIANTS_COLUMNS)
    assert len(var) > 0, "de novo errors give ALT sites on sim100k"
    assert {"swap", "block"} & set(var["alt_type"]), "de novo swaps/blocks expected"


def test_example_1_and_4_vcf_reads_in_pysam_with_sample_run_a(local: Path) -> None:
    pysam = pytest.importorskip("pysam")
    vcf = local / "tables" / "variants.vcf"
    with pysam.VariantFile(str(vcf)) as vf:
        assert list(vf.header.samples) == ["run_A"]
        recs = list(vf)
    assert recs, "ALT sites give VCF records"
    for r in recs:
        assert r.samples["run_A"]["DP"] >= 1


def test_example_1_summary_event_in_the_run_log(local: Path) -> None:
    lines = (local / "logs" / "events.jsonl").read_text(encoding="utf-8").splitlines()
    events = [json.loads(x) for x in lines if x.strip()]
    found = [e for e in events if e.get("event") == "variants.summary"]
    assert len(found) == 1, [e.get("event") for e in events]
    data = found[0]["data"]
    var = _read(local / "tables" / "variants.tsv")
    assert sum(int(v) for v in data["sites_by_call"].values()) == len(var)


# ----------------------------------------------------------------------------- example 7


def test_example_7_tsvs_validate(local: Path) -> None:
    _validate(local / "tables" / "pileup.tsv.gz", "pileup")
    _validate(local / "tables" / "variants.tsv", "variants")


def test_decision_2_datapackage_resolves_against_the_run_root(local: Path) -> None:
    dp = local / "schemas" / "datapackage.json"
    assert dp.is_file(), "OUT/schemas/datapackage.json not written"
    res = {r["name"]: r for r in json.loads(dp.read_text(encoding="utf-8"))["resources"]}
    for name, rel in (("pileup", "tables/pileup.tsv.gz"), ("variants", "tables/variants.tsv")):
        assert name in res, f"datapackage lists no {name} resource"
        assert res[name]["path"] == rel
        assert (local / res[name]["path"]).is_file(), rel


# ----------------------------------------------------------------------------- TE7


def test_te7_second_run_byte_identical(
    sim: Any, run_a: Path, db: Path, local: Path, tmp_path: Path
) -> None:
    again = _run(sim, run_a, tmp_path / "OUT", db, "--executor", "local")
    assert _digests(again) == _digests(local)


# ----------------------------------------------------------------------------- example 3


def test_example_3_same_strain_database_peptides(sim: Any, db: Path, tmp_path: Path) -> None:
    out = _run(sim, Path(sim.peptides), tmp_path / "OUT", db, "--executor", "local")
    for rel in TABLES:
        assert (out / rel).is_file(), rel
    genes = pd.read_parquet(out / "work" / "genes")
    starts = {}
    for g in genes.itertuples():
        first = g.start if g.strand == "+" else g.end - 3
        starts[(g.contig, int(first), g.strand)] = g.start_codon
    for r in _read(out / "tables" / "variants.tsv").itertuples():
        key = (r.contig, int(r.codon_start) - 1, r.strand)
        assert key in starts and starts[key] != "ATG" and r.alt_aa == "M", r
    pysam = pytest.importorskip("pysam")
    with pysam.VariantFile(str(out / "tables" / "variants.vcf")) as vf:
        assert list(vf.header.samples) == ["peptides"]  # auto = the peptide file stem


# ----------------------------------------------------------------------------- step failure


def test_failing_variants_step_exits_4_naming_it(
    sim: Any, run_a: Path, db: Path, tmp_path: Path
) -> None:
    from giap.exitcodes import ExitCode
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get("variants")

    def boom(_ctx: Any) -> None:
        raise ValueError("peptide PEPTIDEK: n_diff 0, recomputed 3")

    out = tmp_path / "OUT"
    with pytest.MonkeyPatch.context() as mp:
        mp.setitem(reg._REGISTRY, "variants", dataclasses.replace(spec, func=boom))
        res = _invoke(sim, run_a, out, db, "--executor", "local")
    assert res.exit_code == ExitCode.STEP_FAILED, res.output
    assert "variants" in res.output
    assert (out / "logs" / "variants.log").is_file()
    assert "variants.log" in res.output


# ----------------------------------------------------------------------------- example 2


@pytest.mark.needs_nextflow
def test_example_2_nextflow_gives_byte_identical_files(
    sim: Any, run_a: Path, db: Path, local: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    if shutil.which("nextflow") is None:
        pytest.skip("nextflow not on PATH")
    nf = _run(sim, run_a, tmp_path_factory.mktemp("NF") / "OUT", db,
              "--executor", "nextflow", "--profile", "farmshare_local")  # fmt: skip
    for rel in TABLES:
        assert (nf / rel).is_file(), f"Nextflow run did not publish {rel}"
    assert _digests(nf) == _digests(local)
