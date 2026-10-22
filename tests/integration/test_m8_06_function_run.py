"""M8-06 examples 1, 2, 3 and 9 through ``giap run`` on sim100k (docs/specs/M8-06-function-step-
and-run.md, 06b; GIAP_SPEC.md §6.3 module switch, §6.4 "byte-identical tables", §9.1).

Inputs: the M6-08 / M2-21 executor-equality inputs (``test_m6_08_operons_run.py``): giapsim
SIM100K with M3-19 proteomics and every precomputed layer, seed 7. ``$GIAP_DB_DIR`` is the
synthetic function DB of ``tests/helpers/function_db.py`` (Q6, lead decision 1), merged with the
M3-17 test DB (``helpers.mapping_steps.write_test_db``: every preflight requirement listed); never
the real databases.

Test-author choice (reported to the lead): SIM100K's Bakta-like JSON carries no ``db_xrefs`` at
all (``layout_from_simulation(sim).xrefs`` is empty), so with the spec's DB alone every
collection would be empty and example 1's ">= 1 row" could not hold. The fixture therefore writes
deterministic synthetic ``db_xrefs`` into the simulated JSON before either run (COG/GO/KEGG IDs by
CDS index; 1 CDS in 6 left without), and the synthetic DB is built from those same IDs. Loci,
coordinates and precomputed layers are untouched.

* Example 1: ``--executor local --modules core,function``: the three ``tables/*.tsv`` exist and
  validate (example 9); summary and detection have >= 1 row; every detection and GSEA row's
  ``flags`` has ``ci_unestimated``; every detection row has ``ora_uncalibrated`` and
  ``tu_recal_skipped`` (lead decision 3); the ``function.run_flags`` event has
  ``tu_recal_skipped`` (constant model, Q10); ``genetable/tu_folds.parquet`` is published and
  header-only; the five resources were staged under their item names.
* Example 3: ``--modules core``: no function table, nothing staged, ``tables/operons.tsv`` still
  written.
* Example 2: ``--executor nextflow --profile farmshare_local``: the three tables and
  ``genetable/tu_folds.parquet`` byte-identical to example 1.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner

pytestmark = [pytest.mark.slow]

SEED = 7
TABLES = ("category_summary.tsv", "enrichment_detection.tsv", "enrichment_gsea.tsv")
RESOURCES = ("go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names", "cog2020_fun")
SCHEMA_DIR = Path(__file__).resolve().parents[2] / "src" / "giap" / "schemas" / "tables"


def _xrefs_for(i: int) -> list[str]:
    """Synthetic Bakta ``db_xrefs`` of the i-th CDS: 3 COGs (J, K, L), 12 GO leaves (every
    namespace/parent of the helper's OBO), 5 KOs (5 synthetic pathways); 1 in 6 has none."""
    if i % 6 == 5:
        return []
    return [f"COG:COG{i % 3:04d}", f"GO:{i % 12 + 1:07d}", f"KEGG:K{i % 5 + 1:05d}"]


def _add_xrefs(bakta_json: Path) -> None:
    data = json.loads(bakta_json.read_text(encoding="utf-8"))
    cds = [f for f in data["features"] if f.get("type") in ("cds", "sorf")]
    for i, f in enumerate(cds):
        f["db_xrefs"] = _xrefs_for(i)
    bakta_json.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim
    from giapsim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        giapsim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    res = giapsim.simulate(cfg, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")
    _add_xrefs(Path(res.bakta_json))
    return res


@pytest.fixture(scope="module")
def db(sim: Any, tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.function_db import ITEMS, PATHS, write_function_db
    from helpers.mapping_steps import write_test_db

    from giapsim.layout import layout_from_simulation

    xrefs = layout_from_simulation(sim).xrefs
    assert {"COG", "GO", "KEGG"} <= set(xrefs["db"]), "the fixture's db_xrefs were not parsed"
    root = tmp_path_factory.mktemp("function_db") / "db"
    write_function_db(root, xrefs)
    write_test_db(root)  # every requirement of every module, per-name manifest form
    doc = json.loads((root / "db_manifest.json").read_text(encoding="utf-8"))
    doc |= {n: {"version": "synthetic", "path": PATHS[n]} for n in ITEMS}
    (root / "db_manifest.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return root


def _run(sim: Any, out: Path, db: Path, modules: str, *extra: str) -> None:
    from giap.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides)]
    assert sim.abundance is not None, "the fixture writes abundance.tsv"
    args += ["--abundance", str(sim.abundance)]
    args += ["--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", modules, "--seed", str(SEED), *extra]
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIAP_DB_DIR", str(db))
        mp.setenv("NXF_ANSI_LOG", "false")
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output


@pytest.fixture(scope="module")
def local(sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("LOCAL") / "OUT"
    _run(sim, out, db, "core,function", "--executor", "local")
    return out


@pytest.fixture(scope="module")
def local_core(sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("LOCAL_CORE") / "OUT"
    _run(sim, out, db, "core", "--executor", "local")
    return out


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _tokens(cell: str) -> set[str]:
    return {t for t in cell.split(";") if t and t != "NA"}


def _validate(path: Path, stem: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    d = json.loads((SCHEMA_DIR / f"{stem}.schema.json").read_text(encoding="utf-8"))
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent),
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, (stem, report.flatten(["rowNumber", "fieldName", "type", "note"]))


def _events(out: Path, name: str) -> list[dict[str, Any]]:
    path = out / "logs" / "events.jsonl"
    assert path.is_file(), f"no {path}"
    evs = [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    return [e for e in evs if e.get("event") == name]


# ---------------------------------------------------------------- example 1 (and 9)


@pytest.mark.parametrize("name", TABLES)
def test_example_1_tables_written_and_valid(local: Path, name: str) -> None:
    path = local / "tables" / name
    assert path.is_file(), f"tables/{name} not published under --executor local"
    _validate(path, name.removesuffix(".tsv"))


def test_example_1_header_rows_equal_the_module_columns(local: Path) -> None:
    from giap.function import gsea, ora, summary

    for name, cols in zip(TABLES, (summary.COLUMNS, ora.COLUMNS, gsea.COLUMNS), strict=True):
        head = (local / "tables" / name).read_text(encoding="utf-8").splitlines()[0]
        assert head == "\t".join(cols), name


def test_example_1_summary_and_detection_have_rows(local: Path) -> None:
    summ = _read(local / "tables" / "category_summary.tsv")
    det = _read(local / "tables" / "enrichment_detection.tsv")
    assert len(summ) >= 1, "category_summary.tsv is header-only"
    assert len(det) >= 1, "enrichment_detection.tsv is header-only"
    assert "cog" in set(summ["collection"]), set(summ["collection"])


def test_example_1_row_flags(local: Path) -> None:
    det = _read(local / "tables" / "enrichment_detection.tsv")
    for cell in det["flags"]:
        assert {"ci_unestimated", "ora_uncalibrated", "tu_recal_skipped"} <= _tokens(cell), cell
    gs = _read(local / "tables" / "enrichment_gsea.tsv")
    for cell in gs["flags"]:
        assert "ci_unestimated" in _tokens(cell), cell


def test_example_1_run_flags_event_has_tu_recal_skipped(local: Path) -> None:
    evs = _events(local, "function.run_flags")
    assert len(evs) == 1, evs
    flags = evs[0]["data"]["flags"]
    assert "tu_recal_skipped" in flags, flags
    for missing in ("layer_missing:cog", "layer_missing:kegg", "layer_missing:go"):
        assert missing not in flags, f"{missing}: the resources were not staged ({flags})"


def test_example_1_resources_staged_under_item_names(local: Path, db: Path) -> None:
    from helpers.function_db import PATHS

    for name in RESOURCES:
        staged = local / "work" / name
        assert staged.is_file(), f"work/{name} not staged"
        assert staged.read_bytes() == (db / PATHS[name]).read_bytes(), name


def test_example_1_tu_folds_published_and_header_only(local: Path) -> None:
    """sim100k is a constant detectability model (Q10): ``tu_folds.parquet`` has no rows."""
    path = local / "genetable" / "tu_folds.parquet"
    assert path.is_file(), "genetable/tu_folds.parquet not published under --executor local"
    df = pd.read_parquet(path)
    assert list(df.columns) == ["tu_id", "fold"]
    assert len(df) == 0


# ---------------------------------------------------------------- example 3


def test_example_3_modules_core_writes_no_function_table(local_core: Path) -> None:
    for name in TABLES:
        assert not (local_core / "tables" / name).exists(), f"tables/{name} with --modules core"
        assert not (local_core / "work" / name.removesuffix(".tsv")).exists(), name
    for name in RESOURCES:
        assert not (local_core / "work" / name).exists(), f"{name} staged with --modules core"
    assert (local_core / "tables" / "operons.tsv").is_file()
    assert not _events(local_core, "function.run_flags")


# ---------------------------------------------------------------- example 2


@pytest.mark.needs_nextflow
def test_example_2_nextflow_gives_byte_identical_tables(
    sim: Any, db: Path, local: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    nf = tmp_path_factory.mktemp("NF") / "OUT"
    _run(sim, nf, db, "core,function", "--executor", "nextflow", "--profile", "farmshare_local")
    for rel in (*(f"tables/{t}" for t in TABLES), "genetable/tu_folds.parquet"):
        got = nf / rel
        assert got.is_file(), f"Nextflow run did not publish {rel}"
        assert got.read_bytes() == (local / rel).read_bytes(), rel
