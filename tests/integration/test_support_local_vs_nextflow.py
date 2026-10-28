"""M9-17 follow-up: ``OUT/plots/support/`` is the same under both executors (GIPPA_SPEC.md §6.4
executor equality; docs/specs/M9-17-plot-support-analysis.md lead decisions 5-6).

Inputs as tests/integration/test_m8_06_function_run.py: gippasim SIM100K with M3-19 proteomics,
every precomputed layer and the deterministic synthetic ``db_xrefs`` (so the function step has
testable sets), ``$GIPPA_DB_DIR`` = the synthetic function DB merged with the M3-17 test DB, seed 7,
``--modules core,function``. One run with ``--executor local``, one with ``--executor nextflow
--profile farmshare_local`` (real Nextflow on the prebuilt GIPPA env; no tool runs, as in
test_m2_21_local_vs_nextflow.py; stub files cannot be compared).

On SIM100K the determinants step fits no model (too_few_rows, test_determinants_run.py), so it
writes no ``determinants_support``: the Nextflow run must succeed with that optional output absent.

* the local run publishes support from at least one analysis step (else the comparison says
  nothing about the M9-17 wiring);
* the two runs publish the same set of ``plots/support/`` paths;
* each of those files is byte-identical between the executors.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

pytestmark = [pytest.mark.slow, pytest.mark.integration]

SEED = 7
#: the support files the three M9-17 analysis steps write (src/gippa/schemas/datapackage.json)
ANALYSIS = {"gsea_running", "determinants_shap", "determinants_pd", "null_determinants",
            "null_regions", "null_islands"}  # fmt: skip


def _xrefs_for(i: int) -> list[str]:
    """As test_m8_06_function_run._xrefs_for: 3 COGs, 12 GO leaves, 5 KOs; 1 CDS in 6 none."""
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
    import gippasim
    from gippasim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        gippasim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    res = gippasim.simulate(cfg, gippasim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")
    _add_xrefs(Path(res.bakta_json))
    return res


@pytest.fixture(scope="module")
def db(sim: Any, tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.function_db import ITEMS, PATHS, write_function_db
    from helpers.mapping_steps import write_test_db

    from gippasim.layout import layout_from_simulation

    root = tmp_path_factory.mktemp("function_db") / "db"
    write_function_db(root, layout_from_simulation(sim).xrefs)
    write_test_db(root)
    doc = json.loads((root / "db_manifest.json").read_text(encoding="utf-8"))
    doc |= {n: {"version": "synthetic", "path": PATHS[n]} for n in ITEMS}
    (root / "db_manifest.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return root


def _run(sim: Any, out: Path, db: Path, *extra: str) -> None:
    from gippa.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides)]
    assert sim.abundance is not None, "the fixture writes abundance.tsv"
    args += ["--abundance", str(sim.abundance)]
    args += ["--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", "core,function", "--seed", str(SEED), *extra]
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIPPA_DB_DIR", str(db))
        mp.setenv("NXF_ANSI_LOG", "false")
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output


def _support(out: Path) -> dict[str, bytes]:
    base = out / "plots" / "support"
    if not base.is_dir():
        return {}
    return {p.relative_to(out).as_posix(): p.read_bytes() for p in sorted(base.rglob("*"))
            if p.is_file()}  # fmt: skip


def _stem(rel: str) -> str:
    return Path(rel).name.split(".")[0]


@pytest.fixture(scope="module")
def local(sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory) -> dict[str, bytes]:
    out = tmp_path_factory.mktemp("LOCAL") / "OUT"
    _run(sim, out, db, "--executor", "local")
    return _support(out)


@pytest.fixture(scope="module")
def flow(sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory) -> dict[str, bytes]:
    out = tmp_path_factory.mktemp("NF") / "OUT"
    _run(sim, out, db, "--executor", "nextflow", "--profile", "farmshare_local")
    return _support(out)


def test_local_run_publishes_analysis_support(local: dict[str, bytes]) -> None:
    got = {_stem(r) for r in local}
    assert got & ANALYSIS, f"no analysis-step support under the local plots/support/: {sorted(got)}"
    assert not got & {"determinants_shap", "determinants_pd", "null_determinants"}, (
        "sim100k fits no determinants model, so no determinants support is expected"
    )


@pytest.mark.needs_nextflow
def test_nextflow_publishes_the_same_support_paths(
    local: dict[str, bytes], flow: dict[str, bytes]
) -> None:
    only_local, only_flow = sorted(set(local) - set(flow)), sorted(set(flow) - set(local))
    assert sorted(flow) == sorted(local), f"only local: {only_local}; only Nextflow: {only_flow}"


@pytest.mark.needs_nextflow
def test_nextflow_support_files_byte_identical(
    local: dict[str, bytes], flow: dict[str, bytes]
) -> None:
    differ = sorted(r for r in set(local) & set(flow) if local[r] != flow[r])
    assert not differ, f"plots/support files differ between the executors: {differ}"
    assert set(local) <= set(flow), "see test_nextflow_publishes_the_same_support_paths"
