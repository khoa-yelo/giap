"""M4-15 lead decision 15: ``gippa run --executor local`` stages ``cog2020_def`` from the setup
manifest (like ``contaminants``), so ``cog_category`` is set under the local executor too; and the
local gene table keeps string NAs null (decision 11).

Inputs: the 50-kb gippasim run of tests/integration/test_mapping_local.py with a test
``$GIPPA_DB_DIR`` (never the real one) whose manifest also lists ``cog2020_def`` and holds the
fixture ``cog/cog-20.def.tab``. The simulated Bakta JSON has no COG xrefs, so a copy gives the first
two CDS ``COG:COG0593`` / ``COG:COG0002``. The local work directory names files by key
(``work/gene_table`` is the ``detectability`` output).
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner
from helpers.mapping_steps import write_test_db

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import gene_table_fixture as gt

SEED = 7


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SimConfig, simulate

    out = tmp_path_factory.mktemp("sim50k_cog")
    cfg = SimConfig(proteomics=ProteomicsConfig(), write_precomputed=True)
    return simulate(cfg, SEED, out / "sim")


def _db(root: Path) -> Path:
    db = write_test_db(root / "db")
    gt.write_cog_def(db / "cog" / "cog-20.def.tab")
    man = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    man["cog2020_def"] = {"version": "2020", "path": "cog/cog-20.def.tab"}
    (db / "db_manifest.json").write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")
    return db


def _annotation(sim: Any, root: Path) -> tuple[Path, dict[str, str]]:
    doc = json.loads(Path(sim.bakta_json).read_text(encoding="utf-8"))
    cogs = iter(["COG0593", "COG0002"])
    given: dict[str, str] = {}
    for f in doc["features"]:
        if f.get("type") == "cds" and f.get("locus") and len(given) < 2:
            cog = next(cogs)
            f["db_xrefs"] = [*f.get("db_xrefs", []), f"COG:{cog}"]
            given[str(f["locus"])] = cog
    assert len(given) == 2, "simulated Bakta JSON has fewer than two CDS"
    path = root / "bakta_cog.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path, given


def test_local_executor_stages_cog_table(
    sim: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.cli import cli

    monkeypatch.setenv("GIPPA_DB_DIR", str(_db(tmp_path)))
    ann, given = _annotation(sim, tmp_path)
    out = tmp_path / "OUT"
    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides),
            "--annotation", str(ann), "--config", str(sim.config_yaml), "--outdir", str(out),
            "--modules", "core", "--seed", str(SEED), "--executor", "local"]  # fmt: skip
    r = CliRunner().invoke(cli, args)
    assert r.exit_code == 0, r.output
    table = pq.read_table(out / "work" / "gene_table")
    by = dict(zip(table.column("locus_tag").to_pylist(),
                  table.column("cog_category").to_pylist(), strict=True))  # fmt: skip
    want = {"COG0593": "L", "COG0002": "KT"}
    for locus, cog in given.items():
        assert by[locus] == want[cog], (locus, by[locus])
    n = table.num_rows
    assert table.column("cog_category").null_count == n - 2
    assert table.column("island_id").null_count == n


def test_local_executor_stages_abundance_and_preflight(
    sim: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """Lead decision 20(b): with ``--abundance`` the local run stages ``abundance`` and
    ``preflight`` (into ``work/preflight/``), the scale is read without a
    ``contract.preflight.abundance_scale`` issue, and the abundance columns are filled."""
    from gippa.cli import cli

    assert sim.abundance is not None and Path(sim.abundance).is_file()
    monkeypatch.setenv("GIPPA_DB_DIR", str(_db(tmp_path)))
    out = tmp_path / "OUT"
    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides),
            "--abundance", str(sim.abundance), "--annotation", str(sim.bakta_json),
            "--config", str(sim.config_yaml), "--outdir", str(out), "--modules", "core",
            "--seed", str(SEED), "--executor", "local"]  # fmt: skip
    with caplog.at_level(logging.WARNING):
        r = CliRunner().invoke(cli, args)
    assert r.exit_code == 0, r.output
    work = out / "work"
    staged = work / "preflight" / "preflight.json"
    assert staged.is_file()
    scale = json.loads(staged.read_text(encoding="utf-8"))["routing"]["abundance_scale"]
    assert scale in ("log2", "linear")
    assert (work / "abundance").is_file()
    msgs = [rec.getMessage() for rec in caplog.records]
    assert not [m for m in msgs if "contract.preflight" in m], msgs
    table = pq.read_table(work / "gene_table")
    types = set(table.column("abundance_type").to_pylist())
    assert types and types != {"none"}
    assert table.column("abundance_log2").null_count < table.num_rows
