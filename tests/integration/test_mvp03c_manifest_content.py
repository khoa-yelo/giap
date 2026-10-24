"""MVP-03c content blocks on real runs under both executors (docs/specs/MVP-03-manifest-finalize.md
acceptance 03c, D-C rows marked c, "Goal": executor-identical stable views; GIAP_SPEC.md §6.4,
§6.7).

The ``test_mvp03_manifest_executors`` setup (giapsim SIM100K with M3-19 proteomics and every
precomputed layer, peptides only, ``--modules core,mobile,bgc --seed 7``, ``--executor local``
and ``--executor nextflow --profile farmshare_local``). The database directory is the test DB of
``helpers.mapping_steps`` (provisional per-name form, never the real ``$GIAP_DB_DIR``); the
``retrieved`` dates are pinned by the unit tests (tests/unit/manifest/test_content_*.py). Run
through Slurm (``scripts/slurm/test_branch.sbatch``), never in the session.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import math
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest
from click.testing import CliRunner
from helpers.manifest_data import validate

pytestmark = [pytest.mark.needs_nextflow, pytest.mark.slow]

SEED = 7
BLOCKS = ("sensitivity", "context", "operons", "mapping", "calibration", "databases")


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim
    from giapsim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        giapsim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    return giapsim.simulate(cfg, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")


def _run(sim: Any, out: Path, *extra: str) -> None:
    from giap.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides)]
    args += ["--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", "core,mobile,bgc", "--seed", str(SEED), *extra]
    res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output


@pytest.fixture(scope="module")
def outs(sim: Any, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    mp = pytest.MonkeyPatch()
    try:
        from helpers.mapping_steps import write_test_db

        mp.setenv("GIAP_DB_DIR", str(write_test_db(tmp_path_factory.mktemp("test_db_mvp03c"))))
        mp.setenv("NXF_ANSI_LOG", "false")
        local = tmp_path_factory.mktemp("LOCAL") / "OUT"
        nf = tmp_path_factory.mktemp("NF") / "OUT"
        _run(sim, local, "--executor", "local")
        _run(sim, nf, "--executor", "nextflow", "--profile", "farmshare_local")
    finally:
        mp.undo()
    return {"local": local, "nextflow": nf}


def _manifest(out: Path) -> dict[str, Any]:
    path = out / "run_manifest.json"
    assert path.is_file(), f"no {path}"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _events(out: Path) -> list[dict[str, Any]]:
    path = out / "logs" / "events.jsonl"
    return [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


EXECUTORS = ["local", "nextflow"]


@pytest.mark.parametrize("executor", EXECUTORS)
def test_content_blocks_present_and_valid(outs: dict[str, Path], executor: str) -> None:
    m = _manifest(outs[executor])
    validate(m)
    assert set(BLOCKS) <= set(m)
    assert "Infinity" not in (outs[executor] / "run_manifest.json").read_text(encoding="utf-8")


@pytest.mark.parametrize("executor", EXECUTORS)
def test_context_replicons_equal_the_diagnostics_file(outs: dict[str, Path], executor: str) -> None:
    out = outs[executor]
    files = sorted(out.rglob("island_diagnostics.parquet"))
    reps = _manifest(out)["context"]["replicons"]
    if not files:  # CONTEXT skipped (precomputed regions layer): nothing to describe
        assert reps == []
        return
    want = [{k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in r.items()}
            for r in pq.read_table(files[0]).to_pylist()]  # fmt: skip
    assert reps == want


@pytest.mark.parametrize("executor", EXECUTORS)
def test_oric_block_equals_the_oriter_events(outs: dict[str, Path], executor: str) -> None:
    out = outs[executor]
    data = [e["data"] for e in _events(out) if e["event"] == "context_stats.oriter.multiple_oric"]
    uniq = {json.dumps(d, sort_keys=True): d for d in data}
    want = sorted(uniq.values(), key=lambda d: str(d["contig"]))
    assert _manifest(out)["context"]["oric"] == want


@pytest.mark.parametrize("executor", EXECUTORS)
def test_sensitivity_block_counts_the_table(outs: dict[str, Path], executor: str) -> None:
    out = outs[executor]
    m = _manifest(out)
    table = out / "tables" / "sensitivity.tsv"
    n = 0
    if table.is_file():
        with table.open(newline="", encoding="utf-8") as fh:
            n = sum(1 for _ in csv.DictReader(fh, delimiter="\t"))
    assert m["sensitivity"]["n_rows"] == n
    assert m["sensitivity"]["level"] == m["config"]["output"]["sensitivity"]


@pytest.mark.parametrize("executor", EXECUTORS)
def test_operons_and_mapping_blocks_follow_their_events(
    outs: dict[str, Path], executor: str
) -> None:
    out = outs[executor]
    m = _manifest(out)
    names = {e["event"] for e in _events(out)}
    assert (m["operons"] is None) == ("operons.missed_operons" not in names)
    assert (m["mapping"] is None) == ("mapping.search_db_resolved" not in names)
    if m["operons"] is not None:
        for b in m["operons"]["bin_counts"]:
            for key in ("lower_bp", "upper_bp"):
                assert b[key] is None or math.isfinite(b[key])


@pytest.mark.parametrize("executor", EXECUTORS)
def test_calibration_block(outs: dict[str, Path], executor: str) -> None:
    m = _manifest(outs[executor])
    # all three calibration.yaml keys, calibrated since ADR-0022 (M6-12 C): no uncalibrated flag
    assert set(m["calibration"]) == {"operon_distance", "uniop_pair_threshold", "terminatornet_yes"}
    assert all(e["calibrated"] is True for e in m["calibration"].values())
    assert not [f for f in m["flags"] if f.startswith("uncalibrated:")]


@pytest.mark.parametrize("executor", EXECUTORS)
def test_databases_listed(outs: dict[str, Path], executor: str) -> None:
    dbs = _manifest(outs[executor])["databases"]
    assert dbs, "the test DB lists every requirement"
    assert [d["name"] for d in dbs] == sorted(d["name"] for d in dbs)
    assert all(d["retrieved"] is None for d in dbs)  # provisional form


def test_content_blocks_equal_across_executors(outs: dict[str, Path]) -> None:
    """§6.4 extended to the 03c blocks (part of ``stable_view`` equality)."""
    a, b = _manifest(outs["local"]), _manifest(outs["nextflow"])
    differ = [k for k in BLOCKS if a[k] != b[k]]
    assert not differ, f"local and Nextflow differ on {differ}"
