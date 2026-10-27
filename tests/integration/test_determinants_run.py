"""M8-10 through ``gippa run``: examples 1 and 2 of docs/specs/M8-10-determinants-step-and-run.md
(GIPPA_SPEC.md §6.2, §6.4, §7.8.3, §9.1, §9.5).

Example 1: ``gippa run --executor local`` on gippasim SIM100K (tests/helpers/m4_real.py, the MVP-01
harness) writes ``tables/determinants_{features,model,attribution}.tsv``; the model table has
univariate rows and exactly one ``too_few_rows`` row per model (sim100k is below the 100-row
minimum); every attribution row (if any gene is a model row) has an NA prediction; exit 0. The
three tables validate against their schemas and the caveats reach ``OUT/logs/events.jsonl``.

Both runs use ``--modules core,mobile,bgc`` and the M4 test database (``m4_chain.test_db``).

Example 2: the same ``gippa run`` with ``--executor nextflow`` gives byte-identical tables. Real
Nextflow on the prebuilt GIPPA env (``--profile farmshare_local``, as
test_m2_21_local_vs_nextflow.py: every annotation layer is precomputed, so no tool runs). The
spec says ``-profile test``; ``farmshare_local`` is the profile the existing local-vs-Nextflow
test uses for a real (non-stub) run on FarmShare, so it is used here too.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import determinants_step_data as s

SEED = 20261014
SETS = ("determinants.permutations=19", "stats.bootstrap=49")


def _sim(root: Path) -> Any:
    import gippasim
    from gippasim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        gippasim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    return gippasim.simulate(cfg, gippasim.SIM100K_SEED, root / "sim")


def _run(sim: Any, out: Path, *extra: str) -> None:
    from gippa.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides),
            "--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml),
            "--outdir", str(out), "--seed", str(SEED), "--modules", "core,mobile,bgc",
            *extra]  # fmt: skip
    for x in SETS:
        args += ["--set", x]
    res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output


@pytest.fixture(scope="module")
def local_out(tmp_path_factory: pytest.TempPathFactory) -> Path:
    s.mod()  # a missing step fails before the run
    root = tmp_path_factory.mktemp("m810_local")
    mp = pytest.MonkeyPatch()
    try:
        mp.setenv("GIPPA_DB_DIR", str(s.mc.test_db(root / "db")))
        out = root / "OUT"
        _run(_sim(root), out, "--executor", "local")
    finally:
        mp.undo()
    return out


def test_example_1_three_tables(local_out: Path) -> None:
    for key, rel in s.PUBLISHED.items():
        assert (local_out / rel).is_file(), f"{rel} not published"
        s.validate(local_out / rel, s.SCHEMAS[key])


def test_example_1_model_table(local_out: Path) -> None:
    tbl = s.read(local_out / s.PUBLISHED["determinants_model"])
    assert (tbl["model"] == "univariate").any(), "no univariate rows"
    for m in s.MODELS:
        r = s.rows(tbl, model=m)
        assert len(r) == 1 and "too_few_rows" in s.flags(r["flags"].iloc[0]), r.to_dict("records")
    assert s.sort_key(tbl) == sorted(s.sort_key(tbl))


def test_example_1_attribution_na_predictions(local_out: Path) -> None:
    path = local_out / s.PUBLISHED["determinants_attribution"]
    assert s.header(path) == list(s.ATTRIBUTION_COLUMNS)
    a = s.read(path)
    assert (a["predicted_log2"] == "NA").all()
    assert all("too_few_rows" in s.flags(f) for f in a["flags"])


def test_example_1_features_rows(local_out: Path) -> None:
    feats = s.read(local_out / s.PUBLISHED["determinants_features"])
    assert len(feats) > 50, f"{len(feats)} feature rows on sim100k"


def test_example_1_caveats_event(local_out: Path) -> None:
    lines = (local_out / "logs" / "events.jsonl").read_text(encoding="utf-8").splitlines()
    ev = [e for e in map(json.loads, lines) if e["event"] == s.CAVEATS_EVENT]
    assert len(ev) == 1 and ev[0]["step"] == "determinants", ev
    assert ev[0]["data"]["caveats"] == list(s.CAVEATS)


@pytest.mark.needs_nextflow
@pytest.mark.slow
def test_example_2_nextflow_byte_identical(
    local_out: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    root = tmp_path_factory.mktemp("m810_nf")
    mp = pytest.MonkeyPatch()
    try:
        mp.setenv("GIPPA_DB_DIR", str(s.mc.test_db(root / "db")))
        mp.setenv("NXF_ANSI_LOG", "false")
        out = root / "OUT"
        _run(_sim(root), out, "--executor", "nextflow", "--profile", "farmshare_local")
    finally:
        mp.undo()
    for rel in s.PUBLISHED.values():
        assert (out / rel).is_file(), f"{rel} not published by Nextflow"
        assert (out / rel).read_bytes() == (local_out / rel).read_bytes(), rel
