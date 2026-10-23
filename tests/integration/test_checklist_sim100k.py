"""M9-24 sweep on a real run: every registered plot rendered from a local-executor ``giap run`` of
sim100k (giapsim seed 20261014, ``--modules core``, as tests/integration/test_report_run.py)
passes every gated item of the §14.4 style checklist; plots that skip carry a declared reason.

This is the evidence before M9-25 freezes the snapshot fixtures (tests/snapshot/
test_checklist_all.py sweeps those).
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.checklist_sweep import sweep

pytestmark = [pytest.mark.slow]

SEED = 7
REASONS = ("module_off:", "layer_missing:", "input_absent:")


@pytest.fixture(scope="module")
def out(tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.mapping_steps import write_test_db

    import giapsim
    from giap.cli import cli
    from giapsim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        giapsim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    sim: Any = giapsim.simulate(
        cfg, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim"
    )
    db = Path(write_test_db(tmp_path_factory.mktemp("test_db_m9_24")))
    res_dir = tmp_path_factory.mktemp("LOCAL") / "OUT"
    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides)]
    args += ["--abundance", str(sim.abundance), "--annotation", str(sim.bakta_json)]
    args += ["--config", str(sim.config_yaml), "--outdir", str(res_dir), "--modules", "core"]
    args += ["--seed", str(SEED), "--executor", "local"]
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIAP_DB_DIR", str(db))
        mp.setenv("NXF_ANSI_LOG", "false")
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output
    return res_dir


def test_sim100k_plots_pass_every_gated_item(out: Path, tmp_path: Path) -> None:
    failures, skipped = sweep(out, tmp_path / "plots")
    assert failures == {}
    assert all(r.startswith(REASONS) for r in skipped.values()), skipped
