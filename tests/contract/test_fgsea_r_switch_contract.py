"""M8-17b on real R (docs/specs/M8-17-r-fgsea-engine.md §9, §11 T7, A2, A6).

Slow and ``requires_r``: runs only on Slurm (``sbatch scripts/slurm/test_branch.sbatch <worktree>
-m requires_r``, at most 2 CPUs) where the ``gippa-fgsea`` env is built.

- ``gsea_table`` on the M8-06 step fixture (four recorded collections, seed 7, ``n_perm_simple``
  200) with the real ``RFgseaEngine`` writes exactly the committed replay table
  ``tests/data/expected/fgsea_replay/sim_main_enrichment_gsea.tsv``: the recording, the replay and
  the live engine agree byte for byte.
- The ``function`` step with the real engine, twice: every ``enrichment_gsea.tsv`` row has
  ``method = fgsea_multilevel`` and both runs are byte-identical (A2; the Python stand-in for the
  ``function_rfgsea.nf.test`` real run, which needs #187's ``function.nf``).
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest
from helpers.fgsea_r_data import real_prefix
from helpers.fgsea_switch_data import EXPECTED_SIM_MAIN, run_sim_gsea, sim_cfg, sim_world
from helpers.function_db import FUNCTION_RESOURCES
from helpers.function_step_data import write_tu_folds

PREFIX = real_prefix()

pytestmark = [
    pytest.mark.slow,
    pytest.mark.requires_r,
    pytest.mark.skipif(
        PREFIX is None,
        reason="the gippa-fgsea env is not built (sbatch scripts/slurm/build_envs.sbatch fgsea); "
        "run on Slurm with scripts/slurm/test_branch.sbatch -m requires_r",
    ),
]


def test_live_engine_reproduces_the_replay_table(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import r_fgsea_engine, resolve_launcher
    from gippa.io.function_tables import write_enrichment_gsea_tsv

    assert PREFIX is not None
    eng = r_fgsea_engine(resolve_launcher({"GIPPA_FGSEA_PREFIX": str(PREFIX)}), tmp_root=tmp_path)
    w = sim_world(tmp_path / "world")
    df = run_sim_gsea(w, eng)
    out = tmp_path / "enrichment_gsea.tsv"
    write_enrichment_gsea_tsv(out, df)
    assert out.read_bytes() == EXPECTED_SIM_MAIN.read_bytes()
    assert eng.provenance is not None and eng.provenance.fgsea == "1.36.2"


def test_function_step_on_real_r_is_byte_identical(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    assert PREFIX is not None
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(PREFIX))
    monkeypatch.chdir(tmp_path)
    load_step_modules()
    w = sim_world(tmp_path / "world")
    folds = write_tu_folds(tmp_path / "tu_folds.parquet", {})
    ins = {"gene_table": w.gene_table, "replicons": w.replicons_path, "tu_folds": folds}
    ins |= dict(w.db.inputs(FUNCTION_RESOURCES))
    texts = []
    for run in ("a", "b"):
        outs = {
            k: tmp_path / run / f"{k}.tsv"
            for k in ("category_summary", "enrichment_detection", "enrichment_gsea")
        }
        reg.run_step("function", sim_cfg(), ins, outs)
        texts.append({k: p.read_bytes() for k, p in outs.items()})
        df = pd.read_csv(outs["enrichment_gsea"], sep="\t", dtype=str, keep_default_na=False)
        assert len(df) > 0
        assert set(df["method"]) == {"fgsea_multilevel"}
    assert texts[0] == texts[1]
    assert not list(Path(os.getcwd()).glob("gippa-fgsea-*"))
