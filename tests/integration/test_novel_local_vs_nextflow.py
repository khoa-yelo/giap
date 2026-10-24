"""The four Module 6 files through ``giap run`` under both executors (issue M7-08 example 2) and
TE7 on them (GIAP_SPEC.md §6.4 executor equality, §11.6 TE7).

Inputs: giapsim ``SIM100K`` with proteomics, planted novel ORFs and N-terminal extensions, and
every precomputed layer including the simulator's ``rbs_novel`` (tests/helpers/
novel_steps_fixture.py). As in tests/integration/test_m2_21_local_vs_nextflow.py, the Nextflow
run is real Nextflow with ``--profile farmshare_local``: every layer is precomputed, so no tool
runs and NOVEL runs only its py steps -- the precomputed ``promoters`` layer serves both rounds
and the precomputed ``rbs_novel`` layer replaces OSTIR_NOVEL + NOVEL_RBS_LAYER (§6.4: a non-null
``annotation.precomputed.<key>`` replaces its tool and parse step under both executors; the
spec's "-profile test, stubs" cannot give byte-identical tables, stub outputs are empty files).

* example 2: ``tables/novel_orfs.tsv``, ``novel_orfs.gff3``, ``novel_orfs.faa``,
  ``annotation_patch.gff3`` byte-identical between the executors;
* TE7: two local runs with the same seed give byte-identical four files.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from helpers import novel_steps_fixture as nf


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    return nf.simulate(tmp_path_factory.mktemp("m7_08_int") / "sim")


def _four(out: Path) -> dict[str, bytes]:
    got = {}
    for name in nf.TABLES:
        path = out / "tables" / name
        assert path.is_file(), f"{path} not published"
        got[name] = path.read_bytes()
    return got


def test_te7_two_local_runs_byte_identical(
    sim: Any, tmp_path_factory: pytest.TempPathFactory
) -> None:
    root = tmp_path_factory.mktemp("te7")
    outs = []
    with nf.db_env(root):
        for tag in ("a", "b"):
            code, log = nf.run_giap(sim, root / tag / "OUT")
            assert code == 0, log
            outs.append(_four(root / tag / "OUT"))
    assert outs[0] == outs[1]
    assert outs[0]["novel_orfs.tsv"].count(b"\n") > 1, "no novel ORF row on the planted genome"


@pytest.mark.needs_nextflow
@pytest.mark.slow
def test_example_2_local_and_nextflow_byte_identical(
    sim: Any, tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path_factory.mktemp("ex2")
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    with nf.db_env(root):
        code, log = nf.run_giap(sim, root / "LOCAL" / "OUT")
        assert code == 0, log
        code, log = nf.run_giap(sim, root / "NF" / "OUT", "--profile", "farmshare_local",
                                executor="nextflow")  # fmt: skip
        assert code == 0, log
    local, flow = _four(root / "LOCAL" / "OUT"), _four(root / "NF" / "OUT")
    differ = sorted(k for k in local if local[k] != flow[k])
    assert not differ, f"local and Nextflow differ: {differ}"


# ---------------------------------------------------------------- M7-11b acceptance criterion 5
#
# The class-FDR wiring (docs/specs/M7-11b-class-fdr-wiring.md criterion 5): a giapsim
# ``novel_orf_power(seed, "100k")`` run, whose database peptide table carries makedb six-frame
# ``rev_GIAPORF_...`` decoys, under both executors. The four tables and the published side table
# ``novel/novel_class_fdr.parquet`` must be the same bytes. The side table must show the picked
# regime, so the parity covers the wired class FDR and not only its ``decoys_not_makedb`` path.

SIDE = "novel/novel_class_fdr.parquet"


@pytest.fixture(scope="module")
def power_sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim
    import giapsim.scenarios as sc

    seed = giapsim.SIM100K_SEED
    root = tmp_path_factory.mktemp("m7_11b_power")
    return giapsim.simulate(sc.novel_orf_power(seed, "100k"), seed, root / "sim")


@pytest.mark.needs_nextflow
@pytest.mark.slow
def test_m7_11b_power_run_local_and_nextflow_byte_identical_with_side_table(
    power_sim: Any, tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    import pandas as pd

    root = tmp_path_factory.mktemp("m7_11b_ex5")
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    with nf.db_env(root):
        code, log = nf.run_giap(power_sim, root / "LOCAL" / "OUT")
        assert code == 0, log
        code, log = nf.run_giap(power_sim, root / "NF" / "OUT", "--profile", "farmshare_local",
                                executor="nextflow")  # fmt: skip
        assert code == 0, log
    lo, fl = root / "LOCAL" / "OUT", root / "NF" / "OUT"
    local, flow = _four(lo), _four(fl)
    for name, d in (("local", lo), ("nextflow", fl)):
        assert (d / SIDE).is_file(), f"{name}: {SIDE} not published"
    local[SIDE], flow[SIDE] = (lo / SIDE).read_bytes(), (fl / SIDE).read_bytes()
    side = pd.read_parquet(lo / SIDE)
    assert "picked_decoy" in set(side["regime"]), sorted(set(side["regime"]))
    assert local["novel_orfs.tsv"].count(b"\n") > 1, "no novel ORF row on the power genome"
    differ = sorted(k for k in local if local[k] != flow[k])
    assert not differ, f"local and Nextflow differ: {differ}"
