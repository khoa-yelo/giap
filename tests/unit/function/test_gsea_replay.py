"""M8-17b T7 and A6: ``gsea_table`` replays real fgsea 1.36.2 output without R
(docs/specs/M8-17-r-fgsea-engine.md §11 T7, A6; lead notes for M8-17b 2 and 3).

``helpers.fgsea_switch_data.ReplayEngine`` checks that every call ``gsea_table`` makes writes
exactly the recorded input TSVs (``write_call`` byte for byte) with the recorded options
(``command.txt``: ``--n-perm-simple 200`` and the seed), then returns ``read_results`` of the
recorded ``results.tsv``. The M8-06 step fixture is rebuilt as ``scripts/dev/record_fgsea.py``
did (seed 7, ``n_perm_simple`` 200) and restricted to the four recorded collections (lead
note 2).

- sim_main: one call; the table written by ``write_enrichment_gsea_tsv`` equals the committed
  ``tests/data/expected/fgsea_replay/sim_main_enrichment_gsea.tsv`` (the test-author's oracle,
  ``helpers.fgsea_switch_data.expected_sim_main``).
- sim_retest: no fixture GO term is significant, so the recording forced significance (each
  go_bp parent minus its lowest-p tested descendant). The replay forces the same significance by
  overriding the main-call p of exactly those four terms (GO:0990100 with GO:0990000, GO:0990101
  with GO:0990002); ``gsea_table``'s own §8 re-test construction must then send exactly the
  recorded sim_retest call, and ``redundant`` follows the recorded re-test p.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.fgsea_switch_data import (
    EXPECTED_SIM_MAIN,
    RECORDED,
    ReplayEngine,
    SimWorld,
    run_sim_gsea,
    sim_world,
)

from giap.stats.multitest import bh

FORCED = {"GO:0990100": 1e-4, "GO:0990000": 2e-4, "GO:0990101": 1e-4, "GO:0990002": 2e-4}
RETESTED = ("GO:0990100", "GO:0990101")


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> SimWorld:
    return sim_world(tmp_path_factory.mktemp("fgsea_replay"))


def test_expected_table_is_committed() -> None:
    assert EXPECTED_SIM_MAIN.is_file()
    text = EXPECTED_SIM_MAIN.read_text(encoding="utf-8")
    assert len(text.splitlines()) == 33  # header + the 32 recorded rows


def test_sim_main_replay_is_one_call(world: SimWorld) -> None:
    eng = ReplayEngine([RECORDED / "sim_main"])
    run_sim_gsea(world, eng)
    assert len(eng.calls) == 1


def test_sim_main_replay_equals_the_committed_table(world: SimWorld, tmp_path: Path) -> None:
    from giap.io.function_tables import write_enrichment_gsea_tsv

    df = run_sim_gsea(world, ReplayEngine([RECORDED / "sim_main"]))
    assert df.attrs["run_flags"] == []
    out = tmp_path / "enrichment_gsea.tsv"
    write_enrichment_gsea_tsv(out, df)
    assert out.read_bytes() == EXPECTED_SIM_MAIN.read_bytes()


def test_sim_main_replay_validates_against_the_schema(world: SimWorld, tmp_path: Path) -> None:
    from helpers.fgsea_switch_data import validate_gsea

    from giap.io.function_tables import write_enrichment_gsea_tsv

    out = tmp_path / "enrichment_gsea.tsv"
    write_enrichment_gsea_tsv(out, run_sim_gsea(world, ReplayEngine([RECORDED / "sim_main"])))
    validate_gsea(out)


def _forced() -> dict[tuple[str, str], dict[str, Any]]:
    return {("main:go_bp", t): {"p": p} for t, p in FORCED.items()}


def _recorded_retest_p() -> dict[str, float]:
    lines = (RECORDED / "sim_retest" / "results.tsv").read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    out = {}
    for line in lines[1:]:
        r = dict(zip(head, line.split("\t"), strict=True))
        out[r["pathway"]] = math.nan if r["pval"] == "NA" else float(r["pval"])
    return out


def test_sim_retest_replay_sends_the_recorded_retest_call(world: SimWorld) -> None:
    eng = ReplayEngine([RECORDED / "sim_main", RECORDED / "sim_retest"], override=_forced())
    run_sim_gsea(world, eng)
    assert len(eng.calls) == 2
    main, re = eng.calls
    assert [j.job_id for j in re.jobs] == [f"retest:go_bp:{t}" for t in RETESTED]
    assert re.seed == main.seed and re.tu_ids == main.tu_ids


def test_sim_retest_replay_redundancy_follows_the_recorded_p(world: SimWorld) -> None:
    eng = ReplayEngine([RECORDED / "sim_main", RECORDED / "sim_retest"], override=_forced())
    df = run_sim_gsea(world, eng)
    bp = df[df["collection"] == "go_bp"].reset_index(drop=True)
    p = bp["p_value"].to_numpy(float)
    for t, t_p in FORCED.items():
        assert bp.loc[bp["term"] == t, "p_value"].item() == t_p
    np.testing.assert_array_equal(bp["q_value"].to_numpy(float), bh(p))
    p_re = _recorded_retest_p()
    alpha = world.cfg.stats.alpha.enrichment
    for t in RETESTED:
        i = int(bp.index[bp["term"] == t][0])
        sub = p.copy()
        sub[i] = p_re[t]
        want = bool(bh(sub)[i] > alpha)
        got = bp.loc[i, "redundant"]
        assert not pd.isna(got) and bool(got) is want, (t, p_re[t])
    for t in set(bp["term"]) - set(RETESTED):
        assert bool(bp.loc[bp["term"] == t, "redundant"].item()) is False, t
