"""``gippa _step placement`` writes the ``search_db: auto`` resolution event (M3-17 decision 14;
docs/backlog.md "m3/mapping-nf: write the search_db auto-resolution event from the Nextflow
``gippa _step`` path too").

Under Nextflow every py step is one ``gippa _step <name>`` task in a fresh interpreter whose
working directory is the task directory; ``--executor local`` already writes the event to
``OUT/logs/events.jsonl`` (test_mapping_steps_review.py::test_d14_...). On the ``_step`` path the
event goes to ``logs/events.jsonl`` under the task's working directory (the same ``logs/`` the
command already uses for ``logs/<step>.log`` on a step error), one §6.7 JSON object per line
(``ts``, ``step``, ``level``, ``event``, ``data``).

Set-up: a 50-kb same-strain gippasim run, every peptide's ``search_db`` set to ``auto``; the chain
up to ``map_tolerant`` (plus ``parse_bakta``) runs in-process through ``run_step``; ``placement``
then runs as a subprocess (fresh interpreter, as under Nextflow).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml
from helpers.mapping_steps import drive, fresh_steps, roots_of, simulate_with_proteomics

from gippa.config import resolve_config

REPO = Path(__file__).resolve().parents[2]
_MAIN = "import sys; from gippa.cli import main; main(sys.argv[1:])"
SEED = 20260917


def _env() -> dict[str, str]:
    env = dict(os.environ)
    paths = [str(REPO / "src"), str(REPO / "tests" / "sim")]
    env["PYTHONPATH"] = os.pathsep.join([*paths, env.get("PYTHONPATH", "")])
    return env


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    from gippasim.proteomics import ProteomicsConfig

    tmp = tmp_path_factory.mktemp("mapping_nf_events")
    run = simulate_with_proteomics(tmp, SEED, proteomics=ProteomicsConfig(search_db="same"))
    df = pd.read_parquet(run.peptides)
    df["search_db"] = pd.Series(["auto"] * len(df), dtype="string")
    auto = tmp / "peptides.auto.parquet"
    df.to_parquet(auto, index=False)
    roots = roots_of(run) | {"peptides": auto}
    cfg = resolve_config(run.config_yaml, [f"seed={SEED}"])
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, tmp / "work", roots, targets=("map_tolerant",))
        spec = reg.get("placement")
    resolved = tmp / "resolved_config.yaml"
    resolved.write_text(yaml.safe_dump(cfg.model_dump(mode="json"), sort_keys=True))
    ins = {k: avail[k] for k in spec.inputs}
    return {"ins": ins, "outputs": list(spec.outputs), "config": resolved}


def _step_placement(chain: dict[str, Any], task: Path) -> subprocess.CompletedProcess[str]:
    task.mkdir(parents=True)
    argv = ["_step", "placement", "--config", str(chain["config"])]
    for k, p in chain["ins"].items():
        argv += ["--in", f"{k}={p}"]
    for k in chain["outputs"]:
        argv += ["--out", f"{k}={k}.parquet"]
    return subprocess.run(
        [sys.executable, "-c", _MAIN, *argv],
        cwd=task, env=_env(), capture_output=True, text=True, timeout=600, check=False,
    )  # fmt: skip


def _search_db_events(path: Path) -> list[dict[str, Any]]:
    found = []
    for ln in path.read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        ev = json.loads(ln)
        data = ev.get("data") if isinstance(ev.get("data"), dict) else {}
        if "search_db" in str(ev.get("event", "")) or "search_db" in data:
            found.append(ev)
    return found


def test_gippa_step_placement_writes_the_search_db_event(
    chain: dict[str, Any], tmp_path: Path
) -> None:
    task = tmp_path / "task"
    res = _step_placement(chain, task)
    assert res.returncode == 0, res.stdout + res.stderr
    for k in chain["outputs"]:
        assert (task / f"{k}.parquet").is_file(), f"gippa _step placement wrote no {k}"
    events = task / "logs" / "events.jsonl"
    assert events.is_file(), (
        "gippa _step placement wrote no logs/events.jsonl in its working directory; files: "
        f"{sorted(str(p.relative_to(task)) for p in task.rglob('*'))}"
    )
    evs = _search_db_events(events)
    assert evs, events.read_text(encoding="utf-8")[-2000:]
    ev = evs[0]
    assert set(ev) >= {"ts", "step", "level", "event", "data"}, ev
    assert ev["step"] == "placement", ev
    assert "same" in {str(v) for v in ev["data"].values()}, ev
