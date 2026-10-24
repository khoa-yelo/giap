"""The recorded runtimes of the M7-15 real-tool arm (docs/specs/M7-15-integration-hardening.md
item 4, example 5: "Promotech/OSTIR second-round runtimes recorded"; acceptance "a recorded run
(job IDs in PROGRESS.md)").

Pinned shape of ``benchmarks/baseline.json`` > ``novel_second_round_k12`` (from
``scripts/measure_runtimes.py`` on the run's ``logs/trace.tsv``; process names are its lowercased
last ``:`` component)::

    {"source": "scripts/measure_runtimes.py ...", "genome": "k12_mg1655",
     "peptides": "<the TE8 in-silico set>", "profile": "farmshare,apptainer",
     "slurm_job": <head job ID>, "date": "YYYY-MM-DD",
     "processes": {"promotech_novel": {"wall_s": <float > 0>, "peak_rss": "<trace text>",
                                       "node_class": "<wheat|barley|...>", "tasks": <int >= 1>},
                   "ostir_novel": {...}}}

The job ID must also appear in PROGRESS.md.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
KEY = "novel_second_round_k12"
DATE = re.compile(r"^20\d\d-\d\d-\d\d$")


def _entry() -> dict[str, object]:
    base = json.loads((REPO / "benchmarks" / "baseline.json").read_text(encoding="utf-8"))
    assert KEY in base, f"benchmarks/baseline.json has no {KEY} (the K-12 real-tool arm)"
    entry = base[KEY]
    assert isinstance(entry, dict)
    return entry


def test_second_round_runtimes_recorded() -> None:
    e = _entry()
    assert "measure_runtimes.py" in str(e["source"])
    assert e["genome"] == "k12_mg1655"
    assert "farmshare" in str(e["profile"]) and "apptainer" in str(e["profile"])
    assert isinstance(e["slurm_job"], int) and e["slurm_job"] > 0
    assert DATE.match(str(e["date"]))
    procs = e["processes"]
    assert isinstance(procs, dict)
    for name in ("promotech_novel", "ostir_novel"):
        p = procs[name]
        assert float(p["wall_s"]) > 0, name
        assert int(p["tasks"]) >= 1, name
        assert p["peak_rss"] and p["node_class"], name


def test_second_round_job_is_in_progress_md() -> None:
    job = str(_entry()["slurm_job"])
    text = (REPO / "PROGRESS.md").read_text(encoding="utf-8")
    assert re.search(rf"\b{job}\b", text), f"Slurm job {job} is not recorded in PROGRESS.md"
