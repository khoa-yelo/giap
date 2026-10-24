"""Runtime and peak RSS of the ``plots`` and ``report`` steps on a 5-Mb giapsim run
(docs/specs/M9-26-scale-and-determinism.md work item 1; GIAP_SPEC.md §10.1, §11.2, §12).

The run is the M4 chain on ``build_m4_fixture(size="sim4500")`` (≈ 4.6 Mb, 4,500 genes,
≈ 45,000 peptide placements, an iBAQ abundance table), prepared once and not timed. Each step
then runs through ``giap.steps.registry.run_step`` in a fresh interpreter, so its wall time and
peak RSS (``VmHWM``) are its own. The measured entries live in ``benchmarks/baseline.json`` as
two flat entries, ``plots_5mb_sim4500`` and ``report_5mb_sim4500`` (found by their ``benchmark``
field; M9-26 lead decision 4). After the timed runs every registered plot is rendered once more
and the M9-24 checklist item "≤ 20,000 drawn elements" (per axes, per V1 track) must pass
on each instance (M9-26 lead decision 2, gap 5). Budgets (lead decision 4): 120 s for
``plots`` and 30 s for ``report``; the spec's 10-minute line stays the K-12 one. The timed span
starts after imports, so interpreter start-up is excluded.

Run with ``pytest benchmarks/test_plots_bench.py --benchmark-only -p no:xdist`` (needs
``tests/sim`` on PYTHONPATH for giapsim).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests" / "helpers"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests"))
import m4_chain as mc
from helpers.viz_landscape import write_preflight

SEED = 20261023
#: per-step wall-time budgets in seconds (M9-26 lead decision 4)
BUDGETS_S = {"plots": 120.0, "report": 30.0}
#: the 5-Mb run's placement lower bound (sim4500 gives ≈ 42k; lead decision 4)
MIN_PLACEMENTS = 35_000
MAX_ELEMENTS = 20_000
ELEMENT_ITEM = "≤ 20,000 drawn elements"

#: one step in a fresh interpreter; prints ``{"secs", "maxrss_mb"}`` as the last line
_CHILD = """
import json, sys, time
from pathlib import Path
from giap.config.model import GiapConfig
from giap.steps import load_step_modules
from giap.steps.registry import run_step
name, seed = sys.argv[1], int(sys.argv[2])
ins, outs = json.loads(sys.argv[3]), json.loads(sys.argv[4])
load_step_modules()
cfg = GiapConfig.model_validate({"seed": seed})
t0 = time.perf_counter()
run_step(name, cfg, {k: Path(v) for k, v in ins.items()}, {k: Path(v) for k, v in outs.items()})
secs = time.perf_counter() - t0
# VmHWM, not ru_maxrss: Linux keeps ru_maxrss across fork + exec (the parent pytest's RSS)
hwm = [ln for ln in Path("/proc/self/status").read_text().splitlines() if ln.startswith("VmHWM:")]
rss = int(hwm[0].split()[1]) / 1024.0
print(json.dumps({"secs": round(secs, 3), "maxrss_mb": round(rss, 1)}))
"""


def _inputs(work: Path) -> dict[str, Path]:
    """The step inputs of the work directory: every published key it holds, the ``mapping``
    directory and the support writers' work tables."""
    from giap.runner.local import PUBLISHED_FILES

    ins = {"preflight": work / "preflight" / "preflight.json", "mapping": work / "mapping",
           "replicons": work / "replicons.parquet", "features": work / "features.parquet",
           "genes": work / "genes.parquet"}  # fmt: skip
    ins |= {k: work / rel for k, rel in PUBLISHED_FILES.items() if (work / rel).is_file()}
    return ins


Prepared = tuple[Path, dict[str, Path], dict[str, int]]


@pytest.fixture(scope="module")
def prepared(tmp_path_factory: pytest.TempPathFactory) -> Prepared:
    import pandas as pd

    root = tmp_path_factory.mktemp("plots_bench")
    work = mc.build(root / "work", seed=SEED, size="sim4500",
                    evidence=mc.evidence_config(abundance_table="ibaq"))  # fmt: skip
    mc.run_m4(work, mc.config(SEED), mc.test_db(root / "db"))
    reps = pd.read_parquet(work / "replicons.parquet")
    write_preflight(work, [(str(r.contig), int(r.length), str(r.topology), str(r.kind))
                           for r in reps.itertuples(index=False)])  # fmt: skip
    size = {
        "n_genes": len(pd.read_parquet(work / "genes.parquet")),
        "n_placements": len(pd.read_csv(work / "mapping" / "peptide_hits.tsv", sep="\t")),
        "genome_bp": int(reps["length"].sum()),
    }
    return work, _inputs(work), size


def _child(name: str, ins: dict[str, Path], outs: dict[str, Path]) -> dict[str, float]:
    args = [json.dumps({k: str(v) for k, v in ins.items()}),
            json.dumps({k: str(v) for k, v in outs.items()})]  # fmt: skip
    cmd = [sys.executable, "-c", _CHILD, name, str(SEED), *args]
    res = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ), check=False)
    assert res.returncode == 0, res.stderr[-4000:]
    rec: dict[str, float] = json.loads(res.stdout.strip().splitlines()[-1])
    return rec


def _element_failures(ins: dict[str, Path], plots: Path, dest: Path) -> list[str]:
    """Every registered plot rendered on the run, checked for the M9-24 element item."""
    from giap.io.plot_files import write_instance
    from giap.io.run_view import RunView
    from giap.steps.plots import load_plot_modules
    from giap.viz import base, checklist, style

    load_plot_modules()
    style.apply()
    cfg = mc.config(SEED)
    view = RunView.from_inputs({**ins, "plots_support": plots / "support"})
    view = view.with_modules(cfg.modules)
    bad: list[str] = []
    for pid, spec in base.registered().items():
        result = spec.fn(view, cfg)
        if isinstance(result, base.Skip):
            continue
        for inst in result:
            files = write_instance(inst, pid, spec.name, dest, ["svg"], 72)
            checks = checklist.check_instance(inst, pid, files)
            bad += [f"{files[0].stem}: {r.detail}" for r in checks
                    if r.item == ELEMENT_ITEM and r.passed is False]  # fmt: skip
    return bad


def test_plots_and_report_5mb_budget(benchmark: Any, prepared: Prepared, tmp_path: Path) -> None:
    work, ins, size = prepared
    plots, report = work / "plots", work / "report"

    def run() -> dict[str, Any]:
        rep_ins = {**ins, "plots": plots}
        return {"plots": _child("plots", ins, {"plots": plots}),
                "report": _child("report", rep_ins, {"report": report})}  # fmt: skip

    rec = benchmark.pedantic(run, rounds=1, iterations=1)
    rec |= size | {"report_bytes": (report / "report.html").stat().st_size}
    benchmark.extra_info.update(rec)
    print("PLOTS_BENCH " + json.dumps(rec, sort_keys=True))
    assert 4_000 <= size["n_genes"] <= 5_500, size
    assert MIN_PLACEMENTS <= size["n_placements"] <= 60_000, size
    assert _element_failures(ins, plots, tmp_path / "plots") == []
    for step in ("plots", "report"):
        budget = BUDGETS_S[step]
        assert rec[step]["secs"] <= budget, f"{step} took {rec[step]['secs']:.1f} s > {budget:g} s"
