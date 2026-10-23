"""Run the M9-24 style checklist over every registered plot rendered on one ``OUT/`` tree.

Shared by ``tests/snapshot/test_checklist_all.py`` (the M9-25 frozen fixtures) and
``tests/integration/test_checklist_sim100k.py`` (a real sim100k run). Each plot is called through
its registry wrapper (so declared skips skip), each instance is written with ``write_instance``
into ``plots_dir`` and checked with ``giap.viz.checklist.check_instance``.
"""

from __future__ import annotations

from pathlib import Path

from giap.config import resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView


def sweep(out: Path, plots_dir: Path) -> tuple[dict[str, list[str]], dict[str, str]]:
    """``(failures, skipped)``: stem -> failing gated items with details; pid -> skip reason."""
    from giap.steps.plots import load_plot_modules
    from giap.viz import base, checklist, style

    load_plot_modules()
    style.apply()
    runs = (out / "resolved_config.yaml", out / "work" / "resolved_config.yaml")
    cfg = resolve_config(next((p for p in runs if p.is_file()), None))
    view = RunView.from_outdir(out)
    failures: dict[str, list[str]] = {}
    skipped: dict[str, str] = {}
    for pid, spec in base.registered().items():
        result = spec.fn(view, cfg)
        if isinstance(result, base.Skip):
            skipped[pid] = result.reason
            continue
        for inst in result:
            files = write_instance(inst, pid, spec.name, plots_dir, ["svg", "png"], 100)
            bad = [
                f"{r.item}: {r.detail}"
                for r in checklist.check_instance(inst, pid, files)
                if r.passed is False
            ]
            if bad:
                failures[files[0].stem] = bad
    return failures, skipped
