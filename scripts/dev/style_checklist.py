"""Run the M9-24 style checklist (GIPPA_SPEC.md §14.4) over every registered plot of one run.

Usage: python scripts/dev/style_checklist.py OUT [NAME]

Renders each plot from ``OUT`` into ``local/review/<NAME>/plots`` (``NAME`` defaults to the
name of ``OUT``; the simulated-deuteranopia PNGs land beside the originals), prints one table per
plot and writes ``local/review/<NAME>/checklist.tsv`` (not committed).
"""

from __future__ import annotations

import sys
from pathlib import Path

from gippa.config import resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.io.tsv import write_tsv
from gippa.steps.plots import load_plot_modules
from gippa.viz import base, checklist, style

REPO = Path(__file__).resolve().parents[2]


def main(argv: list[str]) -> int:
    if not 1 <= len(argv) <= 2:
        print(__doc__, file=sys.stderr)
        return 2
    out = Path(argv[0])
    review = REPO / "local" / "review" / (argv[1] if len(argv) == 2 else out.name)
    load_plot_modules()
    style.apply()
    runs = (out / "resolved_config.yaml", out / "work" / "resolved_config.yaml")
    cfg = resolve_config(next((p for p in runs if p.is_file()), None))
    view = RunView.from_outdir(out)
    rows: list[list[object]] = []
    for pid, spec in base.registered().items():
        result = spec.fn(view, cfg)
        if isinstance(result, base.Skip):
            print(f"{pid} {spec.name}: skipped ({result.reason})\n")
            continue
        for inst in result:
            files = write_instance(inst, pid, spec.name, review / "plots", ["svg", "png"], 100)
            results = checklist.check_instance(inst, pid, files)
            print(f"{results[0].plot}  score {checklist.score(results):.2f}")
            for r in results:
                mark = {True: "pass", False: "FAIL", None: "----"}[r.passed]
                print(f"  {mark}  {r.item:<32} {r.detail}")
                rows.append([r.plot, r.item, r.passed, r.detail])
            print()
    write_tsv(review / "checklist.tsv", ["plot", "item", "passed", "detail"], rows)
    print(f"wrote {review / 'checklist.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
