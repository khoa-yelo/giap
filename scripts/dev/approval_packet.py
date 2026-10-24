"""Write the M9 approval packet for Khoa (docs/specs/M9-25-snapshot-baselines.md item 5).

Usage: python scripts/dev/approval_packet.py [--dest DIR] [--baselines DIR] [--k12 OUT]

``DEST`` (default ``local/review/m9-approval``, not committed) gets: ``baseline/`` (a copy of
every baseline PNG) and ``deuteranopia/`` (their simulations); with ``--k12`` the same plots of
the M9-26 K-12 run under ``k12/plots/`` and ``k12/deuteranopia/``; ``checklist.tsv`` and
``score.json`` (the M9-24 checklist of every instance on each snapshot fixture and on K-12);
``report.html`` and ``compare.html`` built on the rich fixture under ``pages/``; ``index.html``,
which links all of it with relative paths. Regenerate the baselines first if the plots changed
(tests/data/snapshot/PROVENANCE.md): ``pytest tests/snapshot/test_plots.py -k
test_plot_snapshot --mpl-generate-path=tests/snapshot/baseline``.
"""

from __future__ import annotations

import argparse
import html
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / "src"), str(REPO / "tests"), str(REPO / "tests" / "sim")]

DEFAULT_DEST = REPO / "local" / "review" / "m9-approval"
DEFAULT_BASELINES = REPO / "tests" / "snapshot" / "baseline"
DPI = 100
HEADER = ["fixture", "plot", "item", "passed", "detail"]


def deuteranopia(src: Path, dest: Path) -> None:
    """``dest``: the simulated-deuteranopia PNG of ``src`` (alpha kept)."""
    from giap.io.plot_files import read_png, write_png
    from giap.viz.checklist import simulate_deuteranopia

    img = read_png(src)
    img[..., :3] = simulate_deuteranopia(img[..., :3])
    write_png(dest, img)


def check(results: dict[str, Any], plots: Path) -> list[Any]:
    """Write every instance of ``results`` (``pid -> instances | Skip``) into ``plots`` and run
    the checklist on it; returns the ``CheckResult``s."""
    from giap.io.plot_files import write_instance
    from giap.viz import base, checklist

    out = []
    for pid, result in results.items():
        if isinstance(result, base.Skip):
            continue
        for inst in result:
            files = write_instance(inst, pid, base.registered()[pid].name, plots, ["svg", "png"],
                                   DPI)  # fmt: skip
            out += checklist.check_instance(inst, pid, files)
    return out


def k12_results(out: Path) -> dict[str, Any]:
    from giap.config import resolve_config
    from giap.io.run_view import RunView
    from giap.steps.plots import load_plot_modules
    from giap.viz import base

    load_plot_modules()
    runs = (out / "resolved_config.yaml", out / "work" / "resolved_config.yaml")
    cfg = resolve_config(next((p for p in runs if p.is_file()), None))
    view = RunView.from_outdir(out)
    return {pid: spec.fn(view, cfg) for pid, spec in base.registered().items()}


def _cell(v: bool | None) -> str:
    return "NA" if v is None else str(v).lower()


def _imgs(dest: Path, sub: str, sims: str, title: str) -> str:
    rows = []
    for png in sorted((dest / sub).glob("*.png")):
        pair = [f"{sub}/{png.name}", f"{sims}/{png.name}"]
        cells = "".join(f'<td><a href="{p}"><img src="{p}" width="480"></a></td>' for p in pair)
        rows.append(f"<tr><th>{html.escape(png.stem)}</th>{cells}</tr>")
    head = "<tr><th>plot</th><th>as drawn</th><th>deuteranopia</th></tr>"
    return f"<h2>{title}</h2><table>{head}{''.join(rows)}</table>"


def write_index(dest: Path, score: dict[str, float], pages: list[Path], k12: bool) -> None:
    items = "".join(f"<li>{html.escape(f)}: {s:.2f}</li>" for f, s in score.items())
    links = "".join(f'<li><a href="{p.relative_to(dest).as_posix()}">{p.name}</a></li>'
                    for p in pages)  # fmt: skip
    body = [
        "<h1>M9 approval packet: plot baselines and report layout</h1>",
        f"<h2>Checklist score (M9-24)</h2><ul>{items}</ul>",
        '<p><a href="checklist.tsv">checklist.tsv</a> (every item of every instance)</p>',
        f"<h2>Pages (rich snapshot fixture)</h2><ul>{links}</ul>",
        _imgs(dest, "baseline", "deuteranopia", "Baselines (tests/snapshot/baseline, 100 dpi)"),
        _imgs(dest, "k12/plots", "k12/deuteranopia", "K-12 run (M9-26)")
        if k12
        else "<h2>K-12 run (M9-26)</h2><p>The K-12 run was not supplied (--k12).</p>",
    ]
    page = '<!doctype html><html><head><meta charset="utf-8"><title>M9 approval</title></head>'
    (dest / "index.html").write_text(page + "<body>" + "".join(body) + "</body></html>\n",
                                     encoding="utf-8")  # fmt: skip


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--dest", type=Path, default=DEFAULT_DEST)
    ap.add_argument("--baselines", type=Path, default=DEFAULT_BASELINES)
    ap.add_argument("--k12", type=Path, default=None, help="OUT/ of the M9-26 K-12 run")
    ns = ap.parse_args(argv)
    if ns.k12 is not None and not ns.k12.is_dir():
        print(f"approval_packet: --k12 {ns.k12} is not a directory", file=sys.stderr)
        return 2
    from helpers import snapshot_fixture as sf

    from giap.io.tsv import write_tsv
    from giap.viz import checklist

    dest: Path = ns.dest
    for png in sorted(ns.baselines.glob("*.png")):
        (dest / "baseline").mkdir(parents=True, exist_ok=True)
        shutil.copyfile(png, dest / "baseline" / png.name)
        deuteranopia(png, dest / "deuteranopia" / png.name)
    labelled = {f: sf.styled(lambda f=f: sf.raw_results(f)) for f in sf.FIXTURES}
    if ns.k12 is not None:
        labelled["k12"] = sf.styled(lambda: k12_results(ns.k12))
    rows, score = [], {}
    for label, results in labelled.items():
        plots = dest / "k12" / "plots" if label == "k12" else dest / "fixtures" / label
        found = sf.styled(lambda r=results, p=plots: check(r, p))
        score[label] = checklist.score(found)
        rows += [[label, r.plot, r.item, _cell(r.passed), r.detail] for r in found]
    for png in sorted((dest / "k12" / "plots").glob("*.png")):
        deuteranopia(png, dest / "k12" / "deuteranopia" / png.name)
    write_tsv(dest / "checklist.tsv", HEADER, rows)
    (dest / "score.json").write_text(json.dumps(score, indent=1) + "\n", encoding="utf-8")
    shutil.rmtree(dest / "pages", ignore_errors=True)
    pages = [sf.build_report(sf.RICH, dest / "pages" / "run"),
             sf.build_compare(sf.RICH, dest / "pages" / "cmp")]  # fmt: skip
    write_index(dest, score, pages, ns.k12 is not None)
    print(f"wrote {os.path.relpath(dest / 'index.html')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
