"""The ``plots`` and ``report`` steps (GIAP_SPEC.md §6.2, §10.3; M9-02 behaviour 1, 2, 8).

``plots`` renders every registered plot (``giap.viz.base``) from a ``RunView`` of its inputs into
the ``plots`` directory and logs ``plots.rendered`` / ``plots.skipped`` events; ``report`` reads
that directory and writes ``report/report.html`` (``giap.report.build``). Both take every
published key as an optional input, so a subworkflow that did not run contributes nothing and the
plot skips itself. The run's modules come from ``cfg.modules`` (the steps have no resolved-config
input). ``render_plots`` is shared with ``giap report``. ``write_report`` builds the Plotly twins
(``giap.viz.interactive``) from the data TSVs under ``plots/data/`` when ``report.interactive``
is on, and logs ``report.size_exceeded {bytes}`` (WARNING) above ``REPORT_BUDGET`` (M9-14).
"""

from __future__ import annotations

import importlib
import logging
import os
import pkgutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from giap.config import GiapConfig
from giap.io.fasta import read_genome_fasta
from giap.io.plot_files import PlotFiles, read_plot_files, remove_plot_files, write_instance
from giap.io.plot_support import SUPPORT_DIR, write_support
from giap.io.run_view import PATHS, RunView
from giap.runner.local import PUBLISHED, PUBLISHED_FILES
from giap.steps.registry import InputPath, StepContext, step

#: the report size target (GIAP_SPEC.md §10.3: < 25 MB), static SVGs and twins together
REPORT_BUDGET = 25_000_000
_LOG = logging.getLogger("giap.report")

#: every key published under ``OUT/`` (and readable by ``RunView``), except the two outputs here
#: and ``run_manifest`` (written after the run; ``giap report`` reads it from ``OUT/``, M9-16)
ALL_PUBLISHED_KEYS = sorted(
    {*PUBLISHED, *PUBLISHED_FILES, *PATHS} - {"plots", "report", "preflight", "run_manifest"}
)
#: the support writers' inputs (M9-03, M9-18) and the RunView directories; ``report`` takes them
#: too, so REPORT can hand both steps the same [key, path] channel
_SUPPORT_INPUTS = ["genome", "features", "replicons", "genes", "step_events", "plots_support"]


@dataclass(frozen=True)
class PlotEvent:
    """``plots.rendered {plot, files, instances}`` or ``plots.skipped {plot, reason}``."""

    event: str
    data: dict[str, Any]


class PlotFailed(Exception):
    """A plot function raised (``giap report`` exit code 4; a ``StepError`` in ``plots``)."""

    def __init__(self, pid: str, cause: BaseException) -> None:
        super().__init__(f"plot {pid} failed: {type(cause).__name__}: {cause}")
        self.pid = pid


def load_plot_modules() -> None:
    """Import every ``giap.viz`` module, so each plot registers."""
    import giap.viz

    for info in pkgutil.walk_packages(giap.viz.__path__, prefix="giap.viz."):
        importlib.import_module(info.name)


def render_plots(
    view: RunView, cfg: GiapConfig, plots_dir: Path, only: Sequence[str] | None
) -> list[PlotEvent]:
    """Render the registered plots (``only``: normalised IDs; ``None``: all) in ID order into
    ``plots_dir``. Each plot's older files are removed first; other files are left alone."""
    from giap.viz import base, style

    style.apply()
    plots_dir.mkdir(parents=True, exist_ok=True)
    wanted = None if only is None else set(only)
    events: list[PlotEvent] = []
    for pid, spec in base.registered().items():
        if wanted is not None and pid not in wanted:
            continue
        remove_plot_files(plots_dir, pid)
        try:
            result = spec.fn(view, cfg)
        except Exception as exc:
            raise PlotFailed(pid, exc) from exc
        if isinstance(result, base.Skip):
            events.append(PlotEvent("plots.skipped", {"plot": pid, "reason": result.reason}))
            continue
        files: list[str] = []
        for inst in result:
            written = write_instance(
                inst, pid, spec.name, plots_dir, cfg.report.formats, cfg.report.dpi
            )
            files += [p.relative_to(plots_dir).as_posix() for p in written]
        data = {"plot": pid, "files": files, "instances": len(result)}
        events.append(PlotEvent("plots.rendered", data))
    return events


def skipped_plots(view: RunView) -> dict[str, str]:
    """Plot ID -> declared skip reason, for every registered plot that skips on ``view``, and
    ``not_registered`` for a report section's plot that this GIAP version does not have."""
    from giap.report.build import SECTIONS
    from giap.viz import base

    reasons = {pid: base.skip_reason(pid, view) for pid in base.registered()}
    absent = {p: "not_registered" for s in SECTIONS for p in s.plots if p not in reasons}
    return {pid: r for pid, r in reasons.items() if r is not None} | absent


def twin_meta(view: RunView, cfg: GiapConfig) -> dict[str, Any]:
    """The run context of the twins (M9-15 lead decision 5): the ``stats.alpha`` levels and the
    replicon lengths in FASTA order."""
    reps = (view.preflight() or {}).get("replicons", [])
    return {"alpha": cfg.stats.alpha.model_dump(),
            "replicons": {str(r["id"]): int(r["length"]) for r in reps}}  # fmt: skip


#: identifier columns of plot data TSVs, read as text so ``01`` / ``0001`` keep their spelling
TWIN_ID_DTYPES = dict.fromkeys(
    ("contig", "locus_tag", "pair_id", "label", "category", "product", "name"), str
)


def twins(
    plots_dir: Path, files: Mapping[str, PlotFiles], meta: Mapping[str, Any] | None = None
) -> dict[str, tuple[str, str]]:
    """Instance stem -> ``(div id, figure JSON)`` for every instance whose plot has a twin; the
    div id is ``plotly-<pid>[-<suffix>]`` (for an unregistered plot, the stem after the ID).
    The twin's title is the static plot's finding title (its SVG ``<title>``), else the stem.
    A twin that raises is not caught."""
    import pandas as pd

    from giap.io.plot_files import svg_title
    from giap.viz import base
    from giap.viz.interactive import figure_json, has_twin

    out: dict[str, tuple[str, str]] = {}
    for pid, pf in files.items():
        spec = base.registered().get(pid)
        for inst in pf.instances:
            if inst.data is None or not has_twin(pid):
                continue
            data = pd.read_csv(plots_dir / inst.data, sep="\t", na_values=["NA"],
                               keep_default_na=False, dtype=TWIN_ID_DTYPES)  # fmt: skip
            title = (svg_title(inst.svg) if inst.svg else None) or inst.stem
            text = figure_json(pid, data, title, meta)
            suffix = inst.stem.removeprefix(f"{pid}_{spec.name}" if spec else pid).lstrip("_")
            if text is not None:
                out[inst.stem] = (f"plotly-{pid}-{suffix}" if suffix else f"plotly-{pid}", text)
    return out


def write_report(
    view: RunView,
    plots_dir: Path | None,
    report_dir: Path,
    cfg: GiapConfig,
    log: logging.Logger = _LOG,
) -> Path:
    """Write ``report_dir/report.html`` from the plot files under ``plots_dir`` and the tables of
    ``view``; returns its path. Logs ``report.size_exceeded`` when it is over budget."""
    html = _report_html(view, plots_dir, cfg)
    report_dir.mkdir(parents=True, exist_ok=True)
    page = report_dir / "report.html"
    page.write_text(html, encoding="utf-8")
    _check_size(page, log)
    return page


def rerender_report(outdir: Path, cfg: GiapConfig) -> Path:
    """``OUT/report/report.html`` from ``OUT/`` and the existing ``OUT/plots`` (no plot is drawn);
    atomic replace through ``report/.report.html.tmp``; returns the path. Shared by the
    ``giap run`` tail (after ``run_manifest.json``) and ``giap report`` (M9-16b)."""
    load_plot_modules()
    plots_dir = outdir / "plots" if (outdir / "plots").is_dir() else None
    html = _report_html(RunView.from_outdir(outdir), plots_dir, cfg)
    report_dir = outdir / "report"
    report_dir.mkdir(parents=True, exist_ok=True)
    page, tmp = report_dir / "report.html", report_dir / ".report.html.tmp"
    try:
        tmp.write_text(html, encoding="utf-8")
        os.replace(tmp, page)  # replaces a symlinked page, never writes through it
    finally:
        tmp.unlink(missing_ok=True)
    _check_size(page, _LOG)
    _LOG.info("report.rerendered %s", page)
    return page


def _report_html(view: RunView, plots_dir: Path | None, cfg: GiapConfig) -> str:
    from giap.report.build import build_report

    files: Mapping[str, PlotFiles] = {} if plots_dir is None else read_plot_files(plots_dir)
    live = None
    if plots_dir is not None and cfg.report.interactive:
        live = twins(plots_dir, files, twin_meta(view, cfg))
    return build_report(view, files, skipped_plots(view), cfg, live)


def _check_size(page: Path, log: logging.Logger) -> None:
    size = page.stat().st_size
    if size > REPORT_BUDGET:
        data = {"bytes": size}
        log.warning("report.size_exceeded %d bytes > %d", size, REPORT_BUDGET,
                    extra={"event": "report.size_exceeded", "data": data})  # fmt: skip


def _view(ctx: StepContext) -> RunView:
    load_plot_modules()
    return RunView.from_inputs(ctx.inputs).with_modules(ctx.cfg.modules)


def _first(value: InputPath) -> Path:
    return value if isinstance(value, Path) else value[0]


def write_plot_support(inputs: Mapping[str, InputPath], plots_dir: Path) -> list[Path]:
    """``plots_dir/support/``: ``replicon_landmarks`` from ``replicons`` (and ``features``, absent
    = no oriC call) and ``gc_skew`` when ``genome`` is given too (M9-03); ``variant_windows``
    from ``variants`` and ``genes`` (M9-18). Only the ``plots`` step calls this; ``giap report``
    never does (ADR-0043)."""
    from giap.viz.support import WINDOW_BP, gc_skew, replicon_landmarks, variant_windows

    written: list[Path] = []
    if "replicons" in inputs:
        reps = pd.read_parquet(_first(inputs["replicons"]))
        feats = pd.read_parquet(_first(inputs["features"])) if "features" in inputs else None
        written.append(
            write_support("replicon_landmarks", replicon_landmarks(reps, feats), plots_dir)
        )
        if "genome" in inputs:
            genome = {r.id: r.sequence for r in read_genome_fasta(_first(inputs["genome"]))}
            written.append(write_support("gc_skew", gc_skew(genome, reps, WINDOW_BP), plots_dir))
    variants = RunView.from_inputs(inputs).table("variants") if "genes" in inputs else None
    if variants is not None:
        genes = pd.read_parquet(_first(inputs["genes"]), columns=["locus_tag", "aa"])
        written.append(
            write_support("variant_windows", variant_windows(variants, genes), plots_dir)
        )
    return written


@step("plots", inputs=["preflight"], optional_inputs=[*ALL_PUBLISHED_KEYS, *_SUPPORT_INPUTS],
      outputs=["plots"])  # fmt: skip
def plots_step(ctx: StepContext) -> None:
    """Tables -> ``plots/``: ``Vnn_*.{svg,png,pdf}`` and ``data/Vnn_*.tsv`` (behaviour 1), after
    the support tables ``support/*.tsv[.gz]`` from the run inputs (M9-03)."""
    out = ctx.outputs["plots"]
    inputs: dict[str, InputPath] = dict(ctx.inputs)
    if write_plot_support(inputs, out):  # before rendering: V1/V2 read them (M9-03)
        inputs["plots_support"] = out / SUPPORT_DIR
    load_plot_modules()
    view = RunView.from_inputs(inputs).with_modules(ctx.cfg.modules)
    for e in render_plots(view, ctx.cfg, out, None):
        ctx.log.info("%s %s", e.event, e.data["plot"], extra={"event": e.event, "data": e.data})


def _plots_dir(value: InputPath) -> Path | None:
    path = value if isinstance(value, Path) else value[0]
    return path if path.is_dir() else None


@step("report", inputs=["plots", "preflight"],
      optional_inputs=[*ALL_PUBLISHED_KEYS, *_SUPPORT_INPUTS], outputs=["report"])  # fmt: skip
def report_step(ctx: StepContext) -> None:
    """Tables + ``plots/`` -> ``report/report.html`` (behaviour 2)."""
    plots = _plots_dir(ctx.inputs["plots"])
    write_report(_view(ctx), plots, ctx.outputs["report"], ctx.cfg, ctx.log)
