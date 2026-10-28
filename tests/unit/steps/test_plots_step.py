"""The ``plots`` step and ``render_plots`` (M9-02 behaviour 1, acceptance "registry run with two
fake plots (one renders, one skips): files, events, skip reasons; ``--plots`` subset; stale-file
removal"; test idea "for any subset of registered plots, ``render_plots`` writes files only for
that subset"). Render and structure checks only, no ``--mpl`` baselines (M9-25).

Interface assumed beyond the spec text (M9-02 test author):
- ``gippa.steps.plots.render_plots(view, cfg, plots_dir, only)`` returns ``PlotEvent`` objects in
  plot-ID order, one per plot considered; a ``PlotEvent`` has ``event`` (``"plots.rendered"`` or
  ``"plots.skipped"``) and ``data`` (``{plot, files, instances}`` or ``{plot, reason}``);
  ``plot`` is the two-digit ID (``"V90"``), ``files`` the written paths relative to ``plots_dir``
  (no absolute path, §6.6), ``instances`` the number of instances; ``only`` holds normalised IDs
  (``None``: every registered plot).
- The step logs each event through ``ctx.log`` with ``extra={"event": ..., "data": ...}``
  (``gippa.io.events``), like the other steps.
- The ``plots`` output is a directory: ``Vnn_*.<fmt>`` and ``data/Vnn_*.tsv``.
- The module-off skip follows the run's ``cfg.modules`` in the step (the step has no
  ``resolved_config`` input), so a run with ``--modules core`` records ``module_off:<m>``.
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from types import ModuleType

import pytest
from helpers.report_fixtures import (
    cfg,
    register_fake,
    run_out,
    swap_registry,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.io.run_view import RunView
from gippa.runner.local import PUBLISHED, PUBLISHED_FILES
from gippa.steps import load_step_modules
from gippa.steps import registry as reg


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _step_inputs(out: Path, *, islands: bool) -> dict[str, Path]:
    inputs = {"preflight": out / "preflight" / "preflight.json"}
    if islands:
        inputs["islands"] = out / "tables" / "islands.tsv"
    return inputs


def _events(caplog: pytest.LogCaptureFixture) -> list[tuple[str, dict[str, object]]]:
    return [
        (r.event, dict(r.data))  # type: ignore[attr-defined]
        for r in caplog.records
        if isinstance(getattr(r, "event", None), str) and r.event.startswith("plots.")  # type: ignore[attr-defined]
    ]


def _files(plots: Path) -> set[str]:
    return {p.relative_to(plots).as_posix() for p in plots.rglob("*") if p.is_file()}


# ---------------------------------------------------------------- registration


def test_plots_and_report_steps_are_registered() -> None:
    load_step_modules()
    plots, report = reg.get("plots"), reg.get("report")
    assert plots.kind == "py" and report.kind == "py"
    assert plots.module == "core" and report.module == "core"
    assert plots.inputs == ("preflight",)
    assert plots.outputs == ("plots",)
    assert set(report.inputs) == {"plots", "preflight"}
    assert report.outputs == ("report",)


def test_plots_step_takes_every_published_key_as_optional_input() -> None:
    load_step_modules()
    opt = set(reg.get("plots").optional_inputs)
    published = (set(PUBLISHED) | set(PUBLISHED_FILES)) - {"plots", "report"}
    assert published <= opt, sorted(published - opt)
    # the support writers' inputs (M9-03, M9-18) and the M9-01 RunView directories
    assert {"events", "genome", "features", "replicons", "genes"} <= opt
    assert {"step_events", "plots_support"} <= opt
    assert "plots" not in opt and "report" not in opt


def test_report_step_takes_the_published_keys_and_events() -> None:
    load_step_modules()
    opt = set(reg.get("report").optional_inputs)
    published = (set(PUBLISHED) | set(PUBLISHED_FILES)) - {"plots", "report"}
    assert published <= opt, sorted(published - opt)
    assert "events" in opt
    assert "report" not in opt


# ---------------------------------------------------------------- the step


def test_one_plot_renders_one_skips_files_and_events(
    fake_registry: ModuleType, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    register_fake("V90", name="fake", suffixes=(None, "b"))
    register_fake("V91", name="other", needs=["pairs"])
    out = run_out(tmp_path, islands=1)
    plots = tmp_path / "work" / "plots"
    caplog.set_level(logging.INFO, logger="gippa.steps")
    reg.run_step("plots", cfg(), _step_inputs(out, islands=True), {"plots": plots})

    assert plots.is_dir()
    assert _files(plots) == {
        "V90_fake.svg",
        "V90_fake_b.svg",
        "data/V90_fake.tsv",
        "data/V90_fake_b.tsv",
    }
    events = _events(caplog)
    assert [e for e, _ in events] == ["plots.rendered", "plots.skipped"]
    rendered, skipped = events[0][1], events[1][1]
    assert rendered["plot"] == "V90"
    assert rendered["instances"] == 2
    files = [str(f) for f in rendered["files"]]  # type: ignore[attr-defined]
    assert not any(Path(f).is_absolute() for f in files), files
    assert {Path(f).name for f in files} == {
        "V90_fake.svg",
        "V90_fake_b.svg",
        "V90_fake.tsv",
        "V90_fake_b.tsv",
    }
    assert skipped == {"plot": "V91", "reason": "input_absent:pairs"}


def test_formats_follow_report_formats(fake_registry: ModuleType, tmp_path: Path) -> None:
    register_fake("V90", name="fake")
    out = run_out(tmp_path)
    plots = tmp_path / "work" / "plots"
    c = cfg("report.formats=[svg, png, pdf]")
    reg.run_step("plots", c, _step_inputs(out, islands=False), {"plots": plots})
    assert _files(plots) == {"V90_fake.svg", "V90_fake.png", "V90_fake.pdf", "data/V90_fake.tsv"}


def test_module_off_follows_the_run_modules(
    fake_registry: ModuleType, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    calls: list[str] = []
    register_fake("V14", module="function", calls=calls)
    out = run_out(tmp_path)
    plots = tmp_path / "work" / "plots"
    caplog.set_level(logging.INFO, logger="gippa.steps")
    reg.run_step("plots", cfg("modules=[core]"), _step_inputs(out, islands=False), {"plots": plots})
    assert calls == []
    assert _events(caplog) == [("plots.skipped", {"plot": "V14", "reason": "module_off:function"})]
    assert not any(p.name.startswith("V14_") for p in plots.rglob("*"))


def test_no_registered_plot_still_writes_the_plots_directory(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = run_out(tmp_path)
    plots = tmp_path / "work" / "plots"
    reg.run_step("plots", cfg(), _step_inputs(out, islands=False), {"plots": plots})
    assert plots.is_dir()


def test_a_raising_plot_fails_the_step(fake_registry: ModuleType, tmp_path: Path) -> None:
    register_fake("V90", raises=True)
    out = run_out(tmp_path)
    with pytest.raises(reg.StepError, match="V90"):
        reg.run_step(
            "plots", cfg(), _step_inputs(out, islands=False), {"plots": tmp_path / "work" / "plots"}
        )


def test_a_rerun_where_the_plot_now_skips_leaves_no_stale_file(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    register_fake("V90", name="fake", needs=["islands"])
    out = run_out(tmp_path, islands=1)
    plots = tmp_path / "work" / "plots"
    reg.run_step("plots", cfg(), _step_inputs(out, islands=True), {"plots": plots})
    assert (plots / "V90_fake.svg").is_file()
    reg.run_step("plots", cfg(), _step_inputs(out, islands=False), {"plots": plots})
    assert not any(p.name.startswith("V90_") for p in plots.rglob("*"))


# ---------------------------------------------------------------- render_plots


def test_render_plots_subset_and_event_order(fake_registry: ModuleType, tmp_path: Path) -> None:
    from gippa.steps.plots import render_plots

    for pid in ("V93", "V90", "V92", "V91"):
        register_fake(pid, name=f"p{pid}")
    view = RunView.from_outdir(run_out(tmp_path))
    plots = tmp_path / "plots"
    events = render_plots(view, cfg(), plots, ["V92", "V90"])
    assert [e.data["plot"] for e in events] == ["V90", "V92"]
    assert {e.event for e in events} == {"plots.rendered"}
    assert {p.name for p in plots.glob("*.svg")} == {"V90_pV90.svg", "V92_pV92.svg"}

    events = render_plots(view, cfg(), tmp_path / "all", None)
    assert [e.data["plot"] for e in events] == ["V90", "V91", "V92", "V93"]


def test_render_plots_leaves_other_plot_files_alone(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    from gippa.steps.plots import render_plots

    register_fake("V90", name="fake")
    register_fake("V91", name="other")
    plots = tmp_path / "plots"
    (plots / "data").mkdir(parents=True)
    (plots / "V91_other.svg").write_text("<svg/>")
    (plots / "legacy.txt").write_text("keep")
    render_plots(RunView.from_outdir(run_out(tmp_path)), cfg(), plots, ["V90"])
    assert (plots / "V91_other.svg").read_text() == "<svg/>"
    assert (plots / "legacy.txt").read_text() == "keep"
    assert (plots / "V90_fake.svg").is_file()
    assert (plots / "data" / "V90_fake.tsv").is_file()


def test_render_plots_twice_gives_identical_data_tsvs(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    from gippa.steps.plots import render_plots

    register_fake("V90", name="fake", suffixes=(None, "x"))
    view = RunView.from_outdir(run_out(tmp_path))
    for d in ("a", "b"):
        render_plots(view, cfg(), tmp_path / d, None)
    for name in ("V90_fake.tsv", "V90_fake_x.tsv"):
        a = (tmp_path / "a" / "data" / name).read_bytes()
        assert a == (tmp_path / "b" / "data" / name).read_bytes()
        assert a.startswith(b"x\ty\n")


_POOL = ("V90", "V91", "V92", "V93", "V94")


@settings(
    max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(subset=st.sets(st.sampled_from(_POOL)))
def test_property_render_plots_writes_only_the_subset(
    fake_registry: ModuleType, tmp_path: Path, subset: set[str]
) -> None:
    from gippa.steps.plots import render_plots
    from gippa.viz import base

    if not base.registered():
        for pid in _POOL:
            register_fake(pid, name="p")
    view = RunView.from_outdir(run_out(tmp_path))
    with tempfile.TemporaryDirectory(dir=tmp_path) as d:
        plots = Path(d) / "plots"
        events = render_plots(view, cfg(), plots, sorted(subset))
        written = {p.name.split("_")[0] for p in plots.rglob("*") if p.is_file()}
        assert written == subset
        assert [e.data["plot"] for e in events] == sorted(subset)
