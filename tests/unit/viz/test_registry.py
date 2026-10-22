"""The plot registry and its skip wrapper (M9-01 behaviour 6, examples 3, 4 and 8; index
decision 4).

Interface assumed beyond the spec text (the ``giap.steps.registry`` pattern):
- ``giap.viz.base._REGISTRY`` is the module-level ``dict[str, PlotSpec]``; the ``plots`` fixture
  swaps it for an empty dict so fake plots never collide with real ones.
- ``@plot(...)`` returns the *wrapped* function (same ``PlotFn`` signature) that applies the skip
  checks before calling the plot; ``PlotSpec`` has the fields ``pid``, ``name``, ``needs``,
  ``optional``, ``module``, ``layers`` and ``fn`` (the wrapped function).
- Registering an already registered ID raises ``ValueError``.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest
from helpers.viz_out_tree import (
    ISLAND_ROW,
    make_out,
    write_events,
    write_resolved_config,
    write_schema_table,
)

from giap.config import GiapConfig, resolve_config
from giap.io.run_view import RunView
from giap.viz import base
from giap.viz.base import Instance, Skip


@pytest.fixture
def plots(monkeypatch: pytest.MonkeyPatch) -> Iterator[ModuleType]:
    monkeypatch.setattr(base, "_REGISTRY", {})
    yield base


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


def _register(reg: ModuleType, pid: str, name: str = "fake", **kw: object) -> list[str]:
    calls: list[str] = []

    @reg.plot(pid, name, **kw)
    def _fn(view: RunView, cfg: GiapConfig) -> list[Instance] | Skip:
        calls.append(pid)
        return []

    return calls


def _flags_out(tmp_path: Path, *, modules: list[str] | None = None, islands: bool = False) -> Path:
    out = make_out(tmp_path)
    write_events(
        out / "logs" / "steps" / "context_layers.events.jsonl",
        [
            {
                "event": "context_layers.layer_missing",
                "data": {"flag": "layer_missing:genomad", "layer": "genomad"},
                "step": "context_layers",
                "level": "WARNING",
                "ts": "2026-10-22T00:00:00+00:00",
            }
        ],
    )
    write_events(
        out / "logs" / "events.jsonl",
        [{"event": "function.run_flags", "data": {"flags": ["layer_missing:kegg"]}}],
    )
    if modules is not None:
        write_resolved_config(out / "resolved_config.yaml", modules)
    if islands:
        write_schema_table(out / "tables" / "islands.tsv", "islands", [ISLAND_ROW])
    return out


# ---------------------------------------------------------------- IDs


def test_example_3_normalise_ids(plots: ModuleType) -> None:
    for pid in ("V18", "V01", "V03"):
        _register(plots, pid, f"p{pid}")
    assert plots.normalise_ids("V1, v18,V03") == ["V01", "V03", "V18"]


def test_example_4_deferred_id_is_unknown(plots: ModuleType) -> None:
    for pid in ("V01", "V03", "V18"):
        _register(plots, pid, f"p{pid}")
    with pytest.raises(ValueError, match="unknown plot ID V04"):
        plots.normalise_ids("V04")
    with pytest.raises(ValueError, match=r"V0?4"):
        plots.normalise_ids("V01,V4")


def test_registered_is_sorted_by_id(plots: ModuleType) -> None:
    for pid in ("V18", "V01", "V03"):
        _register(plots, pid, f"p{pid}")
    reg = plots.registered()
    assert list(reg) == ["V01", "V03", "V18"]
    assert reg["V03"].pid == "V03"
    assert reg["V03"].name == "pV03"


def test_ids_are_unique(plots: ModuleType) -> None:
    _register(plots, "V90", "one")
    with pytest.raises(ValueError):
        _register(plots, "V90", "two")


def test_spec_records_its_declaration(plots: ModuleType) -> None:
    _register(
        plots, "V90", needs=["islands"], optional=["genes_tsv"], module="mobile", layers=["genomad"]
    )
    spec = plots.registered()["V90"]
    assert list(spec.needs) == ["islands"]
    assert list(spec.optional) == ["genes_tsv"]
    assert spec.module == "mobile"
    assert list(spec.layers) == ["genomad"]


# ---------------------------------------------------------------- skips


def test_example_8_layer_missing_skips_without_calling_the_plot(
    plots: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    calls: list[str] = []

    @plots.plot("V90", "fake", needs=[], layers=["genomad"])
    def fake(view: RunView, cfg: GiapConfig) -> list[Instance] | Skip:
        calls.append("called")
        return []

    view = RunView.from_outdir(_flags_out(tmp_path))
    assert fake(view, cfg) == Skip("layer_missing:genomad")
    assert plots.registered()["V90"].fn(view, cfg) == Skip("layer_missing:genomad")
    assert calls == []


def test_module_off_skips(plots: ModuleType, cfg: GiapConfig, tmp_path: Path) -> None:
    calls = _register(plots, "V90", module="mobile", needs=["islands"])
    view = RunView.from_outdir(_flags_out(tmp_path, modules=["core"], islands=True))
    assert plots.registered()["V90"].fn(view, cfg) == Skip("module_off:mobile")
    assert calls == []


def test_input_absent_skips(plots: ModuleType, cfg: GiapConfig, tmp_path: Path) -> None:
    calls = _register(plots, "V90", needs=["islands"])
    view = RunView.from_outdir(_flags_out(tmp_path))
    result = plots.registered()["V90"].fn(view, cfg)
    assert isinstance(result, Skip)
    assert result.reason.startswith("input_absent:")
    assert "islands" in result.reason
    assert calls == []


def test_skip_order_module_then_layer_then_input(
    plots: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    all_three = _register(plots, "V90", module="mobile", layers=["genomad"], needs=["islands"])
    layer_and_input = _register(plots, "V91", layers=["kegg"], needs=["islands"])
    view = RunView.from_outdir(_flags_out(tmp_path, modules=["core"]))
    reg = plots.registered()
    assert reg["V90"].fn(view, cfg) == Skip("module_off:mobile")
    assert reg["V91"].fn(view, cfg) == Skip("layer_missing:kegg")
    assert all_three == [] and layer_and_input == []


def test_core_is_never_off_and_present_inputs_run_the_plot(
    plots: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    calls = _register(plots, "V90", module="core", layers=["bakta"], needs=["islands"])
    view = RunView.from_outdir(_flags_out(tmp_path, modules=["core"], islands=True))
    assert plots.registered()["V90"].fn(view, cfg) == []
    assert calls == ["V90"]


def test_absent_optional_input_does_not_skip(
    plots: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    seen: list[pd.DataFrame | None] = []

    @plots.plot("V90", "fake", needs=["islands"], optional=["genes_tsv"])
    def fake(view: RunView, cfg: GiapConfig) -> list[Instance] | Skip:
        seen.append(view.table("genes_tsv"))
        return []

    view = RunView.from_outdir(_flags_out(tmp_path, islands=True))
    assert fake(view, cfg) == []
    assert seen == [None]


def test_a_raising_plot_is_not_swallowed(
    plots: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    @plots.plot("V90", "fake", needs=[])
    def fake(view: RunView, cfg: GiapConfig) -> list[Instance] | Skip:
        raise RuntimeError("boom")

    view = RunView.from_outdir(_flags_out(tmp_path))
    with pytest.raises(RuntimeError, match="boom"):
        fake(view, cfg)
