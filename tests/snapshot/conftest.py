"""Keeps ``pytest tests/snapshot --mpl -n 2`` inside ``ci/gates.yaml`` ``budgets_s.snapshot``
(#325 review S5) without changing what any test renders or checks.

* ``helpers.snapshot_fixture.SESSION_CACHE`` points at a directory every xdist worker of this
  session shares (the parent of the workers' base temp dirs; the base temp dir without xdist):
  the report pages' plot files are rendered once per session and the raw plot results once per
  process.
* ``$GIPPA_SNAPSHOT_CASES`` (the same directory) holds the case lists, so the subprocesses of
  ``test_plots_sensitivity.py`` read them instead of rendering both fixtures, and draw only the
  plots they save (still in their own interpreter, with their own ``style`` module).
* In a snapshot-only xdist session the long tests run first, longest first.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from helpers import snapshot_fixture as sf


@pytest.hookimpl(trylast=True)
def pytest_configure(config: pytest.Config) -> None:
    factory = getattr(config, "_tmp_path_factory", None)
    if factory is None:
        return
    base = factory.getbasetemp()
    root = base.parent if hasattr(config, "workerinput") else base
    sf.SESSION_CACHE = root / "gippa-snapshot-cache"
    os.environ[sf.CASES_ENV] = str(sf.SESSION_CACHE)
    # xdist controller: hand out one test at a time per worker (plus the one it holds next), so
    # the long tests that pytest_collection_modifyitems puts first go out longest first
    controller = not hasattr(config, "workerinput")
    unset = getattr(config.option, "maxschedchunk", 0) is None
    if controller and unset and _snapshot_only(config):
        config.option.maxschedchunk = 1


def _snapshot_only(config: pytest.Config) -> bool:
    here = Path(__file__).resolve().parent
    paths = [Path(str(a).split("::")[0]).resolve() for a in config.args]
    return bool(paths) and all(p == here or here in p.parents for p in paths)


#: Approximate call seconds of the long snapshot tests (``--durations`` at -n 2 on FarmShare,
#: with this cache); every other test counts as ``_CHEAP``. Only the run order uses them.
_COST = {
    "test_snapshot_byte_identity.py::test_every_plot_byte_identical_on_the_snapshot_fixture[sim_plots]": 80.0,  # noqa: E501
    "test_checklist_all.py::test_every_registered_plot_passes_every_gated_item[sim_plots]": 34.0,
    "test_report_structure.py::test_page_structure_matches_the_snapshot[sim_plots/report.html]": 27.0,  # noqa: E501
    "test_plots_sensitivity.py::test_detected_colour_change_breaks_v01_v02_v16": 16.0,
    "test_plots_sensitivity.py::test_unchanged_colour_matches_the_baselines": 15.0,
    "test_snapshot_byte_identity.py::test_every_plot_byte_identical_on_the_snapshot_fixture[sim_empty]": 14.0,  # noqa: E501
    "test_checklist_all.py::test_every_registered_plot_passes_every_gated_item[sim_empty]": 3.0,
    "test_report_structure.py::test_page_structure_matches_the_snapshot[sim_empty/report.html]": 2.5,  # noqa: E501
}
_CHEAP = 0.3


def _cost(item: pytest.Item) -> float:
    return _COST.get(f"{item.path.name}::{item.name}", _CHEAP)


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Under xdist, order a snapshot-only session longest first, each long test followed by one
    cheap test (the one its worker holds next), the other cheap tests last: with one test
    scheduled at a time (``pytest_configure``) that is longest-first greedy balancing. The order
    never changes what a test checks."""
    workers = int(os.environ.get("PYTEST_XDIST_WORKER_COUNT", "0") or 0)
    here = Path(__file__).resolve().parent
    if workers < 2 or not items or any(here not in item.path.parents for item in items):
        return
    long = sorted((i for i in items if _cost(i) > _CHEAP), key=_cost, reverse=True)
    cheap = iter([i for i in items if _cost(i) <= _CHEAP])
    order: list[pytest.Item] = []
    for item in long:
        order.append(item)
        order += [c for c in (next(cheap, None),) if c is not None]
    order += list(cheap)
    items[:] = order
