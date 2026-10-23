"""The interactive report (M9-14 behaviour 1-3, 6; examples 1-3; degenerate rows "no twin
registered for any rendered plot" and "a twin raises"; acceptance "the JSON spec of the V17 twin
equals the data TSV values"; GIAP_SPEC.md §10.2 column "Interactive in report", §10.3
"Self-contained ... Plotly.js is inlined once", §5.1 ``--no-interactive``, §5.2
``report.interactive``).

Twins are test-only (lead decision 3), registered into a swapped ``giap.viz.interactive`` registry
over fake plots (tests/helpers/report_fixtures.py). Structure checks with html.parser only; no
browser (lead decision 1). Never touches the real ``$GIAP_DB_DIR``; no network.

Interface assumed beyond the spec text (M9-14 test author):
- ``giap.viz.interactive.interactive(pid)`` returns the decorated function unchanged and raises
  ``ValueError`` when ``pid`` already has a twin; ``figure_json(pid, data, title)`` is
  ``None`` for a plot without a twin, else a JSON object string with ``data`` and ``layout``,
  identical across calls;
- each twin's figure JSON sits in ``<script type="application/json"
  data-plotly="plotly-<pid>[-<suffix>]">`` (the value of ``data-plotly`` is the twin div's id);
  the ``Plotly.newPlot`` call is in an inline ``<script>`` naming that id;
- the suffix in the div id is the instance suffix (``V17_qc_b`` -> ``plotly-V17-b``);
- twin divs and the ``<details>`` holding the static SVG sit inside the plot's
  ``data-plot="Vnn"`` slot; the ``<summary>`` text contains "static figure";
- the footer is a ``<footer>`` element; it names plotly.js and MIT only when plotly.js is inlined
  (example 2: no ``plotly.js`` string at all under ``--no-interactive``).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner, Result
from helpers.interactive_fixtures import (
    BANNER,
    external_loads,
    head,
    page_of,
    plotly_scripts,
    register_twin,
    report_step,
    run_steps,
    swap_twins,
    twin_divs,
)
from helpers.report_fixtures import Page, register_fake, run_out, swap_registry

from giap.exitcodes import ExitCode

ROOT = Path(__file__).resolve().parents[3]
FAST = ["--set", "report.formats=[svg]", "--set", "report.dpi=50"]


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


@pytest.fixture
def twins(monkeypatch: pytest.MonkeyPatch, fake_registry: ModuleType) -> ModuleType:
    return swap_twins(monkeypatch)


def _plotlyjs() -> str:
    import plotly.offline

    text: str = plotly.offline.get_plotlyjs()
    return text


def _specs(page: Page) -> dict[str, dict[str, Any]]:
    """Twin div id -> its parsed figure JSON."""
    out: dict[str, dict[str, Any]] = {}
    for s in page.root.by_tag("script"):
        if s.attrs.get("type") == "application/json" and "data-plotly" in s.attrs:
            out[s.attrs["data-plotly"]] = json.loads(s.text())
    return out


def _newplot_scripts(page: Page) -> list[str]:
    return [s.text() for s in page.root.by_tag("script") if "Plotly.newPlot" in s.text()]


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def _cli(out: Path, *args: str) -> Result:
    from giap.cli import cli

    res = CliRunner().invoke(cli, ["report", str(out), *FAST, *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


# ---------------------------------------------------------------- the registry


def test_interactive_registers_and_figure_json_builds_from_the_frame(twins: ModuleType) -> None:
    calls: list[tuple[str, list[str], str]] = []
    fn = register_twin("V17", calls=calls)
    assert twins._REGISTRY["V17"] is fn
    data = pd.DataFrame({"x": [0, 1, 2], "y": [0.0, 1.0, 0.5]})
    text = twins.figure_json("V17", data, "QC title")
    assert isinstance(text, str)
    spec = json.loads(text)
    assert {"data", "layout"} <= set(spec)
    assert [float(v) for v in spec["data"][0]["x"]] == [0.0, 1.0, 2.0]
    assert [float(v) for v in spec["data"][0]["y"]] == [0.0, 1.0, 0.5]
    assert calls == [("V17", ["x", "y"], "QC title")]
    assert twins.figure_json("V17", data, "QC title") == text


def test_figure_json_is_none_without_a_twin(twins: ModuleType) -> None:
    register_twin("V17")
    assert twins.figure_json("V07", pd.DataFrame({"x": [1], "y": [2.0]}), "t") is None


def test_a_second_twin_for_one_plot_raises(twins: ModuleType) -> None:
    register_twin("V17")
    with pytest.raises(ValueError, match="V17"):
        register_twin("V17")


# ---------------------------------------------------------------- example 1: embedding


def test_example_1_plotly_inlined_once_in_head(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V07")
    register_fake("V17")
    register_twin("V17")
    html = run_steps(tmp_path).read_text(encoding="utf-8")
    assert html.count(BANNER) == 1
    page = page_of(tmp_path / "work" / "report" / "report.html")
    scripts = plotly_scripts(page)
    assert len(scripts) == 1
    assert scripts[0] in head(page).iter(), "plotly.js is not in <head>"
    assert "src" not in scripts[0].attrs
    assert _plotlyjs().strip() in scripts[0].text(), "plotly.js is not inlined verbatim"
    assert len(_newplot_scripts(page)) == 1
    assert html.count("Plotly.newPlot") == 1


def test_twin_div_spec_and_newplot_config(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V17")
    page = page_of(run_steps(tmp_path))
    slot = page.slot("V17")
    divs = twin_divs(slot)
    assert [d.attrs["id"] for d in divs] == ["plotly-V17"]
    assert set(_specs(page)) == {"plotly-V17"}
    (call,) = _newplot_scripts(page)
    assert re.search(r"""["']plotly-V17["']""", call), call[:300]
    for pattern in (
        r"""["']?displaylogo["']?\s*:\s*false""",
        r"""["']?responsive["']?\s*:\s*true""",
        r"""["']?toImageButtonOptions["']?\s*:\s*\{\s*["']?format["']?\s*:\s*["']svg["']""",
    ):
        assert re.search(pattern, call), (pattern, call[:300])


def test_static_svg_stays_in_details_beside_the_twin(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V17")
    slot = page_of(run_steps(tmp_path)).slot("V17")
    details = slot.by_tag("details")
    assert len(details) == 1
    assert details[0].by_tag("svg"), "no static SVG inside <details>"
    summary = details[0].by_tag("summary")
    assert summary and "static figure" in summary[0].text()
    assert not twin_divs(details[0]), "the twin must sit beside the static figure, not in it"


def test_acceptance_v17_spec_equals_data_tsv_values(twins: ModuleType, tmp_path: Path) -> None:
    calls: list[tuple[str, list[str], str]] = []
    register_fake("V17")
    register_twin("V17", calls=calls)
    page = page_of(run_steps(tmp_path))
    data = _tsv(tmp_path / "work" / "plots" / "data" / "V17_qc.tsv")
    trace = _specs(page)["plotly-V17"]["data"][0]
    assert [float(v) for v in trace["x"]] == [float(v) for v in data["x"]]
    assert [float(v) for v in trace["y"]] == [float(v) for v in data["y"]]
    assert [c[:2] for c in calls] == [("V17", ["x", "y"])]
    assert isinstance(calls[0][2], str) and calls[0][2].strip(), "twin got no title"


def test_twin_is_built_from_the_data_tsv_not_the_tables(twins: ModuleType, tmp_path: Path) -> None:
    """Behaviour 1: editing the data TSV between ``plots`` and ``report`` changes the twin."""
    register_fake("V17")
    register_twin("V17")
    out = run_out(tmp_path)
    run_steps(tmp_path, out=out)
    tsv = tmp_path / "work" / "plots" / "data" / "V17_qc.tsv"
    tsv.write_text("x\ty\n5\t7.25\n6\t8.5\n", encoding="utf-8")
    page = page_of(report_step(tmp_path, plots=tmp_path / "work" / "plots", out=out))
    trace = _specs(page)["plotly-V17"]["data"][0]
    assert [float(v) for v in trace["x"]] == [5.0, 6.0]
    assert [float(v) for v in trace["y"]] == [7.25, 8.5]


def test_one_twin_per_instance_with_suffix_ids(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17", suffixes=(None, "b"))
    register_twin("V17")
    page = page_of(run_steps(tmp_path))
    ids = sorted(d.attrs["id"] for d in twin_divs(page.slot("V17")))
    assert ids == ["plotly-V17", "plotly-V17-b"]
    assert set(_specs(page)) == set(ids)
    assert page.html.count("Plotly.newPlot") == 2  # one call per twin (behaviour 2)
    calls = _newplot_scripts(page)
    for i in ids:
        assert [c for c in calls if re.search(rf"""["']{i}["']""", c)], i
    assert len(plotly_scripts(page)) == 1
    all_ids = page.ids()
    assert not sorted({i for i in all_ids if all_ids.count(i) > 1})


def test_plot_without_twin_stays_static_when_another_has_one(
    twins: ModuleType, tmp_path: Path
) -> None:
    register_fake("V07")
    register_fake("V17")
    register_twin("V17")
    page = page_of(run_steps(tmp_path))
    v07 = page.slot("V07")
    assert not twin_divs(v07)
    assert v07.by_tag("svg")


# ---------------------------------------------------------------- degenerate rows


def test_no_twin_for_any_rendered_plot_means_no_plotly(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V07")  # V07 is not rendered in this run
    html = run_steps(tmp_path, "report.interactive=true").read_text(encoding="utf-8")
    assert BANNER not in html
    assert "Plotly.newPlot" not in html
    assert not twin_divs(page_of(tmp_path / "work" / "report" / "report.html").root)


def test_twin_of_a_skipped_plot_means_no_plotly(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V11", needs=["pairs"])  # skips: input_absent:pairs
    register_twin("V11")
    html = run_steps(tmp_path).read_text(encoding="utf-8")
    assert BANNER not in html and "Plotly.newPlot" not in html


def test_a_twin_that_raises_fails_the_report_step(twins: ModuleType, tmp_path: Path) -> None:
    from giap.steps import registry as reg

    register_fake("V17")
    register_twin("V17", raises=True)
    with pytest.raises(reg.StepError, match="twin V17 failed"):
        run_steps(tmp_path)


def test_a_twin_that_raises_is_exit_4_in_giap_report(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V17", raises=True)
    out = run_out(tmp_path)
    res = _cli(out)
    assert res.exit_code == ExitCode.STEP_FAILED, res.output


# ---------------------------------------------------------------- example 2: static only


def _assert_static_only(html: str, page: Page) -> None:
    assert "plotly.js" not in html.lower()
    assert "Plotly" not in html
    assert not twin_divs(page.root)
    assert not _specs(page)
    slot = page.slot("V17")
    assert slot.by_tag("svg"), "V17 is not an inline SVG"


def test_example_2_interactive_false_is_static_only(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V17")
    path = run_steps(tmp_path, "report.interactive=false")
    _assert_static_only(path.read_text(encoding="utf-8"), page_of(path))


def test_example_2_cli_no_interactive_is_static_only(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V17")
    out = run_out(tmp_path)
    res = _cli(out, "--no-interactive")
    assert res.exit_code == ExitCode.SUCCESS, res.output
    path = out / "report" / "report.html"
    _assert_static_only(path.read_text(encoding="utf-8"), page_of(path))


def test_cli_default_is_interactive(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V17")
    out = run_out(tmp_path)
    res = _cli(out)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    html = (out / "report" / "report.html").read_text(encoding="utf-8")
    assert html.count(BANNER) == 1
    assert twin_divs(page_of(out / "report" / "report.html").slot("V17"))


def test_cli_set_report_interactive_false_is_static_only(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V17")
    out = run_out(tmp_path)
    res = _cli(out, "--set", "report.interactive=false")
    assert res.exit_code == ExitCode.SUCCESS, res.output
    path = out / "report" / "report.html"
    _assert_static_only(path.read_text(encoding="utf-8"), page_of(path))


# ---------------------------------------------------------------- example 3: no external loads


def test_external_loads_checker_flags_loads() -> None:
    """The local checker itself (stand-in for M9-13's ``giap.report.checks.external_loads``)."""
    bad = (
        '<html><head><script src="https://cdn.plot.ly/x.js"></script>'
        '<link rel="stylesheet" href="a.css"><style>@import "https://f/x.css";'
        "body{background:url(http://x/y.png)}</style></head><body>"
        '<img src="//cdn/x.png"><a href="https://example.org">ok</a>'
        "<script>var u = 'https://not.a.load/';</script></body></html>"
    )
    found = external_loads(bad)
    assert "https://cdn.plot.ly/x.js" in found and "a.css" in found
    assert "//cdn/x.png" in found and "http://x/y.png" in found
    assert "https://f/x.css" in found
    assert "https://example.org" not in found and "https://not.a.load/" not in found


def test_example_3_no_external_loads(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V07")
    register_fake("V17", suffixes=(None, "b"))
    register_twin("V17")
    html = run_steps(tmp_path, out=run_out(tmp_path, islands=3)).read_text(encoding="utf-8")
    assert BANNER in html
    assert external_loads(html) == []
    page = page_of(tmp_path / "work" / "report" / "report.html")
    assert not [s for s in page.root.by_tag("script") if s.attrs.get("src")]


# ---------------------------------------------------------------- licence (behaviour 6)


def test_footer_names_plotly_js_and_mit(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V17")
    page = page_of(run_steps(tmp_path))
    footers = page.root.by_tag("footer")
    assert len(footers) == 1
    text = " ".join(footers[0].text().split())
    assert "plotly.js" in text and "MIT" in text, text
    assert footers[0] not in head(page).iter()


def test_notice_md_has_a_plotly_js_mit_row() -> None:
    notice = (ROOT / "NOTICE.md").read_text(encoding="utf-8")
    rows = [ln for ln in notice.splitlines() if re.match(r"^\|\s*plotly\.js\b", ln)]
    assert len(rows) == 1, rows
    cells = [c.strip() for c in rows[0].strip().strip("|").split("|")]
    assert cells[1] == "MIT", cells
    assert "4.1.1" in rows[0], rows[0]


# ---------------------------------------------------------------- determinism


def test_building_twice_gives_identical_html(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V07")
    register_fake("V17", suffixes=(None, "b"))
    register_twin("V17")
    out = run_out(tmp_path, islands=3)
    first = run_steps(tmp_path, out=out).read_bytes()
    second = run_steps(tmp_path, out=out).read_bytes()
    assert BANNER.encode() in first
    assert first == second


def test_cli_twice_gives_identical_html(twins: ModuleType, tmp_path: Path) -> None:
    register_fake("V17")
    register_twin("V17")
    out = run_out(tmp_path)
    assert _cli(out).exit_code == 0
    first = (out / "report" / "report.html").read_bytes()
    assert _cli(out).exit_code == 0
    assert (out / "report" / "report.html").read_bytes() == first
