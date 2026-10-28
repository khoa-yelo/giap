"""``gippa report OUT`` (M9-02 behaviour 7, examples 3, 4, 5, 7, 8 and the degenerate cases;
GIPPA_SPEC.md §5.1 exit codes 0, 1, 4; §7.9.4 "Running it twice MUST produce identical
``plots/data/*.tsv``").

Fake plots (tests/helpers/report_fixtures.py) stand in for the real ones, so these tests do not
depend on which plot issues have merged. ``--set report.formats=[svg] --set report.dpi=50`` keeps
them fast unless a test needs the other formats. Never touches the real ``$GIPPA_DB_DIR``.

Interface assumed beyond the spec text: the report path is printed on stdout; error messages go to
the combined CliRunner output; the HTML hooks of tests/helpers/report_fixtures.py.
"""

from __future__ import annotations

import re
from pathlib import Path
from types import ModuleType

import pytest
from click.testing import CliRunner, Result
from helpers.report_fixtures import (
    SECTION_IDS,
    parse,
    register_fake,
    run_out,
    swap_registry,
)
from helpers.viz_out_tree import write_resolved_config, write_schema_table

from gippa.exitcodes import ExitCode

FAST = ["--set", "report.formats=[svg]", "--set", "report.dpi=50"]
EXTERNAL = re.compile(r"^\s*(https?:|//|file:)", re.I)


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _report(out: Path, *args: str, fast: bool = True) -> Result:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["report", str(out), *(FAST if fast else []), *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _pairs(out: Path) -> Path:
    base = {"class": "adjacent", "gene_a": "g1", "gene_b": "g2", "strand": "+", "test": "none"}
    return write_schema_table(out / "tables" / "pairs.tsv", "pairs", [{**base, "pair_id": "P1"}])


def _data(out: Path) -> dict[str, bytes]:
    return {p.name: p.read_bytes() for p in sorted((out / "plots" / "data").glob("*.tsv"))}


def _two_plots(tmp_path: Path) -> Path:
    register_fake("V07")
    register_fake("V17", suffixes=(None, "b"))
    return run_out(tmp_path, islands=3)


# ---------------------------------------------------------------- exit 0, outputs


def test_report_is_implemented_and_writes_plots_and_report(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = _two_plots(tmp_path)
    res = _report(out)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert "not implemented" not in res.output
    assert (out / "report" / "report.html").is_file()
    assert "report.html" in res.stdout
    for name in ("V07_obs_vs_exp.svg", "V17_qc.svg", "V17_qc_b.svg"):
        assert (out / "plots" / name).is_file(), name
    assert set(_data(out)) == {"V07_obs_vs_exp.tsv", "V17_qc.tsv", "V17_qc_b.tsv"}
    page = parse((out / "report" / "report.html").read_text(encoding="utf-8"))
    assert [s.attrs.get("id") for s in page.sections()] == list(SECTION_IDS)


def test_no_interactive_is_accepted(fake_registry: ModuleType, tmp_path: Path) -> None:
    out = _two_plots(tmp_path)
    res = _report(out, "--no-interactive")
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert (out / "report" / "report.html").is_file()


def test_example_3_twice_gives_identical_data_and_report(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = _two_plots(tmp_path)
    assert _report(out).exit_code == 0
    first = (_data(out), (out / "report" / "report.html").read_bytes())
    assert _report(out).exit_code == 0
    second = (_data(out), (out / "report" / "report.html").read_bytes())
    assert first[0] and first[0] == second[0]
    assert first[1] == second[1]


def test_example_7_links_resolve_and_nothing_is_external(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = _two_plots(tmp_path)
    res = _report(out, "--set", "report.formats=[svg, png, pdf]", "--set", "report.dpi=50",
                  fast=False)  # fmt: skip
    assert res.exit_code == 0, res.output
    page = parse((out / "report" / "report.html").read_text(encoding="utf-8"))
    links = page.links()
    assert "../tables/islands.tsv" in links
    assert "../plots/V07_obs_vs_exp.png" in links and "../plots/V07_obs_vs_exp.pdf" in links
    assert "../plots/data/V17_qc_b.tsv" in links
    for url in links:
        assert not EXTERNAL.match(url), url
        if url.startswith("#"):
            continue
        assert (out / "report" / url).resolve().is_file(), url
    assert str(tmp_path) not in page.html


# ---------------------------------------------------------------- exit 1


def test_example_4_unknown_plot_id_exits_1(fake_registry: ModuleType, tmp_path: Path) -> None:
    register_fake("V01")
    out = run_out(tmp_path)
    res = _report(out, "--plots", "V1,V99")
    assert res.exit_code == ExitCode.INVALID_INPUT, res.output
    assert "unknown plot ID V99" in res.output
    assert not (out / "report" / "report.html").exists()


def test_plots_accepts_v1_and_v01(fake_registry: ModuleType, tmp_path: Path) -> None:
    register_fake("V01")
    register_fake("V17")
    out = run_out(tmp_path)
    assert _report(out, "--plots", "V1").exit_code == 0
    assert (out / "plots" / "V01_circos.svg").is_file()
    assert not (out / "plots" / "V17_qc.svg").exists()
    (out / "plots" / "V01_circos.svg").unlink()
    assert _report(out, "--plots", "V01").exit_code == 0
    assert (out / "plots" / "V01_circos.svg").is_file()


def test_example_5_empty_directory_exits_1(tmp_path: Path) -> None:
    empty = tmp_path / "empty_dir"
    empty.mkdir()
    res = _report(empty)
    assert res.exit_code == ExitCode.INVALID_INPUT, res.output
    assert "no tables/ or preflight/" in res.output
    assert not (empty / "report").exists() and not (empty / "plots").exists()


def test_missing_out_directory_exits_1(tmp_path: Path) -> None:
    from gippa.cli import main  # main maps click usage errors to 1 (test_cli_exit_codes.py)

    with pytest.raises(SystemExit) as info:
        main(["report", str(tmp_path / "nope")])
    assert info.value.code == ExitCode.INVALID_INPUT


# ---------------------------------------------------------------- exit 4


def test_a_raising_plot_exits_4_naming_the_plot(fake_registry: ModuleType, tmp_path: Path) -> None:
    register_fake("V07")
    register_fake("V17", raises=True)
    res = _report(run_out(tmp_path))
    assert res.exit_code == ExitCode.STEP_FAILED, res.output
    assert "V17" in res.output


# ---------------------------------------------------------------- subsets and stale files


def test_plots_subset_rebuilds_only_that_plot(fake_registry: ModuleType, tmp_path: Path) -> None:
    out = _two_plots(tmp_path)
    assert _report(out).exit_code == 0
    marker = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"><g id="marker17"/></svg>'
    )
    (out / "plots" / "V17_qc.svg").write_text(marker, encoding="utf-8")
    (out / "plots" / "data" / "V07_obs_vs_exp.tsv").unlink()
    res = _report(out, "--plots", "V07")
    assert res.exit_code == 0, res.output
    assert (out / "plots" / "data" / "V07_obs_vs_exp.tsv").is_file()  # V07 rebuilt
    assert (out / "plots" / "V17_qc.svg").read_text(encoding="utf-8") == marker  # V17 left alone
    html = (out / "report" / "report.html").read_text(encoding="utf-8")
    assert "marker17" in html  # the report is rebuilt from every present plot file


def test_example_8_deleted_svg_not_rebuilt_and_noted(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = _two_plots(tmp_path)
    assert _report(out).exit_code == 0
    (out / "plots" / "V07_obs_vs_exp.svg").unlink()
    res = _report(out, "--plots", "V17")
    assert res.exit_code == 0, res.output
    assert not (out / "plots" / "V07_obs_vs_exp.svg").exists()
    page = parse((out / "report" / "report.html").read_text(encoding="utf-8"))
    slot = page.slot("V07")
    assert slot in page.section("detection").iter()
    assert "plot file absent" in " ".join(slot.text().split())


def test_a_plot_that_now_skips_leaves_no_stale_files(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    register_fake("V11", needs=["pairs"])
    out = run_out(tmp_path)
    _pairs(out)
    assert _report(out).exit_code == 0
    assert (out / "plots" / "V11_operon_pairs.svg").is_file()
    (out / "tables" / "pairs.tsv").unlink()
    assert _report(out).exit_code == 0
    leftovers = [p for p in (out / "plots").rglob("V11_*")]
    assert not leftovers, leftovers
    page = parse((out / "report" / "report.html").read_text(encoding="utf-8"))
    text = " ".join(page.slot("V11").text().split())
    assert "not run" in text and "pairs" in text, text


def test_older_plot_files_of_other_plots_are_kept(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    register_fake("V07")
    out = run_out(tmp_path)
    (out / "plots" / "data").mkdir(parents=True)
    (out / "plots" / "V42_legacy.svg").write_text("<svg/>")
    (out / "plots" / "notes.txt").write_text("keep")
    (out / "plots" / "V07_obs_vs_exp_old.svg").write_text("<svg/>")
    assert _report(out).exit_code == 0
    assert (out / "plots" / "V42_legacy.svg").read_text() == "<svg/>"
    assert (out / "plots" / "notes.txt").read_text() == "keep"
    # V07 is rebuilt: its own older instance files go (remove_plot_files), the new one is written
    assert not (out / "plots" / "V07_obs_vs_exp_old.svg").exists()
    assert (out / "plots" / "V07_obs_vs_exp.svg").is_file()


# ---------------------------------------------------------------- modules off


def test_function_off_notes_v14_and_exits_0(fake_registry: ModuleType, tmp_path: Path) -> None:
    calls: list[str] = []
    register_fake("V14", module="function", calls=calls)
    out = run_out(tmp_path)
    write_resolved_config(out / "resolved_config.yaml", ["core", "mobile", "bgc"])
    res = _report(out)
    assert res.exit_code == 0, res.output
    assert calls == []
    assert not list((out / "plots").glob("V14_*"))
    page = parse((out / "report" / "report.html").read_text(encoding="utf-8"))
    text = " ".join(page.slot("V14").text().split())
    assert "not run" in text and "function" in text and "off" in text, text
