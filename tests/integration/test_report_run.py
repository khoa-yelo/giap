"""M9-02 examples 1, 2, 3, 6 and 7 on a real ``giap run`` of sim100k (giapsim seed 20261014;
docs/specs/M9-02-plots-report-steps-and-cli.md; GIAP_SPEC.md §6.4 executor equality, §7.9.4).

Inputs as tests/integration/test_m6_08_operons_run.py: the giapsim SIM100K genome with M3-19
proteomics and every precomputed layer, ``--modules core``, one seed, the test DB of
tests/helpers/mapping_steps.py. Whatever plots are registered when this runs (none in M9-02
itself, every pass-1 plot by M9-13) must either render or be noted as not run in the report.

* example 1: ``--executor local`` writes ``plots/`` (with ``data/`` for every rendered plot) and
  ``report/report.html``; exit 0;
* examples 6 and 7 on that report: 13 sections in order, no external URL, every relative link
  resolves;
* example 3: ``giap report OUT`` twice gives identical ``plots/data/*.tsv`` and ``report.html``;
* example 2: ``--executor nextflow --profile farmshare_local`` gives byte-identical
  ``plots/data/*.tsv`` (beside the existing executor-equality tests).
"""

from __future__ import annotations

import dataclasses
import re
import shutil
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.report_fixtures import SECTION_IDS, parse

pytestmark = [pytest.mark.slow]

SEED = 7
EXTERNAL = re.compile(r"^\s*(https?:|//|file:)", re.I)


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim
    from giapsim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        giapsim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    assert giapsim.SIM100K_SEED == 20261014
    return giapsim.simulate(cfg, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")


@pytest.fixture(scope="module")
def db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.mapping_steps import write_test_db

    return Path(write_test_db(tmp_path_factory.mktemp("test_db_m9_02")))


def _invoke(args: list[str], db: Path) -> Any:
    from giap.cli import cli

    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIAP_DB_DIR", str(db))
        mp.setenv("NXF_ANSI_LOG", "false")
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _run(sim: Any, out: Path, db: Path, *extra: str) -> None:
    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides)]
    args += ["--abundance", str(sim.abundance), "--annotation", str(sim.bakta_json)]
    args += ["--config", str(sim.config_yaml), "--outdir", str(out), "--modules", "core"]
    args += ["--seed", str(SEED), *extra]
    res = _invoke(args, db)
    assert res.exit_code == 0, res.output


@pytest.fixture(scope="module")
def local(sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("LOCAL") / "OUT"
    _run(sim, out, db, "--executor", "local")
    return out


def _data(out: Path) -> dict[str, bytes]:
    return {p.name: p.read_bytes() for p in sorted((out / "plots" / "data").glob("*.tsv"))}


def test_example_1_plots_and_report_published(local: Path) -> None:
    from giap.viz.base import registered

    assert (local / "plots").is_dir()
    html = local / "report" / "report.html"
    assert html.is_file()
    page = parse(html.read_text(encoding="utf-8"))
    for pid, spec in registered().items():
        svgs = list((local / "plots").glob(f"{pid}_*.svg"))
        if svgs:
            assert (local / "plots" / "data" / f"{svgs[0].stem}.tsv").is_file(), pid
        elif page.root.by_attr("data-plot", pid):
            assert "not run" in page.slot(pid).text(), (pid, spec.name)


def test_example_6_thirteen_sections(local: Path) -> None:
    page = parse((local / "report" / "report.html").read_text(encoding="utf-8"))
    assert [s.attrs.get("id") for s in page.sections()] == list(SECTION_IDS)


def test_example_7_offline_and_links_resolve(local: Path) -> None:
    page = parse((local / "report" / "report.html").read_text(encoding="utf-8"))
    for url in page.links():
        assert not EXTERNAL.match(url), url
        if not url.startswith(("#", "data:")):  # inlined images (V13 colorbar; M9-13: allowed)
            assert (local / "report" / url).resolve().is_file(), url
    assert "../tables/islands.tsv" in page.links()
    assert str(local.parent) not in page.html


def test_example_3_giap_report_twice_identical(
    local: Path, db: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    out = tmp_path_factory.mktemp("REPORT") / "OUT"
    shutil.copytree(local, out, ignore=shutil.ignore_patterns("work"))
    runs = []
    for _ in range(2):
        res = _invoke(["report", str(out)], db)
        assert res.exit_code == 0, res.output
        runs.append((_data(out), (out / "report" / "report.html").read_bytes()))
    assert runs[0] == runs[1]
    # giap report rebuilds the same plot data as the run itself (no recomputed statistic)
    assert runs[0][0] == _data(local)


@pytest.mark.needs_nextflow
def test_example_2_nextflow_gives_identical_plot_data(
    sim: Any, db: Path, local: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    nf = tmp_path_factory.mktemp("NF") / "OUT"
    _run(sim, nf, db, "--executor", "nextflow", "--profile", "farmshare_local")
    assert (nf / "report" / "report.html").is_file(), "Nextflow run did not publish the report"
    assert (nf / "plots").is_dir(), "Nextflow run did not publish plots/"
    assert _data(nf) == _data(local)
