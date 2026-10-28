"""M9-13, the pass-1 closer: every pass-1 plot in ``report.html`` from a real sim100k ``gippa run``
(docs/specs/M9-13-pass1-integration.md examples 1-5 and the zero-detection degenerate row;
GIPPA_SPEC.md §6.4 executor equality, §7.9.4 ``gippa report`` idempotence, §10.3 link and offline
checks, R7/TE7 extended to ``plots/data/``).

Inputs, as tests/integration/test_report_run.py (M9-02) and tests/integration/
test_m2_21_local_vs_nextflow.py (the existing executor-equality test): the gippasim ``SIM100K``
genome at the default seed ``SIM100K_SEED`` with M3-19 proteomics and every precomputed layer
(``config.yaml``), peptides and abundance, ``--modules core``, ``--seed 7``, the test DB of
tests/helpers/mapping_steps.py in ``GIPPA_DB_DIR`` (never the real one); default report formats.
No network.

Interface choices (M9-13 test author, beyond the spec text):

- the pass-1 plot list of example 1 is pinned as ``(plot ID, §10.2 name)``; a plot missing from
  ``gippa.viz.base.registered()`` fails, it is never skipped (lead decision 2). V03 is satisfied by
  its index ``plots/data/V03_locus.tsv`` (M9-10 decision 4); every other plot by at least one
  ``plots/<ID>_<name>[_<suffix>].svg`` whose ``plots/data/<stem>.tsv`` exists;
- the events are read from ``OUT/logs/events.jsonl`` (``{"event": ..., "data": {...}}``, §6.7):
  each registered plot has exactly one ``plots.rendered`` (``data.files`` relative to ``plots/``,
  all on disk) or one ``plots.skipped`` with a non-empty declared reason
  (``module_off:``/``layer_missing:``/``input_absent:``) -- never both;
- examples 2 and 3 use ``gippa.report.checks`` (tests/unit/report/test_checks.py pins it); the
  relative links of ``report/report.html`` resolve against ``OUT/report/``;
- example 4 runs ``--executor nextflow --profile farmshare_local`` (real Nextflow, prebuilt env)
  like every existing executor-equality test: the spec's "-profile test" with stubs cannot give
  byte-identical tables (see tests/integration/test_novel_local_vs_nextflow.py). Compared: the
  recursive file list and bytes of ``plots/data/`` (``diff -r``);
- example 5 copies ``OUT`` (without ``work/``) and runs ``gippa report OUT`` twice;
- zero detections: the same run with every peptide sequence replaced by a distinct
  ``[WCMH]{11}K`` string (absent from a 100-kb genome's six-frame translation) and no abundance;
  precondition: ``tables/genes.tsv`` has no detected gene. The run exits 0, writes
  ``report/report.html``, and V07, V10 and V18 still draw (an empty state is a rendered plot, not a
  skip), so each has an SVG and a data TSV and its report slot is not "not run".
"""

from __future__ import annotations

import csv
import dataclasses
import json
import re
import shutil
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

pytestmark = [pytest.mark.slow, pytest.mark.integration]

SEED = 7
#: example 1: the pass-1 plots (ID, §10.2 name), every one rendered on sim100k
PASS1 = (
    ("V01", "circos"),
    ("V02", "abundance_landscape"),
    ("V03", "locus"),
    ("V05", "mapping_sankey"),
    ("V06", "mapping_fdr"),
    ("V07", "obs_vs_exp"),
    ("V08", "cog_abundance"),
    ("V10", "correlogram"),
    ("V11", "operon_pairs"),
    ("V13", "novel_orfs"),
    ("V17", "qc"),
    ("V18", "stat_diagnostics"),
)
REASON = re.compile(r"^(module_off|layer_missing|input_absent):\S+$")
ZERO_EMPTY = ("V07", "V10", "V18")


# ---------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import gippasim
    from gippasim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        gippasim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    assert gippasim.SIM100K_SEED == 20261014
    return gippasim.simulate(cfg, gippasim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")


@pytest.fixture(scope="module")
def db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.mapping_steps import write_test_db

    return Path(write_test_db(tmp_path_factory.mktemp("test_db_m9_13")))


def _invoke(args: list[str], db: Path) -> Any:
    from gippa.cli import cli

    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIPPA_DB_DIR", str(db))
        mp.setenv("NXF_ANSI_LOG", "false")
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _run(sim: Any, out: Path, db: Path, *extra: str, peptides: Path | None = None,
         abundance: bool = True) -> None:  # fmt: skip
    args = ["run", "--genome", str(sim.genome_fasta)]
    args += ["--peptides", str(peptides or sim.peptides)]
    if abundance:
        args += ["--abundance", str(sim.abundance)]
    args += ["--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", "core", "--seed", str(SEED), *extra]
    res = _invoke(args, db)
    assert res.exit_code == 0, res.output


@pytest.fixture(scope="module")
def local(sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("LOCAL") / "OUT"
    _run(sim, out, db, "--executor", "local")
    return out


def _registered() -> dict[str, Any]:
    from gippa.steps.plots import load_plot_modules
    from gippa.viz.base import registered

    load_plot_modules()
    return dict(registered())


def _events(out: Path) -> list[dict[str, Any]]:
    path = out / "logs" / "events.jsonl"
    assert path.is_file(), f"no {path}"
    lines = path.read_text(encoding="utf-8").splitlines()
    evs = [json.loads(ln) for ln in lines if ln.strip()]
    return [e for e in evs if e.get("event") in ("plots.rendered", "plots.skipped")]


def _data_tree(out: Path) -> dict[str, bytes]:
    root = out / "plots" / "data"
    return {p.relative_to(root).as_posix(): p.read_bytes()
            for p in sorted(root.rglob("*")) if p.is_file()}  # fmt: skip


def _data_tsvs(out: Path) -> dict[str, bytes]:
    return {p.name: p.read_bytes() for p in sorted((out / "plots" / "data").glob("*.tsv"))}


def _svgs(out: Path, pid: str, name: str) -> list[Path]:
    plots = out / "plots"
    return sorted([*plots.glob(f"{pid}_{name}.svg"), *plots.glob(f"{pid}_{name}_*.svg")])


def _html(out: Path) -> str:
    page = out / "report" / "report.html"
    assert page.is_file(), f"no {page}"
    return page.read_text(encoding="utf-8")


# ---------------------------------------------------------------- example 1


@pytest.mark.parametrize("pid, name", PASS1)
def test_example_1_pass1_plot_is_registered(pid: str, name: str) -> None:
    reg = _registered()
    assert pid in reg, f"{pid} ({name}) is not registered: the pass-1 list needs it (decision 2)"
    assert reg[pid].name == name


@pytest.mark.parametrize("pid, name", PASS1)
def test_example_1_pass1_plot_files_written(local: Path, pid: str, name: str) -> None:
    assert pid in _registered(), f"{pid} ({name}) is not registered (decision 2)"
    data = local / "plots" / "data"
    if pid == "V03":
        assert (data / "V03_locus.tsv").is_file(), "the V03 index data/V03_locus.tsv is missing"
        return
    svgs = _svgs(local, pid, name)
    assert svgs, f"no plots/{pid}_{name}[_*].svg"
    for svg in svgs:
        assert (data / f"{svg.stem}.tsv").is_file(), f"no data TSV for {svg.name}"


def test_example_1_report_written(local: Path) -> None:
    assert (local / "report" / "report.html").is_file()


def test_every_registered_plot_rendered_or_skipped_with_reason(local: Path) -> None:
    evs = _events(local)
    for pid in _registered():
        mine = [e for e in evs if e.get("data", {}).get("plot") == pid]
        assert len(mine) == 1, (pid, mine)
        ev = mine[0]
        if ev["event"] == "plots.skipped":
            reason = str(ev["data"].get("reason", ""))
            # V16 is written only by `gippa compare` (M9-23 decision 2): its run-time skip says so
            assert REASON.match(reason), (pid, ev)  # V16: input_absent:compare (M9-23 amended)
        else:
            files = ev["data"]["files"]
            assert files, (pid, ev)
            for rel in files:
                assert (local / "plots" / rel).is_file(), (pid, rel)


def test_no_plot_event_for_an_unregistered_id(local: Path) -> None:
    reg = _registered()
    stray = sorted({e["data"]["plot"] for e in _events(local)} - set(reg))
    assert not stray, stray


@pytest.mark.parametrize("pid, name", PASS1)
def test_example_1_pass1_plots_render_not_skip_on_sim100k(local: Path, pid: str, name: str) -> None:
    rendered = {e["data"]["plot"] for e in _events(local) if e["event"] == "plots.rendered"}
    assert pid in rendered, f"{pid} ({name}) did not render on sim100k"


# ---------------------------------------------------------------- examples 2 and 3


def test_example_2_report_loads_nothing_external(local: Path) -> None:
    from gippa.report.checks import external_loads

    assert external_loads(_html(local)) == []


def test_example_3_every_relative_link_exists(local: Path) -> None:
    from gippa.report.checks import relative_links

    links = relative_links(_html(local))
    assert links, "report.html links to no file under OUT/"
    report_dir = local / "report"
    broken = [u for u in links if not (report_dir / u).resolve().exists()]
    assert not broken, broken
    for u in links:
        assert (report_dir / u).resolve().is_relative_to(local.resolve()), u


def test_example_3_report_links_every_rendered_data_tsv(local: Path) -> None:
    from gippa.report.checks import relative_links

    links = {(local / "report" / u).resolve() for u in relative_links(_html(local))}
    for e in _events(local):
        if e["event"] != "plots.rendered":
            continue
        for rel in e["data"]["files"]:
            if rel.startswith("data/") and rel.endswith(".tsv") and rel != "data/V03_locus.tsv":
                assert (local / "plots" / rel).resolve() in links, rel


# ---------------------------------------------------------------- example 5


def test_example_5_gippa_report_twice_identical(
    local: Path, db: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    out = tmp_path_factory.mktemp("REPORT13") / "OUT"
    shutil.copytree(local, out, ignore=shutil.ignore_patterns("work"))
    runs = []
    for _ in range(2):
        res = _invoke(["report", str(out)], db)
        assert res.exit_code == 0, res.output
        runs.append((_data_tsvs(out), (out / "report" / "report.html").read_bytes()))
    assert runs[0][0], "gippa report wrote no plots/data/*.tsv"
    assert runs[0][0] == runs[1][0], "plots/data/*.tsv differ between two gippa report runs"
    assert runs[0][1] == runs[1][1], "report.html differs between two gippa report runs"
    assert runs[0][0] == _data_tsvs(local), "gippa report rebuilt different plot data than the run"


# ---------------------------------------------------------------- example 4


@pytest.mark.needs_nextflow
def test_example_4_nextflow_plots_data_identical(
    sim: Any, db: Path, local: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    nf = tmp_path_factory.mktemp("NF13") / "OUT"
    _run(sim, nf, db, "--executor", "nextflow", "--profile", "farmshare_local")
    assert (nf / "report" / "report.html").is_file(), "Nextflow run did not publish the report"
    a, b = _data_tree(local), _data_tree(nf)
    assert a, "local run wrote no plots/data"
    assert sorted(b) == sorted(a), "plots/data file sets differ between executors"
    differ = sorted(k for k in a if a[k] != b[k])
    assert not differ, f"plots/data differ between executors: {differ}"


# ---------------------------------------------------------------- zero detections


def _unmappable(i: int) -> str:
    letters = []
    for _ in range(11):
        letters.append("WCMH"[i % 4])
        i //= 4
    return "".join(letters) + "K"


@pytest.fixture(scope="module")
def zero(sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    base = tmp_path_factory.mktemp("ZERO")
    with open(sim.peptides, encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    col = rows[0].index("peptide")
    codes: dict[str, str] = {}
    for row in rows[1:]:
        row[col] = codes.setdefault(row[col], _unmappable(len(codes)))
    peptides = base / "peptides_unmappable.tsv"
    with open(peptides, "w", encoding="utf-8", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    out = base / "OUT"
    _run(sim, out, db, "--executor", "local", peptides=peptides, abundance=False)
    return out


def test_zero_detections_precondition(zero: Path) -> None:
    with open(zero / "tables" / "genes.tsv", encoding="utf-8", newline="") as fh:
        genes = list(csv.DictReader(fh, delimiter="\t"))
    assert genes, "genes.tsv is empty"
    detected = [g for g in genes if g["detected"].strip().lower() in ("true", "1")]
    assert not detected, f"{len(detected)} genes detected in the zero-detection run"


def test_zero_detections_report_written(zero: Path) -> None:
    from gippa.report.checks import external_loads, relative_links

    html = _html(zero)
    assert external_loads(html) == []
    for u in relative_links(html):
        assert (zero / "report" / u).resolve().exists(), u


@pytest.mark.parametrize("pid", ZERO_EMPTY)
def test_zero_detections_empty_states_drawn(zero: Path, pid: str) -> None:
    from helpers.report_fixtures import parse

    names = dict(PASS1)
    assert pid in _registered(), f"{pid} is not registered (decision 2)"
    svgs = _svgs(zero, pid, names[pid])
    assert svgs, f"{pid} drew nothing on the zero-detection run (expected an empty state)"
    for svg in svgs:
        assert (zero / "plots" / "data" / f"{svg.stem}.tsv").is_file(), svg.name
    page = parse(_html(zero))
    assert "not run" not in page.slot(pid).text(), f"{pid} shown as not run"
