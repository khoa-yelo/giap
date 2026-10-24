"""M9-16b examples 1, 2, 3 and 9 on a real ``giap run`` of sim100k (giapsim seed 20261014;
docs/specs/M9-16b-report-after-manifest.md; acceptance "tests/integration/test_report_run.py
(extended)", kept in its own module because the role guard forbids editing existing tests).

Inputs as tests/integration/test_report_run.py: the giapsim SIM100K genome with M3-19 proteomics
and every precomputed layer, ``--modules core``, one seed, the test DB of
tests/helpers/mapping_steps.py (never the real ``$GIAP_DB_DIR``). Slow: run on Slurm
(scripts/slurm/test_branch.sbatch ... -m slow).

* example 1: ``--executor local``: ``OUT/report/report.html`` carries the manifest's provenance
  (version, commit, a command line with ``--executor local``, seed, file digest, provenance
  digest) and not the "not written" sentence; ``OUT/work/report/report.html`` still has it;
* example 2: the same on ``OUT/report/report.html`` with the Nextflow executor;
* example 3: ``giap report OUT`` right after gives the same ``report.html`` bytes;
* example 9: no key or value of ``run_manifest.json`` names ``report/report.html`` or a
  ``plots/`` file.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import re
import shutil
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.report_fixtures import Node, parse

import giap

pytestmark = [pytest.mark.slow]

SEED = 7
NOT_WRITTEN = "run_manifest.json not written by this GIAP version"
DIGEST_LABEL = "manifest digest (SHA-256 of run_manifest.json)"


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim
    from giapsim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        giapsim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    return giapsim.simulate(cfg, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")


@pytest.fixture(scope="module")
def db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.mapping_steps import write_test_db

    return Path(write_test_db(tmp_path_factory.mktemp("test_db_m9_16b")))


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


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _prov(page_path: Path) -> tuple[dict[str, str], str]:
    page = parse(page_path.read_text(encoding="utf-8"))
    text = " ".join(page.section("provenance").text().split())
    rows: dict[str, str] = {}
    for table in page.root.by_attr("data-table-provenance", "run"):
        for tr in table.by_tag("tr"):
            cells = [c for c in tr.children if isinstance(c, Node) and c.tag in ("th", "td")]
            if len(cells) == 2:
                rows[cells[0].text().strip()] = cells[1].text().strip()
    return rows, text


def _check_full_provenance(out: Path, executor: str | None) -> None:
    m = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert m["status"] == "completed"
    rows, text = _prov(out / "report" / "report.html")
    assert NOT_WRITTEN not in text, "the published page predates the manifest"
    assert rows["GIAP version"] == m["giap_version"] == giap.__version__
    if m.get("git_commit"):
        assert rows["commit"].startswith(m["git_commit"])
    if executor is not None:
        assert f"--executor {executor}" in rows["command line"]
    assert rows["seed"] == str(SEED)
    assert rows[DIGEST_LABEL] == _sha(out / "run_manifest.json")
    assert rows["provenance digest"] == m["provenance_digest"]


def test_example_1_local_page_has_the_manifest_provenance(local: Path) -> None:
    _check_full_provenance(local, "local")
    work = local / "work" / "report" / "report.html"
    assert NOT_WRITTEN in _prov(work)[1], "the work copy keeps the in-graph page (behaviour 7)"


def _strings(obj: Any) -> Iterator[str]:
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield str(k)
            yield from _strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _strings(v)
    elif isinstance(obj, str):
        yield obj


def test_example_9_no_report_or_plot_checksum_in_the_manifest(local: Path) -> None:
    m = json.loads((local / "run_manifest.json").read_text(encoding="utf-8"))
    for s in _strings(m):
        assert "report.html" not in s, s
        assert not re.search(r"(^|/)plots/", s), s


def test_example_3_giap_report_after_the_run_is_byte_identical(
    local: Path, db: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    out = tmp_path_factory.mktemp("REPORT_16B") / "OUT"
    shutil.copytree(local, out, symlinks=True)
    tail = (out / "report" / "report.html").read_bytes()
    res = _invoke(["report", str(out)], db)
    assert res.exit_code == 0, res.output
    assert (out / "report" / "report.html").read_bytes() == tail


@pytest.mark.needs_nextflow
def test_example_2_nextflow_page_has_the_manifest_provenance(
    sim: Any, db: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    out = tmp_path_factory.mktemp("NF_16B") / "OUT"
    _run(sim, out, db, "--executor", "nextflow", "--profile", "farmshare_local")
    _check_full_provenance(out, None)
