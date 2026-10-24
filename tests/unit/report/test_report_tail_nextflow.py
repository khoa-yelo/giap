"""M9-16b under the Nextflow executor (stubbed ``nextflow``): the tail calls ``rerender_report``
after ``run_tail`` on exit 0 and not on exit 4 (docs/specs/M9-16b-report-after-manifest.md
behaviours 1, 2; acceptance "test_manifest_cli_nextflow.py (extended)", kept in its own module
because the role guard forbids editing existing tests).

The fake ``nextflow`` (as tests/unit/manifest/test_manifest_cli_nextflow.py) exits ``$FAKE_NF_RC``
and publishes nothing, so ``OUT/report/`` is absent after a completed graph: the tail writes the
page fresh (degenerate row "``OUT/report/report.html`` absent"). ``rerender_report`` is looked up
on ``giap.steps.plots`` at call time. ``$GIAP_ROOT`` and ``$GIAP_DB_DIR`` are fakes from
``helpers.fake_env``.
"""

from __future__ import annotations

import hashlib
import json
import os
import random
import stat
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_giap_root, write_setup_manifest
from helpers.report_fixtures import parse, swap_registry

from giap.exitcodes import ExitCode

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit "${FAKE_NF_RC:-0}"
"""


@pytest.fixture
def nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    monkeypatch.setenv("GIAP_ROOT", str(fake_giap_root(tmp_path / "giap_root")))
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    for var in ("SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "FAKE_NF_RC"):
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))
    rng = random.Random(3)
    genome = tmp_path / "in" / "genome.fa"
    genome.parent.mkdir()
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    genome.write_text(">chr1 [topology=circular]\n" + "\n".join(seq[i:i + 60]
                      for i in range(0, 6000, 60)) + "\n")  # fmt: skip
    aa = "ACDEFGHIKLMNPQRSTVWY"
    proteins = tmp_path / "in" / "proteins.faa"
    proteins.write_text("".join(f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n"
                                for i in range(3)))  # fmt: skip
    swap_registry(monkeypatch)
    return {"genome": genome, "proteins": proteins}


@pytest.fixture
def spy(monkeypatch: pytest.MonkeyPatch) -> list[bytes | None]:
    """Wrap ``giap.steps.plots.rerender_report``; record the manifest bytes at call time."""
    from giap.steps import plots as P

    calls: list[bytes | None] = []
    real = getattr(P, "rerender_report", None)

    def wrapper(outdir: Path, cfg: Any) -> Any:
        man = Path(outdir) / "run_manifest.json"
        calls.append(man.read_bytes() if man.is_file() else None)
        assert real is not None
        return real(outdir, cfg)

    monkeypatch.setattr(P, "rerender_report", wrapper, raising=False)
    return calls


def _run(nf: dict[str, Path], out: Path) -> tuple[int, str]:
    from giap.cli import cli

    argv = ["run", "--genome", str(nf["genome"]), "--proteins", str(nf["proteins"]),
            "--outdir", str(out)]  # fmt: skip
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def test_exit_0_rerenders_after_the_manifest(
    nf: dict[str, Path], tmp_path: Path, spy: list[bytes | None]
) -> None:
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    assert len(spy) == 1, "the Nextflow tail did not call rerender_report"
    assert spy[0] is not None, "rerender_report ran before run_manifest.json was written"
    assert json.loads(spy[0])["status"] == "completed"
    manifest = out / "run_manifest.json"
    assert manifest.read_bytes() == spy[0], "the manifest was rewritten after the re-render"
    page = out / "report" / "report.html"
    assert page.is_file(), "the tail writes the page even with no in-graph copy"
    text = " ".join(parse(page.read_text(encoding="utf-8")).section("provenance").text().split())
    assert hashlib.sha256(manifest.read_bytes()).hexdigest() in text


def test_exit_4_does_not_rerender(
    nf: dict[str, Path], tmp_path: Path, spy: list[bytes | None], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_NF_RC", "1")
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.STEP_FAILED, log
    assert json.loads((out / "run_manifest.json").read_text())["status"] == "failed"
    assert spy == []
    assert not (out / "report" / "report.html").exists()
    assert not (out / "logs" / "report_tail.log").exists()
