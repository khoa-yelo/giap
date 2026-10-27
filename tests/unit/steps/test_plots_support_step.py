"""The ``plots`` step writes ``plots/support/`` from its run inputs; ``gippa report`` never does
(M9-03 acceptance "support written when ``genome`` is supplied, absent otherwise; ``gippa report``
never writes ``plots/support/``", degenerate case "``genome`` input absent: no ``gc_skew.tsv``;
landmarks still written", lead decisions 2026-10-22: GC-skew window 10 kb as a module constant,
``gippa report`` never recomputes GC skew).

A separate file from tests/unit/steps/test_plots_step.py (M9-02's tests are not edited here).

Interface assumed beyond the spec text (M9-03 test author):

- Inputs as the other steps take them: ``genome`` a FASTA path (headers ``>id [topology=...]``),
  ``features`` and ``replicons`` the C15 / C1 parquet files.
- Landmarks need ``replicons`` (``features`` absent behaves like no oriC call); GC skew needs
  ``genome`` and ``replicons``. Without ``replicons`` no support file is written.
- The support tables are written **before** rendering, and the plots rendered by the step see
  them through ``view.support(name)`` (V1/V2 need them during ``gippa run``).
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner
from helpers.report_fixtures import cfg, register_fake, run_out, swap_registry

from gippa.steps import registry as reg

FEATURE_COLUMNS = [
    "feature_id", "type", "contig", "start", "end", "strand",
    "gene", "product", "amino_acid", "anticodon",
]  # fmt: skip


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _inputs_dir(tmp: Path) -> dict[str, Path]:
    """chr1 25,000 bp circular with one oriC call; plasmid_1 4,000 bp circular, none."""
    d = tmp / "inputs"
    d.mkdir(parents=True, exist_ok=True)
    chr1 = ("GGGC" + "A" * 9_996) + ("ACGT" * 2_500) + ("C" * 5_000)
    plasmid = "G" * 4_000
    fasta = d / "genome.fasta"
    fasta.write_text(
        f">chr1 [topology=circular]\n{chr1}\n>plasmid_1 [topology=circular]\n{plasmid}\n",
        encoding="utf-8",
    )
    reps = pd.DataFrame(
        [("chr1", len(chr1), "circular", "chromosome"),
         ("plasmid_1", len(plasmid), "circular", "plasmid")],
        columns=["contig", "length", "topology", "kind"],
    ).astype({"length": "int64"})  # fmt: skip
    reps.to_parquet(d / "replicons.parquet", index=False)
    feats = pd.DataFrame(
        [("F1", "oric", "chr1", 900, 1_100, ".", None, "oriC", None, None)],
        columns=FEATURE_COLUMNS,
    ).astype({"start": "int64", "end": "int64"})
    feats.to_parquet(d / "features.parquet", index=False)
    return {
        "genome": fasta,
        "replicons": d / "replicons.parquet",
        "features": d / "features.parquet",
    }


def _run(tmp: Path, out: Path, extra: dict[str, Path]) -> Path:
    plots = tmp / "work" / "plots"
    inputs = {"preflight": out / "preflight" / "preflight.json", **extra}
    reg.run_step("plots", cfg(), inputs, {"plots": plots})
    return plots


def _support(plots: Path) -> list[str]:
    d = plots / "support"
    return sorted(p.name for p in d.iterdir()) if d.is_dir() else []


def _tsv(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]


def test_step_writes_both_support_tables_with_the_genome(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    plots = _run(tmp_path, run_out(tmp_path), _inputs_dir(tmp_path))
    assert _support(plots) == ["gc_skew.tsv", "replicon_landmarks.tsv"]

    lm = _tsv(plots / "support" / "replicon_landmarks.tsv")
    assert lm[0] == ["contig", "length_bp", "topology", "oric_position", "ter_position", "status"]
    assert lm[1] == ["chr1", "25000", "circular", "1001", "13501", "called"]
    assert lm[2] == ["plasmid_1", "4000", "circular", "NA", "NA", "no_oric"]

    gs = _tsv(plots / "support" / "gc_skew.tsv")
    assert gs[0] == ["contig", "start", "end", "gc_skew", "cumulative_gc_skew"]
    # the 10 kb module window (lead decision 1): chr1 1-10,000 / 10,001-20,000 / 20,001-25,000
    assert [(r[0], r[1], r[2]) for r in gs[1:]] == [
        ("chr1", "1", "10000"),
        ("chr1", "10001", "20000"),
        ("chr1", "20001", "25000"),
        ("plasmid_1", "1", "4000"),
    ]
    assert [float(r[3]) for r in gs[1:]] == pytest.approx([0.5, 0.0, -1.0, 1.0])
    assert [float(r[4]) for r in gs[1:]] == pytest.approx([0.5, 0.5, -0.5, 1.0])


def test_step_without_genome_writes_landmarks_only(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    ins = _inputs_dir(tmp_path)
    del ins["genome"]
    plots = _run(tmp_path, run_out(tmp_path), ins)
    assert _support(plots) == ["replicon_landmarks.tsv"]


def test_step_without_features_marks_every_circular_replicon_no_oric(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    ins = _inputs_dir(tmp_path)
    del ins["features"]
    plots = _run(tmp_path, run_out(tmp_path), ins)
    assert "replicon_landmarks.tsv" in _support(plots)
    lm = _tsv(plots / "support" / "replicon_landmarks.tsv")
    assert [r[5] for r in lm[1:]] == ["no_oric", "no_oric"]


def test_step_without_run_inputs_writes_no_support(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    plots = _run(tmp_path, run_out(tmp_path), {})
    assert _support(plots) == []


def test_step_support_is_byte_identical_across_runs(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = run_out(tmp_path)
    ins = _inputs_dir(tmp_path)
    a = tmp_path / "a"
    b = tmp_path / "b"
    for d in (a, b):
        reg.run_step(
            "plots", cfg(), {"preflight": out / "preflight" / "preflight.json", **ins},
            {"plots": d},
        )  # fmt: skip
    for name in ("gc_skew.tsv", "replicon_landmarks.tsv"):
        assert (a / "support" / name).read_bytes() == (b / "support" / name).read_bytes()


def test_plots_rendered_by_the_step_see_the_support_tables(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    from gippa.viz import base

    seen: dict[str, Any] = {}

    @base.plot("V90", "support_probe")
    def _probe(view: Any, c: Any) -> Any:
        seen["landmarks"] = view.support("replicon_landmarks")
        seen["skew"] = view.support("gc_skew")
        return base.Skip("probe")

    _run(tmp_path, run_out(tmp_path), _inputs_dir(tmp_path))
    assert seen["landmarks"] is not None and list(seen["landmarks"]["contig"]) == [
        "chr1",
        "plasmid_1",
    ]
    assert seen["skew"] is not None and len(seen["skew"]) == 4


# ---------------------------------------------------------------- gippa report


def _report(out: Path) -> Any:
    from gippa.cli import cli

    res = CliRunner().invoke(
        cli, ["report", str(out), "--set", "report.formats=[svg]", "--set", "report.dpi=50"]
    )
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _tree(path: Path) -> dict[str, bytes]:
    return {p.relative_to(path).as_posix(): p.read_bytes() for p in path.rglob("*") if p.is_file()}


def test_gippa_report_without_support_writes_none(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    register_fake("V07")
    out = run_out(tmp_path, islands=1)
    # run inputs lying next to OUT (and a genome in OUT) are never used to recompute (decision 2)
    ins = _inputs_dir(tmp_path)
    (out / "genome.fasta").write_bytes(ins["genome"].read_bytes())
    res = _report(out)
    assert res.exit_code == 0, res.output
    assert not (out / "plots" / "support").exists()


def test_gippa_report_reads_but_never_rewrites_support(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    from gippa.io.plot_support import write_support
    from gippa.viz.support import replicon_landmarks

    register_fake("V07")
    out = run_out(tmp_path, islands=1)
    reps = pd.DataFrame(
        [("chr1", 10_000, "circular", "chromosome")],
        columns=["contig", "length", "topology", "kind"],
    )
    write_support("replicon_landmarks", replicon_landmarks(reps, None), out / "plots")
    before = _tree(out / "plots" / "support")
    stamps = {p: p.stat().st_mtime_ns for p in (out / "plots" / "support").iterdir()}
    res = _report(out)
    assert res.exit_code == 0, res.output
    assert _tree(out / "plots" / "support") == before
    assert {p: p.stat().st_mtime_ns for p in (out / "plots" / "support").iterdir()} == stamps
