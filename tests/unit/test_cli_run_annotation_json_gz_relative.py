# ruff: noqa: F811  (pytest fixtures imported from test_cli_run_wiring)
"""A relative ``--annotation`` Bakta JSON reaches Nextflow as an absolute path (PR #298 review
nit 2). Nextflow resolves a relative ``checkIfExists`` path against its launch directory, not
the user's working directory, so ``annotation.precomputed.bakta_json`` must be absolute and
name the file the user meant. Every other test passes absolute tmp paths.

A fake ``nextflow`` records its argv; no network or real database is touched.
"""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import pytest
from unit.test_cli_run_wiring import (  # noqa: F401  (fixtures)
    _fake_giap_root,
    _run,
    db,
    fake_nf,
    genome,
    proteins,
)

from giap.exitcodes import ExitCode

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "recorded"
    / "makedb"
    / "sim100k_plasmid_1.bakta.json"
)


@pytest.mark.parametrize("name", ["bakta.json", "bakta.json.gz", "BAKTA.JSON.GZ"])
def test_relative_annotation_becomes_absolute_precomputed_bakta_json(
    name: str,
    proteins: Path,
    genome: Path,
    db: Path,
    fake_nf: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ann = tmp_path / "ann" / name
    ann.parent.mkdir(parents=True)
    if name.lower().endswith(".gz"):
        with FIXTURE.open("rb") as src, gzip.open(ann, "wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copyfile(FIXTURE, ann)
    out = tmp_path / "OUT"
    monkeypatch.chdir(tmp_path)
    rel = Path("ann") / name
    code, text = _run(
        ["--genome", str(genome), "--annotation", str(rel), "--outdir", str(out)], proteins
    )
    assert code == ExitCode.SUCCESS, text
    params = json.loads((out / "work" / "params.json").read_text())
    got = params["annotation"]["precomputed"]["bakta_json"]
    assert isinstance(got, str), f"--annotation {rel} left precomputed.bakta_json {got!r}"
    assert Path(got).is_absolute(), f"precomputed.bakta_json {got!r} is relative"
    assert Path(got).resolve() == ann.resolve(), f"{got!r} is not {ann}"


@pytest.mark.parametrize("name", ["bakta.json", "bakta.json.gz", "BAKTA.JSON.GZ"])
def test_bakta_json_helper_returns_cwd_absolute_path(
    name: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``cli._bakta_json`` itself makes a relative Bakta JSON absolute against the cwd (the end
    to end test above also passes if a later stage absolutises the path)."""
    from giap.cli import _bakta_json

    monkeypatch.chdir(tmp_path)
    got = _bakta_json(Path("ann") / name)
    assert got == str(tmp_path / "ann" / name), got
