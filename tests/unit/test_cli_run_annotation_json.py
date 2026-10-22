# ruff: noqa: F811  (pytest fixtures imported from test_cli_run_wiring)
"""``giap run --annotation bakta.json`` is the precomputed Bakta layer (GIAP_SPEC.md §5.3: "Bakta
JSON: used as is, and Bakta is not run"). The MVP real-tool runs (Slurm 1864792, 1865149) re-ran
Bakta because the flag reached only preflight, so locus tags differed from the peptides' IDs."""

from __future__ import annotations

import json
from pathlib import Path

from unit.test_cli_run_wiring import (  # noqa: F401  (fixtures)
    _fake_giap_root,
    _run,
    db,
    fake_nf,
    genome,
    proteins,
)

from giap.exitcodes import ExitCode


def _precomputed(out: Path) -> dict[str, object]:
    params = json.loads((out / "work" / "params.json").read_text())
    return params["annotation"]["precomputed"]


def test_annotation_json_becomes_precomputed_bakta_json(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    ann = tmp_path / "in" / "bakta.json"
    ann.write_text("{}")
    out = tmp_path / "OUT"
    code, text = _run(
        ["--genome", str(genome), "--annotation", str(ann), "--outdir", str(out)], proteins
    )
    assert code == ExitCode.SUCCESS, text
    assert _precomputed(out)["bakta_json"] == str(ann.absolute())


def test_annotation_gff3_is_not_a_bakta_json(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    ann = tmp_path / "in" / "ann.gff3"
    ann.write_text("##gff-version 3\n")
    out = tmp_path / "OUT"
    code, text = _run(
        ["--genome", str(genome), "--annotation", str(ann), "--outdir", str(out)], proteins
    )
    assert code == ExitCode.SUCCESS, text
    assert _precomputed(out)["bakta_json"] is None


def test_set_wins_over_annotation_flag(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    ann = tmp_path / "in" / "bakta.json"
    ann.write_text("{}")
    other = tmp_path / "in" / "other.json"
    other.write_text("{}")
    out = tmp_path / "OUT"
    argv = ["--genome", str(genome), "--annotation", str(ann), "--outdir", str(out),
            "--set", f"annotation.precomputed.bakta_json={other}"]  # fmt: skip
    code, text = _run(argv, proteins)
    assert code == ExitCode.SUCCESS, text
    assert _precomputed(out)["bakta_json"] == str(other.absolute())
