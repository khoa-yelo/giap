"""``giap makedb`` must take the genome (M3-24 lead decision 8; §15.2 step 2).

``build_strain_db`` called ``parse_bakta_json(opts.annotation)`` without ``opts.genome``, so the
recorded Bakta JSON of §15.2 step 1 - whose ``sequences`` carry no ``nt`` (R11: the nucleotides
are not committed twice) - could not be turned into ``validation/pilot/db_a2``. These tests pin
the fix: with ``--genome`` the stripped JSON builds exactly the database the full JSON builds,
and without it the error is a clean ``ValueError``/exit 2, not a traceback.

The last test rebuilds the committed ``db_a2`` from the recorded Nissle JSON and the genome of
decision 8 (``local/data/genomes/nissle1917/…``, uncommitted) and compares it byte for byte; it
is marked ``slow`` and skips when that FASTA is absent.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_kit import ID_MAP, MANIFEST, NISSLE_GENOME, RECORDED_NISSLE, STRAIN_DB

from giap.exitcodes import ExitCode

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from makedb_genomes import (  # tests/helpers is on sys.path, as in test_makedb.py
    ExampleGenome,
    build_example,
    write_bakta_json,
    write_fasta,
)

OUTPUTS = ("strain_db.faa", "db_id_map.tsv")


def _build(annotation: Path, outdir: Path, **kw: Any) -> Path:
    from giap.makedb.build import MakedbOptions, build_strain_db

    genome = kw.pop("genome", None)
    build_strain_db(MakedbOptions(annotation=annotation, genome=genome, outdir=outdir, **kw))
    return outdir


def _cli(argv: list[str]) -> int:
    from giap.cli import main

    with pytest.raises(SystemExit) as info:
        main(argv)
    return 0 if info.value.code is None else int(info.value.code)


def _strip_nucleotides(path: Path, out: Path) -> Path:
    """A recorded Bakta JSON: sequence records without ``nt`` (as under tests/data/recorded)."""
    data = json.loads(path.read_text())
    for record in data["sequences"]:
        record.pop("nt", None)
    out.write_text(json.dumps(data, indent=1) + "\n")
    return out


@pytest.fixture(scope="module")
def example() -> ExampleGenome:
    return build_example()


@pytest.fixture
def full_json(tmp_path: Path, example: ExampleGenome) -> Path:
    return write_bakta_json(tmp_path / "bakta.json", example.contigs, example.features)


@pytest.fixture
def stripped_json(tmp_path: Path, full_json: Path) -> Path:
    return _strip_nucleotides(full_json, tmp_path / "bakta_recorded.json")


@pytest.fixture
def genome(tmp_path: Path, example: ExampleGenome) -> Path:
    return write_fasta(tmp_path / "genome.fasta", example.contigs)


def test_a_json_without_nucleotides_needs_the_genome(stripped_json: Path, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=r"no nucleotides|no genome"):
        _build(stripped_json, tmp_path / "DB")
    assert not (tmp_path / "DB" / "strain_db.faa").exists(), "a partial database was written"


def test_the_cli_reports_it_as_invalid_input(
    stripped_json: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    argv = ["makedb", "--annotation", str(stripped_json), "--outdir", str(tmp_path / "DB"),
            "--no-contaminants"]  # fmt: skip
    assert _cli(argv) == ExitCode.INVALID_INPUT
    assert "giap makedb:" in capsys.readouterr().err


@pytest.mark.parametrize("six_frame_orfs", [False, True], ids=["proteins", "six-frame"])
def test_the_genome_replaces_the_missing_nucleotides(
    full_json: Path, stripped_json: Path, genome: Path, tmp_path: Path, six_frame_orfs: bool
) -> None:
    """Decision 8: ``opts.genome`` reaches ``parse_bakta_json``, so both builds agree exactly."""
    want = _build(full_json, tmp_path / "FULL", six_frame_orfs=six_frame_orfs)
    kw = {"genome": genome, "six_frame_orfs": six_frame_orfs}
    got = _build(stripped_json, tmp_path / "RECORDED", **kw)
    for name in OUTPUTS:
        assert (got / name).read_bytes() == (want / name).read_bytes(), name


def test_the_genome_must_match_the_annotation(
    stripped_json: Path, example: ExampleGenome, tmp_path: Path
) -> None:
    """A FASTA of the wrong length is refused before anything is written (M3-23 behaviour)."""
    contigs = [c for c in example.contigs]
    truncated = tmp_path / "short.fasta"
    text = write_fasta(tmp_path / "g.fasta", contigs).read_text().splitlines()
    truncated.write_text("\n".join([text[0], text[1][:60]]) + "\n")
    with pytest.raises(ValueError):
        _build(stripped_json, tmp_path / "DB", genome=truncated)


def test_the_cli_builds_the_db_a2_command_of_step_2(
    stripped_json: Path, genome: Path, tmp_path: Path
) -> None:
    """§15.2 step 2 as decision 8 rewrites it, with ``--genome`` and ``--outdir``."""
    out = tmp_path / "db_a2"
    argv = ["makedb", "--annotation", str(stripped_json), "--genome", str(genome),
            "--outdir", str(out), "--no-contaminants"]  # fmt: skip
    assert _cli(argv) == ExitCode.SUCCESS
    assert {p.name for p in out.iterdir()} >= {*OUTPUTS, "makedb_manifest.json"}
    manifest = json.loads((out / "makedb_manifest.json").read_text())
    assert set(manifest["inputs"]) == {"annotation", "genome"}, manifest["inputs"]


@pytest.mark.slow
def test_the_committed_db_a2_comes_from_the_recorded_nissle_json(tmp_path: Path) -> None:
    """Acceptance criterion: ``strain_db.faa.gz`` was produced by ``giap makedb`` from that JSON."""
    if not NISSLE_GENOME.is_file():
        pytest.skip(f"{NISSLE_GENOME} is not present (uncommitted, R11)")
    out = _build(RECORDED_NISSLE, tmp_path / "db_a2", genome=NISSLE_GENOME)
    rebuilt = (out / "strain_db.faa").read_bytes()
    assert rebuilt.count(b"\n>") + 1 == sum(
        json.loads((out / "makedb_manifest.json").read_text())["counts"].values()
    )
    assert STRAIN_DB.is_file(), f"{STRAIN_DB} does not exist (M3-24)"
    with gzip.open(STRAIN_DB, "rb") as fh:
        committed = fh.read()
    assert rebuilt == committed, (
        "validation/pilot/db_a2/strain_db.faa.gz differs from a fresh `giap makedb` on "
        f"{RECORDED_NISSLE.name} + {NISSLE_GENOME.name}"
    )
    assert (out / "db_id_map.tsv").read_bytes() == ID_MAP.read_bytes()
    fresh = json.loads((out / "makedb_manifest.json").read_text())
    committed_manifest = json.loads(MANIFEST.read_text())
    assert committed_manifest["counts"] == fresh["counts"]
    assert committed_manifest["options"] == fresh["options"]
