"""``giap makedb`` must take the genome (M3-24 lead decision 8; §15.2 step 2).

``build_strain_db`` called ``parse_bakta_json(opts.annotation)`` without ``opts.genome``, so the
recorded Bakta JSON of §15.2 step 1 - whose ``sequences`` carry no ``nt`` (R11: the nucleotides
are not committed twice) - could not be turned into ``validation/pilot/db_a2``. These tests pin
the fix: with ``--genome`` the stripped JSON builds exactly the database the full JSON builds,
and without it the error is a clean ``ValueError``/exit 2, not a traceback.

Decisions 13 and 14 (code review of PR #122) guard the genome against the annotation:

* the ``build.py`` check ``cid not in fasta or len(fasta[cid]) != rep.length`` is only reached
  when the JSON carries its own nucleotides, so the mismatch cases below use the **full** JSON
  (with a stripped JSON the parser raises first and ``build.py``'s check is never exercised);
* a sequence record with neither ``nt`` nor ``length`` makes ``rep.length`` come from the very
  FASTA it is compared against, so the check is a tautology and a truncated genome would build
  silently: such a record is refused outright, whatever the genome.

The last test rebuilds the committed ``db_a2`` from the recorded Nissle JSON and the genome of
decision 8 (``local/data/genomes/nissle1917/…``, uncommitted) and compares it byte for byte; it
is marked ``slow`` and skips when that FASTA is absent.
"""

from __future__ import annotations

import gzip
import json
import re
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
from makedb_oracle import revcomp

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


def _strip_nucleotides(path: Path, out: Path, *, lengths: bool = True) -> Path:
    """A recorded Bakta JSON: sequence records without ``nt`` (as under tests/data/recorded).

    ``lengths=False`` also drops the declared ``length`` (decision 13: then nothing is left to
    check the FASTA against).
    """
    data = json.loads(path.read_text())
    for record in data["sequences"]:
        record.pop("nt", None)
        if not lengths:
            record.pop("length", None)
    out.write_text(json.dumps(data, indent=1) + "\n")
    return out


def _read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif name:
            out[name] += line.strip()
    return out


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    chunks = []
    for name, seq in records.items():
        chunks.append(f">{name}\n")
        chunks += [seq[i : i + 80] + "\n" for i in range(0, len(seq), 80)]
    path.write_text("".join(chunks))
    return path


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
    """A FASTA of the wrong length is refused before anything is written (M3-23 behaviour).

    With a stripped JSON this is the parser's check (``annotate/bakta.py``); the cases below
    reach ``build.py``'s own check (decision 14).
    """
    contigs = [c for c in example.contigs]
    truncated = tmp_path / "short.fasta"
    text = write_fasta(tmp_path / "g.fasta", contigs).read_text().splitlines()
    truncated.write_text("\n".join([text[0], text[1][:60]]) + "\n")
    with pytest.raises(ValueError):
        _build(stripped_json, tmp_path / "DB", genome=truncated)


# ---------------------------------------------------------------- decisions 13 and 14


def _damaged(genome: Path, out: Path, how: str, cid: str) -> Path:
    """``genome`` with contig ``cid`` shortened, lengthened or renamed."""
    records = _read_fasta(genome)
    assert cid in records, (cid, sorted(records))
    if how == "short":
        records[cid] = records[cid][:-30]
    elif how == "long":
        records[cid] = records[cid] + "A" * 30
    elif how == "renamed":
        records = {("other" if k == cid else k): v for k, v in records.items()}
    else:  # pragma: no cover - guarded by the parametrisation
        raise AssertionError(how)
    return _write_fasta(out, records)


@pytest.mark.parametrize("how", ["short", "long", "renamed"])
def test_a_genome_that_contradicts_the_full_json_is_refused(
    full_json: Path, genome: Path, example: ExampleGenome, tmp_path: Path, how: str
) -> None:
    """Decision 14: the full JSON carries its own ``nt``, so ``build.py``'s check is reached.

    ``long`` kills a ``!=`` weakened to ``<``; ``renamed`` kills a dropped missing-contig half
    (a ``KeyError`` is not a ``ValueError``); ``short`` kills a dropped length half or a skipped
    block.
    """
    cid = example.contigs[0].id
    bad = _damaged(genome, tmp_path / f"{how}.fasta", how, cid)
    with pytest.raises(ValueError, match=re.escape(cid)):
        _build(full_json, tmp_path / "DB", genome=bad)
    assert not (tmp_path / "DB" / "strain_db.faa").exists(), "a partial database was written"


@pytest.mark.parametrize("how", ["match", "short", "long"])
def test_a_record_without_nt_and_without_length_is_refused(
    full_json: Path, genome: Path, example: ExampleGenome, tmp_path: Path, how: str
) -> None:
    """Decision 13: ``rep.length`` would come from the FASTA itself, so there is nothing to check.

    The refusal does not depend on the genome being wrong: with ``match`` the FASTA is correct
    and the build must still stop, because a truncated one could not be told apart from it.
    """
    annotation = _strip_nucleotides(full_json, tmp_path / "no_length.json", lengths=False)
    records = json.loads(annotation.read_text())["sequences"]
    assert all("length" not in r and "nt" not in r for r in records), records
    cid = example.contigs[0].id
    fasta = genome if how == "match" else _damaged(genome, tmp_path / f"{how}.fasta", how, cid)
    with pytest.raises(ValueError) as info:
        _build(annotation, tmp_path / "DB", genome=fasta)
    if how == "match":  # a correct FASTA cannot trip any other check: this is decision 13's
        message = str(info.value)
        assert cid in message, message
        # Decision 17: the fault is in the annotation, so the message must say so. "missing or of
        # the wrong length" would blame the genome file and would also pass a `cid in …` check.
        assert "declares no length" in message, message
        assert str(annotation) in message, message
    assert not (tmp_path / "DB" / "strain_db.faa").exists(), "a partial database was written"


@pytest.mark.parametrize("how", ["short", "long"])
def test_a_declared_length_still_catches_a_wrong_genome(
    stripped_json: Path, genome: Path, example: ExampleGenome, tmp_path: Path, how: str
) -> None:
    """Decision 13: a record with ``length`` but no ``nt`` (the recorded shape) stays checkable."""
    cid = example.contigs[0].id
    bad = _damaged(genome, tmp_path / f"{how}.fasta", how, cid)
    with pytest.raises(ValueError, match=re.escape(cid)):
        _build(stripped_json, tmp_path / "DB", genome=bad)
    assert not (tmp_path / "DB" / "strain_db.faa").exists(), "a partial database was written"


@pytest.mark.parametrize("shape", ["contigs-key", "sequence-key"])
def test_the_declared_length_follows_the_parsers_own_fallbacks(
    full_json: Path, genome: Path, tmp_path: Path, shape: str
) -> None:
    """Decision 17: both shapes ``annotate/bakta.py`` accepts must build, not be refused.

    A Bakta JSON may name the records ``contigs`` instead of ``sequences`` and the nucleotides
    ``sequence`` instead of ``nt``. Whatever reads the declared length (decision 13) has to know
    both: with either fallback missing every record looks length-less and a valid annotation is
    rejected as "declares no length".
    """
    data = json.loads(full_json.read_text())
    if shape == "contigs-key":
        data["contigs"] = data.pop("sequences")
        records = data["contigs"]
        assert all("length" in r and "nt" in r for r in records), records
    else:
        records = data["sequences"]
        for record in records:
            record["sequence"] = record.pop("nt")
            record.pop("length", None)  # a declared length would hide the nucleotide fallback
        assert all("length" not in r and "nt" not in r for r in records), records
    other = tmp_path / f"{shape}.json"
    other.write_text(json.dumps(data, indent=1) + "\n")
    want = _build(full_json, tmp_path / "BASE", genome=genome)
    got = _build(other, tmp_path / "SHAPE", genome=genome)
    for name in OUTPUTS:
        assert (got / name).read_bytes() == (want / name).read_bytes(), name


def test_the_genome_overrides_the_nucleotides_of_the_json(
    full_json: Path, genome: Path, example: ExampleGenome, tmp_path: Path
) -> None:
    """``--genome`` wins where both carry a sequence: the database is the FASTA's, not the JSON's.

    (Not a mismatch case: the variant has the same length, so no length check fires.)
    """
    records = _read_fasta(genome)
    cid = example.contigs[0].id
    variant = dict(records)
    variant[cid] = revcomp(records[cid])
    fasta = _write_fasta(tmp_path / "variant.fasta", variant)
    data = json.loads(full_json.read_text())
    for record in data["sequences"]:
        record["nt"] = variant[record["id"]]
    variant_json = tmp_path / "variant.json"
    variant_json.write_text(json.dumps(data, indent=1) + "\n")

    from_fasta = _build(full_json, tmp_path / "FASTA", genome=fasta, six_frame_orfs=True)
    from_json = _build(variant_json, tmp_path / "JSON", six_frame_orfs=True)
    untouched = _build(full_json, tmp_path / "ORIG", six_frame_orfs=True)
    for name in OUTPUTS:
        assert (from_fasta / name).read_bytes() == (from_json / name).read_bytes(), name
        assert (from_fasta / name).read_bytes() != (untouched / name).read_bytes(), (
            f"precondition: the variant genome must change {name}"
        )


@pytest.mark.parametrize("min_orf_aa", [0, -1])
def test_the_api_guard_on_min_orf_aa(full_json: Path, tmp_path: Path, min_orf_aa: int) -> None:
    """The API guard, not only the CLI one (a surviving mutant of the #122 review)."""
    with pytest.raises(ValueError, match="min_orf_aa"):
        _build(full_json, tmp_path / "DB", min_orf_aa=min_orf_aa)
    assert not (tmp_path / "DB" / "strain_db.faa").exists()


def test_the_recorded_nissle_json_declares_every_length() -> None:
    """Decision 13 cannot break ``db_a2``: all four Nissle records declare a length, none has nt."""
    with gzip.open(RECORDED_NISSLE, "rt") as fh:
        records = json.load(fh)["sequences"]
    assert len(records) == 4, len(records)
    for record in records:
        assert isinstance(record.get("length"), int) and record["length"] > 0, record.get("id")
        assert "nt" not in record, record.get("id")


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
