"""``sample_name: auto`` is resolved by ``giap run`` (docs/specs/M7-04-variants-step-and-run.md
rule 4, examples 4 and 5, lead decision 1; GIAP_SPEC.md §5.2 "auto = stem of the peptide (or
protein) file", §7.5.5 "Sample column", ADR-0005 VCF-safe IDs).

Interface pinned here: ``giap.config.resolve_sample_name(cfg, peptides_path, proteins_path) ->
str``, a pure function. If the implementation lives in ``giap.cli``, ``giap.config`` re-exports
it.

* ``auto`` becomes the peptide file's name without ``.gz`` and its last suffix
  (``run_A.tsv.gz`` -> ``run_A``); with no peptides, the protein file's; with neither,
  ``sample``.
* The value, whether from ``auto`` or given explicitly, then maps to a VCF-safe ID: every
  character outside ``[A-Za-z0-9._-]`` becomes ``_``.
* ``giap run`` writes the resolved value into ``OUT/resolved_config.yaml`` before Nextflow
  starts. This is checked with a fake ``nextflow`` on ``PATH``, as in
  tests/unit/test_cli_run_nextflow.py, so nothing is launched. ``$GIAP_ROOT`` and
  ``$GIAP_DB_DIR`` are fakes there too.
"""

from __future__ import annotations

import os
import random
import re
import stat
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner
from helpers.fake_env import fake_giap_root, write_setup_manifest
from hypothesis import given
from hypothesis import strategies as st

from giap.config import GiapConfig
from giap.exitcodes import ExitCode

SAFE = re.compile(r"[A-Za-z0-9._-]+")


def _resolve(name: str, peptides: str | None, proteins: str | None) -> str:
    from giap.config import resolve_sample_name

    cfg = GiapConfig.model_validate({"sample_name": name})
    pep = Path(peptides) if peptides is not None else None
    prot = Path(proteins) if proteins is not None else None
    return resolve_sample_name(cfg, pep, prot)


# ----------------------------------------------------------------------------- pure function


@pytest.mark.parametrize(
    ("peptides", "proteins", "want"),
    [
        ("/data/run_A.tsv.gz", None, "run_A"),  # example 4
        ("/data/run_A.tsv", None, "run_A"),
        ("run_A.tsv.gz", "/data/prot.faa", "run_A"),  # peptides win over proteins
        ("/data/run.A.tsv", None, "run.A"),  # only the last suffix goes
        ("/data/report.parquet", None, "report"),
        ("/data/runB", None, "runB"),
        (None, "/data/ecoli_proteins.faa.gz", "ecoli_proteins"),
        (None, "/data/ecoli_proteins.faa", "ecoli_proteins"),
        (None, None, "sample"),
        ("/data/my run (2).tsv", None, "my_run__2_"),  # rule 4: VCF-safe
    ],
)
def test_auto_is_the_input_file_stem(peptides: str | None, proteins: str | None, want: str) -> None:
    assert _resolve("auto", peptides, proteins) == want


def test_example_5_spaces_map_to_underscores() -> None:
    assert _resolve("my sample", "/data/run_A.tsv.gz", None) == "my_sample"


@pytest.mark.parametrize("name", ["S1", "run_A", "strain-7.rep2"])
def test_a_safe_explicit_name_passes_unchanged(name: str) -> None:
    assert _resolve(name, "/data/run_A.tsv.gz", "/data/p.faa") == name


@given(st.text(min_size=1, max_size=40).filter(lambda s: s != "auto" and s.strip() == s and s))
def test_resolved_names_are_vcf_safe_and_keep_the_safe_characters(name: str) -> None:
    got = _resolve(name, None, None)
    assert SAFE.fullmatch(got), got
    assert len(got) == len(name)
    for a, b in zip(name, got, strict=True):
        assert b == (a if SAFE.fullmatch(a) else "_")


@given(st.from_regex(r"[A-Za-z0-9._ -]{1,20}", fullmatch=True).filter(lambda s: s != "auto"))
def test_resolving_twice_changes_nothing(name: str) -> None:
    once = _resolve(name, None, None)
    assert _resolve(once, None, None) == once


@given(st.from_regex(r"[A-Za-z0-9_-][A-Za-z0-9._-]{0,15}", fullmatch=True))
def test_auto_from_a_gz_peptide_file(stem: str) -> None:
    assert _resolve("auto", f"/x/{stem}.tsv.gz", None) == stem


def test_auto_is_never_returned() -> None:
    for pep, prot in ((None, None), ("/d/a.tsv", None), (None, "/d/b.faa")):
        assert _resolve("auto", pep, prot) != "auto"


# ----------------------------------------------------------------------------- giap run

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""


@pytest.fixture(autouse=True)
def _env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root = fake_giap_root(tmp_path / "giap_root")
    monkeypatch.setenv("GIAP_ROOT", str(root))
    db = tmp_path / "db"
    write_setup_manifest(db)
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))


@pytest.fixture
def evidence(tmp_path: Path) -> dict[str, Path]:
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    d = tmp_path / "in"
    d.mkdir()
    genome = d / "genome.fa"
    lines = [seq[i : i + 60] for i in range(0, len(seq), 60)]
    genome.write_text(">chr1 [topology=circular]\n" + "\n".join(lines) + "\n")
    aa = "ACDEFGHIKLMNPQRSTVWY"
    proteins = d / "ecoli_proteins.faa"
    proteins.write_text(
        "".join(f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n" for i in range(3))
    )
    import gzip

    peptides = d / "run_A.tsv.gz"
    rows = "".join(f"{''.join(rng.choice(aa) for _ in range(12))}K\t0.001\n" for _ in range(5))
    peptides.write_bytes(gzip.compress(("peptide\tq_value\n" + rows).encode(), mtime=0))
    return {"genome": genome, "proteins": proteins, "peptides": peptides}


def _run(argv: list[str]) -> tuple[int, str]:
    from giap.cli import cli

    result = CliRunner().invoke(cli, ["run", *argv])
    if result.exception is not None and not isinstance(result.exception, SystemExit):
        raise result.exception
    return result.exit_code, result.output


def _resolved(out: Path) -> str:
    data = yaml.safe_load((out / "resolved_config.yaml").read_text(encoding="utf-8"))
    return str(data["sample_name"])


def test_example_4_giap_run_resolves_auto_from_the_peptide_file(
    evidence: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(evidence["genome"]), "--peptides", str(evidence["peptides"]),
                       "--proteins", str(evidence["proteins"]), "--outdir", str(out)])  # fmt: skip
    assert code == ExitCode.SUCCESS, text
    assert _resolved(out) == "run_A"


def test_giap_run_resolves_auto_from_the_protein_file_without_peptides(
    evidence: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(evidence["genome"]), "--proteins",
                       str(evidence["proteins"]), "--outdir", str(out)])  # fmt: skip
    assert code == ExitCode.SUCCESS, text
    assert _resolved(out) == "ecoli_proteins"


def test_example_5_giap_run_maps_an_explicit_name(
    evidence: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(evidence["genome"]), "--peptides", str(evidence["peptides"]),
                       "--proteins", str(evidence["proteins"]), "--outdir", str(out),
                       "--set", "sample_name=my sample"])  # fmt: skip
    assert code == ExitCode.SUCCESS, text
    assert _resolved(out) == "my_sample"
