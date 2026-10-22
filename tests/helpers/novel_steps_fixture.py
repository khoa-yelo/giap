"""Shared inputs for the M7-08 tests (docs/specs/M7-08-novel-steps-and-second-round.md).

* ``SIM_CFG``: giapsim ``SIM100K`` with M3-19 proteomics, planted novel ORFs of the three M3-18
  classes and planted N-terminal extensions (example 1: "sim100k with planted novel ORFs").
* ``run_giap``: ``giap run`` through the click runner (the local executor by default), with a test
  ``$GIAP_DB_DIR`` (never the real one).
* ``library_orfs``: the M7-06 library on the run's own work files (``candidate_rows`` +
  ``reconstruct``), the oracle of what ``novel_candidates`` wraps.

Conventions pinned by the test author (spec rule 2, ADR-0024, M7-07 lead decision 1):

* ``rbs_novel`` is keyed by ``(contig, strand, start)``; ``start`` is the start codon's **first
  base in transcription direction** (on ``-`` its highest coordinate), 0-based in
  ``rbs_novel.parquet`` and in ``load_precomputed('rbs_novel', ...)``, 1-based in the TSV layer --
  the convention of ``start_longest`` / ``start_best_ostir`` in ``novel_orfs.tsv``.
  ``Orf.starts`` (and ``start_best_ostir`` in the M7-07 ``scored`` frame) keep the codon's
  *lowest* coordinate; ``tpos`` converts.
* the simulator's ``rbs_novel`` layer path is read from its ``config.yaml``
  (``annotation.precomputed.rbs_novel``), so ``SimResult.precomputed`` may keep its M2-05 keys.
"""

from __future__ import annotations

import dataclasses
import os
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pandas as pd

SEED = 7
NOVEL_COUNTS = (("intergenic", 3), ("cds_outframe", 2), ("antisense", 2))
N_NTERM = 2
TABLES = ("novel_orfs.tsv", "novel_orfs.gff3", "novel_orfs.faa", "annotation_patch.gff3")
PUBLISHED = {
    "novel_orfs": "tables/novel_orfs.tsv",
    "novel_orfs_gff3": "tables/novel_orfs.gff3",
    "novel_orfs_faa": "tables/novel_orfs.faa",
    "annotation_patch": "tables/annotation_patch.gff3",
}
OSTIR_HEADER = ("name,start_codon,start_position,expression,RBS_distance_bp,dG_total,"
                "dG_rRNA:mRNA,dG_mRNA,dG_spacing,dG_standby,dG_start_codon")  # fmt: skip
OSTIR_INPUT_HEADER = "id,seq,start,end,anti-Shine-Dalgarno"
PROMOTECH_HEADER = "CHROM\tSEQ\tPRED"
RBS_NOVEL_COLUMNS = ("contig", "strand", "start", "tir_ostir", "sd_dg", "mrna_dg", "spacing_dg",
                     "standby_dg", "spacer_bp", "start_codon")  # fmt: skip


def sim_config(**over: Any) -> Any:
    import giapsim
    from giapsim.proteomics import ProteomicsConfig

    kw: dict[str, Any] = {"proteomics": ProteomicsConfig(), "n_novel_orfs": NOVEL_COUNTS,
                          "n_nterm_extensions": N_NTERM}  # fmt: skip
    kw.update(over)
    return dataclasses.replace(giapsim.SIM100K, **kw)


def simulate(out: Path, **over: Any) -> Any:
    import giapsim

    return giapsim.simulate(sim_config(**over), giapsim.SIM100K_SEED, out)


@contextmanager
def db_env(root: Path) -> Iterator[Path]:
    from helpers.mapping_steps import write_test_db

    old = os.environ.get("GIAP_DB_DIR")
    db = write_test_db(root / "test_db")
    os.environ["GIAP_DB_DIR"] = str(db)
    try:
        yield db
    finally:
        if old is None:
            os.environ.pop("GIAP_DB_DIR", None)
        else:
            os.environ["GIAP_DB_DIR"] = old


def run_giap(sim: Any, out: Path, *extra: str, executor: str = "local") -> tuple[int, str]:
    """``giap run`` on the simulated genome, peptides, Bakta JSON and config (exit code, output)."""
    from click.testing import CliRunner

    from giap.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides),
            "--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml),
            "--outdir", str(out), "--modules", "core", "--seed", str(SEED),
            "--executor", executor, *extra]  # fmt: skip
    res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    err = res.stderr if res.stderr_bytes is not None else ""
    return res.exit_code, res.output + err


def rbs_novel_layer(sim: Any) -> Path:
    """The simulator's ``annotation.precomputed.rbs_novel`` (from its ``config.yaml``)."""
    from giap.config import resolve_config

    cfg = resolve_config(Path(sim.config_yaml))
    path = cfg.annotation.precomputed.rbs_novel  # type: ignore[attr-defined]
    assert path is not None, "the simulator's config.yaml sets no annotation.precomputed.rbs_novel"
    p = Path(path)
    return p if p.is_absolute() else Path(sim.config_yaml).parent / p


def genome_of(sim: Any) -> dict[str, str]:
    from giap.io.fasta import read_genome_fasta

    return {r.id: r.sequence.upper() for r in read_genome_fasta(Path(sim.genome_fasta))}


def annotation_of(sim: Any) -> Any:
    from giap.annotate.bakta import parse_bakta_json

    return parse_bakta_json(Path(sim.bakta_json))


def library_orfs(work: Path, sim: Any, mapping: Path | None = None) -> list[Any]:
    """M7-06 ``candidate_rows`` + ``reconstruct`` on the run's peptide hits and C9 peptides."""
    from giap.config import GiapConfig
    from giap.io.work_tables import _read_tsv_table
    from giap.novel.candidates import candidate_rows, reconstruct

    hits = _read_tsv_table("peptide_hits", (mapping or work / "mapping") / "peptide_hits.tsv")
    pep = pd.read_parquet(work / "peptides")
    ann = annotation_of(sim)
    return reconstruct(candidate_rows(hits, pep, None), ann.model, dict(ann.sequences),
                       GiapConfig())  # fmt: skip


def tpos(strand: str, lowest: int) -> int:
    """First base (transcription direction) of the codon whose lowest coordinate is ``lowest``."""
    return lowest if strand == "+" else lowest + 2


_COMP = str.maketrans("ACGT", "TGCA")


def codon_at(seq: str, strand: str, first: int, circular: bool = True) -> str | None:
    """The codon whose first base in transcription direction is ``first`` (0-based)."""
    n = len(seq)
    idx = [first + i for i in range(3)] if strand == "+" else [first - i for i in range(3)]
    if not circular and any(not 0 <= i < n for i in idx):
        return None
    nt = "".join(seq[i % n] for i in idx)
    return nt if strand == "+" else nt.translate(_COMP)


def step_inputs(
    name: str, work: Path, override: Mapping[str, Any] | None = None, skip: tuple[str, ...] = ()
) -> dict[str, Any]:
    """Every declared input of ``name`` found under ``work`` (required ones must exist), with
    ``override`` winning; keys in ``skip`` are left out."""
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get(name)
    override = dict(override or {})
    out: dict[str, Any] = {}
    for key in (*spec.inputs, *spec.optional_inputs):
        if key in skip:
            continue
        if key in override:
            out[key] = override[key]
        elif (work / key).exists():
            out[key] = work / key
        elif key in spec.inputs:
            raise AssertionError(f"{name}: required input {key} not under {work}")
    return out


def run(name: str, cfg: Any, inputs: Mapping[str, Any], outdir: Path) -> dict[str, Path]:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    outs = {k: outdir / k for k in reg.get(name).outputs}
    reg.run_step(name, cfg, dict(inputs), outs)
    return outs


def write_ostir_csv(path: Path, rows: list[list[Any]]) -> Path:
    lines = [OSTIR_HEADER] + [",".join(str(v) for v in r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def ostir_out_row(name: str, codon: str, pos: int, tir: float, spacer: int = 7) -> list[Any]:
    """One OSTIR 1.1.3 output row; dG columns derived from ``tir`` so each row is distinct."""
    return [name, codon, pos, f"{tir:.4f}", spacer, "1.5", f"{-tir / 1000:.4f}", "-10.5",
            "0.25", "-0.5", "-1.194"]  # fmt: skip


def tsv_dicts(path: Path) -> list[dict[str, str | None]]:
    from giap.io.tsv import read_tsv

    header, rows = read_tsv(path)
    return [dict(zip(header, r, strict=True)) for r in rows]


SCHEMA = Path(__file__).resolve().parents[2] / "src/giap/schemas/tables/novel_orfs.schema.json"


def assert_valid_novel_orfs(path: Path) -> None:
    """``novel_orfs.tsv`` validates against its Frictionless schema (M7-07)."""
    import json

    import pytest

    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in json.loads(SCHEMA.read_text(encoding="utf-8")).items()
         if k != "coordinates"}  # fmt: skip
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])
