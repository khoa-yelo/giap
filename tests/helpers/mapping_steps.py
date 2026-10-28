"""Shared helpers for the M3-17 mapping-step tests (docs/specs/M3-17-mapping-steps.md).

Conventions pinned by the tests (the spec names files, not keys; see each test module docstring):

* Root work keys (inputs no mapping step produces): ``genome`` (``genome.norm.fa``, the key
  ``parse_bakta`` already uses), ``peptides`` (``peptides.norm.parquet``, the M4 loader key),
  ``bakta_json`` (until the ext ``bakta`` step supplies it, as in the operon tests); the proteins
  key is whatever the ext ``miniprot`` step declares besides ``genome`` (``proteins.long.faa``).
* Intermediate keys follow the §6.2 file names: ``sixframe`` (``.npz``, ``gippa.mapping.sixframe
  .save_index`` format), ``hits_exact`` and ``hits_tolerant`` (Parquet, ``gippa.mapping.hits``
  columns), ``placement`` (Parquet, one row per placed query with ``peptide``, ``placement_q``).
* The seven §7.2.9 Module 2 files keep their §7.2.9 file names. They are located either as one
  directory (any output key or ``OUT`` subdirectory) holding all seven, or as output keys whose
  names equal the file names without ``.tsv`` and with ``.`` → ``_`` (``mappings_bed``).

``drive`` runs the mapping chain through ``gippa.steps.registry.run_step`` like the local executor,
but with *injection*: an external step listed in ``inject`` is replaced by recorded files for its
output keys (miniprot / tblastn fixtures under ``tests/data/recorded``); any other external step
gets empty files when its outputs are needed (a tool that found nothing). Tests never start a tool.
"""

from __future__ import annotations

import importlib
import json
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner

MAPPING_PY_STEPS: tuple[str, ...] = (
    "sixframe_index",
    "map_exact",
    "map_tolerant",
    "placement",
    "loci",
)
MAPPING_EXT_STEPS: tuple[str, ...] = ("miniprot", "tblastn")
MODULE2_FILES: tuple[str, ...] = (
    "peptide_hits.tsv",
    "protein_hits.tsv",
    "mappings.probed",
    "mappings.bed",
    "mappings.gff3",
    "unmapped.tsv",
    "mapping_fdr.tsv",
)
MAPPING_TSVS: tuple[str, ...] = (
    "peptide_hits.tsv",
    "protein_hits.tsv",
    "unmapped.tsv",
    "mapping_fdr.tsv",
)
REPO = Path(__file__).resolve().parents[2]
CONTAMINANTS_MINI = REPO / "tests/data/contaminants_mini.faa"
RECORDED_MINIPROT = REPO / "tests/data/recorded/miniprot/0.18/sim100k"
RECORDED_TBLASTN = REPO / "tests/data/recorded/blast/2.16.0/sim100k"


# ---------------------------------------------------------------- registry


def fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    """A fresh registry holding every real step (robust to earlier temporary registries)."""
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":  # reloading the registry itself would drop the patched dict
            importlib.reload(importlib.import_module(f"gippa.steps.{n}"))
    return reg


def mapping_module() -> ModuleType:
    return importlib.import_module("gippa.steps.mapping")


# ---------------------------------------------------------------- simulator + preflight


@dataclass(frozen=True)
class SimRun:
    """A simulated genome with proteomics and its preflight work files."""

    res: Any  # gippasim SimResult
    genome: Path  # genome.norm.fa
    peptides: Path  # peptides.norm.parquet
    proteins_long: Path | None  # proteins.long.faa (only with --proteins)
    config_yaml: Path


def write_test_db(db: Path, contaminants: Path | None = None) -> Path:
    """A test ``$GIPPA_DB_DIR`` (M3-17 review decision 12: tests never read the real one).

    ``db_manifest.json`` in the per-name form ``read_manifest`` accepts (``{"<name>": {"version":
    ...}}``; the ``gippa setup`` item form cannot name the ``genomad``/``antismash`` binaries)
    lists every requirement of every module at its pinned version (``1.0`` when unpinned). With
    ``contaminants`` the library is copied to ``contaminants/contaminants.fasta`` (M3-15) and the
    manifest gets a ``contaminants`` entry with that ``path``.
    """
    from gippa.preflight.environment import _ROWS

    db.mkdir(parents=True, exist_ok=True)
    doc: dict[str, dict[str, str]] = {
        name: {"version": version or "1.0"}
        for rows in _ROWS.values()
        for name, _kind, version in rows
    }
    if contaminants is not None:
        rel = "contaminants/contaminants.fasta"
        (db / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(contaminants, db / rel)
        doc["contaminants"] = {"version": "bddf6da9", "path": rel}
    text = json.dumps(dict(sorted(doc.items())), indent=2) + "\n"
    (db / "db_manifest.json").write_text(text, encoding="utf-8")
    return db


def preflight(res: Any, out: Path, proteins: bool = False) -> SimRun:
    """``gippa preflight`` on the simulator's genome, peptides and Bakta-like JSON (with a test
    ``$GIPPA_DB_DIR`` next to ``out``)."""
    from gippa.cli import cli

    args = ["preflight", "--genome", str(res.genome_fasta), "--peptides", str(res.peptides)]
    args += ["--annotation", str(res.bakta_json), "--outdir", str(out)]
    if res.config_yaml is not None:
        args += ["--config", str(res.config_yaml)]
    if proteins:
        args += ["--proteins", str(res.proteins_faa)]
    db = write_test_db(out.with_name(out.name + "_db"))
    r = CliRunner().invoke(cli, args, env={"GIPPA_DB_DIR": str(db)})
    assert r.exit_code == 0, r.output
    work = out / "work" / "preflight"
    long = work / "proteins.long.faa"
    return SimRun(
        res,
        work / "genome.norm.fa",
        work / "peptides.norm.parquet",
        long if long.exists() else None,
        res.config_yaml,
    )


def simulate_with_proteomics(
    out: Path,
    seed: int,
    length: int = 50_000,
    proteomics: Any | None = None,
    proteins: bool = False,
) -> SimRun:
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SimConfig, simulate

    cfg = SimConfig(
        genome_length=length,
        proteomics=proteomics or ProteomicsConfig(),
        write_precomputed=True,
    )
    res = simulate(cfg, seed, out / "sim")
    return preflight(res, out / "pf", proteins=proteins)


# ---------------------------------------------------------------- driver with injection


def _producers(reg: ModuleType) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for s in reg.all_steps():
        for key in s.outputs:
            out.setdefault(key, []).append(s.name)
    return out


def proteins_key(reg: ModuleType) -> str:
    """The proteins input key of the ext ``miniprot`` step (``proteins.long.faa``)."""
    spec = reg.get("miniprot")
    keys = [k for k in (*spec.inputs, *spec.optional_inputs) if k != "genome"]
    assert len(keys) == 1, f"miniprot inputs {spec.inputs}: expected genome + one proteins key"
    return keys[0]


def plan(
    reg: ModuleType, targets: Sequence[str], supplied: set[str], inject: set[str]
) -> tuple[list[str], list[str]]:
    """(py steps in run order, ext steps whose outputs must be faked or injected).

    A step whose required inputs can neither be supplied nor produced is left out (e.g. the
    protein route when no ``proteins.long.faa`` exists), unless it is a target.
    """
    prod = _producers(reg)
    memo: dict[str, bool] = {}

    def ok(name: str, stack: frozenset[str] = frozenset()) -> bool:
        if name in inject:
            return True
        if name not in memo:
            if name in stack:
                return False
            spec = reg.get(name)
            memo[name] = all(
                k in supplied or any(ok(p, stack | {name}) for p in prod.get(k, ()))
                for k in spec.inputs
            )
        return memo[name]

    planned: set[str] = set()
    ext: set[str] = set()
    seen: set[str] = set()
    stack = list(targets)
    while stack:
        name = stack.pop()
        if name in seen:
            continue
        seen.add(name)
        if name not in targets and not ok(name):
            continue
        spec = reg.get(name)
        if spec.kind == "ext":
            ext.add(name)  # injected or faked; its inputs are still produced (e.g. leftovers)
        else:
            planned.add(name)
        for key in (*spec.inputs, *spec.optional_inputs):
            if key in supplied:
                continue
            stack.extend(prod.get(key, ()))
    return reg.topo_order(planned), sorted(ext)


def drive(
    reg: ModuleType,
    cfg: Any,
    work: Path,
    roots: Mapping[str, Path],
    inject: Mapping[str, Mapping[str, Path]] | None = None,
    targets: Sequence[str] = ("loci",),
) -> dict[str, Path]:
    """Run the py steps needed for ``targets`` through ``run_step``; returns every key's path.

    ``roots``: root keys → files (copied to ``work/<key>``). ``inject``: ext step → {output key:
    recorded file}. Other ext steps on the path whose inputs are satisfiable get empty output
    files (a tool that found nothing, e.g. the contaminant blastp); unsatisfiable ones are absent.
    """
    inject = dict(inject or {})
    work.mkdir(parents=True, exist_ok=True)
    avail: dict[str, Path] = {}
    for key, src in roots.items():
        dst = work / key
        shutil.copy(src, dst)
        avail[key] = dst
    order, ext = plan(reg, targets, set(avail), set(inject))
    for name in ext:
        spec = reg.get(name)
        given = inject.get(name)
        for key in spec.outputs:
            dst = work / key
            if given is not None:
                assert key in given, f"inject for ext step {name} lacks output key {key}"
                shutil.copy(given[key], dst)
            else:
                dst.write_text("")
            avail[key] = dst
    for name in order:
        spec = reg.get(name)
        missing = [k for k in spec.inputs if k not in avail]
        assert not missing, f"step {name}: required inputs {missing} not producible in the test"
        ins = {k: avail[k] for k in spec.inputs}
        ins |= {k: avail[k] for k in spec.optional_inputs if k in avail}
        outs = {k: work / k for k in spec.outputs}
        reg.run_step(name, cfg, ins, outs)
        avail.update(outs)
    return avail


def roots_of(run: SimRun, reg: ModuleType | None = None, with_proteins: bool = False) -> dict:
    roots = {"genome": run.genome, "peptides": run.peptides, "bakta_json": run.res.bakta_json}
    if with_proteins:
        assert reg is not None and run.proteins_long is not None
        roots[proteins_key(reg)] = run.proteins_long
    return roots


def key_like(keys: Sequence[str], stem: str) -> str:
    """The one key among ``keys`` naming the §6.2 file ``stem`` (``sixframe`` ↔ ``sixframe``,
    ``work/sixframe.npz``, ``sixframe.npz``)."""
    found = [k for k in keys if Path(k).name.split(".")[0] == stem]
    assert len(found) == 1, f"expected one key for {stem!r} among {list(keys)}"
    return found[0]


# ---------------------------------------------------------------- locating Module 2 outputs


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower().removesuffix(".tsv"))


def module2_files(root: Path, keys: Mapping[str, Path] | None = None) -> dict[str, Path]:
    """The seven §7.2.9 files under ``root`` (one directory) or among output ``keys``."""
    for d in sorted({p.parent for p in root.rglob("peptide_hits.tsv")}):
        if all((d / f).is_file() for f in MODULE2_FILES):
            return {f: d / f for f in MODULE2_FILES}
    if keys:
        by_norm = {_norm(k): p for k, p in keys.items() if p.is_file()}
        found = {f: by_norm[_norm(f)] for f in MODULE2_FILES if _norm(f) in by_norm}
        if len(found) == len(MODULE2_FILES):
            return found
    present = sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file())
    raise AssertionError(f"the seven Module 2 files {MODULE2_FILES} not found; files: {present}")


def read_tsv_df(path: Path) -> pd.DataFrame:
    """A mapping TSV as strings (``NA`` → missing)."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["NA"])
