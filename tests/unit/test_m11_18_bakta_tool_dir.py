"""M11-18 (lead decision on #468): the Bakta JSON and GBFF reach ``annotation/tools/bakta/``
through the M2-20 ``tool_dirs`` input of ``context_layers``, not through a BAKTA ``publishDir``.

A ``publishDir`` on BAKTA cannot work: CONTEXT_LAYERS publishes the whole ``annotation``
directory afterwards and Nextflow's default overwrite deletes the destination first. So the
BAKTA process emits a ``tools/bakta`` directory (``genome.json`` and ``genome.gbff`` only) and
annotate.nf passes it as ``tool_dirs``; ``copy_native`` writes it to ``annotation/tools/bakta/``.

The M2-20 test ``test_layers_step.test_copy_native_byte_copies`` pins ``copy_native`` alone on
geNOMAD / TerminatorNet / Promotech recordings; this file pins the step path (``run_step`` with
the staged ``tool_dirs`` key, the dict keyed by the directory name) on a recorded Bakta 1.12.2
run of gippasim SIM100K, and that ``tools/`` gets nothing but the two files.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from helpers.layers_chain import REC, Chain, build_chain, fresh_steps

from gippa.config import GippaConfig

BAKTA_REC = REC / "bakta" / "1.12.2-light" / "sim100k"


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(tmp_path_factory.mktemp("bakta_tool_dir_chain"))


def _bakta_task(root: Path) -> Path:
    """A BAKTA task dir as the module writes it: bakta/ (every Bakta file) and tools/bakta/."""
    full = root / "bakta"
    full.mkdir(parents=True)
    for ext in ("json", "gbff", "gff3", "faa", "tsv"):
        (full / f"genome.{ext}").write_bytes(
            gzip.decompress((BAKTA_REC / f"sim100k.{ext}.gz").read_bytes())
        )
    tool = root / "tools" / "bakta"
    tool.mkdir(parents=True)
    for ext in ("json", "gbff"):
        (tool / f"genome.{ext}").write_bytes((full / f"genome.{ext}").read_bytes())
    return tool


def _run(out: Path, inputs: dict[str, object]) -> Path:
    outs = {"annotation": out / "annotation", "layers": out / "layers.parquet"}
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        reg.run_step("context_layers", GippaConfig(seed=1), inputs, outs)
    return outs["annotation"]


def _tree(root: Path, *, skip_tools: bool) -> dict[str, bytes]:
    return {
        str(p.relative_to(root)): p.read_bytes()
        for p in sorted(root.rglob("*"))
        if p.is_file() and not (skip_tools and "tools" in p.relative_to(root).parts)
    }


@pytest.mark.parametrize("as_list", [False, True], ids=["path", "list"])
def test_bakta_tool_dir_lands_in_annotation_tools_byte_identical(
    chain: Chain, tmp_path: Path, as_list: bool
) -> None:
    """tool_dirs=<task>/tools/bakta -> annotation/tools/bakta/genome.{json,gbff}, bytes equal."""
    tool = _bakta_task(tmp_path / "task")
    ann = _run(tmp_path / "out", {**chain.inputs, "tool_dirs": [tool] if as_list else tool})
    tools = ann / "tools"
    assert sorted(p.name for p in tools.iterdir()) == ["bakta"], "only bakta/ under tools/"
    got = _tree(tools / "bakta", skip_tools=False)
    assert sorted(got) == ["genome.gbff", "genome.json"], f"tools/bakta holds {sorted(got)}"
    for name, data in got.items():
        assert data == (tool / name).read_bytes(), f"{name} is not a byte copy (no gzip, no rename)"
        assert not (tools / "bakta" / name).is_symlink(), name
    # the uncompressed recorded Bakta outputs, untouched
    for ext in ("json", "gbff"):
        rec = gzip.decompress((BAKTA_REC / f"sim100k.{ext}.gz").read_bytes())
        assert got[f"genome.{ext}"] == rec, ext


def test_bakta_tool_dir_changes_nothing_else_in_annotation(chain: Chain, tmp_path: Path) -> None:
    """The tool dir adds tools/bakta only: all other annotation/* files equal the run without."""
    tool = _bakta_task(tmp_path / "task")
    with_tool = _run(tmp_path / "with", {**chain.inputs, "tool_dirs": tool})
    without = _run(tmp_path / "without", dict(chain.inputs))
    assert _tree(with_tool, skip_tools=True) == _tree(without, skip_tools=True)
    added = sorted(set(_tree(with_tool, skip_tools=False)) - set(_tree(without, skip_tools=False)))
    assert added == ["tools/bakta/genome.gbff", "tools/bakta/genome.json"], added
    leaked = [
        p
        for p in with_tool.rglob("*")
        if p.suffix in {".gff3", ".faa", ".tsv"} and "tools" in p.parts
    ]
    assert not leaked, f"other Bakta files under annotation/tools: {leaked}"
