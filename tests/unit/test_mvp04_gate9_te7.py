"""Gate 9 (TE7 on the inputs) end to end on the 200-kb excerpt
(docs/specs/MVP-04-ecoli-evidence-run.md acceptance 9, §11.6 TE7, rule 3).

The excerpt TE7 test in tests/unit/test_mvp04_real_genome.py compares the tables directly; this
one runs the **summary's** ``gate9`` on two real ``scored_null`` gippasim calls (the evidence
scenario, the decoys that matter for gates 2-7), so the gate cannot pass on a comparison that
misses a table gippasim writes. Same seed: pass with no differing file; another seed: fail, naming
``peptides.tsv`` and the truth tables that changed. Three excerpt simulations, about 3 s in all,
so not marked slow.
"""

from __future__ import annotations

import importlib
import importlib.util
import sys
from dataclasses import replace
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers import mvp04_data as M

from gippasim.expression import IslandSpec

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "dev" / "mvp04_summary.py"
SEED = 20261023
ISLANDS = ((2, 3.0), (18, -3.0), (34, 3.0), (50, -3.0), (66, 3.0), (82, -3.0))
REQUIRED = {"peptides.tsv", "truth/peptides.tsv", "truth/genes.tsv", "truth/islands.tsv",
            "truth/alt_sites.tsv", "truth/novel_orfs.tsv"}  # fmt: skip


def summary() -> ModuleType:
    spec = importlib.util.spec_from_file_location("mvp04_summary_te7", SCRIPT)
    assert spec is not None and spec.loader is not None, SCRIPT
    m = importlib.util.module_from_spec(spec)
    sys.modules["mvp04_summary_te7"] = m
    spec.loader.exec_module(m)
    return m


def cfg() -> Any:
    base = importlib.import_module("gippasim.scenarios").k12_evidence()
    assert base.proteomics.decoys == "scored_null"
    expr = replace(base.expression,
                   islands=tuple(IslandSpec(M.CONTIG, f, 8, e) for f, e in ISLANDS))  # fmt: skip
    return replace(base, expression=expr, n_novel_per_class=20)


def simulate(ex: M.Excerpt, out: Path, seed: int) -> Path:
    importlib.import_module("gippasim.real_genome").simulate_on_genome(
        cfg(), seed, out, genome_fasta=ex.fasta, bakta_json=ex.bakta_json, tu_set=None,
        check_hashes=False,
    )  # fmt: skip
    return out


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    ex = M.make_excerpt(tmp_path_factory.mktemp("k12_excerpt"))
    root = tmp_path_factory.mktemp("te7")
    return {"a": simulate(ex, root / "a", SEED), "b": simulate(ex, root / "b", SEED),
            "other": simulate(ex, root / "other", SEED + 1)}  # fmt: skip


def test_gippasim_writes_every_table_gate9_needs(runs: dict[str, Path]) -> None:
    for rel in REQUIRED:
        assert (runs["a"] / rel).is_file(), rel


def test_same_seed_passes_gate9(runs: dict[str, Path]) -> None:
    res = summary().gate9(runs["a"], runs["b"])
    assert res["value"] == [] and res["status"] == "pass", res


def test_another_seed_fails_gate9(runs: dict[str, Path]) -> None:
    res = summary().gate9(runs["a"], runs["other"])
    assert res["status"] == "fail", res
    assert {"peptides.tsv", "truth/peptides.tsv", "truth/alt_sites.tsv"} <= set(res["value"]), res
    assert set(res["value"]) <= {"peptides.tsv"} | {
        str(p.relative_to(runs["a"])) for p in (runs["a"] / "truth").glob("*.tsv")
    }
