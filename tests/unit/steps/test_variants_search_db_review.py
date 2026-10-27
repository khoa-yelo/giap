"""M7-04 review (PR #164) B1: C9 ``search_db: auto`` (the default) is resolved as placement
resolves it (``resolve_search_db`` against the locus tags), never read as ``same``. Database
peptides searched against another strain's proteins give no ALT evidence under
``genome_relation: same`` (§7.5.1 row 2)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from test_variants_step import _cfg, _chain

from gippa.io.fasta import read_genome_fasta
from gippa.steps.mapping import read_gene_model
from gippa.steps.variants import _peptide_hits, run_variants, support_of


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return _chain(tmp_path_factory.mktemp("variants_sdb"), strain=True)


def _rows(chain: dict[str, Path], search_db: str, protein_ids: str | None) -> int:
    peps = pd.read_parquet(chain["peptides"])
    peps["source"], peps["search_db"] = "database", search_db
    if protein_ids is not None:
        peps["protein_ids"] = protein_ids
    hits = _peptide_hits(Path(chain["mapping"]))
    hits["source"] = "database"
    genome = {r.id: r.sequence for r in read_genome_fasta(Path(chain["genome"]))}
    gm = read_gene_model(Path(chain["genes"]), Path(chain["replicons"]))
    tol = pd.read_parquet(chain["hits_tolerant"])
    _, table, _ = run_variants(hits, tol, peps, gm, genome, _cfg(genome_relation="same"))
    return len(table)


def test_auto_with_other_strain_ids_resolves_to_other(chain: dict[str, Path]) -> None:
    other = _rows(chain, "other", None)
    assert other < _rows(chain, "same", None), "fixture gives no same/other contrast"
    assert _rows(chain, "auto", "OTHERSTRAIN_00001") == other


def test_auto_with_locus_tag_ids_resolves_to_same(chain: dict[str, Path]) -> None:
    gm = read_gene_model(Path(chain["genes"]), Path(chain["replicons"]))
    tag = sorted(g.locus_tag for g in gm.genes)[0]
    assert _rows(chain, "auto", tag) == _rows(chain, "same", None)


def test_support_of_resolves_auto_and_keeps_explicit_values() -> None:
    hits = pd.DataFrame({"peptide": ["AAK", "CCK"], "source": "database", "placement_q": 0.0,
                         "uniqueness": "unique"})  # fmt: skip
    peps = pd.DataFrame({"peptide": ["AAK", "CCK"], "search_db": ["auto", "same"],
                         "protein_ids": ["X_1", "X_2"]})  # fmt: skip
    sup = support_of(hits, peps, {"LT_1"})
    assert (sup["AAK"].search_db, sup["CCK"].search_db) == ("other", "same")
