"""Driver and checks shared by the M4-20 tests (docs/specs/M4-20-integration-golden.md).

The M4 chain on a giapsim fixture, as the spec's "Interface" fixes it:

* ``giapsim.fixtures.build_m4_fixture(work_dir, seed, size, evidence)`` (tooling the implementer
  writes) puts the M4-06 contract tables C1-C14 into ``work_dir`` -- flat (``default_paths``) and
  the three M3 TSVs also under ``mapping/`` (``giap.steps.gene_table.work_inputs`` reads that
  layout, ``run_context_stats`` and ``run_sensitivity`` the flat one, so a fixture that serves the
  whole chain has both);
* ``run_m4`` then calls ``run_gene_table`` -> ``run_detectability`` -> ``run_context_stats`` ->
  ``run_sensitivity`` in the spec's order and returns the four ``tables/`` outputs.

Who writes what in ``tables/``: ``context_stats`` writes ``genes.tsv`` (M5-09 lead decision 3: it
is the only step that holds both the gene table and M5-09's ``island_id`` / ``replichore`` /
``replichore_position``, so writing it earlier would mean writing the file twice),
``regions_stats.tsv`` and ``detection_summary.tsv``; ``sensitivity`` writes ``sensitivity.tsv``.
Known limitation, recorded rather than worked around: a run without C7/C8 skips ``context_stats``
and so produces no ``genes.tsv``; §6.2 gives the file to ``finalize``, which resolves it.

``$GIAP_DB_DIR`` always points at ``test_db`` here (never the real one, M3-17 review decision 12):
the gene table reads ``cog2020_def`` from the manifest, so without a pinned fixture COG table the
golden run would depend on whether the shared database is installed.

``layout=`` is in the spec's Interface line since M4-20 lead decision 4 (example 7 needs
``linear_contigs=(30,)`` and a plasmid, which ``size`` alone cannot express).
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gene_table_fixture as gt
from mapping_steps import write_test_db

REPO = Path(__file__).resolve().parents[2]
SCHEMA_DIR = REPO / "src" / "giap" / "schemas" / "tables"
SIZES = ("sim50k", "sim1000", "sim4500")

#: C1-C9 and C11-C14 of the M4-06 contract, as ``build_m4_fixture`` must write them (C10
#: ``abundance.norm.parquet`` only with ``EvidenceConfig.abundance_table``).
CONTRACT_FILES: tuple[str, ...] = (
    "replicons.parquet", "genes.parquet", "xrefs.parquet", "annotation_meta.json", "tus.parquet",
    "layers.parquet", "regions.parquet", "region_genes.parquet", "peptides.norm.parquet",
    "peptide_hits.tsv", "protein_hits.tsv", "unmapped.tsv", "id_map.tsv",
    "mapping/peptide_hits.tsv", "mapping/protein_hits.tsv", "mapping/unmapped.tsv",
)  # fmt: skip

#: §9.1 ``tables/`` outputs of the M4 chain (spec example 1), file name -> schema stem
OUTPUT_TABLES: Mapping[str, str] = {
    "genes.tsv": "genes",
    "detection_summary.tsv": "detection_summary",
    "regions_stats.tsv": "regions_stats",
    "sensitivity.tsv": "sensitivity",
}


def build(
    work: Path,
    *,
    seed: int,
    size: str = "sim1000",
    evidence: Any | None = None,
    layout: Mapping[str, Any] | None = None,
) -> Path:
    """``build_m4_fixture`` into ``work`` (created); returns ``work``."""
    from giapsim.fixtures import build_m4_fixture

    kw: dict[str, Any] = {}
    if evidence is not None:
        kw["evidence"] = evidence
    if layout is not None:
        kw["layout"] = dict(layout)
    work.mkdir(parents=True, exist_ok=True)
    build_m4_fixture(work, seed=seed, size=size, **kw)
    return work


def evidence_config(**over: Any) -> Any:
    """``giapsim.evidence.EvidenceConfig`` with overrides (M4-09)."""
    from giapsim.evidence import EvidenceConfig

    return EvidenceConfig(**over)


def config(seed: int, **over: Any) -> Any:
    """``GiapConfig`` from the shipped defaults, ``seed`` and nested overrides, e.g.
    ``config(7, stats={"permutations": 199})``."""
    from giap.config.model import GiapConfig

    raw: dict[str, Any] = {"seed": seed}
    for k, v in over.items():
        raw[k] = {**raw[k], **v} if isinstance(v, Mapping) and isinstance(raw.get(k), dict) else v
    return GiapConfig.model_validate(raw)


def test_db(root: Path) -> Path:
    """A test ``$GIAP_DB_DIR`` whose manifest lists the fixture ``cog-20.def.tab``."""
    db = write_test_db(root)
    gt.write_cog_def(db / "cog" / "cog-20.def.tab")
    man = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    man["cog2020_def"] = {"version": "2020", "path": "cog/cog-20.def.tab"}
    (db / "db_manifest.json").write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")
    return db


@contextmanager
def db_env(db: Path) -> Iterator[None]:
    """``$GIAP_DB_DIR`` = ``db`` for the block (restored afterwards)."""
    old = os.environ.get("GIAP_DB_DIR")
    os.environ["GIAP_DB_DIR"] = str(db)
    try:
        yield
    finally:
        if old is None:
            os.environ.pop("GIAP_DB_DIR", None)
        else:
            os.environ["GIAP_DB_DIR"] = old


def run_m4(
    work: Path, cfg: Any, db: Path, *, replicon_depth: Path | None = None
) -> dict[str, Path]:
    """The four M4 steps in spec order; returns the ``tables/`` paths by file name."""
    from giap.steps.context_stats import run_context_stats
    from giap.steps.detectability import run_detectability
    from giap.steps.gene_table import run_gene_table
    from giap.steps.sensitivity import run_sensitivity

    with db_env(db):
        run_gene_table(work, cfg)
        run_detectability(work, cfg)
        run_context_stats(work, cfg, replicon_depth)
        run_sensitivity(work, cfg)
    return {name: work / "tables" / name for name in OUTPUT_TABLES}


def genes(work: Path) -> pd.DataFrame:
    """The final gene table (``genetable/genes.parquet``)."""
    return pd.read_parquet(work / "genetable" / "genes.parquet")


def read(path: Path) -> pd.DataFrame:
    """A written TSV as strings (``NA`` kept verbatim, as the writers emit it)."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def descriptor(stem: str) -> dict[str, Any]:
    """The Frictionless resource of ``src/giap/schemas/tables/<stem>`` (``coordinates`` is a
    GIAP extension Frictionless does not know)."""
    found = [p for p in (SCHEMA_DIR / f"{stem}.json", SCHEMA_DIR / f"{stem}.schema.json")
             if p.is_file()]  # fmt: skip
    assert len(found) == 1, f"expected one {stem} schema in {SCHEMA_DIR}, got {found}"
    d: dict[str, Any] = json.loads(found[0].read_text(encoding="utf-8"))
    return {k: v for k, v in d.items() if k != "coordinates"}


def validate(path: Path, stem: str) -> None:
    """§9.5: the written table validates with ``frictionless`` against its resource."""
    frictionless = pytest.importorskip("frictionless")
    schema = frictionless.Schema.from_descriptor(descriptor(stem))
    resource = frictionless.Resource(path=path.name, basepath=str(path.parent), schema=schema)
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def digests(paths: Mapping[str, Path]) -> dict[str, str]:
    """SHA-256 per written file (TE7: two runs with the same seed give byte-identical tables).

    A path the chain did not write is absent from the result, so the TE7 comparison reports a
    content difference and ``test_example_1_writes_the_four_tables`` the missing table.
    """
    return {
        k: hashlib.sha256(p.read_bytes()).hexdigest()
        for k, p in sorted(paths.items())
        if p.is_file()
    }


def first_diff(want: Path, got: Path) -> list[tuple[int, str, str]]:
    """The first five differing lines as (line number, want, got)."""
    a, b = (
        want.read_text(encoding="utf-8").splitlines(),
        got.read_text(encoding="utf-8").splitlines(),
    )
    out = [(i + 1, x, y) for i, (x, y) in enumerate(zip(a, b, strict=False)) if x != y][:5]
    if not out and len(a) != len(b):
        out = [(min(len(a), len(b)) + 1, f"{len(a)} lines", f"{len(b)} lines")]
    return out
