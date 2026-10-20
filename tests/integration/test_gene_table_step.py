"""M4-15 examples 1 and 10 on the giapsim default fixture (docs/specs/M4-15-gene-table-assembly.md;
GIAP_SPEC.md §6.6/R7, §7.3.3, §9.5).

The M4-07/08/09 simulator writes the contract tables of a 1,000-gene layout to ``tmp_path`` (C1-C8
flat, C9/C10 flat, C11-C13 under ``mapping/``, as the producers do); the registered steps
``gene_table`` and ``detectability`` run through ``run_step``; ``genes.tsv`` is written with
``write_genes_tsv`` and validated with frictionless against ``schemas/tables/genes.json``.

Lead decision 7: giapsim writes ``annotation_source = "bakta_full"``, so its output validates as
written (no rewrite here). The detectability output parquet is the declared output of
``detectability`` that holds a ``p_detect`` column.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import gene_table_fixture as gt

SEED = 20261020
KEYS = {
    "genes": "genes.parquet", "replicons": "replicons.parquet", "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json", "tus": "tus.parquet", "layers": "layers.parquet",
    "regions": "regions.parquet", "region_genes": "region_genes.parquet",
    "peptides": "peptides.norm.parquet", "abundance": "abundance.norm.parquet",
}  # fmt: skip
TABLE_DIR = Path(__file__).resolve().parents[2] / "src" / "giap" / "schemas" / "tables"


def _simulate(work: Path) -> dict[str, Path]:
    import giapsim.evidence as ev
    import giapsim.expression as ex
    import giapsim.layout as la
    from giapsim.contract_tables import write_contract_tables

    rng = np.random.default_rng(SEED)
    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=1000), rng)
    expr = ex.simulate_expression(lay, ex.ExpressionConfig(), rng)
    e = ev.simulate_evidence(lay, expr, ev.EvidenceConfig(), rng)
    write_contract_tables(lay, work)
    ev.write_evidence_tables(e, work)
    keys = {k: work / f for k, f in KEYS.items() if (work / f).exists()}
    return keys | {"mapping": work / "mapping"}


def _pipeline(root: Path) -> tuple[Path, pd.DataFrame, pd.DataFrame]:
    """gene_table -> detectability -> genes.tsv; returns (tsv, genes_pre, genes)."""
    from giap.config import GiapConfig
    from giap.io.genes_tsv import write_genes_tsv
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    cfg = GiapConfig(seed=SEED)
    keys = _simulate(root / "work")
    spec = reg.get("gene_table")
    ins = {k: p for k, p in keys.items() if k in (*spec.inputs, *spec.optional_inputs)}
    outs = {k: root / "out" / k for k in spec.outputs}
    outs["genes_pre"] = root / "out" / "genetable" / "genes_pre.parquet"
    reg.run_step("gene_table", cfg, ins, outs)

    dspec = reg.get("detectability")
    known = keys | {"genes_pre": outs["genes_pre"]}
    dins = {k: p for k, p in known.items() if k in (*dspec.inputs, *dspec.optional_inputs)}
    douts = {k: root / "out" / "det" / k for k in dspec.outputs}
    reg.run_step("detectability", cfg, dins, douts)
    genes = [pd.read_parquet(p) for p in douts.values()
             if p.is_file() and "p_detect" in pd.read_parquet(p).columns]  # fmt: skip
    assert len(genes) == 1, douts
    tsv = root / "genes.tsv"
    write_genes_tsv(tsv, genes[0])
    return tsv, pd.read_parquet(outs["genes_pre"]), genes[0]


@pytest.fixture(scope="module")
def run_a(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, pd.DataFrame, pd.DataFrame]:
    return _pipeline(tmp_path_factory.mktemp("gene_table_a"))


def _schema() -> dict[str, Any]:
    found = [p for p in (TABLE_DIR / "genes.json", TABLE_DIR / "genes.schema.json") if p.exists()]
    assert len(found) == 1, found
    return json.loads(found[0].read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def test_example_1_default_fixture(run_a: tuple[Path, pd.DataFrame, pd.DataFrame]) -> None:
    frictionless = pytest.importorskip("frictionless")
    tsv, _, genes = run_a
    lines = tsv.read_text(encoding="utf-8").splitlines()
    assert tuple(lines[0].split("\t")) == gt.GENES_COLUMNS
    assert len(lines) - 1 == len(genes) == 1000
    d = {k: v for k, v in _schema().items() if k != "coordinates"}
    res = frictionless.Resource(
        path=tsv.name, basepath=str(tsv.parent), schema=frictionless.Schema.from_descriptor(d)
    )
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_p_detect_is_the_cross_fitted_m4_13_prediction(
    run_a: tuple[Path, pd.DataFrame, pd.DataFrame],
) -> None:
    """Methods note: ``p_detect`` is attached only from the cross-fitted M4-13 predictions (the
    ``detectability`` step's own stream ``for_step("detectability", seed)``)."""
    from giap.config import GiapConfig
    from giap.genetable.detectability import fit_detectability
    from giap.rng import for_step

    _, pre, genes = run_a
    pre = pre.set_index("locus_tag")
    detected = pre["detected"].astype("boolean").where(~pre["partial"].astype(bool), pd.NA)
    cov = pre[["length_aa", "n_observable_peptides", "gravy", "pi", "mw_da", "tm_segments"]]
    cfg = GiapConfig(seed=SEED)
    want = fit_detectability(detected, cov.reset_index(), pre["tu_id"], cfg.detectability,
                             for_step("detectability", SEED))  # fmt: skip
    got = genes.set_index("locus_tag")["p_detect"].astype(float)
    exp = want.p_detect.reindex(got.index).astype(float)
    assert np.allclose(got.to_numpy(), exp.to_numpy(), equal_nan=True, rtol=1e-12)
    assert genes["p_detect"].notna().sum() == int((~pre["partial"].astype(bool)).sum())


def test_example_10_same_seed_byte_identical(
    run_a: tuple[Path, pd.DataFrame, pd.DataFrame], tmp_path: Path
) -> None:
    tsv_b, _, _ = _pipeline(tmp_path)
    assert run_a[0].read_bytes() == tsv_b.read_bytes()
    text = tsv_b.read_text(encoding="utf-8")
    assert str(tmp_path) not in text  # no absolute paths in tables (§6.6)
