"""M4-12 on the real input contract: giapsim tables written to disk and read back by M4-06
``load_m4_inputs``, then M4-10 covariates, M4-11 evidence and ``detect``
(docs/specs/M4-12-detection-tiers.md; GIAP_SPEC.md §7.3.2, §9.3).

- decoys ``rev_<locus_tag>`` (M4-09) -> ``decoy_picked``, every decoy paired;
- ``protein_ids = foreign`` (UniProt-like ``UP_<n>``, pilot arm A1 without ``proteins.faa``):
  the peptide-vote map sends every mapped ``UP_<100001 + i>`` to gene i, method ``peptide_vote``
  (M4-12 Q4); M4-14 relies on the map's columns ``protein_id``, ``gene``, ``method``;
- proteins only (C12 rank-1 identical hits) -> ``protein_list``, detected = genes with a hit.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_fixtures as fx

SCORED = frozenset({"q_value", "pep", "score", "is_decoy", "intensity", "protein_ids"})


def _load(tmp_path: Path, **ev_kw: Any) -> tuple[Any, Any]:
    import giapsim.evidence as ev
    import giapsim.expression as ex
    import giapsim.layout as la
    from giap.io.work_tables import FILE_NAMES, load_m4_inputs
    from giapsim.contract_tables import write_contract_tables

    rng = np.random.default_rng(20261015)
    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=200), rng)
    expr = ex.simulate_expression(lay, ex.ExpressionConfig(), rng)
    e = ev.simulate_evidence(lay, expr, ev.EvidenceConfig(**ev_kw), rng)
    write_contract_tables(lay, tmp_path)
    ev.write_evidence_tables(e, tmp_path)
    paths = {}
    for key, name in FILE_NAMES.items():
        p = tmp_path / name
        paths[key] = p if p.exists() else tmp_path / "mapping" / name
    inputs, issues = load_m4_inputs(paths)
    assert not [i for i in issues if i.severity == "error"], issues
    return inputs, e


def _detect(inputs: Any, config: Any = None) -> Any:
    from giap.genetable.covariates import covariate_table
    from giap.genetable.detection import detect
    from giap.genetable.evidence import accepted_peptides, gene_evidence

    c = config if config is not None else fx.cfg()
    acc = accepted_peptides(inputs.peptide_hits, inputs.peptides, c)
    ev = gene_evidence(inputs.genes, inputs.peptide_hits, acc, inputs.replicons)
    cov = covariate_table(inputs.genes, c)
    return detect(inputs.genes, ev, acc, inputs.peptides, inputs.peptide_hits,
                  inputs.protein_hits, cov, c)  # fmt: skip


def test_loaded_tables_with_decoys_give_decoy_picked(tmp_path: Path) -> None:
    inputs, _ = _load(tmp_path, pi_wrong=0.05, emit=SCORED)
    det = _detect(inputs)
    assert det.tier == "decoy_picked"
    t = det.table
    assert t["locus_tag"].astype(str).tolist() == inputs.genes["locus_tag"].astype(str).tolist()
    assert (t["detection_tier"] == "decoy_picked").all()
    assert fx.tier_metric(det, "n_decoys_unpaired") == 0
    assert fx.tier_metric(det, "n_decoys_paired") >= 1
    assert t["detected"].astype(bool).sum() > 0
    s = fx.summary(det)
    assert set(s["section"]) == {"tier"}


def test_foreign_ids_mapped_by_peptide_vote(tmp_path: Path) -> None:
    from giap.genetable.detection import protein_gene_map

    inputs, e = _load(tmp_path, pi_wrong=0.01, emit=SCORED, protein_ids="foreign")
    pmap = protein_gene_map(inputs.genes, inputs.peptides, inputs.peptide_hits,
                            inputs.protein_hits)  # fmt: skip
    assert {"protein_id", "gene", "method"} <= set(pmap.columns)
    tags = inputs.genes["locus_tag"].astype(str).tolist()
    seen = 0
    for i, tag in enumerate(tags):
        gene, method = fx.mapped(pmap, f"UP_{100001 + i}")
        if gene is None:
            continue
        seen += 1
        assert (gene, method) == (tag, "peptide_vote"), i
    n_with_peptides = int((e.gene_truth["n_true_peptides"] > 0).sum())
    assert seen >= n_with_peptides


def test_proteins_only_gives_protein_list(tmp_path: Path) -> None:
    inputs, e = _load(tmp_path, proteins_only=True)
    det = _detect(inputs)
    assert det.tier == "protein_list"
    expected = set(e.gene_truth.loc[e.gene_truth["detected_true"], "locus_tag"].astype(str))
    got = set(det.table.loc[det.table["detected"].astype(bool), "locus_tag"].astype(str))
    assert got == expected
    assert det.table["gene_q_value"].isna().all()


@pytest.mark.parametrize("emit", [frozenset({"protein_ids", "intensity"})])
def test_sequences_only_gives_two_peptide(tmp_path: Path, emit: frozenset[str]) -> None:
    inputs, _ = _load(tmp_path, pi_wrong=0.01, emit=emit)
    det = _detect(inputs)
    assert det.tier == "two_peptide"
    st = set(det.table["detection_status"].astype(str))
    assert st <= {"detected", "not_detected", "single_peptide_unestimated"}
