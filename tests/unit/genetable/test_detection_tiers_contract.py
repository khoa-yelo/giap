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


# ------------------------------------------------------------- lead decisions 5, 6, 8 (PR #92)


def test_pep_only_loaded_tables_detect_genes(tmp_path: Path) -> None:
    """Decision 5: giapsim without q (PEP only, no decoys) -> ``mayu`` with genes detected.

    giapsim's q is the running mean of target PEP, so the PEP-derived q reproduces it and the
    detected set matches the q + PEP run of the same replicate up to the acceptance difference
    (PEP-only rows have no q to filter on)."""
    pep_dir, q_dir = tmp_path / "pep", tmp_path / "q"
    pep_dir.mkdir()
    q_dir.mkdir()
    inputs, e = _load(pep_dir, pi_wrong=0.01, emit=frozenset({"pep", "intensity", "protein_ids"}))
    det = _detect(inputs)
    assert det.tier == "mayu"
    got = set(det.table.loc[det.table["detected"].astype(bool), "locus_tag"].astype(str))
    assert got, "PEP-only input must not detect 0 genes silently"
    q_inputs, _ = _load(q_dir, pi_wrong=0.01,
                        emit=frozenset({"q_value", "pep", "intensity", "protein_ids"}))  # fmt: skip
    ref = _detect(q_inputs)
    want = set(ref.table.loc[ref.table["detected"].astype(bool), "locus_tag"].astype(str))
    assert len(got) >= 0.8 * len(want), (len(got), len(want))
    absent = set(e.gene_truth.loc[~e.gene_truth["present"].astype(bool), "locus_tag"].astype(str))
    assert len(got & absent) <= max(2, 0.05 * len(got))


@pytest.mark.parametrize("flag_contaminant", [False, True])
def test_contaminant_decoys_injected_into_loaded_tables(
    tmp_path: Path, flag_contaminant: bool
) -> None:
    """Decision 6: decoys of contaminants (``rev_CON__...``) added to giapsim's C9 in a number
    that would push the paired share below 0.8 leave tier, counts and gene q unchanged."""
    import dataclasses

    import pandas as pd

    inputs, _ = _load(tmp_path, pi_wrong=0.05, emit=SCORED)
    base = _detect(inputs)
    assert base.tier == "decoy_picked"
    peps = inputs.peptides
    d = peps[peps["is_decoy"].astype(bool) & (peps["q_value"] <= fx.cfg().input.q_value_max)]
    k = len(d) // 2 + 1
    extra = d.iloc[np.arange(k) % len(d)].copy()
    extra["peptide"] = [fx.peptide_name(900_000 + i) for i in range(k)]
    extra["protein_ids"] = [[f"rev_CON__P{10000 + i}"] for i in range(k)]
    extra["is_contaminant"] = flag_contaminant
    inputs2 = dataclasses.replace(inputs, peptides=pd.concat([peps, extra], ignore_index=True))
    det = _detect(inputs2)
    assert det.tier == "decoy_picked"
    assert "decoy_pairing_low" not in det.warnings
    for m in ("n_decoys_paired", "n_decoys_unpaired"):
        assert fx.tier_metric(det, m) == fx.tier_metric(base, m), m
    pd.testing.assert_frame_equal(det.table, base.table)


def test_detection_map_matches_detect_inputs_on_foreign_ids(tmp_path: Path) -> None:
    """Decision 8: ``Detection.protein_gene_map`` is the map ``detect`` built (accepted peptides,
    configured contaminant prefixes), the one M4-14 reuses."""
    from giap.genetable.detection import protein_gene_map
    from giap.genetable.evidence import accepted_peptides

    inputs, _ = _load(tmp_path, pi_wrong=0.05, emit=SCORED, protein_ids="foreign")
    c = fx.cfg()
    det = _detect(inputs, c)
    acc = accepted_peptides(inputs.peptide_hits, inputs.peptides, c)
    ref = protein_gene_map(inputs.genes, inputs.peptides, inputs.peptide_hits, inputs.protein_hits,
                           accepted=acc,
                           contaminant_prefixes=c.input.contaminant_prefixes)  # fmt: skip

    def rows(m: Any) -> set[tuple[str, str, str]]:
        m = m[m["gene"].notna()]
        return set(zip(m["protein_id"].astype(str), m["gene"].astype(str),
                       m["method"].astype(str), strict=True))  # fmt: skip

    assert rows(det.protein_gene_map) == rows(ref)
    assert any(meth == "peptide_vote" for _, _, meth in rows(det.protein_gene_map))
