"""PR #181 review round (M8-06a, docs/specs/M8-06-function-step-and-run.md): two findings.

* **B1.** The step hands C1 (``replicons``) to GSEA, so a term whose TUs sit at the first and the
  last ranked positions of a circular replicon is a ``genomic_cluster`` (decision 9a: they are
  adjacent across the origin), and the same term on a linear replicon is not (they are
  ``n - 1`` TUs apart, beyond the AR(1) kernel's reach). The degenerate case "circular vs
  linear" of the spec and the reason ``replicons`` is a required input.
* **N1.** §6.7 events are JSON lines: with no eligible gene but a resource installed, the
  ``function.collections`` event has an undefined coverage and must write ``null``, not a bare
  ``NaN`` (not RFC 8259 JSON).
"""

from __future__ import annotations

import json
import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.function_db import FunctionDb, write_function_db
from helpers.function_step_data import (
    CONTIG,
    full_frame,
    toy4,
    write_gene_table,
    write_replicons,
    write_tu_folds,
)

from giap.config import GiapConfig

OUTPUTS = {
    "category_summary": "category_summary.tsv",
    "enrichment_detection": "enrichment_detection.tsv",
    "enrichment_gsea": "enrichment_gsea.tsv",
}
COLL_EVENT = "function.collections"
N_RING = 40  # ranked TUs on the ring; the term's two TUs are TU_1 and TU_40


def _cfg(**function: Any) -> GiapConfig:
    raw: dict[str, Any] = {"seed": 7, "stats": {"bootstrap": 100},
                           "function": {"gsea": {"n_perm_simple": 200}, **function}}  # fmt: skip
    return GiapConfig.model_validate(raw)


def _run_step(inputs: Mapping[str, Path], out: Path, cfg: GiapConfig) -> dict[str, Path]:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    outputs = {k: out / "tables" / v for k, v in OUTPUTS.items()}
    reg.run_step("function", cfg, dict(inputs), outputs)
    return outputs


def _empty_xrefs() -> pd.DataFrame:
    return pd.DataFrame({"locus_tag": [], "db": [], "id": []}, dtype=object)


@pytest.fixture(scope="module")
def toy_db(tmp_path_factory: pytest.TempPathFactory) -> FunctionDb:
    return write_function_db(tmp_path_factory.mktemp("review181_db"), _empty_xrefs())


def _tokens(cell: str) -> set[str]:
    return set() if cell == "NA" else set(cell.split(";"))


# ----------------------------------------------------------------------------- B1


def ring_genes(n: int = N_RING) -> pd.DataFrame:
    """``n`` one-gene TUs TU_1..TU_n on ``chr1``, all detected with distinct abundances; only
    the first and the last (TU_1, TU_n: positions 0 and n - 1 of the ranked order) carry COG
    J, so the J set is one origin-spanning pair."""
    cog = ["J" if i in (0, n - 1) else None for i in range(n)]
    cols: dict[str, Any] = {
        "locus_tag": [f"g{i + 1}" for i in range(n)],
        "feature_type": ["cds"] * n,
        "contig": [CONTIG] * n,
        "start": [1000 * i for i in range(n)],
        "end": [1000 * i + 900 for i in range(n)],
        "strand": ["+"] * n,
        "partial": [False] * n,
        "crosses_origin": [False] * n,
        "hypothetical": [False] * n,
        "tu_id": [f"TU_{i + 1}" for i in range(n)],
        "tu_position": [1] * n,
        "tu_size": [1] * n,
        "detected": [True] * n,
        "detection_status": ["detected"] * n,
        "p_detect": [0.9] * n,
        "abundance_log2": [20.0 + 0.1 * ((7 * i) % n) for i in range(n)],
        "abundance_imputed": [False] * n,
        "abundance_type": ["ibaq"] * n,
        "abundance_scale": ["absolute"] * n,
        "cog_category": cog,
        "cog_id": ["COG0001" if c else None for c in cog],
        "annotation_source": ["bakta_full"] * n,
    }
    return full_frame(cols)


@pytest.fixture(scope="module")
def ring(toy_db: FunctionDb, tmp_path_factory: pytest.TempPathFactory) -> dict[str, pd.DataFrame]:
    """The GSEA table of the ring genes on a circular and on a linear C1 (same everything
    else). ``gsea_size`` [2, 500] admits the two-TU J set."""
    root = tmp_path_factory.mktemp("review181_ring")
    base = {
        "gene_table": write_gene_table(root / "genetable" / "genes.parquet", ring_genes()),
        "tu_folds": write_tu_folds(root / "genetable" / "tu_folds.parquet", {}),
    } | toy_db.inputs(["cog2020_fun"])
    topo = {
        "circular": write_replicons(root / "replicons_circ.parquet", topology="circular"),
        "linear": write_replicons(root / "replicons_lin.parquet", topology="linear"),
    }
    cfg = _cfg(gsea_size=[2, 500])
    out: dict[str, pd.DataFrame] = {}
    for name, reps in topo.items():
        paths = _run_step(base | {"replicons": reps}, root / name, cfg)
        out[name] = pd.read_csv(paths["enrichment_gsea"], sep="\t", dtype=str,
                                keep_default_na=False)  # fmt: skip
    return out


def _j_flags(df: pd.DataFrame) -> set[str]:
    rows = df[(df["collection"] == "cog") & (df["term"] == "J")]
    assert len(rows) == 1, df[["collection", "term"]].to_dict("records")
    assert rows.iloc[0]["size_tus"] == "2"
    return _tokens(rows.iloc[0]["flags"])


def test_origin_spanning_term_is_a_genomic_cluster_on_a_circular_replicon(
    ring: dict[str, pd.DataFrame],
) -> None:
    """TU_1 and TU_40 are adjacent across the origin (d = 1, delta ~ 0.5 >= tau 0.15)."""
    assert "genomic_cluster" in _j_flags(ring["circular"])


def test_origin_spanning_term_is_not_a_genomic_cluster_on_a_linear_replicon(
    ring: dict[str, pd.DataFrame],
) -> None:
    """On a linear C1 the same TUs are 39 apart (> d_max 10): no clustering."""
    assert "genomic_cluster" not in _j_flags(ring["linear"])


def test_topology_changes_only_the_genomic_cluster_flag_of_the_ring_term(
    ring: dict[str, pd.DataFrame],
) -> None:
    a, b = ring["circular"], ring["linear"]
    rest = [c for c in a.columns if c != "flags"]
    assert a[rest].equals(b[rest])
    assert _j_flags(a) - {"genomic_cluster"} == _j_flags(b)


# ----------------------------------------------------------------------------- N1


def _no_constant(name: str) -> Any:
    raise ValueError(f"non-JSON constant {name} in events.jsonl")


@pytest.mark.parametrize("case", ["all_partial", "empty"])
def test_no_eligible_gene_events_are_strict_json_with_null_coverage(
    case: str, toy_db: FunctionDb, tmp_path: Path
) -> None:
    """§6.7: every events.jsonl line is RFC 8259 JSON; the undefined coverage (0 eligible genes)
    of the installed COG collection is ``null``."""
    from giap.io.events import EventsHandler

    genes = toy4(partial=[True] * 4) if case == "all_partial" else toy4().iloc[0:0]
    ins = {
        "gene_table": write_gene_table(tmp_path / "genetable" / "genes.parquet", genes),
        "replicons": write_replicons(tmp_path / "replicons.parquet"),
        "tu_folds": write_tu_folds(tmp_path / "genetable" / "tu_folds.parquet", {}),
    } | toy_db.inputs(["cog2020_fun"])
    logs = tmp_path / "logs"
    handler = EventsHandler(logs)
    steps_log = logging.getLogger("giap.steps")
    old = steps_log.level
    steps_log.addHandler(handler)
    steps_log.setLevel(logging.INFO)
    try:
        _run_step(ins, tmp_path / "out", _cfg())
    finally:
        steps_log.removeHandler(handler)
        steps_log.setLevel(old)
    lines = (logs / "events.jsonl").read_text(encoding="utf-8").splitlines()
    records = [json.loads(line, parse_constant=_no_constant) for line in lines]
    coll = [r for r in records if r["event"] == COLL_EVENT]
    assert len(coll) == 1, records
    assert coll[0]["data"]["cog"]["n_terms"] == 0
    assert coll[0]["data"]["cog"]["coverage"] is None
