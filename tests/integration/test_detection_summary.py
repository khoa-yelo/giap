"""M4-16 example 8 on the simulated fixture (docs/specs/M4-16-detection-summary.md).

Inputs come from merged producers: M4-07 layout and M4-08 expression (giapsim), C1-C13 written
to disk and read back with M4-06 ``load_m4_inputs``, M4-11 ``accepted_peptides`` /
``gene_evidence`` with the M4-12 two-peptide rule; ``p_detect`` is the simulator's
``p_detect_true`` until M4-13 merges, and the detectability / tier parts have the M4-13 / M4-12
row shapes (tests/helpers/detection_summary_fixture.py). The written table validates against the
``detection_summary`` schema; sections come in order tier, basic, detectability, qc.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import detection_summary_fixture as fx

TABLE_DIR = Path(__file__).resolve().parents[2] / "src" / "giap" / "schemas" / "tables"
B = 200


@pytest.fixture(scope="module")
def built(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, pd.DataFrame, Path]:
    from giap.context import summary
    from giap.io.summary_tsv import write_detection_summary
    from giap.rng import for_step

    work = tmp_path_factory.mktemp("work")
    sim = fx.sim_fixture(work, seed=7, n_genes=300)
    basic = summary.basic_rows(sim.genes, B, for_step("context_stats", 7))
    qc = summary.qc_rows(sim.inputs.peptides, sim.inputs.peptide_hits, sim.inputs.unmapped)
    table = summary.detection_summary([sim.tier, basic, sim.detectability, qc])
    out = work / "tables" / "detection_summary.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    write_detection_summary(out, table)
    return sim, table, out


def _descriptor() -> dict[str, Any]:
    found = [p for p in (TABLE_DIR / "detection_summary.json",
                         TABLE_DIR / "detection_summary.schema.json") if p.exists()]  # fmt: skip
    assert len(found) == 1, found
    d = json.loads(found[0].read_text(encoding="utf-8"))
    return {k: v for k, v in d.items() if k != "coordinates"}


def test_example_8_validates_against_schema(built: tuple[Any, pd.DataFrame, Path]) -> None:
    frictionless = pytest.importorskip("frictionless")
    _, _, out = built
    resource = frictionless.Resource(
        path=out.name, basepath=str(out.parent),
        schema=frictionless.Schema.from_descriptor(_descriptor()),
    )  # fmt: skip
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_example_8_sections_in_order(built: tuple[Any, pd.DataFrame, Path]) -> None:
    from giap.io.tsv import read_tsv

    _, _, out = built
    header, rows = read_tsv(out)
    assert tuple(header) == fx.SUMMARY_COLUMNS
    sections = [r[0] for r in rows]
    assert set(sections) == set(fx.SECTIONS)
    assert sections == sorted(sections, key=fx.SECTIONS.index)


def test_example_8_basic_rows_match_oracle(built: tuple[Any, pd.DataFrame, Path]) -> None:
    sim, table, _ = built
    basic = table[table["section"] == "basic"]
    groups = fx.group_members(sim.genes)
    keys = set(zip(basic["group_type"].astype(str), basic["group"].astype(str), strict=True))
    assert keys == set(groups)
    assert {"chr1", "p1", "ctg1"} <= {g for t, g in keys if t == "replicon"}
    assert any(t == "cog_category" and g == "none" for t, g in keys)
    for key, mask in groups.items():
        for name, v in fx.point_estimates(sim.genes, mask).items():
            got = float(fx.metric(basic, *key, name)["value"])
            assert (math.isnan(v) and math.isnan(got)) or math.isclose(got, v, abs_tol=1e-9)
    # partial genes (linear contig ends) are excluded
    n_elig = int((sim.genes["detection_status"] != "partial_excluded").sum())
    assert n_elig < len(sim.genes)
    total = sum(float(fx.metric(basic, "replicon", c, "n_genes")["value"])
                for t, c in keys if t == "replicon")  # fmt: skip
    assert total == n_elig


def test_example_8_qc_rows_from_contract_tables(built: tuple[Any, pd.DataFrame, Path]) -> None:
    _, table, _ = built
    qc = table[table["section"] == "qc"]
    val = {m: float(v) for m, v in zip(qc["metric"], qc["value"], strict=True)}
    # PR #87 lead decisions 5-6: id_mixups needs the variant module (NA + requires_variants);
    # every target is same-db, so genome_inconsistent_rate (search_db = other) is NA; the two
    # engine/protein mismatches are reported as engine_protein_mismatch_rate.
    mix = qc[qc["metric"] == "id_mixups"].iloc[0]
    assert math.isnan(val["id_mixups"])
    assert fx.REQUIRES_VARIANTS_FLAG in fx.flag_set(mix["flags"])
    assert math.isnan(val["genome_inconsistent_rate"])
    n_placed = int(built[0].inputs.peptide_hits["peptide"].astype(str).nunique())
    assert val["engine_protein_mismatch_rate"] == pytest.approx(2 / n_placed)
    assert val["contaminant_peptides"] == 1 and val["contaminant_proteins"] == 1
    assert val["unmapped_peptides"] == 1 and val["unmapped_proteins"] == 0


def test_example_8_detectability_rows_carried(built: tuple[Any, pd.DataFrame, Path]) -> None:
    sim, table, _ = built
    det = table[table["section"] == "detectability"].reset_index(drop=True)
    assert len(det) == len(sim.detectability)
    slope = det[det["metric"] == "calibration_slope"].iloc[0]
    assert (slope["ci_low"], slope["ci_high"]) == (0.81, 1.13)


def test_example_8_byte_identical_rerun(
    built: tuple[Any, pd.DataFrame, Path], tmp_path: Path
) -> None:
    from giap.context import summary
    from giap.io.summary_tsv import write_detection_summary
    from giap.rng import for_step

    sim, _, out = built
    basic = summary.basic_rows(sim.genes, B, for_step("context_stats", 7))
    qc = summary.qc_rows(sim.inputs.peptides, sim.inputs.peptide_hits, sim.inputs.unmapped)
    again = tmp_path / "again.tsv"
    write_detection_summary(again, summary.detection_summary([qc, sim.detectability, basic,
                                                              sim.tier]))  # fmt: skip
    assert again.read_bytes() == out.read_bytes()
