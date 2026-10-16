"""M4-16 QC rows and the ``detection_summary`` assembler (spec examples 6-7, sorting, validation).

Spec: GIAP_SPEC.md §9.3 (QC rows ``id_mixups``, ``genome_inconsistent_rate``, contaminant counts),
§7.2.7 (genome relation table), docs/specs/M4-16-detection-summary.md (``qc_rows``,
``detection_summary``; lead question 1 default for ``genome_inconsistent_rate``).

Resolved ambiguities (pinned here):
- QC metrics: ``id_mixups``, ``engine_protein_mismatch_rate``, ``genome_inconsistent_rate``,
  ``contaminant_peptides``, ``contaminant_proteins``, ``unmapped_peptides``,
  ``unmapped_proteins`` (one row each, section ``qc``, ``group_type = run``, ``group = all``);
  C13 counts are rows with that ``label``, by ``kind``.
- Aligned after the statistician's review (PR #87 lead decisions 5-6; detailed tests in
  test_detection_summary_review.py): ``genome_inconsistent_rate`` is over target, non-contaminant
  ``search_db = other`` peptides (NA when C9 is absent or has none; ``same`` rows and decoys are
  ignored); a peptide with an exact C11 row is placed exactly even if it also has tolerant rows.
  ``id_mixups`` is NA with flag ``requires_variants``; the distinct ``engine_protein_mismatch``
  peptides are reported as ``engine_protein_mismatch_rate``.
- ``detection_summary`` accepts parts carrying only the seven core columns (M4-13 diagnostics and
  M4-12 summary rows have no ``n`` / ``flags``): missing ``n`` becomes NA, missing ``flags`` empty.
  An unknown section or a part without a core column raises ``ValueError``.
- Within a section rows sort by ``group_type``, ``group``, ``metric`` as strings.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_summary_fixture as fx
from detection_summary_fixture import gene, qc_hit, qc_peptide

QC_METRICS = ("id_mixups", "engine_protein_mismatch_rate", "genome_inconsistent_rate",
              "contaminant_peptides", "contaminant_proteins", "unmapped_peptides",
              "unmapped_proteins")  # fmt: skip


def _summary() -> Any:
    from giap.context import summary

    return summary


def _qc(peps: pd.DataFrame | None, hits: pd.DataFrame, unmapped: pd.DataFrame | None) -> Any:
    return _summary().qc_rows(peps, hits, unmapped)


def _val(df: pd.DataFrame, name: str) -> float:
    return float(fx.metric(df, "run", "all", name)["value"])


def _unmapped() -> pd.DataFrame:
    return fx.unmapped_frame([
        ("u1", "peptide", "unmapped"), ("u2", "peptide", "unmapped"),
        ("u3", "protein", "unmapped"),
        ("c1", "peptide", "contaminant"), ("c2", "peptide", "contaminant"),
        ("c3", "peptide", "contaminant"), ("c4", "protein", "contaminant"),
    ])  # fmt: skip


# --------------------------------------------------------------------------- example 6


def test_example_6_mismatch_rate_distinct_peptides() -> None:
    """Aligned (PR #87 decision 6): was ``id_mixups == 3``; 3 distinct flagged peptides of the 5
    distinct placed peptides (A, C, D, E, F) give ``engine_protein_mismatch_rate = 3 / 5``."""
    hits = fx.hits_frame([
        qc_hit("AAAAAAAK", flags=["engine_protein_mismatch"]),
        qc_hit("AAAAAAAK", flags=["engine_protein_mismatch", "same_db"]),  # same peptide twice
        qc_hit("CCCCCCCK", flags=["engine_protein_mismatch"]),
        qc_hit("DDDDDDDK", route="tolerant", flags=["engine_protein_mismatch"]),
        qc_hit("EEEEEEEK", flags=["same_db"]),
        qc_hit("FFFFFFFK"),
    ])  # fmt: skip
    peps = fx.peptides_frame(qc_peptide(p) for p in sorted(set(hits["peptide"])))
    out = _qc(peps, hits, _unmapped())
    assert _val(out, "engine_protein_mismatch_rate") == pytest.approx(3 / 5)
    assert math.isnan(_val(out, "id_mixups"))


# --------------------------------------------------------------------------- example 7


def test_example_7_genome_inconsistent_rate() -> None:
    """Aligned (PR #87 decision 5): the 100 peptides are ``search_db = other`` (was ``same``)."""
    names = [f"PEP{i:03d}K" for i in range(100)]
    hit_rows = [qc_hit(p) for p in names[:97]]
    hit_rows += [qc_hit(names[97], route="tolerant"), qc_hit(names[98], route="tolerant")]
    hit_rows.append(qc_hit(names[0], route="tolerant"))  # exact + tolerant: still exact
    pep_rows = [qc_peptide(p, search_db="other") for p in names]  # names[99] is absent from C11
    # ignored: decoys, same-database peptides
    pep_rows += [qc_peptide(f"DEC{i:03d}K", search_db="other", is_decoy=True) for i in range(20)]
    pep_rows += [qc_peptide(f"SAM{i:03d}K") for i in range(30)]
    out = _qc(fx.peptides_frame(pep_rows), fx.hits_frame(hit_rows), _unmapped())
    assert _val(out, "genome_inconsistent_rate") == pytest.approx(0.03)


def test_genome_inconsistent_rate_na_without_other_db() -> None:
    """Aligned (PR #87 decision 5): NA for ``search_db = same`` (was: NA for non-same)."""
    pep_rows = [qc_peptide(f"SAM{i:03d}K") for i in range(5)]
    hits = fx.hits_frame([qc_hit("SAM000K")])
    out = _qc(fx.peptides_frame(pep_rows), hits, _unmapped())
    assert math.isnan(_val(out, "genome_inconsistent_rate"))
    out2 = _qc(None, hits, _unmapped())
    assert math.isnan(_val(out2, "genome_inconsistent_rate"))


# --------------------------------------------------------------------------- C13 counts


def test_contaminant_and_unmapped_counts() -> None:
    hits = fx.hits_frame([qc_hit("AAAAAAAK")])
    out = _qc(fx.peptides_frame([qc_peptide("AAAAAAAK")]), hits, _unmapped())
    assert _val(out, "contaminant_peptides") == 3
    assert _val(out, "contaminant_proteins") == 1
    assert _val(out, "unmapped_peptides") == 2
    assert _val(out, "unmapped_proteins") == 1


def test_missing_c13_gives_na_counts() -> None:
    hits = fx.hits_frame([qc_hit("AAAAAAAK", flags=["engine_protein_mismatch"])])
    out = _qc(fx.peptides_frame([qc_peptide("AAAAAAAK")]), hits, None)
    for name in ("contaminant_peptides", "contaminant_proteins", "unmapped_peptides",
                 "unmapped_proteins"):  # fmt: skip
        assert math.isnan(_val(out, name)), name
    assert _val(out, "engine_protein_mismatch_rate") == 1  # aligned: was id_mixups == 1


def test_qc_rows_shape() -> None:
    hits = fx.hits_frame([qc_hit("AAAAAAAK")])
    out = _qc(fx.peptides_frame([qc_peptide("AAAAAAAK")]), hits, _unmapped())
    assert tuple(out.columns) == fx.SUMMARY_COLUMNS
    assert set(out["section"]) == {"qc"}
    assert set(out["group_type"]) == {"run"} and set(out["group"]) == {"all"}
    assert sorted(out["metric"]) == sorted(QC_METRICS)
    assert _val(out, "engine_protein_mismatch_rate") == 0  # aligned: was id_mixups == 0


def test_empty_hits_ok() -> None:
    """Aligned (PR #87 decisions 5-6): an ``other`` peptide (was ``same``); no C11 peptide gives
    an NA mismatch rate (was ``id_mixups == 0``)."""
    hits = fx.hits_frame([])
    out = _qc(fx.peptides_frame([qc_peptide("AAAAAAAK", search_db="other")]), hits, _unmapped())
    assert math.isnan(_val(out, "engine_protein_mismatch_rate"))
    assert _val(out, "genome_inconsistent_rate") == pytest.approx(1.0)


# --------------------------------------------------------------------------- detection_summary


def _parts() -> list[pd.DataFrame]:
    s = _summary()
    from giap.rng import for_step

    rows = [gene(f"G{i}", f"T{i // 2}", i % 2 == 0, 0.4, cog_category=["K", None][i % 2])
            for i in range(12)]  # fmt: skip
    basic = s.basic_rows(fx.genes_frame(rows), 100, for_step("context_stats", 1))
    hits = fx.hits_frame([qc_hit("AAAAAAAK")])
    qc = s.qc_rows(fx.peptides_frame([qc_peptide("AAAAAAAK")]), hits, _unmapped())
    det = fx.detectability_rows([i % 2 == 0 for i in range(40)], [0.5] * 40)
    tier = fx.tier_rows(6, 6)
    return [qc, det, basic, tier]  # deliberately out of order


def test_detection_summary_sections_in_order_and_sorted() -> None:
    out = _summary().detection_summary(_parts())
    assert tuple(out.columns) == fx.SUMMARY_COLUMNS
    sec = out["section"].astype(str).tolist()
    order = {s: i for i, s in enumerate(fx.SECTIONS)}
    assert [order[s] for s in sec] == sorted(order[s] for s in sec)
    assert set(sec) == set(fx.SECTIONS)
    for s in fx.SECTIONS:
        sub = out[out["section"] == s]
        key = list(zip(sub["group_type"].astype(str), sub["group"].astype(str),
                       sub["metric"].astype(str), strict=True))  # fmt: skip
        assert key == sorted(key), s


def test_detection_summary_keeps_every_row_and_value() -> None:
    parts = _parts()
    out = _summary().detection_summary(parts)
    assert len(out) == sum(len(p) for p in parts)
    det = out[out["section"] == "detectability"]
    slope = det[det["metric"] == "calibration_slope"].iloc[0]
    assert (slope["value"], slope["ci_low"], slope["ci_high"]) == (0.97, 0.81, 1.13)
    assert pd.isna(slope["n"])
    assert fx.flag_set(slope["flags"]) == set()
    basic = out[out["section"] == "basic"]
    assert (basic["n"].notna()).all()


def test_detection_summary_dtypes() -> None:
    out = _summary().detection_summary(_parts())
    for c in ("value", "ci_low", "ci_high"):
        assert pd.api.types.is_float_dtype(out[c]), c
    assert pd.api.types.is_integer_dtype(out["n"]), out["n"].dtype


def test_detection_summary_rejects_unknown_section() -> None:
    bad = fx.tier_rows(1, 1).assign(section="bogus")
    with pytest.raises(ValueError):
        _summary().detection_summary([bad])


def test_detection_summary_rejects_missing_core_column() -> None:
    bad = fx.tier_rows(1, 1).drop(columns=["metric"])
    with pytest.raises(ValueError):
        _summary().detection_summary([bad])


def test_detection_summary_deterministic() -> None:
    a = _summary().detection_summary(_parts())
    b = _summary().detection_summary(list(reversed(_parts())))
    pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True))
    assert np.array_equal(a.index.to_numpy(), np.arange(len(a)))
