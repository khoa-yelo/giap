"""M9-22 ``giap.compare.detection``: the gene x run detection matrix and the exact (UpSet)
intersections (docs/specs/M9-22-compare-tables.md "Tables", example 1, lead decision 4;
GIAP_SPEC.md §7.9.3, §6.5 sort order).

File named ``test_compare_detection.py`` (not ``test_detection.py`` as in the acceptance list):
tests/ has no ``__init__.py`` files, so basenames must be unique across the tree.

Interface assumed (test-author reading of the spec):

- ``detection_matrix(runs)``: ``runs`` maps label -> genes frame (``RunView.table("genes_tsv")``
  shape) in label order; returns the **matrix section** in long form with the
  ``compare_detection`` columns: ``section = "matrix"``, one row per gene x run, ``detected`` and
  ``detection_status`` copied, ``intersection`` / ``n_genes`` NA; sorted by gene in §6.5 order
  (contig as first seen, start, end, strand, ID), then run in label order.
- ``intersections(matrix, labels)``: the intersection section: ``section = "intersection"``,
  ``intersection`` (labels joined by ``&`` in label order, or ``none``), ``n_genes``;
  ``locus_tag``, ``run``, ``detected`` NA; sorted by ``n_genes`` descending then name.
"""

from __future__ import annotations

from collections import Counter

import pandas as pd
import pytest
from helpers.compare_data import DETECTION_COLUMNS, gene, genes_frame
from hypothesis import given, settings
from hypothesis import strategies as st

FAST = settings(max_examples=40, deadline=None)
GENES = ["g1", "g2", "g3", "g4", "g5"]
LABELS = ["A1", "A2", "A3"]


def _runs(detected: dict[str, set[str]], genes: list[str] = GENES) -> dict[str, pd.DataFrame]:
    return {
        label: genes_frame(
            [gene(g, start=1000 * i, detected=g in hits) for i, g in enumerate(genes)]
        )
        for label, hits in detected.items()
    }


def _example1() -> dict[str, pd.DataFrame]:
    return _runs({"A1": {"g1", "g2", "g3"}, "A2": {"g2", "g3", "g4"}, "A3": {"g3"}})


def _inter(runs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    from giap.compare.detection import detection_matrix, intersections

    return intersections(detection_matrix(runs), list(runs))


def _counts(inter: pd.DataFrame) -> dict[str, int]:
    return {str(r.intersection): int(r.n_genes) for r in inter.itertuples()}


# ------------------------------------------------------------------- matrix


def test_matrix_has_one_row_per_gene_and_run() -> None:
    from giap.compare.detection import detection_matrix

    m = detection_matrix(_example1())
    assert len(m) == len(GENES) * len(LABELS)
    assert set(DETECTION_COLUMNS) <= set(m.columns)
    assert set(m["section"]) == {"matrix"}
    assert m["intersection"].isna().all() and m["n_genes"].isna().all()


def test_matrix_copies_detected_and_status() -> None:
    from giap.compare.detection import detection_matrix

    m = detection_matrix(_example1())
    got = {(r.locus_tag, r.run): bool(r.detected) for r in m.itertuples()}
    assert got[("g1", "A1")] and not got[("g1", "A2")] and got[("g3", "A3")]
    assert not any(got[("g5", lab)] for lab in LABELS)
    status = {(r.locus_tag, r.run): r.detection_status for r in m.itertuples()}
    assert status[("g4", "A2")] == "detected" and status[("g4", "A1")] == "not_detected"


def test_matrix_sorted_genes_by_position_then_run_in_label_order() -> None:
    """§6.5: by contig (as first seen), then start; labels in the given order (not sorted)."""
    from giap.compare.detection import detection_matrix

    # positions disagree with the lexicographic order of IDs and labels
    names = ["z9", "a1", "m5"]
    rows = [gene(n, start=s) for n, s in zip(names, [10, 5000, 900], strict=True)]
    runs = {"B": genes_frame(rows), "A": genes_frame(list(reversed(rows)))}
    m = detection_matrix(runs)
    assert list(m["locus_tag"]) == ["z9", "z9", "m5", "m5", "a1", "a1"]
    assert list(m["run"]) == ["B", "A"] * 3


def test_matrix_contig_order_is_first_seen_not_alphabetical() -> None:
    from giap.compare.detection import detection_matrix

    rows = [gene("p1", start=0, contig="plasmid"), gene("c1", start=0, contig="chrA"),
            gene("p0", start=50, contig="plasmid")]  # fmt: skip
    runs = {"A": genes_frame(rows), "B": genes_frame(rows)}
    assert list(detection_matrix(runs)["locus_tag"])[::2] == ["p1", "p0", "c1"]


# ------------------------------------------------------------- intersections


def test_example1_exact_intersections() -> None:
    inter = _inter(_example1())
    assert set(inter["section"]) == {"intersection"}
    assert _counts(inter) == {"A1&A2&A3": 1, "A1": 1, "A1&A2": 1, "A2": 1, "none": 1}
    assert inter["locus_tag"].isna().all()
    assert inter["run"].isna().all()
    assert inter["detected"].isna().all()


def test_example1_ordering_by_count_then_name() -> None:
    """All counts tie at 1, so the order is by name (``none`` placement not pinned)."""
    names = [n for n in _inter(_example1())["intersection"] if n != "none"]
    assert names == sorted(names)


def test_intersections_sorted_by_count_descending() -> None:
    # A only: 3 genes; A&B: 2; B only: 1; none: 0 or absent
    genes = [f"g{i}" for i in range(6)]
    runs = _runs({"A": {"g0", "g1", "g2", "g3", "g4"}, "B": {"g3", "g4", "g5"}}, genes)
    inter = _inter(runs)
    real = inter[inter["intersection"] != "none"]
    assert list(real["intersection"]) == ["A", "A&B", "B"]
    assert list(real["n_genes"]) == [3, 2, 1]


def test_intersection_names_follow_label_order_not_alphabet() -> None:
    runs = _runs({"Zed": {"g1"}, "Alpha": {"g1"}}, ["g1", "g2"])
    assert "Zed&Alpha" in _counts(_inter(runs))
    assert "Alpha&Zed" not in _counts(_inter(runs))


def test_intersections_are_exact_not_inclusive() -> None:
    """Decision 4: a gene detected in A1 and A3 only belongs to ``A1&A3`` and nowhere else."""
    runs = _runs({"A1": {"g1"}, "A2": set(), "A3": {"g1"}}, ["g1", "g2"])
    counts = _counts(_inter(runs))
    assert counts.get("A1&A3") == 1
    assert "A1" not in counts and "A3" not in counts
    assert counts.get("none") == 1


def test_empty_intersections_are_not_listed() -> None:
    counts = _counts(_inter(_example1()))
    assert "A3" not in counts and "A2&A3" not in counts and "A1&A3" not in counts
    assert all(v > 0 for k, v in counts.items() if k != "none")


@FAST
@given(st.lists(st.tuples(st.booleans(), st.booleans(), st.booleans()), min_size=1, max_size=25))
def test_property_intersections_partition_the_genes(hits: list[tuple[bool, ...]]) -> None:
    """Every gene lies in exactly one exact intersection (``none`` included)."""
    genes = [f"g{i:02d}" for i in range(len(hits))]
    runs = _runs({lab: {g for g, h in zip(genes, hits, strict=True) if h[j]}
                  for j, lab in enumerate(LABELS)}, genes)  # fmt: skip
    counts = _counts(_inter(runs))
    assert sum(counts.values()) == len(genes)
    oracle = Counter(
        "&".join(lab for lab, x in zip(LABELS, h, strict=True) if x) or "none" for h in hits
    )
    assert {k: v for k, v in counts.items() if v} == dict(oracle)
    for name in counts:
        if name != "none":
            parts = name.split("&")
            assert parts == [lab for lab in LABELS if lab in parts]  # label order


@FAST
@given(st.lists(st.tuples(st.booleans(), st.booleans()), min_size=1, max_size=20))
def test_property_matrix_counts_match_detected_sets(hits: list[tuple[bool, bool]]) -> None:
    from giap.compare.detection import detection_matrix

    genes = [f"g{i:02d}" for i in range(len(hits))]
    runs = _runs({"A": {g for g, h in zip(genes, hits, strict=True) if h[0]},
                  "B": {g for g, h in zip(genes, hits, strict=True) if h[1]}}, genes)  # fmt: skip
    m = detection_matrix(runs)
    for j, lab in enumerate(["A", "B"]):
        sub = m[m["run"] == lab]
        assert int(sub["detected"].astype(bool).sum()) == sum(h[j] for h in hits)


@pytest.mark.parametrize("n_runs", [2, 4])
def test_matrix_works_for_any_run_count(n_runs: int) -> None:
    labels = [f"R{i}" for i in range(n_runs)]
    runs = _runs({lab: {"g1"} for lab in labels})
    counts = _counts(_inter(runs))
    assert counts["&".join(labels)] == 1
    assert counts["none"] == len(GENES) - 1
