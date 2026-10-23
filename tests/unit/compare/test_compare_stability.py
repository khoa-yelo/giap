"""M9-22 ``giap.compare.stability.stability``: Jaccard overlap of five result families between
runs (docs/specs/M9-22-compare-tables.md "compare_stability.tsv", examples 3 and 4, the
``function``-off degenerate row, lead decision 3; GIAP_SPEC.md §7.9.3).

File named ``test_compare_stability.py`` for unique basenames (tests/ has no ``__init__.py``).

Interface assumed (test-author reading of the spec): ``stability(runs, labels, cfg)`` with
``runs`` label -> ``RunView`` (built here with ``RunView.from_outdir`` on small ``OUT/`` trees),
returns one row per family x unordered pair (pairs in label order) with the
``compare_stability`` columns; ``flags`` is ``;``-joined (NA or empty when none). Jaccard is NA
(not 0) for ``both_empty`` and ``missing_in:<label>``.
"""

from __future__ import annotations

import math
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd
import pytest
from helpers.compare_data import (
    FAMILIES,
    STABILITY_COLUMNS,
    flag_set,
    gene,
    gsea,
    island,
    novel,
    operon_flag,
    ora,
    variant,
    write_run,
)
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig, resolve_config
from giap.io.run_view import RunView

FAST = settings(max_examples=25, deadline=None)
GENES = [gene(f"g{i}", start=100 * i) for i in range(1, 6)]
Tables = Mapping[str, Sequence[Mapping[str, object]]]


def _views(root: Path, tables: Mapping[str, Tables]) -> dict[str, RunView]:
    out = {}
    for label, tabs in tables.items():
        write_run(root / label, GENES, tables=tabs)
        out[label] = RunView.from_outdir(root / label)
    return out


def _stab(root: Path, tables: Mapping[str, Tables], cfg: GiapConfig | None = None) -> pd.DataFrame:
    from giap.compare.stability import stability

    views = _views(root, tables)
    return stability(views, list(tables), cfg or resolve_config())


def _get(frame: pd.DataFrame, family: str, a: str = "A", b: str = "B") -> pd.Series:
    sel = frame[(frame["family"] == family) & (frame["run_a"] == a) & (frame["run_b"] == b)]
    assert len(sel) == 1, (family, a, b, frame)
    return sel.iloc[0]


def _na(x: object) -> bool:
    return x is None or x is pd.NA or (isinstance(x, float) and math.isnan(x))


def _full(**over: Sequence[Mapping[str, object]]) -> dict[str, Sequence[Mapping[str, object]]]:
    """Every family table present (empty unless given)."""
    base: dict[str, Sequence[Mapping[str, object]]] = {
        k: [] for k in ("islands", "operon_flags", "enrichment_detection", "enrichment_gsea",
                        "novel_orfs", "variants")
    }  # fmt: skip
    base.update(over)
    return base


# --------------------------------------------------------------------- shape


def test_one_row_per_family_and_pair(tmp_path: Path) -> None:
    out = _stab(tmp_path, {"A": _full(), "B": _full(), "C": _full()})
    assert set(STABILITY_COLUMNS) <= set(out.columns)
    assert len(out) == len(FAMILIES) * 3
    assert set(out["family"]) == set(FAMILIES)
    pairs = list(dict.fromkeys(zip(out["run_a"], out["run_b"], strict=True)))
    assert pairs == [("A", "B"), ("A", "C"), ("B", "C")]


# ------------------------------------------------------------------- islands


def test_example3_island_gene_jaccard(tmp_path: Path) -> None:
    a = _full(islands=[island(["g1", "g2", "g3"], 0.01)])
    b = _full(islands=[island(["g2", "g3", "g4"], 0.01)])
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "islands")
    assert (int(row["n_a"]), int(row["n_b"]), int(row["n_shared"])) == (3, 3, 2)
    assert row["jaccard"] == pytest.approx(0.5)


def test_islands_union_over_islands_and_alpha_context(tmp_path: Path) -> None:
    """Genes of every island with p <= stats.alpha.context (0.05); the 0.2 island is out."""
    a = _full(islands=[island(["g1"], 0.01, "ISL_1"), island(["g2"], 0.05, "ISL_2"),
                       island(["g5"], 0.2, "ISL_3")])  # fmt: skip
    b = _full(islands=[island(["g1", "g2", "g5"], 0.04)])
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "islands")
    assert (int(row["n_a"]), int(row["n_shared"])) == (2, 2)
    assert row["jaccard"] == pytest.approx(2 / 3)


def test_islands_threshold_follows_config(tmp_path: Path) -> None:
    a = _full(islands=[island(["g1"], 0.01, "ISL_1"), island(["g5"], 0.2, "ISL_3")])
    b = _full(islands=[island(["g1", "g5"], 0.01)])
    cfg = resolve_config(overrides=["stats.alpha.context=0.25"])
    row = _get(_stab(tmp_path, {"A": a, "B": b}, cfg), "islands")
    assert row["jaccard"] == pytest.approx(1.0)


# --------------------------------------------------------------- operon flags


def test_operon_flag_pairs_with_alpha_operon_flags(tmp_path: Path) -> None:
    a = _full(operon_flags=[operon_flag("g1", "g2", 0.05), operon_flag("g3", "g4", 0.10),
                            operon_flag("g4", "g5", 0.5)])  # fmt: skip
    b = _full(operon_flags=[operon_flag("g1", "g2", 0.01), operon_flag("g4", "g5", 0.01)])
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "operon_flags")
    assert (int(row["n_a"]), int(row["n_b"]), int(row["n_shared"])) == (2, 2, 1)
    assert row["jaccard"] == pytest.approx(1 / 3)


# -------------------------------------------------------- enriched categories


def test_enriched_categories_key_includes_table_and_nes_sign(tmp_path: Path) -> None:
    a = _full(enrichment_detection=[ora("cog", "J", "over", 0.01), ora("cog", "K", "over", 0.2)],
              enrichment_gsea=[gsea("kegg", "ko00010", 2.0, 0.01)])  # fmt: skip
    b = _full(enrichment_detection=[ora("cog", "J", "over", 0.04)],
              enrichment_gsea=[gsea("kegg", "ko00010", -2.0, 0.01)])  # fmt: skip
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "enriched_categories")
    assert (int(row["n_a"]), int(row["n_b"]), int(row["n_shared"])) == (2, 2, 1)
    assert row["jaccard"] == pytest.approx(1 / 3)


def test_same_term_in_detection_and_gsea_are_distinct(tmp_path: Path) -> None:
    a = _full(enrichment_detection=[ora("cog", "J", "over", 0.01)])
    b = _full(enrichment_gsea=[gsea("cog", "J", 1.5, 0.01)])
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "enriched_categories")
    assert int(row["n_shared"]) == 0 and row["jaccard"] == pytest.approx(0.0)


def test_ora_direction_is_part_of_the_key(tmp_path: Path) -> None:
    a = _full(enrichment_detection=[ora("cog", "J", "over", 0.01)])
    b = _full(enrichment_detection=[ora("cog", "J", "under", 0.01)])
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "enriched_categories")
    assert int(row["n_shared"]) == 0


def test_enrichment_threshold_is_alpha_enrichment(tmp_path: Path) -> None:
    """Decision 3: ``stats.alpha.enrichment``, not a fixed 0.05."""
    a = _full(enrichment_detection=[ora("cog", "J", "over", 0.01), ora("cog", "K", "over", 0.08)])
    b = _full(enrichment_detection=[ora("cog", "J", "over", 0.01), ora("cog", "K", "over", 0.08)])
    default = _get(_stab(tmp_path / "d", {"A": a, "B": b}), "enriched_categories")
    assert int(default["n_a"]) == 1
    cfg = resolve_config(overrides=["stats.alpha.enrichment=0.1"])
    loose = _get(_stab(tmp_path / "l", {"A": a, "B": b}, cfg), "enriched_categories")
    assert int(loose["n_a"]) == 2 and loose["jaccard"] == pytest.approx(1.0)


def test_function_off_in_one_run_is_missing_in(tmp_path: Path) -> None:
    a = _full(enrichment_detection=[ora("cog", "J", "over", 0.01)])
    b = {k: v for k, v in _full().items() if not k.startswith("enrichment_")}
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "enriched_categories")
    assert _na(row["jaccard"])
    assert "missing_in:B" in flag_set(row["flags"])


# --------------------------------------------------------- novel ORFs, ALT sites


def test_novel_orfs_called_only(tmp_path: Path) -> None:
    a = _full(novel_orfs=[novel("N1", True), novel("N2", True), novel("N3", False)])
    b = _full(novel_orfs=[novel("N1", True), novel("N3", True)])
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "novel_orfs")
    assert (int(row["n_a"]), int(row["n_b"]), int(row["n_shared"])) == (2, 2, 1)
    assert row["jaccard"] == pytest.approx(1 / 3)


def test_example4_no_novel_orf_called_is_both_empty(tmp_path: Path) -> None:
    a = _full(novel_orfs=[novel("N1", False)])
    b = _full(novel_orfs=[])
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "novel_orfs")
    assert _na(row["jaccard"])
    assert "both_empty" in flag_set(row["flags"])
    assert (int(row["n_a"]), int(row["n_b"]), int(row["n_shared"])) == (0, 0, 0)


def test_alt_sites_key_and_call_alt_only(tmp_path: Path) -> None:
    a = _full(variants=[variant(30, "K"), variant(60, "R"), variant(90, "W", call="mixed")])
    b = _full(variants=[variant(30, "K"), variant(60, "Q"), variant(90, "W")])
    row = _get(_stab(tmp_path, {"A": a, "B": b}), "alt_sites")
    assert (int(row["n_a"]), int(row["n_b"]), int(row["n_shared"])) == (2, 3, 1)
    assert row["jaccard"] == pytest.approx(1 / 4)


@pytest.mark.parametrize(
    ("stem", "family"),
    [("islands", "islands"), ("operon_flags", "operon_flags"), ("novel_orfs", "novel_orfs"),
     ("variants", "alt_sites")],
)  # fmt: skip
def test_absent_table_is_missing_in_that_label(tmp_path: Path, stem: str, family: str) -> None:
    b = {k: v for k, v in _full().items() if k != stem}
    out = _stab(tmp_path, {"A": _full(), "B": b, "C": _full()})
    for a_, b_ in [("A", "B"), ("B", "C")]:
        row = _get(out, family, a_, b_)
        assert _na(row["jaccard"])
        assert "missing_in:B" in flag_set(row["flags"])
    assert "missing_in:B" not in flag_set(_get(out, family, "A", "C")["flags"])


def test_identical_runs_give_jaccard_one(tmp_path: Path) -> None:
    tabs = _full(islands=[island(["g1", "g2"], 0.01)],
                 operon_flags=[operon_flag("g1", "g2", 0.01)],
                 enrichment_detection=[ora("cog", "J", "over", 0.01)],
                 novel_orfs=[novel("N1", True)], variants=[variant(30, "K")])  # fmt: skip
    out = _stab(tmp_path, {"A": tabs, "B": tabs})
    assert list(out["jaccard"]) == pytest.approx([1.0] * len(FAMILIES))


# ------------------------------------------------------------------ property


@FAST
@given(
    st.sets(st.sampled_from([f"N{i}" for i in range(8)])),
    st.sets(st.sampled_from([f"N{i}" for i in range(8)])),
)
def test_property_novel_orf_jaccard_matches_set_formula(a: set[str], b: set[str]) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        ta = _full(novel_orfs=[novel(x, True) for x in sorted(a)])
        tb = _full(novel_orfs=[novel(x, True) for x in sorted(b)])
        row = _get(_stab(Path(tmp), {"A": ta, "B": tb}), "novel_orfs")
    assert (int(row["n_a"]), int(row["n_b"]), int(row["n_shared"])) == (
        len(a), len(b), len(a & b),
    )  # fmt: skip
    if not a and not b:
        assert _na(row["jaccard"]) and "both_empty" in flag_set(row["flags"])
    else:
        assert row["jaccard"] == pytest.approx(len(a & b) / len(a | b))
        assert 0.0 <= row["jaccard"] <= 1.0
