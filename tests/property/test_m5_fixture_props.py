"""Hypothesis properties of the M5-11 extension of ``gippasim.fixtures.build_m4_fixture``.

Spec "Interface": ``build_m4_fixture(..., oric_frac: float | None = 0.8, islands:
Sequence[IslandSpec] = ())`` -- the builder writes C15 ``features.parquet`` and can plant
islands; ``expression=`` (``ExpressionConfig`` fields, as ``layout=`` passes ``LayoutConfig``
fields) is the test author's addition for the island case's ``detection`` (tests/helpers/
m5_chain.py).

* the oriC is one 232-bp chr1 ``oric`` row near ``oric_frac * length`` and changes no other
  contract table (placing it draws nothing; the ori-ter gradient stays 0);
* the default is ``oric_frac = 0.8`` on the sequence-free sizes; ``None`` writes a header-only
  ``features.parquet``;
* planting an island changes the expression truth only: the layout tables (C1 replicons, genes,
  C5 tus, C7/C8 regions, C15 features) are byte-identical to the island-free build, the evidence
  is not (unless a down island covers only genes without peptides); an island outside its
  replicon raises ``ValueError``;
* the build is a pure function of its arguments (the fixture half of TE7).
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m5_chain as m5

SETTINGS = settings(max_examples=6, deadline=None, suppress_health_check=[HealthCheck.too_slow])
LAYOUT_FILES = ("replicons.parquet", "genes.parquet", "tus.parquet", "regions.parquet",
                "region_genes.parquet", "features.parquet")  # fmt: skip
ORIC_LEN = 232


def _build(work: Path, seed: int, **kw: Any) -> Path:
    from gippasim.fixtures import build_m4_fixture

    work.mkdir(parents=True, exist_ok=True)
    build_m4_fixture(work, seed=seed, size="sim1000", **kw)
    return work


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _all(work: Path) -> dict[str, str]:
    return {str(p.relative_to(work)): _sha(p) for p in sorted(work.rglob("*")) if p.is_file()}


def _oric(work: Path) -> pd.DataFrame:
    f = pd.read_parquet(work / "features.parquet")
    return f[f["type"].astype(str) == "oric"]


@SETTINGS
@given(seed=st.integers(0, 2**31 - 1), frac=st.sampled_from([0.05, 0.3, 0.5, 0.8, 0.95]))
def test_oric_is_one_chr1_feature_near_its_fraction(
    tmp_path_factory: pytest.TempPathFactory, seed: int, frac: float
) -> None:
    root = tmp_path_factory.mktemp("oric")
    with_oric = _build(root / "with", seed, oric_frac=frac)
    without = _build(root / "without", seed, oric_frac=None)
    o = _oric(with_oric)
    assert len(o) == 1 and str(o["contig"].iloc[0]) == "chr1"
    start, end = int(o["start"].iloc[0]), int(o["end"].iloc[0])
    assert end - start == ORIC_LEN
    length = int(pd.read_parquet(with_oric / "replicons.parquet")
                 .set_index("contig").loc["chr1", "length"])  # fmt: skip
    assert abs((start + end) / 2 / length - frac) < 0.02
    assert len(pd.read_parquet(without / "features.parquet")) == 0
    for name in LAYOUT_FILES:
        if name != "features.parquet":
            assert _sha(with_oric / name) == _sha(without / name), name
    for name in ("peptides.norm.parquet", "mapping/peptide_hits.tsv"):
        assert _sha(with_oric / name) == _sha(without / name), f"the oriC changed {name}"


def test_default_oric_frac_is_0_8(tmp_path: Path) -> None:
    default = _build(tmp_path / "default", 3)
    explicit = _build(tmp_path / "explicit", 3, oric_frac=0.8)
    assert len(_oric(default)) == 1
    assert _all(default) == _all(explicit)


@SETTINGS
@given(
    seed=st.integers(0, 2**31 - 1),
    first=st.integers(0, 400),
    n_tus=st.integers(2, 20),
    effect=st.sampled_from([-3.0, 1.5, 3.0]),
)
@example(seed=1506, first=0, n_tus=2, effect=-3.0)
@example(seed=356, first=0, n_tus=2, effect=1.5)
def test_islands_change_the_expression_not_the_layout(
    tmp_path_factory: pytest.TempPathFactory, seed: int, first: int, n_tus: int, effect: float
) -> None:
    from gippasim.expression import IslandSpec

    root = tmp_path_factory.mktemp("isl")
    plain = _build(root / "plain", seed)
    planted = _build(root / "planted", seed, islands=(IslandSpec("chr1", first, n_tus, effect),))
    for name in LAYOUT_FILES:
        assert _sha(plain / name) == _sha(planted / name), name
    if not _has_peptides(plain, first, n_tus):
        # an island over genes without peptides cannot change the evidence, up (seed 356) or
        # down (seed 1506): the default detection does not draw peptides from abundance
        return
    assert _sha(plain / "peptides.norm.parquet") != _sha(planted / "peptides.norm.parquet")


def _has_peptides(work: Path, first: int, n_tus: int) -> bool:
    tus = pd.read_parquet(work / "tus.parquet")
    chr1 = tus[tus["contig"].astype(str) == "chr1"].iloc[first : first + n_tus]
    genes = {str(g) for gs in chr1["genes"] for g in gs}
    hits = pd.read_csv(work / "mapping" / "peptide_hits.tsv", sep="\t", usecols=["gene"])
    return bool(hits["gene"].astype(str).isin(genes).any())


def test_an_island_outside_its_replicon_raises(tmp_path: Path) -> None:
    from gippasim.expression import IslandSpec

    with pytest.raises(ValueError, match="island"):
        _build(tmp_path / "w", 1, islands=(IslandSpec("chr1", 10_000, 5, 3.0),))


def test_expression_overrides_keep_the_layout(tmp_path: Path) -> None:
    """``expression=`` (the island case's ``detection = "abundance"``) changes evidence only."""
    case = m5.PROPOSED_ISLAND
    plain = _build(tmp_path / "plain", 5)
    expr = _build(tmp_path / "expr", 5, expression=case.expression())
    for name in LAYOUT_FILES:
        assert _sha(plain / name) == _sha(expr / name), name
    assert _sha(plain / "peptides.norm.parquet") != _sha(expr / "peptides.norm.parquet")


@SETTINGS
@given(seed=st.integers(0, 2**31 - 1))
def test_the_island_build_is_a_pure_function(
    tmp_path_factory: pytest.TempPathFactory, seed: int
) -> None:
    root = tmp_path_factory.mktemp("pure")
    case = m5.IslandCase("sim1000_probe", "sim1000", 1000, 20, 3.0, "abundance")
    a = m5.build_island(root / "a", seed, case)
    b = m5.build_island(root / "b", seed, case)
    assert _all(a) == _all(b)
    assert len(m5.planted_tus(a, case)) == 20
