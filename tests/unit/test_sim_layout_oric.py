"""giapsim layout oriC feature and the C15 ``features`` writer (issue M5-04; GIAP_SPEC.md §11.3).

Examples 1, 2 and 7, the oriC edge cases and the determinism criterion of
docs/specs/M5-04-sim-oric-gradient.md, checked against tests/helpers/oric_oracle.py (written from
the spec). Lives in tests/unit like the M4-07/M4-08 simulator tests (tests/sim is implementer
tooling); the spec's ``tests/sim/test_layout_oric.py``.

Interpretations pinned here (reported to the lead):
- "the intergenic gap of chr1 whose midpoint is nearest to ``oric_frac x length``" is taken among
  the gaps the 232-bp feature fits in (width >= 232; most inter-TU gaps are 50-300 bp), the wrap gap
  across the origin included; the feature is centred in that gap (+/- 1 bp for odd widths);
- placing the oriC moves no gene: every other frame is identical to ``oric_frac = None``;
- ``feature_id`` is ``<contig>_oric1``, i.e. ``chr1_oric1`` (lead decision 15);
- no gap of width >= 232 -> ``ValueError`` naming ``oric_frac`` (lead decision 14).
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import layout_oracle as lo
import oric_oracle as oo

# Contract tables of ``simulate_layout(LayoutConfig(), default_rng(1))`` before M5-04 (main at
# ecd3f6f; equal to the M4-08 pins in test_sim_expression.py).
PRE_M5_DEFAULT_SHA256 = {
    "annotation_meta.json": "574ce6a5c6468c793385646f795d8986c317903210e7e7e2f585a14293a30517",
    "genes.parquet": "e0d4341901746a7444965ebe8070a3c382d5edc6f864bff8597bd8e697b65e65",
    "layers.parquet": "c7ae231f33b641c11a6b7b03ce95c36e5daf5c6b7fd5c1a46456d31e6d6c5ba5",
    "region_genes.parquet": "1197b721b0b595beaf182adb375dfbf0707c5cf8085b5d4598106988f1215fd7",
    "regions.parquet": "54ca39d8465b7720f4417364363d63ce4c6aa92cd893cad24c2f311c825f2ca5",
    "replicons.parquet": "baa9f94adf047f02a746031e33777abc727b0eea9e650454e0b7a868562bd85e",
    "tus.parquet": "c564de32b5fac4d9d1dc6a9ff5f885065ccdec0c2d7547b7ce48561f0e1ca825",
    "xrefs.parquet": "0fd3d2531b4d94b2dc6370de0378ecb870858dc2f6b091c37854771e474f592e",
}
OTHER_FRAMES = ("replicons", "genes", "tus", "operon_pairs", "regions", "region_genes", "truth",
                "xrefs")  # fmt: skip
ARROW_TYPES = {
    "feature_id": pa.string(), "type": pa.string(), "contig": pa.string(), "start": pa.int64(),
    "end": pa.int64(), "strand": pa.string(), "gene": pa.string(), "product": pa.string(),
    "amino_acid": pa.string(), "anticodon": pa.string(),
}  # fmt: skip


def simulate(seed: int = 1, **kw: Any) -> Any:
    from giapsim.layout import LayoutConfig, simulate_layout

    return simulate_layout(LayoutConfig(**kw), np.random.default_rng(seed))


def write(layout: Any, work: Path) -> None:
    from giapsim.contract_tables import write_contract_tables

    write_contract_tables(layout, work)


def chr1_length(layout: Any) -> int:
    return int(layout.replicons.set_index("contig").loc["chr1", "length"])


def oric_row(layout: Any) -> pd.Series:
    f = layout.features
    rows = f[f["type"] == "oric"]
    assert len(rows) == 1, f
    return rows.iloc[0]


def digests(work: Path) -> dict[str, str]:
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(work.iterdir())}


# ----------------------------------------------------------------------------- interface


def test_layout_config_oric_frac_defaults_to_none() -> None:
    from giapsim.layout import LayoutConfig

    assert LayoutConfig().oric_frac is None
    cfg = LayoutConfig(oric_frac=0.5)
    assert cfg.oric_frac == 0.5
    hash(cfg)


def test_contract_types_have_c15_with_every_column_written() -> None:
    """giapsim's literal C15 copy: the ten parse_bakta columns, none optional (gene, product,
    amino_acid and anticodon are nullable, not absent)."""
    from giapsim.contract_tables import CONTRACT_COLUMNS, CONTRACT_TYPES, REQUIRED_COLUMNS

    assert tuple(CONTRACT_COLUMNS["features"]) == oo.FEATURE_COLUMNS
    assert tuple(REQUIRED_COLUMNS["features"]) == oo.FEATURE_COLUMNS
    types = dict(CONTRACT_TYPES["features"])
    assert {k: v.rstrip("?") for k, v in types.items()} == {
        "feature_id": "str", "type": "str", "contig": "str", "start": "int", "end": "int",
        "strand": "str", "gene": "str", "product": "str", "amino_acid": "str", "anticodon": "str",
    }  # fmt: skip


# ----------------------------------------------------------------------------- example 1


def test_example1_default_layout_has_empty_features_with_c15_columns() -> None:
    lay = simulate(1)
    assert isinstance(lay.features, pd.DataFrame)
    assert len(lay.features) == 0
    assert tuple(lay.features.columns) == oo.FEATURE_COLUMNS


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="pinned Parquet bytes (Linux)")
def test_example1_pre_m5_outputs_byte_identical(tmp_path: Path) -> None:
    write(simulate(1), tmp_path / "w")
    d = digests(tmp_path / "w")
    assert {k: d.get(k) for k in PRE_M5_DEFAULT_SHA256} == PRE_M5_DEFAULT_SHA256
    assert "features.parquet" in d


def test_example1_features_parquet_empty_with_schema(tmp_path: Path) -> None:
    write(simulate(2, n_genes_chromosome=200), tmp_path / "w")
    t = pq.read_table(tmp_path / "w" / "features.parquet")
    assert t.num_rows == 0
    assert t.column_names == list(oo.FEATURE_COLUMNS)
    assert {f.name: f.type for f in t.schema} == ARROW_TYPES


# ----------------------------------------------------------------------------- example 2


@pytest.fixture(scope="module")
def half() -> Any:
    return simulate(1, n_genes_chromosome=1000, oric_frac=0.5)


def test_example2_one_oric_row(half: Any) -> None:
    f = half.features
    assert tuple(f.columns) == oo.FEATURE_COLUMNS
    assert len(f) == 1
    row = oric_row(half)
    assert row["type"] == "oric" and row["contig"] == "chr1" and row["strand"] == "."
    assert int(row["end"]) - int(row["start"]) == oo.ORIC_LEN
    assert row["feature_id"] == "chr1_oric1"  # lead decision 15
    n = chr1_length(half)
    assert 0 <= int(row["start"]) < n and int(row["end"]) <= 2 * n
    for c in ("amino_acid", "anticodon"):
        assert pd.isna(row[c]), c


def _assert_no_gene_overlap(lay: Any) -> None:
    row = oric_row(lay)
    n = chr1_length(lay)
    iv = (int(row["start"]), int(row["end"]))
    g = lay.genes[lay.genes["contig"] == "chr1"]
    bad = [t for t, s, e in zip(g["locus_tag"], g["start"], g["end"], strict=True)
           if oo.overlaps_circular(iv, (int(s), int(e)), n)]  # fmt: skip
    assert bad == [], (iv, bad)


def _enclosing_gap(lay: Any) -> tuple[int, int]:
    """The chr1 gap (left end, right start; unwrapped) that contains the oriC."""
    row = oric_row(lay)
    n = chr1_length(lay)
    s, e = int(row["start"]), int(row["end"])
    for a, b in oo.chr1_gaps(lay.genes, n):
        for k in (-1, 0, 1):
            if a <= s + k * n and e + k * n <= b:
                return a, b
    raise AssertionError(f"oriC [{s}, {e}) lies in no chr1 gap")


def test_example2_no_gene_overlap(half: Any) -> None:
    _assert_no_gene_overlap(half)


def test_example2_centred_in_its_gap(half: Any) -> None:
    a, b = _enclosing_gap(half)
    row = oric_row(half)
    n = chr1_length(half)
    assert oo.circ_dist((int(row["start"]) + int(row["end"])) / 2, (a + b) / 2, n) <= 1.0


def _check_nearest(lay: Any, frac: float) -> None:
    n = chr1_length(lay)
    target = frac * n
    a, b = _enclosing_gap(lay)
    chosen = oo.circ_dist((a + b) / 2, target, n)
    fits = [(x, y) for x, y in oo.chr1_gaps(lay.genes, n) if y - x >= oo.ORIC_LEN]
    best = min(oo.circ_dist((x + y) / 2, target, n) for x, y in fits)
    assert chosen <= best + 1.0, (chosen, best)


def test_example2_gap_is_the_nearest_that_fits(half: Any) -> None:
    _check_nearest(half, 0.5)


def test_example2_midpoint_near_target(half: Any) -> None:
    """'Within one inter-TU gap of 0.5 L': the nearest fitting gap is a few TUs away at most."""
    *_, o = oo.oric_of(half) or (0, 0, math.nan)
    n = chr1_length(half)
    assert oo.circ_dist(o, 0.5 * n, n) <= 0.02 * n


@pytest.mark.parametrize("frac", [0.0, 0.1, 0.25, 0.8, 0.999])
@pytest.mark.parametrize("seed", [3, 4])
def test_placement_rules_at_any_fraction(seed: int, frac: float) -> None:
    lay = simulate(seed, n_genes_chromosome=600, oric_frac=frac)
    _assert_no_gene_overlap(lay)
    _check_nearest(lay, frac)
    row = oric_row(lay)
    assert int(row["end"]) - int(row["start"]) == oo.ORIC_LEN
    assert 0 <= int(row["start"]) < chr1_length(lay)  # reduced past the origin (B1)


def test_oric_moves_no_gene_and_changes_no_other_frame() -> None:
    """Placing the oriC draws nothing and moves nothing: the other frames equal oric_frac=None."""
    kw: dict[str, Any] = dict(n_genes_chromosome=400, plasmids=(30,), linear_contigs=(20,))
    a, b = simulate(9, **kw), simulate(9, oric_frac=0.3, **kw)
    for name in OTHER_FRAMES:
        pd.testing.assert_frame_equal(getattr(a, name), getattr(b, name), obj=name)


def test_oric_only_on_chr1_with_plasmids_and_contigs() -> None:
    lay = simulate(5, n_genes_chromosome=300, plasmids=(40,), linear_contigs=(25,), oric_frac=0.6)
    assert lay.features["contig"].tolist() == ["chr1"]


# ------------------------------------------------------- _place_oric units (review B1, S3)


def place(frac: float, genes: list[tuple[int, int]], length: int) -> tuple[int, int]:
    """``(start, end)`` of ``_place_oric`` for chr1 genes ``(start, end)`` on a circle."""
    from giapsim.layout import _place_oric

    contig = ["chr1"] * len(genes)
    row = _place_oric(frac, contig, np.array([s for s, _ in genes]),
                      np.array([e for _, e in genes]), length)  # fmt: skip
    assert len(row) == 1 and tuple(row.columns) == oo.FEATURE_COLUMNS
    return int(row["start"].iloc[0]), int(row["end"].iloc[0])


def test_place_oric_reduces_a_start_past_the_origin() -> None:
    """Review B1: the wrap gap (9,990 ... 10,400) centres the oriC at 10,079, which must be
    reduced to the circle: (79, 311), never (10,079, 10,311)."""
    assert place(0.0, [(400, 900), (1000, 9990)], 10_000) == (79, 311)


@pytest.mark.parametrize("first_start", [240, 300, 420, 999])
def test_place_oric_start_is_always_inside_the_replicon(first_start: int) -> None:
    """Whatever the wrap gap's width, ``0 <= start < length`` and the feature stays in the gap."""
    length = 10_000
    genes = [(first_start, first_start + 500), (2000, 9950)]
    start, end = place(0.0, genes, length)
    assert 0 <= start < length
    assert end - start == oo.ORIC_LEN
    assert oo.overlaps_circular((start, end), (9950, first_start + length), length)
    for g in genes:
        assert not oo.overlaps_circular((start, end), g, length)


def test_place_oric_accepts_a_gap_of_exactly_232_bp() -> None:
    """S3: ``>= ORIC_LEN`` (not ``>``): a 232-bp gap holds the oriC exactly, with no slack."""
    length = 20_000
    genes = [(1000, 5000), (5232, 9000), (12_000, 19_900)]  # the 5,000-5,232 gap is exactly 232
    assert place(0.25, genes, length) == (5000, 5232)


def test_place_oric_just_under_232_bp_is_not_used() -> None:
    length = 20_000
    genes = [(1000, 5000), (5231, 9000), (12_000, 19_900)]
    start, _ = place(0.25, genes, length)
    assert start != 5000


def test_place_oric_ties_take_the_first_gap() -> None:
    """S3: two gaps the same circular distance from the target -> the first one (lowest start)."""
    length = 10_000
    # gaps 2,000-2,400 (midpoint 2,200) and 7,600-8,000 (midpoint 7,800): both 2,800 from 5,000
    genes = [(100, 2000), (2400, 7600), (8000, 9900)]
    assert place(0.5, genes, length) == (2084, 2316)


# ----------------------------------------------------------------------------- edge cases


@pytest.mark.parametrize("frac", [-0.1, 1.0, 1.5, math.nan, math.inf])
def test_oric_frac_outside_unit_interval_is_value_error(frac: float) -> None:
    from giapsim.layout import LayoutConfig, simulate_layout

    with pytest.raises(ValueError, match="oric_frac"):
        simulate_layout(LayoutConfig(n_genes_chromosome=100, oric_frac=frac),
                        np.random.default_rng(1))  # fmt: skip


# One-TU, three-gene chromosomes whose every gap (wrap gap included) is < 232 bp (seeds found on
# main at ecd3f6f; the test re-checks the widths so a layout change cannot make it vacuous).
DENSE_KW: dict[str, Any] = {"n_genes_chromosome": 3, "mean_tu_size": 1000.0, "regions": ()}
DENSE_SEEDS = (6, 17, 22)


@pytest.mark.parametrize("seed", DENSE_SEEDS)
def test_no_fitting_gap_is_value_error(seed: int) -> None:
    """Lead decision 14: when no chr1 gap is >= 232 bp, ``ValueError`` naming ``oric_frac``."""
    lay = simulate(seed, **DENSE_KW)
    widths = [b - a for a, b in oo.chr1_gaps(lay.genes, chr1_length(lay))]
    assert max(widths) < oo.ORIC_LEN, widths
    with pytest.raises(ValueError, match="oric_frac"):
        simulate(seed, oric_frac=0.5, **DENSE_KW)


def test_feature_id_is_contig_oric1_for_every_layout() -> None:
    """Lead decision 15: ``<contig>_oric1``, also with plasmids, contigs and origin crossing."""
    for seed in (1, 2):
        lay = simulate(seed, n_genes_chromosome=300, plasmids=(20,), linear_contigs=(10,),
                       origin_crossing=True, oric_frac=0.2)  # fmt: skip
        assert lay.features["feature_id"].tolist() == ["chr1_oric1"]


@pytest.mark.parametrize("frac", [0.0, 0.001, 0.999])
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_origin_crossing_with_oric_near_zero(seed: int, frac: float) -> None:
    """Both options: the oriC still sits in a gap; the crossing gene keeps ``end > length``."""
    lay = simulate(seed, n_genes_chromosome=500, origin_crossing=True, oric_frac=frac)
    n = chr1_length(lay)
    cross = lay.genes[lay.genes["crosses_origin"]]
    assert len(cross) == 1 and int(cross["end"].iloc[0]) > n
    _assert_no_gene_overlap(lay)
    _enclosing_gap(lay)
    row = oric_row(lay)
    assert 0 <= int(row["start"]) < n and int(row["end"]) <= 2 * n


# ----------------------------------------------------------------------------- example 7


def test_example7_features_parquet_written_and_valid(tmp_path: Path, half: Any) -> None:
    from giap.genetable.contract import validate_inputs, validate_table

    work = tmp_path / "w"
    write(half, work)
    assert sorted(p.name for p in work.iterdir()) == sorted(lo.M2_FILES)
    t = pq.read_table(work / "features.parquet")
    assert t.column_names == list(oo.FEATURE_COLUMNS)
    assert {f.name: f.type for f in t.schema} == ARROW_TYPES
    df = t.to_pandas()
    assert len(df) == 1
    row = oric_row(half)
    for c in ("feature_id", "type", "contig", "start", "end", "strand"):
        assert df[c].iloc[0] == row[c], c
    assert validate_table("features", df) == []
    inputs, issues = lo.load_work(work)
    assert lo.errors(issues) == []
    assert validate_inputs(_WithFeatures(inputs, df)) == []


def test_example7_without_features_file_no_issue(tmp_path: Path, half: Any) -> None:
    """Optional: the other tables load and cross-validate the same without features.parquet."""
    from giap.genetable.contract import validate_inputs

    work = tmp_path / "w"
    write(half, work)
    (work / "features.parquet").unlink()
    inputs, issues = lo.load_work(work)
    assert lo.errors(issues) == []
    assert validate_inputs(_WithFeatures(inputs, None)) == []


class _WithFeatures:
    def __init__(self, base: Any, features: pd.DataFrame | None) -> None:
        self._base = base
        self.features = features

    def __getattr__(self, name: str) -> Any:
        return getattr(self._base, name)


# ----------------------------------------------------------------------------- determinism


def test_same_seed_byte_identical_features_parquet(tmp_path: Path) -> None:
    kw: dict[str, Any] = dict(n_genes_chromosome=300, plasmids=(20,), oric_frac=0.7)
    np.random.seed(0)  # noqa: NPY002 (the global numpy state must not matter)
    write(simulate(11, **kw), tmp_path / "a")
    np.random.seed(99)  # noqa: NPY002
    write(simulate(11, **kw), tmp_path / "b")
    da, db = digests(tmp_path / "a"), digests(tmp_path / "b")
    assert "features.parquet" in da
    assert da == db


def test_different_seed_moves_the_oric() -> None:
    a = simulate(11, n_genes_chromosome=300, oric_frac=0.7)
    b = simulate(12, n_genes_chromosome=300, oric_frac=0.7)
    assert int(oric_row(a)["start"]) != int(oric_row(b)["start"])


# ----------------------------------------------------------------------------- M0 adapter


@pytest.fixture(scope="module")
def m0_oric(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, Any]:
    import giapsim
    from giapsim.layout import layout_from_simulation

    sim = giapsim.simulate(giapsim.SimConfig(oric=True, crispr=True), 3,
                           tmp_path_factory.mktemp("m0") / "sim")  # fmt: skip
    return sim, layout_from_simulation(sim)


def test_adapter_fills_features_from_bakta_like_json(m0_oric: tuple[Any, Any]) -> None:
    """Edge case: the M0 simulator's ``oriC`` (and ``crispr``) JSON features become C15 rows
    (0-based half-open; Bakta's strand ``?`` as ``.``)."""
    sim, lay = m0_oric
    data = json.loads(Path(sim.bakta_json).read_text(encoding="utf-8"))
    length = {s["id"]: int(s["length"]) for s in data["sequences"]}
    want = []
    for f in data["features"]:
        if f["type"] not in ("oriC", "crispr"):
            continue
        cid = f.get("sequence", f.get("contig"))
        s0, stop = int(f["start"]) - 1, int(f["stop"])
        end = stop + length[cid] if stop <= s0 else stop
        typ = "oric" if f["type"] == "oriC" else "crispr"
        strand = f["strand"] if f["strand"] in ("+", "-") else "."
        want.append((typ, cid, s0, end, strand))
    assert sorted(t for t, *_ in want) == ["crispr", "oric"]
    f = lay.features
    assert tuple(f.columns) == oo.FEATURE_COLUMNS
    got = sorted(zip(f["type"], f["contig"], f["start"].astype(int), f["end"].astype(int),
                     f["strand"], strict=True))  # fmt: skip
    assert got == sorted(want)
    assert f["feature_id"].is_unique


def test_adapter_features_valid_in_contract(m0_oric: tuple[Any, Any], tmp_path: Path) -> None:
    from giap.genetable.contract import validate_table

    _, lay = m0_oric
    write(lay, tmp_path / "w")
    df = pd.read_parquet(tmp_path / "w" / "features.parquet")
    assert len(df) == len(lay.features) == 2
    assert validate_table("features", df) == []


def _fake_sim(tmp_path: Path, features: list[dict[str, Any]], length: int = 6000) -> Any:
    """The two files ``layout_from_simulation`` reads: a bakta-like JSON and the gene truth."""
    from types import SimpleNamespace

    genes = [
        {"type": "cds", "sequence": "chr1", "start": 101, "stop": 400, "strand": "+",
         "locus": "G_00010", "aa": "MKV", "nt": "ATG" + "GCT" * 98 + "TAA"},
        {"type": "cds", "sequence": "chr1", "start": 1001, "stop": 1300, "strand": "-",
         "locus": "G_00020", "aa": "MKV", "nt": "ATG" + "GCT" * 98 + "TAA"},
    ]  # fmt: skip
    data = {
        "sequences": [{"id": "chr1", "type": "chromosome", "topology": "circular",
                       "length": length}],
        "features": [*genes, *features],
    }  # fmt: skip
    tmp_path.mkdir(parents=True, exist_ok=True)
    js, truth = tmp_path / "bakta_like.json", tmp_path / "genes.tsv"
    js.write_text(json.dumps(data), encoding="utf-8")
    truth.write_text("locus_tag\ttu_id\nG_00010\tTU1\nG_00020\tTU2\n", encoding="utf-8")
    return SimpleNamespace(bakta_json=js, genes_truth=truth)


def test_adapter_unwraps_an_origin_crossing_feature(tmp_path: Path) -> None:
    """S3: a Bakta feature with ``stop <= start`` crosses the origin; C15 keeps ``end``
    unwrapped (``end = stop + length``, §6.5), so ``length < end < 2 length``."""
    from giapsim.layout import layout_from_simulation

    sim = _fake_sim(tmp_path / "m0", [{"type": "oriC", "sequence": "chr1", "start": 5901,
                                       "stop": 132, "strand": "?", "id": "O1"}])  # fmt: skip
    f = layout_from_simulation(sim).features
    assert len(f) == 1
    row = f.iloc[0]
    assert (row["type"], row["strand"]) == ("oric", ".")
    assert (int(row["start"]), int(row["end"])) == (5900, 6132)  # not (5900, 132)
    assert 6000 < int(row["end"]) < 12_000
    assert int(row["end"]) - int(row["start"]) == oo.ORIC_LEN


def test_adapter_keeps_a_non_crossing_feature_unchanged(tmp_path: Path) -> None:
    from giapsim.layout import layout_from_simulation

    sim = _fake_sim(tmp_path / "m0", [{"type": "oriC", "sequence": "chr1", "start": 2001,
                                       "stop": 2232, "strand": "?", "id": "O1"}])  # fmt: skip
    row = layout_from_simulation(sim).features.iloc[0]
    assert (int(row["start"]), int(row["end"])) == (2000, 2232)


def test_adapter_crossing_feature_is_contract_valid(tmp_path: Path) -> None:
    from giap.genetable.contract import validate_table
    from giapsim.layout import layout_from_simulation

    sim = _fake_sim(tmp_path / "m0", [{"type": "oriC", "sequence": "chr1", "start": 5901,
                                       "stop": 132, "strand": "?", "id": "O1"}])  # fmt: skip
    lay = layout_from_simulation(sim)
    write(lay, tmp_path / "w")
    assert validate_table("features", pd.read_parquet(tmp_path / "w" / "features.parquet")) == []


def test_adapter_without_oric_has_empty_features(tmp_path: Path) -> None:
    import giapsim
    from giapsim.layout import layout_from_simulation

    sim = giapsim.simulate(giapsim.SimConfig(), 2, tmp_path / "sim")
    f = layout_from_simulation(sim).features
    assert len(f) == 0 and tuple(f.columns) == oo.FEATURE_COLUMNS
