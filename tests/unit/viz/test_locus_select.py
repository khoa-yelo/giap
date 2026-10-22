"""V3 locus pool, order and cap (M9-10 "Selection and order"; GIAP_SPEC.md §10.2 V3, §10.1;
M9-index decision 18; lead decisions 2026-10-22).

Interface assumed where the spec leaves it open (reported to the lead):

- ``giap.viz.locus`` exports ``Locus``, ``select_loci``, ``render_locus``; importing it registers
  ``V03`` with the name ``locus``. ``needs`` holds none of ``regions_stats``, ``islands``,
  ``operon_flags``, ``peptide_hits`` (each may be absent; the degenerate table).
- ``select_loci`` returns the **whole** ordered pool (the index lists every candidate); the plot
  renders its first ``report.max_locus_plots`` loci.
- ``Locus.start0`` / ``end0`` are the locus *span* (0-based half-open; an island crossing the
  origin keeps ``end0 > length``); ``error`` is ``None`` for NA.
- The plot returns one ``Instance`` per rendered locus (``suffix = locus_id``) plus one index
  ``Instance`` with ``suffix = None`` (written as ``V03_locus.*`` and ``data/V03_locus.tsv``):
  data columns exactly ``rank, locus_id, kind, contig, start, end, error, rendered``, ``start`` in
  ``coords`` (so the TSV is 1-based closed), ``rank`` 1-based. With no loci the index is the
  framework empty state titled ``"no loci to show"``.
- A merged operon locus (flags in one TU) keeps the ``pair_id`` of its best (first-ordered) flag
  and that flag's q as ``error``.
- Numbers in labels print in ``%g`` form (``0.003``, ``0.01``, ``21.5``), as in V11.
"""

from __future__ import annotations

import math
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import viz_locus as vl
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip
from giap.viz.locus import Locus, select_loci

INDEX_COLUMNS = ["rank", "locus_id", "kind", "contig", "start", "end", "error", "rendered"]
KIND_ORDER = ("prophage", "bgc", "island", "operon_flag")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


def _cfg(*overrides: str) -> GiapConfig:
    return resolve_config(overrides=list(overrides))


def _v03() -> base.PlotSpec:
    import giap.viz.locus  # noqa: F401  (registers V03)

    return base.registered()["V03"]


def _run(out: Path, cfg: GiapConfig) -> tuple[Instance, list[Instance]]:
    result = _v03().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    index = [i for i in result if i.suffix is None]
    loci = [i for i in result if i.suffix is not None]
    assert len(index) == 1, [i.suffix for i in result]
    return index[0], loci


def _ids(loci: Sequence[Locus]) -> list[str]:
    return [lo.locus_id for lo in loci]


# --- a TU pair for operon flags -------------------------------------------------------------


def _operon_genes() -> tuple[list[dict[str, object]], list[tuple]]:
    """g3, g4 in TU0001 (chr1 1-based 40,001-42,000); g5 alone in TU0002 (44,001-45,000)."""
    genes = [
        vl.gene("g3", "chr1", 40_001, 41_000, tu_id="TU0001", log2=20.0),
        vl.gene("g4", "chr1", 41_101, 42_000, tu_id="TU0001", log2=22.0),
        vl.gene("g5", "chr1", 44_001, 45_000, tu_id="TU0002", log2=19.0),
    ]
    tus = [
        vl.tu("TU0001", "chr1", 40_001, 42_000, "+", ["g3", "g4"]),
        vl.tu("TU0002", "chr1", 44_001, 45_000, "+", ["g5"]),
    ]
    return genes, tus


def _example1_out(tmp: Path) -> Path:
    genes, tus = _operon_genes()
    return vl.build_out(
        tmp,
        genes=genes,
        tus=tus,
        regions_stats=[
            vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 35_000, 0.003),
            vl.region_stat("prophage_2", "prophage", "chr1", 60_001, 62_000, None),
        ],
        islands=[vl.island("ISL_1", "chr1", 70_001, 75_000, 0.01)],
        flags=[vl.flag("PAIR_7", "g3", "g4", 0.02)],
    )


# --- registration and the pool --------------------------------------------------------------


def test_v03_is_registered_as_locus_without_hard_table_needs() -> None:
    spec = _v03()
    assert spec.pid == "V03"
    assert spec.name == "locus"
    for key in ("regions_stats", "islands", "operon_flags", "peptide_hits", "novel_orfs"):
        assert key not in spec.needs


def test_locus_is_a_frozen_dataclass(tmp_path: Path, cfg: GiapConfig) -> None:
    lo = select_loci(RunView.from_outdir(_example1_out(tmp_path)), cfg)[0]
    assert isinstance(lo, Locus)
    with pytest.raises(AttributeError):
        lo.locus_id = "x"  # type: ignore[misc]


def test_example1_order_across_kinds(tmp_path: Path, cfg: GiapConfig) -> None:
    loci = select_loci(RunView.from_outdir(_example1_out(tmp_path)), cfg)
    assert _ids(loci) == ["prophage_1", "ISL_1", "PAIR_7", "prophage_2"]
    assert [lo.kind for lo in loci] == ["prophage", "island", "operon_flag", "prophage"]
    assert [lo.error for lo in loci[:3]] == pytest.approx([0.003, 0.01, 0.02])
    assert loci[3].error is None


def test_spans_are_zero_based_half_open(tmp_path: Path, cfg: GiapConfig) -> None:
    loci = {
        lo.locus_id: lo for lo in select_loci(RunView.from_outdir(_example1_out(tmp_path)), cfg)
    }
    assert (loci["prophage_1"].contig, loci["prophage_1"].start0, loci["prophage_1"].end0) == (
        "chr1",
        10_000,
        35_000,
    )
    assert (loci["ISL_1"].start0, loci["ISL_1"].end0) == (70_000, 75_000)
    # the flag's span is the TU containing gene_a
    assert (loci["PAIR_7"].start0, loci["PAIR_7"].end0) == (40_000, 42_000)


def test_island_crossing_the_origin_keeps_its_unwrapped_span(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    out = vl.build_out(tmp_path, islands=[vl.island("ISL_1", "chr1", 99_001, 101_500, 0.01)])
    (lo,) = select_loci(RunView.from_outdir(out), cfg)
    assert (lo.kind, lo.start0, lo.end0) == ("island", 99_000, 101_500)


def test_pool_takes_only_all_subset_o_minus_e_prophage_and_bgc_rows(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    rows = [
        vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 20_000, 0.04),
        vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 20_000, 0.0001, stat="hl_shift"),
        vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 20_000, None, subset="integrase"),
        vl.region_stat("bgc_1", "bgc", "chr1", 30_001, 40_000, 0.2),
        vl.region_stat("plasmid_1", "plasmid", "chr1", 1, 100_000, 0.0001),
        vl.region_stat("crispr_1", "crispr", "chr1", 50_001, 50_500, 0.0001),
        vl.region_stat("amr_1", "amr", "chr1", 60_001, 61_000, 0.0001),
    ]
    out = vl.build_out(tmp_path, regions_stats=rows)
    loci = select_loci(RunView.from_outdir(out), cfg)
    assert _ids(loci) == ["prophage_1", "bgc_1"]
    assert loci[0].error == pytest.approx(0.04)
    assert [lo.kind for lo in loci] == ["prophage", "bgc"]


def test_regions_are_not_filtered_by_alpha(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(
        tmp_path,
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 20_000, 0.9)],
        islands=[vl.island("ISL_1", "chr1", 50_001, 52_000, 1.0)],
    )
    assert _ids(select_loci(RunView.from_outdir(out), cfg)) == ["prophage_1", "ISL_1"]


def test_novel_orfs_are_not_in_the_pool(tmp_path: Path, cfg: GiapConfig) -> None:
    orf = {
        "orf_id": "NORF_1",
        "class": "intergenic",
        "contig": "chr1",
        "strand": "+",
        "frame": 1,
        "start_longest": 5_001,
        "start_best_ostir": 5_001,
        "end": 5_300,
        "called": True,
        "rank": 1,
        "class_fdr": "0.001",
    }
    out = vl.build_out(
        tmp_path,
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 20_000, 0.04)],
        novel_orfs=[orf],
    )
    loci = select_loci(RunView.from_outdir(out), cfg)
    assert _ids(loci) == ["prophage_1"]
    assert all(lo.kind != "novel_orf" for lo in loci)


# --- operon flags ------------------------------------------------------------------------------


def test_flags_are_kept_at_q_at_most_alpha(tmp_path: Path) -> None:
    genes, tus = _operon_genes()
    genes += [
        vl.gene("g6", "chr1", 50_001, 51_000, tu_id="TU0003"),
        vl.gene("g7", "chr1", 51_101, 52_000, tu_id="TU0003"),
        vl.gene("g8", "chr1", 60_001, 61_000, tu_id="TU0004"),
        vl.gene("g9", "chr1", 61_101, 62_000, tu_id="TU0004"),
    ]
    tus += [
        vl.tu("TU0003", "chr1", 50_001, 52_000, "+", ["g6", "g7"]),
        vl.tu("TU0004", "chr1", 60_001, 62_000, "+", ["g8", "g9"]),
    ]
    flags = [
        vl.flag("PAIR_1", "g3", "g4", 0.10),  # = alpha (0.10): kept
        vl.flag("PAIR_2", "g6", "g7", 0.11),  # > alpha: dropped
        vl.flag("PAIR_3", "g8", "g9", None),  # NA: dropped
    ]
    out = vl.build_out(tmp_path, genes=genes, tus=tus, flags=flags)
    view = RunView.from_outdir(out)
    assert _ids(select_loci(view, _cfg())) == ["PAIR_1"]
    assert _ids(select_loci(view, _cfg("stats.alpha.operon_flags=0.05"))) == []


def test_flag_genes_in_different_tus_span_both_genes(tmp_path: Path, cfg: GiapConfig) -> None:
    genes, tus = _operon_genes()
    out = vl.build_out(tmp_path, genes=genes, tus=tus, flags=[vl.flag("PAIR_2", "g4", "g5", 0.01)])
    (lo,) = select_loci(RunView.from_outdir(out), cfg)
    assert (lo.locus_id, lo.kind, lo.contig) == ("PAIR_2", "operon_flag", "chr1")
    assert (lo.start0, lo.end0) == (41_100, 45_000)  # g4 start .. g5 end


def test_example7_flags_in_one_tu_merge_into_one_locus(tmp_path: Path, cfg: GiapConfig) -> None:
    genes, tus = _operon_genes()
    flags = [
        vl.flag("PAIR_8", "g4", "g3", 0.05, flag="unexplained_discordance"),
        vl.flag("PAIR_7", "g3", "g4", 0.02),
    ]
    out = vl.build_out(tmp_path, genes=genes, tus=tus, flags=flags)
    loci = select_loci(RunView.from_outdir(out), cfg)
    assert _ids(loci) == ["PAIR_7"]
    (lo,) = loci
    assert lo.error == pytest.approx(0.02)
    assert (lo.start0, lo.end0) == (40_000, 42_000)
    assert "internal_promoter_supported between g3 and g4 (q = 0.02)" in lo.label
    assert "unexplained_discordance between g4 and g3 (q = 0.05)" in lo.label
    _, rendered = _run(out, cfg)
    assert [i.suffix for i in rendered] == ["PAIR_7"]
    assert rendered[0].title == lo.label


# --- order: ties and NA ------------------------------------------------------------------------


def test_ties_break_by_kind_then_id(tmp_path: Path, cfg: GiapConfig) -> None:
    genes, tus = _operon_genes()
    out = vl.build_out(
        tmp_path,
        genes=genes,
        tus=tus,
        regions_stats=[
            vl.region_stat("bgc_1", "bgc", "chr1", 20_001, 22_000, 0.05),
            vl.region_stat("prophage_2", "prophage", "chr1", 12_001, 14_000, 0.05),
            vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 11_000, 0.05),
        ],
        islands=[vl.island("ISL_1", "chr1", 70_001, 75_000, 0.05)],
        flags=[vl.flag("PAIR_1", "g3", "g4", 0.05)],
    )
    loci = select_loci(RunView.from_outdir(out), cfg)
    assert _ids(loci) == ["prophage_1", "prophage_2", "bgc_1", "ISL_1", "PAIR_1"]


def test_na_errors_sort_last_then_by_kind(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(
        tmp_path,
        regions_stats=[
            vl.region_stat("prophage_2", "prophage", "chr1", 12_001, 14_000, None),
            vl.region_stat("bgc_1", "bgc", "chr1", 20_001, 22_000, 0.9),
        ],
        islands=[vl.island("ISL_2", "chr1", 70_001, 75_000, None)],
    )
    loci = select_loci(RunView.from_outdir(out), cfg)
    assert _ids(loci) == ["bgc_1", "prophage_2", "ISL_2"]
    assert [lo.error is None for lo in loci] == [False, True, True]


@st.composite
def _pools(draw: Any) -> list[tuple[str, float | None]]:
    errors = st.sampled_from([None, 0.001, 0.01, 0.05, 0.5, 1.0])
    pool: list[tuple[str, float | None]] = []
    for kind in ("prophage", "bgc", "island"):
        pool += [(kind, draw(errors)) for _ in range(draw(st.integers(0, 3)))]
    flag_q = st.sampled_from([0.001, 0.01, 0.05, 0.1])
    pool += [("operon_flag", draw(flag_q)) for _ in range(draw(st.integers(0, 3)))]
    return pool


@settings(max_examples=30, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(_pools())
def test_property_order_is_error_na_last_then_kind_then_id(pool: list[tuple[str, Any]]) -> None:
    regions, islands, flags, genes, tus = [], [], [], [], []
    expected: list[tuple[bool, float, int, str]] = []
    counters = dict.fromkeys(KIND_ORDER, 0)
    for kind, err in pool:
        counters[kind] += 1
        n = counters[kind]
        base_pos = 7_000 * len(expected) + 1
        if kind in ("prophage", "bgc"):
            lid = f"{kind}_{n}"
            regions.append(vl.region_stat(lid, kind, "chr1", base_pos, base_pos + 999, err))
        elif kind == "island":
            lid = f"ISL_{n}"
            islands.append(vl.island(lid, "chr1", base_pos, base_pos + 999, err))
        else:
            lid = f"PAIR_{n}"
            a, b, tid = f"ga{n}", f"gb{n}", f"TU{n:04d}"
            genes += [
                vl.gene(a, "chr1", base_pos, base_pos + 299, tu_id=tid),
                vl.gene(b, "chr1", base_pos + 400, base_pos + 699, tu_id=tid),
            ]
            tus.append(vl.tu(tid, "chr1", base_pos, base_pos + 699, "+", [a, b]))
            flags.append(vl.flag(lid, a, b, err))
        expected.append((err is None, 0.0 if err is None else err, KIND_ORDER.index(kind), lid))
    with tempfile.TemporaryDirectory() as tmp:
        out = vl.build_out(
            Path(tmp), genes=genes, tus=tus, regions_stats=regions, islands=islands, flags=flags
        )
        loci = select_loci(RunView.from_outdir(out), resolve_config())
    assert _ids(loci) == [e[3] for e in sorted(expected)]


# --- the cap and the index ---------------------------------------------------------------------


def _many_out(tmp: Path, n: int) -> Path:
    """``n`` prophages on chr1, 1 kb each, q = (i + 1) / 1000 so the order is prophage_1 .. n."""
    rows = [
        vl.region_stat(f"prophage_{i + 1}", "prophage", "chr1", i * 1_500 + 1, i * 1_500 + 1_000,
                       (i + 1) / 1000)
        for i in range(n)
    ]  # fmt: skip
    return vl.build_out(tmp, regions_stats=rows)


def test_example2_sixty_candidates_fifty_rendered(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _many_out(tmp_path, 60)
    assert cfg.report.max_locus_plots == 50
    assert len(select_loci(RunView.from_outdir(out), cfg)) == 60
    index, rendered = _run(out, cfg)
    assert [i.suffix for i in rendered] == [f"prophage_{i}" for i in range(1, 51)]
    assert list(index.data.columns) == INDEX_COLUMNS
    assert len(index.data) == 60
    assert list(index.data["rank"]) == list(range(1, 61))
    assert list(index.data["locus_id"]) == [f"prophage_{i}" for i in range(1, 61)]
    assert [bool(v) for v in index.data["rendered"]] == [True] * 50 + [False] * 10


def test_the_cap_reads_report_max_locus_plots(tmp_path: Path) -> None:
    out = _many_out(tmp_path, 4)
    index, rendered = _run(out, _cfg("report.max_locus_plots=2"))
    assert [i.suffix for i in rendered] == ["prophage_1", "prophage_2"]
    assert [bool(v) for v in index.data["rendered"]] == [True, True, False, False]


def test_max_locus_plots_zero_renders_none_and_lists_all(tmp_path: Path) -> None:
    out = _many_out(tmp_path, 3)
    index, rendered = _run(out, _cfg("report.max_locus_plots=0"))
    assert rendered == []
    assert list(index.data["locus_id"]) == ["prophage_1", "prophage_2", "prophage_3"]
    assert [bool(v) for v in index.data["rendered"]] == [False, False, False]


def test_no_loci_writes_a_header_only_index(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(tmp_path, genes=[vl.gene("g1", "chr1", 101, 400, log2=20.0)])
    index, rendered = _run(out, cfg)
    assert rendered == []
    assert list(index.data.columns) == INDEX_COLUMNS
    assert len(index.data) == 0
    assert index.title == "no loci to show"


def test_header_only_tables_are_no_loci(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(tmp_path, regions_stats=[], islands=[], flags=[], tus=[])
    index, rendered = _run(out, cfg)
    assert rendered == []
    assert len(index.data) == 0


def test_example1_index_tsv_is_one_based_with_kind_and_na_error(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    out = _example1_out(tmp_path / "run")
    index, _ = _run(out, cfg)
    paths = write_instance(index, "V03", "locus", tmp_path / "plots", ["svg"], 100)
    tsv = tmp_path / "plots" / "data" / "V03_locus.tsv"
    assert tsv in paths
    lines = tsv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "\t".join(INDEX_COLUMNS)
    assert lines[1:] == [
        "1\tprophage_1\tprophage\tchr1\t10001\t35000\t0.003\ttrue",
        "2\tISL_1\tisland\tchr1\t70001\t75000\t0.01\ttrue",
        "3\tPAIR_7\toperon_flag\tchr1\t40001\t42000\t0.02\ttrue",
        "4\tprophage_2\tprophage\tchr1\t60001\t62000\tNA\ttrue",
    ]


def test_index_rows_match_select_loci(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _example1_out(tmp_path)
    loci = select_loci(RunView.from_outdir(out), cfg)
    index, _ = _run(out, cfg)
    frame: pd.DataFrame = index.data
    assert list(frame["locus_id"]) == _ids(loci)
    assert list(frame["kind"]) == [lo.kind for lo in loci]
    assert [int(v) for v in frame["start"]] == [lo.start0 for lo in loci]
    assert [int(v) for v in frame["end"]] == [lo.end0 for lo in loci]
    errs = [None if pd.isna(v) else float(v) for v in frame["error"]]
    assert errs == [lo.error for lo in loci]


# --- labels and file names -----------------------------------------------------------------------


def test_labels_follow_the_title_templates(tmp_path: Path, cfg: GiapConfig) -> None:
    loci = {
        lo.locus_id: lo for lo in select_loci(RunView.from_outdir(_example1_out(tmp_path)), cfg)
    }
    assert loci["prophage_1"].label == "prophage_1: 12 % detected vs 31 % expected (q = 0.003)"
    assert loci["prophage_2"].label == "prophage_2: 12 % detected vs 31 % expected (q = NA)"
    assert loci["ISL_1"].label == "ISL_1 (high): mean log2 21.5, FWER p 0.01"
    assert loci["PAIR_7"].label == "internal_promoter_supported between g3 and g4 (q = 0.02)"


def test_low_island_label(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(
        tmp_path, islands=[vl.island("ISL_3", "chr1", 5_001, 9_000, 0.002, direction="low",
                                     mean=17.25)]
    )  # fmt: skip
    (lo,) = select_loci(RunView.from_outdir(out), cfg)
    assert lo.label == "ISL_3 (low): mean log2 17.25, FWER p 0.002"


def test_instance_files_are_named_by_sanitised_locus_id(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vl.build_out(
        tmp_path / "run",
        regions_stats=[vl.region_stat("bgc/1 x", "bgc", "chr1", 20_001, 22_000, 0.01)],
    )
    _, rendered = _run(out, cfg)
    (inst,) = rendered
    assert inst.suffix == "bgc/1 x"
    paths = write_instance(inst, "V03", "locus", tmp_path / "plots", ["svg", "png", "pdf"], 100)
    names = sorted(p.relative_to(tmp_path / "plots").as_posix() for p in paths)
    assert names == [
        "V03_locus_bgc_1_x.pdf",
        "V03_locus_bgc_1_x.png",
        "V03_locus_bgc_1_x.svg",
        "data/V03_locus_bgc_1_x.tsv",
    ]


def test_error_is_finite_or_none(tmp_path: Path, cfg: GiapConfig) -> None:
    for lo in select_loci(RunView.from_outdir(_example1_out(tmp_path)), cfg):
        assert lo.error is None or math.isfinite(lo.error)
