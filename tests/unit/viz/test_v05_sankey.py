"""V5 ``mapping_sankey``: where did my sequences go? (M9-05; GIPPA_SPEC.md §10.2 V5, §9.3).

Tests the flow table (the ``Instance`` data frame = ``plots/data/V05_mapping_sankey.tsv``), the
node order, the title, the degenerate cases and determinism. No ``--mpl`` baselines (M9-25).

Interface assumed where the spec leaves it open:

- The plot lives in ``gippa.viz.mapping_sankey`` (§4: one module per plot ID); importing it
  registers ``V05`` (name ``mapping_sankey``, module ``core``). ``needs`` | ``optional`` is
  exactly {``peptide_hits``, ``protein_hits``, ``unmapped``}, and ``protein_hits`` is optional
  (a run without protein input still draws V5). It returns one ``Instance``, ``suffix=None``,
  ``coords=()``.
- Data columns are exactly ``input_class, route, outcome, n``. Outcomes are spelled ``placed``,
  ``contaminant``, ``rejected by FDR``, ``unexplained``. The table is the fixed grid of the ten
  valid flows (peptide x {exact, tolerant} x placed; protein x {miniprot, tblastn} x placed;
  {peptide, protein} x unmapped x {contaminant, rejected by FDR, unexplained}) in node order,
  zero flows included, then the ``all / all / conflicts`` row.
- Every drawn node carries a text label that starts with the node name (``exact``,
  ``exact (1)``, ``exact\\n1`` …); a node with zero flow has no label. Higher on screen = earlier
  in the node order.
- When every input table is header-only, V5 is the framework empty state: a header-only frame.
"""

from __future__ import annotations

import os
import re
import tempfile
from collections import Counter
from pathlib import Path

import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers.viz_mapping import (
    contaminant,
    mapping_out,
    peptide_hit,
    protein_hit,
    unmapped,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.patches import Patch
from matplotlib.text import Text

from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
GOLDEN = REPO / "tests" / "data" / "golden" / "sim50k" / "mapping"
VAL_OUT = Path(os.environ.get("GIPPA_ROOT", REPO / "local")) / "runs" / "val-mvp-2026-10-22" / "OUT"
COLUMNS = ["input_class", "route", "outcome", "n"]
GRID = [
    ("peptide", "exact", "placed"),
    ("peptide", "tolerant", "placed"),
    ("peptide", "unmapped", "contaminant"),
    ("peptide", "unmapped", "rejected by FDR"),
    ("peptide", "unmapped", "unexplained"),
    ("protein", "miniprot", "placed"),
    ("protein", "tblastn", "placed"),
    ("protein", "unmapped", "contaminant"),
    ("protein", "unmapped", "rejected by FDR"),
    ("protein", "unmapped", "unexplained"),
]
CONFLICTS = ("all", "all", "conflicts")
ROUTES = ("exact", "tolerant", "miniprot", "tblastn", "unmapped")
OUTCOMES = ("placed", "contaminant", "rejected by FDR", "unexplained")

#: example 1
EX1_HITS = [
    peptide_hit("PEPAK", 101),
    peptide_hit("PEPAK", 5001),
    peptide_hit("PEPBK", 9001, route="tolerant"),
]
EX1_UNMAPPED = [
    unmapped("PEPCK", "no_hit"),
    contaminant("PEPDK"),
    unmapped("PEPEK", "rejected_placement_fdr"),
]
#: example 2
EX2_PROTEINS = [protein_hit("P1", "miniprot", 1), protein_hit("P1", "tblastn", 2, start=20001)]
EX2_UNMAPPED = [unmapped("P2", "below_identity", kind="protein")]


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def _v05() -> base.PlotSpec:
    import gippa.viz.mapping_sankey  # noqa: F401  (registers V05)

    return base.registered()["V05"]


def _render(out: Path, cfg: GippaConfig) -> Instance:
    result = _v05().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _ex1(tmp: Path, **kw: object) -> Path:
    args: dict[str, object] = {"hits": EX1_HITS, "proteins": [], "unmapped_rows": EX1_UNMAPPED}
    args.update(kw)
    return mapping_out(tmp, **args)  # type: ignore[arg-type]


def _flows(inst: Instance) -> dict[tuple[str, str, str], int]:
    return {
        (str(c), str(r), str(o)): int(float(str(n)))
        for c, r, o, n in zip(
            inst.data["input_class"],
            inst.data["route"],
            inst.data["outcome"],
            inst.data["n"],
            strict=True,
        )
    }


def _nonzero(inst: Instance) -> dict[tuple[str, str, str], int]:
    return {k: v for k, v in _flows(inst).items() if v and k != CONFLICTS}


def _texts(inst: Instance) -> list[Text]:
    return [t for t in inst.figure.findobj(Text) if t.get_text()]


def _label(inst: Instance, node: str) -> list[Text]:
    pat = re.compile(rf"^{re.escape(node)}(\s|\(|$)")
    return [t for t in _texts(inst) if pat.match(t.get_text())]


def _screen_y(t: Text) -> float:
    return float(t.get_transform().transform(t.get_position())[1])


def _patch_colours(inst: Instance) -> set[str]:
    out = set()
    for ax in inst.figure.axes:
        for p in ax.findobj(Patch):
            if p is ax.patch or not p.get_visible():
                continue
            rgba = p.get_facecolor()
            if rgba[3] == 0:
                continue
            out.add(mcolors.to_hex(rgba[:3]).upper())
    return out


# --- registration and skips -------------------------------------------------------------------


def test_v05_is_registered_as_mapping_sankey() -> None:
    spec = _v05()
    assert spec.pid == "V05"
    assert spec.name == "mapping_sankey"
    assert spec.module == "core"
    assert set(spec.needs) | set(spec.optional) == {"peptide_hits", "protein_hits", "unmapped"}
    assert "protein_hits" in spec.optional
    assert "protein_hits" not in spec.needs


def test_no_mapping_tables_is_skipped_as_input_absent(tmp_path: Path, cfg: GippaConfig) -> None:
    out = mapping_out(tmp_path)
    result = _v05().fn(RunView.from_outdir(out), cfg)
    assert isinstance(result, Skip)
    assert result.reason.startswith("input_absent:")
    assert result.reason.split(":", 1)[1] in {"peptide_hits", "unmapped"}


# --- the flow table ---------------------------------------------------------------------------


def test_example1_peptide_flows_and_title(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex1(tmp_path), cfg)
    assert _nonzero(inst) == {
        ("peptide", "exact", "placed"): 1,  # PEPA counted once over its two loci
        ("peptide", "tolerant", "placed"): 1,
        ("peptide", "unmapped", "unexplained"): 1,
        ("peptide", "unmapped", "contaminant"): 1,
        ("peptide", "unmapped", "rejected by FDR"): 1,
    }
    assert _flows(inst)[CONFLICTS] == 0
    assert inst.title == "40 % of 5 input peptides placed in the genome"
    assert inst.title in [t.get_text() for t in _texts(inst)]
    assert "express" not in inst.title.lower()


def test_example2_protein_best_rank_route(tmp_path: Path, cfg: GippaConfig) -> None:
    out = mapping_out(tmp_path, hits=[], proteins=EX2_PROTEINS, unmapped_rows=EX2_UNMAPPED)
    inst = _render(out, cfg)
    assert _nonzero(inst) == {
        ("protein", "miniprot", "placed"): 1,
        ("protein", "unmapped", "unexplained"): 1,
    }
    assert _flows(inst)[("protein", "tblastn", "placed")] == 0
    assert "50 % of 2 proteins" in inst.title


def test_best_rank_wins_whatever_the_row_order(tmp_path: Path, cfg: GippaConfig) -> None:
    proteins = [protein_hit("P1", "tblastn", 2), protein_hit("P1", "miniprot", 1, start=5001)]
    proteins += [protein_hit("P3", "tblastn", 1), protein_hit("P3", "miniprot", 3, start=9001)]
    inst = _render(mapping_out(tmp_path, hits=[], proteins=proteins, unmapped_rows=[]), cfg)
    assert _nonzero(inst) == {
        ("protein", "miniprot", "placed"): 1,
        ("protein", "tblastn", "placed"): 1,
    }


def test_peptides_and_proteins_title(tmp_path: Path, cfg: GippaConfig) -> None:
    out = mapping_out(
        tmp_path,
        hits=EX1_HITS,
        proteins=EX2_PROTEINS,
        unmapped_rows=[*EX1_UNMAPPED, *EX2_UNMAPPED],
    )
    inst = _render(out, cfg)
    assert inst.title == ("40 % of 5 input peptides placed in the genome; 50 % of 2 proteins")
    assert len(_nonzero(inst)) == 7


def test_example3_placed_and_unmapped_counts_once_as_placed(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    out = _ex1(tmp_path, hits=[peptide_hit("PEPAK")], unmapped_rows=[unmapped("PEPAK", "no_hit")])
    inst = _render(out, cfg)
    assert _nonzero(inst) == {("peptide", "exact", "placed"): 1}
    assert _flows(inst)[CONFLICTS] == 1
    assert inst.title == "100 % of 1 input peptides placed in the genome"


def test_exact_wins_over_tolerant_for_one_peptide(tmp_path: Path, cfg: GippaConfig) -> None:
    hits = [
        peptide_hit("MIXEDK", 101, route="tolerant"),
        peptide_hit("MIXEDK", 901, route="exact"),
        peptide_hit("TOLK", 2001, route="tolerant"),
        peptide_hit("TOLK", 3001, route="tolerant"),
    ]
    inst = _render(_ex1(tmp_path, hits=hits, unmapped_rows=[]), cfg)
    assert _nonzero(inst) == {
        ("peptide", "exact", "placed"): 1,
        ("peptide", "tolerant", "placed"): 1,
    }


@pytest.mark.parametrize(
    ("reason", "outcome"),
    [
        ("no_hit", "unexplained"),
        ("below_identity", "unexplained"),
        ("too_short", "unexplained"),
        ("rejected_placement_fdr", "rejected by FDR"),
    ],
)
def test_unmapped_reason_to_outcome(
    tmp_path: Path, cfg: GippaConfig, reason: str, outcome: str
) -> None:
    inst = _render(_ex1(tmp_path, hits=[], unmapped_rows=[unmapped("QQQQK", reason)]), cfg)
    assert _nonzero(inst) == {("peptide", "unmapped", outcome): 1}


def test_contaminant_label_wins_over_reason(tmp_path: Path, cfg: GippaConfig) -> None:
    rows = [contaminant("CONTK"), contaminant("P9", kind="protein")]
    inst = _render(_ex1(tmp_path, hits=[], unmapped_rows=rows), cfg)
    assert _nonzero(inst) == {
        ("peptide", "unmapped", "contaminant"): 1,
        ("protein", "unmapped", "contaminant"): 1,
    }


def test_flow_table_is_the_fixed_grid_in_node_order(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex1(tmp_path), cfg)
    assert list(inst.data.columns) == COLUMNS
    keys = list(
        zip(inst.data["input_class"], inst.data["route"], inst.data["outcome"], strict=True)
    )
    assert keys == [*GRID, CONFLICTS]


def test_zero_routes_stay_in_the_tsv_and_are_not_drawn(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex1(tmp_path), cfg)
    flows = _flows(inst)
    for route in ("miniprot", "tblastn"):
        assert flows[("protein", route, "placed")] == 0
        assert _label(inst, route) == [], route
    assert _label(inst, "protein") == []  # no protein_hits rows: the class is absent
    for node in ("peptide", "exact", "tolerant", "unmapped", *OUTCOMES):
        assert len(_label(inst, node)) == 1, node


def test_no_protein_hits_file_draws_peptides_only(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex1(tmp_path, proteins=None), cfg)
    assert [k for k in _flows(inst) if k[0] == "protein"] == [k for k in GRID if k[0] == "protein"]
    assert all(v == 0 for k, v in _flows(inst).items() if k[0] == "protein")
    assert inst.title == "40 % of 5 input peptides placed in the genome"


def test_protein_only_run_has_only_the_protein_class(tmp_path: Path, cfg: GippaConfig) -> None:
    out = mapping_out(tmp_path, hits=[], proteins=EX2_PROTEINS, unmapped_rows=EX2_UNMAPPED)
    inst = _render(out, cfg)
    assert all(v == 0 for k, v in _flows(inst).items() if k[0] == "peptide")
    assert _label(inst, "peptide") == []
    assert len(_label(inst, "protein")) == 1
    assert "NA" not in inst.title


def test_nodes_are_drawn_in_the_fixed_order(tmp_path: Path, cfg: GippaConfig) -> None:
    out = mapping_out(
        tmp_path,
        hits=EX1_HITS,
        proteins=[protein_hit("P1", "miniprot", 1), protein_hit("P4", "tblastn", 1)],
        unmapped_rows=[*EX1_UNMAPPED, *EX2_UNMAPPED],
    )
    inst = _render(out, cfg)
    for column in (("peptide", "protein"), ROUTES, OUTCOMES):
        labels = [_label(inst, node) for node in column]
        assert all(len(found) == 1 for found in labels), column
        ys = [_screen_y(found[0]) for found in labels]
        assert ys == sorted(ys, reverse=True), (column, ys)


def test_placed_is_blue_and_outcomes_have_distinct_colours(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    inst = _render(_ex1(tmp_path), cfg)
    colours = _patch_colours(inst)
    assert "#0072B2" in colours
    okabe = {c.upper() for c in style.OKABE_ITO.values()} - {"#0072B2"}
    assert len(colours & okabe) >= 3  # contaminant, rejected by FDR, unexplained


def test_every_text_is_at_least_6_pt(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex1(tmp_path), cfg)
    sizes = [t.get_fontsize() for t in _texts(inst)]
    assert sizes
    assert min(sizes) >= 6


def test_width_is_a_style_width(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_ex1(tmp_path), cfg)
    assert inst.suffix is None
    assert inst.coords == ()
    width = round(float(inst.figure.get_size_inches()[0]), 6)
    assert width in (style.SINGLE_IN, style.DOUBLE_IN)


def test_every_table_header_only_is_the_empty_state(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(mapping_out(tmp_path, hits=[], proteins=[], unmapped_rows=[]), cfg)
    assert list(inst.data.columns) == COLUMNS
    assert len(inst.data) == 0
    assert inst.title
    assert "%" not in inst.title


def test_many_sequences_draw_at_most_the_grid(tmp_path: Path, cfg: GippaConfig) -> None:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    def pep(i: int) -> str:
        chars = []
        for _ in range(4):
            chars.append(alphabet[i % 20])
            i //= 20
        return "".join(chars) + "GGGK"

    hits = [peptide_hit(pep(i), 1 + 30 * i) for i in range(20_000)]
    rows = [unmapped(pep(i), "no_hit") for i in range(20_000, 20_500)]
    inst = _render(_ex1(tmp_path, hits=hits, unmapped_rows=rows), cfg)
    assert _nonzero(inst) == {
        ("peptide", "exact", "placed"): 20_000,
        ("peptide", "unmapped", "unexplained"): 500,
    }
    n_patches = sum(len(ax.patches) for ax in inst.figure.axes)
    assert n_patches <= 2 * 5 + 5 + 4 + 2 * 5 * 4  # nodes + ribbons, never per sequence


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    placed=st.lists(st.sampled_from(["exact", "tolerant", "both"]), max_size=8),
    lost=st.lists(
        st.sampled_from(["no_hit", "too_short", "rejected_placement_fdr", "contaminant"]),
        max_size=8,
    ),
    n_conflicts=st.integers(min_value=0, max_value=3),
)
def test_property_flows_partition_distinct_sequences(
    placed: list[str], lost: list[str], n_conflicts: int
) -> None:
    style.apply()
    hits = []
    for i, kind in enumerate(placed):
        name = f"PLC{i:02d}K"
        if kind in ("exact", "both"):
            hits.append(peptide_hit(name, 101 + 300 * i))
        if kind in ("tolerant", "both"):
            hits.append(peptide_hit(name, 50_101 + 300 * i, route="tolerant"))
    rows = [
        contaminant(f"LST{i:02d}K") if r == "contaminant" else unmapped(f"LST{i:02d}K", r)
        for i, r in enumerate(lost)
    ]
    k = min(n_conflicts, len(placed))
    rows += [unmapped(f"PLC{i:02d}K", "no_hit") for i in range(k)]
    if not placed and not rows:
        return
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_ex1(Path(tmp), hits=hits, unmapped_rows=rows), resolve_config())
    flows = _flows(inst)
    assert flows[CONFLICTS] == k
    total = len(placed) + len(lost)
    assert sum(v for key, v in flows.items() if key != CONFLICTS) == total
    n_placed = flows[("peptide", "exact", "placed")] + flows[("peptide", "tolerant", "placed")]
    assert n_placed == len(placed)
    assert flows[("peptide", "tolerant", "placed")] == placed.count("tolerant")
    lost_by = Counter(lost)
    unexplained = lost_by["no_hit"] + lost_by["too_short"]
    assert flows[("peptide", "unmapped", "unexplained")] == unexplained
    assert flows[("peptide", "unmapped", "rejected by FDR")] == lost_by["rejected_placement_fdr"]
    assert flows[("peptide", "unmapped", "contaminant")] == lost_by["contaminant"]
    assert inst.title.startswith(f"{base.title_pct(len(placed), total)} of {total} input peptides")


# --- the golden M3 tables, files and determinism ----------------------------------------------


def _golden_out(tmp: Path) -> Path:
    out = mapping_out(tmp)
    (out / "mapping").mkdir(parents=True, exist_ok=True)
    for name in ("peptide_hits", "protein_hits", "unmapped", "mapping_fdr"):
        (out / "mapping" / f"{name}.tsv").write_bytes((GOLDEN / f"{name}.tsv").read_bytes())
    return out


def test_renders_on_the_sim50k_golden_tables(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_golden_out(tmp_path), cfg)
    assert _nonzero(inst) == {
        ("peptide", "exact", "placed"): 905,
        ("peptide", "unmapped", "unexplained"): 9,  # 5 no_hit + 4 too_short
        ("peptide", "unmapped", "rejected by FDR"): 4,
    }
    assert _flows(inst)[CONFLICTS] == 0
    assert inst.title == "99 % of 918 input peptides placed in the genome"


def test_data_tsv_and_images_are_byte_identical(tmp_path: Path, cfg: GippaConfig) -> None:
    out = _ex1(tmp_path / "run")
    written = []
    for name in ("x", "y"):
        inst = _render(out, cfg)
        paths = write_instance(
            inst, "V05", "mapping_sankey", tmp_path / name, ["svg", "png", "pdf"], 100
        )
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert set(written[0]) == {
        Path("V05_mapping_sankey.svg"),
        Path("V05_mapping_sankey.png"),
        Path("V05_mapping_sankey.pdf"),
        Path("data/V05_mapping_sankey.tsv"),
    }
    assert written[0] == written[1]
    tsv = written[0][Path("data/V05_mapping_sankey.tsv")].decode("utf-8").splitlines()
    assert tsv[0] == "\t".join(COLUMNS)
    assert tsv[1] == "peptide\texact\tplaced\t1"
    assert "peptide\tunmapped\trejected by FDR\t1" in tsv
    assert "protein\tminiprot\tplaced\t0" in tsv
    assert tsv[-1] == "all\tall\tconflicts\t0"


def test_render_does_not_mutate_between_calls(tmp_path: Path, cfg: GippaConfig) -> None:
    view = RunView.from_outdir(_ex1(tmp_path))
    first = _v05().fn(view, cfg)
    second = _v05().fn(view, cfg)
    assert isinstance(first, list) and isinstance(second, list)
    pd.testing.assert_frame_equal(first[0].data, second[0].data)
    assert first[0].title == second[0].title


# --- the validation run (local/, slow): V5 and V6 together -------------------------------------


@pytest.mark.slow
@pytest.mark.skipif(not VAL_OUT.is_dir(), reason=f"validation run absent: {VAL_OUT}")
def test_v05_and_v06_render_on_the_validation_run(tmp_path: Path, cfg: GippaConfig) -> None:
    """1644 exact peptides, 24 unmapped (13 no_hit, 4 too_short, 7 rejected), no protein hits;
    ``mapping_fdr`` all/1 (reverse) = 7 targets, fdr 0.142857; strata 0 and 2 NA."""
    import gippa.viz.mapping_fdr  # noqa: F401  (registers V06)

    view = RunView.from_outdir(VAL_OUT)
    inst = _render(VAL_OUT, cfg)
    assert _nonzero(inst) == {
        ("peptide", "exact", "placed"): 1644,
        ("peptide", "unmapped", "unexplained"): 17,
        ("peptide", "unmapped", "rejected by FDR"): 7,
    }
    assert _flows(inst)[CONFLICTS] == 0
    assert re.fullmatch(r"99 % of 1,?668 input peptides placed in the genome", inst.title)
    v06 = base.registered()["V06"].fn(view, cfg)
    assert isinstance(v06, list) and len(v06) == 1
    b = v06[0].data[v06[0].data["panel"] == "b"]
    fdr = b[(b["metric"] == "fdr_stratum")]
    assert len(fdr) == 3
    assert "unestimated" in v06[0].title
    for pid, name, one in (("V05", "mapping_sankey", inst), ("V06", "mapping_fdr", v06[0])):
        paths = write_instance(one, pid, name, tmp_path / "plots", ["svg", "png"], 100)
        assert all(p.stat().st_size > 0 for p in paths)


def test_protein_rows_without_a_gene_are_not_placed(tmp_path: Path, cfg: GippaConfig) -> None:
    """Review N1 (#198): M3 keeps miniprot rows with no gene; tblastn then assigns the protein at
    its own rank 1. V5 credits the gene-assigned route; a protein with no gene row is not placed."""
    proteins = [
        protein_hit("P1", "miniprot", 1, gene=None),
        protein_hit("P1", "tblastn", 1, gene="g1"),
        protein_hit("P2", "miniprot", 1, gene=None),
    ]
    rows = [unmapped("P2", "below_identity", kind="protein")]
    inst = _render(mapping_out(tmp_path, hits=[], proteins=proteins, unmapped_rows=rows), cfg)
    flows = _flows(inst)
    assert flows[("protein", "tblastn", "placed")] == 1
    assert flows[("protein", "miniprot", "placed")] == 0
    assert flows[("protein", "unmapped", "unexplained")] == 1


@pytest.mark.parametrize(
    ("n", "d", "want"), [(1667, 1668, "99 %"), (1, 300, "1 %"), (0, 5, "0 %"), (5, 5, "100 %")]
)
def test_title_pct_never_rounds_a_partial_share_to_0_or_100(n: int, d: int, want: str) -> None:
    """Review N3 (#198): 1667 of 1668 placed is not "100 % placed"."""
    assert base.title_pct(n, d) == want
