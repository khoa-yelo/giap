"""``giap.viz.support.variant_windows``: the reference protein around each top ALT site (M9-18
"Support: ``plots/support/variant_windows.tsv``", example 3, the ``-`` strand, origin-crossing
and near-start rows of the degenerate-case table; lead decision 2: ±10 aa).

A separate file from tests/unit/viz/test_support.py (M9-03's tests are not edited here; the spec's
"test_support.py extended" lands as this file).

Interface assumed beyond the spec text (M9-18 test author):

- ``variant_windows(variants, genes) -> DataFrame``. ``variants`` is the ``variants.tsv`` frame as
  ``RunView.table("variants")`` reads it (schema types; ``class_error_rate`` a string);
  ``genes`` is a C2 ``genes.parquet`` frame (``locus_tag``, ``aa``, ...). Pure.
- The result is the table **as written** (M9-03 convention): columns exactly ``variant_id, gene,
  aa_pos, window_start_aa, ref_window``, ``aa_pos`` / ``window_start_aa`` 1-based integers,
  ``ref_window`` the protein residues ``[max(1, aa_pos - 10), min(len, aa_pos + 10)]``.
- One row per panel (b) site: the at most 30 ``call = ALT`` sites in the V12 top-site order
  (confidence, ``n_alt`` desc, class error asc NA last, ``variant_id``); ``mixed`` sites get none.
  Row order is not pinned; ``variant_id`` is unique.
"""

from __future__ import annotations

import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import viz_variants as vv
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.io.run_view import RunView

COLUMNS = ["variant_id", "gene", "aa_pos", "window_start_aa", "ref_window"]
SETTINGS = settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)


def _fn() -> Any:
    from giap.viz.support import variant_windows

    return variant_windows


def _variants(rows: Sequence[Mapping[str, object]]) -> pd.DataFrame:
    """``rows`` written as ``variants.tsv`` and read back through ``RunView``."""
    with tempfile.TemporaryDirectory() as tmp:
        out = vv.build_out(Path(tmp), rows)
        frame = RunView.from_outdir(out).table("variants")
    assert frame is not None
    return frame


def _by_id(df: pd.DataFrame) -> dict[str, pd.Series]:
    assert df["variant_id"].is_unique, df
    return {str(r["variant_id"]): r for _, r in df.iterrows()}


def _windows(
    rows: Sequence[Mapping[str, object]], proteins: Mapping[str, str] = vv.PROTEINS
) -> pd.DataFrame:
    return _fn()(_variants(rows), vv.genes_parquet_frame(proteins))


# ---------------------------------------------------------------- examples


def test_example_3_window_204_to_224() -> None:
    out = _windows([vv.variant("VAR0001")])
    assert list(out.columns) == COLUMNS
    r = _by_id(out)["VAR0001"]
    assert str(r["gene"]) == "GSIM_00105"
    assert int(r["aa_pos"]) == 214
    assert int(r["window_start_aa"]) == 204
    assert str(r["ref_window"]) == vv.P105[203:224]
    assert len(str(r["ref_window"])) == 21
    assert str(r["ref_window"])[10] == "Q"


def test_near_the_protein_start_is_clipped_at_1() -> None:
    r = _by_id(_windows([vv.variant("VAR0001", aa_pos=3, ref_aa=vv.P105[2])]))["VAR0001"]
    assert int(r["window_start_aa"]) == 1
    assert str(r["ref_window"]) == vv.P105[0:13]


def test_near_the_protein_end_is_clipped_at_the_last_residue() -> None:
    r = _by_id(_windows([vv.variant("VAR0001", aa_pos=298, ref_aa=vv.P105[297])]))["VAR0001"]
    assert int(r["window_start_aa"]) == 288
    assert str(r["ref_window"]) == vv.P105[287:300]


def test_minus_strand_window_is_the_protein_in_transcription_order() -> None:
    site = vv.variant("VAR0002", gene="GSIM_00200", aa_pos=50, ref_aa="D", alt_aa="N", strand="-")
    r = _by_id(_windows([site]))["VAR0002"]
    assert int(r["window_start_aa"]) == 40
    assert str(r["ref_window"]) == vv.P200[39:60]  # not reversed
    assert str(r["ref_window"])[10] == "D"


def test_origin_crossing_gene_uses_protein_positions() -> None:
    site = vv.variant("VAR0004", gene="GSIM_00900", aa_pos=35, ref_aa="S", alt_aa="T")
    r = _by_id(_windows([site]))["VAR0004"]
    assert int(r["window_start_aa"]) == 25
    assert str(r["ref_window"]) == vv.P900[24:45]


def test_short_protein_window_is_the_whole_protein() -> None:
    prot = {"GSHORT": "MKTAYIAK"}
    site = vv.variant("VAR0001", gene="GSHORT", aa_pos=4, ref_aa="A", contig="chr2",
                      codon_start=10)  # fmt: skip
    r = _by_id(_windows([site], prot))["VAR0001"]
    assert int(r["window_start_aa"]) == 1
    assert str(r["ref_window"]) == "MKTAYIAK"


def test_main_fixture_rows_match_the_oracle() -> None:
    sites = vv.main_variants()
    out = _windows(sites)
    alt = [s for s in sites if s["call"] == "ALT"]
    want = vv.windows_frame(alt)
    got = out.sort_values("variant_id").reset_index(drop=True)
    assert list(got["variant_id"]) == list(want["variant_id"])
    assert [int(v) for v in got["aa_pos"]] == list(want["aa_pos"])
    assert [int(v) for v in got["window_start_aa"]] == list(want["window_start_aa"])
    assert [str(v) for v in got["ref_window"]] == list(want["ref_window"])


def test_mixed_sites_get_no_window() -> None:
    sites = vv.main_variants()
    assert "VAR0005" not in set(_windows(sites)["variant_id"])


def test_only_the_top_30_alt_sites() -> None:
    sites = [
        vv.variant(f"VAR{i:04d}", aa_pos=5 + 3 * i, n_alt=100 - i, ref_aa=vv.P105[4 + 3 * i])
        for i in range(1, 36)
    ]
    sites.append(vv.variant("VAR0099", aa_pos=280, call="mixed", n_alt=500))
    out = _windows(sites)
    assert set(out["variant_id"]) == {f"VAR{i:04d}" for i in range(1, 31)}


def test_top_30_follow_confidence_before_n_alt() -> None:
    high = [vv.variant(f"H{i:03d}", aa_pos=2 + 3 * i, n_alt=2) for i in range(30)]
    cand = [vv.variant("C999", aa_pos=295, confidence="candidate", n_alt=50)]
    assert "C999" not in set(_windows([*cand, *high])["variant_id"])


def test_header_only_variants_give_a_header_only_table() -> None:
    out = _windows([])
    assert list(out.columns) == COLUMNS
    assert len(out) == 0


def test_inputs_are_not_mutated() -> None:
    variants = _variants(vv.main_variants())
    genes = vv.genes_parquet_frame()
    v0, g0 = variants.copy(), genes.copy()
    _fn()(variants, genes)
    pd.testing.assert_frame_equal(variants, v0)
    pd.testing.assert_frame_equal(genes, g0)


def test_deterministic() -> None:
    a = _windows(vv.main_variants())
    b = _windows(vv.main_variants())
    pd.testing.assert_frame_equal(a, b)


# ---------------------------------------------------------------- properties


@SETTINGS
@given(st.integers(1, 400), st.data())
def test_property_window_is_the_clipped_slice_around_the_site(n: int, data: st.DataObject) -> None:
    pos = data.draw(st.integers(1, n))
    aa = vv.protein(n, seed=n)
    site = vv.variant("VAR0001", gene="GPROP", aa_pos=pos, ref_aa=aa[pos - 1], contig="chr2",
                      codon_start=1 + 3 * (pos - 1))  # fmt: skip
    r = _by_id(_windows([site], {"GPROP": aa}))["VAR0001"]
    lo, seq = int(r["window_start_aa"]), str(r["ref_window"])
    assert lo == max(1, pos - 10)
    assert seq == aa[lo - 1 : min(n, pos + 10)]
    assert 1 <= len(seq) <= 21
    assert seq[pos - lo] == aa[pos - 1]


# ---------------------------------------------------------------- V12 reads it back


@pytest.fixture
def _styled() -> None:
    from giap.viz import style

    style.apply()


def test_v12_uses_the_written_support_table(tmp_path: Path, _styled: None) -> None:
    import giap.viz.variants  # noqa: F401
    from giap.config import resolve_config
    from giap.io.plot_support import write_support
    from giap.viz import base

    sites = vv.main_variants()
    out = vv.build_out(tmp_path, sites, vv.main_hits())
    write_support("variant_windows", _windows(sites), out / "plots")
    res = base.registered()["V12"].fn(RunView.from_outdir(out), resolve_config())
    assert not isinstance(res, base.Skip)
    b = res[0].data[res[0].data["panel"] == "b"]
    got = {str(v): str(w) for v, w in zip(b["variant_id"], b["ref_window"], strict=True)}
    assert got["VAR0001"] == vv.P105[203:224]
    assert got["VAR0002"] == vv.P200[39:60]
