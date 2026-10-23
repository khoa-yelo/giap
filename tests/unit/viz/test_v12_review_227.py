"""PR #227 review B1 (lead decision 11): an ALT site no gene frames (``gene`` / ``aa_pos`` NA,
e.g. a sixframe or novel-ORF peptide) neither crashes ``variant_windows`` nor V12."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import viz_variants as vv

from giap.config import resolve_config
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.support import variant_windows


def _v12(view: RunView) -> list[base.Instance] | base.Skip:
    import giap.viz.variants  # noqa: F401  (registers V12)

    return base.registered()["V12"].fn(view, resolve_config())


def _no_gene(vid: str) -> dict[str, object]:
    row = vv.variant(vid, aa_pos=30, codon_start=7001)
    row.update(gene=None, aa_pos=None)
    return row


def test_variant_windows_leaves_out_a_site_without_a_gene(tmp_path: Path) -> None:
    rows = [vv.variant("VAR0001"), _no_gene("VAR0002")]
    out = vv.build_out(tmp_path, rows, [])
    table = RunView.from_outdir(out).table("variants")
    assert table is not None
    got = variant_windows(table, vv.genes_parquet_frame())
    assert list(got["variant_id"]) == ["VAR0001"]


def test_v12_shows_a_site_without_a_gene_unaligned(tmp_path: Path) -> None:
    style.apply()
    hit = vv.plus_hit("MKTAYIAK", 6995, 1)
    out = vv.build_out(tmp_path, [vv.variant("VAR0001"), _no_gene("VAR0002")], [hit])
    result = _v12(RunView.from_outdir(out))
    assert not isinstance(result, base.Skip), result
    inst = result[0]
    b = inst.data[inst.data["panel"] == "b"]
    assert set(b["variant_id"]) == {"VAR0001", "VAR0002"}
    assert pd.isna(b.set_index("variant_id").loc["VAR0002", "ref_window"])
    text = "".join(t.get_text() for ax in inst.figure.axes for t in ax.texts)
    assert "no gene" in text


@pytest.mark.parametrize("n", [1, 3])
def test_only_no_gene_sites_still_render(tmp_path: Path, n: int) -> None:
    style.apply()
    out = vv.build_out(tmp_path, [_no_gene(f"VAR{i:04d}") for i in range(n)], [])
    result = _v12(RunView.from_outdir(out))
    assert not isinstance(result, base.Skip), result
