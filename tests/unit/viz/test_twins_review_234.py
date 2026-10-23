"""PR #234 review B1: data TSVs are read with identifier columns as text (lead fix).

A contig named ``01`` and locus tags such as ``0001`` keep their spelling through
``giap.steps.plots.twins``: the V2 subplot is titled ``01``, the replicon length from ``meta`` is
found (Mb scale), and the hover shows ``0001``.
"""

from __future__ import annotations

import json
from pathlib import Path

from helpers import twins as tw

from giap.io.plot_files import read_plot_files
from giap.steps.plots import load_plot_modules, twins


def _rows() -> list[tw.Row]:
    return [(r[0], r[1], None if r[2] is None else f"000{r[2][1:]}", *r[3:])
            for r in tw.v02_rows(contig="01")]  # fmt: skip


def test_numeric_looking_ids_keep_their_spelling(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    (plots / "data").mkdir(parents=True)
    (plots / "data" / "V02_abundance_landscape.tsv").write_text(
        tw.tsv_text(tw.V02_COLUMNS, _rows()), encoding="utf-8"
    )
    load_plot_modules()
    meta = {"replicons": {"01": 5_000_000}}
    out = twins(plots, read_plot_files(plots), meta)
    ((_, text),) = out.values()
    fig = json.loads(text)
    titles = [a.get("text") for a in fig["layout"].get("annotations", [])]
    blob = json.dumps(fig)
    assert "01" in titles or '"01"' in blob, titles
    assert any("0001" in h for h in tw.all_hovers(fig)), "locus tag lost its leading zeros"
    assert "Mb" in blob, "replicon length from meta not found for contig 01"
