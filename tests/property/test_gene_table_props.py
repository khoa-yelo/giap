"""Hypothesis properties of the M4-15 gene table (docs/specs/M4-15-gene-table-assembly.md;
GIPPA_SPEC.md §6.5 sort order and writer coordinates, §9.4 closed flag vocabulary).

- ``gene_flags``: the union of its parts, unique, in vocabulary order, idempotent; any word outside
  the §9.4 Gene vocabulary raises ``ValueError``.
- ``function_columns``: invariant under row order and duplicate rows of C3; list cells are the
  sorted distinct IDs per gene.
- ``assemble_genes``: the output does not depend on the row order of any input and is §6.5-sorted.
- ``write_genes_tsv`` / ``read_genes_tsv``: TSV ``start = start + 1``, ``end`` unchanged, and the
  internal coordinates come back.
"""

from __future__ import annotations

import dataclasses
import functools
import sys
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import gene_table_fixture as gt
from gene_table_fixture import as_list

FLAGS = st.sampled_from(gt.GENE_FLAGS)
SETTINGS = settings(max_examples=40, deadline=None,
                    suppress_health_check=[HealthCheck.too_slow])  # fmt: skip


@functools.cache
def _world() -> gt.World:
    fr = gt.frames_with()
    g = fr["genes"]
    g.loc[1, ["start", "end", "strand"]] = [450, 750, "-"]  # ties on start/end
    g.loc[2, ["start", "end", "strand"]] = [450, 750, "+"]
    g.loc[5, ["start", "end"]] = [1500, 1800]  # same as gene 4
    return gt.build_world(frames=fr)


@functools.cache
def _table() -> pd.DataFrame:
    return _world().assemble()


@SETTINGS
@given(st.dictionaries(st.sampled_from(["layers", "evidence", "abundance", "detectability"]),
                       st.lists(FLAGS, max_size=8), max_size=4))  # fmt: skip
def test_gene_flags_union_order_idempotent(parts: dict[str, list[str]]) -> None:
    from gippa.genetable.assemble import gene_flags

    out = gene_flags(parts)
    assert set(out) == {f for v in parts.values() for f in v}
    assert len(out) == len(set(out))
    idx = [gt.GENE_FLAGS.index(f) for f in out]
    assert idx == sorted(idx)
    assert gene_flags({"all": out}) == out


@SETTINGS
@given(st.lists(FLAGS, max_size=4), st.text(min_size=0, max_size=12))
def test_gene_flags_closed_vocabulary(ok: list[str], word: str) -> None:
    import pytest

    from gippa.genetable.assemble import gene_flags

    if word in gt.GENE_FLAGS:
        return
    with pytest.raises(ValueError):
        gene_flags({"x": [*ok, word]})


_TAGS = st.sampled_from(["G1", "G2", "G3"])
_XREF = st.one_of(
    st.tuples(_TAGS, st.just("GO"), st.integers(0, 30).map(lambda i: f"GO:{i:07d}")),
    st.tuples(_TAGS, st.just("KEGG"), st.integers(0, 30).map(lambda i: f"K{i:05d}")),
    st.tuples(_TAGS, st.just("EC"), st.integers(0, 9).map(lambda i: f"1.1.1.{i}")),
)


@SETTINGS
@given(st.lists(_XREF, min_size=1, max_size=25), st.randoms(use_true_random=False))
def test_function_columns_order_free(rows: list[tuple[str, str, str]], rnd: Any) -> None:
    from gippa.genetable.assemble import function_columns

    cols = ["locus_tag", "db", "id"]
    a = function_columns(pd.DataFrame(rows, columns=cols), None)
    shuffled = rows + rows[: len(rows) // 2]
    rnd.shuffle(shuffled)
    b = function_columns(pd.DataFrame(shuffled, columns=cols), None)
    key = {"go_terms": "GO", "kegg_ko": "KEGG", "ec": "EC"}
    for t in sorted({r[0] for r in rows}):
        ra = a[a["locus_tag"] == t].iloc[0]
        rb = b[b["locus_tag"] == t].iloc[0]
        for c, db in key.items():
            want = sorted({r[2] for r in rows if r[0] == t and r[1] == db})
            assert as_list(ra[c]) == want == as_list(rb[c]), (t, c)


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.permutations(range(20)))
def test_assemble_independent_of_input_order(perm: list[int]) -> None:
    w = _world()

    def sh(df: pd.DataFrame) -> pd.DataFrame:
        return df.iloc[perm].reset_index(drop=True)

    inputs = w.with_inputs(genes=sh(w.inputs.genes), layers=sh(w.inputs.layers),
                           xrefs=w.inputs.xrefs.iloc[::-1].reset_index(drop=True))  # fmt: skip
    out = w.assemble(
        inputs=inputs, cov=sh(w.cov), evidence=sh(w.evidence),
        detection=dataclasses.replace(w.detection, table=sh(w.detection.table)),
        abundance=dataclasses.replace(w.abundance, table=sh(w.abundance.table)),
    )  # fmt: skip
    base = _table()
    pd.testing.assert_frame_equal(out.reset_index(drop=True), base.reset_index(drop=True))
    keys = gt.sort_key(out, ["chr1"])
    assert keys == sorted(keys)


@SETTINGS
@given(st.lists(st.tuples(st.integers(0, 10**7), st.integers(1, 30_000)), min_size=20,
                max_size=20))  # fmt: skip
def test_writer_converts_coordinates_once(spans: list[tuple[int, int]]) -> None:
    from gippa.io.genes_tsv import read_genes_tsv, write_genes_tsv

    df = _table().copy()
    df["start"] = [s for s, _ in spans]
    df["end"] = [s + n for s, n in spans]
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "genes.tsv"
        write_genes_tsv(path, df)
        cells = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        back = read_genes_tsv(path)
    by_tag = df.set_index("locus_tag")
    for t, s, e in zip(cells["locus_tag"], cells["start"], cells["end"], strict=True):
        assert int(s) == int(by_tag.loc[t, "start"]) + 1
        assert int(e) == int(by_tag.loc[t, "end"])
    bb = back.set_index("locus_tag")
    for t in by_tag.index:
        assert (int(bb.loc[t, "start"]), int(bb.loc[t, "end"])) == (
            int(by_tag.loc[t, "start"]), int(by_tag.loc[t, "end"]),
        )  # fmt: skip
