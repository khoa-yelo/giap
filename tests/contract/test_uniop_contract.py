"""Contract test: ``uniop.pred`` parser on recorded UniOP runs (M2-12; GIAP_SPEC.md §5.4, §7.1.1).

Marked ``pending_sample`` until M2-10's harness records
``tests/data/recorded/uniop/3ed14dfc2a5d/{sim100k,k12}/`` (UniOP pinned at 3ed14dfc, M2-01) on the
Bakta full-DB genes. Each directory holds the input GFF written by ``write_prodigal_gff``
(``*.gff``), ``uniop.pred`` and, when the harness writes it, ``uniop_map.tsv`` (``prodigal_id``,
``locus_tag``); without a map the Prodigal IDs map to themselves.

§7.1.1: UniOP keeps "probability for each adjacent same-strand pair", so every parsed pair is a
consecutive same-strand pair of the GFF (never across contigs, never the wrap pair) and every such
pair is reported.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import pytest

RECORDED = Path(__file__).resolve().parents[1] / "data" / "recorded" / "uniop" / "3ed14dfc2a5d"
EXPECTED = ["sim100k", "k12"]
SAMPLES = [RECORDED / d for d in EXPECTED if (RECORDED / d).is_dir()]

pytestmark = [
    pytest.mark.pending_sample,
    pytest.mark.skipif(not SAMPLES, reason=f"pending_sample: no recorded UniOP run in {RECORDED}"),
]


def _gff(path: Path) -> list[tuple[str, str, str]]:
    out = []
    for ln in path.read_text().splitlines():
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        out.append((f[8].split(";")[0].removeprefix("ID="), f[0], f[6]))
    return out


@pytest.mark.parametrize("sample", SAMPLES, ids=[p.name for p in SAMPLES])
def test_recorded_uniop_pred_parses_to_adjacent_pairs(sample: Path) -> None:
    from giap.annotate.uniop import parse_uniop_pred
    from giap.io.tsv import read_tsv

    gffs = sorted(sample.glob("*.gff"))
    pred = sample / "uniop.pred"
    assert gffs and pred.is_file(), f"{sample}: need *.gff and uniop.pred"
    rows = _gff(gffs[0])
    ids = [r[0] for r in rows]
    if (sample / "uniop_map.tsv").is_file():
        header, mrows = read_tsv(sample / "uniop_map.tsv")
        assert header == ["prodigal_id", "locus_tag"]
        id_map = {str(r[0]): str(r[1]) for r in mrows}
    else:
        id_map = {i: i for i in ids}
    pairs = parse_uniop_pred(pred, id_map)
    assert pairs, "no pair parsed"
    expected = {
        (id_map[a[0]], id_map[b[0]])
        for a, b in itertools.pairwise(rows)
        if a[1] == b[1] and a[2] == b[2]
    }
    got = [(p.gene_a, p.gene_b) for p in pairs]
    assert len(set(got)) == len(got), "duplicate pairs"
    assert set(got) == expected
    assert all(p.source == "uniop" and 0.0 <= p.uniop_p <= 1.0 for p in pairs)
