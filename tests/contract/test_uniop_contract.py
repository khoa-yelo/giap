"""Contract test: ``uniop.pred`` parser on recorded UniOP runs (M2-12; GIAP_SPEC.md §5.4, §7.1.1).

Marked ``pending_sample`` until the implementer records a real UniOP run (pin 3ed14dfc, built
``giap-uniop`` env, run as a process; lead decision 3 after the PR #49 review) under
``tests/data/recorded/uniop/3ed14dfc/sim100k/`` (``k12/`` optional). Each directory holds the exact
input GFF written by ``write_prodigal_gff`` (``*.gff``), UniOP's ``uniop.pred``, ``uniop_map.tsv``
(``prodigal_id``, ``locus_tag``) and ``PROVENANCE.md`` (command, env lock, input sha256).

Real format (lead decision 1): title ``operon prediction scores``, header
``Gene A<TAB>Gene B<TAB>Prediction``, rows ``k<TAB>k+1<TAB>p`` of 1-based row numbers. Row numbers
index UniOP's own gene order (lead decision 2): all GFF genes stable-sorted by start across contigs;
row ``k`` pairs sorted gene ``k`` with the next gene on the same contig; its probability is empty
exactly when that partner is missing or on the other strand.

§7.1.1: UniOP keeps "probability for each adjacent same-strand pair", so the parsed pairs are
exactly the consecutive same-strand genes of each contig in the GFF (never across contigs, never
the wrap).
"""

from __future__ import annotations

import hashlib
import itertools
from pathlib import Path

import pytest

RECORDED = Path(__file__).resolve().parents[1] / "data" / "recorded" / "uniop" / "3ed14dfc"
EXPECTED = ["sim100k", "k12"]
SAMPLES = [RECORDED / d for d in EXPECTED if (RECORDED / d).is_dir()]
TITLE = "operon prediction scores"
HEADER = ["Gene A", "Gene B", "Prediction"]

pytestmark = [
    pytest.mark.pending_sample,
    pytest.mark.skipif(not SAMPLES, reason=f"pending_sample: no recorded UniOP run in {RECORDED}"),
]


def _gff(path: Path) -> list[tuple[str, str, int, str]]:
    """(prodigal ID, contig, 1-based start, strand) per CDS line in file order."""
    out = []
    for ln in path.read_text().splitlines():
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        out.append((f[8].split(";")[0].removeprefix("ID="), f[0], int(f[3]), f[6]))
    return out


def _files(sample: Path) -> tuple[Path, Path, dict[str, str]]:
    from giap.io.tsv import read_tsv

    gffs = sorted(sample.glob("*.gff"))
    pred = sample / "uniop.pred"
    assert len(gffs) == 1 and pred.is_file(), f"{sample}: need one *.gff and uniop.pred"
    header, mrows = read_tsv(sample / "uniop_map.tsv")
    assert header == ["prodigal_id", "locus_tag"]
    return gffs[0], pred, {str(r[0]): str(r[1]) for r in mrows}


@pytest.mark.parametrize("sample", SAMPLES, ids=[p.name for p in SAMPLES])
def test_recorded_run_has_provenance(sample: Path) -> None:
    gff, _, _ = _files(sample)
    prov = (sample / "PROVENANCE.md").read_text()
    assert "3ed14dfc" in prov
    assert hashlib.sha256(gff.read_bytes()).hexdigest() in prov


@pytest.mark.parametrize("sample", SAMPLES, ids=[p.name for p in SAMPLES])
def test_recorded_uniop_pred_layout_follows_uniop_gene_order(sample: Path) -> None:
    gff, pred, _ = _files(sample)
    genes = _gff(gff)
    lines = pred.read_text().splitlines()
    assert lines[0].strip() == TITLE
    assert lines[1].split("\t") == HEADER
    rows = [ln.split("\t") for ln in lines[2:] if ln.strip()]
    assert [(int(r[0]), int(r[1])) for r in rows] == [(k, k + 1) for k in range(1, len(genes))]
    order = sorted(genes, key=lambda g: g[2])  # stable: file order breaks ties
    for k, r in enumerate(rows):
        a = order[k]
        nxt = next((g for g in order[k + 1 :] if g[1] == a[1]), None)
        same_strand = nxt is not None and nxt[3] == a[3]
        cell = r[2].strip()
        assert (cell != "") == same_strand, (r, a, nxt)
        if cell:
            assert 0.0 <= float(cell) <= 1.0


@pytest.mark.parametrize("sample", SAMPLES, ids=[p.name for p in SAMPLES])
def test_recorded_uniop_pred_parses_to_adjacent_pairs(sample: Path) -> None:
    from giap.annotate.uniop_adapter import parse_uniop_pred

    gff, pred, id_map = _files(sample)
    rows = _gff(gff)
    assert sorted(id_map) == sorted(r[0] for r in rows)
    pairs = parse_uniop_pred(pred, gff, id_map)
    assert pairs, "no pair parsed"
    expected = {
        (id_map[a[0]], id_map[b[0]])
        for a, b in itertools.pairwise(rows)
        if a[1] == b[1] and a[3] == b[3]
    }
    got = [(p.gene_a, p.gene_b) for p in pairs]
    assert len(set(got)) == len(got), "duplicate pairs"
    assert set(got) == expected
    assert all(p.source == "uniop" and 0.0 <= p.uniop_p <= 1.0 for p in pairs)
    cells = [ln.split("\t")[2] for ln in pred.read_text().splitlines()[2:] if ln.strip()]
    assert len(pairs) == sum(1 for c in cells if c.strip())
