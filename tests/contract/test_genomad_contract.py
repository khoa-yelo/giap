"""Contract test: geNOMAD 1.12.0 summary parser on recorded runs (M2-18; GIAP_SPEC.md §7.1.1).

Committed recordings (``tests/data/recorded/genomad/PROVENANCE.md``): the complete
``genome_summary/genome_{virus_summary,plasmid_summary,virus_genes}.tsv`` of sim100k, K-12 MG1655
and Nissle 1917 (prefix ``genome``), plus ``replicons.tsv``.

* sim100k: headers only → ``GenomadResult((), (), ())`` (example 5).
* K-12: 5 proviruses, Nissle: 3 proviruses and 2 plasmids (acceptance: ≥ 1 provirus each).
* Q2: every virus gene lies inside its provirus with the coordinates read as absolute 1-based
  closed contig positions (an ``a-1`` offset would move them out), which decides the rule.
* full recordings under ``$GIAP_ROOT/recordings/genomad/`` (skipped when missing): committed
  copies are byte-identical and the full output directory parses to the same result.
"""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path

import pytest

from giap.model.coords import Interval, Replicon, check_interval

pytestmark = pytest.mark.pending_sample

RECORDED = Path(__file__).resolve().parents[1] / "data" / "recorded" / "genomad" / "1.12.0"
SAMPLES = ["sim100k", "k12_mg1655", "nissle1917"]
TSVS = ["virus_summary", "plasmid_summary", "virus_genes"]
FULL = Path(os.environ.get("GIAP_ROOT", "/nonexistent")) / "recordings" / "genomad"
TAX = "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;"

K12_PROVIRUSES = {  # seq_name coordinates (1-based closed) → score, n_hallmarks
    (1412000, 1434224): (0.9628, 5),
    (1196867, 1213107): (0.9605, 8),
    (2463012, 2476510): (0.9555, 3),
    (563848, 584430): (0.8940, 7),
    (1627517, 1656149): (0.8863, 5),
}
NISSLE_PROVIRUSES = {
    (2171337, 2214647): (0.9785, 19),
    (2869019, 2912533): (0.9638, 19),
    (1694618, 1706222): (0.9626, 3),
}


def _replicons(sample: str) -> dict[str, Replicon]:
    with (RECORDED / sample / "replicons.tsv").open(newline="") as fh:
        return {
            r["contig"]: Replicon(r["contig"], int(r["length"]), r["topology"])  # type: ignore[arg-type]
            for r in csv.DictReader(fh, delimiter="\t")
        }


def _rows(sample: str, name: str, base: Path = RECORDED) -> list[dict[str, str]]:
    path = base / sample / "genome_summary" / f"genome_{name}.tsv"
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _parse(sample: str, base: Path = RECORDED):  # type: ignore[no-untyped-def]
    from giap.annotate.genomad import parse_genomad

    return parse_genomad(base / sample, "genome", _replicons(sample))


def test_sim100k_is_empty() -> None:
    from giap.annotate.genomad import GenomadResult

    assert _parse("sim100k") == GenomadResult((), (), ())


@pytest.mark.parametrize(
    ("sample", "contig", "want"),
    [("k12_mg1655", "U00096.3", K12_PROVIRUSES), ("nissle1917", "CP171242.1", NISSLE_PROVIRUSES)],
)
def test_proviruses(
    sample: str, contig: str, want: dict[tuple[int, int], tuple[float, int]]
) -> None:
    res = _parse(sample)
    assert len(res.proviruses) >= 1
    reps = _replicons(sample)
    got = {}
    for p in res.proviruses:
        assert p.contig == contig and p.interval is not None and p.interval.strand == "."
        check_interval(p.interval, reps[contig])
        assert p.fdr is None  # all NA in these runs
        assert p.taxonomy == TAX
        got[(p.interval.start + 1, p.interval.end)] = (p.virus_score, p.n_hallmarks)
    assert got.keys() == want.keys()
    for k, (score, n) in want.items():
        assert got[k][0] == pytest.approx(score) and got[k][1] == n


@pytest.mark.parametrize("sample", ["k12_mg1655", "nissle1917"])
def test_virus_genes_lie_inside_their_provirus_q2(sample: str) -> None:
    res = _parse(sample)
    rows = _rows(sample, "virus_genes")
    assert len(res.virus_genes) == len(rows) > 0
    by_name = {
        f"{p.contig}|provirus_{p.interval.start + 1}_{p.interval.end}": p.interval
        for p in res.proviruses
        if p.interval is not None
    }
    assert all(isinstance(vg.interval, Interval) for vg in res.virus_genes)
    # every gene (as parsed) inside the provirus its name refers to, same contig
    want = {}
    for r in rows:
        prov = r["gene"].rsplit("_", 1)[0]
        iv = by_name[prov]
        strand = {"1": "+", "-1": "-"}[r["strand"]]
        want[(iv.contig, int(r["start"]) - 1, int(r["end"]), strand)] = iv
    got = {
        (g.contig, g.interval.start, g.interval.end, g.interval.strand): g for g in res.virus_genes
    }
    assert got.keys() == want.keys()
    for key, iv in want.items():
        g = got[key]
        assert g.interval.contig == g.contig == iv.contig
        assert iv.start <= g.interval.start < g.interval.end <= iv.end


@pytest.mark.parametrize("sample", ["k12_mg1655", "nissle1917"])
def test_virus_gene_markers_descriptions_and_hallmarks(sample: str) -> None:
    res = _parse(sample)
    rows = _rows(sample, "virus_genes")
    na = lambda v: None if v == "NA" else v  # noqa: E731
    want = sorted(
        (
            int(r["start"]) - 1,
            int(r["end"]),
            na(r["marker"]),
            na(r["annotation_description"]),
            r["virus_hallmark"] == "1",
        )
        for r in rows
    )
    got = sorted(
        (g.interval.start, g.interval.end, g.marker, g.description, g.virus_hallmark)
        for g in res.virus_genes
    )
    assert got == want
    assert any(g.virus_hallmark for g in res.virus_genes)
    assert any(g.marker is None for g in res.virus_genes)


def test_nissle_plasmids() -> None:
    res = _parse("nissle1917")
    got = {p.contig: p for p in res.plasmids}
    assert set(got) == {"CP171244.1", "CP171245.1"}  # CP171243.1 not called (default filter)
    want = {"CP171244.1": (0.9951, ()), "CP171245.1": (0.9932, ("MOBP1",))}
    for contig, (score, conj) in want.items():
        p = got[contig]
        assert p.plasmid_score == pytest.approx(score)
        assert (p.fdr, p.n_hallmarks, p.conjugation_genes, p.amr_genes) == (None, 2, conj, ())


def test_k12_has_no_plasmids() -> None:
    assert _parse("k12_mg1655").plasmids == ()


@pytest.mark.parametrize("sample", SAMPLES)
def test_parse_is_deterministic(sample: str) -> None:
    assert _parse(sample) == _parse(sample)


# ---------------------------------------------------------------- full recordings (not committed)


@pytest.mark.parametrize("sample", SAMPLES)
def test_full_recording_matches_committed_copy(sample: str) -> None:
    full = FULL / sample / "genomad"
    if not (full / "genome_summary").is_dir():
        pytest.skip(f"full geNOMAD recording missing: {full}")
    for name in TSVS:
        rel = Path("genome_summary") / f"genome_{name}.tsv"
        a = hashlib.sha256((full / rel).read_bytes()).hexdigest()
        b = hashlib.sha256((RECORDED / sample / rel).read_bytes()).hexdigest()
        assert a == b, f"{sample}: {rel} differs from the committed copy"
    from giap.annotate.genomad import parse_genomad

    assert parse_genomad(full, "genome", _replicons(sample)) == _parse(sample)
