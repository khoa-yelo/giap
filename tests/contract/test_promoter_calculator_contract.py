"""Contract test: Promoter Calculator 1.2.4 regions, parser and ``tx_layer`` on recorded runs
(M2-14; GIPPA_SPEC.md §5.4, §7.1.1, §7.1.4).

Committed excerpts (``tests/data/recorded/promoter_calculator/PROVENANCE.md``): 5 sim100k TU
regions (both strands, chromosome and plasmid, the shortest region) and 5 K-12 MG1655 regions
(both strands, and two regions for which the tool wrote **no CSV**), each with the recorded
``regions/<TU>.fa``, ``csv/<TU>.csv`` and the region coordinates in ``tus.tsv``.

* ``region_fasta_text`` reproduces the recorded sim100k region files byte for byte from gippasim's
  regenerated SIM100K genome.
* every recorded row parses; ``PcHit.tss`` is the 1-based +1 base: the 20 nt the tool reports as
  ``ITR`` start there (forward region strand for ``+``, reverse complement for ``-``).
* ``tx_layer`` keeps exactly the ``+`` rows (TU strand), maps them to the genome TSS where the
  ``ITR`` is found on the TU's genome strand, and ``tu_tx_features`` gives max/sum/distance.
* the full K-12 run under ``$GIPPA_ROOT/recordings/promoter_calculator/k12_mg1655`` (not committed;
  skipped when missing): all 2,751 CSVs parse with the ITR property, exactly the 20 known short
  regions have no CSV, and with the recorded Bakta JSON the regions and genome TSSs are checked
  genome-wide.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pandas as pd
import pytest
from helpers.promoter_calculator_fixtures import (
    HEADER,
    genome_itr,
    header_only,
    itr_at,
    read_dicts,
    read_tus_tsv,
    recorder_tus,
    region_seq,
)

from gippa.annotate.operons import Tu
from gippa.model.coords import Replicon

DATA = Path(__file__).resolve().parents[1] / "data" / "recorded"
RECORDED = DATA / "promoter_calculator" / "1.2.4"
SAMPLES = ["sim100k", "k12_mg1655"]
NO_CSV = {"sim100k": set(), "k12_mg1655": {"TU_00058", "TU_00318"}}
SIM_GENOME_SHA = "dd10216ee21cb0ac13940b2f7d5a26bf0fdc8c2af54af8baf26b380fc65a8d69"
GIPPA_ROOT = Path(os.environ.get("GIPPA_ROOT", "/nonexistent")) / "recordings"
FULL_K12 = GIPPA_ROOT / "promoter_calculator" / "k12_mg1655"
K12_JSON = GIPPA_ROOT / "bakta_full" / "k12_mg1655" / "bakta" / "k12_mg1655.json"
K12_JSON_SHA = "13c54a8961651053cda596ffa4cc0c611565588112275faf7f3ddc27e28bacec"
K12_NO_CSV = {
    "TU_00058", "TU_00059", "TU_00318", "TU_00319", "TU_00493", "TU_00494", "TU_00768",
    "TU_00769", "TU_01213", "TU_01214", "TU_01217", "TU_01218", "TU_01279", "TU_01280",
    "TU_01865", "TU_01866", "TU_02029", "TU_02030", "TU_02612", "TU_02613",
}  # fmt: skip


def _pc():
    from gippa.annotate import promoter_calculator as pc

    return pc


def _sample(sample: str) -> tuple[list[Tu], dict[str, Replicon], dict[str, int]]:
    return read_tus_tsv(RECORDED / sample / "tus.tsv")


def _genome_tss(tu: Tu, tss1: int, size: int) -> int:
    """Independent oracle of the spec's region → genome rule."""
    a, b = tu.upstream.start, tu.upstream.end
    return (a + tss1 - 1 if tu.strand == "+" else b - tss1) % size


def _oracle(tus: list[Tu], reps: dict[str, Replicon], csv_dir: Path) -> set[tuple[object, ...]]:
    """(contig, strand, tss, tx_rate, tu_id) for every kept recorded row."""
    want = set()
    for tu in tus:
        path = csv_dir / f"{tu.tu_id}.csv"
        if not path.exists():
            continue
        n = len(region_seq(csv_dir.parent / "regions" / f"{tu.tu_id}.fa"))
        for r in read_dicts(path):
            if r["strand"] != "+":
                continue
            tss1 = int(r["TSS"]) + 1
            assert 1 <= tss1 <= n
            g = _genome_tss(tu, tss1, reps[tu.contig].length)
            want.add((tu.contig, tu.strand, g, float(r["Tx_rate"]), tu.tu_id))
    return want


def _rows(df: pd.DataFrame) -> set[tuple[object, ...]]:
    cols = ["contig", "strand", "tss", "tx_rate", "tu_id"]
    return {
        (str(c), str(s), int(t), float(x), str(u))
        for c, s, t, x, u in df[cols].itertuples(index=False, name=None)
    }


def _hits(tus: list[Tu], csv_dir: Path) -> dict[str, list[object]]:
    pc = _pc()
    return {
        t.tu_id: pc.parse_pc_csv(csv_dir / f"{t.tu_id}.csv")
        for t in tus
        if (csv_dir / f"{t.tu_id}.csv").exists()
    }


def _itr_map(tus: list[Tu], csv_dir: Path) -> dict[tuple[str, float], set[str]]:
    """(tu_id, Tx_rate) → recorded ITRs of the + rows."""
    out: dict[tuple[str, float], set[str]] = {}
    for t in tus:
        path = csv_dir / f"{t.tu_id}.csv"
        for r in read_dicts(path) if path.exists() else []:
            if r["strand"] == "+":
                out.setdefault((t.tu_id, float(r["Tx_rate"])), set()).add(r["ITR"])
    return out


def _check_csv(csv_path: Path, region: str) -> int:
    pc = _pc()
    rows = read_dicts(csv_path)
    hits = pc.parse_pc_csv(csv_path)
    assert len(hits) == len(rows), csv_path.name
    for h, r in zip(hits, rows, strict=True):
        assert h.strand == r["strand"] and h.strand in {"+", "-"}
        assert h.tx_rate == float(r["Tx_rate"])
        assert 1 <= h.tss <= len(region)
        assert itr_at(region, h.tss, h.strand) == r["ITR"], (csv_path.name, r["TSS"], h.strand)
    return len(hits)


# ---------------------------------------------------------------- committed excerpts


def test_recorded_header_is_the_27_column_header() -> None:
    hdr = (DATA / "headers" / "promoter_calculator.csv.txt").read_text().strip()
    assert hdr == HEADER
    for sample in SAMPLES:
        for path in sorted((RECORDED / sample / "csv").glob("*.csv")):
            assert path.read_bytes().split(b"\r\n", 1)[0].decode() == hdr, path


def test_excerpt_files_are_complete() -> None:
    for sample in SAMPLES:
        tus, _, _ = _sample(sample)
        assert len(tus) == 5
        csvs = {p.stem for p in (RECORDED / sample / "csv").glob("*.csv")}
        assert csvs == {t.tu_id for t in tus} - NO_CSV[sample]
        assert {t.strand for t in tus} == {"+", "-"}
        for t in tus:
            assert (RECORDED / sample / "regions" / f"{t.tu_id}.fa").is_file()


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_recorded_row_parses_with_the_tss_at_the_itr(sample: str) -> None:
    tus, _, _ = _sample(sample)
    total = 0
    for t in tus:
        csv_path = RECORDED / sample / "csv" / f"{t.tu_id}.csv"
        if csv_path.exists():
            total += _check_csv(
                csv_path, region_seq(RECORDED / sample / "regions" / f"{t.tu_id}.fa")
            )
    assert total > 0


@pytest.mark.parametrize("sample", SAMPLES)
def test_recorded_region_lengths_match_the_coordinates(sample: str) -> None:
    tus, _, _ = _sample(sample)
    for t in tus:
        n = len(region_seq(RECORDED / sample / "regions" / f"{t.tu_id}.fa"))
        assert n == t.upstream.end - t.upstream.start, t.tu_id


def test_no_csv_regions_get_a_header_only_file_with_no_hits(tmp_path: Path) -> None:
    pc = _pc()
    for tu_id in sorted(NO_CSV["k12_mg1655"]):
        assert pc.parse_pc_csv(header_only(tmp_path / f"{tu_id}.csv")) == []


@pytest.mark.parametrize("sample", SAMPLES)
def test_tx_layer_on_recorded_hits(sample: str) -> None:
    pc = _pc()
    tus, reps, _ = _sample(sample)
    csv_dir = RECORDED / sample / "csv"
    tx = pc.tx_layer(tus, _hits(tus, csv_dir), reps)
    want = _oracle(tus, reps, csv_dir)
    assert _rows(tx) == want
    assert len(tx) == len(want)
    # the kept +1 base is where the recorded ITR starts in the region text (TU orientation)
    itr = _itr_map(tus, csv_dir)
    by_id = {t.tu_id: t for t in tus}
    for tu_id, g, x in zip(tx["tu_id"], tx["tss"], tx["tx_rate"], strict=True):
        t = by_id[str(tu_id)]
        text = region_seq(RECORDED / sample / "regions" / f"{tu_id}.fa")
        idx = int(g) - t.upstream.start if t.strand == "+" else t.upstream.end - 1 - int(g)
        assert text[idx : idx + 20] in itr[(str(tu_id), float(x))]


@pytest.mark.parametrize("sample", SAMPLES)
def test_tu_features_on_recorded_hits(sample: str) -> None:
    pc = _pc()
    tus, reps, sc = _sample(sample)
    csv_dir = RECORDED / sample / "csv"
    tx = pc.tx_layer(tus, _hits(tus, csv_dir), reps)
    feats = pc.tu_tx_features(tus, tx, replicons=reps)
    assert list(feats["tu_id"]) == [t.tu_id for t in tus]
    want = _oracle(tus, reps, csv_dir)
    for t in tus:
        f = feats[feats["tu_id"] == t.tu_id].iloc[0]
        mine = [(g, x) for _, _, g, x, u in want if u == t.tu_id]
        if not mine:
            assert f["tx_rate_max"] == 0 and f["tx_rate_sum"] == 0
            assert pd.isna(f["tss_distance_bp"])
            continue
        best = max(x for _, x in mine)
        assert f["tx_rate_max"] == best
        assert f["tx_rate_sum"] == pytest.approx(sum(x for _, x in mine), rel=1e-12)
        dist = {
            sc[t.tu_id] - g if t.strand == "+" else g - sc[t.tu_id] for g, x in mine if x == best
        }
        assert f["tss_distance_bp"] in dist
        assert -50 <= int(f["tss_distance_bp"]) <= 300


def test_sim100k_regions_byte_for_byte_and_genome_tss(tmp_path: Path) -> None:
    import gippasim
    from gippa.io.fasta import read_genome_fasta

    pc = _pc()
    res = gippasim.simulate(gippasim.SIM100K, gippasim.SIM100K_SEED, tmp_path / "sim")
    fa = Path(res.genome_fasta)
    if hashlib.sha256(fa.read_bytes()).hexdigest() != SIM_GENOME_SHA:
        pytest.fail("gippasim SIM100K genome changed: the recording no longer matches it")
    seq = {r.id: r.sequence for r in read_genome_fasta(fa)}
    tus, reps, _ = _sample("sim100k")
    assert {c: r.length for c, r in reps.items()} == {c: len(s) for c, s in seq.items()}
    for t in tus:
        want = (RECORDED / "sim100k" / "regions" / f"{t.tu_id}.fa").read_text()
        assert pc.region_fasta_text(t, seq) == want, t.tu_id
    csv_dir = RECORDED / "sim100k" / "csv"
    tx = pc.tx_layer(tus, _hits(tus, csv_dir), reps)
    itr = _itr_map(tus, csv_dir)
    for contig, strand, g, x, tu_id in _rows(tx):
        assert genome_itr(seq[str(contig)], int(g), str(strand)) in itr[(str(tu_id), float(x))]


# ---------------------------------------------------------------- full K-12 run (not committed)


@pytest.mark.skipif(not (FULL_K12 / "csv").is_dir(), reason="full K-12 recording not present")
def test_full_k12_every_csv_parses_and_the_missing_ones_are_the_short_regions() -> None:
    regions = {p.stem: p for p in (FULL_K12 / "regions").glob("*.fa")}
    csvs = {p.stem: p for p in (FULL_K12 / "csv").glob("*.csv")}
    assert len(regions) == 2771 and len(csvs) == 2751
    assert set(regions) - set(csvs) == K12_NO_CSV
    total = sum(_check_csv(p, region_seq(regions[s])) for s, p in sorted(csvs.items()))
    assert total == 238_387
    assert all(len(region_seq(regions[s])) < 83 for s in K12_NO_CSV)  # shorter than any promoter


@pytest.mark.slow  # ~75 s: 2,771 regions and ~120k genome TSSs
@pytest.mark.skipif(
    not (FULL_K12 / "csv").is_dir() or not K12_JSON.is_file(),
    reason="full K-12 recording or its Bakta JSON not present",
)
def test_full_k12_regions_and_genome_tss() -> None:
    if hashlib.sha256(K12_JSON.read_bytes()).hexdigest() != K12_JSON_SHA:
        pytest.fail("recorded K-12 Bakta JSON changed")
    pc = _pc()
    tus, reps, seqs = recorder_tus(K12_JSON)
    assert len(tus) == 2771
    for t in tus:
        assert pc.region_fasta_text(t, seqs) == (FULL_K12 / "regions" / f"{t.tu_id}.fa").read_text()
    csv_dir = FULL_K12 / "csv"
    tx = pc.tx_layer(tus, _hits(tus, csv_dir), reps)
    want = _oracle(tus, reps, csv_dir)
    assert _rows(tx) == want
    assert len(tx) == sum(
        1 for t in tus if (csv_dir / f"{t.tu_id}.csv").exists()
        for r in read_dicts(csv_dir / f"{t.tu_id}.csv") if r["strand"] == "+"
    )  # fmt: skip
    itr = _itr_map(tus, csv_dir)
    for contig, strand, g, x, tu_id in _rows(tx):
        assert genome_itr(seqs[str(contig)], int(g), str(strand)) in itr[(str(tu_id), float(x))]
