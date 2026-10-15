"""Contract test: OSTIR 1.1.3 input writer and parser on recorded runs (M2-16; GIAP_SPEC.md §5.4,
§7.1.1, §7.1.6).

Committed recordings (``tests/data/recorded/ostir/PROVENANCE.md``): the complete sim100k run and
excerpts of the K-12 MG1655 and Nissle 1917 runs (both strands, sORFs, pseudogenes, all Nissle
plasmid genes, the last K-12 gene, and the one gene per genome for which OSTIR wrote no row).

* sim100k: the M2-16 writer reproduces the recorded ``ostir_input.csv`` **byte for byte** from
  the Bakta gene coordinates (``genes.tsv``) and giapsim's regenerated SIM100K genome.
* the parser keeps every recorded row (all at ``start_position == start == 51``) with exactly the
  recorded values; genes without a row are absent.
* full K-12 / Nissle runs under ``$GIAP_ROOT/recordings/ostir/`` (not committed; skipped when
  missing): every row parses, exactly one gene each is missing, and (when the recorded Bakta JSON is
  there too) the writer reproduces the full recorded inputs byte for byte.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
from pathlib import Path

import pytest

from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

pytestmark = pytest.mark.pending_sample

RECORDED = Path(__file__).resolve().parents[1] / "data" / "recorded" / "ostir" / "1.1.3"
SAMPLES = ["sim100k", "k12_mg1655", "nissle1917"]
NO_ROW = {"sim100k": set(), "k12_mg1655": {"KOIAIB_02000"}, "nissle1917": {"JMJJOG_02078"}}
ASD = "ACCTCCTTA"
SIM_GENOME_SHA = "dd10216ee21cb0ac13940b2f7d5a26bf0fdc8c2af54af8baf26b380fc65a8d69"
GIAP_ROOT = Path(os.environ.get("GIAP_ROOT", "/nonexistent")) / "recordings"
FULL = {
    "k12_mg1655": GIAP_ROOT / "ostir" / "k12_mg1655",
    "nissle1917": GIAP_ROOT / "ostir" / "nissle1917",
}


def _dicts(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _genes(sample: str) -> tuple[list[Gene], dict[str, Replicon]]:
    genes, reps = [], {}
    with (RECORDED / sample / "genes.tsv").open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            reps[r["contig"]] = Replicon(r["contig"], int(r["contig_length"]), r["topology"])  # type: ignore[arg-type]
            genes.append(
                Gene(
                    locus_tag=r["locus_tag"],
                    feature_type=r["feature_type"],  # type: ignore[arg-type]
                    interval=Interval(r["contig"], int(r["start"]), int(r["end"]), r["strand"]),  # type: ignore[arg-type]
                    aa="",
                    nt=None,
                    start_codon=r["start_codon"],
                )
            )
    return genes, reps


def _requested(inp: Path) -> dict[str, int]:
    return {r["id"]: int(r["start"]) for r in _dicts(inp)}


def _check_parsed(inp: Path, out: Path, no_row: set[str]) -> None:
    from giap.annotate.ostir import RbsFeatures, parse_ostir

    requested = _requested(inp)
    got = parse_ostir(out, requested)
    assert set(got) == set(requested) - no_row
    for r in _dicts(out):
        assert int(r["start_position"]) == requested[r["name"]]
        f = got[r["name"]]
        assert isinstance(f, RbsFeatures)
        assert f.tir_ostir == float(r["expression"])
        assert f.sd_dg == float(r["dG_rRNA:mRNA"])
        assert f.mrna_dg == float(r["dG_mRNA"])
        assert f.spacing_dg == float(r["dG_spacing"])
        assert f.standby_dg == float(r["dG_standby"])
        assert f.spacer_bp == int(r["RBS_distance_bp"]) and isinstance(f.spacer_bp, int)
        assert all(math.isfinite(x) for x in (f.tir_ostir, f.sd_dg, f.mrna_dg))


# ---------------------------------------------------------------- committed recordings


def test_recorded_files_present() -> None:
    for s in SAMPLES:
        for name in ("ostir_input.csv", "ostir.csv", "genes.tsv", "command.txt", "versions.txt"):
            assert (RECORDED / s / name).is_file(), f"{s}/{name}"
    assert "ostir 1.1.3" in (RECORDED / "sim100k" / "versions.txt").read_text()
    head = (RECORDED / "sim100k" / "ostir.csv").read_text().splitlines()[0]
    hdr = Path(__file__).resolve().parents[1] / "data" / "recorded" / "headers"
    assert head == (hdr / "ostir.ostir.csv.txt").read_text().strip()


def test_sim100k_writer_reproduces_the_recorded_input_byte_for_byte(tmp_path: Path) -> None:
    import giapsim
    from giap.annotate.ostir import ostir_row, write_ostir_input
    from giap.io.fasta import read_genome_fasta

    res = giapsim.simulate(giapsim.SIM100K, giapsim.SIM100K_SEED, tmp_path / "sim")
    fa = Path(res.genome_fasta)
    if hashlib.sha256(fa.read_bytes()).hexdigest() != SIM_GENOME_SHA:
        pytest.fail("giapsim SIM100K genome changed: the recording no longer matches it")
    seq = {r.id: r.sequence for r in read_genome_fasta(fa)}
    genes, reps = _genes("sim100k")
    rows = [ostir_row(g, reps[g.interval.contig], seq[g.interval.contig], ASD) for g in genes]
    out = tmp_path / "ostir_input.csv"
    write_ostir_input(out, rows)
    assert out.read_bytes() == (RECORDED / "sim100k" / "ostir_input.csv").read_bytes()


@pytest.mark.parametrize("sample", SAMPLES)
def test_recorded_inputs_follow_the_geometry(sample: str) -> None:
    genes, _ = _genes(sample)
    by_id = {g.locus_tag: g for g in genes}
    rows = _dicts(RECORDED / sample / "ostir_input.csv")
    assert [r["id"] for r in rows] == [g.locus_tag for g in genes]
    for r in rows:
        g = by_id[r["id"]]
        s = int(r["start"])
        assert (r["end"], r["anti-Shine-Dalgarno"]) == (r["start"], ASD)
        assert r["seq"][s - 1 : s + 2] == g.start_codon
        assert (len(r["seq"]), s) == (100, 51)  # no recorded gene is within 50 nt of a contig end
    assert {g.interval.strand for g in genes} == {"+", "-"}


@pytest.mark.parametrize("sample", SAMPLES)
def test_recorded_outputs_parse(sample: str) -> None:
    d = RECORDED / sample
    _check_parsed(d / "ostir_input.csv", d / "ostir.csv", NO_ROW[sample])


def test_sim100k_first_row_values() -> None:
    from giap.annotate.ostir import RbsFeatures, parse_ostir

    d = RECORDED / "sim100k"
    got = parse_ostir(d / "ostir.csv", _requested(d / "ostir_input.csv"))
    assert len(got) == 86
    assert got["KABHIH_0001"] == RbsFeatures(1129.1028, -12.781, -14.6, 0.0, 0.0, 5, "ATG")
    assert got["KABHIH_0002"].spacing_dg == 0.3399 and got["KABHIH_0002"].standby_dg == -0.2


def test_excerpts_cover_the_listed_cases() -> None:
    k12, _ = _genes("k12_mg1655")
    nis, nreps = _genes("nissle1917")
    types = {g.feature_type for g in k12} | {g.feature_type for g in nis}
    assert types == {"cds", "sorf", "pseudogene"}
    assert {g.start_codon for g in k12} == {"ATG", "GTG", "TTG"}
    plasmids = {c for c in nreps if c != "CP171242.1"}
    assert {g.interval.contig for g in nis} >= plasmids and len(plasmids) == 3
    last = max(k12, key=lambda g: g.interval.end)
    assert last.locus_tag == "KOIAIB_04715" and last.interval.end > 4641652 - 50


# ---------------------------------------------------------------- full runs ($GIAP_ROOT)


def _full(sample: str) -> Path:
    d = FULL[sample]
    if not (d / "ostir.csv").is_file() or not (d / "ostir_input.csv").is_file():
        pytest.skip(f"pending_sample: full recording not under {d}")
    return d


@pytest.mark.parametrize("sample", sorted(FULL))
def test_full_recording_parses(sample: str) -> None:
    d = _full(sample)
    _check_parsed(d / "ostir_input.csv", d / "ostir.csv", NO_ROW[sample])
    assert len(_dicts(d / "ostir_input.csv")) > 4000


@pytest.mark.parametrize("sample", sorted(FULL))
def test_full_writer_reproduces_the_recorded_input(sample: str, tmp_path: Path) -> None:
    from giap.annotate.bakta import parse_bakta_json
    from giap.annotate.ostir import ostir_row, write_ostir_input
    from giap.model.coords import interval_sort_key

    d = _full(sample)
    js = sorted((GIAP_ROOT / "bakta_full" / sample / "bakta").glob("*.json"))
    if not js:
        pytest.skip(f"pending_sample: recorded Bakta JSON for {sample} missing")
    ann = parse_bakta_json(js[0])
    reps = ann.model.replicons
    rank = {c: k for k, c in enumerate(reps)}
    genes = sorted(ann.model.genes, key=lambda g: interval_sort_key(g.interval, rank, g.locus_tag))
    rows = [
        ostir_row(g, reps[g.interval.contig], ann.sequences[g.interval.contig], ASD) for g in genes
    ]
    out = tmp_path / "ostir_input.csv"
    write_ostir_input(out, rows)
    assert out.read_bytes() == (d / "ostir_input.csv").read_bytes()
