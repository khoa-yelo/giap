"""giapsim M0 skeleton: one circular chromosome, genes, Bakta-like JSON (GIAP_SPEC.md §11.3)."""

from __future__ import annotations

import csv
import itertools
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

BASES = "TCAG"
AA_11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: AA_11[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}
STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def translate11(nt: str) -> str:
    """Table-11 translation of a CDS (stop omitted, initiator forced to M)."""
    aa = "".join(CODON_TABLE[nt[i : i + 3]] for i in range(0, len(nt) - 2, 3))
    if aa.endswith("*"):
        aa = aa[:-1]
    return "M" + aa[1:]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def cds_features(data: dict[str, Any]) -> list[dict[str, Any]]:
    return [f for f in data["features"] if f["type"] == "cds"]


def at_frac(s: str) -> float:
    return sum(1 for c in s if c in "AT") / len(s)


# ---------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def small(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim

    out = tmp_path_factory.mktemp("sim50k")
    return giapsim.simulate(giapsim.SimConfig(genome_length=50_000), seed=1, outdir=out)


@pytest.fixture(scope="module")
def big(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim

    out = tmp_path_factory.mktemp("sim1m")
    return giapsim.simulate(giapsim.SimConfig(genome_length=1_000_000), seed=2, outdir=out)


@pytest.fixture(scope="module")
def small_json(small: Any) -> dict[str, Any]:
    data: dict[str, Any] = json.loads(Path(small.bakta_json).read_text())
    return data


@pytest.fixture(scope="module")
def big_json(big: Any) -> dict[str, Any]:
    data: dict[str, Any] = json.loads(Path(big.bakta_json).read_text())
    return data


@pytest.fixture(scope="module")
def small_genome(small: Any) -> str:
    return read_fasta(Path(small.genome_fasta))[0][1]


def read_truth(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- files


def test_default_genome_length_is_50kb() -> None:
    import giapsim

    assert giapsim.SimConfig().genome_length == 50_000


def test_writes_expected_files(small: Any) -> None:
    for attr, name in [
        ("genome_fasta", "genome.fasta"),
        ("bakta_json", "bakta_like.json"),
        ("genes_truth", "genes.tsv"),
    ]:
        p = Path(getattr(small, attr))
        assert p.is_file()
        assert p.name == name
    assert Path(small.genes_truth).parent.name == "truth"


def test_same_seed_gives_identical_files(tmp_path: Path) -> None:
    import giapsim

    cfg = giapsim.SimConfig(genome_length=30_000)
    a = giapsim.simulate(cfg, seed=11, outdir=tmp_path / "a")
    b = giapsim.simulate(cfg, seed=11, outdir=tmp_path / "b")
    for attr in ("genome_fasta", "bakta_json", "genes_truth"):
        assert Path(getattr(a, attr)).read_bytes() == Path(getattr(b, attr)).read_bytes()


def test_different_seed_gives_different_genome(tmp_path: Path) -> None:
    import giapsim

    cfg = giapsim.SimConfig(genome_length=30_000)
    a = giapsim.simulate(cfg, seed=11, outdir=tmp_path / "a")
    b = giapsim.simulate(cfg, seed=12, outdir=tmp_path / "b")
    assert read_fasta(Path(a.genome_fasta))[0][1] != read_fasta(Path(b.genome_fasta))[0][1]


# ---------------------------------------------------------------- genome


def test_single_replicon_of_configured_length(small: Any) -> None:
    records = read_fasta(Path(small.genome_fasta))
    assert len(records) == 1
    assert len(records[0][1]) == 50_000
    assert set(records[0][1]) <= set("ACGT")


def test_fasta_header_marks_circular(small: Any) -> None:
    header = read_fasta(Path(small.genome_fasta))[0][0]
    assert "[topology=circular]" in header


# ---------------------------------------------------------------- CDS invariants


def test_cds_length_multiple_of_three(small_json: dict[str, Any]) -> None:
    for f in cds_features(small_json):
        assert len(f["nt"]) % 3 == 0, f["locus"]
        assert len(f["nt"]) == f["stop"] - f["start"] + 1, f["locus"]


def test_cds_start_codons(small_json: dict[str, Any]) -> None:
    for f in cds_features(small_json):
        assert f["nt"][:3] in STARTS, f["locus"]
        assert f["start_type"] == f["nt"][:3], f["locus"]


def test_cds_ends_with_stop(small_json: dict[str, Any]) -> None:
    for f in cds_features(small_json):
        assert f["nt"][-3:] in STOPS, f["locus"]


def test_cds_no_internal_stop(small_json: dict[str, Any]) -> None:
    for f in cds_features(small_json):
        nt = f["nt"]
        internal = [nt[i : i + 3] for i in range(3, len(nt) - 3, 3)]
        assert not STOPS.intersection(internal), f["locus"]


def test_cds_within_replicon(small_json: dict[str, Any]) -> None:
    lengths = {s["id"]: s["length"] for s in small_json["sequences"]}
    for f in cds_features(small_json):
        assert 1 <= f["start"] <= f["stop"] <= lengths[f["contig"]], f["locus"]


def test_both_strands_present(small_json: dict[str, Any]) -> None:
    strands = {f["strand"] for f in cds_features(small_json)}
    assert strands == {"+", "-"}


def test_nt_matches_genome_on_both_strands(small_json: dict[str, Any], small_genome: str) -> None:
    for f in cds_features(small_json):
        seg = small_genome[f["start"] - 1 : f["stop"]]
        expected = seg if f["strand"] == "+" else revcomp(seg)
        assert f["nt"] == expected, f["locus"]


def test_aa_is_table11_translation(small_json: dict[str, Any]) -> None:
    for f in cds_features(small_json):
        assert f["aa"] == translate11(f["nt"]), f["locus"]
        assert f["aa"].startswith("M")
        assert "*" not in f["aa"]


def test_start_codon_frequencies(big_json: dict[str, Any]) -> None:
    counts = Counter(f["nt"][:3] for f in cds_features(big_json))
    n = sum(counts.values())
    assert n > 500
    assert 0.85 <= counts["ATG"] / n <= 0.95
    assert 0.04 <= counts["GTG"] / n <= 0.12
    assert 0.002 <= counts["TTG"] / n <= 0.05


def test_genes_overlap_at_most_20bp(big_json: dict[str, Any]) -> None:
    feats = sorted(cds_features(big_json), key=lambda f: (f["contig"], f["start"]))
    for a, b in itertools.pairwise(feats):
        if a["contig"] != b["contig"]:
            continue
        overlap = a["stop"] - b["start"] + 1
        assert overlap <= 20, (a["locus"], b["locus"], overlap)


# ---------------------------------------------------------------- Bakta-like JSON


def test_json_top_level(small_json: dict[str, Any], small_genome: str) -> None:
    assert isinstance(small_json["features"], list)
    seqs = small_json["sequences"]
    assert isinstance(seqs, list) and len(seqs) == 1
    s = seqs[0]
    assert s["length"] == 50_000
    assert s["topology"] == "circular"
    assert s["sequence"] == small_genome
    assert isinstance(s["id"], str) and s["id"]


def test_json_cds_fields(small_json: dict[str, Any]) -> None:
    feats = cds_features(small_json)
    assert len(feats) > 10
    seq_ids = {s["id"] for s in small_json["sequences"]}
    for f in feats:
        for key in ("contig", "start", "stop", "strand", "locus", "aa", "nt", "start_type"):
            assert key in f, key
        assert f["contig"] in seq_ids
        assert f["strand"] in {"+", "-"}
        assert isinstance(f["start"], int) and isinstance(f["stop"], int)


def test_json_locus_tags_unique(small_json: dict[str, Any]) -> None:
    loci = [f["locus"] for f in cds_features(small_json)]
    assert len(loci) == len(set(loci))


# ---------------------------------------------------------------- truth / TUs


def test_truth_columns(small: Any) -> None:
    with Path(small.genes_truth).open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
    assert header[:5] == ["locus_tag", "tu_id", "start", "end", "strand"]


def test_truth_matches_json(small: Any, small_json: dict[str, Any]) -> None:
    rows = {r["locus_tag"]: r for r in read_truth(Path(small.genes_truth))}
    feats = cds_features(small_json)
    assert set(rows) == {f["locus"] for f in feats}
    for f in feats:
        r = rows[f["locus"]]
        assert (int(r["start"]), int(r["end"]), r["strand"]) == (f["start"], f["stop"], f["strand"])


def test_tu_genes_share_strand_and_are_adjacent(big: Any) -> None:
    rows = read_truth(Path(big.genes_truth))
    rows.sort(key=lambda r: int(r["start"]))
    order = {r["locus_tag"]: i for i, r in enumerate(rows)}
    by_tu: dict[str, list[dict[str, str]]] = {}
    for r in rows:
        by_tu.setdefault(r["tu_id"], []).append(r)
    for tu, members in by_tu.items():
        assert len({m["strand"] for m in members}) == 1, tu
        idx = sorted(order[m["locus_tag"]] for m in members)
        assert idx == list(range(idx[0], idx[0] + len(idx))), tu


def test_intra_tu_gaps_within_spec(big: Any) -> None:
    rows = sorted(read_truth(Path(big.genes_truth)), key=lambda r: int(r["start"]))
    for a, b in itertools.pairwise(rows):
        if a["tu_id"] == b["tu_id"]:
            gap = int(b["start"]) - int(a["end"]) - 1
            assert -20 <= gap <= 40, (a["locus_tag"], b["locus_tag"], gap)


def test_mean_tu_size_near_1_8(big: Any) -> None:
    rows = read_truth(Path(big.genes_truth))
    sizes = Counter(r["tu_id"] for r in rows)
    mean = sum(sizes.values()) / len(sizes)
    assert 1.5 <= mean <= 2.1


# ---------------------------------------------------------------- composition


def test_intergenic_more_at_rich_than_coding(small_json: dict[str, Any], small_genome: str) -> None:
    coding = bytearray(len(small_genome))
    for f in cds_features(small_json):
        coding[f["start"] - 1 : f["stop"]] = b"\x01" * (f["stop"] - f["start"] + 1)
    cod = "".join(c for c, m in zip(small_genome, coding, strict=True) if m)
    inter = "".join(c for c, m in zip(small_genome, coding, strict=True) if not m)
    assert len(inter) > 500
    assert at_frac(inter) > at_frac(cod) + 0.05


# ---------------------------------------------------------------- speed


def test_100kb_genome_is_fast(tmp_path: Path) -> None:
    import giapsim

    t0 = time.perf_counter()
    giapsim.simulate(giapsim.SimConfig(genome_length=100_000), seed=5, outdir=tmp_path)
    # loose 10x guard; the 2 s budget is enforced in benchmarks/test_giapsim.py
    assert time.perf_counter() - t0 < 20.0
