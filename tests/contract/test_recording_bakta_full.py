"""Full-DB Bakta fixtures with circular replicons (M2-10 G17; GIAP_SPEC.md §11.4, R11).

G17: "`scripts/slurm/record_bakta_full.sbatch` (16 CPU, 64 GB) runs Bakta 1.12.2 with the full DB
v6.0 on: K-12 MG1655 (GCA_000005845.2); Nissle 1917 (GCA_043228065.1); the 100-kb synthetic
genome. … commits the gzipped JSON, GFF3 and GBFF … under `tests/data/recorded/bakta/1.12.2-full/`."
M2-10: each genome is preflight-normalised first and Bakta gets preflight's ``replicons.tsv``
(``--replicons``); ``contig_map.tsv`` is kept beside the fixture.

Backlog (M3-23 review): the first recordings (2026-10-14) were made from NCBI FASTA without a
circular flag, so preflight wrote ``l`` for every replicon and Bakta annotated the K-12 and Nissle
chromosomes and Nissle's plasmids as linear ``contig``s. These tests pin the re-recording.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import pytest

FULL = Path(__file__).resolve().parents[1] / "data" / "recorded" / "bakta" / "1.12.2-full"
MAX_FILE = 5 * 1024 * 1024

# sample -> {sequence id: (Bakta type, length)}; every replicon circular
EXPECTED: dict[str, dict[str, tuple[str, int]]] = {
    "k12": {"U00096.3": ("chromosome", 4_641_652)},
    "nissle": {
        "CP171242.1": ("chromosome", 5_051_730),
        "CP171243.1": ("plasmid", 4_523),
        "CP171244.1": ("plasmid", 3_173),
        "CP171245.1": ("plasmid", 5_514),
    },
    "sim100k": {"chr1": ("chromosome", 100_000), "plasmid_1": ("plasmid", 8_000)},
}
SAMPLES = sorted(EXPECTED)


def _one(d: Path, suffix: str) -> Path:
    hits = sorted(d.glob(f"*{suffix}"))
    assert len(hits) == 1, f"{d}: expected one *{suffix}, found {[h.name for h in hits]}"
    return hits[0]


def _json(sample: str) -> dict[str, Any]:
    loaded: dict[str, Any] = json.loads(
        gzip.decompress(_one(FULL / sample, ".json.gz").read_bytes())
    )
    return loaded


@pytest.mark.parametrize("sample", SAMPLES)
def test_gzipped_json_gff3_gbff_committed(sample: str) -> None:
    d = FULL / sample
    assert d.is_dir(), f"missing {d} (G17)"
    for suffix in (".json.gz", ".gff3.gz", ".gbff.gz"):
        f = _one(d, suffix)
        raw = f.read_bytes()
        assert raw[:2] == b"\x1f\x8b" and raw[4:8] == b"\0\0\0\0", f"{f.name}: not `gzip -n` output"
        assert len(raw) <= MAX_FILE, f"{f.name} is {len(raw) / 2**20:.1f} MB > 5 MB (R11)"
    assert (
        not list(d.glob("*.json")) and not list(d.glob("*.gff3")) and not list(d.glob("*.gbff"))
    ), "R11: text fixtures are gzipped"


@pytest.mark.parametrize("sample", SAMPLES)
def test_json_is_bakta_1_12_2_full_db_6_0(sample: str) -> None:
    js = _json(sample)
    assert js["version"]["bakta"] == "1.12.2"
    assert js["version"]["db"] == {"version": "6.0", "type": "full"}


@pytest.mark.parametrize("sample", SAMPLES)
def test_json_replicons_are_circular_with_preflight_types(sample: str) -> None:
    seqs = {s["id"]: s for s in _json(sample)["sequences"]}
    assert set(seqs) == set(EXPECTED[sample])
    for sid, (kind, length) in EXPECTED[sample].items():
        s = seqs[sid]
        assert (s["type"], s["topology"], s["length"]) == (kind, "circular", length), sid


@pytest.mark.parametrize("sample", SAMPLES)
def test_gbff_locus_lines_are_circular(sample: str) -> None:
    text = gzip.decompress(_one(FULL / sample, ".gbff.gz").read_bytes()).decode()
    loci = [ln.split() for ln in text.splitlines() if ln.startswith("LOCUS")]
    assert {f[1] for f in loci} == set(EXPECTED[sample])
    assert all("circular" in f and "linear" not in f for f in loci), loci


@pytest.mark.parametrize("sample", SAMPLES)
def test_gff3_sequence_regions_cover_every_replicon(sample: str) -> None:
    text = gzip.decompress(_one(FULL / sample, ".gff3.gz").read_bytes()).decode()
    regions = {
        ln.split()[1]: int(ln.split()[3])
        for ln in text.splitlines()
        if ln.startswith("##sequence-region")
    }
    assert regions == {sid: length for sid, (_, length) in EXPECTED[sample].items()}


@pytest.mark.parametrize("sample", SAMPLES)
def test_preflight_tables_kept_and_command_uses_replicons(sample: str) -> None:
    d = FULL / sample
    reps = d / "replicons.tsv"
    cmap = d / "contig_map.tsv"
    assert reps.is_file() and cmap.is_file(), (
        f"{d}: keep preflight's replicons.tsv and contig_map.tsv"
    )
    rows = [ln.split("\t") for ln in reps.read_text().splitlines() if ln and not ln.startswith("#")]
    table = {r[0]: (r[2], r[3]) for r in rows}
    want_type = {"chromosome": "chromosome", "plasmid": "plasmid"}
    assert table == {sid: (want_type[k], "c") for sid, (k, _) in EXPECTED[sample].items()}
    prov = (d / "PROVENANCE.md").read_text()
    assert "--replicons" in prov and "replicons.tsv" in prov, "recorded command lacks --replicons"


@pytest.mark.parametrize("sample", SAMPLES)
def test_stripped_json_keys_are_documented(sample: str) -> None:
    """Q1 (lead decision 4): `sequences[].nt` may be stripped to meet R11; PROVENANCE says so and
    the genome is committed under tests/data/genomes/."""
    js = _json(sample)
    if any("nt" not in s for s in js["sequences"]):
        prov = (FULL / sample / "PROVENANCE.md").read_text()
        assert "sequences.nt" in prov or "sequences[].nt" in prov


def test_k12_has_an_oric_call() -> None:
    """TE6 needs Bakta's oriC call on K-12 (full DB, DoriC 12)."""
    assert any(f["type"] == "oriC" for f in _json("k12")["features"])
