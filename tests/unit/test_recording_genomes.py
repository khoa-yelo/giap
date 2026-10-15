"""Committed test genomes (M2-10; GIAP_SPEC.md §11.4, R11).

§11.4: "The K-12 and Nissle FASTA files (gzipped, < 2 MB each) are committed under
`tests/data/genomes/`, with a provenance note." M2-10 `scripts/fetch_test_genomes.py` writes
``{k12_mg1655,nissle1917,sim100k}.fna.gz`` (NCBI headers kept, ``[topology=circular]`` added) and
``PROVENANCE.md``. The circular flag is what makes preflight's ``replicons.tsv`` (and therefore
Bakta ``--replicons``) mark the chromosomes and Nissle's plasmids circular; the first full-DB
recordings lacked it and annotated every replicon as linear (backlog, M3-23 review).

Sequence SHA-256s (upper-case sequence, no header/newlines) were computed by the test author from
the NCBI Datasets downloads GCA_000005845.2 and GCA_043228065.1 (``$GIAP_DATA_DIR/genomes``,
fetched 2026-10-14 by ``fetch_resources.sbatch``).
"""

from __future__ import annotations

import gzip
import hashlib
import re
from pathlib import Path

import pytest

from giap.config.model import PreflightConfig
from giap.io.fasta import read_genome_fasta
from giap.preflight.genome import classify_replicons, resolve_topology, summarize_genome

GENOMES = Path(__file__).resolve().parents[1] / "data" / "genomes"
FILES = ("k12_mg1655.fna.gz", "nissle1917.fna.gz", "sim100k.fna.gz")
MAX_REAL_BYTES = 2 * 1024 * 1024  # §11.4 "< 2 MB each"

# id -> (NCBI header text, length, sha256 of the sequence)
NCBI = {
    "k12_mg1655.fna.gz": {
        "U00096.3": (
            "U00096.3 Escherichia coli str. K-12 substr. MG1655, complete genome",
            4_641_652,
            "e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0",
        ),
    },
    "nissle1917.fna.gz": {
        "CP171242.1": (
            "CP171242.1 Escherichia coli Nissle 1917 chromosome, complete genome",
            5_051_730,
            "14d7914bdd9073ccc04beab477fc4b0341d17ce3180be59fb46b9e4ea370d3fb",
        ),
        "CP171243.1": (
            "CP171243.1 Escherichia coli Nissle 1917 plasmid p15A-ILA-dapA, complete sequence",
            4_523,
            "40c11dbad5ab475d2b9f87cccf56190c705c018840ea3eb807064e90bc8d5c34",
        ),
        "CP171244.1": (
            "CP171244.1 Escherichia coli Nissle 1917 plasmid pMUT1, complete sequence",
            3_173,
            "7278bc2cf0fd2a8af41f3d7f3e6751a5d0038651526e2fcc1a133a6c11ae7ee8",
        ),
        "CP171245.1": (
            "CP171245.1 Escherichia coli Nissle 1917 plasmid pMUT2, complete sequence",
            5_514,
            "89b32fa4e38e2127d0c664bc37b3e1c041522fa060995617275d1aa44dc677ce",
        ),
    },
}
SHA_LINE = re.compile(r"^\s*([0-9a-f]{64})\s+\*?(\S+)\s*$", re.MULTILINE)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _replicons(path: Path) -> dict[str, tuple[str, str]]:
    """id -> (topology, kind) exactly as `giap preflight` classifies them."""
    cfg = PreflightConfig()
    summary = summarize_genome(read_genome_fasta(path), cfg.n_run_min_bp)
    reps = classify_replicons(
        summary,
        resolve_topology(summary, "auto"),
        small_bp=cfg.small_replicon_bp,
        chromosome_min_bp=cfg.min_largest_contig_warn,
    )
    return {r.id: (r.topology, r.kind) for r in reps}


@pytest.mark.parametrize("name", FILES)
def test_genome_committed_gzipped_reproducibly(name: str) -> None:
    path = GENOMES / name
    assert path.is_file(), f"missing {path} (M2-10 scripts/fetch_test_genomes.py)"
    raw = path.read_bytes()
    assert raw[:2] == b"\x1f\x8b", f"{name} is not gzip"
    # `gzip -n`: no stored file name (FLG.FNAME) and MTIME 0, so re-fetching gives identical bytes
    assert raw[3] & 0x08 == 0, f"{name}: gzip header stores a file name (use gzip -n)"
    assert raw[4:8] == b"\x00\x00\x00\x00", f"{name}: gzip header stores an mtime (use gzip -n)"
    gzip.decompress(raw)
    if name != "sim100k.fna.gz":
        assert len(raw) < MAX_REAL_BYTES, f"{name} is {len(raw)} bytes >= 2 MB (§11.4)"


@pytest.mark.parametrize("name", sorted(NCBI))
def test_real_genome_is_the_ncbi_assembly_with_circular_flags(name: str) -> None:
    path = GENOMES / name
    assert path.is_file(), f"missing {path}"
    recs = read_genome_fasta(path)
    assert [r.original_id for r in recs] == list(NCBI[name])
    for r in recs:
        header, length, sha = NCBI[name][r.original_id]
        assert r.header.startswith(header), f"NCBI header changed: {r.header!r}"
        assert "[topology=circular]" in r.header, r.header
        assert r.circular_tag
        assert len(r.sequence) == length
        assert _sha(r.sequence.encode()) == sha, f"{r.original_id}: not the NCBI sequence"


def test_preflight_makes_k12_chromosome_circular() -> None:
    assert _replicons(GENOMES / "k12_mg1655.fna.gz") == {"U00096.3": ("circular", "chromosome")}


def test_preflight_lists_nissle_plasmids_as_circular_plasmids() -> None:
    """Edge case: "Nissle has two plasmids (pMUT1, pMUT2); the replicon table must list them as
    plasmids" (the assembly has a third, p15A-ILA-dapA; it is a plasmid too)."""
    assert _replicons(GENOMES / "nissle1917.fna.gz") == {
        "CP171242.1": ("circular", "chromosome"),
        "CP171243.1": ("circular", "plasmid"),
        "CP171244.1": ("circular", "plasmid"),
        "CP171245.1": ("circular", "plasmid"),
    }


def test_sim100k_genome_is_the_simulator_genome(tmp_path: Path) -> None:
    from giapsim import SIM100K, SIM100K_SEED, simulate

    path = GENOMES / "sim100k.fna.gz"
    assert path.is_file(), f"missing {path}"
    want = read_genome_fasta(simulate(SIM100K, SIM100K_SEED, tmp_path / "sim").genome_fasta)
    got = read_genome_fasta(path)
    assert [(r.id, r.sequence, r.circular_tag) for r in got] == [
        (r.id, r.sequence, r.circular_tag) for r in want
    ]
    assert _replicons(path) == {
        "chr1": ("circular", "chromosome"),
        "plasmid_1": ("circular", "plasmid"),
    }


def test_genomes_provenance_is_complete_and_checksums_match() -> None:
    prov = GENOMES / "PROVENANCE.md"
    assert prov.is_file(), f"missing {prov} (§11.4 'with a provenance note')"
    text = prov.read_text()
    for needle in (
        "GCA_000005845.2",
        "GCA_043228065.1",
        "api.ncbi.nlm.nih.gov",
        "SIM100K",
        "20261014",
    ):
        assert needle in text, f"PROVENANCE.md does not mention {needle}"
    assert "[topology=circular]" in text, "PROVENANCE.md must say the circular flag was added"
    assert re.search(r"\b20\d\d-\d\d-\d\d\b", text), "PROVENANCE.md has no date"
    listed = {Path(p).name: h for h, p in SHA_LINE.findall(text)}
    for name in FILES:
        assert name in listed, f"PROVENANCE.md lists no `sha256  {name}` line"
        assert listed[name] == _sha((GENOMES / name).read_bytes()), f"{name}: sha256 mismatch"
