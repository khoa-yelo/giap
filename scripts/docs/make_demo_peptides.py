"""Make the quickstart demo peptide table from giapsim (backlog 2026-10-25: quickstart demo table).

Usage: python scripts/docs/make_demo_peptides.py [--check]

Simulates SIM100K (seed SIM100K_SEED) with simulated proteomics in one run and writes its
``peptides.tsv``, gzipped without name or mtime, to docs/tutorials/data/sim100k_peptides.tsv.gz.
The simulated genome must equal the committed tests/data/genomes/sim100k.fna.gz, so the table and
that genome belong together. ``--check`` writes nothing and exits 1 when the committed table
differs from what giapsim generates now (compared uncompressed). No network access.
"""

from __future__ import annotations

import argparse
import dataclasses
import gzip
import hashlib
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / "src"), str(REPO / "tests" / "sim")]

from giapsim import SIM100K, SIM100K_SEED, simulate  # noqa: E402
from giapsim.proteomics import ProteomicsConfig  # noqa: E402

GENOME = REPO / "tests" / "data" / "genomes" / "sim100k.fna.gz"
OUT = REPO / "docs" / "tutorials" / "data" / "sim100k_peptides.tsv.gz"
CONFIG = dataclasses.replace(SIM100K, proteomics=ProteomicsConfig(n_runs=1))
SEED = SIM100K_SEED
MAX_BYTES = 100_000


def generate() -> bytes:
    """The uncompressed peptide table; fails if the simulated genome is not the committed one."""
    with tempfile.TemporaryDirectory(prefix="giap-demo-") as tmp:
        res = simulate(CONFIG, SEED, tmp)
        if res.genome_fasta.read_bytes() != gzip.decompress(GENOME.read_bytes()):
            raise SystemExit(f"simulated genome differs from {GENOME.relative_to(REPO)}")
        assert res.peptides is not None
        return res.peptides.read_bytes()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="compare with the committed table")
    args = ap.parse_args(argv)
    table = generate()
    if args.check:
        if OUT.is_file() and OUT.stat().st_size > MAX_BYTES:
            raise SystemExit(f"demo table is {OUT.stat().st_size} bytes (limit {MAX_BYTES})")
        same = OUT.is_file() and gzip.decompress(OUT.read_bytes()) == table
        print(f"{OUT.relative_to(REPO)}: {'up to date' if same else 'DIFFERS'}")
        return 0 if same else 1
    packed = gzip.compress(table, compresslevel=9, mtime=0)
    if len(packed) > MAX_BYTES:
        raise SystemExit(f"demo table is {len(packed)} bytes (limit {MAX_BYTES})")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(packed)
    n_rows = table.count(b"\n") - 1
    print(f"wrote {OUT.relative_to(REPO)}: {n_rows} rows, {len(packed)} bytes")
    print(f"sha256 (tsv)    {hashlib.sha256(table).hexdigest()}")
    print(f"sha256 (tsv.gz) {hashlib.sha256(packed).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
