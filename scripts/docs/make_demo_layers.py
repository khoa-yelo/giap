"""Make the ``gippa demo`` bundle under src/gippa/data/demo/ from gippasim (M11-15 item 2).

Usage: python scripts/docs/make_demo_layers.py [--check]

Simulates SIM100K (seed SIM100K_SEED) and writes, gzipped without name or mtime, its Bakta-like
JSON as ``layers/bakta.json.gz`` and every other ``annotation.precomputed`` layer as
``layers/<key>.tsv.gz`` (operon pairs, promoters, Tx rates, terminators, RBS, regions and the
novel-ORF RBS layer ``rbs_novel``, ADR-0024). It also copies the SIM100K genome
(tests/data/genomes/sim100k.fna.gz) and the tutorial peptide table
(docs/tutorials/data/sim100k_peptides.tsv.gz, written by make_demo_peptides.py) byte for byte.
The simulated genome must equal the committed one. ``--check`` writes nothing and exits 1 when a
committed file differs (layers compared uncompressed, their gzip mtime must be 0). No network.
``PROVENANCE.md`` is written by hand.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / "src"), str(REPO / "tests" / "sim")]

from gippasim import SIM100K, SIM100K_SEED, simulate  # noqa: E402

DEMO = REPO / "src" / "gippa" / "data" / "demo"
GENOME = REPO / "tests" / "data" / "genomes" / "sim100k.fna.gz"
PEPTIDES = REPO / "docs" / "tutorials" / "data" / "sim100k_peptides.tsv.gz"
LAYER_KEYS = ("operon_pairs", "promoters", "tx_rates", "terminators", "rbs", "regions",
              "rbs_novel")  # fmt: skip
MAX_BYTES = 1024 * 1024


def generate() -> dict[str, bytes]:
    """Uncompressed layer bytes by path relative to the demo directory."""
    with tempfile.TemporaryDirectory(prefix="gippa-demo-") as tmp:
        res = simulate(SIM100K, SIM100K_SEED, Path(tmp) / "sim")
        if res.genome_fasta.read_bytes() != gzip.decompress(GENOME.read_bytes()):
            raise SystemExit(f"simulated genome differs from {GENOME.relative_to(REPO)}")
        pre = dict(res.precomputed)
        pre["rbs_novel"] = res.bakta_json.parent / "precomputed" / "rbs_novel.tsv"
        out = {"layers/bakta.json.gz": res.bakta_json.read_bytes()}
        out |= {f"layers/{k}.tsv.gz": Path(pre[k]).read_bytes() for k in LAYER_KEYS}
        return out


def _same_layer(path: Path, want: bytes) -> bool:
    if not path.is_file():
        return False
    data = path.read_bytes()
    return data[4:8] == b"\x00\x00\x00\x00" and gzip.decompress(data) == want


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--check", action="store_true", help="compare with the committed bundle")
    args = ap.parse_args(argv)
    layers = generate()
    copies = {
        "sim100k.fna.gz": GENOME.read_bytes(),
        "sim100k_peptides.tsv.gz": PEPTIDES.read_bytes(),
    }
    if args.check:
        bad = [n for n, b in layers.items() if not _same_layer(DEMO / n, b)]
        bad += [n for n, b in copies.items() if not (DEMO / n).is_file()
                or (DEMO / n).read_bytes() != b]  # fmt: skip
        size = sum(p.stat().st_size for p in DEMO.rglob("*") if p.is_file())
        for name in bad:
            print(f"{(DEMO / name).relative_to(REPO)}: DIFFERS")
        if size > MAX_BYTES:
            print(f"{DEMO.relative_to(REPO)}: {size} bytes (limit {MAX_BYTES})")
        print(f"{DEMO.relative_to(REPO)}: {'up to date' if not bad else 'DIFFERS'}")
        return 0 if not bad and size <= MAX_BYTES else 1
    (DEMO / "layers").mkdir(parents=True, exist_ok=True)
    for name, data in layers.items():
        (DEMO / name).write_bytes(gzip.compress(data, compresslevel=9, mtime=0))
    for name, data in copies.items():
        (DEMO / name).write_bytes(data)
    size = sum(p.stat().st_size for p in DEMO.rglob("*") if p.is_file())
    print(f"wrote {DEMO.relative_to(REPO)}: {len(layers) + len(copies)} files, {size} bytes")
    return 0 if size <= MAX_BYTES else 1


if __name__ == "__main__":
    raise SystemExit(main())
