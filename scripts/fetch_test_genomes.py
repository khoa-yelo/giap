#!/usr/bin/env python3
"""Write the committed test genomes tests/data/genomes/*.fna.gz (M2-10; GIPPA_SPEC.md §11.4, R11).

K-12 MG1655 (GCA_000005845.2) and Nissle 1917 (GCA_043228065.1) come from NCBI Datasets
(api.ncbi.nlm.nih.gov); headers are kept as NCBI gives them and ``[topology=circular]`` is appended
to every replicon of these complete assemblies (lead decision 1). With ``--compare DIR`` the
sequences must equal those of the NCBI FTP downloads of ``fetch_resources.sbatch``
(``DIR/<name>/*_genomic.fna.gz``). ``sim100k.fna.gz`` is the gippasim SIM100K/SIM100K_SEED genome.
Files are gzipped reproducibly (no name, mtime 0); PROVENANCE.md lists their sha256 sums.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import socket
import tempfile
import urllib.request
import zipfile
from datetime import UTC, datetime
from pathlib import Path

API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"
GENOMES = {"k12_mg1655": "GCA_000005845.2", "nissle1917": "GCA_043228065.1"}
TAG = "[topology=circular]"


def gz(data: bytes) -> bytes:
    """``gzip -n -9`` equivalent (no name, mtime 0): identical input gives identical bytes."""
    return gzip.compress(data, compresslevel=9, mtime=0)


def records(text: str) -> list[tuple[str, str]]:
    """(header, upper-case sequence) for every FASTA record of ``text``."""
    out: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        if line.startswith(">"):
            out.append((line[1:].strip(), []))
        elif line.strip():
            out[-1][1].append(line.strip().upper())
    return [(h, "".join(s)) for h, s in out]


def tag_circular(text: str) -> str:
    """Append ``[topology=circular]`` to every header that lacks it; sequence lines unchanged."""
    lines = [
        f"{ln.rstrip()} {TAG}" if ln.startswith(">") and TAG not in ln else ln
        for ln in text.splitlines()
    ]
    return "\n".join(lines) + "\n"


def download(acc: str) -> str:
    """The genomic FASTA of ``acc`` from the NCBI Datasets v2 API."""
    with urllib.request.urlopen(API.format(acc=acc), timeout=300) as r:
        zf = zipfile.ZipFile(io.BytesIO(r.read()))
    names = [n for n in zf.namelist() if n.endswith(("_genomic.fna", ".fna"))]
    if len(names) != 1:
        raise SystemExit(f"{acc}: expected one FASTA in the Datasets zip, got {names}")
    return zf.read(names[0]).decode()


def seq_sha(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=Path("tests/data/genomes"))
    ap.add_argument("--compare", type=Path, help="dir with <name>/*_genomic.fna.gz (NCBI FTP)")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    now = datetime.now(UTC)
    notes: list[str] = []
    for name, acc in GENOMES.items():
        text = download(acc)
        recs = records(text)
        if args.compare:
            (ftp,) = sorted((args.compare / name).glob("*_genomic.fna.gz"))
            want = records(gzip.decompress(ftp.read_bytes()).decode())
            if [(h.split()[0], s) for h, s in want] != [(h.split()[0], s) for h, s in recs]:
                raise SystemExit(f"{name}: Datasets and FTP sequences differ")
            notes.append(f"- `{name}`: sequences identical to the NCBI FTP file `{ftp.name}`.")
        (args.out / f"{name}.fna.gz").write_bytes(gz(tag_circular(text).encode()))
        for h, s in recs:
            notes.append(f"  - `{h.split()[0]}` {len(s):,} bp, sequence sha256 `{seq_sha(s)}`")

    from gippasim import SIM100K, SIM100K_SEED, simulate

    with tempfile.TemporaryDirectory() as tmp:
        sim = simulate(SIM100K, SIM100K_SEED, Path(tmp) / "sim").genome_fasta.read_bytes()
    (args.out / "sim100k.fna.gz").write_bytes(gz(sim))

    sums = "\n".join(
        f"{hashlib.sha256((args.out / f).read_bytes()).hexdigest()}  {f}"
        for f in ("k12_mg1655.fna.gz", "nissle1917.fna.gz", "sim100k.fna.gz")
    )
    acc_lines = "\n".join(f"- `{n}.fna.gz`: {a}" for n, a in GENOMES.items())
    (args.out / "PROVENANCE.md").write_text(
        "# Test genomes (M2-10; GIPPA_SPEC.md §11.4)\n\n"
        f"- Command: `python scripts/fetch_test_genomes.py --out {args.out}"
        f"{' --compare ' + str(args.compare) if args.compare else ''}`\n"
        f"- Date: {now:%Y-%m-%dT%H:%M:%SZ}\n"
        f"- Host: {socket.gethostname()}\n"
        f"- Source: NCBI Datasets v2 (`{API.split('/datasets')[0]}`), assemblies:\n{acc_lines}\n"
        "- Cross-check: NCBI FTP downloads of `scripts/slurm/fetch_resources.sbatch`, fetched\n"
        "  20261014 (2026-10-14) into `$GIPPA_DATA_DIR/genomes/<name>/` (md5 verified).\n"
        f"- Headers are NCBI's with ` {TAG}` appended to every replicon (chromosomes and plasmids\n"
        "  of these complete assemblies are circular; lead decision 1). Sequences are unchanged.\n"
        "- `sim100k.fna.gz`: gippasim `simulate(SIM100K, SIM100K_SEED)` genome FASTA, verbatim.\n"
        "- Gzipped without name or mtime (`gzip -n -9` equivalent).\n\n"
        "## Sequences (sha256 of the upper-case sequence, no header or newlines)\n\n"
        + "\n".join(notes)
        + f"\n\n## Files (sha256sum)\n\n```\n{sums}\n```\n"
    )
    print(sums)


if __name__ == "__main__":
    main()
