"""NCBI genome package → topology-tagged FASTA for Bakta (campaign SPEC §6).

The package (``ncbi_dataset/data/{acc}/*_genomic.fna`` and ``genomic.gbff``) is fetched and
extracted by the driver; these helpers pick its files, read each replicon's topology from the
GenBank ``LOCUS`` lines and tag circular FASTA headers with `` [topology=circular]``.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import IO, cast

__all__ = ["CIRCULAR_TAG", "find_genome_files", "gbff_topologies", "genome_zip_path", "tag_fasta"]

CIRCULAR_TAG = b" [topology=circular]"
_GZIP_MAGIC = b"\x1f\x8b"


def genome_zip_path(accession: str) -> str:
    """Datasets v2 path (below the API root) of the FASTA + GBFF package of ``accession``."""
    return (
        f"/genome/accession/{accession}/download"
        "?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GBFF"
    )


def _open_bytes(path: Path) -> IO[bytes]:
    """``path`` opened for binary reading, gunzipped when it starts with the gzip magic."""
    with path.open("rb") as fh:
        gzipped = fh.read(2) == _GZIP_MAGIC
    return cast(IO[bytes], gzip.open(path, "rb")) if gzipped else path.open("rb")


def find_genome_files(paths: Iterable[Path]) -> tuple[Path, Path]:
    """The one ``*_genomic.fna`` and the one ``genomic.gbff`` among the files in ``paths``."""
    files = [Path(p) for p in paths if Path(p).is_file()]
    fna = [p for p in files if p.name.endswith("_genomic.fna")]
    gbff = [p for p in files if p.name == "genomic.gbff"]
    if len(fna) != 1 or len(gbff) != 1:
        raise ValueError(
            f"expected one *_genomic.fna and one genomic.gbff, found {len(fna)} and {len(gbff)}"
        )
    return fna[0], gbff[0]


def gbff_topologies(gbff: Path) -> dict[str, str]:
    """``VERSION`` accession and ``LOCUS`` name → ``circular`` / ``linear`` (default ``linear``)."""
    topologies: dict[str, str] = {}
    topology = "linear"
    with _open_bytes(Path(gbff)) as fh:
        for raw in fh:
            line = raw.decode("utf-8", errors="replace")
            if line.startswith("LOCUS"):
                tokens = line.split()
                topology = "circular" if "circular" in tokens[2:] else "linear"
                if len(tokens) > 1:
                    topologies[tokens[1]] = topology
            elif line.startswith("VERSION"):
                tokens = line.split()
                if len(tokens) > 1:
                    topologies[tokens[1]] = topology
    return topologies


def tag_fasta(src: Path, dst: Path, topologies: Mapping[str, str]) -> dict[str, int]:
    """Copy ``src`` to plain ``dst`` byte for byte; headers with a circular first token tagged."""
    counts = {"circular": 0, "linear": 0, "unknown": 0}
    with _open_bytes(Path(src)) as fin, Path(dst).open("wb") as fout:
        for line in fin:
            if line.startswith(b">"):
                tokens = line[1:].split(maxsplit=1)
                seq_id = tokens[0].decode("utf-8", errors="replace") if tokens else ""
                topology = topologies.get(seq_id)
                if topology not in ("circular", "linear"):
                    counts["unknown"] += 1
                else:
                    counts[topology] += 1
                if topology == "circular":
                    body = line.rstrip(b"\r\n")
                    line = body + CIRCULAR_TAG + line[len(body) :]
            fout.write(line)
    return counts
