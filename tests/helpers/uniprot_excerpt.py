"""Build the recorded UniProt K-12 excerpt for the M3-22 TE8 CI test (deterministic).

Usage: ``python -P tests/helpers/uniprot_excerpt.py <UP000000625.fasta> <out.fasta>``

Selection: entries sorted by accession, every ``STEP``-th one (starting at index 0), plus the
named ``EXTRA`` accessions (EF-Tu, beta-galactosidase); records are written verbatim (header and
sequence, 60 residues per line), in accession order.
"""

from __future__ import annotations

import sys
from pathlib import Path

STEP = 110
EXTRA = ("P0CE47", "P00722")


def _records(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    head, seq = "", []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if head:
                out.append((head, "".join(seq)))
            head, seq = line, []
        elif line.strip():
            seq.append(line.strip())
    if head:
        out.append((head, "".join(seq)))
    return out


def accession(header: str) -> str:
    return header[1:].split()[0].split("|")[1]


def build(src: Path, dst: Path) -> int:
    recs = sorted(_records(src), key=lambda r: accession(r[0]))
    keep = [r for i, r in enumerate(recs) if i % STEP == 0 or accession(r[0]) in EXTRA]
    lines = []
    for head, seq in keep:
        lines.append(head)
        lines += [seq[i : i + 60] for i in range(0, len(seq), 60)]
    dst.write_text("\n".join(lines) + "\n")
    return len(keep)


if __name__ == "__main__":
    print(build(Path(sys.argv[1]), Path(sys.argv[2])))
