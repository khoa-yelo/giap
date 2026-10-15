"""Regenerate proteins.faa of the recorded blastp run (see PROVENANCE.md).

Usage: ``PYTHONPATH=src python make_proteins.py LIBRARY OUT`` with LIBRARY =
``tests/data/contaminants_mini.faa``. Query IDs encode the truth (tests/contract/
test_blastp_contract.py): per library entry two ``HIT_<ACC>_<k>`` (exact substitutions to an
identity drawn uniformly in [0.93, 0.99], full length) and one ``LOWCOV_<ACC>_<k>`` (the first
40 % of the entry fused to 60 % random residues), plus four ``RAND_<k>`` random sequences.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

from giap.io.protein_fasta import read_protein_fasta

AA = "ACDEFGHIKLMNPQRSTVWY"
SEED = 15


def diverge(seq: str, identity: float, rng: random.Random) -> str:
    out = list(seq)
    for i in rng.sample(range(len(seq)), round((1.0 - identity) * len(seq))):
        out[i] = rng.choice([a for a in AA if a != seq[i]])
    return "".join(out)


def main(library: Path, out: Path) -> None:
    rng = random.Random(SEED)
    recs: list[tuple[str, str]] = []
    for r in sorted(read_protein_fasta(library), key=lambda r: r.id):
        for k in (1, 2):
            recs.append((f"HIT_{r.id}_{k}", diverge(r.sequence, rng.uniform(0.93, 0.99), rng)))
        n_keep = int(0.4 * len(r.sequence))
        tail = "".join(rng.choice(AA) for _ in range(len(r.sequence) - n_keep))
        recs.append((f"LOWCOV_{r.id}_1", r.sequence[:n_keep] + tail))
    for k in range(1, 5):
        recs.append(
            (f"RAND_{k}", "M" + "".join(rng.choice(AA) for _ in range(rng.randint(80, 400))))
        )
    lines: list[str] = []
    for rid, seq in recs:
        lines.append(f">{rid}")
        lines.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
    out.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
