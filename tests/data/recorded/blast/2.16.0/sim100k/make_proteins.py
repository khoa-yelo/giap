"""Regenerate proteins.faa of the recorded tblastn run (see PROVENANCE.md).

Usage: ``PYTHONPATH=src:tests/sim python make_proteins.py OUTDIR``. Writes ``OUTDIR/genome.fa``
(the giapsim SIM100K genome) and ``OUTDIR/proteins.faa``: every CDS/sORF translation > 40 aa of
the simulated gene model (ID = locus tag, terminal ``*`` removed) with an exact number of random
substitutions (never the initial residue) for a target identity drawn uniformly in [0.82, 0.93].
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

from giap.annotate.bakta import parse_bakta_json
from giapsim import SIM100K, SIM100K_SEED, simulate

AA = "ACDEFGHIKLMNPQRSTVWY"
SEED = 13


def diverge(seq: str, identity: float, rng: random.Random) -> str:
    n_sub = round((1.0 - identity) * len(seq))
    out = list(seq)
    for i in rng.sample(range(1, len(seq)), n_sub):
        out[i] = rng.choice([a for a in AA if a != seq[i]])
    return "".join(out)


def main(outdir: Path) -> None:
    sim = simulate(SIM100K, SIM100K_SEED, outdir / "sim")
    (outdir / "genome.fa").write_bytes(sim.genome_fasta.read_bytes())
    model = parse_bakta_json(sim.bakta_json).model
    rng = random.Random(SEED)
    lines: list[str] = []
    for g in model.genes:
        aa = (g.aa or "").rstrip("*")
        if g.feature_type not in ("cds", "sorf") or len(aa) <= 40:
            continue
        seq = diverge(aa, rng.uniform(0.82, 0.93), rng)
        lines.append(f">{g.locus_tag}")
        lines.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
    (outdir / "proteins.faa").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
