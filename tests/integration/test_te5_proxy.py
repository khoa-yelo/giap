"""TE5 proxy (spec §11.6, M3-15): trypsin and keratin peptides land in the contaminant bucket.

In silico: a giapsim genome stands in for the sample; its six-frame tryptic peptides are the
"simulated" genome peptides (placed and accepted). Fully tryptic digests (7-40 aa) of porcine
trypsin and human keratins 1 and 10 from tests/data/contaminants_mini.faa are added to the table.
Genome hits are decided by an independent oracle (tests/helpers/contaminant_oracle.py). Every
added peptide that does not map to the genome must come out ``contaminant`` with a library ID
whose sequence contains it (I = L); genome peptides never appear.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
from contaminant_oracle import hits_genome, il, six_frames, tryptic_digest

from giap.io.protein_fasta import read_protein_fasta
from giap.mapping.placement_fdr import PER_PEPTIDE_COLUMNS
from giapsim.simulate import SimConfig, simulate

MINI = Path(__file__).resolve().parents[1] / "data" / "contaminants_mini.faa"
CONTAMINANT_SOURCES = ("P00761", "P04264", "P13645")  # trypsin, KRT1, KRT10
SEED = 20261015


def _read_genome(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif line.strip():
            seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


@pytest.fixture(scope="module")
def setup(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    res = simulate(SimConfig(genome_length=30_000), SEED, tmp_path_factory.mktemp("te5"))
    genome = _read_genome(res.genome_fasta)
    frames = [il(f) for seq in genome.values() for f in six_frames(seq)]
    genome_peps: set[str] = set()
    for f in frames:
        for orf in f.split("*"):
            genome_peps.update(tryptic_digest(orf, 7, 30))
    genome_peps = set(sorted(genome_peps)[:400])
    library = read_protein_fasta(MINI)
    contam = sorted(
        {p for r in library if r.id in CONTAMINANT_SOURCES for p in tryptic_digest(r.sequence)}
    )
    return {
        "frames": frames,
        "genome_peps": genome_peps,
        "contam": contam,
        "library": library,
        "lib_seq": {r.id: r.sequence for r in library},
    }


def test_te5_contaminant_digests_land_in_bucket(setup: dict[str, Any]) -> None:
    c = importlib.import_module("giap.mapping.contaminants")
    frames, library = setup["frames"], setup["library"]
    contam: list[str] = setup["contam"]
    assert len(contam) >= 20  # the digests are not trivially small
    assert "VATVSLPR" in contam and "SLNNQFASFIDK" in contam

    all_peps = sorted(set(setup["genome_peps"]) | set(contam))
    mapped = {p for p in all_peps if hits_genome(p, frames)}
    assert setup["genome_peps"] <= mapped

    peptides = pd.DataFrame(
        {
            "peptide": all_peps,
            "protein_ids": pd.Series([[] for _ in all_peps], dtype=object),
            "is_decoy": [False] * len(all_peps),
            "q_value": [0.0] * len(all_peps),
            "is_contaminant": [False] * len(all_peps),
        }
    )
    placed = sorted(mapped)
    placements = pd.DataFrame(
        {
            "peptide": placed,
            "source": ["database"] * len(placed),
            "decoy_type": ["reverse"] * len(placed),
            "length": [len(p) for p in placed],
            "m_target": [0] * len(placed),
            "e_target": [1e-6] * len(placed),
            "placement_q": [0.0] * len(placed),
            "accepted": pd.Series([True] * len(placed), dtype=bool),
            "flags": [""] * len(placed),
        },
        columns=list(PER_PEPTIDE_COLUMNS),
    )
    proteins = pd.DataFrame(
        {
            "protein_id": pd.Series([], dtype=str),
            "sequence": pd.Series([], dtype=str),
            "assigned": pd.Series([], dtype=bool),
            "has_alignment": pd.Series([], dtype=bool),
        }
    )
    unplaced = [p for p in all_peps if p not in mapped]
    pep_matches = c.match_contaminant_peptides(unplaced, library)
    rows = c.unmapped_table(peptides, placements, proteins, pep_matches, {}, ["CON__"])
    by_id = {r.sequence_id: r for r in rows}

    expected = [p for p in contam if p not in mapped]
    assert expected, "every contaminant digest peptide mapped to the genome; test is vacuous"
    for p in expected:
        r = by_id[p]
        assert r.kind == "peptide" and r.label == "contaminant", p
        assert r.reason in ("no_hit", "too_short"), p
        assert r.contaminant_match in CONTAMINANT_SOURCES, p
        assert il(p) in il(setup["lib_seq"][r.contaminant_match]), p
    assert not (set(by_id) & mapped)
    assert set(by_id) == set(unplaced)
