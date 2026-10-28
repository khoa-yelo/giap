"""M11-16c review item M2: ``gippasim.real_genome`` at genetic code 4 (docs/specs/
M11-16-translation-tables-4-25.md item 8: "``real_genome.py`` (peptides from a real genome)
gains the same ``genetic_code`` argument").

Interface assumed: ``RealGenomeConfig.genetic_code: int = 11`` (implementer, e7cc97f2);
``provenance.json`` has ``"genetic_code": t`` only when ``t != 11`` (lead decision on #474).

The fixture genome is a small gippasim genome simulated at code 4 (Bakta-1.12-style JSON from
the v2 writer), run through ``simulate_on_genome(..., check_hashes=False, tu_set=None)`` (gap
rule TUs, no planted islands). Checks: the sample strain's correct peptides are translations of
the sample genome under the table (initiator M over GTG/TTG excepted) with no ``*``; novel ORFs
end in TAA/TAG and their ``aa`` is the translation; alt-site residues are read under the table.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from helpers.m11_16c import STARTS, STOPS_4_25, codon_at, codons, fasta, translate, tsv
from helpers.m11_16c import slice_ as _slice

SEED = 3


def _cfg(table: int) -> Any:
    from gippasim import scenarios

    base = scenarios.k12_evidence()
    return replace(
        base,
        n_novel_per_class=5,
        genetic_code=table,
        expression=replace(base.expression, islands=()),
    )


@pytest.fixture(scope="module")
def genome(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from gippasim.simulate import SimConfig, simulate

    d = tmp_path_factory.mktemp("code4_genome")
    cfg = SimConfig(genetic_code=4, genome_length=100_000, write_precomputed=True,
                    tga_sense_share=0.9)  # fmt: skip
    return simulate(cfg, 1, d)


def _run(genome: Any, table: int, out: Path) -> Any:
    from gippasim.real_genome import simulate_on_genome

    return simulate_on_genome(
        _cfg(table), SEED, out, genome_fasta=genome.genome_fasta, bakta_json=genome.bakta_json,
        tu_set=None, check_hashes=False,
    )  # fmt: skip


def test_real_genome_config_genetic_code() -> None:
    from gippasim.real_genome import RealGenomeConfig

    assert RealGenomeConfig().genetic_code == 11
    assert RealGenomeConfig(genetic_code=4).genetic_code == 4
    for bad in (0, 1, 12, 26):
        with pytest.raises(ValueError):
            RealGenomeConfig(genetic_code=bad)


def test_real_genome_code4_outputs_follow_table(genome: Any, tmp_path: Path) -> None:
    from gippasim.codons import code

    t = code(4)
    res = _run(genome, 4, tmp_path / "o")
    out = res.outdir
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["genetic_code"] == 4
    sample = fasta(out / "truth" / "sample_genome.fna")
    ref = fasta(genome.genome_fasta)

    truth = tsv(out / "truth" / "peptides.tsv")
    correct = [r for r in truth if r["is_wrong"] == "false" and r["contig"] in sample]
    assert correct
    n_tga = 0
    for r in correct:
        pep = r["peptide"]
        assert "*" not in pep
        cod = codons(_slice(sample[r["contig"]], int(r["start"]) - 1, int(r["end"]), r["strand"]))
        assert len(cod) == len(pep), r
        special = False
        for i, (c, a) in enumerate(zip(cod, pep, strict=True)):
            if t[c] == a or {t[c], a} == {"I", "L"}:
                continue
            assert i == 0 and a == "M" and c in STARTS, (r, i, c, a)
            special = True
        if not special:
            assert int(r["n_loci_genome"]) >= 1, r
        n_tga += "TGA" in cod and r["locus_tag"].startswith("GSIM_")
    assert n_tga > 0, "no annotated-gene peptide over a TGA codon"

    novel = tsv(out / "truth" / "novel_orfs.tsv")
    assert novel
    for r in novel:
        nt = _slice(sample[r["contig"]], int(r["start"]) - 1, int(r["end"]), r["strand"])
        cod = codons(nt)
        assert cod[-1] in STOPS_4_25, r
        assert not STOPS_4_25.intersection(cod[1:-1]), r
        assert r["aa"] == translate(nt, 4), r

    alt = tsv(out / "truth" / "alt_sites.tsv")
    assert alt
    for r in alt:
        c0 = int(r["codon_start"]) - 1
        sc = codon_at(sample[r["contig"]], c0, r["strand"])
        rc = codon_at(ref[r["contig"]], c0, r["strand"])
        assert r["alt_aa"] == t[sc] and r["ref_aa"] == t[rc], (r, sc, rc)


def test_real_genome_provenance_omits_code_11(genome: Any, tmp_path: Path) -> None:
    res = _run(genome, 11, tmp_path / "o")
    prov = json.loads((res.outdir / "provenance.json").read_text())
    assert "genetic_code" not in prov
