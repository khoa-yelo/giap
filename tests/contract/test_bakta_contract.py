"""Contract test: Bakta 1.12.2 JSON parser on recorded runs (M2-03; GIAP_SPEC.md §5.4, §7.1.1).

Active since M2-10 (no longer ``pending_sample``): the four directories
``tests/data/recorded/bakta/1.12.2-light/sim100k/`` and ``1.12.2-full/{k12,nissle,sim100k}/`` must
exist. Each holds the Bakta JSON (``*.json`` or ``*.json.gz``), GFF3 and FAA. Recorded JSON may
have ``sequences[].nt`` dropped (lead decision M2 Q4, R11 size limit); then a genome FASTA
(``*.fa``, ``*.fasta``, ``*.fna``, optionally gzipped) in the same directory, or else the committed
genome ``tests/data/genomes/<genome>.fna.gz`` (§11.4), supplies the nucleotides.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

RECORDED = Path(__file__).resolve().parents[1] / "data" / "recorded" / "bakta"
GENOMES = Path(__file__).resolve().parents[1] / "data" / "genomes"
GENOME_FILE = {
    "k12": "k12_mg1655.fna.gz",
    "nissle": "nissle1917.fna.gz",
    "sim100k": "sim100k.fna.gz",
}
EXPECTED = ["1.12.2-light/sim100k", "1.12.2-full/k12", "1.12.2-full/nissle", "1.12.2-full/sim100k"]
SAMPLES = [RECORDED / d for d in EXPECTED]


def _text(path: Path) -> str:
    raw = path.read_bytes()
    return (gzip.decompress(raw) if path.suffix == ".gz" else raw).decode()


def _one(d: Path, *patterns: str) -> Path | None:
    hits = sorted(p for pat in patterns for p in d.glob(pat))
    return hits[0] if hits else None


def _fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur = ""
    for line in _text(path).splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif line:
            out[cur].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}


@pytest.mark.parametrize("sample", SAMPLES, ids=[str(p.relative_to(RECORDED)) for p in SAMPLES])
def test_recorded_bakta_json_parses(sample: Path) -> None:
    from giap.annotate.bakta import parse_bakta_json

    js = _one(sample, "*.json", "*.json.gz")
    gff = _one(sample, "*.gff3", "*.gff3.gz")
    faa = _one(sample, "*.faa", "*.faa.gz")
    assert sample.is_dir(), f"missing {sample} (M2-10 wave A)"
    assert js and gff and faa, f"{sample}: need JSON, GFF3 and FAA"
    fa = _one(sample, "*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz")
    if fa is None and (GENOMES / GENOME_FILE[sample.name]).is_file():
        fa = GENOMES / GENOME_FILE[sample.name]
    genome = _fasta(fa) if fa else None

    ann = parse_bakta_json(js, genome=genome)
    assert ann.bakta_version == "1.12.2"
    assert ann.db_type == ("light" if "light" in sample.parent.name else "full")
    assert ann.db_version

    # one gene per GFF3 CDS locus tag (a split origin-crossing CDS or a pseudogene row counts once)
    gff_loci: set[str] = set()
    for line in _text(gff).splitlines():
        if line.startswith("##FASTA"):
            break
        cols = line.split("\t")
        if not line.startswith("#") and len(cols) == 9 and cols[2] in ("CDS", "pseudogene"):
            attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
            gff_loci.add(attrs["locus_tag"])
    genes = ann.model.genes
    assert len(genes) == len(gff_loci)
    assert {g.locus_tag for g in genes} == gff_loci

    loci = {g.locus_tag for g in genes}
    faa_ids = {ln[1:].split()[0] for ln in _text(faa).splitlines() if ln.startswith(">")}
    assert faa_ids <= loci

    reps = ann.model.replicons
    for g in genes:
        if g.feature_type != "pseudogene":
            assert g.aa, g.locus_tag
        assert g.nt and len(g.nt) == g.interval.length_bp, g.locus_tag
        rep = reps[g.interval.contig]
        if g.interval.end > rep.length:
            assert rep.topology == "circular", g.locus_tag
