"""TE6: Bakta's oriC call on E. coli K-12 sits next to mnmG (gidA) (M2-04 example 10; §11.6).

"On K-12, Bakta's oriC call sits next to mnmG (gidA): mnmG or mioC is the nearest gene, and mnmG
is within 2 kb." Runs on the recorded full-DB JSON ``tests/data/recorded/bakta/1.12.2-full/k12/``
(active since M2-10; the chromosome is circular).
Recorded JSON may lack ``sequences[].nt`` (M2 lead decision 4); a genome FASTA in the same
directory, or else ``tests/data/genomes/k12_mg1655.fna.gz``, then supplies the nucleotides.
"""

from __future__ import annotations

import gzip
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
K12 = DATA / "recorded" / "bakta" / "1.12.2-full" / "k12"
K12_GENOME = DATA / "genomes" / "k12_mg1655.fna.gz"


def _one(*patterns: str) -> Path | None:
    hits = sorted(p for pat in patterns for p in K12.glob(pat))
    return hits[0] if hits else None


def _fasta(path: Path) -> dict[str, str]:
    raw = path.read_bytes()
    text = (gzip.decompress(raw) if path.suffix == ".gz" else raw).decode()
    out: dict[str, list[str]] = {}
    cur = ""
    for line in text.splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif line:
            out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def test_te6_oric_next_to_mnmg() -> None:
    from giap.annotate.bakta import parse_bakta_json
    from giap.annotate.bakta_layers import nearest_genes

    js = _one("*.json", "*.json.gz")
    assert js is not None, f"no Bakta JSON in {K12}"
    fa = _one("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz")
    if fa is None and K12_GENOME.is_file():
        fa = K12_GENOME
    ann = parse_bakta_json(js, genome=_fasta(fa) if fa else None)
    assert ann.db_type == "full"
    assert ann.model.replicons["U00096.3"].topology == "circular"
    orics = [f for f in ann.model.features if f.type == "oric"]
    assert orics, "Bakta called no oriC on K-12"
    genes = list(ann.model.genes)
    gene_name = {g.locus_tag: g.gene for g in genes}
    for oric in orics:
        rep = ann.model.replicons[oric.interval.contig]
        near = nearest_genes(oric, genes, rep, n=20)
        assert gene_name[near[0][0]] in ("mnmG", "mioC"), near[:3]
        mnmg = [d for loc, d in near if gene_name[loc] == "mnmG"]
        assert mnmg and mnmg[0] <= 2_000, near[:5]
