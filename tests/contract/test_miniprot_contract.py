"""Contract test: miniprot 0.18 GFF parser on a recorded run (issue M3-12 example 10; §11.2).

Marked ``pending_sample`` until ``tests/data/recorded/miniprot/0.18/sim100k/`` exists. The
directory holds the native miniprot ``--gff`` output (``*.gff``, optionally gzipped), the query
``proteins.faa`` and ``PROVENANCE.md`` (command, env lock or image digest, date, host). The gene
model is the giapsim 100-kb genome (``SIM100K``, ``SIM100K_SEED``), regenerated here; a Bakta-like
``*.json`` in the directory is used instead when present. Protein IDs are the locus tags.

Example 10: the output parses and >= 99 % of the simulated proteins > 40 aa are assigned
``identical`` to their own gene.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pytest

from giap.config.model import ProteinsConfig

RECORDED = (
    Path(__file__).resolve().parents[1] / "data" / "recorded" / "miniprot" / "0.18" / "sim100k"
)

pytestmark = [
    pytest.mark.pending_sample,
    pytest.mark.skipif(
        not RECORDED.is_dir(), reason=f"pending_sample: no recorded miniprot run in {RECORDED}"
    ),
]


def _one(*patterns: str) -> Path | None:
    hits = sorted(p for pat in patterns for p in RECORDED.glob(pat))
    return hits[0] if hits else None


def _plain(path: Path, tmp: Path) -> Path:
    if path.suffix != ".gz":
        return path
    out = tmp / path.name.removesuffix(".gz")
    out.write_bytes(gzip.decompress(path.read_bytes()))
    return out


@pytest.fixture(scope="module")
def recorded(tmp_path_factory: pytest.TempPathFactory):  # type: ignore[no-untyped-def]
    from giap.annotate.bakta import parse_bakta_json
    from giap.io.protein_fasta import read_protein_fasta

    tmp = tmp_path_factory.mktemp("miniprot_contract")
    gff = _one("*.gff", "*.gff.gz", "*.gff3", "*.gff3.gz")
    faa = _one("proteins.faa", "proteins.faa.gz", "*.faa", "*.faa.gz")
    assert gff and faa, f"{RECORDED}: need the miniprot GFF and proteins.faa"
    js = _one("*.json", "*.json.gz")
    if js is None:
        from giapsim import SIM100K, SIM100K_SEED, simulate

        js = simulate(SIM100K, SIM100K_SEED, tmp / "sim").bakta_json
    model = parse_bakta_json(_plain(js, tmp)).model
    proteins = {r.id: r.sequence for r in read_protein_fasta(_plain(faa, tmp))}
    return _plain(gff, tmp), model, proteins


def test_provenance_recorded() -> None:
    prov = RECORDED / "PROVENANCE.md"
    assert prov.is_file(), "PROVENANCE.md (command, env lock or digest, date, host) missing"
    text = prov.read_text()
    assert "miniprot" in text and "0.18" in text
    assert re.search(r"-T\s+11", text), "indexing command with -T 11 not recorded"
    assert "--gff" in text
    assert re.search(r"\d{4}-\d{2}-\d{2}", text), "date missing"


def test_recorded_output_parses(recorded) -> None:  # type: ignore[no-untyped-def]
    from giap.mapping.miniprot import parse_miniprot_gff

    gff, model, proteins = recorded
    lens = {k: len(v) for k, v in proteins.items()}
    alns = parse_miniprot_gff(gff, model.replicons, lens)
    assert alns, "no alignments parsed from the recorded GFF"
    n_mrna = sum(
        1
        for line in gff.read_text().splitlines()
        if not line.startswith("#") and line.count("\t") >= 8 and line.split("\t")[2] == "mRNA"
    )
    assert len(alns) == n_mrna
    for a in alns:
        rep = model.replicons[a.interval.contig]
        assert 0 <= a.interval.start < a.interval.end <= rep.length
        assert 0 <= a.q_start < a.q_end <= a.protein_len_aa == lens[a.protein_id]
        assert 0.0 <= a.identity <= 1.0
        assert a.rank >= 1
        assert a.route == "miniprot"


def test_example10_own_gene_identical(recorded) -> None:  # type: ignore[no-untyped-def]
    from giap.mapping.miniprot import parse_miniprot_gff
    from giap.mapping.protein_assign import assign_alignments

    gff, model, proteins = recorded
    lens = {k: len(v) for k, v in proteins.items()}
    df, _left = assign_alignments(
        parse_miniprot_gff(gff, model.replicons, lens), model, ProteinsConfig()
    )
    tags = {g.locus_tag for g in model.genes if g.feature_type in ("cds", "sorf")}
    eligible = {pid for pid, seq in proteins.items() if len(seq.rstrip("*")) > 40 and pid in tags}
    assert eligible, "no simulated proteins > 40 aa among the recorded queries"
    ok = df[(df["class"] == "identical") & (df["gene"] == df["protein_id"])]
    hit = set(ok["protein_id"]) & eligible
    frac = len(hit) / len(eligible)
    assert frac >= 0.99, (
        f"{len(hit)}/{len(eligible)} = {frac:.4f} identical to own gene; "
        f"missing: {sorted(eligible - hit)[:10]}"
    )
