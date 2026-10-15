"""Contract test: tblastn (BLAST+ 2.16.0) parser and HSP merge on a recorded run (M3-13 ex. 8).

Marked ``pending_sample`` until ``tests/data/recorded/blast/2.16.0/sim100k/`` exists. The
directory holds the native tabular output (``*.tsv``/``*.tsv.gz`` in ``OUTFMT``), the query
``proteins.faa`` (leftover proteins: simulated CDS/sORF translations diverged to 80-95 %
identity, IDs = locus tags of their source gene) and ``PROVENANCE.md`` (both commands, image
digest or env lock, input sha256, date, host). The gene model is the giapsim 100-kb genome
(``SIM100K``, ``SIM100K_SEED``), regenerated here; a Bakta-like ``*.json`` in the directory is
used instead when present.

Example 8: the output parses, every merged interval is canonical, and the diverged proteins are
assigned ``divergent`` to their own gene (>= 95 % of the queries; a few may fall below the
80 % identity/coverage thresholds by chance).
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pytest

from giap.config.model import ProteinsConfig

RECORDED = (
    Path(__file__).resolve().parents[1] / "data" / "recorded" / "blast" / "2.16.0" / "sim100k"
)
OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
    "qlen nident positive"
)

pytestmark = [
    pytest.mark.pending_sample,
    pytest.mark.skipif(
        not RECORDED.is_dir(), reason=f"pending_sample: no recorded tblastn run in {RECORDED}"
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

    tmp = tmp_path_factory.mktemp("tblastn_contract")
    tsv = _one("*.tsv", "*.tsv.gz", "*.tab", "*.tab.gz")
    faa = _one("proteins.faa", "proteins.faa.gz", "*.faa", "*.faa.gz")
    assert tsv and faa, f"{RECORDED}: need the tblastn tabular output and proteins.faa"
    js = _one("*.json", "*.json.gz")
    if js is None:
        from giapsim import SIM100K, SIM100K_SEED, simulate

        js = simulate(SIM100K, SIM100K_SEED, tmp / "sim").bakta_json
    model = parse_bakta_json(_plain(js, tmp)).model
    proteins = {r.id: r.sequence for r in read_protein_fasta(_plain(faa, tmp))}
    return _plain(tsv, tmp), model, proteins


def test_provenance_recorded() -> None:
    prov = RECORDED / "PROVENANCE.md"
    assert prov.is_file(), "PROVENANCE.md (commands, digest or env lock, date, host) missing"
    text = prov.read_text()
    assert "tblastn" in text and "makeblastdb" in text and "2.16.0" in text
    # -parse_seqids is no longer required (M3-13 review decision 1; sim IDs are unaffected by it)
    for flag in ("-db_gencode 11", "-seg no", "-evalue 1e-5", "-dbtype nucl"):
        assert flag in text, f"command flag {flag!r} not recorded"
    assert OUTFMT in text, "outfmt string not recorded verbatim"
    assert re.search(r"sha256", text), "input checksums missing"
    assert re.search(r"\d{4}-\d{2}-\d{2}", text), "date missing"


def test_recorded_output_parses(recorded) -> None:  # type: ignore[no-untyped-def]
    from giap.mapping.tblastn import merge_hsps, parse_tblastn

    tsv, model, proteins = recorded
    hsps = parse_tblastn(tsv)
    n_lines = sum(1 for line in tsv.read_text().splitlines() if line.strip())
    assert hsps and len(hsps) == n_lines
    for h in hsps:
        assert h.qseqid in proteins, f"query {h.qseqid!r} not in proteins.faa"
        assert h.sseqid in model.replicons, f"subject {h.sseqid!r} not a sim replicon"
        assert h.qlen == len(proteins[h.qseqid].rstrip("*"))
        assert 1 <= h.qstart <= h.qend <= h.qlen
        assert 0 <= h.nident <= h.length
    alns = merge_hsps(hsps, model.replicons)
    assert alns
    for a in alns:
        rep = model.replicons[a.interval.contig]
        assert 0 <= a.interval.start < rep.length
        assert a.interval.end <= a.interval.start + rep.length
        if rep.topology == "linear":
            assert a.interval.end <= rep.length
        assert 0 < a.query_coverage <= 1.0
        assert 0.0 <= a.identity <= 1.0
        assert a.route == "tblastn"


def test_example8_diverged_proteins_divergent(recorded) -> None:  # type: ignore[no-untyped-def]
    from giap.mapping.protein_assign import assign_alignments
    from giap.mapping.tblastn import merge_hsps, parse_tblastn

    tsv, model, proteins = recorded
    tags = {g.locus_tag for g in model.genes if g.feature_type in ("cds", "sorf")}
    queries = set(proteins) & tags
    assert queries, "query IDs must be locus tags of the simulated CDS/sORF"
    df, _left = assign_alignments(
        merge_hsps(parse_tblastn(tsv), model.replicons), model, ProteinsConfig()
    )
    own = df[(df["gene"] == df["protein_id"]) & (df["class"] == "divergent")]
    hit = set(own["protein_id"]) & queries
    frac = len(hit) / len(queries)
    assert frac >= 0.95, (
        f"{len(hit)}/{len(queries)} = {frac:.3f} divergent on own gene; "
        f"missing: {sorted(queries - hit)[:10]}"
    )
    wrong = df[df["gene"].notna() & (df["gene"] != df["protein_id"])]
    assert wrong.empty, f"assigned to a foreign gene: {wrong[['protein_id', 'gene']].head()}"
