"""``parse_bakta_json`` on gippasim output (issue M2-03 acceptance; GIPPA_SPEC.md §11.3).

The M0 simulator JSON (``contig``, ``start``, ``stop``, ``locus``, ``start_type``, ``nt``, ``aa``,
no ``version``) parses with ``db_type="unknown"`` and reproduces ``truth/genes.tsv``. The M3-18
options (plasmids, a linear contig, rRNA/tRNA, pseudogenes, sORFs, an ``edge: true``
origin-crossing gene, ``reference_bakta_like.json``) are exercised when the simulator has them.
"""

from __future__ import annotations

import csv
import dataclasses
from pathlib import Path
from typing import Any

import pytest

from gippa.model import Interval
from gippasim import SimConfig, simulate

_FIELDS = {f.name for f in dataclasses.fields(SimConfig)}
needs_m318 = pytest.mark.skipif(
    "strain" not in _FIELDS, reason="gippasim lacks the M3-18 options (not merged yet)"
)


def _parse(path: Path, **kw: Any) -> Any:
    from gippa.annotate.bakta import parse_bakta_json

    return parse_bakta_json(path, **kw)


def _truth(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as fh:
        return {r["locus_tag"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif line:
            out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


_COMP = str.maketrans("ACGTN", "TGCAN")


def _slice(seq: str, iv: Interval) -> str:
    doubled = seq + seq
    s = doubled[iv.start : iv.end] if iv.end > len(seq) else seq[iv.start : iv.end]
    return s.translate(_COMP)[::-1] if iv.strand == "-" else s


def _assert_matches_truth(ann: Any, truth_path: Path) -> None:
    truth = _truth(truth_path)
    cds = {g.locus_tag: g for g in ann.model.genes if g.feature_type in ("cds", "pseudogene")}
    assert set(cds) == set(truth)
    for locus, row in truth.items():
        g = cds[locus]
        assert (g.interval.start + 1, g.interval.end, g.interval.strand) == (
            int(row["start"]),
            int(row["end"]),
            row["strand"],
        ), locus
        assert g.start_codon == row["start_codon"], locus


def test_m0_simulator_json(tmp_path: Path) -> None:
    res = simulate(SimConfig(genome_length=50_000), seed=1, outdir=tmp_path)
    ann = _parse(res.bakta_json)
    assert ann.db_type == "unknown"
    assert dict(ann.skipped_types) == {}
    rep = ann.model.replicons["chr1"]
    assert (rep.length, rep.topology, rep.kind) == (50_000, "circular", "chromosome")
    assert dict(ann.sequences) == _fasta(res.genome_fasta)
    assert len(ann.model.genes) == res.n_genes
    _assert_matches_truth(ann, res.genes_truth)


@needs_m318
def test_m318_simulator_json_all_features(tmp_path: Path) -> None:
    cfg = SimConfig(
        genome_length=100_000,
        n_plasmids=1,
        linear_contig_length=20_000,
        origin_crossing_gene=True,
        n_rrna_operons=1,
        n_trna=3,
        n_pseudogenes=2,
        n_sorfs=3,
    )
    res = simulate(cfg, seed=1, outdir=tmp_path)
    ann = _parse(res.bakta_json)
    assert ann.db_type == "unknown"
    reps = ann.model.replicons
    kinds = {c: (r.kind, r.topology) for c, r in reps.items()}
    assert kinds["chr1"] == ("chromosome", "circular")
    assert kinds["plasmid_1"] == ("plasmid", "circular")
    assert kinds["linear_1"] == ("unknown", "linear")

    types = [g.feature_type for g in ann.model.genes]
    assert types.count("pseudogene") == 2
    assert types.count("sorf") == 3
    crossing = [g for g in ann.model.genes if g.interval.crosses_origin(reps[g.interval.contig])]
    assert len(crossing) == 1
    assert crossing[0].interval.end > reps["chr1"].length

    nc = [f.type for f in ann.model.features]
    assert nc.count("rrna") == 3
    assert nc.count("trna") == 3
    assert all(
        f.anticodon and set(f.anticodon) <= set("ACGT")
        for f in ann.model.features
        if f.type == "trna"
    )
    assert dict(ann.skipped_types) == {}

    seqs = _fasta(res.genome_fasta)
    assert dict(ann.sequences) == seqs
    for g in ann.model.genes:
        assert g.nt == _slice(seqs[g.interval.contig], g.interval), g.locus_tag
    _assert_matches_truth(ann, res.genes_truth)


@needs_m318
def test_m318_reference_bakta_like_json(tmp_path: Path) -> None:
    from gippasim import StrainConfig

    cfg = SimConfig(
        genome_length=100_000,
        origin_crossing_gene=True,
        n_pseudogenes=1,
        n_sorfs=1,
        strain=StrainConfig(snv_rate_per_kb=1.0, indel_rate_per_kb=0.2),
    )
    res = simulate(cfg, seed=3, outdir=tmp_path)
    assert res.reference_json is not None and res.reference_fasta is not None
    ref = _parse(res.reference_json)
    sample = _parse(res.bakta_json)
    ref_seqs = _fasta(res.reference_fasta)
    assert dict(ref.sequences) == ref_seqs
    assert {g.locus_tag for g in ref.model.genes} == {g.locus_tag for g in sample.model.genes}
    for g in ref.model.genes:
        assert g.nt == _slice(ref_seqs[g.interval.contig], g.interval), g.locus_tag
    reps = ref.model.replicons
    assert sum(g.interval.crosses_origin(reps[g.interval.contig]) for g in ref.model.genes) == 1


@needs_m318
def test_m318_simulator_json_with_preflight_replicons(tmp_path: Path) -> None:
    from gippa.model import Replicon

    cfg = SimConfig(genome_length=60_000, n_plasmids=1, origin_crossing_gene=True)
    res = simulate(cfg, seed=2, outdir=tmp_path)
    seqs = _fasta(res.genome_fasta)
    reps = {
        "chr1": Replicon("chr1", len(seqs["chr1"]), "circular", "chromosome"),
        "plasmid_1": Replicon("plasmid_1", len(seqs["plasmid_1"]), "circular", "plasmid"),
    }
    ann = _parse(res.bakta_json, replicons=reps, genome=seqs)
    assert dict(ann.model.replicons) == reps
    _assert_matches_truth(ann, res.genes_truth)
