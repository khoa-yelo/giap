"""geNOMAD review tests (M2-18 lead decision 5 after the PR #84 code review).

Kills the surviving mutant that accepted position 0: geNOMAD coordinates are 1-based closed, so
``<contig>|provirus_0_<b>`` (and a virus gene starting at 0) is corrupt input → ``ValueError``,
while ``provirus_1_<b>`` is the first valid start (``Interval(contig, 0, b, ".")``). The error
is geNOMAD's own range check naming tool and contig, not ``Interval``'s generic one.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from giap.model.coords import Interval, Replicon

VS_HEADER = (
    "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\t"
    "n_hallmarks\tmarker_enrichment\ttaxonomy"
)
PS_HEADER = (
    "seq_name\tlength\ttopology\tn_genes\tgenetic_code\tplasmid_score\tfdr\tn_hallmarks\t"
    "marker_enrichment\tconjugation_genes\tamr_genes"
)
VG_HEADER = (
    "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\tmarker\tevalue\t"
    "bitscore\tuscg\tplasmid_hallmark\tvirus_hallmark\ttaxid\ttaxname\tannotation_conjscan\t"
    "annotation_amr\tannotation_accessions\tannotation_description"
)
TAX = "Viruses;Duplodnaviria;;"
REPS = {"chr1": Replicon("chr1", 100_000, "circular", "chromosome")}


def _write(tmp: Path, virus: list[str], genes: list[str]) -> Path:
    d = tmp / "genome_summary"
    d.mkdir(parents=True, exist_ok=True)
    for name, header, rows in (
        ("virus_summary", VS_HEADER, virus),
        ("plasmid_summary", PS_HEADER, []),
        ("virus_genes", VG_HEADER, genes),
    ):
        (d / f"genome_{name}.tsv").write_text("\n".join([header, *rows]) + "\n")
    return tmp


def _vs(a: int, b: int) -> str:
    return f"chr1|provirus_{a}_{b}\t{b - a + 1}\tProvirus\t{a}-{b}\t5\t11\t0.9\tNA\t2\t3.0\t{TAX}"


def _vg(gene: str, start: int, end: int) -> str:
    return (
        f"{gene}\t{start}\t{end}\t{end - start + 1}\t1\t0.5\t11\tNone\tNA\t"
        "NA\tNA\t0\t0\t0\t1\tNA\tNA\tNA\tNA\tNA"
    )


def test_provirus_at_position_zero_raises(tmp_path: Path) -> None:
    from giap.annotate.genomad import parse_genomad

    with pytest.raises(ValueError, match=r"geNOMAD.*chr1"):
        parse_genomad(_write(tmp_path, [_vs(0, 5000)], []), "genome", REPS)


def test_virus_gene_at_position_zero_raises(tmp_path: Path) -> None:
    from giap.annotate.genomad import parse_genomad

    out = _write(tmp_path, [_vs(1, 5000)], [_vg("chr1|provirus_1_5000_1", 0, 300)])
    with pytest.raises(ValueError, match=r"geNOMAD.*chr1"):
        parse_genomad(out, "genome", REPS)


def test_provirus_at_position_one_is_the_first_base(tmp_path: Path) -> None:
    from giap.annotate.genomad import parse_genomad

    out = _write(tmp_path, [_vs(1, 5000)], [_vg("chr1|provirus_1_5000_1", 1, 300)])
    res = parse_genomad(out, "genome", REPS)
    (p,) = res.proviruses
    assert p.interval == Interval("chr1", 0, 5000, ".")
    (g,) = res.virus_genes
    assert g.interval == Interval("chr1", 0, 300, "+")
