"""M3 gate, Route 1 share: >= 99.9 % recall of planted peptides (GIPPA_SPEC.md §13 row M3, §7.2.2).

gippasim genomes (circular, 50 kb and 100 kb); every protein is digested with trypsin (cut after
K/R, not before P, no missed cleavages); peptides of 7-40 aa are "planted" at the known gene
coordinates. N-terminal peptides of genes with a non-ATG start are planted too (M11-20 AC 10,
ADR-0068: their initiator is written as M over a GTG/TTG codon and ``match_exact`` credits it at
the annotated start cells, passed with ``starts=`` as ``map_exact`` does). Recall = share of
planted peptides whose origin interval is among the target rows; threshold ``mapping_recall_min``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from gippa.io.fasta import read_genome_fasta
from gippa.mapping import initiator as ini
from gippa.mapping.decoys import build_decoys
from gippa.mapping.exact import match_exact
from gippa.mapping.sixframe import build_sixframe
from gippa.model import Interval, Replicon
from gippa.model.genes import gene_model_from_bakta_json
from gippa.rng import for_step
from gippasim.codons import CODE_11
from gippasim.simulate import SimConfig, revcomp, simulate

pytestmark = pytest.mark.calibration


def _tryptic(aa: str) -> list[tuple[int, str]]:
    out, a = [], 0
    for m in re.finditer(r"[KR](?!P)", aa):
        out.append((a, aa[a : m.end()]))
        a = m.end()
    if a < len(aa):
        out.append((a, aa[a:]))
    return out


def translate(nt: str) -> str:
    return "".join(CODE_11[nt[i : i + 3]] for i in range(0, len(nt), 3))


@pytest.mark.parametrize(("length", "seed"), [(50_000, 101), (100_000, 102)])
def test_planted_peptide_recall(
    length: int, seed: int, tmp_path: Path, tolerances: dict[str, object]
) -> None:
    res = simulate(SimConfig(genome_length=length), seed=seed, outdir=tmp_path)
    (rec,) = read_genome_fasta(res.genome_fasta)
    rep = Replicon(rec.id, len(rec.sequence), "circular")
    model = gene_model_from_bakta_json(res.bakta_json, {rec.id: rep})
    idx = build_sixframe({rec.id: rec.sequence}, {rec.id: rep})
    ext = rec.sequence + rec.sequence[:120]
    planted: dict[str, set[Interval]] = {}
    for g in model.genes:
        iv = g.interval
        for off, pep in _tryptic(g.aa):
            if not 7 <= len(pep) <= 40:
                continue
            # minus strand: protein position 0 sits at the 3' end of the interval
            s = iv.start + 3 * off if iv.strand == "+" else iv.end - 3 * (off + len(pep))
            s %= rep.length
            want = Interval(rec.id, s, s + 3 * len(pep), iv.strand)
            nt = ext[want.start : want.end]
            aa = translate(nt if iv.strand == "+" else revcomp(nt))
            if off == 0 and g.start_codon != "ATG":  # the initiator M over GTG/TTG (M11-20)
                aa = "M" + aa[1:]
            assert aa == pep  # oracle check
            planted.setdefault(pep, set()).add(want)
    assert len(planted) > 100
    queries = sorted(planted)
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", seed))}
    genome = {rec.id: rec.sequence}
    sites, _ = ini.annotated_starts(model, genome, ("ATG", "GTG", "TTG"), table=11)
    tc = ini.start_frame_cells(idx, sites)
    starts = ini.StartCells(
        sites, tc, {"reverse": ini.mirror_start_cells("reverse", decoys["reverse"], tc)}
    )
    df = match_exact(queries, idx, decoys, dict.fromkeys(queries, "reverse"), starts=starts)
    tgt = df[df["search"] == "target"]
    found: dict[str, set[Interval]] = {}
    for r in tgt.itertuples(index=False):
        found.setdefault(r.peptide, set()).add(
            Interval(r.contig, int(r.start), int(r.end), r.strand)
        )
    n_planted = sum(len(v) for v in planted.values())
    n_found = sum(len(v & found.get(p, set())) for p, v in planted.items())
    recall = n_found / n_planted
    assert recall >= float(tolerances["mapping_recall_min"]), (recall, n_found, n_planted)  # type: ignore[arg-type]
