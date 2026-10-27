"""gippasim stop-codon read-through (M7-16 rule 14; PR 16b-3, lead decision in the #335 review
round; ADR-0058; acceptance "16b": sim100k with ``n_readthrough = 10``).

``SimConfig.n_readthrough: int = 0``. With 0 every existing scenario and seed is byte-identical
and no ``truth/readthrough.tsv`` is written. With ``n > 0`` the simulator picks ``n`` annotated
genes whose next in-frame stop is 5-60 codons after the annotated stop, inside the replicon
(wrapping on circular ones). The expressed protein reads the stop as ``W`` (TGA) or ``Q``
(TAA/TAG) and continues to the next stop; the annotation keeps the original stop. Truth:
``truth/readthrough.tsv`` with ``locus_tag, contig, strand, codon_start, codon, residue,
ext_stop, n_true_peptides_spanning``.

The oracle (``helpers/readthrough_sim_oracle.py``) pins the interpretations of the spec's open
points: eligible = non-pseudo CDS and sORFs; ``k`` (codons from the annotated stop to the new
one) in ``[5, 60]``; ``ext_stop`` the new stop's lowest coordinate, wrapped; the selection draws
after the genome is built, so genome and annotation are byte-identical to ``n = 0``; fewer
eligible genes than ``n`` raise ``ValueError`` (no silent shortfall, the ``n_stop_loss``
convention of rule 20); ``n_true_peptides_spanning`` counts detected correct peptides.

``MIXED`` (seed 39) holds every geometry: a circular chromosome whose origin-crossing gene is
eligible with its stop codon past the origin (wrapped ``codon_start``), a circular plasmid, a
linear contig with one gene whose scan leaves the contig before the next stop (not eligible),
eligible sORFs, genes whose next stop is < 5 codons away and genes with none within 60.
"""

from __future__ import annotations

import dataclasses
import hashlib
from pathlib import Path
from typing import Any

import pytest
from helpers.readthrough_sim_oracle import (
    K_MAX,
    K_MIN,
    READTHROUGH_COLUMNS,
    RESIDUE,
    codon_at,
    eligible,
    expressed,
    extension,
    holds,
    read_genome,
    readthrough_table,
    spanning_peptides,
    step,
    translate_interval,
    tsv,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import gippasim
from gippasim.codons import STOPS
from gippasim.denovo import DenovoConfig
from gippasim.proteomics import ProteomicsConfig
from gippasim.simulate import SimConfig, simulate

MIXED = SimConfig(
    genome_length=30_000,
    n_plasmids=1,
    plasmid_length=(5000, 5000),
    linear_contig_length=4000,
    origin_crossing_gene=True,
    n_sorfs=3,
    n_pseudogenes=2,
    proteomics=ProteomicsConfig(),
)
MIXED_SEED = 39
#: measured on main 1549a930 with the oracle: 26 eligible genes (3 sORFs, the origin-crossing
#: GSIM_00155), GSIM_00180 on the linear contig runs off its end
MIXED_ELIGIBLE = 26
SIM100K_DENOVO = dataclasses.replace(
    gippasim.SIM100K, proteomics=ProteomicsConfig(), denovo=DenovoConfig(), write_precomputed=True
)
#: float-free files of a build, pinned to main below. ``peptides.tsv``, ``abundance.tsv`` and the
#: de novo table carry formatted floats and are not pinned: the n = 0 vs default comparison only
#: checks the explicit 0 against the (also 0) default. The #350 review checked byte-identity of
#: whole trees against main by hand (three configs, float files included).
KEEP = (
    "genome.fasta",
    "bakta_like.json",
    "annotation.gff3",
    "proteins.faa",
    "truth/features.tsv",
    "truth/genes.tsv",
    "truth/peptides.tsv",
    "truth/denovo.tsv",
)
#: the annotation side, which read-through never changes
GENOME_SIDE = (
    "genome.fasta",
    "bakta_like.json",
    "annotation.gff3",
    "proteins.faa",
    "truth/features.tsv",
    "truth/genes.tsv",
    "truth/novel_orfs.tsv",
)
#: digests of KEEP recorded on main 1549a930, before ``n_readthrough`` existed
PINNED = {
    "mixed": "28307def866b2d9bdeafb296fcd6d6b45bf63cc4001734d8499ec80305642ce5",
    "sim100k_denovo": "779151d6c41e27f59826cba423a07110f7746129c68e2bf36a8a1a642953e7a8",
}
CASES: dict[str, tuple[SimConfig, int]] = {
    "mixed": (MIXED, MIXED_SEED),
    "sim100k_denovo": (SIM100K_DENOVO, gippasim.SIM100K_SEED),
}


def _rt(cfg: SimConfig, n: int) -> SimConfig:
    return dataclasses.replace(cfg, n_readthrough=n)  # type: ignore[call-arg]


def _digest(out: Path) -> str:
    h = hashlib.sha256()
    for rel in sorted(KEEP):
        p = out / rel
        if p.exists():
            h.update(rel.encode() + b"\0" + hashlib.sha256(p.read_bytes()).digest())
    return h.hexdigest()


def _files(out: Path) -> dict[str, bytes]:
    return {
        str(p.relative_to(out)): p.read_bytes()
        for p in sorted(out.rglob("*"))
        if p.is_file() and p.name != "config.yaml"
    }


@pytest.fixture(scope="module")
def mixed(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """MIXED at n = 0, 10 and every eligible gene (26)."""
    base = tmp_path_factory.mktemp("rt_mixed")
    return {
        str(n): simulate(_rt(MIXED, n), MIXED_SEED, base / f"n{n}").outdir
        for n in (0, 10, MIXED_ELIGIBLE)
    }


@pytest.fixture(scope="module")
def sim100k(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """The e2e build (sim100k, de novo) at n = 0 and n = 10."""
    base = tmp_path_factory.mktemp("rt_100k")
    return {
        str(n): simulate(_rt(SIM100K_DENOVO, n), gippasim.SIM100K_SEED, base / f"n{n}").outdir
        for n in (0, 10)
    }


# ---------------------------------------------------------------- default off


def test_the_knob_defaults_to_zero() -> None:
    assert SimConfig().n_readthrough == 0  # type: ignore[attr-defined]
    assert gippasim.SIM100K.n_readthrough == 0  # type: ignore[attr-defined]


@pytest.mark.parametrize("case", sorted(CASES))
def test_default_build_is_byte_identical_to_main(tmp_path: Path, case: str) -> None:
    cfg, seed = CASES[case]
    out = simulate(cfg, seed, tmp_path / "sim").outdir
    assert _digest(out) == PINNED[case]
    assert not (out / "truth" / "readthrough.tsv").exists()


@pytest.mark.parametrize("case", sorted(CASES))
def test_n_readthrough_zero_equals_the_default(tmp_path: Path, case: str) -> None:
    cfg, seed = CASES[case]
    a = simulate(_rt(cfg, 0), seed, tmp_path / "zero").outdir
    b = simulate(cfg, seed, tmp_path / "default").outdir
    assert _files(a) == _files(b)
    assert _digest(a) == PINNED[case]
    assert not (a / "truth" / "readthrough.tsv").exists()


def test_negative_n_readthrough_raises() -> None:
    with pytest.raises(ValueError):
        SimConfig(n_readthrough=-1)  # type: ignore[call-arg]


# ---------------------------------------------------------------- the truth table


@pytest.mark.parametrize("n", [10, MIXED_ELIGIBLE])
def test_truth_table_columns_and_count(mixed: dict[str, Path], n: int) -> None:
    out = mixed[str(n)]
    df = tsv(out / "truth" / "readthrough.tsv")
    assert list(df.columns) == READTHROUGH_COLUMNS
    assert len(df) == n
    assert df["locus_tag"].is_unique


@pytest.mark.parametrize("n", [10, MIXED_ELIGIBLE])
def test_every_row_is_an_eligible_gene_with_its_own_stop(mixed: dict[str, Path], n: int) -> None:
    """codon_start is the gene's terminal codon, codon its gene-strand bases, residue W/Q, and
    ext_stop the first in-frame stop after it, 5-60 codons on, inside the replicon."""
    out = mixed[str(n)]
    g = read_genome(out)
    feats = tsv(out / "truth" / "features.tsv").set_index("feature_id")
    elig = eligible(mixed["0"])
    for x in readthrough_table(out).itertuples(index=False):
        f = feats.loc[x.locus_tag]
        assert f["type"] in ("cds", "sorf") and f["pseudo_reason"] == "", x.locus_tag
        assert (f["contig"], f["strand"]) == (x.contig, x.strand)
        nn = g.n(x.contig)
        last = (int(f["end"]) - 3) % nn if x.strand == "+" else int(f["start"]) - 1
        assert x.codon_start == last, "codon_start is the gene's terminal codon (0-based, wrapped)"
        assert codon_at(g, x.contig, x.codon_start, x.strand) == x.codon
        assert x.codon in STOPS and x.residue == RESIDUE[x.codon]
        assert 0 <= x.codon_start < nn and 0 <= x.ext_stop < nn
        k = ((x.ext_stop - x.codon_start) if x.strand == "+" else (x.codon_start - x.ext_stop)) % nn
        assert k % 3 == 0 and K_MIN <= k // 3 <= K_MAX, (x.locus_tag, k)
        assert codon_at(g, x.contig, x.ext_stop, x.strand) in STOPS
        for j in range(1, k // 3):
            c = codon_at(g, x.contig, x.codon_start + step(x.strand, j), x.strand)
            assert c is not None and c not in STOPS, (x.locus_tag, j)
        e = elig[x.locus_tag]
        assert (e.codon_start, e.codon, e.ext_stop) == (x.codon_start, x.codon, x.ext_stop)


def test_selection_is_a_subset_of_the_eligible_genes(mixed: dict[str, Path]) -> None:
    elig = set(eligible(mixed["0"]))
    assert len(elig) == MIXED_ELIGIBLE
    assert set(readthrough_table(mixed["10"])["locus_tag"]) <= elig


def test_every_eligible_gene_is_selected_when_n_equals_the_supply(mixed: dict[str, Path]) -> None:
    """Exhaustive: the wrapped stop of the origin-crossing gene and the sORFs are in, the linear
    contig's edge gene, the genes with a stop < 5 codons on and those with none within 60 out."""
    elig = eligible(mixed["0"])
    got = readthrough_table(mixed[str(MIXED_ELIGIBLE)])
    assert set(got["locus_tag"]) == set(elig)
    feats = tsv(mixed["0"] / "truth" / "features.tsv").set_index("feature_id")
    lengths = {c: len(s) for c, s in read_genome(mixed["0"]).seqs.items()}
    crossing = [t for t in elig if int(feats.loc[t, "end"]) > lengths[feats.loc[t, "contig"]]]
    assert crossing, "MIXED's origin-crossing gene is eligible (its stop lies past the origin)"
    assert {e.kind for e in elig.values()} == {"cds", "sorf"}
    linear = feats[feats["contig"].str.startswith("linear")]
    assert set(linear.index) - set(elig), "a linear-contig gene runs off the edge"


def test_too_few_eligible_genes_raise(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        simulate(_rt(MIXED, MIXED_ELIGIBLE + 1), MIXED_SEED, tmp_path / "sim")


# ---------------------------------------------------------------- the annotation keeps the stop


@pytest.mark.parametrize("n", ["10", str(MIXED_ELIGIBLE)])
def test_genome_and_annotation_are_unchanged(mixed: dict[str, Path], n: str) -> None:
    for rel in GENOME_SIDE:
        a, b = mixed["0"] / rel, mixed[n] / rel
        assert a.read_bytes() == b.read_bytes(), rel
    pre = mixed["0"] / "precomputed"
    for p in sorted(pre.rglob("*")) if pre.exists() else []:
        if p.is_file():
            assert (mixed[n] / p.relative_to(mixed["0"])).read_bytes() == p.read_bytes(), p.name


def test_sim100k_genome_and_annotation_are_unchanged(sim100k: dict[str, Path]) -> None:
    for rel in GENOME_SIDE:
        assert (sim100k["0"] / rel).read_bytes() == (sim100k["10"] / rel).read_bytes(), rel
    assert len(readthrough_table(sim100k["10"])) == 10


# ---------------------------------------------------------------- the expressed protein


def _check_truth_peptides(out: Path) -> None:
    """Every correct truth row of a planted gene translates from the genome with the stop read
    as the planted residue (so no row reaches the new stop); no other gene has a row over its
    own terminal stop."""
    g = read_genome(out)
    rt = readthrough_table(out).set_index("locus_tag")
    feats = tsv(out / "truth" / "features.tsv").set_index("feature_id")
    t = tsv(out / "truth" / "peptides.tsv")
    t = t[(t["is_wrong"] == "false") & (t["locus_tag"].isin(feats.index))]
    for r in t.itertuples(index=False):
        f = feats.loc[r.locus_tag]
        if f["type"] not in ("cds", "sorf"):
            continue
        c, n = f["contig"], g.n(f["contig"])
        s0, e0 = int(r.start) - 1, int(r.end)
        cs = (int(f["end"]) - 3) % n if f["strand"] == "+" else int(f["start"]) - 1
        over = holds(s0, e0, cs, n, g.circular[c])
        if r.locus_tag not in rt.index:
            assert not over, ("a non-planted gene's peptide over its stop", r.peptide)
            continue
        res = rt.loc[r.locus_tag, "residue"]
        aa = translate_interval(g, c, s0, e0, f["strand"], cs, res)
        want = r.peptide[1:] if r.alt_start_nterm == "true" else r.peptide
        got = aa[1:] if r.alt_start_nterm == "true" else aa
        assert got == want, (r.locus_tag, r.peptide, aa)
        assert "*" not in aa


def test_truth_peptides_follow_the_expressed_protein(mixed: dict[str, Path]) -> None:
    _check_truth_peptides(mixed[str(MIXED_ELIGIBLE)])


def test_sim100k_truth_peptides_follow_the_expressed_protein(sim100k: dict[str, Path]) -> None:
    _check_truth_peptides(sim100k["10"])


@pytest.mark.parametrize("n", ["10", str(MIXED_ELIGIBLE)])
def test_n_true_peptides_spanning_counts_detected_peptides_over_the_stop(
    mixed: dict[str, Path], n: str
) -> None:
    out = mixed[n]
    for x in readthrough_table(out).itertuples(index=False):
        peps = spanning_peptides(out, x.locus_tag, x.codon_start, x.contig)
        assert x.n_true_peptides_spanning == len(peps), (x.locus_tag, sorted(peps))


def test_spanning_peptides_are_products_of_the_extended_protein(mixed: dict[str, Path]) -> None:
    out = mixed[str(MIXED_ELIGIBLE)]
    g = read_genome(out)
    elig = eligible(mixed["0"])
    total = 0
    for x in readthrough_table(out).itertuples(index=False):
        e = elig[x.locus_tag]
        prot = expressed(out, e, g)
        assert len(extension(g, e)) == e.k - 1
        for pep in spanning_peptides(out, x.locus_tag, x.codon_start, x.contig):
            assert pep in prot and x.residue in pep, (x.locus_tag, pep)
            total += 1
    assert total >= 1, "MIXED has detected read-through peptides"


def test_sim100k_has_spanning_peptides(sim100k: dict[str, Path]) -> None:
    rt = readthrough_table(sim100k["10"])
    for x in rt.itertuples(index=False):
        peps = spanning_peptides(sim100k["10"], x.locus_tag, x.codon_start, x.contig)
        assert x.n_true_peptides_spanning == len(peps), x.locus_tag
    assert int((rt["n_true_peptides_spanning"] >= 1).sum()) >= 3
    assert int((rt["n_true_peptides_spanning"] >= 2).sum()) >= 1


def test_planted_build_is_deterministic(tmp_path: Path, mixed: dict[str, Path]) -> None:
    again = simulate(_rt(MIXED, 10), MIXED_SEED, tmp_path / "again").outdir
    assert _files(again) == _files(mixed["10"])


# ---------------------------------------------------------------- property over seeds


@settings(
    max_examples=12,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture, HealthCheck.too_slow],
)
@given(seed=st.integers(0, 10_000), want=st.integers(1, 4))
def test_property_selection_and_truth_agree_with_the_oracle(
    tmp_path_factory: pytest.TempPathFactory, seed: int, want: int
) -> None:
    cfg = SimConfig(
        genome_length=12_000,
        n_plasmids=1,
        plasmid_length=(3000, 3000),
        origin_crossing_gene=True,
        n_sorfs=1,
        proteomics=ProteomicsConfig(),
    )
    base = tmp_path_factory.mktemp("rt_prop")
    zero = simulate(cfg, seed, base / "n0").outdir
    elig = eligible(zero)
    n = min(want, len(elig))
    if n == 0:
        with pytest.raises(ValueError):
            simulate(_rt(cfg, want), seed, base / "nx")
        return
    out = simulate(_rt(cfg, n), seed, base / "n").outdir
    for rel in GENOME_SIDE:
        assert (zero / rel).read_bytes() == (out / rel).read_bytes(), rel
    rt = readthrough_table(out)
    assert len(rt) == n and rt["locus_tag"].is_unique
    for x in rt.itertuples(index=False):
        e = elig[x.locus_tag]
        got: tuple[Any, ...] = (x.contig, x.strand, x.codon_start, x.codon, x.residue, x.ext_stop)
        assert got == (e.contig, e.strand, e.codon_start, e.codon, RESIDUE[e.codon], e.ext_stop)
        peps = spanning_peptides(out, x.locus_tag, x.codon_start, x.contig)
        assert x.n_true_peptides_spanning == len(peps)
    _check_truth_peptides(out)
