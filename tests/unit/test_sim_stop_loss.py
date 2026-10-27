"""gippasim strain stop-loss (M7-16b rule 20, D5; statistician review of the amendment, item 1;
acceptance "With ``n_stop_loss = 0`` every existing gippasim strain build is byte-identical").

``StrainConfig.n_stop_loss: int = 0``. With ``k > 0`` ``derive_reference`` also places ``k`` SNVs,
from a generator separate from the existing draws, each in the terminal stop codon of a distinct
non-pseudo sample CDS, turning it into a sense codon in the reference, ``_SPACING`` from every other
event, written as an ordinary ``snv`` row of ``strain_differences.tsv``. A gene is eligible only if
the reference protein then has a C-terminal digest product (``digest_len``, ``missed_cleavages``)
that contains the new residue, is >= 10 aa and holds no other strain event in its codons. Fewer
eligible genes than ``k``, or fewer wrong-peptide slots than ``k``, raise ``ValueError``.
``truth/stop_loss.tsv`` lists **every** stop-loss of the build (``planted`` true for the knob's,
false for the uniform draws', ``cterm_peptide`` empty when none qualifies). With ``wrong_source =
"other_strain"`` every planted ``cterm_peptide`` absent from the sample (I = L) is a wrong peptide.

The gippasim strain model is SNV-only in these builds, so sample and reference coordinates agree.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from gippasim.codons import CODE_11, STOPS
from gippasim.proteomics import ProteomicsConfig
from gippasim.simulate import SimConfig, simulate
from gippasim.strain import StrainConfig

COLUMNS = [
    "locus_tag",
    "contig",
    "strand",
    "codon_start",
    "sample_codon",
    "reference_codon",
    "reference_residue",
    "cterm_peptide",
    "planted",
]
SPACING = 12  # gippasim.strain._SPACING
PROT = ProteomicsConfig(search_db="other", pi_wrong=0.05, wrong_source="other_strain")
#: digests of the float-free files of two strain builds, recorded on main (c523c1c0 + the 16b
#: spec commits) before ``n_stop_loss`` existed
KEEP = (
    "genome.fasta",
    "reference.fasta",
    "bakta_like.json",
    "reference_bakta_like.json",
    "annotation.gff3",
    "proteins.faa",
    "truth/strain_differences.tsv",
    "truth/features.tsv",
    "truth/genes.tsv",
    "truth/peptides.tsv",
)
PINNED = {
    "other_strain_100k": "699ed56e520786433c62dd2dd8bb65444308d77111b2f82b00286a236a6be81d",
    "indels_50k": "dedb61f2805eda642212f691e3fd350ed1552b11b75f8df63ec73afc26b40a86",
}


def _config(case: str, **strain: Any) -> tuple[SimConfig, int]:
    if case == "other_strain_100k":
        return SimConfig(
            genome_length=100_000, strain=StrainConfig(snv_rate_per_kb=5, **strain), proteomics=PROT
        ), 11
    return SimConfig(
        genome_length=50_000,
        strain=StrainConfig(snv_rate_per_kb=2, indel_rate_per_kb=0.2, **strain),
        n_plasmids=1,
        origin_crossing_gene=True,
    ), 3


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


# ---------------------------------------------------------------- readers and oracle


def fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    name = ""
    for ln in path.read_text().splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            out[name] = []
        elif ln:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def codon(seq: str, a: int, strand: str) -> str:
    c = "".join(seq[(a + i) % len(seq)] for i in range(3))
    return c if strand == "+" else rc(c)


def tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def stop_loss(out: Path) -> pd.DataFrame:
    df = tsv(out / "truth" / "stop_loss.tsv")
    df["codon_start"] = df["codon_start"].astype(int)
    df["planted"] = df["planted"].str.lower().map({"true": True, "false": False})
    return df


def proteins(json_path: Path) -> dict[str, str]:
    doc = json.loads(json_path.read_text())
    return {f["locus"]: f.get("aa") or "" for f in doc["features"] if f.get("locus")}


def oracle_stop_losses(out: Path) -> dict[str, tuple[str, int, str, str, str]]:
    """Every non-pseudo sample CDS whose terminal stop is a sense codon in the reference (same
    coordinates: SNV-only builds): locus -> (contig, codon_start, strand, sample, reference)."""
    g, r = fasta(out / "genome.fasta"), fasta(out / "reference.fasta")
    f = tsv(out / "truth" / "features.tsv")
    found = {}
    for x in f.itertuples(index=False):
        if x.type != "cds" or x.pseudo_reason:
            continue
        a = int(x.end) - 3 if x.strand == "+" else int(x.start) - 1
        a %= len(g[x.contig])
        sc, rf = codon(g[x.contig], a, x.strand), codon(r[x.contig], a, x.strand)
        if sc in STOPS and rf not in STOPS:
            found[x.feature_id] = (x.contig, a, x.strand, sc, rf)
    return found


def _il(s: str) -> str:
    return s.replace("I", "L")


@pytest.fixture(scope="module")
def planted(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    """(k = 5, k = 0) builds of one seed: 100 kb, 5 SNV/kb, other-strain wrong peptides."""
    base = tmp_path_factory.mktemp("stoploss")
    cfg0, seed = _config("other_strain_100k")
    cfg5, _ = _config("other_strain_100k", n_stop_loss=5)
    return simulate(cfg5, seed, base / "k5").outdir, simulate(cfg0, seed, base / "k0").outdir


# ---------------------------------------------------------------- default off


@pytest.mark.parametrize("case", sorted(PINNED))
def test_default_build_is_byte_identical_to_main(tmp_path: Path, case: str) -> None:
    cfg, seed = _config(case)
    out = simulate(cfg, seed, tmp_path / "sim").outdir
    assert _digest(out) == PINNED[case]
    assert not (out / "truth" / "stop_loss.tsv").exists()


@pytest.mark.parametrize("case", sorted(PINNED))
def test_n_stop_loss_zero_is_the_default(tmp_path: Path, case: str) -> None:
    assert StrainConfig().n_stop_loss == 0  # type: ignore[attr-defined]
    cfg0, seed = _config(case, n_stop_loss=0)
    cfg, _ = _config(case)
    a = simulate(cfg0, seed, tmp_path / "zero").outdir
    b = simulate(cfg, seed, tmp_path / "default").outdir
    assert _files(a) == _files(b)
    assert _digest(a) == PINNED[case]
    assert not (a / "truth" / "stop_loss.tsv").exists()


# ---------------------------------------------------------------- planted stop-losses


def test_truth_table_columns_and_planted_count(planted: tuple[Path, Path]) -> None:
    k5, _ = planted
    df = tsv(k5 / "truth" / "stop_loss.tsv")
    assert list(df.columns) == COLUMNS
    sl = stop_loss(k5)
    assert int(sl["planted"].sum()) == 5
    assert sl.loc[sl["planted"], "locus_tag"].is_unique


def test_every_row_is_a_stop_loss_of_the_genomes(planted: tuple[Path, Path]) -> None:
    k5, _ = planted
    g, r = fasta(k5 / "genome.fasta"), fasta(k5 / "reference.fasta")
    feats = tsv(k5 / "truth" / "features.tsv").set_index("feature_id")
    for x in stop_loss(k5).itertuples(index=False):
        f = feats.loc[x.locus_tag]
        assert (f["type"], f["contig"], f["strand"], f["pseudo_reason"]) == (
            "cds",
            x.contig,
            x.strand,
            "",
        )
        n = len(g[x.contig])
        last = (int(f["end"]) - 3) % n if x.strand == "+" else int(f["start"]) - 1
        assert x.codon_start == last, "codon_start is the sample CDS's terminal codon (0-based)"
        assert codon(g[x.contig], x.codon_start, x.strand) == x.sample_codon
        assert codon(r[x.contig], x.codon_start, x.strand) == x.reference_codon
        assert x.sample_codon in STOPS and x.reference_codon not in STOPS
        assert sum(a != b for a, b in zip(x.sample_codon, x.reference_codon, strict=True)) == 1
        assert x.reference_residue == CODE_11[x.reference_codon]


def test_planted_cterm_peptide_is_the_reference_c_terminal_product(
    planted: tuple[Path, Path],
) -> None:
    k5, _ = planted
    ref, smp = proteins(k5 / "reference_bakta_like.json"), proteins(k5 / "bakta_like.json")
    lo, hi = PROT.digest_len
    for x in stop_loss(k5)[lambda d: d["planted"]].itertuples(index=False):
        pep, prot = x.cterm_peptide, ref[x.locus_tag]
        assert prot == smp[x.locus_tag] + x.reference_residue, x.locus_tag
        assert max(10, lo) <= len(pep) <= hi, pep
        assert prot.endswith(pep) and pep.endswith(x.reference_residue)
        before = prot[: len(prot) - len(pep)]
        assert before == "" or before[-1] in "KR", (x.locus_tag, pep)


def test_planted_snvs_are_ordinary_strain_differences(planted: tuple[Path, Path]) -> None:
    k5, _ = planted
    sd = tsv(k5 / "truth" / "strain_differences.tsv")
    sd["p0"] = sd["pos_sample"].astype(int) - 1
    for x in stop_loss(k5)[lambda d: d["planted"]].itertuples(index=False):
        rows = sd[
            (sd["contig"] == x.contig)
            & (sd["p0"] >= x.codon_start)
            & (sd["p0"] < x.codon_start + 3)
        ]
        assert len(rows) == 1, (x.locus_tag, rows)
        (row,) = rows.to_dict("records")
        assert (row["type"], row["gene"]) == ("snv", x.locus_tag)


def test_planted_snvs_use_a_separate_generator(planted: tuple[Path, Path]) -> None:
    """The k = 0 draws are unchanged: same sample genome, every k = 0 strain row kept, exactly
    the k planted rows added, the references differing at exactly those k bases."""
    k5, k0 = planted
    assert (k5 / "genome.fasta").read_bytes() == (k0 / "genome.fasta").read_bytes()
    a, b = (
        tsv(k0 / "truth" / "strain_differences.tsv"),
        tsv(k5 / "truth" / "strain_differences.tsv"),
    )
    key = ["contig", "pos_sample", "pos_reference", "ref", "alt", "type", "gene"]
    rows0 = set(map(tuple, a[key].to_numpy().tolist()))
    rows5 = set(map(tuple, b[key].to_numpy().tolist()))
    assert rows0 <= rows5 and len(rows5 - rows0) == 5
    r0, r5 = fasta(k0 / "reference.fasta"), fasta(k5 / "reference.fasta")
    diffs = sum(x != y for c in r0 for x, y in zip(r0[c], r5[c], strict=True))
    assert diffs == 5


def test_planted_snvs_keep_spacing_and_their_peptides_are_clean(
    planted: tuple[Path, Path],
) -> None:
    k5, _ = planted
    sd = tsv(k5 / "truth" / "strain_differences.tsv")
    sd["p0"] = sd["pos_sample"].astype(int) - 1
    n = {c: len(s) for c, s in fasta(k5 / "genome.fasta").items()}
    for x in stop_loss(k5)[lambda d: d["planted"]].itertuples(index=False):
        on = sd[sd["contig"] == x.contig]
        mine = on[(on["p0"] >= x.codon_start) & (on["p0"] < x.codon_start + 3)]
        others = on.drop(mine.index)
        p = int(mine["p0"].iloc[0])
        dist = [min(abs(p - q), n[x.contig] - abs(p - q)) for q in others["p0"]]
        assert min(dist, default=SPACING) >= SPACING, (x.locus_tag, min(dist))
        span = 3 * len(x.cterm_peptide)
        lo, hi = (
            (x.codon_start + 3 - span, x.codon_start + 3)
            if x.strand == "+"
            else (x.codon_start, x.codon_start + span)
        )
        inside = [q for q in others["p0"] if lo <= q < hi or lo <= q + n[x.contig] < hi]
        assert not inside, (x.locus_tag, inside)


def test_planted_cterm_peptides_join_the_wrong_peptides(planted: tuple[Path, Path]) -> None:
    k5, _ = planted
    truth = tsv(k5 / "truth" / "peptides.tsv")
    wrong = set(
        truth[(truth["is_wrong"] == "true") & (truth["origin"] == "other_strain")]["peptide"]
    )
    sample = "|".join(_il(a) for a in proteins(k5 / "bakta_like.json").values())
    n_absent = 0
    for x in stop_loss(k5)[lambda d: d["planted"]].itertuples(index=False):
        if _il(x.cterm_peptide) in sample:
            continue
        n_absent += 1
        assert x.cterm_peptide in wrong, x.cterm_peptide
    assert n_absent >= 1


def test_planted_build_is_deterministic(tmp_path: Path, planted: tuple[Path, Path]) -> None:
    k5, _ = planted
    cfg, seed = _config("other_strain_100k", n_stop_loss=5)
    again = simulate(cfg, seed, tmp_path / "again").outdir
    assert _files(again) == _files(k5)


# ---------------------------------------------------------------- every stop-loss is listed


def test_natural_stop_losses_are_listed_unplanted(tmp_path: Path) -> None:
    """20 SNV/kb at seed 6 gives 6 uniform-draw stop-losses (measured at k = 0 on main)."""
    cfg = SimConfig(genome_length=100_000, strain=StrainConfig(snv_rate_per_kb=20, n_stop_loss=3))
    out = simulate(cfg, 6, tmp_path / "sim").outdir
    want = oracle_stop_losses(out)
    sl = stop_loss(out)
    assert set(sl["locus_tag"]) == set(want)
    assert int(sl["planted"].sum()) == 3
    natural = sl[~sl["planted"]]
    assert len(natural) >= 6
    ref = proteins(out / "reference_bakta_like.json")
    for x in sl.itertuples(index=False):
        assert want[x.locus_tag] == (
            x.contig,
            x.codon_start,
            x.strand,
            x.sample_codon,
            x.reference_codon,
        )
    for x in natural.itertuples(index=False):
        pep = x.cterm_peptide
        assert pep == "" or (ref[x.locus_tag].endswith(pep) and pep.endswith(x.reference_residue))


# ---------------------------------------------------------------- no silent shortfall


def test_too_few_eligible_genes_raise(tmp_path: Path) -> None:
    cfg = SimConfig(
        genome_length=50_000, strain=StrainConfig(snv_rate_per_kb=5, n_stop_loss=10_000)
    )
    with pytest.raises(ValueError):
        simulate(cfg, 4, tmp_path / "sim")


def test_too_few_wrong_peptide_slots_raise(tmp_path: Path) -> None:
    prot = ProteomicsConfig(search_db="other", pi_wrong=0.0, wrong_source="other_strain")
    cfg = SimConfig(
        genome_length=100_000,
        strain=StrainConfig(snv_rate_per_kb=5, n_stop_loss=2),
        proteomics=prot,
    )
    with pytest.raises(ValueError):
        simulate(cfg, 11, tmp_path / "sim")
