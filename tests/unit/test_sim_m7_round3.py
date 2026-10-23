"""giapsim M7-12 round 3: lead decisions 21 and 22 after the statistician's round-2 block on
PR #225 (comment 6358601933; docs/specs/M7-12-sim-variant-novel-scenarios.md; GIAP_SPEC.md §7.6.6,
§7.9.2, §8.4).

The oracle is tests/helpers/m7_round3_oracle.py (no giapsim import). Decoy IDs are checked against
the database GIAP's own ``giap makedb --six-frame-orfs`` builds from the simulator's
``bakta_like.json`` and genome (``min_orf_aa = 7``), because that is what M7-11's picked class FDR
pairs against. The 1-Mb case is ``slow``; the pooled distribution checks are in
tests/calibration/test_sim_m7_decoy_calibration.py.

Interface pinned here (``wrong_source = "sixframe_noncoding"`` only; every other config keeps its
bytes):

* **Decision 21.**
  - A six-frame decoy row of ``peptides.tsv`` has ``is_decoy = true`` and ``protein_ids`` =
    ``rev_GIAPORF_{class}_{contig}_{start1}_{end1}_{p|m}``. With ``rev_`` stripped, that is the
    ``db_id`` of a makedb entry: the stop-to-stop six-frame ORF the product lies in, or for
    ``nterm_extension`` the gene's makedb N-terminal entry (upstream part + first 30 codons). The
    entry's protein contains the reversed decoy peptide, and the ID's class is the entry's makedb
    class.
  - Each six-frame wrong target (its ``truth/peptides.tsv`` locus) lies in exactly one such
    entry, so its decoy ID is defined.
  - Decoy ``score`` values are finite draws from the null N(mu0, sd0) of
    ``ProteomicsConfig.score_mixture``, the null the wrong targets use. ``pep`` and ``q_value``
    come from the same mixture with the target rows' ``pi`` (``truth/proteomics.json``
    ``pi_wrong_realised``).
  - Each null PSM is a wrong target row or a six-frame decoy row with probability 1/2 (pooled
    check in the calibration file).
  - The ``rev_<locus_tag>`` protein decoys stay (reversed products of annotated proteins), so
    GIAP's gene-level ``pair_decoys`` / ``choose_tier`` still give ``decoy_picked`` with no
    ``decoy_pairing_low``.
* **Decision 22.** The six-frame wrong peptides come from four territories, drawn in turn:
  intergenic, cds_outframe, antisense (round 2) and ``nterm_extension``, the products in frame
  with an annotated CDS/sORF upstream of its true start with no in-frame stop in between. The test
  classes peptides with the round-2 ``wrong_class`` rule, after an ``nterm_extension`` check placed
  first (oracle ``wrong_class4``). GIAP's ``LocusClassifier`` must class the nterm-territory
  products ``nterm_extension``.
* Configs that do not use the six-frame source (``variant_f1``, ``isobaric_attack``,
  ``novel_orf_power`` with ``wrong_source = "shuffled"``) are byte-identical to 9bd0b7a6. So are
  the non-proteomics files of ``novel_orf_power``.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import m7_round2_oracle as R
from helpers import m7_round3_oracle as Q
from helpers import m7_sim_oracle as O
from scipy import stats

import giapsim

SEED = 20261022
KS_MIN_P = 0.001  # tests/calibration/tolerances.yaml ks_min_p

#: 9bd0b7a6 outputs (seed SEED, 100 kb) of the configs decisions 21-22 must not touch.
VARIANT_F1_SHA256 = {
    "abundance.tsv": "c7576c8c6902c1a74892ddf5edd453fa73d67315b993243485d51a9960e7fba6",
    "annotation.gff3": "239698807121fd27c9da3a01ce03d6ceb353ef2f3657f3e357ca43bbdafaaf73",
    "bakta_like.json": "7cf8b02cb7ada6536cca689db577882b28d47a7ec5ac4c1daa1b27c6831e1d9b",
    "genome.fasta": "fc8f62a258ffc2b22cd0e81cb833083101942c0ba91d89dfc1d1a822afcbfe08",
    "peptides.tsv": "45bfcd9b2bbaeee0a0b74e896de3200490d28eb7b90aa4d20aaa5ed8445dbf17",
    "peptides_denovo.tsv": "0e103bbf4f7052e40c6dd287436010b2d8da4605d78494032ef24b8b9cc397a0",
    "proteins.faa": "adb46b9a094ac2376d021de57e4c7a91be21d460f5b53fdf2f0b453f449d6fb6",
    "reference.fasta": "617be93bdb838d30616c867006dd978e5ce42570460247a2a9fb5e8f310ff2ea",
    "reference_bakta_like.json": "be7cb0e9cc68b0ba43145bfe05d636f8e94a9aed390cc1d73a726fbad19b06f3",
    "reference_peptides.tsv": "52ffc6e46341519a6b3e14b46466a53ea9cf6ee157eb9261db42766d72daa0c6",
    "truth/alt_sites.tsv": "bb2a0b018202b96e50f8aa31898311fb7ae074b09d98750e65aa1095723898ec",
    "truth/denovo.tsv": "449fa9e3de34c6785454bbb9ec52a827e1a94dd3362e5b8234ebffb0b4645b6e",
    "truth/denovo_precision.tsv": "688c31e3929158135c30bcccbfd7599aeb271f8aafe13ac25f3e4970eb91129a",  # noqa: E501
    "truth/features.tsv": "11f87f300857c69fbaae6a0a5d524a2ddcc3846d759019e02ec8c91cb3449b86",
    "truth/genes.tsv": "4c581641292d96e35d0f1fea5d4658199efc4c987a1f3d7c71b1b9ecfb7e2b3f",
    "truth/indel_sites.tsv": "f885380ebd17b8689b5cbcf04769790e297613ce128b807544c047d6a6b533a3",
    "truth/novel_orfs.tsv": "ca3502b76760e37c9a78f8ba2b0bae532d6eb9b68ef1b57c80460bbc160a83a4",
    "truth/peptides.tsv": "24f899375b38a69aaab97dd98cff853825d67cd203ee7c20c0dd528901f08897",
    "truth/proteomics.json": "876c25f36509a08793a6218c1050525e95b14e9659744f6f31152a5c8354007c",
    "truth/strain_differences.tsv": "d22c76d35d4b4990bf68d2a75237b8b4a9a6e47e9ac8570954abba5347d95978",  # noqa: E501
}  # fmt: skip
ISOBARIC_SHA256 = {
    "abundance.tsv": "c7576c8c6902c1a74892ddf5edd453fa73d67315b993243485d51a9960e7fba6",
    "annotation.gff3": "239698807121fd27c9da3a01ce03d6ceb353ef2f3657f3e357ca43bbdafaaf73",
    "bakta_like.json": "7cf8b02cb7ada6536cca689db577882b28d47a7ec5ac4c1daa1b27c6831e1d9b",
    "genome.fasta": "fc8f62a258ffc2b22cd0e81cb833083101942c0ba91d89dfc1d1a822afcbfe08",
    "peptides.tsv": "45bfcd9b2bbaeee0a0b74e896de3200490d28eb7b90aa4d20aaa5ed8445dbf17",
    "peptides_denovo.tsv": "23e6b4be2412e75eb8b6aa6b941eb92e4e43c2d1b3e070f12678d6879747180f",
    "proteins.faa": "adb46b9a094ac2376d021de57e4c7a91be21d460f5b53fdf2f0b453f449d6fb6",
    "reference_peptides.tsv": "52ffc6e46341519a6b3e14b46466a53ea9cf6ee157eb9261db42766d72daa0c6",
    "truth/alt_sites.tsv": "074935cea3c4e0fca14c95b0375855716dc8ba17ecde7ff770605519a398b8c6",
    "truth/denovo.tsv": "1226ab8a1ef2e30d33b3990fc9a6404dec362fa925676c1fc9125ad06d1bfe02",
    "truth/denovo_precision.tsv": "688c31e3929158135c30bcccbfd7599aeb271f8aafe13ac25f3e4970eb91129a",  # noqa: E501
    "truth/features.tsv": "e0b807008a6f0f71a51a07de62e3a386152d3b48f67582071a6dc52e8855f87b",
    "truth/genes.tsv": "4c581641292d96e35d0f1fea5d4658199efc4c987a1f3d7c71b1b9ecfb7e2b3f",
    "truth/novel_orfs.tsv": "ca3502b76760e37c9a78f8ba2b0bae532d6eb9b68ef1b57c80460bbc160a83a4",
    "truth/peptides.tsv": "24f899375b38a69aaab97dd98cff853825d67cd203ee7c20c0dd528901f08897",
    "truth/proteomics.json": "876c25f36509a08793a6218c1050525e95b14e9659744f6f31152a5c8354007c",
}  # fmt: skip
#: ``novel_orf_power`` with ``wrong_source = "shuffled"``: the proteomics files (9bd0b7a6).
NOVEL_SHUFFLED_SHA256 = {
    "peptides.tsv": "e339c1b02b3fefa6fc253480aeaae4c6925879bba53ab85d5e561feaff723df9",
    "truth/peptides.tsv": "6eb1e00271ded401771d41ceaa9e653cc7d4617086643e305292893ade75d06e",
    "truth/proteomics.json": "c0350df4568aadb68fced37c970c071481e75b7d9c7b5093e91924a2f5118ad2",
    "peptides_denovo.tsv": "d7a1d869a631305183bf2b73ebc56dc4d48e52dd54561d3fa8056f97a1d15269",
}  # fmt: skip
#: ``novel_orf_power`` (six-frame): files outside the proteomics step (9bd0b7a6).
NOVEL_NONPROTEOMICS_SHA256 = {
    "annotation.gff3": "b6a1d46fa9548287611fd48eaf38ab37a829477446e711c7c6fea20a95eb815e",
    "bakta_like.json": "2ef3d4b7faf94bcc864a9743ecad4470d2f17724769e1585afc445c9737c4503",
    "genome.fasta": "d5480626271116e811665a0e06b35739a70ce063657bd277fc65e611582e5bde",
    "precomputed/operon_pairs.tsv": "1edbfd598b1332fcdb8fb908203049da0cb303aa0fe2b4383db39ab51200a2f2",  # noqa: E501
    "precomputed/promoters.tsv": "1014b150c5134312f6dcf55ac5fcc9bdf90dcbd6d8f6e342565c07ab3928aabe",
    "precomputed/rbs.tsv": "5c95152e264c554ae266ab47658897f49b1e34091ae6b910280b7510541ff823",
    "precomputed/rbs_novel.tsv": "6d3a61763eeca93724b3171e84556afa2426bc359e2d7a0abb9e452774cfa63a",
    "precomputed/regions.tsv": "f5718b0001a9bc220c9151b1ec3b193091530fcd472731d9659925e937733cb1",
    "precomputed/terminators.tsv": "ff5cd2b0becc62ad90c2fa09745af7fa63374d8cff60a89c5cc319b08897f4b9",  # noqa: E501
    "precomputed/tx_rates.tsv": "9590a6f0d19aa5f3f820b12e1c8b800ea0bdfc68f603943721efa8e4ad0cdc29",
    "truth/features.tsv": "6f8c6ccea5fcc6e59c713dfdd581c24ad8da2006dcf1a35ae64adc8c877b93f7",
    "truth/genes.tsv": "8be1b09fd2932c65366378b0201c45cbafa64b20c452731cbd0a9e3931791932",
    "truth/novel_orfs.tsv": "1b44f1d5f647a9bcf63df499a150921291d720abc96d08c02ccbea43a8c0020b",
    "truth/promoters.tsv": "9d7f3aec059d1b9e962703025ce810ccb29f2fc137dc72850e30e8b068e4087b",
    "truth/rbs.tsv": "ddbf0cc999430ba5d2c6859131ced89f2c1ecd8a2deb12fb275405e581286d37",
    "truth/regions.tsv": "8c31320cff1a63dbf46cc5bed7dff2d354035ba8f869175f3142036e51039fc1",
    "truth/terminators.tsv": "fa66d450a9944ed6e44845f6068347720c665d29e4e24cfb5a0edb6b69ca146c",
    "truth/tus.tsv": "a5d19dc08acab6976d9ff0483f86304e8d93f6e727f4ac930a53d8ff69368a2f",
}  # fmt: skip


def sc() -> Any:
    import giapsim.scenarios

    return giapsim.scenarios


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def files(d: Path) -> dict[str, str]:
    return {str(p.relative_to(d)): sha(p) for p in sorted(d.rglob("*"))
            if p.is_file() and p.name != "config.yaml"}  # fmt: skip


@dataclasses.dataclass
class Out:
    dir: Path
    entries: dict[str, Q.Entry]
    index: Q.EntryIndex


@pytest.fixture(scope="module")
def outs(tmp_path_factory: pytest.TempPathFactory) -> Any:
    cache: dict[tuple[int, str], Out] = {}

    def make(seed: int = SEED, size: str = "100k") -> Out:
        if (seed, size) not in cache:
            d = tmp_path_factory.mktemp(f"r3_{seed}_{size}")
            giapsim.simulate(sc().novel_orf_power(seed, size), seed, d)
            ent = Q.makedb_entries(d, tmp_path_factory.mktemp(f"r3db_{seed}_{size}"))
            cache[(seed, size)] = Out(d, ent, Q.EntryIndex(ent, Q.lengths(d)))
        return cache[(seed, size)]

    return make


def classes4(out: Path) -> dict[str, str]:
    """Decision 22 territory of every six-frame wrong peptide (oracle ``wrong_class4``)."""
    genome = O.receiving_genome(out)
    cds, genes = R.cds_spans(out), Q.true_starts(out)
    return {p: Q.wrong_class4(s, cds, genes, genome) for p, s in Q.wrong_loci(out).items()}


# ---------------------------------------------------------------- decision 21: decoy IDs


def test_sixframe_decoys_carry_makedb_class_tagged_ids(outs: Any) -> None:
    o = outs()
    _, d6, dp = Q.split_rows(o.dir)
    assert d6, "no rev_GIAPORF_ decoy rows in six-frame mode"
    for r in d6:
        for pid in Q.ids(r["protein_ids"]):
            assert Q.DECOY_ID.match(pid), r
    assert not any("rev_sixframe" in r["protein_ids"] for r in [*d6, *dp])


def test_sixframe_decoys_name_their_makedb_entry(outs: Any) -> None:
    """Stripped of ``rev_``, the ID is a makedb ``db_id`` whose protein holds the reversed
    peptide; the class in the ID is the entry's makedb class."""
    o = outs()
    _, d6, _ = Q.split_rows(o.dir)
    assert d6
    for r in d6:
        for pid in Q.ids(r["protein_ids"]):
            db_id = pid[len("rev_") :]
            assert db_id in o.entries, r
            e = o.entries[db_id]
            m = Q.DECOY_ID.match(pid)
            assert m is not None and m.group(1) == e.cls, (r, e)
            assert r["peptide"][::-1] in e.aa, (r, e)


def test_every_wrong_target_lies_in_one_makedb_entry(outs: Any) -> None:
    o = outs()
    loci = Q.wrong_loci(o.dir)
    assert len(loci) >= 20
    for pep, s in loci.items():
        hits = o.index.entry_of(s.contig, s.s1, s.e1, s.strand)
        assert len(hits) == 1, (pep, s, hits)
        assert hits[0].cls in Q.CLASSES4, (pep, hits[0])


def test_decoy_classes_cover_all_four(outs: Any) -> None:
    o = outs()
    _, d6, _ = Q.split_rows(o.dir)
    got = {m.group(1) for r in d6 for p in Q.ids(r["protein_ids"])
           if (m := Q.DECOY_ID.match(p))}  # fmt: skip
    assert got == set(Q.CLASSES4), got


# ---------------------------------------------------------------- decision 21: decoy scores


def _mixture(o: Out) -> tuple[float, tuple[float, float, float, float]]:
    pj = json.loads((o.dir / "truth" / "proteomics.json").read_text())
    mix = tuple(float(x) for x in sc().novel_orf_power(SEED).proteomics.score_mixture)
    return float(pj["pi_wrong_realised"]), mix  # type: ignore[return-value]


def test_sixframe_decoy_scores_finite_and_from_the_mixture(outs: Any) -> None:
    o = outs()
    pi, mix = _mixture(o)
    t, d6, _ = Q.split_rows(o.dir)
    assert d6 and t
    for r in [*t, *d6]:  # targets: the oracle's own guard
        s, q, p = (float(r[c]) for c in ("score", "q_value", "pep"))
        assert np.isfinite([s, q, p]).all(), r
        assert 0.0 <= q <= 1.0 and 0.0 <= p <= 1.0, r
        assert p == pytest.approx(Q.mixture_pep(s, pi, mix), rel=1e-6, abs=1e-12), r
        assert q == pytest.approx(Q.mixture_fdr(s, pi, mix), rel=1e-6, abs=1e-12), r


def test_sixframe_decoy_scores_follow_the_null(outs: Any) -> None:
    o = outs()
    _, _, mu0, sd0 = _mixture(o)[1]
    t, d6, _ = Q.split_rows(o.dir)
    ds = np.array([float(r["score"]) for r in d6])
    ts = np.array([float(r["score"]) for r in t])
    assert len(ds) >= 30, len(ds)
    assert stats.kstest(ds, stats.norm(mu0, sd0).cdf).pvalue >= KS_MIN_P
    assert stats.ks_2samp(ds, ts).pvalue >= KS_MIN_P


def test_null_psms_split_between_targets_and_decoys(outs: Any) -> None:
    """Single replicate: both sides present and the split within a two-sided binomial 0.001
    (the pooled calibration check is the real test)."""
    t, d6, _ = Q.split_rows(outs().dir)
    n = len(t) + len(d6)
    assert t and d6
    assert stats.binomtest(len(d6), n, 0.5).pvalue >= 0.001, (len(d6), n)


# ---------------------------------------------------------------- decision 21: gene level


def test_protein_decoys_stay(outs: Any) -> None:
    o = outs()
    prot = {k: f["aa"] for k, f in O.bakta_features(o.dir / "bakta_like.json").items()}
    _, _, dp = Q.split_rows(o.dir)
    assert dp, "no rev_<locus> protein decoys"
    for r in dp:
        for pid in Q.ids(r["protein_ids"]):
            assert pid.startswith("rev_") and pid[4:] in prot, r
            assert r["peptide"][::-1] in prot[pid[4:]], r


def test_gene_tier_is_mayu_while_giapsim_decoys_are_unscored(outs: Any) -> None:
    """GIAP's own gene-level pairing on ``peptides.tsv``: the ``rev_<locus>`` decoys pair, but
    giapsim writes them with no q, PEP or score, so M4-24 rule 4 gives ``mayu`` +
    ``decoy_unscored``. The giapsim scored-null-decoys backlog item restores ``decoy_picked``."""
    from giap.config.model import InputConfig
    from giap.genetable.detection import choose_tier, pair_decoys

    o = outs()
    pep = pd.read_csv(o.dir / "peptides.tsv", sep="\t", na_values=["NA"], keep_default_na=False)
    pep["is_decoy"] = pep["is_decoy"].astype(str).str.lower() == "true"
    tags = list(O.bakta_features(o.dir / "bakta_like.json"))
    pmap = pd.DataFrame({"protein_id": tags, "gene": tags, "method": "locus_tag"}).astype(object)
    ic = InputConfig()
    pairing = pair_decoys(pep, pmap, ic.decoy_prefixes, q_value_max=ic.q_value_max,
                          contaminant_prefixes=ic.contaminant_prefixes)  # fmt: skip
    empty = pd.DataFrame()
    tier, warns = choose_tier(pep, empty, empty, pairing)
    assert (tier, warns) == ("mayu", ["decoy_unscored"]), (tier, warns)


# ---------------------------------------------------------------- decision 22: nterm territory


def test_four_territories_split_evenly_100k(outs: Any) -> None:
    by: dict[str, int] = {c: 0 for c in Q.CLASSES4}
    for c in classes4(outs().dir).values():
        by[c] += 1
    n = sum(by.values())
    assert n >= 20, by
    for c in Q.CLASSES4:
        assert by[c] >= n // 4 - 2, by


def test_nterm_territory_is_giap_nterm_extension(outs: Any) -> None:
    """GIAP's M3-14 classifier calls the nterm-territory wrong peptides ``nterm_extension``."""
    from giap.annotate.bakta import parse_bakta_json
    from giap.mapping.loci import LocusClassifier
    from giap.model.coords import Interval

    o = outs()
    genome = O.receiving_genome(o.dir)
    lc = LocusClassifier(parse_bakta_json(o.dir / "bakta_like.json", genome=genome).model, genome)
    loci = Q.wrong_loci(o.dir)
    nterm = [p for p, c in classes4(o.dir).items() if c == "nterm_extension"]
    assert nterm
    for p in nterm:
        s = loci[p]
        cls, _ = lc.classify(Interval(s.contig, s.s1 - 1, s.e1, s.strand))  # type: ignore[arg-type]
        assert cls == "nterm_extension", (p, s, cls)


def test_nterm_territory_lies_upstream_of_the_true_start(outs: Any) -> None:
    """Choice 2's exclusions hold except the in-frame upstream extension itself: no nterm product
    overlaps its gene's true protein span (true start to stop) in frame."""
    o = outs()
    loci = Q.wrong_loci(o.dir)
    coding = R.coding_spans(o.dir)  # annotated CDS/sORFs and every novel_orfs row (true starts)
    for p, c in classes4(o.dir).items():
        if c == "nterm_extension":
            assert not any(R.in_frame_with(loci[p], g) for g in coding), p


@pytest.mark.slow
def test_sixframe_false_candidates_four_classes_1mb(outs: Any) -> None:
    """The round-2 per-class check extended to the fourth class (decision 22)."""
    by: dict[str, int] = {}
    for c in classes4(outs(SEED, "1M").dir).values():
        by[c] = by.get(c, 0) + 1
    for c in Q.CLASSES4:
        assert by.get(c, 0) >= 5, by


# ---------------------------------------------------------------- determinism


@pytest.mark.parametrize(("name", "want"), [("variant_f1", VARIANT_F1_SHA256),
                                           ("isobaric_attack", ISOBARIC_SHA256)])  # fmt: skip
def test_non_sixframe_scenarios_byte_identical(name: str, want: dict[str, str],
                                               tmp_path: Path) -> None:  # fmt: skip
    giapsim.simulate(getattr(sc(), name)(SEED), SEED, tmp_path)
    assert files(tmp_path) == want


def test_novel_orf_power_shuffled_source_byte_identical(tmp_path: Path) -> None:
    cfg = sc().novel_orf_power(SEED)
    prot = dataclasses.replace(cfg.proteomics, wrong_source="shuffled")
    cfg = dataclasses.replace(cfg, proteomics=prot)
    giapsim.simulate(cfg, SEED, tmp_path)
    got = files(tmp_path)
    assert {k: got.get(k) for k in NOVEL_SHUFFLED_SHA256} == NOVEL_SHUFFLED_SHA256


def test_novel_orf_power_non_proteomics_files_unchanged(outs: Any) -> None:
    got = files(outs().dir)
    assert {k: got.get(k) for k in NOVEL_NONPROTEOMICS_SHA256} == NOVEL_NONPROTEOMICS_SHA256
