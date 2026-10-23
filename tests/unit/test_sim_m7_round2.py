"""giapsim M7-12 round 2: lead decisions 11-16 after the statistician's block on PR #225
(docs/specs/M7-12-sim-variant-novel-scenarios.md; GIAP_SPEC.md §7.5.6, §7.6.6, §8.4, §11.3).

The oracle is tests/helpers/m7_round2_oracle.py (no giapsim import). The 1-Mb cases are ``slow``.

Interface choices pinned here (names follow the existing giapsim API):

* **Decision 11.** ``truth/alt_sites.tsv`` keeps its seven columns first and adds
  ``frameshift_upstream`` (``true``/``false``). Each row is the SNV codon in the *sample* gene's
  frame lifted to ``reference.fasta`` (``codon_start = cs_sample + pos_reference - pos_sample``);
  ``ref_aa`` / ``alt_aa`` are the reference / sample translations of that codon (code 11); rows
  synonymous or nonsense in that frame are dropped. ``frameshift_upstream`` is false when no
  frameshifting indel of the same gene lies upstream (transcription direction) of the SNV and true
  when the net length change of the gene's upstream indels is not a multiple of 3 (the
  net-zero case with frameshifts is left open). ``truth/indel_sites.tsv`` now lists **every**
  ins/del inside a gene (in-frame and frameshifting), adding ``frameshift`` (bool) and
  ``junction`` (int: 1-based reference coordinate of the first reference base after the VCF
  anchor, ``pos_reference + 1``). ``giapsim.truth.read_alt_sites`` / ``read_indel_sites`` parse
  the new columns as bool / int.
* **Decision 12.** ``ProteomicsConfig.wrong_source = "sixframe_noncoding"`` (opt-in; the default
  stays ``"shuffled"``): wrong peptides are Trypsin/P products of six-frame translations of
  ``genome.fasta`` that are in frame with no annotated CDS/sORF and no ``truth/novel_orfs.tsv``
  ORF and occur in no true protein (I = L). Their ``truth/peptides.tsv`` rows have
  ``is_wrong = true``, ``origin = sixframe_noncoding`` and the source locus in ``contig, start,
  end, strand`` (1-based closed), ``n_loci_genome >= 1``; they get null-component scores, and
  ``n_decoy_rows = "match_wrong"`` gives one decoy row per wrong target row. ``novel_orf_power``
  turns it on; ``variant_f1`` and ``isobaric_attack`` keep ``"shuffled"``.
* **Decision 13.** ``giapsim.scenarios.SIZES["300k"] = 300_000`` and
  ``NOVEL_PER_CLASS["300k"]``; ``novel_orf_power(seed, "300k")`` has >= 200 annotated
  non-pseudo CDSs (so GIAP's decoy-picked detection floor 1/n is below q = 0.01).
* **Decision 14.** ``SimConfig.planted_tir: tuple[float, float]`` = (mu, sigma) of the planted
  starts' log TIR in ``rbs_novel``; ``novel_orf_power`` sets (8.0, 0.8) at every size. Changing it
  changes only those TIRs.
* **Decision 15.** ``giapsim.features.novel_orfs`` uses ``min_gap = 3 * (hi + 1) + 2 * margin + 1``
  (still at least 160) and, with fewer usable gaps than requested intergenic ORFs, plants one per
  gap and logs a WARNING mentioning "intergenic"; ``truth/novel_orfs.tsv`` then holds the realised
  rows (power denominators read the truth).
* **Decision 16.** ``truth/novel_orfs.tsv`` ``stop`` is the 1-based lowest coordinate of the stop
  codon (``end - 2`` on ``+``, ``start`` on ``-``).
"""

from __future__ import annotations

import dataclasses
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers import m7_round2_oracle as R
from helpers import m7_sim_oracle as O

import giapsim

SEED = 20261022


def sc() -> Any:
    import giapsim.scenarios

    return giapsim.scenarios


@pytest.fixture(scope="module")
def outs(tmp_path_factory: pytest.TempPathFactory) -> Any:
    cache: dict[tuple[str, int, str], Path] = {}

    def make(name: str, seed: int = SEED, size: str = "100k") -> Path:
        key = (name, seed, size)
        if key not in cache:
            d = tmp_path_factory.mktemp(f"r2_{name}_{seed}_{size}")
            fn = getattr(sc(), name)
            cfg = fn(seed, size) if name in ("variant_f1", "novel_orf_power") else fn(seed)
            giapsim.simulate(cfg, seed, d)
            cache[key] = d
        return cache[key]

    return make


def alt_key(r: dict[str, str]) -> R.AltKey:
    return (r["contig"], int(r["codon_start"]), r["strand"], r["gene"], r["ref_aa"], r["alt_aa"])


def check_sample_frame(out: Path) -> R.FrameTruth:
    """Decision 11 oracle: alt_sites (outside ambiguous codons) == the sample-frame lift."""
    ft = R.sample_frame_alt_sites(out)
    got = [r for r in O.rows(out / "truth" / "alt_sites.tsv")
           if (r["contig"], int(r["codon_start"]), r["gene"]) not in ft.ambiguous]  # fmt: skip
    want = {k for k in ft.want if (k[0], k[1], k[3]) not in ft.ambiguous}
    keys = {alt_key(r) for r in got}
    assert len(keys) == len(got), "one row per codon"
    assert keys == want, (sorted(keys - want)[:5], sorted(want - keys)[:5])
    return ft


# ---------------------------------------------------------------- decision 11


@pytest.mark.parametrize("seed", [SEED, 7])
def test_alt_sites_follow_sample_frame_lift(outs: Any, seed: int) -> None:
    """Seed 7 holds the statistician's example (truth 42592 N>V; peptides support 42593 M>V)."""
    check_sample_frame(outs("variant_f1", seed))


def test_statistician_example_seed7(outs: Any) -> None:
    out = outs("variant_f1", 7)
    rows = {alt_key(r): r for r in O.rows(out / "truth" / "alt_sites.tsv")}
    assert ("chr1", 42593, "+", "GSIM_00205", "M", "V") in rows
    assert ("chr1", 42592, "+", "GSIM_00205", "N", "V") not in rows
    assert O.truthy(rows[("chr1", 42593, "+", "GSIM_00205", "M", "V")]["frameshift_upstream"])


def test_alt_sites_header_gains_frameshift_upstream(outs: Any) -> None:
    for name in ("variant_f1", "isobaric_attack"):
        h = O.header(outs(name) / "truth" / "alt_sites.tsv")
        assert h[:7] == O.ALT_COLUMNS and "frameshift_upstream" in h[7:], (name, h)


@pytest.mark.parametrize("seed", [SEED, 7])
def test_frameshift_upstream_flag(outs: Any, seed: int) -> None:
    out = outs("variant_f1", seed)
    ft = R.sample_frame_alt_sites(out)
    n_true = 0
    for r in O.rows(out / "truth" / "alt_sites.tsv"):
        assert r["frameshift_upstream"] in ("true", "false"), r
        key = (r["contig"], int(r["codon_start"]), r["gene"])
        if key in ft.ambiguous or key not in ft.shift:
            continue
        any_fs, net = ft.shift[key]
        flag = O.truthy(r["frameshift_upstream"])
        n_true += flag
        if not any_fs:
            assert not flag, r
        if net != 0:
            assert flag, r
    if seed == 7:
        assert n_true >= 3, "seed 7 has three SNV codons downstream of a frameshift"


def test_alt_sites_counts_after_lift(outs: Any) -> None:
    """Rule 2 still holds on the lifted rows: peptides covering the SNV's sample position."""
    out = outs("variant_f1", 7)
    ft = R.sample_frame_alt_sites(out)
    sample = O.read_fasta(out / "genome.fasta")
    correct = O.correct_peptide_rows(out)
    for r in O.rows(out / "truth" / "alt_sites.tsv"):
        key = (r["contig"], int(r["codon_start"]), r["gene"])
        if key in ft.ambiguous or key not in ft.snv:
            continue
        s = ft.snv[key]
        n = O.recount_site_peptides(out, s["gene"], s["contig"], int(s["pos_sample"]),
                                    len(sample[s["contig"]]), correct)  # fmt: skip
        assert int(r["n_true_peptides"]) == n, (r, n)


@pytest.mark.parametrize("seed", [SEED, 7])
def test_indel_sites_list_frameshifts_with_junction(outs: Any, seed: int) -> None:
    out = outs("variant_f1", seed)
    path = out / "truth" / "indel_sites.tsv"
    h = O.header(path)
    assert {"frameshift", "junction"} <= set(h), h
    got = {(r["contig"], r["pos_reference"], r["ref"], r["alt"], r["gene"],
            O.truthy(r["frameshift"]), int(r["junction"])) for r in O.rows(path)}  # fmt: skip
    want = {(r["contig"], r["pos_reference"], r["ref"], r["alt"], r["gene"], r["frameshift"],
             r["junction"]) for r in R.indels_in_genes(out)}  # fmt: skip
    assert got == want, sorted(got ^ want)[:5]
    assert any(w[5] for w in want), "variant_f1 plants frameshifting indels in genes"


def test_truth_readers_parse_new_columns(outs: Any) -> None:
    from giapsim.truth import read_alt_sites, read_indel_sites

    out = outs("variant_f1", 7)
    alt = read_alt_sites(out)
    assert alt["frameshift_upstream"].dtype == bool
    assert alt["frameshift_upstream"].any()
    ind = read_indel_sites(out)
    assert ind["frameshift"].dtype == bool
    assert np.issubdtype(ind["junction"].dtype, np.integer)
    empty = read_alt_sites(outs("isobaric_attack"))
    assert empty.empty and "frameshift_upstream" in empty.columns


@pytest.mark.slow
def test_alt_sites_sample_frame_1mb(outs: Any) -> None:
    out = outs("variant_f1", 1, "1M")
    check_sample_frame(out)
    alt = O.rows(out / "truth" / "alt_sites.tsv")
    assert sum(O.truthy(r["frameshift_upstream"]) for r in alt) >= 10


# ---------------------------------------------------------------- decision 12


def wrong_rows(out: Path) -> list[dict[str, str]]:
    return [r for r in O.rows(out / "truth" / "peptides.tsv")
            if O.truthy(r["is_wrong"]) and r["origin"] != "decoy"]  # fmt: skip


def test_sixframe_source_is_opt_in() -> None:
    from giapsim.proteomics import ProteomicsConfig

    assert ProteomicsConfig().wrong_source == "shuffled"
    assert ProteomicsConfig(wrong_source="sixframe_noncoding").wrong_source == "sixframe_noncoding"
    with pytest.raises(ValueError):
        ProteomicsConfig(wrong_source="nonsense")  # type: ignore[arg-type]
    assert sc().novel_orf_power(SEED).proteomics.wrong_source == "sixframe_noncoding"
    assert sc().novel_orf_power(SEED, "1M").proteomics.wrong_source == "sixframe_noncoding"
    assert sc().variant_f1(SEED).proteomics.wrong_source == "shuffled"
    assert sc().isobaric_attack(SEED).proteomics.wrong_source == "shuffled"


def test_sixframe_wrong_peptides_map_to_noncoding_frames(outs: Any) -> None:
    out = outs("novel_orf_power")
    seqs = O.read_fasta(out / "genome.fasta")
    wrong = wrong_rows(out)
    assert len(wrong) >= 20, len(wrong)
    assert {r["origin"] for r in wrong} == {"sixframe_noncoding"}
    coding = R.coding_spans(out)
    hay = "|".join(O.il(p) for p in [*(f["aa"] for f in
                   O.bakta_features(out / "bakta_like.json").values() if f.get("aa")),
                   *(r["aa"] for r in O.rows(out / "truth" / "novel_orfs.tsv"))])  # fmt: skip
    for r in wrong:
        pep, c, st = r["peptide"], r["contig"], r["strand"]
        s1, e1 = int(r["start"]), int(r["end"])
        assert 7 <= len(pep) <= 30 and e1 - s1 + 1 == 3 * len(pep), r
        assert R.translate_nt(O.strand_slice(seqs[c], s1 - 1, e1, st)) == pep, r
        assert int(r["n_loci_genome"]) >= 1, r
        assert O.il(pep) not in hay, r  # in no true or annotated protein
        span = R.Span(c, s1, e1, st)
        assert not any(R.in_frame_with(span, g) for g in coding), r
        # Trypsin/P: cut after K/R, or the frame's stop
        up0, dn0 = (s1 - 4, e1) if st == "+" else (e1, s1 - 4)
        before, after = R.codon_aa(seqs[c], up0, st), R.codon_aa(seqs[c], dn0, st)
        if before is not None:
            assert before in "KR*", r
        if after is not None:
            assert pep[-1] in "KR" or after == "*", r


def test_sixframe_no_shortfall_and_matched_decoys(outs: Any) -> None:
    import json

    out = outs("novel_orf_power")
    pj = json.loads((out / "truth" / "proteomics.json").read_text())
    assert pj["n_wrong"] == pj["n_wrong_requested"] > 0, pj
    peps = O.rows(out / "peptides.tsv")
    wrong = {r["peptide"] for r in wrong_rows(out)}
    assert {r["origin"] for r in wrong_rows(out)} == {"sixframe_noncoding"}
    n_wrong_rows = sum(r["is_decoy"] == "false" and r["peptide"] in wrong for r in peps)
    n_decoy_rows = sum(r["is_decoy"] == "true" for r in peps)
    assert n_decoy_rows == n_wrong_rows > 0
    sw = [float(r["score"]) for r in peps if r["is_decoy"] == "false" and r["peptide"] in wrong]
    sc_ = [
        float(r["score"]) for r in peps if r["is_decoy"] == "false" and r["peptide"] not in wrong
    ]
    assert np.mean(sc_) - np.mean(sw) > 2.0  # wrong-ID (null) scores


@pytest.mark.slow
def test_sixframe_false_candidates_per_class_1mb(outs: Any) -> None:
    """The statistician's check: at 1 Mb, >= 5 wrong peptides in each novel class's territory."""
    out = outs("novel_orf_power", SEED, "1M")
    cds = R.cds_spans(out)
    by: dict[str, int] = {}
    for r in wrong_rows(out):
        cls = R.wrong_class(R.Span(r["contig"], int(r["start"]), int(r["end"]), r["strand"]), cds)
        by[cls] = by.get(cls, 0) + 1
    for cls in ("intergenic", "cds_outframe", "antisense"):
        assert by.get(cls, 0) >= 5, by


# ---------------------------------------------------------------- decision 13


def test_sim300k_size() -> None:
    assert sc().SIZES["300k"] == 300_000
    assert sc().novel_orf_power(SEED, "300k").genome_length == 300_000
    assert (
        sc().NOVEL_PER_CLASS["100k"] <= sc().NOVEL_PER_CLASS["300k"] <= sc().NOVEL_PER_CLASS["1M"]
    )


def test_sim300k_has_200_genes(outs: Any) -> None:
    out = outs("novel_orf_power", SEED, "300k")
    cds = [f for f in O.bakta_features(out / "bakta_like.json").values()
           if f["type"] == "cds" and not f.get("pseudo")]  # fmt: skip
    assert len(cds) >= 200, len(cds)
    assert 1 / len(cds) < 0.01
    rows = O.rows(out / "truth" / "novel_orfs.tsv")
    want = sc().NOVEL_PER_CLASS["300k"]
    for cls in O.NOVEL_CLASSES_M7:
        assert sum(r["class"] == cls for r in rows) == want, cls


# ---------------------------------------------------------------- decision 14


def test_planted_tir_is_a_config_field() -> None:
    fields = {f.name for f in dataclasses.fields(giapsim.SimConfig)}
    assert "planted_tir" in fields
    for size in ("100k", "300k", "1M"):
        assert tuple(sc().novel_orf_power(SEED, size).planted_tir) == (8.0, 0.8)


def test_planted_tir_changes_only_planted_tirs(outs: Any, tmp_path: Path) -> None:
    base = outs("novel_orf_power")
    cfg = dataclasses.replace(sc().novel_orf_power(SEED), planted_tir=(9.0, 0.5))
    giapsim.simulate(cfg, SEED, tmp_path)
    for f in ("genome.fasta", "truth/novel_orfs.tsv", "precomputed/rbs.tsv", "peptides.tsv"):
        assert (tmp_path / f).read_bytes() == (base / f).read_bytes(), f
    planted = {
        (r["contig"], r["strand"], int(r["true_start"]))
        for r in O.rows(base / "truth" / "novel_orfs.tsv")
        if O.truthy(r["rbs_planted"])
    }

    def layer(d: Path) -> dict[tuple[str, str, int], float]:
        return {(r["contig"], r["strand"], int(r["start"])): float(r["tir_ostir"])
                for r in O.rows(d / "precomputed" / "rbs_novel.tsv")}  # fmt: skip

    a, b = layer(base), layer(tmp_path)
    assert a.keys() == b.keys()
    diff = {k for k in a if a[k] != b[k]}
    assert diff and diff <= planted
    # stronger arm: higher mean log TIR at the planted starts
    assert np.mean([np.log(b[k]) for k in planted]) > np.mean([np.log(a[k]) for k in planted])


# ---------------------------------------------------------------- decision 15


def test_novel_orf_power_100k_seed14_does_not_crash(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.WARNING):
        giapsim.simulate(sc().novel_orf_power(14), 14, tmp_path)
    rows = O.rows(tmp_path / "truth" / "novel_orfs.tsv")
    n_inter = sum(r["class"] == "intergenic" for r in rows)
    assert 1 <= n_inter <= sc().NOVEL_PER_CLASS["100k"]
    if n_inter < sc().NOVEL_PER_CLASS["100k"]:
        assert any("intergenic" in m for m in caplog.messages), caplog.messages
    genome = O.receiving_genome(tmp_path)
    for r in rows:
        nt = O.strand_slice(genome[r["contig"]], int(r["start"]) - 1, int(r["end"]), r["strand"])
        assert O.translate_orf(nt) == r["aa"], r["orf_id"]


@pytest.mark.slow
def test_novel_orf_power_1mb_seed103_does_not_crash(tmp_path: Path) -> None:
    giapsim.simulate(sc().novel_orf_power(103, "1M"), 103, tmp_path)
    rows = O.rows(tmp_path / "truth" / "novel_orfs.tsv")
    assert sum(r["class"] == "intergenic" for r in rows) == sc().NOVEL_PER_CLASS["1M"]


def _gap_case(gap: int, hi: int = 60, margin: int = 20) -> tuple[list[list[Any]], int, int]:
    from giapsim.features import Feat, novel_orfs
    from giapsim.simulate import SimConfig, _Sim

    a = 300
    b = a + gap
    n = b + 300
    seqs = {"chr1": bytearray(b"A" * n)}
    feats = [Feat("cds", "chr1", 0, a, "+", "A" * a, "M" + "K" * 99),
             Feat("cds", "chr1", b, b + 300, "+", "A" * 300, "M" + "R" * 99)]  # fmt: skip
    sim = _Sim(SimConfig(), np.random.default_rng(0))
    rows = novel_orfs(sim, {"intergenic": 1}, seqs, feats, (hi, hi), margin)
    return rows, a, b


def test_min_gap_off_by_one_exact_gap_does_not_crash(caplog: pytest.LogCaptureFixture) -> None:
    """A gap of exactly 3 (hi + 1) + 2 margin made ``rng.integers(margin, margin)`` raise."""
    gap = 3 * 61 + 2 * 20
    with caplog.at_level(logging.WARNING):
        rows, a, b = _gap_case(gap)
    for r in rows:
        assert a + 20 <= int(r[3]) - 1 and int(r[4]) <= b - 20, r
    if not rows:
        assert any("intergenic" in m for m in caplog.messages), caplog.messages


def test_min_gap_one_more_plants_with_margins() -> None:
    gap = 3 * 61 + 2 * 20 + 1
    rows, a, b = _gap_case(gap)
    assert len(rows) == 1
    s1, e0 = int(rows[0][3]), int(rows[0][4])
    assert e0 - (s1 - 1) == 3 * 61
    assert a + 20 <= s1 - 1 and e0 <= b - 20


# ---------------------------------------------------------------- decision 16


def test_novel_stop_is_one_based(outs: Any) -> None:
    out = outs("novel_orf_power")
    genome = O.receiving_genome(out)
    rows = O.rows(out / "truth" / "novel_orfs.tsv")
    assert rows
    for r in rows:
        s1, e1, st, stop = int(r["start"]), int(r["end"]), r["strand"], int(r["stop"])
        assert stop == (e1 - 2 if st == "+" else s1), r
        assert O.strand_slice(genome[r["contig"]], stop - 1, stop + 2, st) in O.STOPS, r


def test_novel_reader_stop_one_based(outs: Any) -> None:
    from giapsim.truth import read_novel_orfs

    df = read_novel_orfs(outs("novel_orf_power"))
    plus = df[df["strand"] == "+"]
    minus = df[df["strand"] == "-"]
    assert (plus["stop"] == plus["end"] - 2).all()
    assert (minus["stop"] == minus["start"]).all()
