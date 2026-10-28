"""M11-16e review round (PR #478): behavioural pins for the run's genetic code at every call
site the white-box keyword check (``test_m11_16e_downstream.test_call_sites_pass_the_run_table``)
cannot tell from a hard-wired ``table=11`` (docs/specs/M11-16-translation-tables-4-25.md
section 6, AC 10-12).

Each test is built so that reading the site at table 11 instead of the run's 4 or 25 changes an
observable output:

* 6.1 codon usage: rare codons 2-10 where TGG is a minority Trp codon at 4 (it is the only one at
  11); the step's ENC/tAI move with the table;
* 6.3 novel ORFs / class FDR: ``makedb_id_map`` tags the ``nterm_extension`` whose upstream walk
  crosses a TGA; ``class_fdr`` pairs it (picked regime); a ``novel_orfs`` step run at 4 writes W
  at the TGA, walks a position-less ORF up past a TGA and classifies a de novo decoy upstream of a
  gene across a TGA as ``nterm_extension`` (not in the intergenic list);
* 6.4 variants: TTA -> W is one base (TGA) at 4 and two (TGG) at 11; a left 2->1 block
  ``TGA GGT`` reads ``GG`` ~ N at 25; a frameshift whose prefix or suffix covers a TGA aligns;
* 6.5 makedb: a six-frame ORF upstream of a gene, X-delimited, with a TGA between it and the
  gene, is ``nterm_extension`` (left out) at 4 and ``intergenic`` at 11.

Reading pinned by the test author: a six-frame ORF the classifier calls ``nterm_extension`` is
not written (makedb keeps ``ORF_CLASSES`` only; the nterm records come from ``nterm_entries``).
"""

from __future__ import annotations

import math
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
from makedb_genomes import Contig, feature, stretch, write_bakta_json

from gippa.config import GippaConfig
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene
from gippa.model.genetic_code import codon_table

_SENSE_TABLES = (4, 25)


def tr(nt: str, table: int) -> str:
    ct = codon_table(table)
    return "".join(ct.get(nt[i : i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def il(p: str) -> str:
    return p.replace("I", "L")


def cfg_at(table: int, **kw: Any) -> GippaConfig:
    return GippaConfig.model_validate({"seed": 1, "annotation": {"translation_table": table}, **kw})


# ======================================================================== 6.1 codon usage


def cu() -> Any:
    import gippa.annotate.codon_usage as m

    return m


def _gene(tag: str, nt: str, **kw: Any) -> Gene:
    return Gene(tag, "cds", Interval("c1", 0, len(nt), "+"), "M", nt, nt[:3], **kw)


# TGA-rich genes: at 4 Trp is TGA 120 + TGG 9 (TGG share 0.07, rare); at 11 TGG is all of Trp.
_TGA_RICH = [_gene(f"T{k}", "ATG" + "TGAAAAGCC" * 40 + "CTG" + "TAA") for k in range(3)]
_RIBO = [
    _gene(f"R{k}", "ATG" + "GCCAAACTG" * 5 + "TAA", gene=f"rpl{chr(65 + k)}",
          product=f"50S ribosomal protein L{k + 1}")
    for k in range(3)
]  # fmt: skip
_TGG_EARLY = _gene("E1", "ATG" + "TGGAAAGCC" * 3 + "AAAGCC" * 4 + "TAA")  # TGG at codons 2, 5, 8


@pytest.mark.parametrize(("table", "want"), [(4, 3), (11, 0), (25, 0)])
def test_rare_codons_2_10_follow_the_table(table: int, want: int) -> None:
    """N2 (M27): TGG is rare at 4 (minority of the 2-codon Trp family), never at 11 or 25."""
    shares = cu().synonymous_shares(
        [g.nt[3:-3] for g in (*_TGA_RICH, *_RIBO, _TGG_EARLY)], table
    )  # the oracle: same eligible set as codon_usage
    assert ("TGG" in cu().rare_codons(shares)) is (table == 4)
    feats, _ = cu().codon_usage([*_TGA_RICH, *_RIBO, _TGG_EARLY], {"CCA": 1, "TCA": 1}, table)
    assert feats["E1"].rare_codons_2_10 == want


def test_ac11_tga_only_gene_reads_its_tga_weight() -> None:
    """B2 (M26): a TGA-only gene has tAI = w(TGA) > 0 at 4 (0.84 with tGCN TCA:1, CCA:1, measured
    with codon-bias 0.3.x on 2026-10-28), and no tAI at 11 (TGA is a stop, no sense codon)."""
    x = _gene("X", "ATG" + "TGA" * 12 + "TAA")
    trna = {"TCA": 1, "CCA": 1}
    f4, _ = cu().codon_usage([*_RIBO, x], trna, 4)
    f11, _ = cu().codon_usage([*_RIBO, x], trna, 11)
    assert f4["X"].tai is not None and f4["X"].tai == pytest.approx(0.8409, abs=1e-3)
    assert f11["X"].tai is None


# ================================================================== 6.3 / 6.5 shared world

G1_NT = "ATG" + stretch(8, 1) + "TGA" + stretch(8, 5) + "TAA"  # CDS [300, 357), TGA [327, 330)
# ORF A (+, [30, 93)): ATG stretch(8) TGA stretch(10, 3) TAA; TGA at [57, 60)
ORF_A_NT = "ATG" + stretch(8) + "TGA" + stretch(10, 3) + "TAA"
# ORF B (+, no position): TAA [120, 123) stretch(5, 1) TGA [138, 141) stretch(4, 2) TAA [153, 156)
ORF_B_NT = "TAA" + stretch(5, 1) + "TGA" + stretch(4, 2) + "TAA"
# Upstream of G1 in its frame: TAA [216, 219), stretch(20) [219, 279), NNN [279, 282),
# stretch(2) [282, 288), TGA [288, 291), stretch(3) [291, 300). The 20-aa six-frame ORF [219, 279)
# is X-delimited: the classifier's walk to G1 is blocked by the TGA at 11 only. G1's own upstream
# walk meets the X at 4/25 and the TGA at 11, so no nterm_entries record at any table.
UP_NT = "TAA" + stretch(20) + "NNN" + stretch(2, 4) + "TGA" + stretch(3, 6)


def _world_contig() -> Contig:
    c = Contig("c1", 600, "linear")
    c.put(30, ORF_A_NT)
    c.put(120, ORF_B_NT)
    c.put(216, UP_NT)
    c.put(300, G1_NT)
    return c


def _world_json(path: Path) -> tuple[Path, str]:
    c = _world_contig()
    js = write_bakta_json(path, [c], [feature("cds", "c1", 300, 357, "+", locus="G1", product="p")])
    return js, c.text


def _nterm_json(path: Path) -> Path:
    """``test_m11_16e_downstream._makedb_json``: upstream of G1 TAA [264, 267) + 7 codons + TGA
    [288, 291) + 3 codons, so the nterm_extension is [267, 354) at 4/25 and absent at 11."""
    c = Contig("c1", 600, "linear")
    c.put(99, stretch(10) + "TGA" + stretch(10, 3))
    c.put(264, "TAA" + stretch(7, 2) + "TGA" + stretch(3, 6))
    c.put(300, G1_NT)
    return write_bakta_json(path, [c], [feature("cds", "c1", 300, 357, "+", locus="G1",
                                                product="p")])  # fmt: skip


def _ann(js: Path, table: int) -> tuple[Any, dict[str, str]]:
    from gippa.annotate.bakta import parse_bakta_json

    ann = parse_bakta_json(js, translation_table=table)
    return ann.model, dict(ann.sequences)


NTERM_ID = "GIPPAORF_nterm_extension_c1_268_354_p"


def _pep_rows(target: str, ids: list[str], decoy: list[str] | None = None) -> pd.DataFrame:
    from helpers.novel_fdr_fixture import pep_row, peptides_frame

    rows = [pep_row(target, ids, 1, is_decoy=False)]
    if decoy is not None:
        rows.append(pep_row("DECYAAAK", decoy, 2, is_decoy=True))
    return peptides_frame(rows)


# ======================================================================== 6.3 class FDR


@pytest.mark.parametrize(("table", "tag"), [(4, "G1"), (25, "G1"), (11, None)])
def test_makedb_id_map_tags_the_nterm_entry_of_the_run_table(
    table: int, tag: str | None, tmp_path: Path
) -> None:
    """B3a (M12): the nterm_extension at [268, 354] exists only when TGA is sense."""
    from gippa.novel.class_call import makedb_id_map

    gm, genome = _ann(_nterm_json(tmp_path / "b.json"), table)
    m = makedb_id_map(_pep_rows("PEPAAAAK", [NTERM_ID]), gm, genome, ["rev_"], table=table)
    (r,) = m.to_dict("records")
    assert (r["db_id"], r["class"], int(r["start"]), int(r["end"])) == (NTERM_ID,
                                                                        "nterm_extension", 268,
                                                                        354)  # fmt: skip
    got = None if pd.isna(r["locus_tag"]) else r["locus_tag"]
    assert got == tag


@pytest.mark.parametrize("table", (4, 25, 11))
def test_class_fdr_pairs_the_nterm_orf_through_the_run_table(table: int, tmp_path: Path) -> None:
    """B3a (M12, M13): ``class_fdr`` at the run's table pairs an nterm_extension ORF with its
    makedb entry by locus tag; at 11 the entry has no tag, so the ORF is ``no_makedb_pair``."""
    from helpers.class_call_fixture import orf_at

    from gippa.novel.class_call import class_fdr

    gm, genome = _ann(_nterm_json(tmp_path / "b.json"), table)
    peps = _pep_rows("PEPAAAAK", [NTERM_ID], ["rev_" + NTERM_ID])
    orf = orf_at("NOVEL_c1_268_p", "nterm_extension", "c1", "+", 354, ["PEPAAAAK"], start=267,
                 overlapping_gene="G1")  # fmt: skip
    res = class_fdr([orf], peps, gm=gm, genome=genome, cfg=cfg_at(table), denovo=None)
    r = res.orfs[orf.orf_id]
    assert r.regime == "picked_decoy"
    if table == 11:
        assert r.reason == "no_makedb_pair"
    else:
        assert r.reason != "no_makedb_pair"
        assert [s.cls for s in res.classes] == ["nterm_extension"]


# ---------------------------------------------------------------- novel_orfs step at 4

_PEP_A = il(tr(stretch(8)[-12:] + "TGA" + stretch(10, 3)[:12], 4))  # 9 aa, W at index 4
_DECOY_IV = Interval("c1", 222, 249, "+")  # in G1's frame, upstream of the X and the TGA


def _pep_b(table: int) -> str:
    """ORF B's peptide: 5 aa with the table's TGA residue at index 2."""
    return il(tr(stretch(5, 1)[-6:] + "TGA" + stretch(4, 2)[:6], table))


def _orfs(table: int = 4) -> list[Any]:
    from gippa.novel.candidates import Orf

    a = Orf(orf_id="NOVEL_c1_31_p", cls="intergenic", contig="c1", strand="+", frame=0,
            stop=90, starts=(30,), start_longest=30, end_5p_open=False, end_3p_open=False,
            peptides=(_PEP_A,), best_q=None, overlapping_gene=None,
            first_peptide=45)  # fmt: skip
    b = Orf(orf_id="NOVEL_c1_124_p", cls="intergenic", contig="c1", strand="+", frame=0,
            stop=153, starts=(), start_longest=None, end_5p_open=False, end_3p_open=False,
            peptides=(_pep_b(table),), best_q=None, overlapping_gene=None,
            first_peptide=None)  # fmt: skip
    return [a, b]


def _denovo_frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    from helpers.class_call_fixture import costs_frame, decoy_hits, target_hits
    from helpers.novel_fdr_fixture import pep_row, peptides_frame

    peps = peptides_frame([pep_row(_PEP_A, [], 0, is_decoy=False, source="denovo", q=math.nan)])
    costs = costs_frame([_PEP_A], e_target=1e-4, e_decoy=1e-5, pi=0.01)
    hits = pd.concat([target_hits([_PEP_A], [Interval("c1", 45, 72, "+")]),
                      decoy_hits([_PEP_A], [_DECOY_IV])], ignore_index=True)  # fmt: skip
    return peps, costs, hits


@pytest.fixture(scope="module")
def novel_world(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """parse_bakta work files of the shared world, the candidates, a gene table and the de novo
    peptides / hits / costs of ORF A."""
    from helpers.class_call_fixture import gene_table

    from gippa.io.parquet import write_parquet_deterministic
    from gippa.steps import load_step_modules, registry
    from gippa.steps.novel import _CAND_COLS, _candidate_records

    load_step_modules()
    tmp = tmp_path_factory.mktemp("novel16e")
    js, seq = _world_json(tmp / "bakta.json")
    fa = tmp / "genome.fa"
    fa.write_text(f">c1\n{seq}\n")
    names = {"genes": "genes.parquet", "features": "features.parquet",
             "replicons": "replicons.parquet", "rrna16s": "rrna16s.fa",
             "trna_counts": "trna_counts.tsv", "xrefs": "xrefs.parquet",
             "annotation_meta": "annotation_meta.json"}  # fmt: skip
    work = {k: tmp / "work" / v for k, v in names.items()}
    registry.run_step("parse_bakta", cfg_at(4), {"bakta_json": js, "genome": fa}, work)
    reps = {"c1": Replicon("c1", 600, "linear", "chromosome")}
    cands = {}
    for t in (*_SENSE_TABLES, 11):
        cands[t] = tmp / f"novel_candidates_{t}.parquet"
        write_parquet_deterministic(cands[t], _CAND_COLS, _candidate_records(_orfs(t), reps, {}))
    gt = tmp / "gene_table.parquet"
    gene_table().to_parquet(gt, index=False)
    peps, costs, hits = _denovo_frames()
    paths = {"peptides": tmp / "peptides.parquet", "placement_costs": tmp / "costs.parquet",
             "hits_exact": tmp / "hits_exact.parquet",
             "hits_tolerant": tmp / "hits_tolerant.parquet"}  # fmt: skip
    peps.to_parquet(paths["peptides"], index=False)
    costs.to_parquet(paths["placement_costs"], index=False)
    hits.to_parquet(paths["hits_exact"], index=False)
    hits.iloc[0:0].to_parquet(paths["hits_tolerant"], index=False)
    ins = {"gene_table": gt, "genes": work["genes"],
           "replicons": work["replicons"], "genome": fa}  # fmt: skip
    return {"ins": ins, "cands": cands, "denovo": paths, "seq": seq, "js": js, "fa": fa}


def _novel_run(world: dict[str, Any], table: int, out: Path, denovo: bool) -> dict[str, Path]:
    from gippa.steps import registry

    ins = {**world["ins"], "novel_candidates": world["cands"][table]}
    if denovo:
        ins.update(world["denovo"])
    outs = {k: out / k for k in registry.get("novel_orfs").outputs}
    registry.run_step("novel_orfs", cfg_at(table), ins, outs)
    return outs


def _faa(path: Path) -> dict[str, str]:
    from gippa.io.protein_fasta import read_protein_fasta

    return {r.id: r.sequence for r in read_protein_fasta(path)}


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_novel_orfs_step_writes_the_run_table_residue(
    table: int, novel_world: dict[str, Any], tmp_path: Path
) -> None:
    """B3b (M15, M16, M18): ORF A's protein reads TGA as the table's residue (W / G); the
    position-less ORF B walks up past its TGA to the TAA at [120, 123), so its peptide is
    found 3 codons into the 10-codon body: 7 aa (at 11 the walk stops at the TGA and misses it)."""
    outs = _novel_run(novel_world, table, tmp_path, denovo=False)
    want = tr(ORF_A_NT[:-3], table)
    assert want[9] == codon_table(table)["TGA"]
    faa = _faa(outs["novel_orfs_faa"])
    assert "M" + want[1:] in faa.values(), faa
    t = pd.read_csv(outs["novel_orfs"], sep="\t")
    by = {r["orf_id"]: r for r in t.to_dict("records")}
    assert int(by["NOVEL_c1_124_p"]["length_aa_longest"]) == 7
    assert int(by["NOVEL_c1_31_p"]["length_aa_longest"]) == 20


def test_novel_orfs_step_classifies_de_novo_decoys_at_the_run_table(
    novel_world: dict[str, Any], tmp_path: Path
) -> None:
    """B3b (M14): the de novo decoy at [222, 249) is upstream of G1 across a TGA, so at 4 it is
    ``nterm_extension`` (no intergenic decoy: ORF A's q is the no-decoy q) and at 11 it would be
    an intergenic decoy that outscores ORF A's target. Oracle: ``class_fdr`` on the same inputs
    with ``LocusClassifier(table=4)``; the ``table=11`` classifier must give another q."""
    from gippa.mapping.loci import LocusClassifier
    from gippa.novel import class_call as cc

    outs = _novel_run(novel_world, 4, tmp_path, denovo=True)
    side = pd.read_parquet(outs["novel_class_fdr"])
    got = {r["orf_id"]: r for r in side.to_dict("records")}["NOVEL_c1_31_p"]
    assert got["regime"] == "denovo_subgroup"
    gm, genome = _ann(novel_world["js"], 4)
    assert LocusClassifier(gm, genome, table=4).classify(_DECOY_IV) == ("nterm_extension", "G1")
    assert LocusClassifier(gm, genome, table=11).classify(_DECOY_IV)[0] == "intergenic"
    peps, costs, hits = _denovo_frames()

    def q_with(classifier_table: int) -> float:
        dn = cc.denovo_inputs(costs, hits, LocusClassifier(gm, genome, table=classifier_table))
        res = cc.class_fdr(_orfs(), peps, gm=gm, genome=genome, cfg=cfg_at(4), denovo=dn)
        return float(res.orfs["NOVEL_c1_31_p"].q)

    assert q_with(4) != q_with(11)
    assert got["q"] == pytest.approx(q_with(4))


# ======================================================================== 6.5 makedb


@pytest.mark.parametrize(("table", "present"), [(4, False), (25, False), (11, True)])
def test_makedb_classifier_reads_the_run_table(table: int, present: bool, tmp_path: Path) -> None:
    """B3c (M22): the X-delimited six-frame ORF [219, 279) is ``nterm_extension`` of G1 at 4/25
    (TGA [288, 291) is sense) and so left out; at 11 the TGA blocks the walk: ``intergenic``."""
    from gippa.io.tsv import read_tsv
    from gippa.makedb.build import MakedbOptions, build_strain_db

    js, _ = _world_json(tmp_path / "bakta.json")
    out = tmp_path / "DB"
    build_strain_db(MakedbOptions(annotation=js, genome=None, outdir=out, six_frame_orfs=True,
                                  translation_table=table))  # fmt: skip
    header, rows = read_tsv(out / "db_id_map.tsv")
    recs = [dict(zip(header, r, strict=True)) for r in rows]
    at = [r for r in recs if r["kind"] == "orf" and (r["start"], r["end"]) == ("220", "279")]
    if present:
        assert [(r["class"], r["strand"]) for r in at] == [("intergenic", "+")]
    else:
        assert at == []
    assert not [r for r in recs if r["class"] == "nterm_extension"]  # no nterm_entries record


# ======================================================================== 6.4 variants

BODY, FLANK = "MSDKPLTEAGRWQLNEFK", "VPAGHTLDMRK"
N_TOY = 3 * (len(BODY) + 1 + len(FLANK) + 1)  # codon 18 [54, 57), terminal TAA [90, 93)


def _toy(codon18: str, letter: str, edits: dict[int, str]) -> tuple[dict[str, str], Any]:
    """The M7-16c toy with codon 18 = ``codon18`` (gene aa ``letter`` there) and codon edits."""
    from helpers.readthrough_toys import back_translate, gene, model

    nt = back_translate(BODY + "*" + FLANK + "*", {18: codon18, 30: "TAA"})
    for k, c in edits.items():
        nt = nt[: 3 * k] + c + nt[3 * k + 3 :]
    aa = tr(nt, 11 if letter == "U" else 4)  # only the letter at 18 depends on the table
    aa = aa[:18] + letter + aa[19:30]
    reps = {"c1": Replicon("c1", N_TOY, "linear")}
    g1 = replace(gene("G1", "c1", 0, N_TOY, "+"), aa=aa)
    return {"c1": nt}, model(reps, [g1])


def _run(genome: dict[str, str], gm: Any, pep: str, start: int, n_diff: int, table: int) -> Any:
    from gippa.steps.variants import run_variants

    hits = pd.DataFrame([{
        "peptide": pep, "peptide_len": len(pep), "source": "denovo", "contig": "c1",
        "start": start, "end": start + 3 * len(pep), "strand": "+", "n_diff": n_diff,
        "diff_types": "substitution" if n_diff else pd.NA, "placement_q": 0.0,
        "uniqueness": "unique",
    }])  # fmt: skip
    peptides = pd.DataFrame({"peptide": [pep], "search_db": "same"})
    cfg = cfg_at(table, genome_relation="same", sample_name="S1")
    return run_variants(hits, None, peptides, gm, genome, cfg)


@pytest.mark.parametrize(
    ("table", "codon18", "alt_codon", "n", "flag"),
    [(4, "TGA", "TGA", 1, "snv"), (25, "TGA", "TGG", 2, "multi_nt"),
     (11, "TGG", "TGG", 2, "multi_nt")],
)  # fmt: skip
def test_tta_to_w_minimal_codon_follows_the_table(
    table: int, codon18: str, alt_codon: str, n: int, flag: str
) -> None:
    """N1 (M5, M8): ref TTA (L, codon 13 [39, 42)) with ALT W: TGA (1 nt) at 4, TGG (2 nt) at
    11 and 25; ``nt_class`` agrees (``snv`` / ``multi_nt``)."""
    from gippa.variants.flags import nt_class
    from gippa.variants.vcf_records import minimal_alt_codon

    assert minimal_alt_codon("TTA", "W", table=table) == (alt_codon, n)
    assert nt_class("TTA", "W", table=table) == flag
    letter = codon_table(table)[codon18]
    genome, gm = _toy(codon18, letter, {13: "TTA"})
    pep = "WQWNEFK" + letter + "VPAGHTL"
    _pile, t, records = _run(genome, gm, pep, 33, 1, table)
    (r,) = t.to_dict("records")
    assert (r["codon_start"], r["ref_codon"], r["ref_aa"], r["alt_aa"]) == (39, "TTA", "L", "W")
    assert r["alt_codon_vcf"] == alt_codon and int(r["nt_changes"]) == n
    assert flag in set(r["flags"]) and ({"snv", "multi_nt"} - {flag}).isdisjoint(r["flags"])
    (rec,) = records
    assert (rec.pos0, rec.ref, rec.alts) == (39, "TTA", (alt_codon,))


@pytest.mark.parametrize(("table", "letter", "mass_eq"), [(25, "G", True), (11, "U", False)])
def test_left_block_reads_the_tga_before_the_site_with_the_table(
    table: int, letter: str, mass_eq: bool
) -> None:
    """N1 (M1, M9): site codon 19 GGT (G) with ALT Q, the codon before it TGA. At 25 the left
    2->1 block ``TGA GGT`` = GG + Methyl (128.058 Da) ~ Q: ``mass_equivalent``; at 11 it is
    ``*G``. G -> Q itself (+71.037 Da) is no B.5 delta, so only the block read flags it.
    (GG ~ N would not do: N - G is the Carbamidomethyl delta, so G -> N is mass_equivalent at
    any table.)"""
    from gippa.variants.flags import left_block_read, mass_equivalent

    tol = GippaConfig().variants.mass_equivalent_tol_da
    assert not mass_equivalent("G", "Q", tol)
    assert left_block_read("QPAGHTL", 0, "GG", 1, tol)
    assert not left_block_read("QPAGHTL", 0, "*G", 1, tol)
    genome, gm = _toy("TGA", letter, {19: "GGT"})
    _pile, t, _records = _run(genome, gm, "QPAGHTLDMRK", 57, 1, table)
    (r,) = t.to_dict("records")
    assert (r["codon_start"], r["ref_codon"], r["ref_aa"], r["alt_aa"]) == (57, "GGT", "G", "Q")
    assert ("mass_equivalent" in set(r["flags"])) is mass_eq


@pytest.mark.parametrize("table", _SENSE_TABLES)
@pytest.mark.parametrize("side", ("prefix", "suffix"))
def test_frameshift_block_covering_a_tga_aligns_at_4_and_25(table: int, side: str) -> None:
    """N1 (M10): a +1 frameshift whose prefix (or suffix) block covers a TGA reads it as the
    table's residue; at 11 the ``*`` outside ``stop_sites`` is refused."""
    from gippa.variants.align import align_hit

    a, b = "ATGAAATGAGCC", "GAAAAACTGCGT"  # MKWA / MKGA ; EKLR
    if side == "suffix":
        a, b = b, a
    seq = a + "C" + b + "CCCCCC"
    rep = Replicon("c", len(seq), "linear", "chromosome")
    pep = tr(a, table) + tr(b, table)
    hit = align_hit(pep, "c", 0, 25, "+", 1, "0-12;13-25", seq, rep, table=table)
    assert hit.frameshift is not None and hit.frameshift.shift == 1
    assert all(o.kind in ("match", "il") for o in hit.ops)
    assert "".join(o.ref for o in hit.ops) == pep
    with pytest.raises(ValueError, match=r"X or \*"):
        align_hit(pep, "c", 0, 25, "+", 1, "0-12;13-25", seq, rep, table=11)
