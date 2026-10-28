"""Genetic codes 4 and 25 downstream of mapping (M11-16e;
docs/specs/M11-16-translation-tables-4-25.md items 6.1-6.5, examples 10-12, AC 10 (makedb), 11
and 12).

Readings pinned by the test author where the spec names no signature:

* the run's table reaches pure functions as a keyword: ``table=`` (as 16b's ``annotated_stops``,
  ``LocusClassifier``) except where the first argument is already a ``table`` DataFrame
  (``novel.outputs.orf_table``, ``protein_records``, ``io.novel_tables.write_novel_faa``), which
  take ``translation_table=``. Functions that already take a ``GippaConfig`` (``apply_evidence``,
  ``flags_for``, ``pileup``, ``run_variants``) read ``cfg.annotation.translation_table``;
* ``MakedbOptions.translation_table: int = 11``; ``makedb_manifest.json`` gains a top-level
  ``"translation_table"``; the CLI fills it from ``--set annotation.translation_table=N``;
* example 10 "TGG is a 1-codon Trp family with no share" under 25 is not asserted (today a
  1-codon family gets share 1.0 at 11 too); ``scored_sequence`` has no internal-stop rule
  today, so "eligible" is tested through ``synonymous_shares`` and ``codon_usage`` scores;
* AC 11: codon-bias 0.3.x gives TGA a tAI weight of 0.84 from anticodon ``TCA`` at code 4 with
  tGCN ``{TCA: 1, CCA: 1}`` (0.37 with ``TRNA_4``; measured 2026-10-28), so the AC holds and D10
  needs no finding.

Table 11 stays byte-identical: ``MAKEDB_11`` pins today's ``strain_db.faa`` / ``db_id_map.tsv``
of the hand-built genome, computed on 16b (3cf0b509) before this change.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
from makedb_genomes import Contig, feature, stretch, write_bakta_json

from gippa.config import GippaConfig
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene
from gippa.model.genetic_code import codon_table

SRC = Path(__file__).resolve().parents[2] / "src" / "gippa"
_COMP = str.maketrans("ACGT", "TGCA")
_SENSE_TABLES = (4, 25)


def tr(nt: str, table: int) -> str:
    """Brute-force reference: codon by codon, ``X`` for a non-ACGT codon, partial codon dropped."""
    ct = codon_table(table)
    return "".join(ct.get(nt[i : i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def rc(nt: str) -> str:
    return nt.translate(_COMP)[::-1]


def cfg_at(table: int, **kw: Any) -> GippaConfig:
    return GippaConfig.model_validate({"seed": 1, "annotation": {"translation_table": table}, **kw})


# ======================================================================== 6.1 codon usage


def cu() -> Any:
    import gippa.annotate.codon_usage as m

    return m


def test_example10_table4_trp_is_a_two_codon_family() -> None:
    shares = cu().synonymous_shares(["TGGTGATGAAAA"], 4)
    assert shares["TGA"] == pytest.approx(2 / 3)
    assert shares["TGG"] == pytest.approx(1 / 3)
    assert shares["AAA"] == 1.0


def test_example10_table25_tga_counts_in_the_gly_family() -> None:
    shares = cu().synonymous_shares(["TGGTGATGAAAAGGC"], 25)
    assert shares["TGA"] == pytest.approx(2 / 3)
    assert shares["GGC"] == pytest.approx(1 / 3)
    assert shares["GGA"] == shares["GGG"] == shares["GGT"] == 0.0


def test_example10_table11_tga_is_a_stop_unchanged() -> None:
    """Pin (passes today): at 11 TGA never enters a family."""
    shares = cu().synonymous_shares(["TGGTGATGAAAA"], 11)
    assert "TGA" not in shares and shares["TGG"] == 1.0


def test_table4_trp_codon_can_be_rare() -> None:
    shares = cu().synonymous_shares(["TGA" * 19 + "TGG"], 4)
    assert cu().rare_codons(shares) >= {"TGG"}
    assert "TGA" not in cu().rare_codons(shares)


@given(st.lists(st.text("ACGT", min_size=0, max_size=60), max_size=5), st.sampled_from((4, 11, 25)))
@settings(max_examples=60, deadline=None)
def test_shares_match_a_brute_force_reference(seqs: list[str], table: int) -> None:
    ct = codon_table(table)
    counts: dict[str, int] = {}
    for s in seqs:
        for i in range(0, len(s) - len(s) % 3, 3):
            c = s[i : i + 3]
            if ct[c] != "*":
                counts[c] = counts.get(c, 0) + 1
    want: dict[str, float] = {}
    for aa in set(ct.values()) - {"*"}:
        fam = [c for c in ct if ct[c] == aa]
        tot = sum(counts.get(c, 0) for c in fam)
        if tot:
            want.update({c: counts.get(c, 0) / tot for c in fam})
    assert cu().synonymous_shares(seqs, table) == pytest.approx(want)


def _gene(tag: str, nt: str, **kw: Any) -> Gene:
    return Gene(tag, "cds", Interval("c1", 0, len(nt), "+"), "M", nt, nt[:3], **kw)


TGA_GENES = [
    _gene(f"T{k}", "ATG" + ("TGAAAAGCC" * (4 + k)) + "TGGCTG" + "TAA") for k in range(6)
] + [
    _gene(f"R{k}", "ATG" + ("GCCTGAAAA" * 5) + "CTG" + "TAA", gene=f"rpl{chr(65 + k)}",
          product=f"50S ribosomal protein L{k + 1}")
    for k in range(3)
]  # fmt: skip
TRNA_4 = {"TCA": 1, "CCA": 1, "TTT": 2, "GGC": 2, "CAG": 3, "CAT": 3}


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_codon_usage_accepts_tables_4_and_25(table: int) -> None:
    feats, flags = cu().codon_usage(TGA_GENES, TRNA_4, table)
    assert set(feats) == {g.locus_tag for g in TGA_GENES}
    assert flags == ()
    for f in feats.values():  # a code-4/25 CDS with internal TGA is scored
        assert f.gc3 is not None and f.enc is not None


def test_ac11_tga_tai_weight_is_positive_at_code4() -> None:
    """AC 11 / D10: tRNA-Trp(UCA) reads TGA at code 4; a TGA-only gene's tAI is w(TGA) > 0
    (0.84), while at 11 it has no sense codon and no tAI (#478 review B2: a TGG codon would
    score 1.0 at either table)."""
    g = _gene("W1", "ATG" + "TGA" * 12 + "TAA")
    feats, _ = cu().codon_usage([*TGA_GENES, g], {"TCA": 1, "CCA": 1}, 4)
    assert feats["W1"].tai is not None and feats["W1"].tai == pytest.approx(0.8409, abs=1e-3)
    f11, _ = cu().codon_usage([*TGA_GENES, g], {"TCA": 1, "CCA": 1}, 11)
    assert f11["W1"].tai is None


@pytest.fixture(scope="module")
def codon_work(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    from helpers.codon_excerpt import write_inputs

    from gippa.steps import load_step_modules, registry

    load_step_modules()
    tmp = tmp_path_factory.mktemp("codon16e")
    js, fa = write_inputs(tmp / "in")
    names = {"genes": "genes.parquet", "features": "features.parquet",
             "replicons": "replicons.parquet", "rrna16s": "rrna16s.fa",
             "trna_counts": "trna_counts.tsv", "xrefs": "xrefs.parquet",
             "annotation_meta": "annotation_meta.json"}  # fmt: skip
    outs = {k: tmp / "work" / v for k, v in names.items()}
    registry.run_step("parse_bakta", GippaConfig(seed=1), {"bakta_json": js, "genome": fa}, outs)
    return outs


def _codon_step(work: dict[str, Path], out: Path, cfg: GippaConfig) -> Path:
    from gippa.steps import registry

    ins = {"genes": work["genes"], "trna_counts": work["trna_counts"]}
    registry.run_step("codon_usage", cfg, ins, {"codon": out})
    return out


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_codon_usage_step_runs_at_4_and_25(
    table: int, codon_work: dict[str, Path], tmp_path: Path
) -> None:
    import pyarrow.parquet as pq

    a = _codon_step(codon_work, tmp_path / "t.parquet", cfg_at(table))
    b = _codon_step(codon_work, tmp_path / "d.parquet", GippaConfig(seed=1))
    ra, rb = pq.read_table(a).to_pylist(), pq.read_table(b).to_pylist()
    assert [r["locus_tag"] for r in ra] == [r["locus_tag"] for r in rb]
    assert any(r["enc"] is not None for r in ra)
    # the table reaches the scores (#478 review N2): ENC and tAI move with the Trp/Gly families
    for col in ("enc", "tai"):
        assert any(x[col] != y[col] for x, y in zip(ra, rb, strict=True)), col


def test_codon_usage_step_table11_is_byte_identical(
    codon_work: dict[str, Path], tmp_path: Path
) -> None:
    """Pin (passes today): ``translation_table: 11`` gives the default run's bytes."""
    a = _codon_step(codon_work, tmp_path / "a.parquet", cfg_at(11))
    b = _codon_step(codon_work, tmp_path / "b.parquet", GippaConfig(seed=1))
    assert a.read_bytes() == b.read_bytes()


@pytest.mark.parametrize("table", (1, 2, 12, 26))
def test_codon_usage_still_rejects_unsupported_tables(table: int) -> None:
    with pytest.raises(ValueError, match="4, 11 and 25"):
        cu().codon_usage(TGA_GENES, TRNA_4, table)


# ================================================================ 6.2 covariates (tests only)


@pytest.mark.parametrize(("table", "want"), [(4, "MKWA"), (25, "MKGA"), (11, "MK")])
def test_covariates_protein_sequence_follows_the_table(table: int, want: str) -> None:
    from gippa.genetable.covariates import covariate_table, protein_sequence

    nt = "ATGAAATGAGCCTAA"
    assert protein_sequence(None, nt, "pseudogene", table) == want
    genes = pd.DataFrame(
        {"locus_tag": ["P1"], "aa": [None], "nt": [nt], "feature_type": ["pseudogene"]}
    )
    assert covariate_table(genes, cfg_at(table))["length_aa"].tolist() == [len(want)]


# ======================================================================== 6.3 novel ORFs

# Example 11: ORF nt ATG AAA TGA GCC TAA at [3, 18) of a linear contig.
NOVEL_SEQ = "CCC" + "ATGAAATGAGCCTAA" + "CCCCCC"
NOVEL_REPS = {"c": Replicon("c", len(NOVEL_SEQ), "linear", "chromosome")}


def _novel_orf(**kw: Any) -> Any:
    from gippa.novel.candidates import Orf

    base: dict[str, Any] = dict(
        orf_id="NOVEL_c_4_p", cls="intergenic", contig="c", strand="+", frame=0, stop=15,
        starts=(3,), start_longest=3, end_5p_open=False, end_3p_open=False,
        peptides=("MKWA",), best_q=0.001, overlapping_gene=None, first_peptide=3,
    )  # fmt: skip
    base.update(kw)
    return Orf(**base)


@pytest.mark.parametrize(("table", "want"), [(4, "MKWA"), (25, "MKGA")])
def test_example11_novel_fasta_protein_reads_tga_as_sense(
    table: int, want: str, tmp_path: Path
) -> None:
    from gippa.io.novel_tables import write_novel_faa
    from gippa.novel.outputs import orf_table, protein_records

    genome = {"c": NOVEL_SEQ}
    t = orf_table([_novel_orf()], None, genome=genome, replicons=NOVEL_REPS)
    recs = protein_records(t, genome, NOVEL_REPS, translation_table=table)
    assert [(r.header, r.sequence) for r in recs] == [("NOVEL_c_4_p|longest", want)]
    path = tmp_path / "novel_orfs.faa"
    write_novel_faa(path, t, genome, NOVEL_REPS, translation_table=table)
    assert want in path.read_text().replace("\n", "")


def test_example11_novel_fasta_table11_unchanged() -> None:
    """Pin: at 11 the ORF ends at the TGA (stop [9, 12)): ``MK``, the default and ``=11`` agree."""
    from gippa.novel.outputs import orf_table, protein_records

    genome = {"c": NOVEL_SEQ}
    orf = _novel_orf(stop=9, peptides=("MK",))
    t = orf_table([orf], None, genome=genome, replicons=NOVEL_REPS)
    got = protein_records(t, genome, NOVEL_REPS, translation_table=11)
    assert [r.sequence for r in got] == ["MK"]
    assert got == protein_records(t, genome, NOVEL_REPS)


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_orf_without_position_walks_up_past_tga(table: int) -> None:
    """``_walk_up`` with the run's stops: TAA GCC TGA AAA AAA TAA, the peptide spans the TGA."""
    from gippa.novel.outputs import orf_table

    seq = "TAAGCCTGAAAAAAATAACCCCCC"
    reps = {"c": Replicon("c", len(seq), "linear", "chromosome")}
    pep = "A" + codon_table(table)["TGA"] + "KK"
    orf = _novel_orf(orf_id="NOVEL_c_4_p", stop=15, starts=(), start_longest=None,
                     first_peptide=None, peptides=(pep,))  # fmt: skip
    t = orf_table([orf], None, genome={"c": seq}, replicons=reps, translation_table=table)
    assert t["length_aa_longest"].tolist() == [4]


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_pseudogene_junction_ignores_tga_under_4_and_25(table: int) -> None:
    from helpers.novel_evidence_fixture import circ, pseudo_case
    from helpers.novel_oracle import cds

    import gippa.novel.evidence as ev

    seq = pseudo_case(None).genome["c1"]  # internal TGA [600, 603) is its only premature stop
    gene = cds("P", "c1", 540, 900, "+", "pseudogene")
    assert ev.pseudogene_junction(gene, seq, circ(), table=11) == 600
    assert ev.pseudogene_junction(gene, seq, circ(), table=table) is None
    minus = cds("P", "c1", 100, 460, "-", "pseudogene")
    assert ev.pseudogene_junction(minus, rc(seq), circ(), table=table) is None


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_apply_evidence_reads_the_configured_table(table: int) -> None:
    from helpers.novel_evidence_fixture import evidence_table, flags_of, pseudo_case

    cfg = cfg_at(table, annotation={"translation_table": table,
                                    "terminators": {"terminatornet_yes": 0.5}})  # fmt: skip
    f11 = flags_of(evidence_table(pseudo_case(None)).iloc[0])
    f = flags_of(evidence_table(pseudo_case(None), config=cfg).iloc[0])
    assert "junction_unknown" not in f11
    assert "junction_unknown" in f


def _calls(path: Path, name: str) -> list[ast.Call]:
    tree = ast.parse(path.read_text())
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            fn = node.func
            if (isinstance(fn, ast.Name) and fn.id == name) or (
                isinstance(fn, ast.Attribute) and fn.attr == name
            ):
                out.append(node)
    return out


@pytest.mark.parametrize(
    ("module", "name", "kw"),
    [
        ("steps/novel.py", "LocusClassifier", "table"),  # 6.3, lead (#476 review)
        ("steps/novel.py", "orf_table", "translation_table"),
        ("steps/novel.py", "write_novel_faa", "translation_table"),
        ("makedb/build.py", "LocusClassifier", "table"),  # 6.5
        ("makedb/build.py", "build_sixframe", "table"),
        ("novel/class_call.py", "makedb_id_map", "table"),  # lead (#478 impl report)
        ("novel/class_call.py", "nterm_entries", "table"),
    ],
)
def test_call_sites_pass_the_run_table(module: str, name: str, kw: str) -> None:
    """White-box: every call site passes the run's table by keyword (no silent table 11)."""
    calls = _calls(SRC / module, name)
    assert calls, f"{module}: no call to {name}"
    for c in calls:
        assert kw in {k.arg for k in c.keywords}, f"{module}:{c.lineno} {name}() without {kw}="


# ======================================================================== 6.4 variants


def vr() -> Any:
    import gippa.variants.vcf_records as m

    return m


def test_example12_codon_sets() -> None:
    f = vr()._codons_of
    assert f("W", 4) == ("TGA", "TGG")
    assert f("W", 11) == ("TGG",)
    assert f("G", 25) == ("GGA", "GGC", "GGG", "GGT", "TGA")
    assert f("W", 25) == ("TGG",)


@pytest.mark.parametrize("table", (4, 11, 25))
def test_codons_of_matches_the_table(table: int) -> None:
    ct = codon_table(table)
    for aa in sorted(set(ct.values()) - {"*"}):
        assert vr()._codons_of(aa, table) == tuple(sorted(c for c in ct if ct[c] == aa))


@settings(deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.text("ACGT", min_size=3, max_size=3), st.sampled_from("ACDEFGHKLMNPQRSTVWY"),
       st.sampled_from((4, 11, 25)))  # fmt: skip
def test_minimal_alt_codon_property(ref: str, alt: str, table: int) -> None:
    ct = codon_table(table)
    best = min((sum(a != b for a, b in zip(c, ref, strict=True)), c) for c in ct if ct[c] == alt)
    assert vr().minimal_alt_codon(ref, alt, table=table) == (best[1], best[0])


def test_u_has_no_zero_change_codon_under_4_and_25() -> None:
    assert vr().minimal_alt_codon("TGA", "U", table=11) == ("TGA", 0)  # M7-16 rule 11, unchanged
    assert vr().minimal_alt_codon("TAG", "O", table=4) == ("TAG", 0)  # Pyl in every table
    for table in _SENSE_TABLES:
        with pytest.raises(ValueError):
            vr().minimal_alt_codon("TGA", "U", table=table)


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_translate_locus_uses_the_table(table: int) -> None:
    from gippa.variants.align import translate_locus

    seq = "ATGTGAGCC"
    rep = Replicon("c", 9, "linear", "chromosome")
    sense = codon_table(table)["TGA"]
    assert translate_locus(seq, rep, 0, 9, "+", table=table) == "M" + sense + "A"
    assert translate_locus(rc(seq), rep, 0, 9, "-", table=table) == "M" + sense + "A"
    assert translate_locus(seq, rep, 0, 9, "+") == "M*A"


@pytest.mark.parametrize(("table", "alt"), [(4, "EK"), (25, "Q")])
def test_block_mass_equivalence_translates_ref_with_the_table(table: int, alt: str) -> None:
    """``TGA GCC`` is ``WA`` (257.116 Da ~ E+K 257.138) under 4, ``GA`` (128.059 ~ Q) under 25."""
    from gippa.variants.flags import flags_for

    site = {"alt_type": "block", "alt_aa": alt, "ref_aa": "*A", "ref_codon": "TGAGCC",
            "contig": "c", "strand": "+", "codon_start": 0}  # fmt: skip
    assert "mass_equivalent" not in flags_for(site, (), cfg_at(11), set())
    assert "mass_equivalent" in flags_for(site, (), cfg_at(table), set())


# --- run_variants on the M7-16c toy (tests/unit/test_variants_recoded.py world)

BODY, FLANK = "MSDKPLTEAGRWQLNEFK", "VPAGHTLDMRK"
N_TOY = 3 * (len(BODY) + 1 + len(FLANK) + 1)  # TGA codon [54, 57), terminal TAA [90, 93)


def _world(letter: str) -> tuple[dict[str, str], Any]:
    from helpers.readthrough_toys import back_translate, gene, model

    nt = back_translate(BODY + "*" + FLANK + "*", {18: "TGA", 30: "TAA"})
    reps = {"c1": Replicon("c1", N_TOY, "linear")}
    g1 = replace(gene("G1", "c1", 0, N_TOY, "+"), aa=BODY + letter + FLANK)
    return {"c1": nt}, model(reps, [g1])


def _variants(letter: str, alt: str, n_diff: int, table: int) -> tuple[Any, Any, list[Any]]:
    from gippa.steps.variants import run_variants

    genome, gm = _world(letter)
    pep = "WQLNEFK" + alt + "VPAGHTL"
    hits = pd.DataFrame([{
        "peptide": pep, "peptide_len": len(pep), "source": "denovo", "contig": "c1",
        "start": 33, "end": 78, "strand": "+", "n_diff": n_diff,
        "diff_types": "substitution" if n_diff else pd.NA, "placement_q": 0.0,
        "uniqueness": "unique",
    }])  # fmt: skip
    peptides = pd.DataFrame({"peptide": [pep], "search_db": "same"})
    cfg = cfg_at(table, genome_relation="same", sample_name="S1")
    return run_variants(hits, None, peptides, gm, genome, cfg)


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_peptide_with_the_table_residue_at_tga_is_a_match(table: int) -> None:
    """AC 12: what Bakta writes at TGA (W at 4, G at 25) is REF, never a ``stop_codon`` site."""
    sense = codon_table(table)["TGA"]
    pile, table_, records = _variants(sense, sense, 0, table)
    assert table_.empty and records == []
    row = pile[pile["codon_start"] == 54].to_dict("records")
    assert [(r["ref_aa"], r["n_ref"], r["call"]) for r in row] == [(sense, 1, "REF")]


@pytest.mark.parametrize(("table", "alt", "alt_codon"), [(4, "G", "GGA"), (25, "W", "TGG")])
def test_substitution_at_tga_is_ordinary_under_4_and_25(
    table: int, alt: str, alt_codon: str
) -> None:
    sense = codon_table(table)["TGA"]
    _pile, t, records = _variants(sense, alt, 1, table)
    (r,) = t.to_dict("records")
    assert (r["ref_codon"], r["ref_aa"], r["alt_aa"], r["call"]) == ("TGA", sense, alt, "ALT")
    assert set(r["flags"]) == {"snv"}  # no stop_codon
    assert r["alt_codon_vcf"] == alt_codon and int(r["nt_changes"]) == 1
    (rec,) = records
    assert (rec.pos0, rec.ref, rec.alts) == (54, "TGA", (alt_codon,))
    assert rec.info["AA_REF"] == sense


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_example12_u_over_tga_is_rare_skipped(table: int) -> None:
    """Example 12 / D6: a precomputed ``U`` under 4/25 gives no site, no ALT, one rare skip."""
    pile, t, records = _variants("U", "U", 1, table)
    assert t.empty and records == []
    assert pile.attrs["nonstandard"]["by_rule"]["rare"] == 1
    assert pile.attrs["nonstandard"]["by_residue"] == {"U": 1}


def test_u_over_tga_at_11_is_unchanged() -> None:
    """Pin (M7-16): ALT U over a recoded TGA, ``stop_codon;vcf_unrepresentable``, no record."""
    _pile, t, records = _variants("U", "U", 1, 11)
    (r,) = t.to_dict("records")
    assert (r["ref_aa"], r["alt_aa"]) == ("*", "U")
    assert set(r["flags"]) == {"stop_codon", "vcf_unrepresentable"} and records == []


@pytest.mark.parametrize("table", (4, 11, 25))
def test_variants_step_passes_the_table_to_annotated_stops(
    table: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    """6.4 (lead, #476 review): ``annotated_stops(gm, genome, recoded=True, table=)``."""
    import gippa.steps.variants as sv

    seen: list[Any] = []
    real = sv.annotated_stops

    def spy(*a: Any, **kw: Any) -> Any:
        seen.append(kw.get("table"))
        return real(*a, **kw)

    monkeypatch.setattr(sv, "annotated_stops", spy)
    sense = codon_table(table)["TGA"] if table != 11 else "U"
    _variants(sense, sense, 0 if table != 11 else 1, table)
    assert seen == [table]


# ======================================================================== 6.5 makedb

MAKEDB_11 = {
    "strain_db.faa": "5df59cec973dae6ec76b5f8c77a4f23505364a191d86e00ff61e62a09e72dcce",
    "db_id_map.tsv": "38393d60b4d235a1a4d0e63a49f1c71dbdcd05d3490e7fa1aa836dc6f0ba4a6b",
}
G1_NT = "ATG" + stretch(8, 1) + "TGA" + stretch(8, 5) + "TAA"  # CDS [300, 357), no Bakta aa


def _makedb_json(path: Path, recorded: int | None = None) -> tuple[Path, str]:
    """Linear 600 nt: an intergenic stretch(10) TGA stretch(10) at [99, 162) (two < 20 aa ORFs at
    11, one at 4/25); G1 with an internal TGA; upstream TAA [264, 267) + 7 codons + TGA [288,
    291) + 3 codons (no nterm_extension at 11: 3 codons; 11 codons at 4/25)."""
    c = Contig("c1", 600, "linear")
    c.put(99, stretch(10) + "TGA" + stretch(10, 3))
    c.put(264, "TAA" + stretch(7, 2) + "TGA" + stretch(3, 6))
    c.put(300, G1_NT)
    js = write_bakta_json(path, [c], [feature("cds", "c1", 300, 357, "+", locus="G1", product="p")])
    if recorded is not None:
        data = json.loads(js.read_text())
        data["genome"] = {"translation_table": recorded}
        js.write_text(json.dumps(data, indent=1) + "\n")
    return js, c.text


def _makedb(tmp: Path, **kw: Any) -> tuple[Path, str]:
    from gippa.makedb.build import MakedbOptions, build_strain_db

    tmp.mkdir(parents=True, exist_ok=True)
    js, seq = _makedb_json(tmp / "bakta.json")
    out = tmp / "DB"
    build_strain_db(
        MakedbOptions(annotation=js, genome=None, outdir=out, six_frame_orfs=True, **kw)
    )
    return out, seq


def _db(out: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    from gippa.io.protein_fasta import read_protein_fasta
    from gippa.io.tsv import read_tsv

    fa = {r.id: r.sequence for r in read_protein_fasta(out / "strain_db.faa")}
    header, rows = read_tsv(out / "db_id_map.tsv")
    return fa, [dict(zip(header, r, strict=True)) for r in rows]


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_makedb_translates_with_the_configured_table(table: int, tmp_path: Path) -> None:
    out, seq = _makedb(tmp_path, translation_table=table)
    fa, rows = _db(out)
    sense = codon_table(table)["TGA"]
    assert fa["G1"] == tr(G1_NT, table).rstrip("*") and sense in fa["G1"]  # _protein
    orfs = [r for r in rows if r["kind"] == "orf"]
    for r in orfs:  # every ORF record is the table's translation of its id_map interval
        s, e = int(r["start"]) - 1, int(r["end"])
        nt = seq[s:e] if r["strand"] == "+" else rc(seq[s:e])
        assert fa[r["db_id"]] == tr(nt, table) and "*" not in fa[r["db_id"]], r
    planted = tr(stretch(10) + "TGA" + stretch(10, 3), table)
    assert any(planted in fa[r["db_id"]] for r in orfs if r["class"] == "intergenic")
    (nt_ext,) = [r for r in orfs if r["class"] == "nterm_extension"]
    assert (nt_ext["locus_tag"], nt_ext["start"], nt_ext["end"]) == ("G1", "268", "354")
    manifest = json.loads((out / "makedb_manifest.json").read_text())
    assert manifest["translation_table"] == table


def test_makedb_table11_is_byte_identical(tmp_path: Path) -> None:
    """Pin: default and ``translation_table=11`` give today's bytes; the manifest gains the key."""
    for name, kw in (("default", {}), ("t11", {"translation_table": 11})):
        out, _ = _makedb(tmp_path / name, **kw)
        for f, digest in MAKEDB_11.items():
            assert hashlib.sha256((out / f).read_bytes()).hexdigest() == digest, (name, f)
        assert json.loads((out / "makedb_manifest.json").read_text())["translation_table"] == 11


def _patch(path: Path) -> Path:
    path.write_text("##gff-version 3\nc1\tGIPPA\tCDS\t301\t357\t.\t+\t0\tID=PX;locus_tag=G1\n")
    return path


@pytest.mark.parametrize("table", _SENSE_TABLES)
def test_makedb_patch_with_tga_is_accepted_under_4_and_25(table: int, tmp_path: Path) -> None:
    out, _ = _makedb(tmp_path, translation_table=table, patch=_patch(tmp_path / "p.gff3"))
    fa, _ = _db(out)
    assert fa["PX"] == "M" + tr(G1_NT[3:], table).rstrip("*")


def test_makedb_patch_with_tga_is_rejected_at_11(tmp_path: Path) -> None:
    """Pin (M3-23 review decision 7): an internal TGA is a stop at 11."""
    with pytest.raises(ValueError, match="internal stop"):
        _makedb(tmp_path, translation_table=11, patch=_patch(tmp_path / "p.gff3"))


def _cli(argv: list[str]) -> int:
    from gippa.cli import main

    with pytest.raises(SystemExit) as info:
        main(argv)
    return 0 if info.value.code is None else int(info.value.code)


@pytest.fixture
def no_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "nodb"))


@pytest.mark.usefixtures("no_home")
def test_cli_makedb_reads_the_table_from_set(tmp_path: Path) -> None:
    from gippa.exitcodes import ExitCode

    js, _ = _makedb_json(tmp_path / "bakta.json")
    out = tmp_path / "DB"
    argv = ["makedb", "--annotation", str(js), "--outdir", str(out), "--no-contaminants",
            "--six-frame-orfs", "--set", "annotation.translation_table=4"]  # fmt: skip
    assert _cli(argv) == ExitCode.SUCCESS
    assert json.loads((out / "makedb_manifest.json").read_text())["translation_table"] == 4
    assert _db(out)[0]["G1"] == tr(G1_NT, 4).rstrip("*")


@pytest.mark.usefixtures("no_home")
def test_cli_makedb_refuses_a_json_recorded_at_another_table(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """AC 10 / 3.4: a JSON with ``genome.translation_table: 11`` and a run at 4 exits
    ``INVALID_INPUT``."""
    from gippa.exitcodes import ExitCode

    js, _ = _makedb_json(tmp_path / "bakta.json", recorded=11)
    argv = ["makedb", "--annotation", str(js), "--outdir", str(tmp_path / "DB"),
            "--no-contaminants", "--set", "annotation.translation_table=4"]  # fmt: skip
    assert _cli(argv) == ExitCode.INVALID_INPUT
    err = capsys.readouterr().err
    assert "annotation JSON was made with translation table 11; this run uses 4" in err
    assert not (tmp_path / "DB" / "strain_db.faa").exists()
    ok = ["makedb", "--annotation", str(js), "--outdir", str(tmp_path / "DB11"),
          "--no-contaminants"]  # fmt: skip
    assert _cli(ok) == ExitCode.SUCCESS
