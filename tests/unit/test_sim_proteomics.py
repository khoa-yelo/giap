"""giapsim proteins, peptides and database-search errors (issue M3-19; GIAP_SPEC.md §11.3).

Examples 1-8 of docs/specs/M3-19-sim-peptides-db-errors.md, the truth-agrees-with-the-genome
checks the M3 consumers rely on (M3-10 placement FDR, M3-14 locus classes incl. ADR-0013
``cterm_extension``, M3-21 calibration with real protein-derived wrong peptides), determinism and
backward compatibility. Statistical claims are in
tests/calibration/test_sim_proteomics_calibration.py, invariants over random seeds in
tests/property/test_sim_proteomics_props.py; the independent oracle is
tests/helpers/peptide_oracle.py.

Conventions these tests pin where the spec is silent:

* ``ProteomicsConfig`` / ``ProteinSimConfig`` and ``digest(protein, digest_len=..,
  missed_cleavages=..)`` live in ``giapsim.proteomics``; ``SimConfig.proteomics`` (default ``None``)
  switches the module on and ``SimResult.truth["peptides"]`` is ``truth/peptides.tsv``.
* Expressed proteins are the *true* proteins: non-pseudo CDS and sORFs translated from their true
  start (N-terminal extensions included) and the unannotated expressed ORFs of
  ``truth/novel_orfs.tsv`` (§11.3 "unannotated expressed ORFs"; truth ``locus_tag`` = ``orf_id``).
  Pseudogenes are not expressed. No readthrough is planted, so no correct peptide crosses an
  annotated stop (M3-14 never calls ``cterm_extension`` on them).
* Peptide detection is one logistic draw per tryptic product (not per run); every row of
  ``peptides.tsv`` is one identification of one run. The share of wrong identifications is counted
  over distinct non-decoy peptides; PEP and q use π = ``pi_wrong``.
* Wrong peptides (both sources) never occur, up to I/L, in any true or annotated protein of the
  sample genome: otherwise GIAP (``i_equals_l``) would place them correctly and the truth label
  would lie. ``other_strain`` wrong peptides are Trypsin/P products of the reference proteins
  (``reference_bakta_like.json``), i.e. real protein-derived (M3-index lead decision 11).
* ``n_loci_genome`` is the number of distinct genomic loci whose six-frame translation is the
  peptide, exact or I/L-equivalent (either is accepted); ``NA``/``0`` for decoy rows. The
  review decisions (GIAP's exact matcher with I = L, 0 allowed, wrong rows counted too,
  ``alt_start_nterm``) are pinned in tests/unit/test_sim_proteomics_review.py.
* ``abundance.tsv`` has ``locus_tag`` and ``log_abundance`` (natural log) or ``abundance``.
* ``wrong_source="other_strain"`` without ``SimConfig.strain`` raises ``ValueError``.
"""

from __future__ import annotations

import dataclasses
import hashlib
import itertools
import math
import re
import time
from pathlib import Path
from typing import Any

import pytest
from helpers.peptide_oracle import (
    Run,
    empty,
    fnum,
    il,
    maps_to,
    mixture_fdr,
    mixture_pep,
    read_fasta,
    rows,
    sixframe_loci,
    truthy,
    trypsin_p,
    trypsin_p_peptides,
)

import giapsim
from giap.config.model import InputConfig
from giap.io.peptides import aggregate_peptides, read_peptide_table, validate_peptides
from giap.mapping.loci import classify_locus
from giap.model.coords import Interval, Replicon
from giap.model.genes import gene_model_from_bakta_json

TRUTH_COLUMNS = [
    "peptide", "is_wrong", "origin", "locus_tag", "contig", "start", "end", "strand",
    "n_loci_genome",
]  # fmt: skip
ORIGINS = {"gene", "shuffled", "other_strain", "decoy"}
MIX = (3.0, 1.0, 0.0, 1.0)
ALL_DETECTED = (50.0, 0.0, 0.0)  # logit >= 50: every tryptic product is detected
FEATURES: dict[str, Any] = {
    "genome_length": 100_000,
    "n_plasmids": 1,
    "plasmid_length": (6_000, 9_000),
    "linear_contig_length": 20_000,
    "origin_crossing_gene": True,
    "n_rrna_operons": 1,
    "n_trna": 4,
    "n_pseudogenes": 3,
    "n_sorfs": 4,
    "n_nterm_extensions": 4,
    "n_novel_orfs": (("intergenic", 3), ("cds_outframe", 3), ("antisense", 3)),
}
#: SHA-256 of the M0 output for ``SimConfig()``, seed 1 (as in tests/unit/test_sim_features.py).
M0_SHA256 = {
    "genome.fasta": "698373a45c787caa64c01d7050dd1d2a7b31e03771e4e389056ca510c6df1971",
    "bakta_like.json": "d312420eb3d404cb1abbb6d07844cccbc83c0ce5e9c912d3a02aa17da6c49064",
    "truth/genes.tsv": "9c2c98152eed218d3855ccddf478bf49ea42bc69727fce466111f704e0387b0f",
    "truth/features.tsv": "14c5e4866fe3beb761ec7652455e4cc90e7496ac3370b6def7eb1c0fc03c6075",
}
NEW_FILES = ("peptides.tsv", "proteins.faa", "abundance.tsv", "truth/peptides.tsv")
#: UniProt accession (https://www.uniprot.org/help/accession_numbers), bare or in ``db|ACC|NAME``.
UNIPROT_ACC = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")


def proteomics(**kw: Any) -> Any:
    from giapsim.proteomics import ProteomicsConfig

    return ProteomicsConfig(**kw)


def strain(**kw: Any) -> Any:
    return giapsim.StrainConfig(**kw)


def simulate(out: Path, seed: int, prot: Any = None, **kw: Any) -> Run:
    cfg = giapsim.SimConfig(**kw) if prot is None else giapsim.SimConfig(proteomics=prot, **kw)
    res = giapsim.simulate(cfg, seed, out)
    return Run(out, {"res": res, "cfg": cfg, "seed": seed})


@pytest.fixture(scope="module")
def mk(tmp_path_factory: pytest.TempPathFactory) -> Any:
    cache: dict[str, Run] = {}

    def make(name: str, seed: int, prot: Any = None, **kw: Any) -> Run:
        if name not in cache:
            cache[name] = simulate(tmp_path_factory.mktemp(name), seed, prot, **kw)
        return cache[name]

    return make


@pytest.fixture(scope="module")
def complete(mk: Any) -> Run:
    """Every tryptic product detected, no wrong identifications (example 2)."""
    run: Run = mk("complete", 11, proteomics(pi_wrong=0.0, detect_logit=ALL_DETECTED), **FEATURES)
    return run


@pytest.fixture(scope="module")
def default(mk: Any) -> Run:
    run: Run = mk("default", 12, proteomics(), **FEATURES)
    return run


@pytest.fixture(scope="module")
def wrong5(mk: Any) -> Run:
    """About 10^4 distinct peptides with pi_wrong = 0.05 (examples 3-6)."""
    prot = proteomics(pi_wrong=0.05, detect_logit=ALL_DETECTED)
    run: Run = mk("wrong5", 13, prot, genome_length=250_000, n_sorfs=3, n_nterm_extensions=2)
    return run


@pytest.fixture(scope="module")
def other(mk: Any) -> Run:
    """Wrong peptides from the related strain's proteome (example 8)."""
    prot = proteomics(pi_wrong=0.05, wrong_source="other_strain")
    st = strain(snv_rate_per_kb=10.0, indel_rate_per_kb=0.5)
    run: Run = mk("other", 14, prot, **{**FEATURES, "strain": st})
    return run


def replicons(run: Run) -> dict[str, Replicon]:
    out = {}
    for s in run.bakta["sequences"]:
        kind = s.get("type") if s.get("type") in ("chromosome", "plasmid") else "unknown"
        out[s["id"]] = Replicon(s["id"], int(s["length"]), s["topology"], kind)
    return out


def correct_rows(run: Run) -> list[dict[str, str]]:
    return [r for r in run.truth if r["origin"] == "gene"]


def wrong_rows(run: Run) -> list[dict[str, str]]:
    return [r for r in run.truth if r["origin"] in ("shuffled", "other_strain")]


def protein_offset(run: Run, r: dict[str, str]) -> int:
    """aa offset of a correct peptide inside its true protein, from the truth coordinates."""
    p = run.proteins[r["locus_tag"]]
    n = len(run.genome[p.contig])
    s0, e0 = int(r["start"]) - 1, int(r["end"])
    d = s0 - p.start if p.strand == "+" else p.end - e0
    if d < 0:
        d += n
    if d >= 3 * len(p.aa) + 3:
        d -= n
    assert d % 3 == 0, (r, p)
    return d // 3


# ---------------------------------------------------------------- config and outputs


def test_config_defaults_match_spec() -> None:
    from giapsim.proteomics import ProteomicsConfig

    c = ProteomicsConfig()
    assert tuple(c.digest_len) == (7, 30)
    assert c.missed_cleavages == 1
    assert tuple(c.detect_logit) == (-2.0, 0.08, 0.9)
    assert c.pi_wrong == 0.01
    assert c.wrong_source == "shuffled"
    assert tuple(c.score_mixture) == MIX
    assert c.n_decoy_rows == "match_wrong"
    assert c.search_db == "same"
    assert c.n_runs == 3
    assert c.proteins is None
    assert giapsim.SimConfig().proteomics is None


def test_simconfig_with_proteomics_stays_hashable() -> None:
    a = giapsim.SimConfig(proteomics=proteomics(pi_wrong=0.05))
    b = giapsim.SimConfig(proteomics=proteomics(pi_wrong=0.05))
    assert a == b and hash(a) == hash(b)
    assert a != giapsim.SimConfig(proteomics=proteomics(pi_wrong=0.02))


def test_outputs_written(default: Run) -> None:
    for name in NEW_FILES:
        assert (default.out / name).is_file(), name
    res = default.meta["res"]
    assert Path(res.truth["peptides"]) == default.out / "truth" / "peptides.tsv"
    header = (default.out / "truth" / "peptides.tsv").read_text().splitlines()[0].split("\t")
    assert header[: len(TRUTH_COLUMNS)] == TRUTH_COLUMNS
    assert {r["origin"] for r in default.truth} <= ORIGINS
    assert {r["origin"] for r in default.truth} >= {"gene", "shuffled", "decoy"}


def test_every_identification_has_truth(default: Run) -> None:
    truth = default.truth_by_peptide
    peps = {r["peptide"] for r in default.peptides}
    assert peps <= set(truth), sorted(peps - set(truth))[:5]
    targets = {r["peptide"] for r in default.target_rows()}
    listed = {r["peptide"] for r in default.truth if r["origin"] != "decoy"}
    assert listed == targets  # truth lists exactly the identified target peptides
    decoys = {r["peptide"] for r in default.decoy_rows()}
    assert decoys <= {r["peptide"] for r in default.truth if r["origin"] == "decoy"}


# ---------------------------------------------------------------- example 1: digest


def test_example1_digest_trypsin_p() -> None:
    from giapsim.proteomics import digest

    protein = "MKAAKPRGGREEDDLLSTAMKPEPTIDEWYGR"
    # fragments MK|AAK|PR|GGR|EEDDLLSTAMK|PEPTIDEWYGR (Trypsin/P also cuts K|P);
    # 7-30 aa with <= 1 missed cleavage:
    expected = {"EEDDLLSTAMK", "PEPTIDEWYGR", "GGREEDDLLSTAMK", "EEDDLLSTAMKPEPTIDEWYGR"}
    assert set(digest(protein, digest_len=(7, 30), missed_cleavages=1)) == expected
    assert "AAKPRGGR" not in set(digest(protein, digest_len=(7, 30), missed_cleavages=1))
    assert set(digest(protein, digest_len=(7, 30), missed_cleavages=0)) == {
        "EEDDLLSTAMK",
        "PEPTIDEWYGR",
    }


def test_example1_digest_c_terminal_peptide_without_k_r() -> None:
    from giapsim.proteomics import digest

    got = set(digest("MAAAAAAKGGGSSSTTTW", digest_len=(7, 30), missed_cleavages=1))
    assert got == {"MAAAAAAK", "GGGSSSTTTW", "MAAAAAAKGGGSSSTTTW"}


def test_example1_simulated_peptides_are_tryptic(default: Run) -> None:
    for r in default.target_rows():
        p = r["peptide"]
        assert 7 <= len(p) <= 30, p
        assert sum(c in "KR" for c in p[:-1]) <= 1, p  # Trypsin/P: every inner K/R is missed


# ---------------------------------------------------------------- example 2: truth = genome


def test_example2_no_wrong_peptides(complete: Run) -> None:
    assert complete.truth and not wrong_rows(complete)
    assert all(not truthy(r["is_wrong"]) for r in complete.truth if r["origin"] != "decoy")
    assert not complete.decoy_rows()  # match_wrong: no wrong rows, no decoy rows


def test_example2_every_peptide_maps_to_its_truth_locus(complete: Run) -> None:
    bad = []
    for r in correct_rows(complete):
        p = complete.proteins.get(r["locus_tag"])
        assert p is not None, f"{r['locus_tag']} is not an expressed true protein"
        ok = r["contig"] == p.contig and r["strand"] == p.strand
        ok = ok and maps_to(
            r["peptide"], complete.genome, r["contig"], int(r["start"]), int(r["end"]), r["strand"]
        )
        if ok:
            a = protein_offset(complete, r)
            ok = complete.proteins[r["locus_tag"]].aa[a : a + len(r["peptide"])] == r["peptide"]
        if not ok:
            bad.append(r)
    assert not bad, bad[:5]


def test_example2_all_tryptic_products_of_true_proteins_detected(complete: Run) -> None:
    expected: dict[str, set[str]] = {
        loc: trypsin_p_peptides(p.aa) for loc, p in complete.proteins.items()
    }
    got = {r["peptide"] for r in correct_rows(complete)}
    union = set().union(*expected.values())
    assert got == union, (sorted(union - got)[:5], sorted(got - union)[:5])
    for r in correct_rows(complete):
        assert r["peptide"] in expected[r["locus_tag"]], r


def test_nterm_extension_peptides_present(complete: Run) -> None:
    """Peptides from the upstream part of N-terminal extensions (true start, §11.3)."""
    ext = {
        loc: (len(p.aa) - len(complete.json_aa[loc]))
        for loc, p in complete.proteins.items()
        if loc in complete.json_aa and len(p.aa) > len(complete.json_aa[loc])
    }
    assert len(ext) == FEATURES["n_nterm_extensions"]
    upstream = [
        r for r in correct_rows(complete)
        if r["locus_tag"] in ext and protein_offset(complete, r) < ext[r["locus_tag"]]
    ]  # fmt: skip
    expected = {
        complete.proteins[loc].aa[a:b]
        for loc, k in ext.items()
        for a, b in trypsin_p(complete.proteins[loc].aa)
        if a < k
    }
    assert expected, "fixture has no tryptic product inside an N-terminal extension"
    assert {r["peptide"] for r in upstream} == expected


def test_no_peptides_from_pseudogenes(complete: Run) -> None:
    assert complete.pseudo_loci
    assert not {r["locus_tag"] for r in correct_rows(complete)} & complete.pseudo_loci


def test_novel_orf_peptides_present(complete: Run) -> None:
    novel = {r["orf_id"] for r in complete.novel}
    assert len(novel) == 9
    have = {r["locus_tag"] for r in correct_rows(complete)} & novel
    want = {o for o in novel if trypsin_p(complete.proteins[o].aa)}
    assert want and have == want


def test_locus_classes_of_correct_peptides(complete: Run) -> None:
    """M3-14 (incl. ADR-0013) classifies every planted peptide at its truth locus as planted."""
    gm = gene_model_from_bakta_json(complete.out / "bakta_like.json", replicons(complete))
    seen: dict[str, int] = {}
    bad = []
    for r in correct_rows(complete):
        p = complete.proteins[r["locus_tag"]]
        n = len(complete.genome[r["contig"]])
        s0, e0 = int(r["start"]) - 1, int(r["end"])
        if s0 >= n:
            s0, e0 = s0 - n, e0 - n
        cls, tag = classify_locus(Interval(r["contig"], s0, e0, r["strand"]), gm, complete.genome)
        if p.kind.startswith("novel:"):
            want, want_tag = p.kind.split(":")[1], None
        elif p.kind == "sorf":
            want, want_tag = "sorf", p.locus
        else:
            k = len(p.aa) - len(complete.json_aa[p.locus])
            want = "nterm_extension" if protein_offset(complete, r) < k else "cds_inframe"
            want_tag = p.locus
        seen[cls] = seen.get(cls, 0) + 1
        if (cls, tag) != (want, want_tag):
            bad.append((r["peptide"], r["locus_tag"], cls, tag, want))
    assert not bad, bad[:5]
    assert "cterm_extension" not in seen
    assert {"cds_inframe", "sorf", "nterm_extension", "intergenic", "antisense"} <= set(seen)


def test_n_loci_genome(complete: Run) -> None:
    topo = complete.topology
    for r in correct_rows(complete)[::7]:
        n = int(r["n_loci_genome"])
        if r["peptide"].startswith("M"):
            continue  # an initiator M may come from GTG/TTG: 0 loci allowed (review decision 3)
        assert n >= 1, r
        exact = len(sixframe_loci(r["peptide"], complete.genome, topo, il_equal=False))
        il_eq = len(sixframe_loci(r["peptide"], complete.genome, topo, il_equal=True))
        assert n in (exact, il_eq), (r, exact, il_eq)


def test_wrong_and_decoy_truth_rows_have_no_coordinates(default: Run) -> None:
    for r in default.truth:
        if r["origin"] == "gene":
            assert not truthy(r["is_wrong"])
            assert all(not empty(r[c]) for c in ("locus_tag", "contig", "start", "end", "strand"))
            continue
        if r["origin"] != "decoy":
            assert truthy(r["is_wrong"]), r
        for c in ("contig", "start", "end", "strand"):
            assert empty(r[c]), r
        if r["origin"] == "decoy":
            assert empty(r["n_loci_genome"]) or int(r["n_loci_genome"]) == 0, r


# ---------------------------------------------------------------- example 3: wrong share


def test_example3_wrong_share(wrong5: Run) -> None:
    labels = {r["peptide"]: truthy(r["is_wrong"]) for r in wrong5.truth if r["origin"] != "decoy"}
    n = len(labels)
    assert n >= 9_000, n
    share = sum(labels.values()) / n
    assert abs(share - 0.05) <= 0.007, share
    assert {r["origin"] for r in wrong_rows(wrong5)} == {"shuffled"}


def test_shuffled_wrong_peptides_do_not_occur_in_the_genome_proteins(wrong5: Run) -> None:
    hay = wrong5.proteome_and_annotation_il
    bad = [r["peptide"] for r in wrong_rows(wrong5) if il(r["peptide"]) in hay]
    assert not bad, bad[:5]
    for r in wrong_rows(wrong5):
        p = r["peptide"]
        assert 7 <= len(p) <= 30 and sum(c in "KR" for c in p[:-1]) <= 1, p


# ---------------------------------------------------------------- examples 4-5: PEP and q


def test_example4_pep_formula_value() -> None:
    assert mixture_pep(1.5, 0.05, MIX) == pytest.approx(0.05, abs=1e-12)


def test_example4_pep_of_every_target_row(wrong5: Run) -> None:
    rs = [r for r in wrong5.target_rows()]
    near = min(rs, key=lambda r: abs(float(r["score"]) - 1.5))
    assert abs(float(near["score"]) - 1.5) < 0.05
    assert float(near["pep"]) == pytest.approx(0.05, abs=0.01)
    for r in rs:
        s, pep = float(r["score"]), float(r["pep"])
        assert pep == pytest.approx(mixture_pep(s, 0.05, MIX), rel=1e-3, abs=1e-9), r


def test_example5_q_values(wrong5: Run) -> None:
    rs = sorted(wrong5.target_rows(), key=lambda r: -float(r["score"]))
    qs = [float(r["q_value"]) for r in rs]
    assert all(0.0 <= q <= 1.0 for q in qs)
    assert all(a <= b + 1e-12 for a, b in itertools.pairwise(qs))  # non-decreasing in -score
    for r in rs:
        s, q = float(r["score"]), float(r["q_value"])
        assert q == pytest.approx(mixture_fdr(s, 0.05, MIX), rel=1e-3, abs=1e-9), r
    lowest_correct = [r for r in rs if not wrong5.is_wrong(r["peptide"])][-1]
    assert float(lowest_correct["q_value"]) <= 1.0


def test_pep_and_q_use_configured_mixture(mk: Any) -> None:
    mix = (4.0, 0.8, -0.5, 1.2)
    prot = proteomics(pi_wrong=0.1, score_mixture=mix, detect_logit=ALL_DETECTED)
    run: Run = mk("mix", 15, prot, genome_length=60_000)
    for r in run.target_rows():
        s = float(r["score"])
        assert float(r["pep"]) == pytest.approx(mixture_pep(s, 0.1, mix), rel=1e-3, abs=1e-9)
        assert float(r["q_value"]) == pytest.approx(mixture_fdr(s, 0.1, mix), rel=1e-3, abs=1e-9)


# ---------------------------------------------------------------- example 6: decoys


def test_example6_decoy_rows(wrong5: Run) -> None:
    decoys = wrong5.decoy_rows()
    assert decoys
    for r in decoys:
        assert truthy(r["is_decoy"])
        ids = r["protein_ids"].split(";")
        assert all(i.startswith("rev_") for i in ids), r
        for c in ("score", "pep", "q_value"):
            assert empty(r.get(c)), r
        prots = [
            wrong5.proteins[i[4:]].aa if i[4:] in wrong5.proteins else wrong5.json_aa.get(i[4:], "")
            for i in ids
        ]
        assert any(r["peptide"][::-1] in trypsin_p_peptides(p) for p in prots), r
    for r in wrong5.target_rows():
        assert not any(i.startswith("rev_") for i in r["protein_ids"].split(";")), r


def test_decoy_rows_match_wrong(wrong5: Run) -> None:
    n_wrong = sum(wrong5.is_wrong(r["peptide"]) for r in wrong5.target_rows())
    assert n_wrong > 0 and len(wrong5.decoy_rows()) == n_wrong


def test_decoy_rows_fixed_count(mk: Any) -> None:
    run: Run = mk("ndecoy", 16, proteomics(pi_wrong=0.02, n_decoy_rows=37), genome_length=60_000)
    assert len(run.decoy_rows()) == 37


# ---------------------------------------------------------------- example 7: M1-06 format


@pytest.mark.parametrize("name", ["complete", "default", "wrong5", "other"])
def test_example7_peptide_table_passes_m106(name: str, request: pytest.FixtureRequest) -> None:
    run: Run = request.getfixturevalue(name)
    df_raw, issues = read_peptide_table(run.out / "peptides.tsv")
    df, issues2 = validate_peptides(df_raw, InputConfig())
    errors = [i for i in [*issues, *issues2] if i.severity == "error"]
    assert not errors, errors
    assert not [i for i in [*issues, *issues2] if i.code in ("short_peptide", "long_peptide")]
    assert len(df) == len(run.peptides)
    agg, _ = aggregate_peptides(df)
    assert set(agg["peptide"]) == {r["peptide"] for r in run.peptides}


def test_per_run_columns(default: Run) -> None:
    rs = default.peptides
    assert {r["run"] for r in rs} and len({r["run"] for r in rs}) == 3
    by_run: dict[str, list[str]] = {}
    for r in rs:
        assert r["source"] == "database" and r["search_db"] == "same", r
        assert r["charge"] in ("2", "3"), r
        if not truthy(r["is_decoy"]):
            x = fnum(r["intensity"])
            assert x is not None and x > 0 and math.isfinite(x), r
        assert not empty(r["spectrum_id"]) and not empty(r["protein_ids"]), r
        by_run.setdefault(r["run"], []).append(r["spectrum_id"])
    for run, ids in by_run.items():
        assert len(ids) == len(set(ids)), run
    for r in default.target_rows():
        if default.is_wrong(r["peptide"]):
            continue
        tags = {
            t["locus_tag"]
            for t in default.truth_by_peptide[r["peptide"]]
            if r["peptide"] in default.json_aa.get(t["locus_tag"], "")
        }
        if tags:  # a peptide of an annotated protein names that protein
            assert tags & set(r["protein_ids"].split(";")), r


def test_search_db_other_and_n_runs(mk: Any) -> None:
    run: Run = mk("other_db", 17, proteomics(search_db="other", n_runs=5), genome_length=60_000)
    assert {r["search_db"] for r in run.peptides} == {"other"}
    assert len({r["run"] for r in run.peptides}) == 5


# ---------------------------------------------------------------- example 8: other strain


def test_example8_other_strain_wrong_peptides(other: Run) -> None:
    wrong = wrong_rows(other)
    assert wrong and {r["origin"] for r in wrong} == {"other_strain"}
    ref = rows_json_aa(other.out / "reference_bakta_like.json")
    ref_products = set().union(*(trypsin_p_peptides(a) for a in ref))
    hay = other.proteome_and_annotation_il
    for r in wrong:
        p = r["peptide"]
        assert p in ref_products, f"{p}: not a tryptic product of a reference protein"
        assert il(p) not in hay, f"{p}: occurs in the sample genome's proteins"


def test_example8_correct_peptides_still_map_to_the_sample_genome(other: Run) -> None:
    for r in correct_rows(other)[::5]:
        assert maps_to(
            r["peptide"], other.genome, r["contig"], int(r["start"]), int(r["end"]), r["strand"]
        ), r


def test_other_strain_needs_a_strain(tmp_path: Path) -> None:
    cfg = giapsim.SimConfig(proteomics=proteomics(wrong_source="other_strain"))
    with pytest.raises(ValueError):
        giapsim.simulate(cfg, 1, tmp_path)


def rows_json_aa(path: Path) -> list[str]:
    import json

    data = json.loads(path.read_text())
    return [f["aa"] for f in data["features"] if f["type"] in ("cds", "sorf") and f.get("aa")]


# ---------------------------------------------------------------- proteins.faa


def test_proteins_faa_holds_bakta_proteins(default: Run) -> None:
    faa = read_fasta(default.out / "proteins.faa")
    for loc, aa in default.json_aa.items():
        if loc in default.pseudo_loci:
            continue
        assert faa.get(loc) == aa, loc
    assert set(faa) <= set(default.json_aa)  # no cross-strain copies unless configured


def test_proteins_faa_cross_strain_copies(mk: Any) -> None:
    from Bio.Align import PairwiseAligner

    from giapsim.proteomics import ProteinSimConfig

    run: Run = mk("copies", 18, proteomics(proteins=ProteinSimConfig()), **FEATURES)
    faa = read_fasta(run.out / "proteins.faa")
    copies = {k: v for k, v in faa.items() if k not in run.json_aa}
    assert copies
    for loc, aa in run.json_aa.items():
        if loc not in run.pseudo_loci:
            assert faa.get(loc) == aa, loc
    aligner = PairwiseAligner(mode="global", match_score=1, mismatch_score=0)
    aligner.open_gap_score, aligner.extend_gap_score = -1, -0.5
    sources = list(run.json_aa.values())
    for name, aa in list(copies.items())[:15]:
        assert UNIPROT_ACC.fullmatch(name) or (
            "|" in name and UNIPROT_ACC.fullmatch(name.split("|")[1])
        ), name
        kmers = {aa[i : i + 5] for i in range(len(aa) - 4)}
        src = max(sources, key=lambda s: sum(s[i : i + 5] in kmers for i in range(len(s) - 4)))
        aln = aligner.align(src, aa)[0]
        ident = sum(
            a == b for a, b in zip(aln[0], aln[1], strict=True) if a != "-" and b != "-"
        ) / max(len(src), len(aa))
        assert 0.79 <= ident <= 0.995, (name, ident)


# ---------------------------------------------------------------- determinism and compatibility


def tree(d: Path) -> dict[str, str]:
    return {
        str(p.relative_to(d)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(d.rglob("*"))
        if p.is_file() and p.name != "config.yaml"
    }


def test_same_seed_byte_identical(tmp_path: Path) -> None:
    prot = proteomics(pi_wrong=0.05)
    a = tree(simulate(tmp_path / "a", 21, prot, **FEATURES).out)
    b = tree(simulate(tmp_path / "b", 21, prot, **FEATURES).out)
    assert a == b
    assert all(n in a for n in NEW_FILES)


def test_different_seed_differs(tmp_path: Path) -> None:
    prot = proteomics(pi_wrong=0.05)
    a = simulate(tmp_path / "a", 22, prot, genome_length=60_000).out
    b = simulate(tmp_path / "b", 23, prot, genome_length=60_000).out
    for name in NEW_FILES:
        assert (a / name).read_bytes() != (b / name).read_bytes(), name


def test_proteomics_options_change_only_new_outputs(tmp_path: Path) -> None:
    """Peptides come from their own random stream: the genome and its truth stay the same."""
    base = tree(simulate(tmp_path / "off", 24, None, **FEATURES).out)
    for i, prot in enumerate([proteomics(), proteomics(pi_wrong=0.2, n_runs=1)]):
        on = tree(simulate(tmp_path / f"on{i}", 24, prot, **FEATURES).out)
        assert {k: v for k, v in on.items() if k in base} == base
        assert set(on) - set(base) >= set(NEW_FILES)


def test_proteomics_off_writes_nothing_new(tmp_path: Path) -> None:
    simulate(tmp_path / "a", 25, None, **FEATURES)
    for name in (*NEW_FILES, "reference_peptides.tsv"):
        assert not (tmp_path / "a" / name).exists(), name
    a = tree(tmp_path / "a")
    giapsim.simulate(giapsim.SimConfig(proteomics=None, **FEATURES), 25, tmp_path / "b")
    assert tree(tmp_path / "b") == a


def test_m0_output_unchanged_with_proteomics_none(tmp_path: Path) -> None:
    giapsim.simulate(giapsim.SimConfig(proteomics=None), 1, tmp_path)
    for name, digest in M0_SHA256.items():
        assert hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() == digest, name


def test_sim100k_has_no_proteomics() -> None:
    assert giapsim.SIM100K.proteomics is None
    assert dataclasses.replace(giapsim.SIM100K, proteomics=proteomics()).proteomics is not None


# ---------------------------------------------------------------- speed


def test_100kb_with_peptides_under_1s(tmp_path: Path) -> None:
    cfg = giapsim.SimConfig(proteomics=proteomics(), **FEATURES)
    best = float("inf")
    for i in range(3):
        t0 = time.perf_counter()
        giapsim.simulate(cfg, 31, tmp_path / str(i))
        best = min(best, time.perf_counter() - t0)
    # loose 10x guard as in test_sim_features.py; the 1 s budget belongs in benchmarks/
    assert best < 10.0, best


@pytest.mark.slow
def test_5mb_with_peptides_under_30s(tmp_path: Path) -> None:
    cfg = giapsim.SimConfig(
        proteomics=proteomics(), **{**FEATURES, "genome_length": 5_000_000, "n_sorfs": 40}
    )
    t0 = time.perf_counter()
    giapsim.simulate(cfg, 32, tmp_path)
    assert time.perf_counter() - t0 < 300.0  # loose 10x guard of the 30 s budget
    assert len(rows(tmp_path / "truth" / "peptides.tsv")) > 10_000
