"""gippasim on the real K-12 genome, part 1 (docs/specs/MVP-04-ecoli-evidence-run.md).

Covers everything that does not depend on ``decoys = "scored_null"`` (that mode, examples 2-4 and
rule 6, comes in part 2): the input hash check (rule 1, example 1, D1), reading the committed
genome and recorded Bakta JSON, the SNV-only sample strain with swapped roles (step 2, Q4), TUs
from RegulonDB or the gap rule (step 3, Q3, example 9), the TU-level AR(1) abundance with planted
islands, centred to ``LOG_ABUNDANCE[0]``, the novel-ORF rule (step 4, Q5, Q7), the outputs (step 8),
TE7 on the inputs (acceptance 9) and the empty-peptide refusal (D9).

All runs use a test-time 200-kb excerpt of the committed FASTA and recorded JSON (spec rule 3),
written to pytest temporary directories only. The RegulonDB tables are hand-made,
RegulonDB-*shaped* fixtures over real Bakta loci (no RegulonDB bytes, R11).

Assumed interface (spec "Interface", names the implementer may refine within gippasim):

* ``gippasim.real_genome.simulate_on_genome(cfg, seed, outdir, *, genome_fasta, bakta_json,
  tu_set, check_hashes=True)``; ``check_hashes=False`` lets an excerpt through rule 1.
* ``gippasim.real_genome.verify_inputs(genome_fasta, bakta_json)`` raises a ``ValueError`` whose
  text contains ``input_hash_mismatch:<file name>``; ``python -m gippasim.real_genome --genome G
  --annotation J --seed S --outdir D`` exits non-zero with that text on stderr.
* ``GeneProductSet.txt`` (gene coordinates) is read from the directory of ``tu_set``.
* ``truth/genes.tsv`` ``tu_id`` values are the layout's ``TU_<n>`` ids, so ``truth/islands.tsv``
  ``first_tu_id`` / ``last_tu_id`` name them.
* novel ORFs appear in ``abundance.tsv`` under their ``orf_id``.
"""

from __future__ import annotations

import hashlib
import importlib
import itertools
import json
import math
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import yaml
from helpers import mvp04_data as M

from gippasim.expression import IslandSpec

SEED = 20261023
N_NOVEL = 20
#: Six excerpt islands of 8 TUs (the 200-kb excerpt has ~110 gap-rule TUs).
ISLANDS = ((2, 3.0), (18, -3.0), (34, 3.0), (50, -3.0), (66, 3.0), (82, -3.0))


def rg() -> ModuleType:
    return importlib.import_module("gippasim.real_genome")


def scenario() -> Any:
    return importlib.import_module("gippasim.scenarios").k12_evidence()


def excerpt_cfg(
    n_novel: int = N_NOVEL, islands: tuple[tuple[int, float], ...] = ISLANDS, **prot: Any
) -> Any:
    base = scenario()
    expr = replace(
        base.expression, islands=tuple(IslandSpec(M.CONTIG, f, 8, e) for f, e in islands)
    )
    cfg = replace(base, expression=expr, n_novel_per_class=n_novel)
    if prot:
        cfg = replace(cfg, proteomics=replace(base.proteomics, **prot))
    return cfg


def run(
    ex: M.Excerpt, outdir: Path, cfg: Any = None, seed: int = SEED, tu_set: Path | None = None
) -> list[Any]:
    """One excerpt run; returns the WARNING records it logged."""
    with M.capture_warnings() as recs:
        rg().simulate_on_genome(
            cfg if cfg is not None else excerpt_cfg(), seed, outdir,
            genome_fasta=ex.fasta, bakta_json=ex.bakta_json, tu_set=tu_set, check_hashes=False,
        )  # fmt: skip
    return recs


@pytest.fixture(scope="module")
def ex(tmp_path_factory: pytest.TempPathFactory) -> M.Excerpt:
    return M.make_excerpt(tmp_path_factory.mktemp("k12_excerpt"))


@pytest.fixture(scope="module")
def base(ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, list[Any], str]:
    """The reference excerpt run (gap rule: no TUSet); (outdir, warnings, input FASTA sha256)."""
    out = tmp_path_factory.mktemp("mvp04_base")
    before = M.sha256(ex.fasta)
    recs = run(ex, out)
    return out, recs, before


def sample_seq(out: Path) -> str:
    seqs = M.read_fasta(out / "truth" / "sample_genome.fna")
    assert list(seqs) == [M.CONTIG]
    return seqs[M.CONTIG]


# ---------------------------------------------------------------- scenario (config only)


def test_k12_evidence_scenario_pins_the_spec_choices() -> None:
    sc = importlib.import_module("gippasim.scenarios")
    cfg = sc.k12_evidence()
    assert isinstance(cfg, rg().RealGenomeConfig)
    assert isinstance(sc.K12_EVIDENCE_SEED, int)
    # Q4: SNVs only, 0.5 per kb
    assert cfg.strain.snv_rate_per_kb == 0.5
    assert cfg.strain.indel_rate_per_kb == 0.0
    # step 3: TU-level AR(1) at rho 0.5, RegulonDB TUs first
    assert cfg.expression.ar1_rho == 0.5
    assert cfg.tu_source == "regulondb"
    # step 4 / Q5: 300 per class (M7-12 NOVEL_PER_CLASS_GATE), 30-150 aa
    assert cfg.n_novel_per_class == 300
    assert tuple(cfg.novel_len_aa) == (30, 150)
    # step 5: the related sample searched, pi_wrong 0.01, 3 runs
    p = cfg.proteomics
    assert p.search_db == "other"
    assert p.pi_wrong == 0.01
    assert p.n_runs == 3


def test_tuned_intercept_changes_b0_only_and_leaves_the_default_alone() -> None:
    """Rule 7: only b0 is tuned (the default -2.0 saturates, example 8); the other detect_logit
    terms and every existing default stay as they are (no pin changes)."""
    from gippasim.proteomics import ProteomicsConfig

    b0, b1, b2 = scenario().proteomics.detect_logit
    assert (b1, b2) == (0.08, 0.9)
    assert math.isfinite(b0) and b0 < -2.0
    assert ProteomicsConfig().detect_logit == (-2.0, 0.08, 0.9)


def test_b0_and_share_are_documented() -> None:
    """Rule 7: b0 (and the measured share) are recorded in the scenario docstring."""
    sc = importlib.import_module("gippasim.scenarios")
    b0 = sc.k12_evidence().proteomics.detect_logit[0]
    doc = (sc.k12_evidence.__doc__ or "") + (sc.__doc__ or "")
    assert f"{b0:g}" in doc or str(b0) in doc
    assert not math.isnan(b0)


def test_recorded_inputs_have_the_spec_counts() -> None:
    genes = M.protein_genes(M.k12_json())
    assert len(M.k12_text()[1]) == M.K12_LENGTH
    assert len(genes) == 4373
    assert sum(g.strand == "+" for g in genes) == 2140
    assert sum(g.truncated for g in genes) == 6
    assert all(g.end <= M.K12_LENGTH and g.start < g.end for g in genes)  # none crosses the origin


# ---------------------------------------------------------------- rule 1: input hashes


def test_committed_inputs_match_the_spec_hashes() -> None:
    assert M.sha256(M.K12_FASTA) == M.SHA_FASTA_FILE
    assert hashlib.sha256(M.k12_text()[1].encode()).hexdigest() == M.SHA_SEQUENCE_UPPER
    assert M.sha256(M.K12_JSON) == M.SHA_JSON_FILE
    rg().verify_inputs(M.K12_FASTA, M.K12_JSON)  # the committed pair passes


def _k12_copy(dst: Path, *, flip_at: int | None = None, header: str | None = None) -> Path:
    head, seq = M.k12_text()
    if flip_at is not None:
        seq = seq[:flip_at] + ("A" if seq[flip_at] != "A" else "C") + seq[flip_at + 1 :]
    dst.mkdir(parents=True, exist_ok=True)
    return M.write_fasta_gz(dst / "k12_mg1655.fna.gz", header or head, seq)


def test_one_changed_base_is_an_input_hash_mismatch(tmp_path: Path) -> None:
    """Example 1: one byte changed in the sequence; the step fails and writes nothing."""
    fasta = _k12_copy(tmp_path / "in", flip_at=1_000_000)
    with pytest.raises(ValueError, match=r"input_hash_mismatch:k12_mg1655\.fna\.gz"):
        rg().verify_inputs(fasta, M.K12_JSON)
    out = tmp_path / "out"
    with pytest.raises(ValueError, match=r"input_hash_mismatch:k12_mg1655\.fna\.gz"):
        rg().simulate_on_genome(
            scenario(), SEED, out, genome_fasta=fasta, bakta_json=M.K12_JSON, tu_set=None
        )
    assert not (out / "peptides.tsv").exists()
    assert not (out / "truth").exists()


def test_lost_topology_tag_fails_the_hash_check_first(tmp_path: Path) -> None:
    """D1: the file hash covers the header, so a FASTA without ``[topology=circular]`` (same
    sequence) is refused before gippa could read the replicon as linear."""
    head, _ = M.k12_text()
    assert head.endswith("[topology=circular]")
    fasta = _k12_copy(tmp_path, header=head.replace(" [topology=circular]", ""))
    with pytest.raises(ValueError, match=r"input_hash_mismatch:k12_mg1655\.fna\.gz"):
        rg().verify_inputs(fasta, M.K12_JSON)


def test_changed_annotation_json_is_an_input_hash_mismatch(tmp_path: Path) -> None:
    data = M.k12_json()
    data["features"] = data["features"][:-1]
    js = M.write_json_gz(tmp_path / "k12_mg1655.json.gz", data)
    with pytest.raises(ValueError, match=r"input_hash_mismatch:k12_mg1655\.json\.gz"):
        rg().verify_inputs(M.K12_FASTA, js)


def test_input_step_cli_exits_non_zero_on_a_mismatch(tmp_path: Path) -> None:
    fasta = _k12_copy(tmp_path / "in", flip_at=17)
    out = tmp_path / "out"
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(M.REPO / "tests" / "sim"), str(M.REPO / "src"), env.get("PYTHONPATH", "")]
    )
    proc = subprocess.run(
        [sys.executable, "-m", "gippasim.real_genome", "--genome", str(fasta), "--annotation",
         str(M.K12_JSON), "--seed", str(SEED), "--outdir", str(out)],
        capture_output=True, text=True, env=env, timeout=600, check=False,
    )  # fmt: skip
    assert proc.returncode != 0
    assert "input_hash_mismatch:k12_mg1655.fna.gz" in proc.stderr + proc.stdout
    assert not (out / "peptides.tsv").exists()


# ---------------------------------------------------------------- outputs (step 8)


def test_outputs_are_written(base: tuple[Path, list[Any], str]) -> None:
    out, _, _ = base
    expected = ["peptides.tsv", "abundance.tsv", "config.yaml", "provenance.json"]
    expected += [
        f"truth/{n}.tsv" for n in ("peptides", "genes", "islands", "alt_sites", "novel_orfs")
    ]
    expected.append("truth/sample_genome.fna")
    for rel in expected:
        assert (out / rel).is_file(), rel
    assert not (out / "sample_genome.fna").exists()  # never next to the gippa inputs


def test_config_yaml_holds_only_seed_and_sample_name(base: tuple[Path, list[Any], str]) -> None:
    """No ``annotation.precomputed`` layers: the real tools run (step 8)."""
    cfg = yaml.safe_load((base[0] / "config.yaml").read_text())
    assert set(cfg) == {"seed", "sample_name"}
    assert cfg["seed"] == SEED


def test_provenance_records_hashes_seed_b0_and_tu_source(
    ex: M.Excerpt, base: tuple[Path, list[Any], str]
) -> None:
    prov = json.loads((base[0] / "provenance.json").read_text())
    assert prov["seed"] == SEED
    assert prov["b0"] == scenario().proteomics.detect_logit[0]
    assert prov["tu_source"] == "gap_rule"
    flat = {str(v) for v in M.values(prov)}
    assert M.sha256(ex.fasta) in flat
    assert M.sha256(ex.bakta_json) in flat


def test_input_genome_is_never_modified(ex: M.Excerpt, base: tuple[Path, list[Any], str]) -> None:
    assert M.sha256(ex.fasta) == base[2]
    assert M.read_fasta(ex.fasta)[M.CONTIG] == ex.seq


def test_peptides_table_has_target_rows(base: tuple[Path, list[Any], str]) -> None:
    rows = M.read_tsv(base[0] / "peptides.tsv")
    assert any(r["is_decoy"] == "false" for r in rows)
    assert {r["run"] for r in rows} == {"run1", "run2", "run3"}


# ---------------------------------------------------------------- step 2: sample strain


def test_sample_differs_from_k12_by_snvs_only(
    ex: M.Excerpt, base: tuple[Path, list[Any], str]
) -> None:
    """Q4: SNVs only, so coordinates are identical; rate 0.5/kb (Poisson mean ~100 on 200 kb,
    the band is about +-4 SD)."""
    s = sample_seq(base[0])
    assert len(s) == ex.length
    diff = [i for i, (a, b) in enumerate(zip(ex.seq, s, strict=True)) if a != b]
    assert 60 <= len(diff) <= 145, len(diff)
    assert set(s) <= set("ACGT")


def test_alt_sites_have_k12_as_ref_and_the_sample_as_alt(
    ex: M.Excerpt, base: tuple[Path, list[Any], str]
) -> None:
    """Swapped roles (step 2, example 6): ``ref`` = the K-12 residue, ``alt`` = the sample's."""
    out = base[0]
    s = sample_seq(out)
    rows = M.read_tsv(out / "truth" / "alt_sites.tsv")
    assert rows, "0.5 SNV/kb leaves non-synonymous sites in 200 kb"
    genes = {g.locus: g for g in ex.genes}
    for r in rows:
        g = genes[r["gene"]]
        c1, strand = int(r["codon_start"]), r["strand"]
        assert strand == g.strand
        assert g.start < c1 and c1 + 2 <= g.end
        assert (c1 - 1 - g.start) % 3 == 0 if strand == "+" else (g.end - (c1 + 2)) % 3 == 0
        ref = M.translate(M.segment(ex.seq, c1, c1 + 2, strand))
        alt = M.translate(M.segment(s, c1, c1 + 2, strand))
        assert (r["ref_aa"], r["alt_aa"]) == (ref, alt), r
        assert ref != alt


def test_true_peptides_lie_on_the_sample_at_k12_coordinates(
    base: tuple[Path, list[Any], str],
) -> None:
    """With SNVs only, sample and K-12 coordinates agree: every correct peptide is the sample's
    translation of its stated span (the N-terminal residue may be the initiator M)."""
    out = base[0]
    s = sample_seq(out)
    n = 0
    for r in M.read_tsv(out / "truth" / "peptides.tsv"):
        if r["is_wrong"] != "false" or r["origin"] == "decoy":
            continue
        aa = M.translate(M.segment(s, int(r["start"]), int(r["end"]), r["strand"]))
        pep = r["peptide"]
        assert aa[1:] == pep[1:] and (aa[0] == pep[0] or pep[0] == "M"), r
        n += 1
    assert n > 0


# ---------------------------------------------------------------- step 3: TUs and abundance


def test_genes_tsv_has_one_row_per_recorded_protein_gene(
    ex: M.Excerpt, base: tuple[Path, list[Any], str]
) -> None:
    out = base[0]
    assert M.header_of(out / "truth" / "genes.tsv")[:2] == ["locus_tag", "tu_id"]
    tus = M.tu_of(out)
    assert set(tus) == {g.locus for g in ex.genes}


def test_gap_rule_when_tuset_is_absent(ex: M.Excerpt, base: tuple[Path, list[Any], str]) -> None:
    """Example 9: no TUSet -> adjacent same-strand genes with gap <= 40 bp share a TU, all
    other neighbours do not; ``tu_source = gap_rule`` with a WARNING."""
    out, recs, _ = base
    tus = M.tu_of(out)
    genes = ex.genes
    for a, b in itertools.pairwise(genes):
        joined = a.strand == b.strand and b.start - a.end <= 40
        assert (tus[a.locus] == tus[b.locus]) == joined, (a.locus, b.locus, b.start - a.end)
    assert any("gap_rule" in r.getMessage() for r in recs)


def test_regulondb_tus_match_by_strand_and_stop(ex: M.Excerpt, tmp_path: Path) -> None:
    """Step 3: a Bakta gene matches a RegulonDB gene on the same strand with the same stop
    coordinate (the right end on +, the left end on -); TUs are the connected components of
    same-strand genes sharing a RegulonDB TU; unmatched genes and genes in no TU are singletons."""
    genes = ex.genes
    ok = [g for g in genes if g.eligible and g.kind == "cds"]
    i = next(k for k in range(len(ok) - 2) if all(g.strand == "+" for g in ok[k : k + 3]))
    a, b, c = ok[i : i + 3]
    j = next(k for k in range(len(ok) - 1) if ok[k].strand == ok[k + 1].strand == "-")
    d, e = ok[j : j + 2]
    rest = [g for g in ok if g not in (a, b, c, d, e)]
    f, h = rest[len(rest) // 3], rest[2 * len(rest) // 3]

    def rec(name: str, g: M.Gene, left_shift: int = 0, right_shift: int = 0,
            strand: str | None = None) -> tuple[str, int, int, str]:  # fmt: skip
        return name, g.start + 1 + left_shift, g.end + right_shift, strand or g.strand

    rdb_genes = [
        rec("fixA", a, left_shift=-30),  # + strand, other start, same stop: still matches
        rec("fixB", b),
        rec("fixC", c),
        rec("fixD", d, right_shift=30),  # - strand, other start (right end), same stop
        rec("fixE", e),
        rec("fixF", f, strand="-" if f.strand == "+" else "+"),  # opposite strand: no match
        rec("fixH", h),
        ("fixGhost", 1_000_001, 1_000_300, "+"),  # RegulonDB-only
    ]
    tuset = M.write_regulondb_fixture(
        tmp_path / "rdb",
        rdb_genes,
        [("FIXTU1", ["fixA", "fixB"]), ("FIXTU2", ["fixB", "fixC"]), ("FIXTU3", ["fixD", "fixE"]),
         ("FIXTU4", ["fixF", "fixH"]), ("FIXTU5", ["fixGhost", "fixH"])],
    )  # fmt: skip
    out = tmp_path / "out"
    run(ex, out, tu_set=tuset)
    tus = M.tu_of(out)
    members: dict[str, set[str]] = {}
    for g, t in tus.items():
        members.setdefault(t, set()).add(g)
    assert members[tus[a.locus]] == {a.locus, b.locus, c.locus}
    assert members[tus[d.locus]] == {d.locus, e.locus}
    assert members[tus[f.locus]] == {f.locus}
    assert members[tus[h.locus]] == {h.locus}
    multi = {frozenset(m) for m in members.values() if len(m) > 1}
    assert multi == {frozenset({a.locus, b.locus, c.locus}), frozenset({d.locus, e.locus})}
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["tu_source"] == "regulondb"
    assert M.sha256(tuset) in {str(v) for v in M.values(prov)}


def test_abundance_is_centred_on_log_abundance_mean(
    ex: M.Excerpt, base: tuple[Path, list[Any], str]
) -> None:
    from gippasim.proteomics import LOG_ABUNDANCE

    ab = M.abundance_of(base[0])
    vals = [ab[g.locus] for g in ex.genes if g.locus in ab]
    assert len(vals) >= 0.95 * len(ex.genes)
    assert abs(float(np.mean(vals)) - LOG_ABUNDANCE[0]) < 0.1


def test_island_truth_mirrors_the_planted_islands(base: tuple[Path, list[Any], str]) -> None:
    out = base[0]
    rows = M.read_tsv(out / "truth" / "islands.tsv")
    assert len(rows) == 6
    ecfg = excerpt_cfg().expression
    stationary = ecfg.tu_sd / math.sqrt(1 - ecfg.ar1_rho**2)
    for r, (first, eff) in zip(rows, ISLANDS, strict=True):
        assert int(r["first_tu"]) == first and int(r["n_tus"]) == 8
        assert float(r["effect_sd"]) == eff
        assert float(r["effect_log2"]) == pytest.approx(eff * stationary)
    sets = list(M.island_genes(out).values())
    assert all(sets)
    for k in range(len(sets)):
        for m in range(k + 1, len(sets)):
            assert not sets[k] & sets[m]  # not overlapping


def test_islands_shift_tu_abundance_exactly(
    ex: M.Excerpt, base: tuple[Path, list[Any], str], tmp_path: Path
) -> None:
    """The islands are deterministic TU shifts: the same seed with zero-effect islands differs,
    after the common centring constant, by exactly ``effect_log2 * ln 2`` on island genes and 0
    elsewhere (expression_for_simulation is natural log)."""
    out0 = tmp_path / "flat"
    run(ex, out0, cfg=excerpt_cfg(islands=tuple((f, 0.0) for f, _ in ISLANDS)))
    a, z = M.abundance_of(base[0]), M.abundance_of(out0)
    loci = [g.locus for g in ex.genes if g.locus in a]
    shift = {}
    for r in M.read_tsv(base[0] / "truth" / "islands.tsv"):
        shift[r["island_id"]] = float(r["effect_log2"]) * math.log(2.0)
    isl = M.island_genes(base[0])
    expect = {g: 0.0 for g in loci}
    for iid, gs in isl.items():
        for g in gs:
            expect[g] = shift[iid]
    resid = np.array([a[g] - z[g] - expect[g] for g in loci])
    assert np.ptp(resid) < 1e-4, (resid.min(), resid.max())


# ---------------------------------------------------------------- step 4: novel ORFs


def test_novel_orfs_follow_the_rule(ex: M.Excerpt, base: tuple[Path, list[Any], str]) -> None:
    out = base[0]
    rows = M.read_tsv(out / "truth" / "novel_orfs.tsv")
    assert M.header_of(out / "truth" / "novel_orfs.tsv")[:7] == [
        "orf_id", "class", "contig", "start", "end", "strand", "aa"]  # fmt: skip
    bad = M.novel_rule_violations(rows, ex.seq, sample_seq(out), ex.genes)
    assert not bad, bad[:10]


def test_each_class_gets_n_novel_per_class(base: tuple[Path, list[Any], str]) -> None:
    rows = M.read_tsv(base[0] / "truth" / "novel_orfs.tsv")
    n = {c: sum(r["class"] == c for r in rows) for c in (*M.NOVEL_CLASSES, "nterm_extension")}
    assert n["cds_outframe"] == N_NOVEL and n["antisense"] == N_NOVEL
    assert 0 < n["intergenic"] <= N_NOVEL
    assert 0 < n["nterm_extension"] <= N_NOVEL
    keys = [
        (r["contig"], r["strand"], r["end"] if r["strand"] == "+" else r["start"]) for r in rows
    ]
    assert len(keys) == len(set(keys)), "an ORF was drawn twice"


def test_nterm_extensions_on_both_strands(ex: M.Excerpt, base: tuple[Path, list[Any], str]) -> None:
    """Example 7 / D2: the extension starts at an in-frame start >= 5 codons upstream (higher
    coordinates on -) with no stop between; the recorded annotation stays unchanged."""
    rows = M.read_tsv(base[0] / "truth" / "novel_orfs.tsv")
    nterm = [r for r in rows if r["class"] == "nterm_extension"]
    assert {r["strand"] for r in nterm} == {"+", "-"}
    bad = M.nterm_violations(nterm, ex.seq, ex.genes)
    assert not bad, bad[:10]
    assert M.read_fasta(ex.fasta)[M.CONTIG] == ex.seq


def test_novel_orf_abundance_comes_from_the_upper_tercile(
    ex: M.Excerpt, base: tuple[Path, list[Any], str]
) -> None:
    """Q5: novel ORFs (not N-terminal extensions, which keep their gene) draw their abundance
    from the upper tercile of the gene distribution."""
    out = base[0]
    ab = M.abundance_of(out)
    gene_vals = np.array([ab[g.locus] for g in ex.genes if g.locus in ab])
    cut = float(np.quantile(gene_vals, 2 / 3))
    novel = [
        r for r in M.read_tsv(out / "truth" / "novel_orfs.tsv") if r["class"] in M.NOVEL_CLASSES
    ]
    assert novel
    for r in novel:
        assert ab[r["orf_id"]] >= cut - 1e-3, (r["orf_id"], ab[r["orf_id"]], cut)
    assert max(ab[r["orf_id"]] for r in novel) <= gene_vals.max() + 1e-6


def test_short_supply_takes_what_exists_and_origin_crossers_are_excluded(tmp_path: Path) -> None:
    """Q7 / D1 and the short-supply rule: the excerpt is rotated so that a feature-free intergenic
    ORF crosses the origin. Asking for far more ORFs than exist realises fewer (with a WARNING),
    and no realised ORF, that one included, crosses the origin."""
    plain = M.make_excerpt(tmp_path / "plain")
    first, last, stop = M.free_intergenic_orf(plain.seq, plain.data)
    k = last + 3  # every start of that ORF lies before the new origin, its stop after it
    ex = M.make_excerpt(tmp_path / "rot", rotate=k)
    assert ex.seq[(first - k) % ex.length :][:3] in M.STARTS
    new_stop_end = (stop - k) % ex.length + 3  # 1-based last base of the stop codon
    assert M.translate(ex.seq[new_stop_end - 3 : new_stop_end]) == "*"
    out = tmp_path / "out"
    recs = run(ex, out, cfg=excerpt_cfg(n_novel=100_000))
    rows = M.read_tsv(out / "truth" / "novel_orfs.tsv")
    n = {c: sum(r["class"] == c for r in rows) for c in M.NOVEL_CLASSES}
    assert all(0 < v < 100_000 for v in n.values()), n
    assert any(r.levelname == "WARNING" for r in recs)
    for r in rows:
        assert 1 <= int(r["start"]) < int(r["end"]) <= ex.length, r
        assert not (r["strand"] == "+" and int(r["end"]) == new_stop_end), r
    bad = M.novel_rule_violations(rows, ex.seq, sample_seq(out), ex.genes)
    assert not bad, bad[:10]


# ---------------------------------------------------------------- TE7 and degenerate cases


def _tables(out: Path) -> dict[str, bytes]:
    files = [out / "peptides.tsv", out / "abundance.tsv", *sorted((out / "truth").glob("*.tsv"))]
    return {str(p.relative_to(out)): p.read_bytes() for p in files}


def test_same_seed_gives_byte_identical_inputs(
    ex: M.Excerpt, base: tuple[Path, list[Any], str], tmp_path: Path
) -> None:
    """Acceptance 9 (TE7 on the inputs)."""
    again = tmp_path / "again"
    run(ex, again)
    first, second = _tables(base[0]), _tables(again)
    assert set(first) == set(second)
    assert {"truth/peptides.tsv", "truth/genes.tsv", "truth/islands.tsv", "truth/alt_sites.tsv",
            "truth/novel_orfs.tsv"} <= set(first)  # fmt: skip
    for name in first:
        assert first[name] == second[name], name
    assert (base[0] / "config.yaml").read_bytes() == (again / "config.yaml").read_bytes()


def test_another_seed_gives_other_inputs(
    ex: M.Excerpt, base: tuple[Path, list[Any], str], tmp_path: Path
) -> None:
    other = tmp_path / "other"
    run(ex, other, seed=SEED + 1)
    assert (other / "peptides.tsv").read_bytes() != (base[0] / "peptides.tsv").read_bytes()
    assert sample_seq(other) != sample_seq(base[0])


def test_no_target_rows_is_refused(ex: M.Excerpt, tmp_path: Path) -> None:
    """D9: the input step never writes a ``peptides.tsv`` without target rows."""
    out = tmp_path / "empty"
    cfg = excerpt_cfg(detect_logit=(-1000.0, 0.08, 0.9), pi_wrong=0.0)
    with pytest.raises(ValueError):
        run(ex, out, cfg=cfg)
    assert not (out / "peptides.tsv").exists()
