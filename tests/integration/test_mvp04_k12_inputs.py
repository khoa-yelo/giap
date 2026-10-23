"""MVP-04 giapsim inputs on the full K-12 genome (docs/specs/MVP-04-ecoli-evidence-run.md), part 1.

Slow: each test simulates the whole 4.6-Mb genome at ``K12_EVIDENCE_SEED`` from the committed
FASTA and recorded Bakta JSON (hash-checked, rule 1). Outputs go to pytest temporary directories
only (rule 2). The RegulonDB TUs are used when ``$GIAP_DATA_DIR/regulondb/14.5.0/TUSet.txt``
exists (read-only); the RegulonDB-specific tests skip without it.

Statistical claims checked here (single fixed seed, as the acceptance gates):

* acceptance 8 / rule 7: the tuned intercept puts the truth-detected share of eligible genes in
  [0.55, 0.70]; example 8: the M3-19 default b0 = -2.0 saturates (>= 0.90);
* Q4: 0.5 SNV/kb -> about 2,300 SNVs (Poisson SD ~48; band [2,100, 2,500]);
* step 3: AR(1) rho = 0.5 over TUs -> lag-1 correlation of adjacent non-island TU means near
  0.5 * 5.33 / (5.33 + ~0.5) ~ 0.46 (band [0.30, 0.60], SE ~0.02); genes of one TU share the TU
  value (correlation > 0.70, expected ~0.88); the four islands (one 20-TU high, three 8-TU
  lows, lead decision 20) shift by +-3 stationary SDs
  (+-4.8 natural-log units), checked at more than half that;
* Q5 / M7-12 decision 26: upper-tercile novel ORFs let at least 2 classes reach >= 100 ORFs with a
  true peptide.
"""

from __future__ import annotations

import importlib
import itertools
import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers import mvp04_data as M

pytestmark = pytest.mark.slow


def _sc() -> Any:
    return importlib.import_module("giapsim.scenarios")


def _simulate(out: Path, cfg: Any = None) -> list[Any]:
    sc = _sc()
    tu = M.TUSET if M.TUSET.is_file() else None
    with M.capture_warnings() as recs:
        importlib.import_module("giapsim.real_genome").simulate_on_genome(
            cfg if cfg is not None else sc.k12_evidence(), sc.K12_EVIDENCE_SEED, out,
            genome_fasta=M.K12_FASTA, bakta_json=M.K12_JSON, tu_set=tu,
        )  # fmt: skip
    return recs


@pytest.fixture(scope="module")
def k12(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, list[Any]]:
    out = tmp_path_factory.mktemp("mvp04_k12")
    return out, _simulate(out)


@pytest.fixture(scope="module")
def genes() -> list[M.Gene]:
    return M.protein_genes(M.k12_json())


@pytest.fixture(scope="module")
def k12_seq() -> str:
    return M.k12_text()[1]


def _share(out: Path, genes: list[M.Gene]) -> float:
    eligible = {g.locus for g in genes if g.eligible}
    return len(M.true_gene_loci(out) & eligible) / len(eligible)


def test_truth_detected_share_is_in_the_tuned_band(
    k12: tuple[Path, list[Any]], genes: list[M.Gene]
) -> None:
    """Acceptance 8 (rule 7): eligible genes (D8: truncated ones excluded) with >= 1 true
    peptide, at the pinned b0 and K12_EVIDENCE_SEED."""
    share = _share(k12[0], genes)
    assert 0.55 <= share <= 0.70, share
    prov = json.loads((k12[0] / "provenance.json").read_text())
    assert prov["b0"] == _sc().k12_evidence().proteomics.detect_logit[0]


def test_default_intercept_saturates(tmp_path: Path, genes: list[M.Gene]) -> None:
    """Example 8: b0 = -2.0 (the M3-19 default) gives a share ~0.98, which is why it is tuned."""
    cfg = _sc().k12_evidence()
    _, b1, b2 = cfg.proteomics.detect_logit
    cfg = replace(cfg, proteomics=replace(cfg.proteomics, detect_logit=(-2.0, b1, b2)))
    _simulate(tmp_path, cfg)
    assert _share(tmp_path, genes) >= 0.90


def test_inputs_untouched_and_hashes_recorded(k12: tuple[Path, list[Any]]) -> None:
    assert M.sha256(M.K12_FASTA) == M.SHA_FASTA_FILE
    assert M.sha256(M.K12_JSON) == M.SHA_JSON_FILE
    flat = {str(v) for v in M.values(json.loads((k12[0] / "provenance.json").read_text()))}
    assert {M.SHA_FASTA_FILE, M.SHA_SEQUENCE_UPPER, M.SHA_JSON_FILE} <= flat


def test_snv_count_matches_half_per_kb(k12: tuple[Path, list[Any]], k12_seq: str) -> None:
    sample = M.read_fasta(k12[0] / "truth" / "sample_genome.fna")[M.CONTIG]
    assert len(sample) == len(k12_seq)
    n = sum(a != b for a, b in zip(k12_seq, sample, strict=True))
    assert 2_100 <= n <= 2_500, n
    assert len(M.read_tsv(k12[0] / "truth" / "alt_sites.tsv")) >= 1_000


def test_one_high_and_three_low_islands_spread_and_not_overlapping(
    k12: tuple[Path, list[Any]], genes: list[M.Gene]
) -> None:
    """Lead decision 20: one 20-TU high island at +3 SD from TU 970 and three 8-TU lows at -3 SD
    from TUs 560, 1380 and 2200 (0-based ``first_tu``), disjoint and spread."""
    out = k12[0]
    rows = M.read_tsv(out / "truth" / "islands.tsv")
    got = sorted((int(r["first_tu"]), int(r["n_tus"]), float(r["effect_sd"])) for r in rows)
    assert got == [(560, 8, -3.0), (970, 20, 3.0), (1380, 8, -3.0), (2200, 8, -3.0)]
    sets = M.island_genes(out)
    by_id = {r["island_id"]: int(r["n_tus"]) for r in rows}
    assert all(len(v) >= by_id[iid] for iid, v in sets.items())
    vals = list(sets.values())
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            assert not vals[i] & vals[j]
    pos = {g.locus: (g.start + g.end) / 2 for g in genes}
    centres = sorted(float(np.mean([pos[g] for g in v])) for v in vals)
    assert centres[-1] - centres[0] >= M.K12_LENGTH / 2  # spread over the chromosome
    assert min(b - a for a, b in itertools.pairwise(centres)) >= M.K12_LENGTH / 20


def test_islands_shift_abundance(k12: tuple[Path, list[Any]], genes: list[M.Gene]) -> None:
    """Example 10's input: the down islands sit well below the background, the up ones above."""
    out = k12[0]
    ab = M.abundance_of(out)
    rows = {
        r["island_id"]: float(r["effect_sd"]) for r in M.read_tsv(out / "truth" / "islands.tsv")
    }
    sets = M.island_genes(out)
    inside = set().union(*sets.values())
    bg = float(np.mean([ab[g.locus] for g in genes if g.locus in ab and g.locus not in inside]))
    for sign in (3.0, -3.0):
        loci = [g for iid, gs in sets.items() if rows[iid] == sign for g in gs if g in ab]
        delta = float(np.mean([ab[g] for g in loci])) - bg
        assert sign * delta > 2.4, (sign, delta)


def test_ar1_over_tus_and_shared_tu_values(
    k12: tuple[Path, list[Any]], genes: list[M.Gene]
) -> None:
    out = k12[0]
    ab = M.abundance_of(out)
    tus = M.tu_of(out)
    inside = set().union(*M.island_genes(out).values())
    by_tu: dict[int, list[float]] = {}
    for g in genes:
        if g.locus in ab:
            by_tu.setdefault(M.tu_number(tus[g.locus]), []).append(ab[g.locus])
    island_tu = {M.tu_number(tus[g]) for g in inside}
    order = sorted(by_tu)
    x, y = [], []
    for a, b in itertools.pairwise(order):
        if b == a + 1 and a not in island_tu and b not in island_tu:
            x.append(np.mean(by_tu[a]))
            y.append(np.mean(by_tu[b]))
    r1 = float(np.corrcoef(x, y)[0, 1])
    assert 0.30 <= r1 <= 0.60, r1
    pairs = [(v[0], v[1]) for t, v in by_tu.items() if len(v) >= 2 and t not in island_tu]
    assert len(pairs) >= 300
    rw = float(np.corrcoef(np.array(pairs).T)[0, 1])
    assert rw > 0.70, rw


def test_novel_orfs_rule_counts_and_supply(
    k12: tuple[Path, list[Any]], genes: list[M.Gene], k12_seq: str
) -> None:
    out, recs = k12
    rows = M.read_tsv(out / "truth" / "novel_orfs.tsv")
    sample = M.read_fasta(out / "truth" / "sample_genome.fna")[M.CONTIG]
    bad = M.novel_rule_violations(rows, k12_seq, sample, genes)
    assert not bad, bad[:10]
    bad = M.nterm_violations(rows, k12_seq, genes)
    assert not bad, bad[:10]
    n = {c: sum(r["class"] == c for r in rows) for c in (*M.NOVEL_CLASSES, "nterm_extension")}
    assert n["cds_outframe"] == 300 and n["antisense"] == 300, n
    assert 0 < n["intergenic"] <= 300 and 0 < n["nterm_extension"] <= 300, n
    if min(n.values()) < 300:
        assert any(r.levelname == "WARNING" for r in recs)


def test_upper_tercile_lets_two_classes_reach_100_with_a_peptide(
    k12: tuple[Path, list[Any]],
) -> None:
    """Q5 / M7-12 decision 26: a class needs >= 100 truth ORFs with a true peptide."""
    out = k12[0]
    hit = M.true_gene_loci(out)
    rows = M.read_tsv(out / "truth" / "novel_orfs.tsv")
    with_pep = {
        c: sum(r["class"] == c and r["orf_id"] in hit for r in rows) for c in M.NOVEL_CLASSES
    }
    assert sum(v >= 100 for v in with_pep.values()) >= 2, with_pep


def test_an_nterm_peptide_crosses_the_annotated_start(
    k12: tuple[Path, list[Any]], genes: list[M.Gene]
) -> None:
    """Example 7: the extended proteins are expressed from the upstream start, so some correct
    peptide of an extended locus lies (partly) upstream of the annotated start."""
    out = k12[0]
    by_end = {(g.end, "+"): g for g in genes} | {(g.start, "-"): g for g in genes}
    ext = {}
    for r in M.read_tsv(out / "truth" / "novel_orfs.tsv"):
        if r["class"] == "nterm_extension":
            key = (int(r["end"]), "+") if r["strand"] == "+" else (int(r["start"]) - 1, "-")
            ext[by_end[key].locus] = by_end[key]
    assert ext
    n = 0
    for r in M.read_tsv(out / "truth" / "peptides.tsv"):
        g = ext.get(r["locus_tag"])
        if g is None or r["is_wrong"] != "false":
            continue
        s0, e0 = int(r["start"]) - 1, int(r["end"])
        n += s0 < g.start if g.strand == "+" else e0 > g.end
    assert n >= 1


def test_regulondb_tus(k12: tuple[Path, list[Any]], genes: list[M.Gene]) -> None:
    if not M.TUSET.is_file():
        pytest.skip(f"RegulonDB 14.5.0 TUSet.txt absent under {M.REGULONDB_DIR}")
    out = k12[0]
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["tu_source"] == "regulondb"
    assert M.sha256(M.TUSET) in {str(v) for v in M.values(prov)}
    tus = M.tu_of(out)
    strand = {g.locus: g.strand for g in genes}
    members: dict[str, set[str]] = {}
    for g, t in tus.items():
        members.setdefault(t, set()).add(g)
    assert all(len({strand[g] for g in m}) == 1 for m in members.values())
    assert sum(len(m) > 1 for m in members.values()) >= 300
    lac = {g.gene: g.locus for g in genes if g.gene in ("lacZ", "lacY", "lacA")}
    assert len({tus[x] for x in lac.values()}) == 1  # the lac operon is one TU


def test_same_seed_full_genome_is_byte_identical(
    k12: tuple[Path, list[Any]], tmp_path: Path
) -> None:
    """Acceptance 9: two giapsim calls with the same seed give byte-identical tables."""
    _simulate(tmp_path)
    for rel in (
        "peptides.tsv",
        *(f"truth/{p.name}" for p in sorted((k12[0] / "truth").glob("*.tsv"))),
    ):
        assert (k12[0] / rel).read_bytes() == (tmp_path / rel).read_bytes(), rel
