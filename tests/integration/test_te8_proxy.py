"""TE8 in-silico proxy on a recorded UniProt K-12 excerpt (issue M3-22; GIAP_SPEC.md §11.6 TE8).

§11.6 TE8: "UniProt K-12 tryptic peptides on the K-12 genome map >= 99.5 % exactly; on the
Nissle genome, the decoy-estimated placement FDR at defaults is <= 1 % and tolerant plus unmapped
shares are reported". Lead decision 3 (2026-10-20): the real UniProt comparison runs only in the
weekly job (tests/e2e/test_te8_uniprot_weekly.py); CI uses the 43-protein excerpt
``tests/data/recorded/uniprot/2026_03/UP000000625/k12_excerpt.fasta`` (PROVENANCE.md there) and
never the network.

Interface assumed (spec "Public interface"; field names pinned here by the test author)::

    validation.proxies.te_proxies.te8_proxy(k12_genome, k12_json, nissle_genome, nissle_json,
                                            uniprot_fasta) -> ProxyResult

* ``counts``: ``peptides`` (distinct sequences of the UniProt digest, missed <= 1, 7-30 aa, exact
  string, no I = L collapsing; lead decision 4), ``non_atg_nterm``, ``k12_exact``,
  ``nissle_exact`` (accepted exact placements; decision 7), ``nissle_tolerant``,
  ``nissle_unmapped``.
* ``non_atg_nterm`` (decision 6): distinct peptides that are N-terminal (start 0) in a UniProt
  protein and equal an N-terminal digest product of a K-12 gene (``k12_json``) whose first
  genome codon is not ATG and whose Bakta protein starts with M (the TE1 strict rule, decision
  5). They leave the K-12 denominator: ``shares["k12_exact"] = k12_exact / (peptides -
  non_atg_nterm)``; ``k12_exact`` counts the remaining peptides with an exact K-12 hit.
* ``shares``: ``k12_exact``, ``nissle_tolerant``, ``nissle_unmapped`` (both over ``peptides``),
  ``nissle_placement_fdr`` (pooled decoy-estimated FDR among accepted Nissle placements).
* ``passed`` = ``k12_exact >= 0.995 and nissle_placement_fdr <= 0.01``; ``exceptions`` <= 50.

Pinned on the excerpt (independent oracle: I = L substring search in the six-frame translation
with the circular wrap): 1,373 peptides, 4 ``non_atg_nterm``, 1,364 ``k12_exact``, share
1364/1369 = 0.99635 >= 0.995 (the 5 remaining misses are UniProt N-termini where Bakta chose a
different start). Full proteome (release 2026_03, weekly): 151,424 / 354 / 150,918, share 0.99899.
"""

from __future__ import annotations

import gzip
import importlib
import json
import sys
from pathlib import Path
from typing import Any

import pytest
from helpers.contaminant_oracle import il, six_frames
from helpers.mapping_steps import write_test_db

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "tests" / "data"
EXCERPT = DATA / "recorded" / "uniprot" / "2026_03" / "UP000000625" / "k12_excerpt.fasta"
BAKTA = DATA / "recorded" / "bakta" / "1.12.2-full"
K12_GENOME, NISSLE_GENOME = (
    DATA / "genomes" / "k12_mg1655.fna.gz",
    DATA / "genomes" / "nissle1917.fna.gz",
)
K12_JSON, NISSLE_JSON = (
    BAKTA / "k12" / "k12_mg1655.json.gz",
    BAKTA / "nissle" / "nissle1917.json.gz",
)
K12_MIN, FDR_MAX = 0.995, 0.01
EXCERPT_PEPTIDES, EXCERPT_NON_ATG_NTERM, EXCERPT_K12_EXACT = 1373, 4, 1364  # decisions 4, 6

pytestmark = pytest.mark.slow


def _te_proxies() -> Any:
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return importlib.import_module("validation.proxies.te_proxies")


@pytest.fixture(autouse=True)
def _test_db(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GIAP_DB_DIR", str(write_test_db(tmp_path_factory.mktemp("test_db"))))


def read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    out: dict[str, list[str]] = {}
    cur = ""
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                cur = line[1:].split()[0]
                out[cur] = []
            elif line:
                out[cur].append(line)
    return {k: "".join(v) for k, v in out.items()}


def _digest_products(fasta: Path) -> list[tuple[int, str]]:
    from giap.mapping.digest import digest

    return [(st, p) for seq in read_fasta(fasta).values() for st, _, p in digest(seq)]


def exact_hits(peptides: set[str], genome: Path) -> set[str]:
    """Peptides with an I = L substring hit in a six-frame translation (circular wrap)."""
    import ahocorasick

    auto = ahocorasick.Automaton()
    for p in {il(p) for p in peptides}:
        auto.add_word(p, p)
    auto.make_automaton()
    hit: set[str] = set()
    for seq in read_fasta(genome).values():
        for frame in six_frames(seq):
            hit.update(p for _, p in auto.iter(il(frame)))
    return {p for p in peptides if il(p) in hit}


def k12_non_atg_nterm(fasta: Path) -> set[str]:
    """Decision 6: UniProt N-terminal peptides equal to a strict non-ATG K-12 N-terminal product."""
    from giap.mapping.digest import digest

    with gzip.open(K12_JSON, "rt") as fh:
        feats = json.load(fh)["features"]
    k12 = {
        p
        for f in feats
        if f["type"] in ("cds", "sorf") and f.get("aa") and not f.get("pseudo")
        and f["nt"][:3] != "ATG" and f["aa"][0] == "M"
        for st, _, p in digest(f["aa"])
        if st == 0
    }  # fmt: skip
    return {p for st, p in _digest_products(fasta) if st == 0} & k12


@pytest.fixture(scope="module")
def result() -> Any:
    for f in (EXCERPT, K12_GENOME, NISSLE_GENOME, K12_JSON, NISSLE_JSON):
        assert f.is_file(), f
    return _te_proxies().te8_proxy(K12_GENOME, K12_JSON, NISSLE_GENOME, NISSLE_JSON, EXCERPT)


def test_excerpt_oracle_matches_pinned_numbers() -> None:
    peptides = {p for _, p in _digest_products(EXCERPT)}
    excluded = k12_non_atg_nterm(EXCERPT)
    assert len(peptides) == EXCERPT_PEPTIDES
    assert len(excluded) == EXCERPT_NON_ATG_NTERM
    assert len(exact_hits(peptides - excluded, K12_GENOME)) == EXCERPT_K12_EXACT
    assert not exact_hits(excluded, K12_GENOME)  # the excluded peptides never map exactly


def test_te8_excerpt_k12_exact_is_pinned(result: Any) -> None:
    c, sh = result.counts, result.shares
    assert c["peptides"] == EXCERPT_PEPTIDES
    assert c["non_atg_nterm"] == EXCERPT_NON_ATG_NTERM
    assert c["k12_exact"] == EXCERPT_K12_EXACT
    denom = EXCERPT_PEPTIDES - EXCERPT_NON_ATG_NTERM
    assert sh["k12_exact"] == pytest.approx(EXCERPT_K12_EXACT / denom)  # 0.99635
    assert sh["k12_exact"] >= K12_MIN


def test_te8_excerpt_nissle_shares_reported_and_bounded(result: Any) -> None:
    c, sh = result.counts, result.shares
    n = c["peptides"]
    peptides = {p for _, p in _digest_products(EXCERPT)}
    n_exact = len(exact_hits(peptides, NISSLE_GENOME))
    assert c["nissle_exact"] <= n_exact  # accepted exact placements only (decision 7)
    assert c["nissle_exact"] + c["nissle_tolerant"] + c["nissle_unmapped"] <= n
    assert c["nissle_unmapped"] <= n - n_exact  # a peptide with an exact hit is not unmapped
    assert c["nissle_unmapped"] > 0  # Nissle lacks some K-12 genes (oracle exact share 0.74)
    for key in ("nissle_tolerant", "nissle_unmapped"):
        assert sh[key] == pytest.approx(c[key] / n)
    assert 0.0 <= sh["nissle_placement_fdr"] <= 1.0


def test_te8_excerpt_pass_flag_follows_the_rule(result: Any) -> None:
    sh = result.shares
    rule = sh["k12_exact"] >= K12_MIN and sh["nissle_placement_fdr"] <= FDR_MAX
    assert result.passed is rule
    assert len(result.exceptions) <= 50
