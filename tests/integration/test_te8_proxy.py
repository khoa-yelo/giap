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

* ``counts``: ``peptides`` (the UniProt digest, missed <= 1, 7-30 aa; either one per
  (protein, start) product or one per distinct sequence), ``k12_exact``, ``nissle_exact``,
  ``nissle_tolerant``, ``nissle_unmapped``.
* ``shares``: ``k12_exact`` (= counts / peptides), ``nissle_tolerant``, ``nissle_unmapped``,
  ``nissle_placement_fdr`` (pooled decoy-estimated FDR among accepted Nissle placements).
* ``passed`` = ``k12_exact >= 0.995 and nissle_placement_fdr <= 0.01``; ``exceptions`` <= 50.

The K-12 exact count is checked against an independent oracle (I = L substring search in the
six-frame translation with the circular wrap). On this excerpt the oracle share is 1364/1373 =
0.9934 (nine N-terminal peptides of non-ATG-start proteins, which UniProt writes with M), so the
excerpt is expected NOT to pass; the full proteome (oracle share 0.9967) is the weekly gate.
"""

from __future__ import annotations

import gzip
import importlib
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


def exact_oracle(fasta: Path, genome: Path) -> dict[int, int]:
    """{n_peptides: n_exact} for the three counting conventions (products, distinct, distinct
    up to I = L); exact = an I = L substring of a six-frame translation (circular wrap)."""
    import ahocorasick

    from giap.mapping.digest import digest

    units = [p for seq in read_fasta(fasta).values() for _, _, p in digest(seq)]
    auto = ahocorasick.Automaton()
    for p in {il(p) for p in units}:
        auto.add_word(p, p)
    auto.make_automaton()
    hit: set[str] = set()
    for seq in read_fasta(genome).values():
        for frame in six_frames(seq):
            hit.update(p for _, p in auto.iter(il(frame)))
    distinct = set(units)
    return {
        len(units): sum(il(p) in hit for p in units),
        len(distinct): sum(il(p) in hit for p in distinct),
        len({il(p) for p in units}): len(hit),
    }


@pytest.fixture(scope="module")
def result() -> Any:
    for f in (EXCERPT, K12_GENOME, NISSLE_GENOME, K12_JSON, NISSLE_JSON):
        assert f.is_file(), f
    return _te_proxies().te8_proxy(K12_GENOME, K12_JSON, NISSLE_GENOME, NISSLE_JSON, EXCERPT)


def test_te8_excerpt_k12_exact_count_equals_oracle(result: Any) -> None:
    c, sh = result.counts, result.shares
    oracle = exact_oracle(EXCERPT, K12_GENOME)
    assert c["peptides"] in oracle, (c["peptides"], oracle)
    assert c["k12_exact"] == oracle[c["peptides"]]
    assert sh["k12_exact"] == pytest.approx(c["k12_exact"] / c["peptides"])
    assert sh["k12_exact"] < K12_MIN  # 0.9934 on this excerpt (module docstring)


def test_te8_excerpt_nissle_shares_reported_and_bounded(result: Any) -> None:
    c, sh = result.counts, result.shares
    n = c["peptides"]
    oracle = exact_oracle(EXCERPT, NISSLE_GENOME)
    assert c["nissle_exact"] <= oracle[n]  # never an exact hit the genome lacks
    assert c["nissle_exact"] + c["nissle_tolerant"] + c["nissle_unmapped"] <= n
    assert c["nissle_unmapped"] <= n - oracle[n]  # a peptide with an exact hit is not unmapped
    assert c["nissle_unmapped"] > 0  # Nissle lacks some K-12 genes (oracle share 0.74)
    for key in ("nissle_tolerant", "nissle_unmapped"):
        assert sh[key] == pytest.approx(c[key] / n)
    assert 0.0 <= sh["nissle_placement_fdr"] <= 1.0


def test_te8_excerpt_pass_flag_follows_the_rule(result: Any) -> None:
    sh = result.shares
    rule = sh["k12_exact"] >= K12_MIN and sh["nissle_placement_fdr"] <= FDR_MAX
    assert result.passed is rule
    assert result.passed is False  # the excerpt's K-12 share is below 0.995
    assert len(result.exceptions) <= 50
