"""TE1 in-silico proxy (issue M3-22 examples 4-5; GIAP_SPEC.md §11.6 TE1, §13 gate M3).

§11.6: "In-silico proxy: in-silico tryptic peptides of Nissle Bakta proteins map >= 99.9 %
exactly in frame, the exceptions being N-terminal peptides at non-ATG starts".

Interface assumed (spec "Public interface"; field names pinned here by the test author)::

    validation.proxies.te_proxies.te1_proxy(genome: Path, bakta_json: Path) -> ProxyResult

``ProxyResult`` "holds counts, shares, the pass/fail flag, and up to 50 example exceptions":

* ``counts: Mapping[str, int]`` with ``peptides`` (distinct digest sequences, exact string, no
  I = L collapsing; lead decision 4), ``non_atg_nterm`` (distinct N-terminal products of genes
  whose first genome codon is not ATG **and** whose Bakta protein starts with M; decision 5;
  excluded from the share),
  ``exact_inframe`` (``route = exact`` with ``locus_class`` in {cds_inframe, sorf} at the
  peptide's own gene) and ``misses`` (the rest); ``peptides == non_atg_nterm + exact_inframe +
  misses``.
* ``shares: Mapping[str, float]`` with ``exact_inframe = exact_inframe / (peptides -
  non_atg_nterm)``.
* ``passed: bool`` (``shares["exact_inframe"] >= 0.999``).
* ``exceptions: Sequence[Mapping[str, Any]]``, at most 50, each with ``locus_tag`` and
  ``peptide``.

Pinned on the Nissle recording (oracle: this file's digest of the recorded proteins): 164,181
peptides, 431 ``non_atg_nterm``; an in-frame translation oracle finds 4 other misses (share
0.99998).

Inputs may be gzip-compressed (§5.3). ``validation/`` is imported from the repository root.
Every run gets a test ``$GIAP_DB_DIR`` (never the real one).
"""

from __future__ import annotations

import gzip
import importlib
import json
import sys
from pathlib import Path
from typing import Any

import pytest
from helpers.mapping_steps import write_test_db

REPO = Path(__file__).resolve().parents[2]
RECORDED = REPO / "tests" / "data" / "recorded" / "bakta" / "1.12.2-full" / "nissle"
NISSLE_JSON = RECORDED / "nissle1917.json.gz"
NISSLE_GENOME = REPO / "tests" / "data" / "genomes" / "nissle1917.fna.gz"
THRESHOLD = 0.999  # §11.6 TE1 proxy
NISSLE_PEPTIDES, NISSLE_NON_ATG_NTERM = 164_181, 431  # lead decisions 4-5

pytestmark = pytest.mark.slow


def _te_proxies() -> Any:
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return importlib.import_module("validation.proxies.te_proxies")


@pytest.fixture(autouse=True)
def _test_db(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GIAP_DB_DIR", str(write_test_db(tmp_path_factory.mktemp("test_db"))))


def _load_json(path: Path) -> dict[str, Any]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        data: dict[str, Any] = json.load(fh)
    return data


def _genes(path: Path) -> list[dict[str, Any]]:
    """CDS/sORF features with a protein, from Bakta 1.12 or giapsim JSON."""
    return [
        f
        for f in _load_json(path)["features"]
        if f["type"] in ("cds", "sorf") and f.get("aa") and not f.get("pseudo")
    ]


def _nterm(aa: str) -> list[str]:
    from giap.mapping.digest import digest

    return [p for s, _, p in digest(aa) if s == 0]


def _check_shape(r: Any) -> None:
    c, sh = r.counts, r.shares
    assert c["peptides"] == c["non_atg_nterm"] + c["exact_inframe"] + c["misses"]
    denom = c["peptides"] - c["non_atg_nterm"]
    assert denom > 0
    assert sh["exact_inframe"] == pytest.approx(c["exact_inframe"] / denom)
    assert r.passed is (sh["exact_inframe"] >= THRESHOLD)
    assert len(r.exceptions) <= 50
    for ex in r.exceptions:
        assert {"locus_tag", "peptide"} <= set(ex)


def _sim(tmp: Path, seed: int, **kw: Any) -> Any:
    from giapsim.simulate import SimConfig, simulate

    return simulate(SimConfig(genome_length=100_000, write_precomputed=True, **kw), seed, tmp)


def _n_peptides(genes: list[dict[str, Any]]) -> int:
    """Distinct digest sequences over every gene (lead decision 4)."""
    from giap.mapping.digest import digest

    return len({p for g in genes for _, _, p in digest(g["aa"])})


def _strict_nterm(genes: list[dict[str, Any]]) -> set[str]:
    """Distinct N-terminal products of genes with a non-ATG first codon and M (decision 5)."""
    return {p for g in genes if g["nt"][:3] != "ATG" and g["aa"][0] == "M" for p in _nterm(g["aa"])}


# ------------------------------------------------------------------ example 4


def test_example_4_all_atg_sim_share_is_one(tmp_path: Path) -> None:
    res = _sim(tmp_path / "sim", 2210, start_codons=(("ATG", 1.0),))
    genes = _genes(res.bakta_json)
    assert genes and all(g["nt"][:3] == "ATG" for g in genes)
    r = _te_proxies().te1_proxy(res.genome_fasta, res.bakta_json)
    _check_shape(r)
    assert r.counts["peptides"] == _n_peptides(genes)
    assert r.counts["non_atg_nterm"] == 0
    assert r.counts["misses"] == 0
    assert r.shares["exact_inframe"] == 1.0
    assert r.passed is True


# ------------------------------------------------------------------ example 5


def test_example_5_gtg_ttg_starts_are_the_only_misses(tmp_path: Path) -> None:
    res = _sim(tmp_path / "sim", 2211)  # default start codons: 90 % ATG, 8 % GTG, 2 % TTG
    genes = _genes(res.bakta_json)
    non_atg = [g for g in genes if g["nt"][:3] != "ATG"]
    assert non_atg, "the simulated genome has no GTG/TTG genes; test is vacuous"
    assert all(g["nt"][:3] in ("GTG", "TTG") and g["aa"][0] == "M" for g in non_atg)
    expected_nterm = {(g["locus"], p) for g in non_atg for p in _nterm(g["aa"])}
    assert expected_nterm

    r = _te_proxies().te1_proxy(res.genome_fasta, res.bakta_json)
    _check_shape(r)
    assert r.counts["misses"] == 0
    assert r.shares["exact_inframe"] == 1.0
    assert r.passed is True
    assert r.counts["peptides"] == _n_peptides(genes)
    assert r.counts["non_atg_nterm"] == len({p for _, p in expected_nterm})
    assert r.counts["non_atg_nterm"] == len(_strict_nterm(genes))
    assert r.exceptions, "the non-ATG N-terminal peptides are listed"
    assert len(r.exceptions) == min(50, r.counts["non_atg_nterm"])  # all listed, up to 50
    for ex in r.exceptions:
        assert (ex["locus_tag"], ex["peptide"]) in expected_nterm, ex


def test_exclusion_is_what_keeps_the_share_at_one(tmp_path: Path) -> None:
    """With 50 % GTG starts the N-terminal products are a visible share of all products: the
    raw exact share is below 1, the share after the §11.6 exclusion is exactly 1."""
    res = _sim(tmp_path / "sim", 2212, start_codons=(("ATG", 0.5), ("GTG", 0.5)))
    r = _te_proxies().te1_proxy(res.genome_fasta, res.bakta_json)
    _check_shape(r)
    c = r.counts
    assert c["non_atg_nterm"] > 0
    assert c["exact_inframe"] / c["peptides"] < 1.0
    assert r.shares["exact_inframe"] == 1.0 and r.passed is True


# ------------------------------------------------------------------ real Nissle (gate M3)


@pytest.mark.pending_sample
@pytest.mark.skipif(
    not (NISSLE_JSON.is_file() and NISSLE_GENOME.is_file()),
    reason=f"pending_sample: no recorded Nissle full-DB Bakta JSON at {NISSLE_JSON}",
)
def test_te1_nissle_full_db_recording_passes() -> None:
    genes = _genes(NISSLE_JSON)
    strict = _strict_nterm(genes)  # decision 5; 5'-truncated "Edge" genes not starting M excluded
    assert _n_peptides(genes) == NISSLE_PEPTIDES  # guards the pinned numbers against the fixture
    assert len(strict) == NISSLE_NON_ATG_NTERM

    r = _te_proxies().te1_proxy(NISSLE_GENOME, NISSLE_JSON)
    _check_shape(r)
    assert r.counts["peptides"] == NISSLE_PEPTIDES
    assert r.counts["non_atg_nterm"] == NISSLE_NON_ATG_NTERM
    assert r.shares["exact_inframe"] >= THRESHOLD
    assert r.passed is True
    assert r.exceptions
    gene_ids = {g["locus"] for g in genes}
    for ex in r.exceptions:
        assert ex["locus_tag"] in gene_ids, ex
    # the listed exceptions include the §11.6 kind: N-terminal peptides at non-ATG starts
    assert any(ex["peptide"] in strict for ex in r.exceptions)
