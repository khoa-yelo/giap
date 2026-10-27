"""M7-10 in the ``variants`` step: ``class_error_rate`` / ``class_error_n`` on every
``variants.tsv`` ALT row and ``CLASS_ERR`` in the VCF, from C17 ``peptide_spectra``
(docs/specs/M7-10-class-error-rates.md rules 1-4, examples 1-6 at step level, degenerate rows;
lead decision 1; GIPPA_SPEC.md §7.5.4, §7.5.5).

Interface pinned here (test-author reading; see the report to the lead):

* The step declares ``peptide_spectra`` among its ``optional_inputs`` (lead decision 1); its
  required inputs are unchanged. Without it every ALT row is ``unestimated`` / ``NA``.
* With it, the rates are pooled per class ``(alt_type, mass_equivalent)``, where a site's
  ``mass_equivalent`` is its M7-05 flag, and every row of a class carries the same value.
  The expected values come from ``helpers.class_error_chain.oracle_rates``: a brute-force count
  of rules 1-3 on the raw gippasim tables (de novo reads, reference peptides with q and decoy
  status), with the sites' supporting peptides from the merged M7-02 ``sites``.
* ``reference_peptides`` without ``spectrum_id``: every row ``unestimated`` and a WARNING naming
  ``no_spectrum_pairing``, either as a preflight check (rule id or message) or in the step's log.
* A C17 table without de novo rows (database-only input): every row ``unestimated``.

The fixture is a 50-kb gippasim genome, M3-19 database peptides (pi_wrong 0.05) and their M3-20
de novo reads as the input, with gippasim's paired ``reference_peptides.tsv`` (183 of ~3,000
reference rows have q > 0.01 at this size, so rule 1's threshold is exercised).
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.class_error_chain import (
    ClassErrorChain,
    build_chain,
    library_preflight,
    measured_rates,
    oracle_rates,
    oracle_sites,
    peptide_classes,
    read_variants,
    run_step,
    table_classes,
)

SEED = 23
UNEST = "unestimated"


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> ClassErrorChain:
    return build_chain(Path(tmp_path_factory.mktemp("m7_10_chain")), SEED)


@pytest.fixture(scope="module")
def oracle(chain: ClassErrorChain) -> dict[str, Any]:
    st, peps = oracle_sites(chain.avail, chain.cfg)
    classes = peptide_classes(st, peps)
    return {
        "sites": st,
        "classes": classes,
        "rates": oracle_rates(classes, chain.denovo_tsv, chain.reference_tsv),
        "rates_no_q": oracle_rates(classes, chain.denovo_tsv, chain.reference_tsv, q_max=1.0),
    }


@pytest.fixture(scope="module")
def with_spectra(chain: ClassErrorChain, tmp_path_factory: pytest.TempPathFactory) -> Any:
    out = Path(tmp_path_factory.mktemp("m7_10_step"))
    paths, recs = run_step(chain.step_inputs(), out, chain.cfg)
    return paths, recs


def _all_unestimated(variants: pd.DataFrame) -> None:
    assert len(variants) > 0
    assert set(variants["class_error_rate"]) == {UNEST}
    assert set(variants["class_error_n"]) <= {"NA", ""}


def test_step_declares_peptide_spectra_optional() -> None:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    spec = reg.get("variants")
    assert "peptide_spectra" in spec.optional_inputs
    assert "peptide_spectra" not in spec.inputs
    assert set(spec.inputs) == {"mapping", "genome", "genes", "replicons", "peptides"}


def test_fixture_has_several_classes_and_q_matters(oracle: dict[str, Any]) -> None:
    """Fixture sanity (independent of the code under test)."""
    rates = oracle["rates"]
    assert len(rates) >= 3, rates
    assert sum(n for n, _e in rates.values()) >= 100
    assert rates != oracle["rates_no_q"], "rule 1's q threshold changes nothing here"


def test_rates_equal_the_brute_force_oracle(with_spectra: Any, oracle: dict[str, Any]) -> None:
    paths, _ = with_spectra
    v = read_variants(paths["variants"])
    assert v["variant_id"].tolist() == oracle["sites"]["variant_id"].tolist()
    want = oracle["rates"]
    got = measured_rates(v)
    assert set(got) == set(want), (sorted(got), sorted(want))
    for c, (n, e) in want.items():
        assert got[c][1] == n, c
        assert got[c][0] == pytest.approx(e / n, rel=1e-5, abs=1e-6), c
    for c, r in zip(table_classes(v), v["class_error_rate"], strict=True):
        assert (r == UNEST) == (c not in want), c


def test_q_above_threshold_is_not_a_pair(with_spectra: Any, oracle: dict[str, Any]) -> None:
    """Degenerate row "reference q above 0.01": the step's ``n`` is the q <= 0.01 count, not the
    count over every reference row."""
    paths, _ = with_spectra
    got = measured_rates(read_variants(paths["variants"]))
    diff = [c for c in oracle["rates"] if oracle["rates"][c][0] != oracle["rates_no_q"][c][0]]
    assert diff, "fixture: no class where the q filter matters"
    for c in diff:
        assert got[c][1] == oracle["rates"][c][0] != oracle["rates_no_q"][c][0]


def test_vcf_class_err_matches_the_table(with_spectra: Any) -> None:
    """Example 6 at step level: each VCF allele's ``CLASS_ERR`` is its row's value."""
    paths, _ = with_spectra
    v = read_variants(paths["variants"])
    vcf_rows = v[v["alt_codon_vcf"] != "NA"]
    vals: list[str] = []
    for line in paths["variants_vcf"].read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        info = dict(kv.split("=", 1) for kv in line.split("\t")[7].split(";") if "=" in kv)
        vals += info["CLASS_ERR"].split(",")
    assert sorted(vals) == sorted(vcf_rows["class_error_rate"].tolist())
    assert any(x != UNEST for x in vals)


def test_without_peptide_spectra_every_row_unestimated(
    chain: ClassErrorChain, tmp_path: Path
) -> None:
    paths, _ = run_step(chain.step_inputs(with_spectra=False), tmp_path, chain.cfg)
    _all_unestimated(read_variants(paths["variants"]))


def _strip_spectrum_ids(src: Path, dst: Path) -> Path:
    with src.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cols = [c for c in rows[0] if c != "spectrum_id"]
    with dst.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")  # fmt: skip
        w.writeheader()
        w.writerows(rows)
    return dst


def test_reference_without_spectrum_id_is_unestimated_with_warning(
    chain: ClassErrorChain, tmp_path: Path
) -> None:
    ref = _strip_spectrum_ids(chain.reference_tsv, tmp_path / "ref_nosid.tsv")
    res = chain.res
    report, work = library_preflight(res.genome_fasta, chain.denovo_tsv, res.bakta_json, ref,
                                     tmp_path / "pf", res.config_yaml)  # fmt: skip
    inputs = chain.step_inputs(with_spectra=False)
    inputs["reference_peptides"] = ref
    if (work / "peptide_spectra.parquet").is_file():
        inputs["peptide_spectra"] = work / "peptide_spectra.parquet"
    paths, recs = run_step(inputs, tmp_path / "out", chain.cfg)
    _all_unestimated(read_variants(paths["variants"]))
    texts = [f"{c.rule_id} {c.message}" for c in report.checks if c.level in ("warn", "error")]
    texts += [r.getMessage() for r in recs if r.levelno >= 30]
    assert any("no_spectrum_pairing" in t for t in texts), texts


def test_database_only_spectra_are_unestimated(chain: ClassErrorChain, tmp_path: Path) -> None:
    """Degenerate row "database-only input": a C17 table with no de novo row (preflight on the
    database peptides of the same spectra) leaves every row ``unestimated``."""
    res = chain.res
    _report, work = library_preflight(res.genome_fasta, Path(res.peptides), res.bakta_json,
                                      chain.reference_tsv, tmp_path / "pf",
                                      res.config_yaml)  # fmt: skip
    spectra = work / "peptide_spectra.parquet"
    assert spectra.is_file()
    assert not (pd.read_parquet(spectra)["source"].astype(str) == "denovo").any()
    inputs = chain.step_inputs(with_spectra=False)
    inputs["peptide_spectra"] = spectra
    paths, _ = run_step(inputs, tmp_path / "out", chain.cfg)
    _all_unestimated(read_variants(paths["variants"]))


def test_rates_do_not_change_calls_or_counts(with_spectra: Any, chain: ClassErrorChain,
                                            tmp_path: Path) -> None:  # fmt: skip
    """The rates annotate rows only (§7.5.4: a rate qualifies a call, it never changes it)."""
    paths, _ = with_spectra
    base, _ = run_step(chain.step_inputs(with_spectra=False), tmp_path, chain.cfg)
    a, b = read_variants(paths["variants"]), read_variants(base["variants"])
    keep = [c for c in a.columns if c not in ("class_error_rate", "class_error_n")]
    pd.testing.assert_frame_equal(a[keep], b[keep])
    assert paths["pileup"].read_bytes() == base["pileup"].read_bytes()
