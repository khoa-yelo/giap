"""gippasim de novo errors and reference peptides (issue M3-20; GIPPA_SPEC.md §11.3, §5.3, §8.4).

Examples 1-8 of docs/specs/M3-20-sim-denovo-errors.md, truth-label checks against the independent
oracle tests/helpers/denovo_oracle.py, the M3-11 consumer contract, determinism and backward
compatibility. Invariants over random peptides are in tests/property/test_sim_denovo_props.py and
the statistical claims in tests/calibration/test_sim_denovo_calibration.py.

Conventions these tests pin where the spec is silent:

* API in ``gippasim.denovo``: ``DenovoConfig``, ``ERROR_TYPES``, ``DEFAULT_ERROR_WEIGHTS``,
  ``BLOCKS``/``NEAR_BLOCKS`` (gippasim's own copy of §7.2.4, pairs of strings), the pure error
  operators ``swap_adjacent(pep, i)``, ``reverse_stretch(pep, start, length)``,
  ``exchange_block(pep, start, old, new)``, ``misread_modification(pep, i)``,
  ``redraw_il(pep, rng)`` and the random operator ``apply_error(pep, error_type, rng)`` (returns
  the erroneous peptide; may raise ``ValueError`` only when the type has no site in ``pep``).
  Positions are 0-based.
* ``SimConfig.denovo`` (default ``None`` = M3-19 behaviour) switches the module on and needs
  ``SimConfig.proteomics`` (else ``ValueError``). ``SimResult.truth`` gains ``"denovo"`` and
  ``"denovo_precision"``.
* ``DenovoConfig(precision_by_decile=<10 floats>, error_weights=<mapping type -> weight>)``:
  ``precision_by_decile[b]`` is the configured precision of score decile b (b = 0 lowest);
  ``truth/denovo_precision.tsv`` states the model's precision for each bin (within 0.02 of the
  configured value). Default weights are the spec's proposal (Questions for the lead, 1).
* Score deciles are pooled over every row of ``peptides_denovo.tsv`` (all runs): edges at the
  k/10 quantiles, right-closed, as ``gippa.mapping.denovo_precision.score_bins``.
* One de novo row per non-decoy row of ``peptides.tsv`` (same ``run``, ``spectrum_id``,
  ``charge``); decoy rows get none. ``true_peptide`` is the database peptide of that spectrum
  whenever the database identification is correct (truth ``is_wrong`` false).
* ``reference_peptides.tsv`` holds every non-decoy row of ``peptides.tsv`` with the same
  ``run``, ``spectrum_id``, ``peptide`` and ``q_value``, ``source = database``,
  ``search_db = same``.
* The de novo stream is separate: switching it on leaves every M3-19 output byte-identical.
"""

from __future__ import annotations

import csv
import hashlib
import math
from collections import Counter
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers.denovo_oracle import (
    AA20,
    ERROR_TYPES,
    MOD_READS,
    SPEC_BLOCKS,
    SPEC_NEAR_BLOCKS,
    decile_bins,
    explains,
    il,
    mass,
    norm_pairs,
)

import gippasim
from gippa.config.model import InputConfig
from gippa.io.peptides import aggregate_peptides, read_peptide_table, validate_peptides

DENOVO_FILES = (
    "peptides_denovo.tsv",
    "reference_peptides.tsv",
    "truth/denovo.tsv",
    "truth/denovo_precision.tsv",
)
TRUTH_COLUMNS = ["run", "spectrum_id", "peptide", "true_peptide", "error_type", "correct"]
PRECISION_COLUMNS = ["score_bin", "precision_true", "n"]
#: M3-19 outputs of SimConfig(genome_length=60_000, proteomics=ProteomicsConfig(pi_wrong=0.05)),
#: seed 11, recorded on main at 482d28a before M3-20 (backward compatibility).
M319_SHA256 = {
    "genome.fasta": "be761b0905ac293822a08f0adb63c7438a254d173ea687cf25ab9785159bdf45",
    "bakta_like.json": "707b8e5f3ca79c722438ae28543aa0992a9b217f7b35028e67e080e82d0f0c70",
    "peptides.tsv": "5e8afae04c8b06a82f0850e88ae9b2e6bd5cfb36a044abd40fc6b68dd466ee63",
    "proteins.faa": "ce1b758fea08ea4af66aa4c98a797df74b2fc54e7905fbb5dfa2e70f4b3445a8",
    "abundance.tsv": "830a9fdcc550b3e99ec492da9a350cd2bc849112a6757054584d260959583504",
    "truth/peptides.tsv": "bc76394203ce2e54db82fdfdc1e9326624b7bb77de04d995ea0959cfb9d55dc1",
    "truth/proteomics.json": "0fd12a1ae51936ebabf933761e8d0d4ab85b1baa8a71260703fbb8f0a8bf10fe",
    "truth/features.tsv": "18f2da469969eadcc322494f3bc082749c0c65ec9e2ef649e4bd59193a68fd4a",
}
M319_FILES = {
    "abundance.tsv", "annotation.gff3", "bakta_like.json", "genome.fasta", "peptides.tsv",
    "proteins.faa", "truth/features.tsv", "truth/genes.tsv", "truth/novel_orfs.tsv",
    "truth/peptides.tsv", "truth/proteomics.json",
}  # fmt: skip
SEED = 11
BASE = {"genome_length": 60_000}
PROPOSED_WEIGHTS = {
    "substitution": 0.3, "block": 0.2, "swap": 0.15, "near_block": 0.05,
    "reversed_stretch": 0.1, "modification": 0.1, "random": 0.1,
}  # fmt: skip
EX6_PRECISION = tuple(float(x) for x in np.round(np.linspace(0.45, 0.9, 10), 6))


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def truthy(v: str) -> bool:
    return v.strip().lower() in ("true", "1")


def dn(**kw: Any) -> Any:
    from gippasim.denovo import DenovoConfig

    return DenovoConfig(**kw)


def prot(**kw: Any) -> Any:
    from gippasim.proteomics import ProteomicsConfig

    return ProteomicsConfig(**kw)


@dataclass
class DnRun:
    out: Path
    res: Any = None
    meta: dict[str, Any] = field(default_factory=dict)

    @cached_property
    def db(self) -> list[dict[str, str]]:
        return rows(self.out / "peptides.tsv")

    @cached_property
    def db_targets(self) -> list[dict[str, str]]:
        return [r for r in self.db if not truthy(r["is_decoy"])]

    @cached_property
    def db_truth(self) -> dict[str, bool]:
        """peptide -> is_wrong (non-decoy truth rows)."""
        out: dict[str, bool] = {}
        for r in rows(self.out / "truth" / "peptides.tsv"):
            if r["origin"] != "decoy":
                out[r["peptide"]] = out.get(r["peptide"], False) or truthy(r["is_wrong"])
        return out

    @cached_property
    def denovo(self) -> list[dict[str, str]]:
        return rows(self.out / "peptides_denovo.tsv")

    @cached_property
    def reference(self) -> list[dict[str, str]]:
        return rows(self.out / "reference_peptides.tsv")

    @cached_property
    def truth(self) -> list[dict[str, str]]:
        return rows(self.out / "truth" / "denovo.tsv")

    @cached_property
    def precision(self) -> list[dict[str, str]]:
        return rows(self.out / "truth" / "denovo_precision.tsv")

    @cached_property
    def truth_by_key(self) -> dict[tuple[str, str], dict[str, str]]:
        return {(r["run"], r["spectrum_id"]): r for r in self.truth}

    def scores_and_correct(self) -> tuple[np.ndarray, np.ndarray]:
        s = np.array([float(r["score"]) for r in self.denovo])
        c = np.array(
            [
                truthy(self.truth_by_key[(r["run"], r["spectrum_id"])]["correct"])
                for r in self.denovo
            ]
        )
        return s, c


def simulate(out: Path, seed: int, **kw: Any) -> DnRun:
    cfg = gippasim.SimConfig(**kw)
    res = gippasim.simulate(cfg, seed, out)
    return DnRun(out, res, {"cfg": cfg, "seed": seed})


@pytest.fixture(scope="module")
def mk(tmp_path_factory: pytest.TempPathFactory) -> Any:
    cache: dict[str, DnRun] = {}

    def make(name: str, seed: int, **kw: Any) -> DnRun:
        if name not in cache:
            cache[name] = simulate(tmp_path_factory.mktemp(name), seed, **kw)
        return cache[name]

    return make


@pytest.fixture(scope="module")
def default(mk: Any) -> DnRun:
    run: DnRun = mk("default", SEED, proteomics=prot(pi_wrong=0.05), denovo=dn(), **BASE)
    return run


@pytest.fixture(scope="module")
def big(mk: Any) -> DnRun:
    """Example 6: >= 10^4 spectra, configured precision 0.9 in the top decile."""
    run: DnRun = mk(
        "big",
        2026,
        proteomics=prot(pi_wrong=0.05, n_runs=5),
        denovo=dn(precision_by_decile=EX6_PRECISION),
        genome_length=100_000,
    )
    return run


# ---------------------------------------------------------------- examples 1-5: error operators


def test_example1_block_n_to_gg() -> None:
    from gippasim.denovo import exchange_block

    out = exchange_block("PEPTNDEKAR", 4, "N", "GG")
    assert out == "PEPTGGDEKAR"
    assert explains("PEPTNDEKAR", out, "block")
    assert abs(mass(out) - mass("PEPTNDEKAR")) < 0.001


def test_block_reverse_direction_and_two_residue_blocks() -> None:
    from gippasim.denovo import exchange_block

    assert exchange_block("PEPTGGDEKAR", 4, "GG", "N") == "PEPTNDEKAR"
    assert exchange_block("PEGTDEKAR", 2, "GT", "SA") == "PESADEKAR"
    with pytest.raises(ValueError):
        exchange_block("PEPTNDEKAR", 3, "N", "GG")  # no N at position 3


def test_example2_adjacent_swap() -> None:
    from gippasim.denovo import swap_adjacent

    assert swap_adjacent("AKLMNPQRSTVW", 4) == "AKLMPNQRSTVW"
    assert explains("AKLMNPQRSTVW", "AKLMPNQRSTVW", "swap")


def test_example3_il_redraw_counted_correct() -> None:
    from gippasim.denovo import redraw_il

    seen = set()
    for s in range(64):
        out = redraw_il("PEPTIDEKAR", np.random.default_rng(s))
        assert out in ("PEPTIDEKAR", "PEPTLDEKAR"), out
        seen.add(out)
    assert seen == {"PEPTIDEKAR", "PEPTLDEKAR"}  # I and L are re-drawn at random
    assert explains("PEPTIDEKAR", "PEPTLDEKAR", "none")


def test_example3_truth_marks_il_variants_correct(default: DnRun) -> None:
    none = [r for r in default.truth if r["error_type"] == "none"]
    assert none
    for r in none:
        assert il(r["peptide"]) == il(r["true_peptide"]), r
        assert truthy(r["correct"]), r
    with_il = [r for r in none if set(r["true_peptide"]) & {"I", "L"}]
    assert any(r["peptide"] != r["true_peptide"] for r in with_il)  # the redraw happens


def test_example4_reversed_stretch() -> None:
    from gippasim.denovo import reverse_stretch

    assert reverse_stretch("MKAESTVLLR", 2, 4) == "MKTSEAVLLR"
    assert explains("MKAESTVLLR", "MKTSEAVLLR", "reversed_stretch")
    for bad in (2, 7):  # stretches are 3-6 residues
        with pytest.raises(ValueError):
            reverse_stretch("MKAESTVLLR", 0, bad)


def test_example5_modification_misread() -> None:
    from gippasim.denovo import misread_modification

    assert misread_modification("PEPTADEKAR", 4) == "PEPTSDEKAR"  # A + Oxidation = S
    assert misread_modification("PEPTNDEKAR", 4) == "PEPTDDEKAR"  # N + Deamidated = D
    assert misread_modification("PEPTGDEKAR", 4) == "PEPTADEKAR"  # G + Methyl = A
    with pytest.raises(ValueError):
        misread_modification("PEPTKDEKAR", 4)  # K has no modelled misread
    assert explains("PEPTADEKAR", "PEPTSDEKAR", "modification")


def test_modification_reads_are_mass_equivalent() -> None:
    from helpers.denovo_oracle import DEAMIDATED, METHYL, OXIDATION, RESIDUE_MASS

    delta = {("A", "S"): OXIDATION, ("N", "D"): DEAMIDATED, ("G", "A"): METHYL}
    for (a, b), d in delta.items():
        assert abs(RESIDUE_MASS[a] + d - RESIDUE_MASS[b]) < 0.001
    from gippasim.denovo import misread_modification

    for a, b in MOD_READS:
        assert misread_modification(f"PEPT{a}DEKAR", 4) == f"PEPT{b}DEKAR"


# ---------------------------------------------------------------- block list (§7.2.4)


def test_block_lists_equal_spec_and_diffops() -> None:
    from gippa.mapping import diffops
    from gippasim.denovo import BLOCKS, NEAR_BLOCKS

    spec_b = norm_pairs(tuple(p) for p in SPEC_BLOCKS)
    spec_n = norm_pairs(tuple(p) for p in SPEC_NEAR_BLOCKS)
    assert norm_pairs(BLOCKS) == spec_b and len(spec_b) == 5
    assert norm_pairs(NEAR_BLOCKS) == spec_n and len(spec_n) == 3
    assert norm_pairs(diffops.BLOCKS) == spec_b
    assert norm_pairs(diffops.NEAR_BLOCKS) == spec_n


def test_block_lists_are_gippasims_own_copy() -> None:
    import gippasim.denovo as d
    from gippa.mapping import diffops

    src = Path(d.__file__).read_text()
    assert "gippa.mapping" not in src  # independent of the code under test
    assert d.BLOCKS is not diffops.BLOCKS


# ---------------------------------------------------------------- config


def test_error_types_and_default_weights() -> None:
    from gippasim.denovo import DEFAULT_ERROR_WEIGHTS
    from gippasim.denovo import ERROR_TYPES as SIM_TYPES

    assert set(SIM_TYPES) == set(ERROR_TYPES)
    assert dict(DEFAULT_ERROR_WEIGHTS) == pytest.approx(PROPOSED_WEIGHTS)
    assert dict(dn().error_weights) == pytest.approx(PROPOSED_WEIGHTS)
    p = dn().precision_by_decile
    assert len(p) == 10 and all(0.0 <= x <= 1.0 for x in p)
    assert list(p) == sorted(p)  # scores correlate with correctness


@pytest.mark.parametrize(
    "kw",
    [
        {"precision_by_decile": (0.5,) * 9},
        {"precision_by_decile": (0.5,) * 9 + (1.2,)},
        {"precision_by_decile": (-0.1,) + (0.5,) * 9},
        {"error_weights": {"substitution": -1.0, "swap": 2.0}},
        {"error_weights": {"substitution": 0.0}},
        {"error_weights": {"frameshift": 1.0}},
    ],
)
def test_config_rejects_invalid(kw: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        dn(**kw)


def test_denovo_needs_proteomics(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        gippasim.simulate(gippasim.SimConfig(denovo=dn(), **BASE), SEED, tmp_path)


# ---------------------------------------------------------------- outputs and truth labels


def test_outputs_and_truth_keys(default: DnRun) -> None:
    for f in DENOVO_FILES:
        assert (default.out / f).is_file(), f
    assert Path(default.res.truth["denovo"]) == default.out / "truth" / "denovo.tsv"
    assert (
        Path(default.res.truth["denovo_precision"])
        == default.out / "truth" / "denovo_precision.tsv"
    )
    with (default.out / "truth" / "denovo.tsv").open() as fh:
        assert fh.readline().rstrip("\n").split("\t")[: len(TRUTH_COLUMNS)] == TRUTH_COLUMNS
    with (default.out / "truth" / "denovo_precision.tsv").open() as fh:
        assert fh.readline().rstrip("\n").split("\t")[:3] == PRECISION_COLUMNS


def test_one_denovo_row_per_target_spectrum(default: DnRun) -> None:
    targets = {(r["run"], r["spectrum_id"]): r for r in default.db_targets}
    decoys = {(r["run"], r["spectrum_id"]) for r in default.db if truthy(r["is_decoy"])}
    keys = [(r["run"], r["spectrum_id"]) for r in default.denovo]
    assert len(keys) == len(set(keys))
    assert set(keys) == set(targets)
    assert not set(keys) & decoys
    for r in default.denovo:
        assert r["source"] == "denovo", r
        assert r["charge"] == targets[(r["run"], r["spectrum_id"])]["charge"], r
        assert not truthy(r.get("is_decoy", "false")), r
        s = float(r["score"])
        assert math.isfinite(s), r


def test_truth_rows_match_denovo_rows(default: DnRun) -> None:
    assert len(default.truth) == len(default.denovo)
    assert set(default.truth_by_key) == {(r["run"], r["spectrum_id"]) for r in default.denovo}
    for r in default.denovo:
        t = default.truth_by_key[(r["run"], r["spectrum_id"])]
        assert t["peptide"] == r["peptide"], (r, t)


def test_true_peptide_is_the_correct_database_id(default: DnRun) -> None:
    db = {(r["run"], r["spectrum_id"]): r["peptide"] for r in default.db_targets}
    n = 0
    for t in default.truth:
        pep = db[(t["run"], t["spectrum_id"])]
        if not default.db_truth[pep]:
            assert t["true_peptide"] == pep, t
            n += 1
    assert n > 0


def test_truth_labels_explained_by_oracle(default: DnRun) -> None:
    types = Counter()
    for t in default.truth:
        et, pep, true = t["error_type"], t["peptide"], t["true_peptide"]
        assert et in {"none", *ERROR_TYPES}, t
        types[et] += 1
        assert set(pep) <= set(AA20), t
        assert 7 <= len(pep) <= 40, t
        assert truthy(t["correct"]) == (il(pep) == il(true)), t
        assert truthy(t["correct"]) == (et == "none"), t
        assert explains(true, pep, et), t
        if et in ("block", "near_block"):
            assert abs(mass(pep) - mass(true)) <= 0.05, t
    assert types["none"] > 0 and sum(v for k, v in types.items() if k != "none") > 0


def test_every_error_type_occurs(big: DnRun) -> None:
    types = Counter(t["error_type"] for t in big.truth)
    for et in ERROR_TYPES:
        assert types[et] > 0, (et, types)


@pytest.mark.parametrize("et", ["substitution", "swap", "reversed_stretch", "random", "block"])
def test_single_error_type_weight(et: str, mk: Any) -> None:
    w = {k: (1.0 if k == et else 0.0) for k in ERROR_TYPES}
    run: DnRun = mk(
        f"only_{et}", 5, proteomics=prot(pi_wrong=0.05), denovo=dn(error_weights=w), **BASE
    )
    errs = [t for t in run.truth if t["error_type"] != "none"]
    assert errs
    if et == "block":  # peptides without a block site may fall back to another type
        share = sum(t["error_type"] == "block" for t in errs) / len(errs)
        assert share >= 0.5, share
    else:
        assert {t["error_type"] for t in errs} == {et}
    for t in errs:
        assert explains(t["true_peptide"], t["peptide"], t["error_type"]), t


# ---------------------------------------------------------------- precision bins


def test_precision_table_shape(default: DnRun) -> None:
    bins = [int(r["score_bin"]) for r in default.precision]
    assert bins == list(range(10))
    n = [int(r["n"]) for r in default.precision]
    assert sum(n) == len(default.denovo)
    mine = np.bincount(decile_bins(default.scores_and_correct()[0]), minlength=10)
    tol = max(2, int(0.01 * len(default.denovo)))
    assert all(abs(a - b) <= tol for a, b in zip(n, mine.tolist(), strict=True)), (n, mine)
    for r, p in zip(default.precision, dn().precision_by_decile, strict=True):
        x = float(r["precision_true"])
        assert 0.0 <= x <= 1.0
        assert abs(x - p) <= 0.02, (r, p)


def test_bins_agree_with_m311_score_bins(default: DnRun) -> None:
    from gippa.mapping.denovo_precision import score_bins

    s, _ = default.scores_and_correct()
    np.testing.assert_array_equal(score_bins(s), decile_bins(s))


def test_example6_top_decile_precision(big: DnRun) -> None:
    s, c = big.scores_and_correct()
    assert len(s) >= 10_000
    top = decile_bins(s) == 9
    assert abs(float(c[top].mean()) - 0.9) <= 0.03, c[top].mean()
    row = next(r for r in big.precision if int(r["score_bin"]) == 9)
    assert abs(float(row["precision_true"]) - 0.9) <= 0.02, row


def test_scores_correlate_with_correctness(big: DnRun) -> None:
    from scipy import stats

    s, c = big.scores_and_correct()
    assert s[c].mean() > s[~c].mean()
    assert stats.mannwhitneyu(s[c], s[~c], alternative="greater").pvalue < 1e-6
    b = decile_bins(s)
    prec = [float(c[b == k].mean()) for k in range(10)]
    assert prec[9] > prec[0] + 0.3, prec


# ---------------------------------------------------------------- example 7: reference peptides


def test_example7_reference_peptides_pair_spectra(default: DnRun) -> None:
    ref = {(r["run"], r["spectrum_id"]): r for r in default.reference}
    for r in default.reference:
        assert r["source"] == "database" and r["search_db"] == "same", r
    for r in default.db_targets:
        k = (r["run"], r["spectrum_id"])
        assert k in ref, k
        assert ref[k]["peptide"] == r["peptide"] and ref[k]["q_value"] == r["q_value"], k
    db = {(r["run"], r["spectrum_id"]): r["peptide"] for r in default.db_targets}
    for r in default.denovo:
        k = (r["run"], r["spectrum_id"])
        t = default.truth_by_key[k]
        if not default.db_truth[db[k]]:  # the spectrum had a correct database ID
            assert k in ref and ref[k]["peptide"] == t["true_peptide"], k


# ---------------------------------------------------------------- example 8: M1-06 validity


@pytest.mark.parametrize("name", ["peptides_denovo.tsv", "reference_peptides.tsv"])
def test_example8_tables_pass_m106(name: str, default: DnRun) -> None:
    df_raw, issues = read_peptide_table(default.out / name)
    df, issues2 = validate_peptides(df_raw, InputConfig())
    errors = [i for i in [*issues, *issues2] if i.severity == "error"]
    assert not errors, errors
    assert len(df) == len(rows(default.out / name))
    want = "denovo" if name == "peptides_denovo.tsv" else "database"
    assert set(df["source"]) == {want}
    agg, _ = aggregate_peptides(df)
    assert len(agg) > 0


# ---------------------------------------------------------------- determinism, compatibility


def test_deterministic(tmp_path: Path, default: DnRun) -> None:
    again = simulate(tmp_path, SEED, proteomics=prot(pi_wrong=0.05), denovo=dn(), **BASE)
    for f in DENOVO_FILES:
        assert sha(again.out / f) == sha(default.out / f), f


def test_seed_changes_output(tmp_path: Path, default: DnRun) -> None:
    other = simulate(tmp_path, SEED + 1, proteomics=prot(pi_wrong=0.05), denovo=dn(), **BASE)
    assert sha(other.out / "truth" / "denovo.tsv") != sha(default.out / "truth" / "denovo.tsv")


def test_off_by_default_unchanged(tmp_path: Path) -> None:
    run = simulate(tmp_path, SEED, proteomics=prot(pi_wrong=0.05), **BASE)
    files = {p.relative_to(tmp_path).as_posix() for p in tmp_path.rglob("*") if p.is_file()}
    assert files == M319_FILES
    for f, h in M319_SHA256.items():
        assert sha(run.out / f) == h, f
    assert gippasim.SimConfig().denovo is None


def test_switching_on_leaves_m319_outputs_unchanged(default: DnRun) -> None:
    for f, h in M319_SHA256.items():
        assert sha(default.out / f) == h, f
