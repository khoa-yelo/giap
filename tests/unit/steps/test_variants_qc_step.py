"""M7-14 rules 3-4 at step level: the ``variants`` step's ``variants_qc`` output and the
``id_mixups`` row ``context_stats`` writes from it (docs/specs/M7-14-te2-proxy-and-id-mixups.md
examples 1 and 5, the degenerate rows, acceptance "the M4-16 ``requires_variants`` path is still
covered for runs without peptides"; GIAP_SPEC.md §7.5.1 row 1, §9.3 QC rows).

Interface pinned here (test-author reading):

* ``@step("variants")`` gains the output key ``variants_qc`` (one Parquet file, the local
  executor's ``work/variants_qc``); its columns are the implementer's, only ``context_stats``
  reads it. ``context_stats`` gains the optional input ``variants_qc``, so the registry orders
  ``variants`` before ``context_stats``.
* With ``variants_qc``, ``detection_summary.tsv`` has one ``qc``/``run``/``all``/``id_mixups`` row
  whose ``value`` is the count and ``n`` the number of same-db ALT sites (rule 3), with no
  ``requires_variants`` flag; without it the M4-16 row stays (``NA``, flag
  ``requires_variants``).
* Oracle on the fixture: every peptide is ``database``, ``search_db = same`` (giapsim's default),
  so every ``variants.tsv`` row is a same-db ALT site: ``n`` = rows, count = rows flagged neither
  ``start_codon`` nor ``stop_codon``.

Fixtures: the 50-kb giapsim producer chain of tests/unit/steps/test_variants_step.py (strain
layer: sample-strain peptides on the reference genome under the default ``genome_relation:
related``, so strain differences are same-db ALT sites that ``id_mixups`` must count, degenerate
row "related | still counted"); and the M4-20 ``sim1000`` contract-table fixture for
``context_stats``. A test ``$GIAP_DB_DIR`` only. Example 5 on a real ``giap run`` of sim100k is in
tests/integration/test_id_mixups_sim100k.py.
"""

from __future__ import annotations

import dataclasses
import hashlib
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.mapping_steps import drive, fresh_steps, preflight, roots_of

from giap.config import GiapConfig, resolve_config

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import m4_chain as mc

SEED = 11
M4_SEED = 20261021
M4_STATS = {"permutations": 199, "bootstrap": 200, "alpha": {"identification": 0.05}}
STEP_KEYS = ("mapping", "genome", "genes", "replicons", "peptides", "hits_tolerant")
OUTPUTS = ("pileup", "variants", "variants_vcf", "variants_qc")
EXEMPT = {"start_codon", "stop_codon"}


def _cfg(**over: Any) -> GiapConfig:
    return GiapConfig.model_validate({"seed": SEED, "sample_name": "S1", **over})


def _chain(root: Path, *, strain: bool) -> dict[str, Path]:
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SimConfig, simulate
    from giapsim.strain import StrainConfig

    st = StrainConfig(snv_rate_per_kb=5.0, indel_rate_per_kb=0.5, indel_len=(1, 1))
    cfg_sim = SimConfig(genome_length=50_000, proteomics=ProteomicsConfig(),
                        strain=st if strain else None, write_precomputed=False)  # fmt: skip
    res = simulate(cfg_sim, SEED, root / "sim")
    if strain:
        res = dataclasses.replace(
            res, genome_fasta=res.reference_fasta, bakta_json=res.reference_json
        )
    run = preflight(res, root / "pf")
    cfg = resolve_config(run.config_yaml, [f"seed={SEED}"])
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, root / "work", roots_of(run))
    return {k: avail[k] for k in STEP_KEYS}


def _run_variants(inputs: Mapping[str, Path], out: Path, cfg: GiapConfig) -> dict[str, Path]:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    out.mkdir(parents=True, exist_ok=True)
    paths = {k: out / k for k in OUTPUTS}
    reg.run_step("variants", cfg, dict(inputs), paths)
    return paths


def _variants(path: Path) -> pd.DataFrame:
    from giap.io.variant_tables import read_variants

    return read_variants(path)


def _oracle(variants: pd.DataFrame) -> tuple[int, int]:
    """All peptides are database/same: every row is a same-db ALT site."""
    n = len(variants)
    count = sum(not (EXEMPT & set(f)) for f in variants["flags"])
    return count, n


@pytest.fixture(scope="module")
def strain_chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return _chain(Path(tmp_path_factory.mktemp("qc_strain")), strain=True)


@pytest.fixture(scope="module")
def same_chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return _chain(Path(tmp_path_factory.mktemp("qc_same")), strain=False)


@pytest.fixture(scope="module")
def strain_run(strain_chain: dict[str, Path], tmp_path_factory: pytest.TempPathFactory
               ) -> dict[str, Path]:  # fmt: skip
    return _run_variants(strain_chain, Path(tmp_path_factory.mktemp("qc_strain_out")), _cfg())


@pytest.fixture(scope="module")
def other_run(strain_chain: dict[str, Path], tmp_path_factory: pytest.TempPathFactory
              ) -> dict[str, Path]:  # fmt: skip
    """The same chain with every C9 peptide relabelled ``search_db = other``."""
    root = Path(tmp_path_factory.mktemp("qc_other"))
    peps = pd.read_parquet(strain_chain["peptides"])
    assert set(peps["source"]) == {"database"}
    peps["search_db"] = "other"
    path = root / "peptides"
    peps.to_parquet(path, index=False)
    return _run_variants({**strain_chain, "peptides": path}, root / "out", _cfg())


@pytest.fixture(scope="module")
def m4_work(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    """M4-20 ``sim1000`` contract tables with the gene table and detectability parts written."""
    from giap.steps.detectability import run_detectability
    from giap.steps.gene_table import run_gene_table

    root = Path(tmp_path_factory.mktemp("qc_m4"))
    db = mc.test_db(root / "db")
    work = mc.build(root / "work", seed=M4_SEED, size="sim1000")
    with mc.db_env(db):
        run_gene_table(work, _m4_cfg())
        run_detectability(work, _m4_cfg())
    return work, db


def _m4_cfg() -> Any:
    return mc.config(M4_SEED, stats=M4_STATS)


def _context_inputs(work: Path) -> dict[str, Path]:
    """The step inputs ``run_context_stats`` would build on ``work`` (M4 half only)."""
    from giap.io.work_tables import OPTIONAL_KEYS, default_paths
    from giap.steps.gene_table import contract_paths

    inputs = dict(contract_paths(work))
    inputs |= {k: p for k, p in default_paths(work).items()
               if k in OPTIONAL_KEYS and k not in inputs and p.exists()}  # fmt: skip
    inputs["gene_table"] = work / "genetable" / "genes.parquet"
    parts = work / "genetable" / "detection_summary_parts"
    for key, name in (("detectability_summary", "detectability"), ("tier_summary", "tier")):
        if (parts / f"{name}.parquet").exists():
            inputs[key] = parts / f"{name}.parquet"
    return inputs


def _context_stats(m4: tuple[Path, Path], out: Path, variants_qc: Path | None) -> pd.DataFrame:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    work, db = m4
    load_step_modules()
    inputs = _context_inputs(work)
    if variants_qc is not None:
        inputs["variants_qc"] = variants_qc
    out.mkdir(parents=True, exist_ok=True)
    outputs = {"regions_stats": out / "regions_stats", "detection_summary": out / "ds"}
    with mc.db_env(db):
        reg.run_step("context_stats", _m4_cfg(), inputs, outputs)
    return mc.read(outputs["detection_summary"])


def _id_mixups_row(ds: pd.DataFrame) -> pd.Series:
    sel = ds[(ds["section"] == "qc") & (ds["group_type"] == "run") & (ds["group"] == "all")
             & (ds["metric"] == "id_mixups")]  # fmt: skip
    assert len(sel) == 1, sel
    return sel.iloc[0]


def _flag_set(v: object) -> set[str]:
    return {x for x in str(v).split(";") if x} if isinstance(v, str) else set()


# ------------------------------------------------------------------ registry (rule 4)


def test_variants_step_gains_the_variants_qc_output() -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get("variants")
    assert "variants_qc" in spec.outputs
    assert "variants_qc" not in spec.optional_outputs


def test_context_stats_takes_variants_qc_as_an_optional_input() -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get("context_stats")
    assert "variants_qc" in spec.optional_inputs
    assert "variants_qc" not in spec.inputs


def test_context_stats_now_runs_after_variants() -> None:
    """Rule 4 / lead decision 2: CONTEXT waits for VARIANTS."""
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    order = reg.topo_order([s.name for s in reg.all_steps()])
    assert order.index("variants") < order.index("context_stats")


# ------------------------------------------------------------------ the variants_qc file


def test_variants_qc_is_written_as_parquet(strain_run: dict[str, Path]) -> None:
    path = strain_run["variants_qc"]
    assert path.is_file()
    pd.read_parquet(path)  # one Parquet file


def test_variants_qc_is_byte_identical_for_the_same_inputs(
    strain_chain: dict[str, Path], strain_run: dict[str, Path], tmp_path: Path
) -> None:
    """TE7: one pass, no randomness."""
    again = _run_variants(strain_chain, tmp_path / "again", _cfg())
    digest = [hashlib.sha256(p["variants_qc"].read_bytes()).hexdigest()
              for p in (strain_run, again)]  # fmt: skip
    assert digest[0] == digest[1]


def test_variants_qc_leaves_the_other_outputs_unchanged(
    strain_chain: dict[str, Path], strain_run: dict[str, Path], tmp_path: Path
) -> None:
    """The extra output adds no row and no flag to ``variants.tsv``."""
    from giap.steps import registry as reg

    paths = {k: tmp_path / k for k in OUTPUTS}
    reg.run_step("variants", _cfg(), dict(strain_chain), paths)
    assert paths["variants"].read_bytes() == strain_run["variants"].read_bytes()


# ------------------------------------------------------------------ example 5 at step level


def test_example_5_id_mixups_row_is_numeric_and_matches_the_oracle(
    strain_run: dict[str, Path], m4_work: tuple[Path, Path], tmp_path: Path
) -> None:
    variants = _variants(strain_run["variants"])
    count, n = _oracle(variants)
    assert n > 0, "the strain fixture is meant to give same-db ALT sites"
    assert count > 0, "strain differences searched as search_db = same are mix-ups"
    row = _id_mixups_row(_context_stats(m4_work, tmp_path / "cs", strain_run["variants_qc"]))
    assert row["value"] != "NA" and row["n"] != "NA", row.to_dict()
    assert float(row["value"]) == count
    assert int(row["n"]) == n
    assert "requires_variants" not in _flag_set(row["flags"])


def test_id_mixups_row_keeps_count_and_n_apart(
    strain_run: dict[str, Path], m4_work: tuple[Path, Path], tmp_path: Path
) -> None:
    """Review #248 B2: count < n (exempt start/stop sites) publishes value = count, n = n."""
    qc = pd.read_parquet(strain_run["variants_qc"])
    mask = qc["metric"].astype(str) == "id_mixups"
    qc.loc[mask, "value"] = 1
    qc.loc[mask, "n"] = 3
    path = tmp_path / "variants_qc.parquet"
    qc.to_parquet(path, index=False)
    row = _id_mixups_row(_context_stats(m4_work, tmp_path / "cs", path))
    assert float(row["value"]) == 1.0
    assert int(row["n"]) == 3


def test_related_relation_is_still_counted(strain_run: dict[str, Path]) -> None:
    """Degenerate row "``genome_relation = related`` | still counted": the fixture runs under
    the default ``related`` (every site flagged ``strain_difference``), and the sites the
    oracle counts are exactly those, so the flag exempts none."""
    assert _cfg().genome_relation == "related"
    variants = _variants(strain_run["variants"])
    assert len(variants) and all("strain_difference" in f for f in variants["flags"])
    assert _oracle(variants)[0] > 0


def test_no_same_db_peptides_gives_zero_with_n_zero(
    other_run: dict[str, Path], m4_work: tuple[Path, Path], tmp_path: Path
) -> None:
    """Spec: "no ``search_db = same`` peptides | ``id_mixups`` 0 with ``n`` 0 (not NA: the check
    ran)"."""
    assert len(_variants(other_run["variants"])) > 0, "ALT sites exist, none same-db"
    row = _id_mixups_row(_context_stats(m4_work, tmp_path / "cs", other_run["variants_qc"]))
    assert row["value"] != "NA", row.to_dict()
    assert float(row["value"]) == 0.0
    assert int(row["n"]) == 0
    assert "requires_variants" not in _flag_set(row["flags"])


def test_without_variants_qc_the_m4_16_row_stays(
    m4_work: tuple[Path, Path], tmp_path: Path
) -> None:
    """Spec: "Without ``variants_qc`` (no peptides), the M4-16 behaviour stays"."""
    row = _id_mixups_row(_context_stats(m4_work, tmp_path / "cs", None))
    assert row["value"] == "NA"
    assert "requires_variants" in _flag_set(row["flags"])


def test_the_other_qc_rows_do_not_change_with_variants_qc(
    strain_run: dict[str, Path], m4_work: tuple[Path, Path], tmp_path: Path
) -> None:
    with_qc = _context_stats(m4_work, tmp_path / "a", strain_run["variants_qc"])
    without = _context_stats(m4_work, tmp_path / "b", None)
    keep = lambda d: d[d["metric"] != "id_mixups"].reset_index(drop=True)  # noqa: E731
    pd.testing.assert_frame_equal(keep(with_qc), keep(without))


# ------------------------------------------------------------------ example 1 (TE2 on sim)


def test_example_1_same_strain_sim_has_no_high_violation(
    same_chain: dict[str, Path], tmp_path: Path
) -> None:
    """Same-strain database peptides through ``run_variants``: zero ``high`` rows without
    ``start_codon`` (M7-12's ``same_strain_insilico`` version is in
    tests/integration/test_te2_proxy.py)."""
    import importlib

    repo = Path(__file__).resolve().parents[3]
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    te = importlib.import_module("validation.proxies.te_proxies")
    paths = _run_variants(same_chain, tmp_path / "out", _cfg())
    r = te.te2_result(_variants(paths["variants"]))
    assert r.counts["violations"] == 0
    assert r.passed is True
