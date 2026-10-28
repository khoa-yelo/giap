"""M5-10 acceptance "Integration": the ``sensitivity`` step with the ``islands`` family on a
gippasim
work directory (docs/specs/M5-10-sensitivity-islands.md examples 4-5, edge case 1, lead decision
38; M4-19 lead decisions 2, 5, 12; GIPPA_SPEC.md §8.6, G9).

Work directory: the M5-04 ``island_power`` scenario (1,000 genes, AR(1) rho = 0.5, one planted
15-TU island at +4 SD), 60 genes stripped of peptide intensities, the M4-15 gene table
(``tests/helpers/sensitivity_islands_data.sim_work_dir``). Expected sets come from the merged
M5-06 producers at B = 199 (``stats.bootstrap`` and ``stats.permutations``) on the family stream
``for_step("sensitivity/islands", seed)`` (M4-19 decision 2: common random numbers, one stream
per family), each variant's genes from ``variant_genes``.

Interface assumed: ``run_sensitivity(work_dir, cfg)`` binds C5 ``tus`` (``tus.parquet``) and C1
``replicons`` to ``islands_family``; the registered step ``sensitivity`` declares ``tus`` as an
input (it is required: islands cannot run without C5, M4-19 decision 12).
"""

from __future__ import annotations

import json
import logging
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import sensitivity_islands_data as sid


def _run(work: Path, cfg: Any) -> Path:
    from gippa.steps.sensitivity import run_sensitivity

    run_sensitivity(work, cfg)
    return work / "tables" / "sensitivity.tsv"


def _copy(src: Path, dst: Path) -> Path:
    shutil.copytree(src, dst)
    (dst / "tables" / "sensitivity.tsv").unlink(missing_ok=True)
    return dst


@pytest.fixture(scope="module")
def work(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return sid.sim_work_dir(tmp_path_factory.mktemp("sens_isl") / "work")


@pytest.fixture(scope="module")
def work_observed(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Every detected gene keeps its peptide intensities (nothing imputed)."""
    return sid.sim_work_dir(tmp_path_factory.mktemp("sens_isl_obs") / "work",
                            blank_intensity_genes=0)  # fmt: skip


@pytest.fixture(scope="module")
def expected(work: Path) -> tuple[frozenset[str], dict[str, frozenset[str]]]:
    return sid.oracle_sets(work, sid.work_config())


def _islands(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["family"] == "islands"].set_index("variant")


def _schema() -> dict[str, Any]:
    import gippa

    path = Path(gippa.__file__).resolve().parent / "schemas" / "tables" / "sensitivity.schema.json"
    return json.loads(path.read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def _validate(tsv: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    res = frictionless.Resource(path=tsv.name, basepath=str(tsv.parent),
                                schema=frictionless.Schema.from_descriptor(_schema()))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


# --------------------------------------------------------------------------- the table


def test_fixture_has_an_island_and_imputed_genes(
    work: Path, expected: tuple[frozenset[str], dict[str, frozenset[str]]]
) -> None:
    g = sid.sim_genes(work)
    assert g["abundance_imputed"].astype(bool).any(), "fixture: nothing imputed"
    assert len(expected[0]) >= 1, "fixture: no island at B = 199"
    assert set(expected[1]) == {"missing=min", "missing=exclude"}


def test_writes_schema_valid_islands_rows(work: Path) -> None:
    tsv = _run(work, sid.work_config())
    _validate(tsv)
    df = sid.read_sensitivity(tsv)
    pairs = set(zip(df["family"], df["variant"], strict=True))
    assert {("islands", "missing=min"), ("islands", "missing=exclude")} <= pairs
    assert {("regions", "missing=min"), ("regions", "missing=exclude")} <= pairs
    isl = df[df["family"] == "islands"]
    assert len(isl) == 2
    assert (isl["agreement"] == "NA").all()


def test_islands_rows_match_the_m5_06_sets_at_b_199(
    work: Path, expected: tuple[frozenset[str], dict[str, frozenset[str]]]
) -> None:
    """Primary recomputed at B = 199 and every variant at B = 199, all on the one family
    stream; counts and Jaccard as M4-19 writes them (6 significant digits)."""
    primary, variants = expected
    df = _islands(sid.read_sensitivity(_run(work, sid.work_config())))
    for name, v in variants.items():
        r = df.loc[name]
        assert int(r["n_primary"]) == len(primary), name
        assert int(r["n_variant"]) == len(v), name
        assert int(r["n_shared"]) == len(primary & v), name
        u = primary | v
        if u:
            assert float(r["jaccard"]) == pytest.approx(len(primary & v) / len(u), abs=1e-6)
        else:
            assert r["jaccard"] == "NA"


def test_islands_run_at_b_199_for_primary_and_every_variant(
    work: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Spec: reruns use ``stats.bootstrap = 199`` and ``stats.permutations = 199`` for the
    primary and each variant (the user's ``stats.bootstrap`` is 50 here)."""
    w = _copy(work, tmp_path / "w")
    cfgs = sid.record_scans(monkeypatch)
    cfg = sid.work_config()
    assert cfg.stats.bootstrap != sid.SENS_B
    _run(w, cfg)
    assert len(cfgs) == 3  # primary, missing=min, missing=exclude
    assert {(c.stats.bootstrap, c.stats.permutations) for c in cfgs} == {(199, 199)}
    assert all(c.context.islands == cfg.context.islands for c in cfgs)


def test_example_5_exclude_row_is_produced(
    work: Path, expected: tuple[frozenset[str], dict[str, frozenset[str]]]
) -> None:
    """Example 5: under ``exclude`` the detected genes without abundance are unscored; the row
    is still written, with the M5-06 set of the ``exclude`` gene table."""
    df = _islands(sid.read_sensitivity(_run(work, sid.work_config())))
    assert "missing=exclude" in df.index
    assert int(df.loc["missing=exclude", "n_variant"]) == len(expected[1]["missing=exclude"])


def test_exclude_with_nothing_imputed_is_identical(work_observed: Path) -> None:
    """Edge case 1: every detected gene observed -> ``exclude`` = primary (jaccard 1 or NA)."""
    g = sid.sim_genes(work_observed)
    assert not g["abundance_imputed"].astype(bool).any(), "fixture: a gene is imputed"
    df = _islands(sid.read_sensitivity(_run(work_observed, sid.work_config())))
    r = df.loc["missing=exclude"]
    assert r["n_primary"] == r["n_variant"] == r["n_shared"]
    assert r["jaccard"] in ("1", "NA")


def test_same_seed_byte_identical(work: Path) -> None:
    cfg = sid.work_config()
    a = _run(work, cfg).read_bytes()
    assert _run(work, cfg).read_bytes() == a


def test_islands_rows_do_not_depend_on_other_families(
    work: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """M4-19 decision 2 / §6.6: the islands family has its own stream, so running it alone
    gives the same rows as running it after ``regions``."""
    from gippa.context import sensitivity as S

    full = _islands(sid.read_sensitivity(_run(work, sid.work_config())))
    w = _copy(work, tmp_path / "w")
    monkeypatch.setattr(S, "FAMILIES", {"islands": S.FAMILIES["islands"]})
    alone = sid.read_sensitivity(_run(w, sid.work_config()))
    assert set(alone["family"]) == {"islands"}
    pd.testing.assert_frame_equal(_islands(alone), full)


def test_islands_rows_depend_on_the_seed_stream(
    work: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The family's generator is ``for_step("sensitivity/islands", seed)``, in the same state for
    the primary and each variant (common random numbers)."""
    from gippa.context import sensitivity as S
    from gippa.rng import for_step

    states: list[Any] = []
    real = S.FAMILIES["islands"]

    def family(genes: pd.DataFrame, cfg: Any, rng: np.random.Generator, *,
               tus: pd.DataFrame, replicons: pd.DataFrame) -> frozenset[str]:  # fmt: skip
        states.append(rng.bit_generator.state)
        return real(genes, cfg, rng, tus=tus, replicons=replicons)  # type: ignore[no-any-return]

    w = _copy(work, tmp_path / "w")
    monkeypatch.setattr(S, "FAMILIES", {"islands": family})
    cfg = sid.work_config()
    _run(w, cfg)
    want = for_step(sid.STREAM, cfg.seed).bit_generator.state
    assert len(states) == 3
    assert all(s == want for s in states)


def test_step_binds_c5_tus_and_c1_replicons(
    work: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Decision 38: the family gets the C5 ``tus`` and C1 ``replicons`` files' tables."""
    from gippa.context import sensitivity as S

    seen: list[tuple[pd.DataFrame, pd.DataFrame]] = []

    def family(genes: pd.DataFrame, cfg: Any, rng: np.random.Generator, *,
               tus: pd.DataFrame, replicons: pd.DataFrame) -> frozenset[str]:  # fmt: skip
        seen.append((tus, replicons))
        return frozenset()

    w = _copy(work, tmp_path / "w")
    monkeypatch.setattr(S, "FAMILIES", {"islands": family})
    _run(w, sid.work_config())
    m4 = sid.sim_m4(w)
    assert len(seen) == 3
    for tus, reps in seen:
        assert sorted(tus["tu_id"].astype(str)) == sorted(m4.tus["tu_id"].astype(str))
        assert list(reps["contig"].astype(str)) == list(m4.replicons["contig"].astype(str))
        assert list(reps["topology"].astype(str)) == list(m4.replicons["topology"].astype(str))


def test_missing_tus_raises(work: Path, tmp_path: Path) -> None:
    """M4-19 decision 12 applied to C5: a missing required table is an error, never rows."""
    w = _copy(work, tmp_path / "w")
    (w / "tus.parquet").unlink()
    with pytest.raises((ValueError, OSError)):
        _run(w, sid.work_config())
    assert not (w / "tables" / "sensitivity.tsv").exists()


# --------------------------------------------------------------------------- off (example 4)


def test_example_4_off_writes_no_file(work: Path, tmp_path: Path) -> None:
    w = _copy(work, tmp_path / "w")
    from gippa.steps.sensitivity import run_sensitivity

    assert run_sensitivity(w, sid.work_config(output={"sensitivity": "off"})) is None
    assert not (w / "tables" / "sensitivity.tsv").exists()


# --------------------------------------------------------------------------- registered step


def _ctx(work: Path, cfg: Any, out: Path, *, tus: bool = True) -> Any:
    from gippa.rng import for_step
    from gippa.steps.gene_table import contract_paths
    from gippa.steps.registry import StepContext

    paths = contract_paths(work)
    keys = ["replicons", "regions", "region_genes"] + (["tus"] if tus else [])
    inputs: dict[str, Any] = {k: paths[k] for k in keys}
    inputs["gene_table"] = work / "genetable" / "genes.parquet"
    return StepContext(cfg=cfg, inputs=inputs, outputs={"sensitivity": out},
                       rng=for_step("sensitivity", cfg.seed),
                       log=logging.getLogger("test.sensitivity"))  # fmt: skip


def test_registered_step_declares_tus() -> None:
    from gippa.steps import load_step_modules
    from gippa.steps.registry import get

    load_step_modules()
    spec = get("sensitivity")
    assert "tus" in spec.inputs, "islands need C5: tus is a required input (decision 38)"
    assert {"gene_table", "replicons", "regions", "region_genes"} <= set(spec.inputs)
    assert spec.outputs == ("sensitivity",)


def test_registered_step_writes_the_same_table(work: Path, tmp_path: Path) -> None:
    from gippa.steps.sensitivity import sensitivity

    cfg = sid.work_config()
    want = _run(work, cfg).read_bytes()
    out = tmp_path / "tables" / "sensitivity.tsv"
    sensitivity(_ctx(work, cfg, out))
    assert out.read_bytes() == want
    assert "islands" in set(sid.read_sensitivity(out)["family"])


def test_registered_step_off_has_no_islands_rows(work: Path, tmp_path: Path) -> None:
    from gippa.steps.sensitivity import sensitivity

    out = tmp_path / "tables" / "sensitivity.tsv"
    sensitivity(_ctx(work, sid.work_config(output={"sensitivity": "off"}), out))
    assert len(sid.read_sensitivity(out)) == 0
