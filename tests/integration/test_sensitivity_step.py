"""M4-19 ``run_sensitivity`` on a giapsim work directory (docs/specs/M4-19-sensitivity-regions.md
acceptance "Integration", examples 6-7, edge cases; lead decisions 2-4; GIAP_SPEC.md §8.6, G9).

Work directory: M4-07/08/09 contract tables of a 1,000-gene layout (prophage detection shifted,
so ``prophage_1`` O - E is flagged), 60 genes stripped of peptide intensities (imputed under the
primary ``median``), and the M4-15 gene table at ``genetable/genes.parquet``
(``tests/helpers/sensitivity_data.sim_work_dir``). Output: ``tables/sensitivity.tsv``, validated
against ``src/giap/schemas/tables/sensitivity(.schema).json``.

Common random numbers (lead decision 2) are checked with a recording family registered in
``FAMILIES``: every call of one family (the primary recomputation and each variant) gets a
generator in the same state, different families different states, and all run with
``stats.permutations = 199``.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import sensitivity_data as sd


def _run(work: Path, cfg: Any) -> None:
    from giap.steps.sensitivity import run_sensitivity

    run_sensitivity(work, cfg)


@pytest.fixture(scope="module")
def work(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return sd.sim_work_dir(tmp_path_factory.mktemp("sens_default") / "work")


@pytest.fixture(scope="module")
def work_2pep(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """No q-values, PEPs or decoys: the ``two_peptide`` tier (single-peptide genes present)."""
    return sd.sim_work_dir(tmp_path_factory.mktemp("sens_2pep") / "work",
                           evidence={"emit": frozenset({"intensity"})})  # fmt: skip


@pytest.fixture(scope="module")
def work_protein(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """No peptide intensities; abundance from the C10 protein table (``top3`` impossible)."""
    return sd.sim_work_dir(
        tmp_path_factory.mktemp("sens_protein") / "work",
        evidence={"emit": frozenset({"q_value", "pep", "is_decoy", "protein_ids"}),
                  "abundance_table": "ibaq"},
        blank_intensity_genes=0,
        cfg_over={"abundance": {"scale": "linear"}},
    )  # fmt: skip


def _fixture_ok(work: Path) -> pd.DataFrame:
    g = sd.sim_genes(work)
    assert g["abundance_imputed"].astype(bool).any(), "fixture: no imputed gene"
    return g


def _schema() -> dict[str, Any]:
    import giap

    table_dir = Path(giap.__file__).resolve().parent / "schemas" / "tables"
    found = [p for p in (table_dir / "sensitivity.schema.json", table_dir / "sensitivity.json")
             if p.is_file()]  # fmt: skip
    assert len(found) == 1, f"expected one sensitivity schema in {table_dir}, got {found}"
    return json.loads(found[0].read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def _validate(tsv: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in _schema().items() if k != "coordinates"}
    res = frictionless.Resource(
        path=tsv.name, basepath=str(tsv.parent), schema=frictionless.Schema.from_descriptor(d)
    )
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def _reported_n(work: Path, cfg: Any) -> int:
    """Lead decision 4 on the reported M4-17 tests (B = 1,000), what regions_stats.tsv holds."""
    from giap.context.regions import test_all_regions
    from giap.rng import for_step

    m4 = sd.sim_m4(work)
    regions = m4.regions[m4.regions["region_type"].isin(cfg.context.region_types)]
    res = test_all_regions(sd.sim_genes(work), regions, m4.region_genes, m4.replicons, cfg,
                           for_step("context_stats", cfg.seed))  # fmt: skip
    return len(sd.flagged_from(res, cfg.stats.alpha.context))


# --------------------------------------------------------------------------- the table


def test_writes_schema_valid_regions_rows(work: Path) -> None:
    _fixture_ok(work)
    cfg = sd.config()
    _run(work, cfg)
    tsv = sd.sensitivity_tsv(work)
    assert tsv.is_file()
    header = tsv.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert tuple(header[:6]) == sd.COLUMNS
    assert header[6:] in ([], ["agreement"])
    _validate(tsv)
    df = sd.read_sensitivity(tsv)
    assert sorted(zip(df["family"], df["variant"], strict=True)) == [
        ("regions", "missing=exclude"), ("regions", "missing=min"),
    ]  # fmt: skip
    assert str(work) not in tsv.read_text(encoding="utf-8")


def test_n_primary_equals_the_reported_count(work: Path) -> None:
    """Lead decision 2: below 20,000 genes the primary recomputed at B = 199 is the reported one,
    so n_primary is the regions_stats.tsv count of q <= alpha.context shift tests (>= 1 here)."""
    cfg = sd.config()
    _run(work, cfg)
    df = sd.read_sensitivity(sd.sensitivity_tsv(work))
    want = _reported_n(work, cfg)
    assert want >= 1, "fixture: no flagged region"
    assert set(df["n_primary"].astype(int)) == {want}


def test_example_7_exclude_keeps_every_o_minus_e_id(
    work: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Example 7 in the step: the ``missing=exclude`` genes keep the primary detection calls and
    p_detect, so the O - E IDs equal the primary's."""
    from giap.context import sensitivity as S

    rec = sd.Recorder()
    real = S.regions_family

    def family(genes: pd.DataFrame, cfg: Any, rng: np.random.Generator) -> frozenset[str]:
        rec(genes, cfg, rng)
        return real(genes, cfg, rng)

    monkeypatch.setattr(S, "FAMILIES", {"regions": family})
    _run(work, sd.config())
    prim = sd.sim_genes(work)
    excl = [c for c in rec.calls if not c.genes["abundance_imputed"].astype(bool).any()]
    assert len(excl) == 1
    g = excl[0].genes.set_index("locus_tag").loc[prim["locus_tag"]]
    assert (g["detected"].astype(bool).to_numpy() == prim["detected"].astype(bool).to_numpy()).all()
    assert np.allclose(g["p_detect"].to_numpy(float), prim["p_detect"].to_numpy(float),
                       equal_nan=True)  # fmt: skip
    imp = prim["abundance_imputed"].astype(bool).to_numpy()
    assert np.isnan(g["abundance_log2"].to_numpy(float)[imp]).all()


def test_same_seed_byte_identical(work: Path) -> None:
    cfg = sd.config()
    _run(work, cfg)
    a = sd.sensitivity_tsv(work).read_bytes()
    _run(work, cfg)
    assert sd.sensitivity_tsv(work).read_bytes() == a


def test_empty_sets_written_na(work: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Example 2 in the file: jaccard(empty, empty) is written ``NA``, counts 0."""
    from giap.context import sensitivity as S

    monkeypatch.setattr(S, "FAMILIES", {"regions": sd.Recorder(frozenset())})
    _run(work, sd.config())
    df = sd.read_sensitivity(sd.sensitivity_tsv(work))
    assert len(df) == 2
    assert (df["jaccard"] == "NA").all()
    assert (df[["n_primary", "n_variant", "n_shared"]] == "0").all(axis=None)


def test_example_6_off_writes_no_file(work: Path) -> None:
    tsv = sd.sensitivity_tsv(work)
    if tsv.exists():
        tsv.unlink()
    _run(work, sd.config(output={"sensitivity": "off"}))
    assert not tsv.exists()


# --------------------------------------------------------------------------- CRN, B = 199


def _record(monkeypatch: pytest.MonkeyPatch, **families: sd.Recorder) -> None:
    from giap.context import sensitivity as S

    monkeypatch.setattr(S, "FAMILIES", dict(families))


def test_primary_recomputed_and_every_variant_at_b_199(
    work: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rec = sd.Recorder(frozenset({"x:chr1:o_minus_e"}))
    _record(monkeypatch, regions=rec)
    cfg = sd.config()
    _run(work, cfg)
    assert len(rec.calls) == 3  # primary + missing=min + missing=exclude
    assert all(c.cfg.stats.permutations == sd.SENS_B for c in rec.calls)
    assert all(c.cfg.seed == cfg.seed for c in rec.calls)
    prim = sd.sim_genes(work)
    same = [c for c in rec.calls
            if np.allclose(c.genes.set_index("locus_tag").loc[prim["locus_tag"],
                                                               "abundance_log2"].to_numpy(float),
                           prim["abundance_log2"].to_numpy(float), equal_nan=True)]  # fmt: skip
    assert same, "no call on the reported (primary) gene table"


def test_common_random_numbers_one_stream_per_family(
    work: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Lead decision 2: one child stream keyed by family (not family x variant): every call of a
    family starts from the same generator state; two families get different states."""
    a, b = sd.Recorder(), sd.Recorder()
    _record(monkeypatch, regions=a, other=b)
    _run(work, sd.config())
    for rec in (a, b):
        assert len(rec.calls) == 3
        assert len({c.first_draw for c in rec.calls}) == 1, [c.first_draw for c in rec.calls]
        assert all(c.state == rec.calls[0].state for c in rec.calls)
    assert a.calls[0].first_draw != b.calls[0].first_draw


def test_family_stream_depends_on_seed_only(work: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The family stream does not depend on which variants run (``for_step("sensitivity", seed)``
    child keyed by family name) but does depend on the seed."""
    a, b, c = sd.Recorder(), sd.Recorder(), sd.Recorder()
    _record(monkeypatch, regions=a)
    _run(work, sd.config())
    _record(monkeypatch, regions=b)
    _run(work, sd.config(abundance={"sensitivity": ["exclude"]}))
    _record(monkeypatch, regions=c)
    _run(work, sd.config(seed=sd.SIM_SEED + 1))
    assert len(b.calls) == 2
    assert a.calls[0].first_draw == b.calls[0].first_draw
    assert a.calls[0].first_draw != c.calls[0].first_draw


def test_variant_genes_differ_only_in_abundance(
    work: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``missing=min`` / ``exclude`` rebuild the abundance only (M4-14 with ``policy``): detection
    and p_detect are the primary's; under ``min`` the imputed genes get the minimum observed
    log2 abundance of the detected non-partial genes."""
    rec = sd.Recorder()
    _record(monkeypatch, regions=rec)
    _run(work, sd.config())
    prim = _fixture_ok(work)
    imp = prim["abundance_imputed"].astype(bool).to_numpy()
    obs = prim.loc[prim["detected"].astype(bool) & ~prim["partial"].astype(bool) & ~imp,
                   "abundance_log2"].astype(float)  # fmt: skip
    seen_min = False
    for c in rec.calls:
        g = c.genes.set_index("locus_tag").loc[prim["locus_tag"]]
        assert (g["detected"].astype(bool).to_numpy() == prim["detected"].astype(bool)).all()
        assert np.allclose(g["p_detect"].to_numpy(float), prim["p_detect"].to_numpy(float),
                           equal_nan=True)  # fmt: skip
        vals = g["abundance_log2"].to_numpy(float)[imp]
        if np.allclose(vals, obs.min()):
            seen_min = True
    assert seen_min, "no missing=min call with the minimum observed abundance imputed"


# --------------------------------------------------------------------------- two_peptide tier


def test_two_peptide_adds_single_peptide_row(work_2pep: Path) -> None:
    g = sd.sim_genes(work_2pep)
    assert set(g["detection_tier"].astype(str)) == {"two_peptide"}
    assert (g["detection_status"] == "single_peptide_unestimated").any()
    _run(work_2pep, sd.config())
    df = sd.read_sensitivity(sd.sensitivity_tsv(work_2pep))
    assert set(df["variant"]) == {"missing=min", "missing=exclude", "single_peptide_detected"}
    _validate(sd.sensitivity_tsv(work_2pep))


def test_single_peptide_variant_detects_and_refits(
    work_2pep: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """§7.3.2: single-peptide genes count as detected in the rerun, and ``p_detect`` is refitted
    (M4-13) on the new detection calls; the family still gets the shared stream."""
    rec = sd.Recorder()
    _record(monkeypatch, regions=rec)
    _run(work_2pep, sd.config())
    prim = sd.sim_genes(work_2pep)
    single = (prim["detection_status"] == "single_peptide_unestimated").to_numpy()
    want = prim["detected"].astype(bool).to_numpy() | single
    sp = [c for c in rec.calls
          if c.genes.set_index("locus_tag").loc[prim["locus_tag"], "detected"].astype(bool)
          .to_numpy().sum() == want.sum()]  # fmt: skip
    assert len(sp) == 1
    g = sp[0].genes.set_index("locus_tag").loc[prim["locus_tag"]]
    assert (g["detected"].astype(bool).to_numpy() == want).all()
    p_new, p_old = g["p_detect"].to_numpy(float), prim["p_detect"].to_numpy(float)
    assert (np.isnan(p_new) == np.isnan(p_old)).all()
    assert not np.allclose(p_new, p_old, equal_nan=True)
    assert len({c.first_draw for c in rec.calls}) == 1


def test_single_peptide_refit_uses_common_folds(
    work_2pep: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Lead decision 2: the M4-13 folds of the refit are common random numbers too. Reading
    (see report): the refit's generator is the one the primary fit used,
    ``for_step("detectability", seed)``, so the TU folds are the primary's and only the
    detection calls differ; a per-variant fold stream fails this."""
    import giap.genetable.detectability as D
    from giap.rng import for_step

    real = D.fit_detectability
    states: list[dict[str, Any]] = []

    def spy(*args: Any, **kwargs: Any) -> Any:
        rng = kwargs.get("rng", args[4] if len(args) > 4 else None)
        states.append(dict(rng.bit_generator.state))
        return real(*args, **kwargs)

    monkeypatch.setattr(D, "fit_detectability", spy)
    import giap.steps.sensitivity as SS

    monkeypatch.setattr(SS, "fit_detectability", spy, raising=False)
    cfg = sd.config()
    _run(work_2pep, cfg)
    assert len(states) == 1
    assert states[0] == for_step("detectability", cfg.seed).bit_generator.state


# --------------------------------------------------------------------------- full / top3


def test_full_adds_top3_with_peptide_intensities(work: Path) -> None:
    _run(work, sd.config(output={"sensitivity": "full"}))
    df = sd.read_sensitivity(sd.sensitivity_tsv(work))
    assert set(df["variant"]) == {"missing=min", "missing=exclude", "top3"}
    _validate(sd.sensitivity_tsv(work))


def test_top3_skipped_without_peptide_intensities(
    work_protein: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Edge case: ``top3`` without peptide intensities is skipped with a log line."""
    cfg = sd.config(output={"sensitivity": "full"}, abundance={"scale": "linear"})
    with caplog.at_level(logging.DEBUG):
        _run(work_protein, cfg)
    df = sd.read_sensitivity(sd.sensitivity_tsv(work_protein))
    assert "top3" not in set(df["variant"])
    assert {"missing=min", "missing=exclude"} <= set(df["variant"])
    assert any("top3" in r.getMessage() for r in caplog.records)
