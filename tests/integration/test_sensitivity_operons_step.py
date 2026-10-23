"""M6-10 at step level: the registered ``sensitivity`` step with the ``operon_flags`` family on the
SIM100K producer chain (docs/specs/M6-10-sensitivity-operon-flags.md Behaviour 1, 3-7, examples
6-7, degenerate cases, lead decisions 1-2; M4-19 decisions 2, 5, 16; GIAP_SPEC.md §8.6, G9).

Interface assumed (test author): the step declares the pair inputs (C2 ``genes``, ``pairs_raw``,
C3 ``features``; ``promoter_windows`` / ``terminators`` / ``tx_rates`` optional) and binds them to
``operon_flags_family`` by keyword (C2 ``genes`` -> ``gene_models``). For each variant it builds
the config as ``with_overrides(base, overrides | config_patch)`` with ``base`` at
``stats.permutations = 199``; a patch-only variant hands the family the primary gene table.
Expected rows: the merged M6-08 chain (``run_operons``) on each variant's gene table and config,
written with the M4-19 writer. Families other than ``operon_flags`` are removed from
``FAMILIES`` here so that only this family's rows are compared.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import sensitivity_operons_data as so
from helpers.layers_chain import Chain, build_chain

REPO = Path(__file__).resolve().parents[2]
CONTEXT_NF = REPO / "workflow" / "subworkflows" / "local" / "context.nf"
PAIR_INPUTS = {"genes", "pairs_raw", "features"}


def _S() -> Any:
    from giap.context import sensitivity as S

    return S


def _spec() -> Any:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    return reg.get("sensitivity")


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(Path(tmp_path_factory.mktemp("sens_operons_step")))


@pytest.fixture
def only_operons(monkeypatch: pytest.MonkeyPatch) -> None:
    S = _S()
    monkeypatch.setattr(S, "FAMILIES", {"operon_flags": S.FAMILIES["operon_flags"]})


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _run(chain: Chain, root: Path, gt: pd.DataFrame, cfg: Any, **kw: Any) -> Path:
    from giap.steps.registry import run_step

    _spec()
    ins = so.write_inputs(chain, root / "in", gt, **kw)
    out = root / "tables" / "sensitivity.tsv"
    run_step("sensitivity", cfg, ins, {"sensitivity": out})
    return out


def _expected(chain: Chain, gt_path: Path, cfg: Any, root: Path) -> pd.DataFrame:
    """M4-19 rows of the oracle sets for every variant the step reruns (top3 skipped without
    peptide intensities, as M4-19's step does), written and read back."""
    from giap.io.genes_tsv import read_genes_parquet
    from giap.io.sensitivity_tsv import write_sensitivity_tsv
    from giap.steps.sensitivity import with_overrides

    S = _S()
    gt = read_genes_parquet(gt_path)
    t = so.tables(chain)
    base = with_overrides(cfg, {"stats.permutations": so.SENS_B})
    primary = {so.FAMILY: so.oracle_flagged(gt, base, t)}
    reruns: dict[tuple[str, str], frozenset[str]] = {}
    for v in S.variants(cfg, str(gt["detection_tier"].iloc[0])):
        if v.name == "top3" and not gt["top3_log2"].notna().any():
            continue
        patch = dict(getattr(v, "config_patch", None) or {})
        vcfg = with_overrides(base, {**dict(v.overrides), **patch})
        g = gt if patch and not v.overrides else S.variant_genes(gt, v, cfg)
        reruns[(so.FAMILY, v.name)] = so.oracle_flagged(g, vcfg, t)
    out = root / "expected.tsv"
    write_sensitivity_tsv(out, S.sensitivity_rows(primary, reruns))
    return _read(out)


# --------------------------------------------------------------------------- wiring


def test_step_declares_the_pair_inputs() -> None:
    """Behaviour 7: the merged step gains the pair-related inputs."""
    spec = _spec()
    declared = set(spec.inputs) | set(spec.optional_inputs)
    assert declared >= PAIR_INPUTS
    assert set(so.LAYER_KEYS) <= set(spec.optional_inputs)
    assert {"gene_table", "replicons", "regions", "region_genes", "tus"} <= set(spec.inputs)
    assert spec.outputs == ("sensitivity",)


def test_context_nf_passes_the_pair_channels() -> None:
    """Behaviour 7: ``context.nf`` passes the pair channels to the ``sensitivity`` step (static
    check of its ``stepIn('sensitivity', ...)`` call; the run-time check is the e2e)."""
    text = re.sub(r"(?m)//.*$", "", CONTEXT_NF.read_text(encoding="utf-8"))
    m = re.search(r"stepIn\(\s*'sensitivity'\s*,(.*?)\)\s*\n\s*SENSITIVITY\(", text, re.S)
    assert m, "no stepIn('sensitivity', ...) followed by SENSITIVITY(...)"
    call = m.group(1)
    for key in ("pairs_raw", "features"):
        assert key in call, f"context.nf: sensitivity does not receive {key}"
    assert re.search(r"\bgenes\b|kv\('genes'", text), "context.nf: no C2 genes channel"


# --------------------------------------------------------------------------- rows


@pytest.mark.usefixtures("only_operons")
@pytest.mark.parametrize("level", ["basic", "full"])
def test_rows_match_the_chain_per_variant(chain: Chain, tmp_path: Path, level: str) -> None:
    """Behaviour 3 and 6: one row per variant, n / shared / J from the chain's sets, agreement
    NA; ``full`` has the two ``pair_threshold`` rows."""
    cfg = so.config(output={"sensitivity": level})
    gt = so.gene_table(so.tables(chain)["genes"], imputed=(1, 2, 3, 10, 11))
    out = _run(chain, tmp_path, gt, cfg)
    got = _read(out)
    assert set(got["family"]) == {so.FAMILY}
    want = _expected(chain, tmp_path / "in" / "genes_table.parquet", cfg, tmp_path)
    pd.testing.assert_frame_equal(got, want)
    names = list(got["variant"])
    if level == "full":
        assert names[-2:] == [so.variant_name(t) for t in so.UNIOP_THRESHOLDS]
        assert "top3" in names
    else:
        assert names == ["missing=exclude"]
    assert set(got["agreement"]) == {"NA"}
    assert got["n_primary"].nunique() == 1 and int(got["n_primary"].iloc[0]) > 0


@pytest.mark.usefixtures("only_operons")
def test_example_6_no_abundance_rows_all_na(chain: Chain, tmp_path: Path) -> None:
    """Example 6 / degenerate case 1 (sim100k, no abundance): every row has n 0 and J NA, no
    error; top3 is skipped (no peptide intensities)."""
    cfg = so.config(output={"sensitivity": "full"})
    gt = so.gene_table(so.tables(chain)["genes"], abundance=False)
    got = _read(_run(chain, tmp_path, gt, cfg))
    assert list(got["variant"]) == ["missing=exclude", "pair_threshold=0.3", "pair_threshold=0.7"]
    assert set(got["family"]) == {so.FAMILY}
    for c in ("n_primary", "n_variant", "n_shared"):
        assert set(got[c]) == {"0"}, c
    assert set(got["jaccard"]) == {"NA"}, "lead decision 2: J of two empty sets is NA"


@pytest.mark.usefixtures("only_operons")
def test_relative_scale_rows_na_with_reason(chain: Chain, tmp_path: Path,
                                            caplog: pytest.LogCaptureFixture) -> None:  # fmt: skip
    """Degenerate case 2: on a relative scale the row is written with NA J and the log says
    why (``relative``)."""
    cfg = so.config(output={"sensitivity": "full"})
    gt = so.gene_table(so.tables(chain)["genes"], scale="relative")
    with caplog.at_level(logging.DEBUG):
        got = _read(_run(chain, tmp_path, gt, cfg))
    assert len(got) and set(got["family"]) == {so.FAMILY}
    assert set(got["jaccard"]) == {"NA"}
    assert any("relative" in r.getMessage() for r in caplog.records
               if r.name.startswith("giap.")), "no log record names the relative scale"  # fmt: skip


@pytest.mark.usefixtures("only_operons")
def test_missing_layers_still_run(chain: Chain, tmp_path: Path) -> None:
    """The evidence layers are optional step inputs: without them the family still runs (as a
    missing layer) and the rows match the chain without layers."""
    cfg = so.config()
    gt = so.gene_table(so.tables(chain)["genes"])
    got = _read(_run(chain, tmp_path, gt, cfg, drop=so.LAYER_KEYS))
    assert list(got["variant"]) == ["missing=exclude"]
    from giap.io.genes_tsv import read_genes_parquet

    g = read_genes_parquet(tmp_path / "in" / "genes_table.parquet")
    t = {k: v for k, v in so.tables(chain).items() if k not in so.LAYER_KEYS}
    base = so.config(stats={"permutations": so.SENS_B})
    assert int(got["n_primary"].iloc[0]) == len(so.oracle_flagged(g, base, t))


@pytest.mark.usefixtures("only_operons")
def test_example_7_same_seed_identical_rows(chain: Chain, tmp_path: Path) -> None:
    cfg = so.config(output={"sensitivity": "full"})
    gt = so.gene_table(so.tables(chain)["genes"], imputed=(1, 2, 3))
    a = _run(chain, tmp_path / "a", gt, cfg).read_bytes()
    b = _run(chain, tmp_path / "b", gt, cfg).read_bytes()
    assert a == b


# --------------------------------------------------------------------------- patch at the step


class _Recorder:
    """A ``(genes, cfg, rng)`` family that records its calls and returns no flag."""

    def __init__(self) -> None:
        self.calls: list[tuple[pd.DataFrame, Any]] = []

    def __call__(self, genes: pd.DataFrame, cfg: Any, rng: np.random.Generator
                 ) -> frozenset[str]:  # fmt: skip
        self.calls.append((genes.copy(), cfg))
        return frozenset()


def _record(chain: Chain, tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch
            ) -> tuple[_Recorder, pd.DataFrame]:  # fmt: skip
    """Run the step function directly (every input given) with only the recorder registered."""
    from giap.io.genes_tsv import read_genes_parquet
    from giap.rng import for_step
    from giap.steps.registry import StepContext
    from giap.steps.sensitivity import sensitivity

    rec = _Recorder()
    monkeypatch.setattr(_S(), "FAMILIES", {"recorder": rec})
    gt = so.gene_table(so.tables(chain)["genes"], imputed=(1, 2, 3))
    ins = so.write_inputs(chain, tmp_path / "in", gt)
    out = tmp_path / "tables" / "sensitivity.tsv"
    sensitivity(StepContext(cfg=cfg, inputs=dict(ins), outputs={"sensitivity": out},
                            rng=for_step("sensitivity", cfg.seed),
                            log=logging.getLogger("giap.steps.sensitivity")))  # fmt: skip
    return rec, read_genes_parquet(ins["gene_table"])


def test_patch_variants_reach_the_family(chain: Chain, tmp_path: Path,
                                         monkeypatch: pytest.MonkeyPatch) -> None:  # fmt: skip
    """Lead decision 1: a ``pair_threshold`` variant gets the primary gene table and a config
    with the patched threshold, ``stats.permutations = 199`` and the primary abundance choices."""
    cfg = so.config(output={"sensitivity": "full"})
    rec, _ = _record(chain, tmp_path, cfg, monkeypatch)
    assert len(rec.calls) == 1 + len(_S().variants(cfg, "mayu"))
    primary_genes, primary_cfg = rec.calls[0]
    seen = []
    for genes, vcfg in rec.calls[1:]:
        t = vcfg.annotation.uniop.pair_threshold
        assert vcfg.stats.permutations == so.SENS_B
        if t != cfg.annotation.uniop.pair_threshold:
            seen.append(t)
            pd.testing.assert_frame_equal(genes, primary_genes)
            assert vcfg.abundance == cfg.abundance
            assert vcfg.model_dump() == so.with_threshold(primary_cfg, t).model_dump()
    assert sorted(seen) == list(so.UNIOP_THRESHOLDS)


def test_unpatched_variants_unchanged(chain: Chain, tmp_path: Path,
                                      monkeypatch: pytest.MonkeyPatch) -> None:  # fmt: skip
    """Lead decision 1, "existing variants are byte-identical" (a regression pin that holds
    before M6-10 too): a variant without patch gets ``variant_genes`` and
    ``with_overrides(base, overrides)``, with the primary threshold."""
    from giap.steps.sensitivity import with_overrides

    S = _S()
    cfg = so.config(output={"sensitivity": "full"})
    rec, gt = _record(chain, tmp_path, cfg, monkeypatch)
    base = with_overrides(cfg, {"stats.permutations": so.SENS_B})
    plain = [v for v in S.variants(cfg, "mayu") if not getattr(v, "config_patch", None)]
    plain = [v for v in plain if v.name != "top3" or gt["top3_log2"].notna().any()]
    assert [c for c in rec.calls[1 : 1 + len(plain)]]
    for v, (genes, vcfg) in zip(plain, rec.calls[1 : 1 + len(plain)], strict=True):
        assert vcfg.model_dump() == with_overrides(base, v.overrides).model_dump(), v.name
        pd.testing.assert_frame_equal(genes, S.variant_genes(gt, v, cfg), obj=v.name)
