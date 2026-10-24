"""M6-13: the ``operons.discordance_fit`` event of the ``operons`` step
(docs/specs/M6-13-discordance-abundance-sigma.md "Public interface", steps/operons.py; §6.7).

One INFO event per run, ``operons.discordance_fit``, whose ``data`` is strict JSON (NaN -> null)
with the keys below; a WARNING when the scaled fit fell back for ``not_converged``, ``collapsed``
or ``worse_than_constant``. ``gamma_per_log2 = gamma_per_sd / abundance_sd_log2`` (S-1). The
fields must equal ``Discordance.scaled`` of the same chain (``discordance`` on the step's own pair
frame). No table gains a column: the existing step tests pin the tables and schemas.

Inputs: the sim100k chain of ``test_operons_step`` (36 both-detected OP pairs with abundance,
none without).
"""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import scipy.optimize
from helpers.layers_chain import Chain
from test_operons_step import _cfg, chain, run_operons_step, step_inputs, write_gene_table

import giap.stats.tdist as tdist

__all__ = ["chain"]  # the module-scoped fixture, re-exported for pytest

EVENT = "operons.discordance_fit"
KEYS = {
    "sigma_model", "sigma_model_used", "reason", "n_fit", "loc_log2", "scale0_log2",
    "gamma_per_sd", "gamma_per_log2", "df", "df_at_bound", "gamma_at_bound",
    "abundance_center_log2", "abundance_sd_log2", "z_min", "z_max", "lr_stat", "lr_p",
}  # fmt: skip
_REAL_MINIMIZE = scipy.optimize.minimize


def _events(records: list[logging.LogRecord]) -> list[logging.LogRecord]:
    return [r for r in records if getattr(r, "event", None) == EVENT]


def _data(records: list[logging.LogRecord]) -> dict[str, Any]:
    found = _events(records)
    assert len(found) == 1, f"expected one {EVENT} event, got {len(found)}"
    assert found[0].levelno == logging.INFO
    data = found[0].data  # type: ignore[attr-defined]
    assert isinstance(data, dict)
    return data


def _scaled(inputs: dict[str, Path]) -> Any:
    """``Discordance.scaled`` of the step's own chain (as ``test_operons_step.oracle``)."""
    import pandas as pd

    from giap.annotate.layers import genes_of, replicons_of, terminators_of
    from giap.io.genes_tsv import read_genes_parquet
    from giap.operons.discordance import discordance
    from giap.operons.pairs import RNA_TYPES, Evidence, pair_frame

    rd = pd.read_parquet
    gt = read_genes_parquet(inputs["gene_table"])
    feats = rd(inputs["features"])
    ev = Evidence(
        windows=rd(inputs["promoter_windows"]),
        terminators=tuple(terminators_of(rd(inputs["terminators"]))),
        tx=rd(inputs["tx_rates"]),
        rna=feats[feats["type"].isin(RNA_TYPES)].reset_index(drop=True),
    )
    pairs = pair_frame(gt, genes_of(rd(inputs["genes"])), replicons_of(rd(inputs["replicons"])),
                       rd(inputs["pairs_raw"]), ev, _cfg())  # fmt: skip
    return discordance(pairs, gt, _cfg(), None).scaled


def _same(event: Any, value: float) -> None:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        assert event is None
    else:
        assert event == pytest.approx(float(value), rel=1e-12, abs=1e-15)


@pytest.fixture(scope="module")
def with_abundance(chain: Chain, tmp_path_factory: pytest.TempPathFactory) -> Any:
    root = tmp_path_factory.mktemp("operons_sigma")
    gt = write_gene_table(chain, root / "genetable" / "genes.parquet", abundance=True)
    ins = step_inputs(chain, gt)
    return run_operons_step(ins, root / "out", _cfg()), ins


def test_one_event_with_the_spec_keys(with_abundance: Any) -> None:
    run, _ = with_abundance
    data = _data(run.records)
    assert set(data) >= KEYS
    json.dumps(data, allow_nan=False)  # strict JSON: no NaN, no inf, no numpy scalars
    assert data["sigma_model"] == "loglinear"
    assert data["sigma_model_used"] in {"loglinear", "constant", "none"}


def test_event_equals_the_scaled_fit_of_the_chain(with_abundance: Any) -> None:
    run, ins = with_abundance
    data = _data(run.records)
    s = _scaled(ins)
    assert data["sigma_model_used"] == s.model
    assert data["reason"] == s.reason
    assert data["n_fit"] == s.n
    for key, value in (("loc_log2", s.loc), ("scale0_log2", s.scale0), ("gamma_per_sd", s.gamma),
                       ("df", s.df), ("abundance_center_log2", s.center),
                       ("abundance_sd_log2", s.sd), ("z_min", s.z_min), ("z_max", s.z_max),
                       ("lr_stat", s.lr_stat), ("lr_p", s.lr_p)):  # fmt: skip
        _same(data[key], value)


def test_gamma_per_log2_is_the_physical_slope(with_abundance: Any) -> None:
    run, _ = with_abundance
    data = _data(run.records)
    if data["sigma_model_used"] == "none":
        assert data["gamma_per_log2"] is None
        return
    assert data["n_fit"] >= 20
    assert data["gamma_per_log2"] == pytest.approx(
        data["gamma_per_sd"] / data["abundance_sd_log2"], rel=1e-12)  # fmt: skip


def test_bound_flags_are_booleans(with_abundance: Any) -> None:
    run, _ = with_abundance
    data = _data(run.records)
    assert isinstance(data["df_at_bound"], bool) and isinstance(data["gamma_at_bound"], bool)
    if data["sigma_model_used"] == "loglinear":
        assert data["gamma_at_bound"] == (abs(data["gamma_per_sd"]) == 2.0)
        assert data["df_at_bound"] == (data["df"] in (2.0, 50.0))


def test_chain_with_abundance_fits_the_loglinear_null(with_abundance: Any) -> None:
    """36 both-detected OP pairs with varied abundance (12 + sin(i)): the scaled fit runs."""
    run, _ = with_abundance
    data = _data(run.records)
    assert data["sigma_model_used"] == "loglinear" and data["reason"] is None
    assert not [r for r in run.records if r.levelno >= logging.WARNING
                and getattr(r, "event", None) == EVENT]  # fmt: skip


def test_no_abundance_event_is_none(chain: Chain, tmp_path: Path) -> None:
    gt = write_gene_table(chain, tmp_path / "genetable" / "genes.parquet", abundance=False)
    run = run_operons_step(step_inputs(chain, gt), tmp_path / "out", _cfg())
    data = _data(run.records)
    assert data["sigma_model_used"] == "none"
    assert data["n_fit"] == 0
    json.dumps(data, allow_nan=False)
    assert data["gamma_per_sd"] is None and data["loc_log2"] is None


def test_fallback_logs_a_warning(chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                 caplog: pytest.LogCaptureFixture) -> None:  # fmt: skip
    """A scaled fit that does not converge: constant null, ``reason = not_converged``, WARNING."""

    def fake(fun: Any, x0: Any, *a: Any, **k: Any) -> Any:
        res = _REAL_MINIMIZE(fun, x0, *a, **k)
        if np.asarray(x0).size == 4:
            res = scipy.optimize.OptimizeResult(dict(res))
            res.update(success=False, status=2)
        return res

    monkeypatch.setattr(scipy.optimize, "minimize", fake)
    if hasattr(tdist, "minimize"):
        monkeypatch.setattr(tdist, "minimize", fake)
    gt = write_gene_table(chain, tmp_path / "genetable" / "genes.parquet", abundance=True)
    with caplog.at_level(logging.WARNING, logger="giap"):
        run = run_operons_step(step_inputs(chain, gt), tmp_path / "out", _cfg())
    data = _data(run.records)
    assert data["sigma_model_used"] == "constant" and data["reason"] == "not_converged"
    assert data["gamma_per_sd"] == 0.0
    warned = [r for r in [*run.records, *caplog.records]
              if r.levelno >= logging.WARNING and "not_converged" in r.getMessage()]  # fmt: skip
    assert warned
