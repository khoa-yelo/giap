"""M6-13: the ``operons.discordance_fit`` event of the ``operons`` step
(docs/specs/M6-13-discordance-abundance-sigma.md "Public interface", steps/operons.py; §6.7).

One event per run, ``operons.discordance_fit``, whose ``data`` is strict JSON (NaN -> null)
with the keys below. It is INFO, or WARNING when the scaled fit fell back for ``not_converged``,
``collapsed`` or ``worse_than_constant`` (review CR S-1: the WARNING-level event is what reaches
``logs/events.jsonl`` and the manifest warnings). ``gamma_per_log2 = gamma_per_sd /
abundance_sd_log2`` (S-1). The fields must equal ``Discordance.scaled`` of the same chain
(``discordance`` on the step's own pair frame). No table gains a column: the existing step tests
pin the tables and schemas.

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

import gippa.stats.tdist as tdist

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


def _data(records: list[logging.LogRecord], level: int = logging.INFO) -> dict[str, Any]:
    found = _events(records)
    assert len(found) == 1, f"expected one {EVENT} event, got {len(found)}"
    assert found[0].levelno == level, logging.getLevelName(found[0].levelno)
    data = found[0].data  # type: ignore[attr-defined]
    assert isinstance(data, dict)
    return data


def _scaled(inputs: dict[str, Path]) -> Any:
    """``Discordance.scaled`` of the step's own chain (as ``test_operons_step.oracle``)."""
    import pandas as pd

    from gippa.annotate.layers import genes_of, replicons_of, terminators_of
    from gippa.io.genes_tsv import read_genes_parquet
    from gippa.operons.discordance import discordance
    from gippa.operons.pairs import RNA_TYPES, Evidence, pair_frame

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
    if data["sigma_model_used"] != "none":
        assert data["df_at_bound"] == (data["df"] in tdist.DF_BOUNDS)
    if data["sigma_model_used"] == "loglinear":
        assert data["gamma_at_bound"] == (data["gamma_per_sd"] in tdist.GAMMA_BOUNDS)


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


def _patch_scaled(monkeypatch: pytest.MonkeyPatch, edit: Any) -> None:
    """Patch the four-parameter optimiser only; the constant fit runs unchanged."""

    def fake(fun: Any, x0: Any, *a: Any, **k: Any) -> Any:
        res = _REAL_MINIMIZE(fun, x0, *a, **k)
        if np.asarray(x0).size != 4:
            return res
        out = scipy.optimize.OptimizeResult(dict(res))
        out.update(edit(res, np.asarray(x0, dtype=float)))
        return out

    monkeypatch.setattr(scipy.optimize, "minimize", fake)
    if hasattr(tdist, "minimize"):
        monkeypatch.setattr(tdist, "minimize", fake)


def _not_converged(res: Any, x0: np.ndarray) -> dict[str, Any]:
    return {"success": False, "status": 2}


def _collapsed(res: Any, x0: np.ndarray) -> dict[str, Any]:
    x = x0.copy()
    x[1] -= 40.0  # sigma0 at e^-40 of the constant scale: under the collapse floor
    return {"x": x, "success": True}


def _worse(res: Any, x0: np.ndarray) -> dict[str, Any]:
    x = x0.copy()
    x[0] += 5.0 * math.exp(x0[1])  # mu five scales off: fits worse than the constant start
    return {"x": x, "success": True}


FALLBACKS = {"not_converged": _not_converged, "collapsed": _collapsed,
             "worse_than_constant": _worse}  # fmt: skip


@pytest.mark.parametrize("reason", sorted(FALLBACKS))
def test_fallback_logs_a_warning(chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                 reason: str) -> None:  # fmt: skip
    """CR S-1, T-2: a scaled fit that falls back (not converged, collapsed, worse than constant)
    gives the constant null with that reason, and the ``operons.discordance_fit`` event itself is
    WARNING in the step's event records (the records that become ``logs/events.jsonl`` and the
    manifest warnings), not only a stderr line on another logger."""
    _patch_scaled(monkeypatch, FALLBACKS[reason])
    gt = write_gene_table(chain, tmp_path / "genetable" / "genes.parquet", abundance=True)
    run = run_operons_step(step_inputs(chain, gt), tmp_path / "out", _cfg())
    data = _data(run.records, logging.WARNING)
    assert data["sigma_model_used"] == "constant" and data["reason"] == reason
    assert data["gamma_per_sd"] == 0.0


def _relative_gene_table(chain: Chain, path: Path) -> Path:
    from gippa.io.genes_tsv import read_genes_parquet, write_genes_parquet

    write_gene_table(chain, path, abundance=True)
    df = read_genes_parquet(path)
    df["abundance_type"] = "lfq"
    df["abundance_scale"] = "relative"
    write_genes_parquet(path, df)
    return path


def test_relative_scale_run_reports_no_null(chain: Chain, tmp_path: Path) -> None:
    """Stat S2: on a relative-scale run no p is computed (every OP pair is ``relative_scale``), so
    the event must not report a fitted log-linear null with a gamma-hat. Choice (the spec leaves
    it open; ``sigma_model_used`` is ``loglinear | constant | none`` in the spec): the null used
    is ``none`` and the reason is the existing ``relative_scale`` flag; no fit values."""
    gt = _relative_gene_table(chain, tmp_path / "genetable" / "genes.parquet")
    run = run_operons_step(step_inputs(chain, gt), tmp_path / "out", _cfg())
    data = _data(run.records)
    assert data["sigma_model_used"] == "none"
    assert data["reason"] == "relative_scale"
    for key in ("gamma_per_sd", "gamma_per_log2", "loc_log2", "scale0_log2", "df", "lr_stat",
                "lr_p"):  # fmt: skip
        assert data[key] is None, key
    assert data["df_at_bound"] is False and data["gamma_at_bound"] is False
    json.dumps(data, allow_nan=False)


@pytest.mark.xfail(
    strict=True,
    reason="stat N3 deferred to the sigma_model config key (M6-13 Q2): the step has no way to"
    " pass a model today, so the event reports the module constant (lead decision)",
)
def test_event_sigma_model_is_the_value_passed(
    chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Stat N3: the event's ``sigma_model`` is the model ``discordance`` was actually called with,
    not the module constant ``SIGMA_MODEL``. The step has no ``sigma_model`` option yet (M6-13
    Q2), so the call is forced to ``constant`` at the ``discordance`` the step uses."""
    import gippa.steps.operons as step_mod

    real = step_mod.discordance
    seen: list[str] = []

    def forced(*a: Any, **k: Any) -> Any:
        k["sigma_model"] = "constant"
        seen.append(k["sigma_model"])
        return real(*a, **k)

    monkeypatch.setattr(step_mod, "discordance", forced)
    gt = write_gene_table(chain, tmp_path / "genetable" / "genes.parquet", abundance=True)
    run = run_operons_step(step_inputs(chain, gt), tmp_path / "out", _cfg())
    assert seen == ["constant"]
    data = _data(run.records)
    assert data["sigma_model"] == "constant"
    assert data["sigma_model_used"] == "constant" and data["reason"] == "disabled"


# --------------------------------------------------------------------------- bound flags (S-3)


def _fit(model: str, *, gamma: float, df: float, reason: str | None = None) -> Any:
    from gippa.stats.tdist import ScaledStudentT

    nan = math.nan
    return ScaledStudentT(
        loc=0.0, scale0=0.5, gamma=gamma, df=df, n=40, converged=model != "none", model=model,
        reason=reason, center=12.0, sd=0.7, z_min=-1.5, z_max=1.5, nll=nan, nll_constant=nan,
        lr_stat=nan, lr_p=nan,
    )  # type: ignore[arg-type]  # fmt: skip


@pytest.mark.parametrize("df", [2.0, 50.0])
def test_df_at_bound_is_reported_for_the_constant_null(df: float) -> None:
    """CR S-3 / stat N1: the unchanged nu-on-a-bound rule applies to the constant fallback too."""
    from gippa.steps.operons import discordance_fit_data

    data = discordance_fit_data(_fit("constant", gamma=0.0, df=df, reason="flat_abundance"))
    assert data["df_at_bound"] is True
    assert data["gamma_at_bound"] is False
    inner = discordance_fit_data(_fit("constant", gamma=0.0, df=7.5, reason="flat_abundance"))
    assert inner["df_at_bound"] is False


def test_bound_flags_read_the_module_constants(monkeypatch: pytest.MonkeyPatch) -> None:
    """CR S-3 / stat N2: the flags come from ``GAMMA_BOUNDS`` and ``DF_BOUNDS`` (gippa.stats.tdist),
    not hard-coded 2 / 50 and +-2: with the constants moved, the flags move with them."""
    import gippa.steps.operons as step_mod

    monkeypatch.setattr(tdist, "GAMMA_BOUNDS", (-1.5, 1.5))
    monkeypatch.setattr(tdist, "DF_BOUNDS", (3.0, 40.0))
    for mod in (step_mod,):
        if hasattr(mod, "GAMMA_BOUNDS"):
            monkeypatch.setattr(mod, "GAMMA_BOUNDS", (-1.5, 1.5))
        if hasattr(mod, "DF_BOUNDS"):
            monkeypatch.setattr(mod, "DF_BOUNDS", (3.0, 40.0))
    on = step_mod.discordance_fit_data(_fit("loglinear", gamma=-1.5, df=40.0))
    assert on["gamma_at_bound"] is True and on["df_at_bound"] is True
    off = step_mod.discordance_fit_data(_fit("loglinear", gamma=-2.0, df=50.0))
    assert off["gamma_at_bound"] is False and off["df_at_bound"] is False


# --------------------------------------------------------------------------- CR S-2


def test_run_operons_returns_the_null_and_run_operons_fit_is_gone(with_abundance: Any) -> None:
    """CR S-2 (lead decision): ``run_operons`` itself returns the five values (pairs, operons,
    operon_flags, bin_counts, ``Discordance.scaled``); the ``run_operons_fit`` duplicate is merged
    back."""
    import pandas as pd

    import gippa.steps.operons as step_mod
    from gippa.io.genes_tsv import read_genes_parquet
    from gippa.stats.tdist import ScaledStudentT

    assert not hasattr(step_mod, "run_operons_fit")
    _, ins = with_abundance
    rd = pd.read_parquet
    layers = {k: rd(ins[k]) for k in ("promoter_windows", "terminators", "tx_rates")}
    out = step_mod.run_operons(read_genes_parquet(ins["gene_table"]), rd(ins["genes"]),
                               rd(ins["replicons"]), rd(ins["tus"]), rd(ins["pairs_raw"]),
                               rd(ins["features"]), layers, _cfg())  # fmt: skip
    assert len(out) == 5
    assert isinstance(out[4], ScaledStudentT)
    assert out[4] == _scaled(ins)
