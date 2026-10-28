"""M6-10 review #243 B2 and B3 at step level: the registered ``sensitivity`` step on the SIM100K
producer chain (docs/specs/M6-10-sensitivity-operon-flags.md, lead decisions 4, 7 and 8;
GIPPA_SPEC.md §8.6, G9).

* B2 (lead decision 8, option a): a variant runs only for the families in ``Variant.families``
  (``None`` = every family). Under ``--sensitivity full`` the ``pair_threshold=*`` variants give
  ``operon_flags`` rows only; ``regions`` and ``islands`` get none, and are not rerun for them.
* B3 (lead decisions 4 and 8): when any one of C2 ``genes``, ``pairs_raw`` or C3 ``features`` is
  absent, ``operon_flags`` is skipped with exactly one ``sensitivity.operon_flags_skipped``
  event (``data.missing`` = the absent keys, sorted), never a KeyError, and the other families
  keep their rows.

``regions`` / ``islands`` are replaced by recording stubs (cheap, and they show which configs a
family is rerun with); ``operon_flags`` is the registered family unless stated otherwise.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import sensitivity_operons_data as so
from helpers.layers_chain import Chain, build_chain

PAIR_KEYS = ("genes", "pairs_raw", "features")
EVENT = "sensitivity.operon_flags_skipped"
PT = tuple(so.variant_name(t) for t in so.UNIOP_THRESHOLDS)


def _S() -> Any:
    from gippa.context import sensitivity as S

    return S


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(Path(tmp_path_factory.mktemp("sens_review_243")))


class _Recorder:
    """A ``(genes, cfg, rng)`` family that records the config of each call and flags one ID."""

    def __init__(self) -> None:
        self.cfgs: list[Any] = []

    def __call__(self, genes: pd.DataFrame, cfg: Any, rng: np.random.Generator
                 ) -> frozenset[str]:  # fmt: skip
        self.cfgs.append(cfg)
        return frozenset({"x"})


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _run(chain: Chain, root: Path, cfg: Any, *, missing: tuple[str, ...] = ()) -> pd.DataFrame:
    from gippa.steps import load_step_modules
    from gippa.steps.registry import run_step

    load_step_modules()
    gt = so.gene_table(so.tables(chain)["genes"], imputed=(1, 2, 3))
    ins = so.write_inputs(chain, root / "in", gt)
    for k in missing:
        ins.pop(k)
    out = root / "tables" / "sensitivity.tsv"
    run_step("sensitivity", cfg, ins, {"sensitivity": out})
    return _read(out)


def _variants_of(rows: pd.DataFrame, family: str) -> list[str]:
    return list(rows.loc[rows["family"] == family, "variant"])


def _stubbed(monkeypatch: pytest.MonkeyPatch, *, real_operons: bool = True
             ) -> dict[str, Any]:  # fmt: skip
    S = _S()
    fams: dict[str, Any] = {"regions": _Recorder(), "islands": _Recorder()}
    fams["operon_flags"] = S.FAMILIES["operon_flags"] if real_operons else _Recorder()
    monkeypatch.setattr(S, "FAMILIES", fams)
    return fams


# --------------------------------------------------------------------------- B2


def test_b2_full_pair_threshold_rows_only_for_operon_flags(
    chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Under ``full`` the ``pair_threshold`` rows exist for ``operon_flags`` and for no other
    family; ``regions`` / ``islands`` are never rerun with a patched UniOP threshold."""
    fams = _stubbed(monkeypatch)
    cfg = so.config(output={"sensitivity": "full"})
    rows = _run(chain, tmp_path, cfg)
    assert _variants_of(rows, "operon_flags") == ["missing=exclude", "top3", *PT]
    for name in ("regions", "islands"):
        assert _variants_of(rows, name) == ["missing=exclude", "top3"], name
        thresholds = {c.annotation.uniop.pair_threshold for c in fams[name].cfgs}
        assert thresholds == {cfg.annotation.uniop.pair_threshold}, name
        assert len(fams[name].cfgs) == 1 + 2, f"{name}: primary + missing=exclude + top3 only"


def test_b2_families_none_runs_for_every_family(
    chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``families = None`` (also on a patch-only variant) reruns every family; a tuple restricts
    the variant to the families it names, whatever they are (operons are not special-cased)."""
    S = _S()
    _stubbed(monkeypatch, real_operons=False)
    vs = [
        S.Variant("all", {"abundance.missing": "exclude"}),
        S.Variant("all_none", {"abundance.missing": "exclude"}, families=None),
        S.Variant("only_regions", {"abundance.missing": "exclude"}, families=("regions",)),
        S.Variant("patch_all", {}, config_patch={so.PATCH_KEY: 0.3}),
        S.Variant("patch_islands", {}, config_patch={so.PATCH_KEY: 0.7}, families=("islands",)),
    ]
    monkeypatch.setattr(S, "variants", lambda cfg, tier: list(vs))
    rows = _run(chain, tmp_path, so.config(output={"sensitivity": "full"}))
    got = {(f, v) for f, v in zip(rows["family"], rows["variant"], strict=True)}
    fams = ("regions", "islands", "operon_flags")
    want = {(f, v) for f in fams for v in ("all", "all_none", "patch_all")}
    want |= {("regions", "only_regions"), ("islands", "patch_islands")}
    assert got == want
    assert len(rows) == len(want)


# --------------------------------------------------------------------------- B3


@pytest.mark.parametrize(
    "missing", [("genes",), ("pairs_raw",), ("features",), PAIR_KEYS],
    ids=["genes", "pairs_raw", "features", "all"],
)  # fmt: skip
def test_b3_missing_pair_table_skips_operon_flags(
    chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture, missing: tuple[str, ...],
) -> None:  # fmt: skip
    """Any one absent pair table skips the family (no KeyError), logs exactly one skip event
    naming the absent keys, writes no ``operon_flags`` / ``pair_threshold`` row and leaves the
    other families' rows as they are."""
    fams = _stubbed(monkeypatch)
    cfg = so.config(output={"sensitivity": "full"})
    with caplog.at_level(logging.DEBUG):
        rows = _run(chain, tmp_path, cfg, missing=missing)
    events = [r for r in caplog.records if getattr(r, "event", None) == EVENT]
    assert len(events) == 1, [r.getMessage() for r in events]
    assert list(events[0].data["missing"]) == sorted(missing)
    assert "operon_flags" not in set(rows["family"])
    assert not rows["variant"].isin(PT).any()
    for name in ("regions", "islands"):
        assert _variants_of(rows, name) == ["missing=exclude", "top3"], name
        assert len(fams[name].cfgs) == 1 + 2, name
    assert set(rows["jaccard"]) == {"1"}


def test_b3_no_skip_event_with_every_pair_table(
    chain: Chain, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:  # fmt: skip
    """Control: with all three pair tables present the family runs and no skip event is logged."""
    _stubbed(monkeypatch)
    with caplog.at_level(logging.DEBUG):
        rows = _run(chain, tmp_path, so.config())
    assert not [r for r in caplog.records if getattr(r, "event", None) == EVENT]
    assert _variants_of(rows, "operon_flags") == ["missing=exclude"]
