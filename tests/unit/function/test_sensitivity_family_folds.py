"""M8-16 review round, code review B2a (lead decision 4 gap 2, lead decision 5): the
``enriched_categories`` family passes the M4-13 TU ``folds`` to the M8-04 detection ORA, so its
``ora_*`` IDs follow the function step's recalibrated q (docs/specs/M8-16-sensitivity-function-
determinants.md item 1; GIPPA_SPEC.md §7.7.3).

The world: 200 one-gene TUs whose ``p_detect`` (0.1 / 0.9) says nothing about detection (half of
each class detected). COG K holds 12 genes at p_detect 0.1, 6 of them detected. Without folds the
raw q makes 6 detections surprising (``ora_over``, q ~ 0.001); the cross-fitted recalibration
flattens q to the base rate and K is no longer enriched. GSEA is switched off by an engine that
returns p = 0.9 for every term.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import sensitivity_function_data as sf

N = 200


def _rng() -> np.random.Generator:
    from gippa.rng import for_step

    return for_step(f"sensitivity/{sf.EC}", sf.SEED)


def _world() -> tuple[pd.DataFrame, list[Any], dict[str, int]]:
    p = [0.1] * 12 + [0.1 if i % 2 else 0.9 for i in range(12, N)]
    det = [i % 2 == 0 for i in range(12)] + [(i // 2) % 2 == 0 for i in range(12, N)]
    ab = [0.01 * i if d else math.nan for i, d in enumerate(det)]
    genes = sf.gene_table(ab, det, p)
    colls = [sf.gene_sets("cog", {"K": sf.tags(range(12))})]
    folds = {f"TU_{i + 1}": i % 5 for i in range(N)}
    return genes, colls, folds


def _call(genes: pd.DataFrame, colls: list[Any], folds: dict[str, int] | None) -> frozenset[str]:
    res = sf.ec_family()(genes, sf.base_config(), _rng(), collections=colls,
                         tus=sf.tus_frame(genes), folds=folds,
                         engine=sf.constant_engine(0.9))  # fmt: skip
    return frozenset(res)


def _oracle(genes: pd.DataFrame, colls: list[Any], folds: dict[str, int] | None) -> frozenset[str]:
    from gippa.function.ora import detection_enrichment

    cfg = sf.base_config()
    lo, hi = cfg.function.ora_size
    alpha = cfg.stats.alpha.enrichment
    res = detection_enrichment(genes, colls, (lo, hi), alpha, folds=folds)
    q = res["q_value"].astype(float)
    keep = res[q.notna() & (q <= alpha)]
    cols = zip(keep["collection"], keep["term"], keep["direction"], strict=True)
    return frozenset(f"{c}:{t}:ora_{d}" for c, t, d in cols)


def test_fixture_folds_change_the_ora() -> None:
    genes, colls, folds = _world()
    assert _oracle(genes, colls, None) == {"cog:K:ora_over"}
    assert _oracle(genes, colls, folds) == frozenset()


def test_folds_recalibrate_the_family_ora() -> None:
    """With folds the family's ORA IDs are the recalibrated ORA's (K not enriched); without them
    the raw-q ORA's (K over)."""
    genes, colls, folds = _world()
    assert _call(genes, colls, folds) == _oracle(genes, colls, folds) == frozenset()
    assert _call(genes, colls, None) == _oracle(genes, colls, None) == {"cog:K:ora_over"}


@pytest.mark.parametrize("folds", [None, "given"])
def test_folds_reach_detection_enrichment(
    folds: str | None, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The family hands ``detection_enrichment`` exactly the folds it was given (``None`` stays
    ``None``)."""
    import gippa.function.sensitivity as FS

    genes, colls, f = _world()
    want = f if folds else None
    seen: list[Any] = []
    real = FS.detection_enrichment

    def spy(*a: Any, **k: Any) -> Any:
        seen.append(k.get("folds"))
        return real(*a, **k)

    monkeypatch.setattr(FS, "detection_enrichment", spy)
    _call(genes, colls, want)
    assert seen, "the family never ran the detection ORA"
    assert all((s is None) if want is None else (dict(s) == want) for s in seen), seen
