"""#217 review N3: M5-11 fixture edges mutation testing left open (decision 41, place_oric)."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest


def _oric_mid(work: Path) -> float:
    f = pd.read_parquet(work / "features.parquet")
    o = f[f["type"].astype(str) == "oric"]
    length = int(
        pd.read_parquet(work / "replicons.parquet").set_index("contig").loc["chr1", "length"]
    )
    return float((o["start"].iloc[0] + o["end"].iloc[0]) / 2 / length)


def test_the_oric_frac_keyword_wins_over_layout(tmp_path: Path) -> None:
    from gippasim.fixtures import build_m4_fixture

    build_m4_fixture(tmp_path, seed=3, size="sim1000", layout={"oric_frac": 0.3}, oric_frac=0.8)
    assert abs(_oric_mid(tmp_path) - 0.8) < 0.02


def _with_other_rows() -> object:
    from gippasim.fixtures import fixture_layout
    from gippasim.layout import place_oric

    lay = place_oric(fixture_layout(3, "sim1000", None), 0.5)
    o = lay.features.iloc[0]
    late = {**o.to_dict(), "type": "crispr", "feature_id": "CR_2", "start": 900_000, "end": 900_100}
    early = {**o.to_dict(), "type": "crispr", "feature_id": "CR_1", "start": 100, "end": 200}
    return replace(lay, features=pd.DataFrame([late, early]).astype(lay.features.dtypes.to_dict()))


def test_place_oric_keeps_other_rows_sorted_and_replaces_old_oric() -> None:
    from gippasim.layout import place_oric

    lay = place_oric(place_oric(_with_other_rows(), 0.2), 0.8)  # type: ignore[arg-type]
    f = lay.features
    assert list(f["type"].astype(str)).count("oric") == 1
    assert sorted(f.loc[f["type"].astype(str) == "crispr", "feature_id"]) == ["CR_1", "CR_2"]
    assert list(f["start"]) == sorted(f["start"])


@pytest.mark.parametrize("frac", [1.0, -0.1, float("nan")])
def test_place_oric_rejects_a_fraction_outside_0_1(frac: float) -> None:
    from gippasim.fixtures import fixture_layout
    from gippasim.layout import place_oric

    with pytest.raises(ValueError, match="oric_frac"):
        place_oric(fixture_layout(3, "sim1000", None), frac)
