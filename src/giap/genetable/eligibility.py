"""§7.3 gene eligibility shared by the region (M6) and function (M8) tests."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pandas as pd

B8 = npt.NDArray[np.bool_]


def eligible(genes: pd.DataFrame) -> B8:
    """Non-partial rows (NA partial counts as not partial) with a finite ``p_detect``."""
    p = pd.to_numeric(genes["p_detect"], errors="coerce").to_numpy(float)
    part = genes["partial"].astype("boolean").fillna(False).to_numpy(bool)
    return np.asarray(~part & np.isfinite(p), dtype=bool)
