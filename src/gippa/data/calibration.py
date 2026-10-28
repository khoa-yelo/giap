"""Loader for the bundled ``calibration.yaml`` (GIPPA_SPEC.md §7.1.3, §15.4; M2-12)."""

from __future__ import annotations

from importlib import resources
from typing import Any

import yaml


def load() -> dict[str, Any]:
    """The calibrated constants; the single reader of ``gippa/data/calibration.yaml``."""
    text = resources.files("gippa.data").joinpath("calibration.yaml").read_text(encoding="utf-8")
    data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError("calibration.yaml: expected a mapping at the top level")
    return data
