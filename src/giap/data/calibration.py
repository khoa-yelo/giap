"""Loader for the bundled ``calibration.yaml`` (GIAP_SPEC.md §7.1.3, §15.4; M2-12)."""

from __future__ import annotations

from importlib import resources
from typing import Any

import yaml


def load() -> dict[str, Any]:
    """The calibrated constants; the single reader of ``giap/data/calibration.yaml``."""
    text = resources.files("giap.data").joinpath("calibration.yaml").read_text(encoding="utf-8")
    data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError("calibration.yaml: expected a mapping at the top level")
    return data
