"""Search-engine output converters to the generic peptide table (`giap convert`, §5.4)."""

from __future__ import annotations

from collections.abc import Callable, Mapping

from giap.io.converters.base import ConvertResult
from giap.io.converters.diann import convert_diann
from giap.io.converters.generic import convert_generic

__all__ = ["CONVERTERS", "ConvertResult", "convert_diann", "convert_generic"]

CONVERTERS: Mapping[str, Callable[..., ConvertResult]] = {
    "diann": convert_diann,
    "generic": convert_generic,
}
