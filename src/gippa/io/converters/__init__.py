"""Search-engine output converters to the generic peptide table (`gippa convert`, §5.4)."""

from __future__ import annotations

from collections.abc import Callable, Mapping

from gippa.io.converters.base import ConvertResult
from gippa.io.converters.cascadia import convert_cascadia
from gippa.io.converters.diann import convert_diann
from gippa.io.converters.generic import convert_generic

__all__ = [
    "CONVERTERS",
    "ConvertResult",
    "convert_cascadia",
    "convert_diann",
    "convert_generic",
]

CONVERTERS: Mapping[str, Callable[..., ConvertResult]] = {
    "cascadia": convert_cascadia,
    "diann": convert_diann,
    "generic": convert_generic,
}
