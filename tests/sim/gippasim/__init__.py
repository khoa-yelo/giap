"""gippasim: synthetic-truth generator for GIPPA tests (GIPPA_SPEC.md §11.3).

M0 skeleton: one circular chromosome, protein-coding genes in transcription units, a Bakta-like
JSON annotation and a gene truth table.

The top-level names load lazily, so the sequence-free modules (``layout``, ``expression``, ...)
import without ``gippa`` (the nucleotide simulator uses ``gippa`` writers).
"""

import importlib
import sys
import types
from typing import Any

__all__ = [
    "SIM5M",
    "SIM5M_SEED",
    "SIM100K",
    "SIM100K_SEED",
    "SimConfig",
    "SimResult",
    "StrainConfig",
    "simulate",
]
_LAZY = {name: "gippasim.simulate" for name in __all__} | {"StrainConfig": "gippasim.strain"}


class _Package(types.ModuleType):
    def __setattr__(self, name: str, value: Any) -> None:
        # ``gippasim.simulate`` is the function, as with the former eager import, even after the
        # import system binds the ``gippasim.simulate`` submodule onto the package.
        if name == "simulate" and isinstance(value, types.ModuleType):
            value = value.simulate
        super().__setattr__(name, value)


def __getattr__(name: str) -> Any:
    if name in _LAZY:
        return getattr(importlib.import_module(_LAZY[name]), name)
    raise AttributeError(f"module 'gippasim' has no attribute {name!r}")


sys.modules[__name__].__class__ = _Package
