"""giapsim: synthetic-truth generator for GIAP tests (GIAP_SPEC.md §11.3).

M0 skeleton: one circular chromosome, protein-coding genes in transcription units, a Bakta-like
JSON annotation and a gene truth table.

The top-level names load lazily, so the sequence-free modules (``layout``, ``expression``, ...)
import without ``giap`` (the nucleotide simulator uses ``giap`` writers).
"""

import importlib
import sys
import types
from typing import Any

__all__ = ["SIM100K", "SIM100K_SEED", "SimConfig", "SimResult", "StrainConfig", "simulate"]
_LAZY = {name: "giapsim.simulate" for name in __all__} | {"StrainConfig": "giapsim.strain"}


class _Package(types.ModuleType):
    def __setattr__(self, name: str, value: Any) -> None:
        # ``giapsim.simulate`` is the function, as with the former eager import, even after the
        # import system binds the ``giapsim.simulate`` submodule onto the package.
        if name == "simulate" and isinstance(value, types.ModuleType):
            value = value.simulate
        super().__setattr__(name, value)


def __getattr__(name: str) -> Any:
    if name in _LAZY:
        return getattr(importlib.import_module(_LAZY[name]), name)
    raise AttributeError(f"module 'giapsim' has no attribute {name!r}")


sys.modules[__name__].__class__ = _Package
