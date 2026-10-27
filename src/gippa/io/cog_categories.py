"""The bundled COG category table ``gippa/data/cog_categories.tsv`` (M9-index decision 8).

Package data, not user input: read once per process. Columns ``letter``, ``super_category``,
``description``; ``#`` lines are comments.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import cache
from importlib.resources import files
from types import MappingProxyType


@cache
def cog_descriptions() -> Mapping[str, str]:
    """COG letter -> its NCBI description (``"J"`` -> ``"Translation, ribosomal ..."``)."""
    text = files("gippa.data").joinpath("cog_categories.tsv").read_text(encoding="utf-8")
    rows = [ln.split("\t") for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    header = rows[0]
    li, di = header.index("letter"), header.index("description")
    return MappingProxyType({r[li]: r[di] for r in rows[1:]})
