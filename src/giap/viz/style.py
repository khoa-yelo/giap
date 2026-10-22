"""The one style module for every GIAP plot (GIAP_SPEC.md §10.1; M9-01, M9-index decisions 6-10).

Pure constants and functions, no file I/O. ``apply()`` sets the rcParams once per process (the
``plots`` step and every test fixture call it); ``new_figure`` makes a figure at one of the two
§10.1 widths. Figures are created without pyplot so no global figure manager holds them.
"""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from types import MappingProxyType
from typing import Any, Literal, cast

import matplotlib
from cycler import cycler
from matplotlib.figure import Figure

SINGLE_IN: float = 3.5
DOUBLE_IN: float = 7.2

OKABE_ITO: Mapping[str, str] = MappingProxyType(
    {
        "black": "#000000",
        "orange": "#E69F00",
        "sky_blue": "#56B4E9",
        "bluish_green": "#009E73",
        "yellow": "#F0E442",
        "blue": "#0072B2",
        "vermillion": "#D55E00",
        "reddish_purple": "#CC79A7",
        "grey": "#999999",
    }
)
DETECTION: Mapping[str, str] = MappingProxyType(
    {"detected": "#0072B2", "not_detected": "#BBBBBB", "unestimated": "#E69F00"}
)
#: ``genes.tsv`` ``detection_status`` -> colour (decision 9; ``partial_excluded`` is hatched)
STATUS_COLOUR: Mapping[str, str] = MappingProxyType(
    {
        "detected": DETECTION["detected"],
        "not_detected": DETECTION["not_detected"],
        "partial_excluded": DETECTION["not_detected"],
        "single_peptide_unestimated": DETECTION["unestimated"],
    }
)
REGION: Mapping[str, str] = MappingProxyType(
    {
        "prophage": "#D55E00",
        "plasmid": "#009E73",
        "bgc": "#CC79A7",
        "crispr": "#56B4E9",
        "amr": "#F0E442",
    }
)
#: island direction -> (colour, hatch) (decision 10)
ISLAND: Mapping[str, tuple[str, str]] = MappingProxyType(
    {"high": ("#0072B2", "///"), "low": ("#D55E00", "\\\\\\")}
)

_GROUPS = {
    "information": "JAKLB",
    "cellular": "DYVTMNZWUOX",
    "metabolism": "CGEFHIPQ",
    "poorly": "RS",
}
#: COG letter -> super-category (the four NCBI groups; ``data/cog_categories.tsv``, decision 8)
COG_SUPER: Mapping[str, str] = MappingProxyType(
    {letter: group for group, letters in _GROUPS.items() for letter in letters}
)
COG_SUPER_COLOUR: Mapping[str, str] = MappingProxyType(
    {"information": "#0072B2", "cellular": "#009E73", "metabolism": "#E69F00", "poorly": "#999999"}
)
CONTINUOUS: str = "viridis"
LOCUS_ABUNDANCE: str = "cividis"

#: the category cycle: Okabe-Ito order without black and grey
_CYCLE = [v for k, v in OKABE_ITO.items() if k not in ("black", "grey")]

RC: Mapping[str, Any] = MappingProxyType(
    {
        "font.family": ["DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "figure.titlesize": 8,
        "svg.fonttype": "none",
        "svg.hashsalt": "giap",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "path.simplify": True,
    }
)


def apply() -> None:
    """Select the Agg backend and set the GIAP rcParams (behaviour 1)."""
    matplotlib.use("Agg", force=True)
    # RcParams is typed with literal keys; RC is checked by the style tests instead
    cast("MutableMapping[str, Any]", matplotlib.rcParams).update(RC)
    matplotlib.rcParams["axes.prop_cycle"] = cycler(color=_CYCLE)


def new_figure(width: Literal["single", "double"], height_in: float, **kw: Any) -> Figure:
    """A figure ``SINGLE_IN`` or ``DOUBLE_IN`` wide and ``height_in`` tall."""
    if width not in ("single", "double"):
        raise ValueError(f"width must be 'single' or 'double', not {width!r}")
    w = SINGLE_IN if width == "single" else DOUBLE_IN
    return Figure(figsize=(w, height_in), **kw)


def cog_colour(letter: str) -> str:
    """The super-category colour of one COG letter; unknown letters are poorly characterised."""
    return COG_SUPER_COLOUR[COG_SUPER.get(letter, "poorly")]
