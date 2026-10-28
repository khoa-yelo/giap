"""``gippa.viz.style``: the §10.1 palettes, widths and rcParams (M9-01 behaviour 1-2, examples 1-2).

Interface assumed beyond the spec text: ``COG_SUPER`` maps every letter A-Z to a super-category
key (the key names are the implementer's); ``src/gippa/data/cog_categories.tsv`` has a ``#``
provenance line naming NCBI COG, then a header with at least ``letter`` and ``description``.
"""

from __future__ import annotations

import string
from pathlib import Path

import matplotlib
import pytest

from gippa.viz import style

ROOT = Path(__file__).resolve().parents[3]
GENES_SCHEMA_STATUS = {"detected", "not_detected", "single_peptide_unestimated", "partial_excluded"}


def test_example_1_unestimated_status_colour() -> None:
    assert style.STATUS_COLOUR["single_peptide_unestimated"] == "#E69F00"


def test_status_colours_cover_the_detection_status_enum() -> None:
    assert set(style.STATUS_COLOUR) == GENES_SCHEMA_STATUS
    assert style.STATUS_COLOUR["detected"] == "#0072B2"
    assert style.STATUS_COLOUR["not_detected"] == "#BBBBBB"
    assert style.STATUS_COLOUR["partial_excluded"] == "#BBBBBB"


@pytest.mark.parametrize(
    ("letter", "colour"), [("J", "#0072B2"), ("C", "#E69F00"), ("S", "#999999"), ("D", "#009E73")]
)
def test_example_2_cog_colours(letter: str, colour: str) -> None:
    assert style.cog_colour(letter) == colour


def test_widths() -> None:
    assert style.SINGLE_IN == 3.5
    assert style.DOUBLE_IN == 7.2


def test_okabe_ito_palette_exact() -> None:
    assert dict(style.OKABE_ITO) == {
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


def test_section_10_1_colour_table_exact() -> None:
    assert dict(style.DETECTION) == {
        "detected": "#0072B2",
        "not_detected": "#BBBBBB",
        "unestimated": "#E69F00",
    }
    assert dict(style.REGION) == {
        "prophage": "#D55E00",
        "plasmid": "#009E73",
        "bgc": "#CC79A7",
        "crispr": "#56B4E9",
        "amr": "#F0E442",
    }
    assert dict(style.ISLAND) == {"high": ("#0072B2", "///"), "low": ("#D55E00", "\\\\\\")}
    assert style.CONTINUOUS == "viridis"
    assert style.LOCUS_ABUNDANCE == "cividis"


GROUPS = {
    "information": "JAKLB",
    "cellular": "DYVTMNZWUO",
    "metabolism": "CGEFHIPQ",
    "poorly": "RS",
}
GROUP_COLOUR = {
    "information": "#0072B2",
    "cellular": "#009E73",
    "metabolism": "#E69F00",
    "poorly": "#999999",
}


def test_cog_super_categories_follow_the_ncbi_groups() -> None:
    assert set(style.COG_SUPER) == set(string.ascii_uppercase)
    assert set(style.COG_SUPER_COLOUR.values()) == set(GROUP_COLOUR.values())
    for group, letters in GROUPS.items():
        keys = {style.COG_SUPER[x] for x in letters}
        assert len(keys) == 1, (group, keys)
        for x in letters:
            assert style.cog_colour(x) == GROUP_COLOUR[group], x
            assert style.COG_SUPER_COLOUR[style.COG_SUPER[x]] == GROUP_COLOUR[group]
    assert len({style.COG_SUPER[g[0]] for g in GROUPS.values()}) == 4


@pytest.mark.parametrize("letter", ["?", "none", "", "1"])
def test_unknown_cog_letter_is_poorly_characterised_grey(letter: str) -> None:
    assert style.cog_colour(letter) == "#999999"


def test_bundled_cog_categories_table() -> None:
    path = ROOT / "src" / "gippa" / "data" / "cog_categories.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    comments = [x for x in lines if x.startswith("#")]
    assert comments and any("NCBI" in x and "COG" in x for x in comments)
    body = [x.split("\t") for x in lines if x and not x.startswith("#")]
    header, rows = body[0], body[1:]
    assert {"letter", "description"} <= set(header)
    li, di = header.index("letter"), header.index("description")
    assert sorted(r[li] for r in rows) == list(string.ascii_uppercase)
    assert all(r[di].strip() for r in rows)


def test_apply_sets_agg_fonts_sizes_and_deterministic_svg() -> None:
    style.apply()
    rc = matplotlib.rcParams
    assert matplotlib.get_backend().lower() == "agg"
    assert list(rc["font.family"]) == ["DejaVu Sans"]
    assert rc["font.size"] == 7
    assert rc["axes.titlesize"] == 8
    assert rc["xtick.labelsize"] == 6
    assert rc["ytick.labelsize"] == 6
    assert rc["svg.fonttype"] == "none"
    assert rc["svg.hashsalt"] == "gippa"
    assert rc["pdf.fonttype"] == 42


def test_apply_sets_the_okabe_ito_cycle() -> None:
    style.apply()
    cycle = [c.upper() for c in matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]]
    order = [c.upper() for c in style.OKABE_ITO.values()]
    assert len(cycle) >= 7
    # the cycle is the Okabe-Ito order (a subsequence of it, e.g. without black or grey)
    positions = [order.index(c) for c in cycle]
    assert positions == sorted(positions)


@pytest.mark.parametrize(("width", "inches"), [("single", 3.5), ("double", 7.2)])
def test_new_figure_sizes(width: str, inches: float) -> None:
    style.apply()
    fig = style.new_figure(width, 2.5)  # type: ignore[arg-type]
    w, h = fig.get_size_inches()
    assert (round(w, 6), round(h, 6)) == (inches, 2.5)
