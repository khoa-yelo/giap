"""``helpers.byte_identity.svg_diff_report``: the failure message of the R7 byte-identity tests
shows a unified diff of each differing SVG (macOS V18_stat_diagnostics.svg investigation)."""

from __future__ import annotations

from helpers.byte_identity import svg_diff_report


def test_report_names_paths_and_diffs_svg_lines() -> None:
    a = {"V18_x.svg": b'<svg>\n<clipPath id="p0123456789"/>\n</svg>\n', "V18_x.png": b"a"}
    b = {"V18_x.svg": b'<svg>\n<clipPath id="p9876543210"/>\n</svg>\n', "V18_x.png": b"b"}
    msg = svg_diff_report(a, b)
    assert "V18_x.png" in msg and "V18_x.svg" in msg
    assert '-<clipPath id="p0123456789"/>' in msg
    assert '+<clipPath id="p9876543210"/>' in msg


def test_report_truncates_long_diffs() -> None:
    a = {"x.svg": "".join(f"a{i}\n" for i in range(200)).encode()}
    b = {"x.svg": "".join(f"b{i}\n" for i in range(200)).encode()}
    msg = svg_diff_report(a, b, max_lines=40)
    assert "more diff lines" in msg
    assert "a150" not in msg


def test_report_on_identical_trees_lists_nothing() -> None:
    assert svg_diff_report({"x.svg": b"1"}, {"x.svg": b"1"}) == "differing: []"
