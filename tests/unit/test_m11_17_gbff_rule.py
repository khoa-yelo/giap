"""M11-17 ``inputs.bakta_gbff``: ``expected_gbff`` examples 1-4 and the AC 7 property,
``bgc_gbff_check`` examples 5-14 and AC 1 (docs/specs/M11-17-preflight-gbff-warning.md "Rule",
"Examples"; GIPPA_SPEC.md §5.3, §7.0, §7.1.1).

Interface assumed (spec "Rule"):
- ``gippa.preflight.gbff.expected_gbff(bakta_json: Path) -> tuple[Path, Path]``;
- ``gippa.preflight.gbff.bgc_gbff_check(cfg: GippaConfig, annotation: Path | None)
  -> Check | None``;
- ``gippa.preflight.gbff.RULE_ID == "inputs.bakta_gbff"``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest
from hypothesis import given
from hypothesis import strategies as st

from gippa.config import GippaConfig, resolve_config
from gippa.preflight.checks import Check

RULE_ID = "inputs.bakta_gbff"
RULE_TEXT = (
    "bgc module on and a precomputed Bakta JSON <stem>.json[.gz] -> <stem>.gbff or "
    "<stem>.gbff.gz must exist next to it, else warn (antiSMASH is skipped: layer_missing)"
)
MSG_BGC_OFF = "not evaluated: the bgc module is off"
MSG_REGIONS = "not evaluated: annotation.precomputed.regions replaces antiSMASH"


def _warn_message(js: Path, gbff: Path, gz: Path) -> str:
    return (
        f"no GBFF next to the precomputed Bakta JSON {js}: expected {gbff} (or {gz}); antiSMASH "
        "will be skipped (layer_missing:antismash) and the bgc module will report no "
        "biosynthetic gene clusters. Put Bakta's GBFF output there, or leave bgc out of --modules."
    )


def expected(p: Path) -> tuple[Path, Path]:
    from gippa.preflight.gbff import expected_gbff

    out = expected_gbff(p)
    assert isinstance(out, tuple) and len(out) == 2
    return out


def check(cfg: GippaConfig, annotation: Path | None = None) -> Check | None:
    from gippa.preflight.gbff import bgc_gbff_check

    c = bgc_gbff_check(cfg, annotation)
    assert c is None or isinstance(c, Check)
    return c


def cfg_with(*overrides: str) -> GippaConfig:
    return resolve_config(overrides=list(overrides))


def bakta_cfg(js: Path, *extra: str) -> GippaConfig:
    return cfg_with(f"annotation.precomputed.bakta_json={json.dumps(str(js))}", *extra)


def touch(p: Path) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(b"LOCUS\n")
    return p


def assert_shape(c: Check) -> None:
    assert c.rule_id == RULE_ID
    assert c.input == "inputs"
    assert c.level in {"ok", "info", "warn"}
    assert c.exit_code == 0
    assert c.rule == RULE_TEXT


def assert_value(c: Check, js: Path, gbff: str | None) -> None:
    v: Any = c.value
    assert isinstance(v, dict)
    assert set(v) == {"bakta_json", "expected", "gbff"}
    assert v["bakta_json"] == str(js)
    stem = re.sub(r"(?i)\.json(\.gz)?$", "", str(js))
    assert v["expected"] == [stem + ".gbff", stem + ".gbff.gz"]
    assert v["gbff"] == gbff


# ---------------------------------------------------------------- expected_gbff


def test_rule_id_constant() -> None:
    from gippa.preflight.gbff import RULE_ID as rid

    assert rid == RULE_ID
    assert re.fullmatch(r"[a-z]+\.[a-z_]+", rid)


@pytest.mark.parametrize(
    ("given_path", "want"),
    [
        ("/d/x.json", ("/d/x.gbff", "/d/x.gbff.gz")),  # example 1
        ("/d/x.json.gz", ("/d/x.gbff", "/d/x.gbff.gz")),  # example 2
        ("/d/sample.bakta.JSON.GZ", ("/d/sample.bakta.gbff", "/d/sample.bakta.gbff.gz")),  # 3
        ("rel/x.Json", ("rel/x.gbff", "rel/x.gbff.gz")),  # example 4
    ],
)
def test_examples_1_to_4_expected_gbff(given_path: str, want: tuple[str, str]) -> None:
    got = expected(Path(given_path))
    assert (str(got[0]), str(got[1])) == want
    assert not got[0].is_absolute() or Path(given_path).is_absolute()


def test_example_4_relative_stays_relative() -> None:
    got = expected(Path("rel/x.Json"))
    assert not got[0].is_absolute() and not got[1].is_absolute()


_STEM = st.text(
    alphabet=st.characters(whitelist_categories=("Ll", "Lu", "Nd"), whitelist_characters="._-"),
    min_size=1,
    max_size=20,
).filter(lambda s: s not in {".", ".."})
_SUFFIX = st.sampled_from([".json", ".json.gz"]).flatmap(
    lambda s: st.lists(st.booleans(), min_size=len(s), max_size=len(s)).map(
        lambda ups: "".join(ch.upper() if u else ch for ch, u in zip(s, ups, strict=True))
    )
)


@given(stem=_STEM, suffix=_SUFFIX, absolute=st.booleans())
def test_ac7_expected_gbff_keeps_dir_and_stem(stem: str, suffix: str, absolute: bool) -> None:
    parent = Path("/data/run 1") if absolute else Path("sub")
    gbff, gz = expected(parent / (stem + suffix))
    assert gbff.parent == parent and gz.parent == parent
    assert gbff.name == stem + ".gbff"
    assert gz.name == stem + ".gbff.gz"


def test_expected_gbff_is_pure(tmp_path: Path) -> None:
    js = tmp_path / "none" / "x.json"  # the directory does not exist: no I/O, no error
    assert expected(js) == (tmp_path / "none" / "x.gbff", tmp_path / "none" / "x.gbff.gz")


# ---------------------------------------------------------------- bgc_gbff_check


def test_example_5_gbff_present_is_ok(tmp_path: Path) -> None:
    js = touch(tmp_path / "x.json")
    touch(tmp_path / "x.gbff")
    c = check(bakta_cfg(js))
    assert c is not None
    assert_shape(c)
    assert c.level == "ok"
    assert c.message == ""
    assert_value(c, js, str(tmp_path / "x.gbff"))


def test_example_6_gbff_gz_present_is_ok(tmp_path: Path) -> None:
    js = touch(tmp_path / "x.json")
    touch(tmp_path / "x.gbff.gz")
    c = check(bakta_cfg(js))
    assert c is not None
    assert_shape(c)
    assert c.level == "ok"
    assert c.message == ""
    assert_value(c, js, str(tmp_path / "x.gbff.gz"))


def test_example_6_json_gz_with_gbff_gz(tmp_path: Path) -> None:
    js = touch(tmp_path / "bakta_like.JSON.GZ")
    touch(tmp_path / "bakta_like.gbff.gz")
    c = check(bakta_cfg(js))
    assert c is not None and c.level == "ok"
    assert_value(c, js, str(tmp_path / "bakta_like.gbff.gz"))


def test_example_7_both_present_prefers_plain_gbff(tmp_path: Path) -> None:
    js = touch(tmp_path / "x.json")
    touch(tmp_path / "x.gbff")
    touch(tmp_path / "x.gbff.gz")
    c = check(bakta_cfg(js))
    assert c is not None and c.level == "ok"
    assert_value(c, js, str(tmp_path / "x.gbff"))


def test_example_8_missing_warns_with_the_expected_path(tmp_path: Path) -> None:
    js = touch(tmp_path / "x.json")
    c = check(bakta_cfg(js))
    assert c is not None
    assert_shape(c)
    assert c.level == "warn"
    assert_value(c, js, None)
    gbff, gz = tmp_path / "x.gbff", tmp_path / "x.gbff.gz"
    assert str(gbff) in c.message and str(gz) in c.message
    assert "layer_missing:antismash" in c.message
    assert c.message == _warn_message(js, gbff, gz)


def test_example_8_json_need_not_exist(tmp_path: Path) -> None:
    js = tmp_path / "missing.json"  # existence of the JSON is Nextflow's checkIfExists, not ours
    c = check(bakta_cfg(js))
    assert c is not None and c.level == "warn"
    assert str(tmp_path / "missing.gbff") in c.message


def test_example_9_bgc_off_is_info(tmp_path: Path) -> None:
    js = touch(tmp_path / "x.json")
    c = check(bakta_cfg(js, 'modules=["core","function"]'))
    assert c is not None
    assert_shape(c)
    assert c.level == "info"
    assert c.message == MSG_BGC_OFF
    assert c.value is None


def test_example_9_bgc_off_wins_over_regions(tmp_path: Path) -> None:
    js = touch(tmp_path / "x.json")
    regions = touch(tmp_path / "regions.tsv")
    c = check(
        bakta_cfg(
            js,
            'modules=["core","mobile"]',
            f"annotation.precomputed.regions={json.dumps(str(regions))}",
        )
    )
    assert c is not None and c.level == "info" and c.message == MSG_BGC_OFF


def test_example_10_precomputed_regions_is_info(tmp_path: Path) -> None:
    js = touch(tmp_path / "x.json")
    regions = touch(tmp_path / "regions.tsv")
    c = check(bakta_cfg(js, f"annotation.precomputed.regions={json.dumps(str(regions))}"))
    assert c is not None
    assert_shape(c)
    assert c.level == "info"
    assert c.message == MSG_REGIONS
    assert c.value is None


@pytest.mark.parametrize("name", [None, "a.gff3", "a.gff", "a.gbff", "a.gbk", "a.gff3.gz"])
def test_example_11_no_bakta_json_is_not_applicable(tmp_path: Path, name: str | None) -> None:
    ann = None if name is None else touch(tmp_path / name)
    assert check(cfg_with(), ann) is None


def test_example_11_not_applicable_even_with_bgc_off(tmp_path: Path) -> None:
    assert check(cfg_with('modules=["core"]'), None) is None


def test_example_12_annotation_json_gz_without_config(tmp_path: Path) -> None:
    ann = touch(tmp_path / "x.json.gz")
    c = check(cfg_with(), ann)
    assert c is not None
    assert_shape(c)
    assert c.level == "warn"
    assert_value(c, ann, None)
    assert str(tmp_path / "x.gbff") in c.message


@pytest.mark.parametrize("name", ["x.json", "x.JSON", "x.Json.Gz"])
def test_example_12_annotation_json_with_gbff_is_ok(tmp_path: Path, name: str) -> None:
    ann = touch(tmp_path / name)
    touch(tmp_path / "x.gbff")
    c = check(cfg_with(), ann)
    assert c is not None and c.level == "ok"
    assert_value(c, ann, str(tmp_path / "x.gbff"))


def test_example_13_configured_json_wins(tmp_path: Path) -> None:
    a = touch(tmp_path / "a.json")
    touch(tmp_path / "a.gbff")
    b = touch(tmp_path / "b.json")
    c = check(bakta_cfg(a), b)
    assert c is not None and c.level == "ok"
    assert_value(c, a, str(tmp_path / "a.gbff"))
    # and the other way round: the configured JSON without a GBFF warns although b has one
    touch(tmp_path / "b.gbff")
    (tmp_path / "a.gbff").unlink()
    c2 = check(bakta_cfg(a), b)
    assert c2 is not None and c2.level == "warn"
    assert str(tmp_path / "a.gbff") in c2.message


def test_example_14_other_case_gbff_does_not_count(tmp_path: Path) -> None:
    probe = tmp_path / "CaseProbe"
    probe.write_text("")
    if (tmp_path / "caseprobe").exists():
        pytest.skip("case-insensitive file system")
    js = touch(tmp_path / "x.json")
    touch(tmp_path / "x.GBFF")
    c = check(bakta_cfg(js))
    assert c is not None and c.level == "warn"


def test_gbff_in_another_directory_does_not_count(tmp_path: Path) -> None:
    js = touch(tmp_path / "a" / "x.json")
    touch(tmp_path / "b" / "x.gbff")
    c = check(bakta_cfg(js))
    assert c is not None and c.level == "warn"


@pytest.mark.parametrize(
    ("modules", "files", "regions"),
    [
        (m, f, r)
        for m in (None, '["core"]', '["core","bgc"]')
        for f in ((), ("x.gbff",), ("x.gbff.gz",))
        for r in (False, True)
    ],
)
def test_ac1_never_error(
    tmp_path: Path, modules: str | None, files: tuple[str, ...], regions: bool
) -> None:
    js = touch(tmp_path / "x.json")
    for f in files:
        touch(tmp_path / f)
    extra = [] if modules is None else [f"modules={modules}"]
    if regions:
        extra.append(f"annotation.precomputed.regions={json.dumps(str(touch(tmp_path / 'r.tsv')))}")
    c = check(bakta_cfg(js, *extra))
    assert c is not None
    assert_shape(c)
    assert c.level != "error" and c.exit_code == 0


def test_check_is_deterministic(tmp_path: Path) -> None:
    js = touch(tmp_path / "x.json")
    cfg = bakta_cfg(js)
    assert check(cfg) == check(cfg)
