"""``validation/pilot/px_record.py``: the official record -> the run selection (M10-01).

Spec: docs/specs/M10-01-pxd078459-download.md, "Parsing the record", "Selection", examples 1-4,
acceptance criteria 1-2; lead decision 1 (the first three wild-type runs by natural sort).

Interface settled by the test author (record vocabulary: ``tests/helpers/pxd_fixtures.py``):

* ``parse_record(xml: bytes) -> list[RecordFile]`` in record order; a ``RecordFile`` has
  ``file_name`` (the record's file name), ``uri`` (verbatim), ``kind`` (``"raw"`` or
  ``"not_raw"``), ``bytes`` (``int | None``), ``checksum`` (``"md5:<lowercase hex>" | None``) and
  ``sample`` (the annotation text the selection reads, ``None`` without one).
* ``select_wild_type(files) -> Selection``; ``Selection.selected`` is in selection order, each item
  with ``run`` (``WT_1``..``WT_3``), ``file_name`` and ``selected_because``;
  ``Selection.rejected`` lists every other file once, each with ``file_name`` and ``reason``
  (``not_wild_type`` | ``beyond_three`` | ``not_raw``).
* Fewer than three wild-type runs raise ``SelectionError`` (a ``ValueError``) whose message is
  ``fewer than three wild-type runs: <names, comma-separated, natural order>``.
* CLI ``python -I px_record.py RECORD.xml --out MANIFEST.tsv``: exit 0 and the manifest (columns
  of the spec, ``NA`` for missing, ``sha256`` = ``NA``); exit 2 and no manifest for fewer than
  three runs or a record that names no timsTOF instrument.

No statistical claim is made in this module, so it has no calibration test (spec §8.4).
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest
from helpers.pxd_fixtures import (
    BASE,
    MANIFEST_COLUMNS,
    PX_RECORD,
    FileSpec,
    example1_files,
    example2_files,
    example3_files,
    fixture_xml,
    md5_of,
    px_xml,
    read_tsv,
)
from hypothesis import given, settings
from hypothesis import strategies as st


@pytest.fixture(scope="module")
def px() -> ModuleType:
    from helpers.pxd_fixtures import load

    return load(PX_RECORD)


def _sel(px: ModuleType, files: list[FileSpec]):  # type: ignore[no-untyped-def]
    return px.select_wild_type(px.parse_record(px_xml(files)))


def _selected(sel) -> list[tuple[str, str, str]]:  # type: ignore[no-untyped-def]
    return [(s.run, s.file_name, s.selected_because) for s in sel.selected]


def _rejected(sel) -> dict[str, str]:  # type: ignore[no-untyped-def]
    out: dict[str, str] = {}
    for r in sel.rejected:
        assert r.file_name not in out, f"{r.file_name} rejected twice"
        out[r.file_name] = r.reason
    return out


def _wt(*names: str) -> list[FileSpec]:
    return [FileSpec(n, f"{BASE}/{n}", size=100 + i, md5=md5_of(n)) for i, n in enumerate(names)]


# ------------------------------------------------------------------ fixtures are committed


@pytest.mark.parametrize("n", [1, 2, 3])
def test_xml_fixtures_are_committed_and_small(n: int) -> None:
    """Acceptance 1: synthetic fixtures under tests/data/pilot/, a few kB each."""
    path = fixture_xml(n)
    assert path.is_file(), path
    assert path.stat().st_size < 20_000, path.stat().st_size


# ------------------------------------------------------------------ parse_record


def test_parse_example1_types_sizes_checksums(px: ModuleType) -> None:
    files = px.parse_record(fixture_xml(1).read_bytes())
    by = {f.file_name: f for f in files}
    assert set(by) == {s.name for s in example1_files()}, sorted(by)
    for spec in example1_files():
        f = by[spec.name]
        assert f.uri == spec.uri
        assert f.bytes == spec.size
        if spec.raw:
            assert f.kind == "raw", (spec.name, f.kind)
            assert f.checksum == f"md5:{spec.md5}", f.checksum
        else:
            assert f.kind == "not_raw", "the Spectronaut .sne is a search file, never downloaded"
            assert f.checksum is None


def test_parse_keeps_record_order(px: ModuleType) -> None:
    files = px.parse_record(fixture_xml(3).read_bytes())
    assert [f.file_name for f in files] == [s.name for s in example3_files()]


def test_parse_example2_reads_the_sample_annotation(px: ModuleType) -> None:
    files = px.parse_record(fixture_xml(2).read_bytes())
    samples = {f.file_name: f.sample or "" for f in files}
    for spec in example2_files():
        assert spec.genotype and spec.genotype in samples[spec.name], samples


def test_parse_missing_size_and_checksum_are_none(px: ModuleType) -> None:
    files = px.parse_record(px_xml([FileSpec("WT_1.d.zip", f"{BASE}/WT_1.d.zip")]))
    assert len(files) == 1
    assert files[0].bytes is None and files[0].checksum is None


def test_parse_keeps_ftp_uris_verbatim(px: ModuleType) -> None:
    uri = "ftp://download.iprox.cn/IPX0012345000/IPX0012345001/WT_1.d.zip"
    files = px.parse_record(px_xml([FileSpec("WT_1.d.zip", uri, size=7)]))
    assert files[0].uri == uri


# ------------------------------------------------------------------ select_wild_type


def test_example1_selects_wt_by_token(px: ModuleType) -> None:
    sel = px.select_wild_type(px.parse_record(fixture_xml(1).read_bytes()))
    assert _selected(sel) == [
        ("WT_1", "WT_1.d.zip", "token:WT"),
        ("WT_2", "WT_2.d.zip", "token:WT"),
        ("WT_3", "WT_3.d.zip", "token:WT"),
    ]
    assert _rejected(sel) == {
        "dLon_1.d.zip": "not_wild_type",
        "dLon_2.d.zip": "not_wild_type",
        "dLon_3.d.zip": "not_wild_type",
        "Report.sne": "not_raw",
    }


def test_example2_annotation_selects_and_labels_runs(px: ModuleType) -> None:
    sel = px.select_wild_type(px.parse_record(fixture_xml(2).read_bytes()))
    assert _selected(sel) == [
        ("WT_1", "S01.d.zip", "sdrf:genotype=wild type"),
        ("WT_2", "S02.d.zip", "sdrf:genotype=wild type"),
        ("WT_3", "S03.d.zip", "sdrf:genotype=wild type"),
    ]
    assert _rejected(sel) == {f"S0{i}.d.zip": "not_wild_type" for i in (4, 5, 6)}


def test_annotation_wins_over_the_file_name(px: ModuleType) -> None:
    """Example 2: a WT-named file annotated as a mutant is not wild type."""
    files = [
        FileSpec(f"S0{i}.d.zip", f"{BASE}/S0{i}.d.zip", genotype="wild type") for i in (1, 2, 3)
    ] + [FileSpec("WT_9.d.zip", f"{BASE}/WT_9.d.zip", genotype="ΔfliC")]
    sel = _sel(px, files)
    assert [s.file_name for s in sel.selected] == ["S01.d.zip", "S02.d.zip", "S03.d.zip"]
    assert _rejected(sel) == {"WT_9.d.zip": "not_wild_type"}


def test_example3_natural_sort_and_beyond_three(px: ModuleType) -> None:
    sel = px.select_wild_type(px.parse_record(fixture_xml(3).read_bytes()))
    assert [(s.run, s.file_name) for s in sel.selected] == [
        ("WT_1", "WT_1.d.zip"),
        ("WT_2", "WT_2.d.zip"),
        ("WT_3", "WT_3.d.zip"),
    ]
    assert _rejected(sel) == {"WT_10.d.zip": "beyond_three"}


def test_example4_fewer_than_three_raises_with_the_candidates(px: ModuleType) -> None:
    files = px.parse_record(px_xml(_wt("WT_2.d.zip", "WT_1.d.zip", "KO_1.d.zip")))
    with pytest.raises(px.SelectionError) as err:
        px.select_wild_type(files)
    assert issubclass(px.SelectionError, ValueError)
    assert "fewer than three wild-type runs: WT_1.d.zip, WT_2.d.zip" in str(err.value)


@pytest.mark.parametrize(
    "name",
    [
        "wildtype_1.d.zip",
        "Wild-Type_1.d.zip",
        "wild_type_1.d.zip",
        "EcN_wt_1.d.zip",
        "20260101_EcN_WT-1.d.zip",
    ],
)
def test_wild_type_name_tokens(px: ModuleType, name: str) -> None:
    sel = _sel(px, _wt(name, "WT_7.d.zip", "WT_8.d.zip"))
    assert name in [s.file_name for s in sel.selected], _rejected(sel)


@pytest.mark.parametrize(
    "name",
    [
        "SWTX_1.d.zip",  # WT not token-bounded
        "WTF_1.d.zip",
        "WT_KO_lon_1.d.zip",
        "WT_delta_lon_1.d.zip",
        "WT_ΔfliC_1.d.zip",
        "WT_mut_1.d.zip",
        "WT_pMUT-cured_1.d.zip",
        "WT_cured_1.d.zip",
        "WT_lon::kan_1.d.zip",
    ],
)
def test_names_that_are_not_wild_type(px: ModuleType, name: str) -> None:
    sel = _sel(px, _wt(name, "WT_7.d.zip", "WT_8.d.zip", "WT_9.d.zip"))
    assert [s.file_name for s in sel.selected] == ["WT_7.d.zip", "WT_8.d.zip", "WT_9.d.zip"]
    assert _rejected(sel) == {name: "not_wild_type"}


def test_multi_file_runs_share_a_run_label(px: ModuleType) -> None:
    """Rule 4: a .d stored as analysis.tdf + analysis.tdf_bin under a folder is one run."""
    files = []
    for n in (2, 1, 3):
        for part in ("analysis.tdf", "analysis.tdf_bin"):
            files.append(FileSpec(f"WT_{n}.d/{part}", f"{BASE}/WT_{n}.d/{part}", size=10))
    sel = _sel(px, files)
    runs = [s.run for s in sel.selected]
    assert runs == ["WT_1", "WT_1", "WT_2", "WT_2", "WT_3", "WT_3"], _selected(sel)
    assert not _rejected(sel)


# ------------------------------------------------------------------ properties


@settings(max_examples=60, deadline=None)
@given(st.sets(st.integers(min_value=1, max_value=300), min_size=3, max_size=12), st.randoms())
def test_natural_sort_picks_the_three_smallest_numbers(
    px: ModuleType, numbers: set[int], rnd
) -> None:  # type: ignore[no-untyped-def]
    order = sorted(numbers)
    rnd.shuffle(order)
    sel = _sel(px, _wt(*(f"WT_{n}.d.zip" for n in order)))
    smallest = sorted(numbers)[:3]
    assert [s.file_name for s in sel.selected] == [f"WT_{n}.d.zip" for n in smallest]
    assert [s.run for s in sel.selected] == ["WT_1", "WT_2", "WT_3"]
    rest = {f"WT_{n}.d.zip": "beyond_three" for n in sorted(numbers)[3:]}
    assert _rejected(sel) == rest


_SAFE = "bfghjnpqsxyz"  # letters that cannot spell WT or a modification token


@settings(max_examples=60, deadline=None)
@given(
    st.text(_SAFE, min_size=1, max_size=4),
    st.text(_SAFE, min_size=1, max_size=4),
    st.sampled_from(["WT", "wt", "Wt"]),
)
def test_wt_token_must_be_bounded(px: ModuleType, left: str, right: str, wt: str) -> None:
    glued = f"{left}{wt}{right}_1.d.zip"
    split = f"{left}_{wt}_{right}_1.d.zip"
    sel = _sel(px, _wt(glued, split, "WT_8.d.zip", "WT_9.d.zip"))
    chosen = [s.file_name for s in sel.selected]
    assert split in chosen and glued not in chosen, (chosen, _rejected(sel))
    assert _rejected(sel)[glued] == "not_wild_type"


@settings(max_examples=40, deadline=None)
@given(
    st.lists(
        st.tuples(st.sampled_from(["WT", "KO", "dLon", "wildtype"]), st.integers(1, 50)),
        min_size=3,
        max_size=10,
        unique=True,
    )
)
def test_every_file_is_selected_or_rejected_exactly_once(
    px: ModuleType, stems: list[tuple[str, int]]
) -> None:
    names = [f"{s}_{n}.d.zip" for s, n in stems]
    files = px.parse_record(px_xml(_wt(*names)))
    try:
        sel = px.select_wild_type(files)
    except px.SelectionError:
        wt = [n for n in names if n.startswith(("WT_", "wildtype_"))]
        assert len(wt) < 3
        return
    chosen = [s.file_name for s in sel.selected]
    rejected = _rejected(sel)
    assert len(chosen) == 3 and not set(chosen) & set(rejected)
    assert sorted([*chosen, *rejected]) == sorted(names)


# ------------------------------------------------------------------ CLI


def _cli(*args: str | Path) -> subprocess.CompletedProcess[str]:
    cmd = [sys.executable, "-I", str(PX_RECORD), *map(str, args)]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=False)


def test_cli_writes_the_example1_manifest(tmp_path: Path) -> None:
    out = tmp_path / "draft.tsv"
    res = _cli(fixture_xml(1), "--out", out)
    assert res.returncode == 0, res.stdout + res.stderr
    header, rows = read_tsv(out)
    assert tuple(header) == MANIFEST_COLUMNS
    specs = {s.name: s for s in example1_files()}
    assert [r["run"] for r in rows] == ["WT_1", "WT_2", "WT_3"]
    for r in rows:
        spec = specs[r["file_name"]]
        assert r["url_record"] == spec.uri
        assert r["url"] == spec.uri, "an https record URI is downloaded as is"
        assert r["bytes"] == str(spec.size)
        assert r["checksum_record"] == f"md5:{spec.md5}"
        assert r["sha256"] == "NA", "the draft has no sha256 yet (two-stage commit, Q5)"
        assert r["selected_because"] == "token:WT"
        assert r["sample"] != ""


def test_cli_writes_na_for_missing_values(tmp_path: Path) -> None:
    xml = tmp_path / "rec.xml"
    xml.write_bytes(px_xml([FileSpec(f"WT_{n}.d.zip", f"{BASE}/WT_{n}.d.zip") for n in (1, 2, 3)]))
    out = tmp_path / "draft.tsv"
    res = _cli(xml, "--out", out)
    assert res.returncode == 0, res.stdout + res.stderr
    _, rows = read_tsv(out)
    assert {r["bytes"] for r in rows} == {"NA"}
    assert {r["checksum_record"] for r in rows} == {"NA"}


def test_cli_fewer_than_three_exits_2_without_a_manifest(tmp_path: Path) -> None:
    xml = tmp_path / "rec.xml"
    xml.write_bytes(px_xml(_wt("WT_1.d.zip", "WT_2.d.zip", "dLon_1.d.zip")))
    out = tmp_path / "draft.tsv"
    res = _cli(xml, "--out", out)
    assert res.returncode == 2, res.stdout + res.stderr
    assert "fewer than three wild-type runs: WT_1.d.zip, WT_2.d.zip" in res.stdout + res.stderr
    assert not out.exists(), "example 4: no manifest written"


def test_cli_refuses_a_record_without_timstof(tmp_path: Path) -> None:
    xml = tmp_path / "rec.xml"
    xml.write_bytes(px_xml(example1_files(), instrument="Orbitrap Exploris 480"))
    out = tmp_path / "draft.tsv"
    res = _cli(xml, "--out", out)
    assert res.returncode == 2, res.stdout + res.stderr
    assert re.search(r"does not describe timsTOF data", res.stdout + res.stderr)
    assert not out.exists()


def test_module_is_stdlib_only() -> None:
    """Deliverable 1: stdlib only, so ``python -I`` runs it on any node."""
    assert PX_RECORD.is_file(), PX_RECORD
    imports = re.findall(r"^\s*(?:from|import)\s+([\w.]+)", PX_RECORD.read_text(), re.M)
    roots = {i.split(".")[0] for i in imports}
    non_std = sorted(r for r in roots if r not in sys.stdlib_module_names and r != "__future__")
    assert not non_std, non_std
