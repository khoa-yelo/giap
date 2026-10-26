"""Raw-file selection and acquisition mode from the data (SPEC §7)."""

from __future__ import annotations

import gzip
import sqlite3
from pathlib import Path

import pytest
from giap_campaign.pride import PrideFile
from giap_campaign.spectra import (
    choose_raw,
    classify_raw_file,
    mzml_acquisition,
    resolve_acquisition,
    select_runs,
    tdf_acquisition,
)
from hypothesis import given, settings
from hypothesis import strategies as st

DATA = Path(__file__).parent / "data" / "spectra"
BASE = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/03/PXD000001/"
GB = 1_000_000_000


def _pf(name: str, size: int | None = 1000, url: str | None = "default") -> PrideFile:
    return PrideFile(
        name=name,
        category="RAW",
        size_bytes=size,
        url=BASE + name if url == "default" else url,
        checksum=None,
    )


# ---------------------------------------------------------------- classify_raw_file


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("Ecoli_WT_1.raw", "thermo_raw"),
        ("ECOLI_WT_1.RAW", "thermo_raw"),
        ("sub/dir/Run1.Raw", "thermo_raw"),
        ("run.d.zip", "bruker_d"),
        ("RUN.D.ZIP", "bruker_d"),
        ("run.d.tar", "bruker_d"),
        ("run.d.tar.gz", "bruker_d"),
        ("Run_Slot1-1_1_123.d.TAR.GZ", "bruker_d"),
        ("run.d.tgz", "bruker_d"),
        ("run.mzML", "mzml"),
        ("run.mzml", "mzml"),
        ("RUN.MZML", "mzml"),
        ("run.mzML.gz", "mzml"),
        ("run.MZML.GZ", "mzml"),
        ("run.wiff", "sciex_wiff"),
        ("RUN.WIFF", "sciex_wiff"),
        ("run.wiff2", "sciex_wiff"),
        ("run.wiff.scan", None),
        ("RUN.WIFF.SCAN", None),
        ("run.mgf", None),
        ("run.mzXML", None),
        ("run.raw.zip", None),
        ("run.mzML.zip", None),
        ("run.zip", None),
        ("run.tar.gz", None),
        ("run.mzid", None),
        ("run.mzid.gz", None),
        ("msms.txt", None),
        ("proteinGroups.txt", None),
        ("report.tsv", None),
        ("results.sage.tsv", None),
        ("README", None),
    ],
)
def test_classify_raw_file(name: str, expected: str | None) -> None:
    assert classify_raw_file(name) == expected


@settings(max_examples=100, deadline=None)
@given(
    stem=st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789_-", min_size=1, max_size=12),
    ext=st.sampled_from(
        [
            *(".raw", ".d.zip", ".d.tar", ".d.tar.gz", ".d.tgz", ".mzml", ".mzml.gz"),
            *(".wiff", ".wiff2", ".wiff.scan", ".mgf", ".txt", ".mzid", ".zip"),
        ]
    ),
    upper=st.lists(st.booleans(), min_size=30, max_size=30),
)
def test_classify_raw_file_is_case_insensitive(stem: str, ext: str, upper: list[bool]) -> None:
    name = stem + ext
    mixed = "".join(c.upper() if u else c for c, u in zip(name, upper, strict=False))
    assert len(mixed) == len(name)
    assert classify_raw_file(mixed) == classify_raw_file(name.lower())
    assert classify_raw_file("PXD000001/" + mixed) == classify_raw_file(name)


# ---------------------------------------------------------------- choose_raw


def test_choose_raw_most_files_wins_and_wiff_is_never_chosen() -> None:
    files = [
        _pf("c.raw"),
        _pf("a.raw"),
        _pf("b.raw"),
        _pf("x.d.zip"),
        _pf("y.d.zip"),
        _pf("m.mzML"),
        *[_pf(f"w{i}.wiff") for i in range(5)],
        _pf("w0.wiff.scan"),
        _pf("search.mzid"),
        _pf("peaks.mgf"),
    ]
    fmt, chosen, seen = choose_raw(files)
    assert fmt == "thermo_raw"
    assert [f.name for f in chosen] == ["a.raw", "b.raw", "c.raw"]
    assert set(seen) == {"thermo_raw", "bruker_d", "mzml", "sciex_wiff"}


def test_choose_raw_returns_pride_files_unchanged() -> None:
    files = [_pf("b.raw", 7), _pf("a.raw", 5)]
    _, chosen, _ = choose_raw(files)
    assert chosen == [files[1], files[0]]


@pytest.mark.parametrize(
    ("names", "expected"),
    [
        (["a.raw", "b.raw", "x.d.zip", "y.d.tar.gz"], "thermo_raw"),
        (["x.d.zip", "y.d.tgz", "a.mzML", "b.mzML.gz"], "bruker_d"),
        (["a.mzML", "b.raw", "x.d.zip"], "thermo_raw"),
        (["a.mzML", "x.d.zip"], "bruker_d"),
        (["a.mzML", "b.mzML", "x.d.zip"], "mzml"),
        (["a.mzML.gz"], "mzml"),
    ],
)
def test_choose_raw_tie_order(names: list[str], expected: str) -> None:
    fmt, chosen, _ = choose_raw([_pf(n) for n in names])
    assert fmt == expected
    assert all(classify_raw_file(f.name) == expected for f in chosen)


def test_choose_raw_ignores_files_without_url() -> None:
    files = [
        _pf("a.raw", url=None),
        _pf("b.raw", url=None),
        _pf("c.raw", url=None),
        _pf("d.raw"),
        _pf("m1.mzML"),
        _pf("m2.mzML"),
    ]
    fmt, chosen, _ = choose_raw(files)
    assert fmt == "mzml"
    assert [f.name for f in chosen] == ["m1.mzML", "m2.mzML"]


def test_choose_raw_chosen_files_exclude_those_without_url() -> None:
    files = [_pf("b.raw"), _pf("a.raw", url=None), _pf("c.raw")]
    fmt, chosen, _ = choose_raw(files)
    assert fmt == "thermo_raw"
    assert [f.name for f in chosen] == ["b.raw", "c.raw"]


def test_choose_raw_only_wiff_gives_none() -> None:
    fmt, chosen, seen = choose_raw([_pf("a.wiff"), _pf("a.wiff.scan"), _pf("b.wiff2")])
    assert fmt is None
    assert list(chosen) == []
    assert set(seen) == {"sciex_wiff"}


def test_choose_raw_no_raw_files() -> None:
    fmt, chosen, seen = choose_raw([_pf("msms.txt"), _pf("peaks.mgf")])
    assert fmt is None
    assert list(chosen) == []
    assert set(seen) == set()
    fmt, chosen, seen = choose_raw([])
    assert (fmt, list(chosen), set(seen)) == (None, [], set())


def test_choose_raw_sorts_by_name() -> None:
    names = ["run10.raw", "run02.raw", "run1.raw", "Run3.raw"]
    _, chosen, _ = choose_raw([_pf(n) for n in names])
    assert [f.name for f in chosen] == sorted(names)


# ---------------------------------------------------------------- select_runs


def _names(files: list[PrideFile]) -> list[str]:
    return [f.name for f in files]


def test_select_runs_evenly_spaced_n10_k4() -> None:
    files = [_pf(f"run{i:02d}.raw") for i in range(10)]
    got = select_runs(files[::-1], max_runs=4)
    assert _names(got) == ["run00.raw", "run03.raw", "run06.raw", "run09.raw"]


def test_select_runs_k1_takes_first() -> None:
    files = [_pf(f"run{i:02d}.raw") for i in range(10)]
    assert _names(select_runs(files, max_runs=1)) == ["run00.raw"]


def test_select_runs_fewer_files_than_max() -> None:
    files = [_pf("c.raw"), _pf("a.raw"), _pf("b.raw")]
    assert _names(select_runs(files)) == ["a.raw", "b.raw", "c.raw"]


def test_select_runs_default_max_runs_is_12() -> None:
    files = [_pf(f"run{i:02d}.raw") for i in range(30)]
    assert len(select_runs(files)) == 12


def test_select_runs_n_equals_k_takes_all() -> None:
    files = [_pf(f"run{i:02d}.raw") for i in range(12)]
    assert _names(select_runs(files, max_runs=12)) == [f"run{i:02d}.raw" for i in range(12)]


def test_select_runs_size_cap_lowers_k() -> None:
    # n = 7: k = 4 -> 0,2,4,6 (160 GB); k = 3 -> 0,3,6 (120 GB); k = 2 -> 0,6 (80 GB)
    files = [_pf(f"run{i}.raw", 40 * GB) for i in range(7)]
    assert _names(select_runs(files, max_runs=4, max_gb=100.0)) == ["run0.raw", "run6.raw"]
    assert _names(select_runs(files, max_runs=4, max_gb=120.0)) == [
        "run0.raw",
        "run3.raw",
        "run6.raw",
    ]
    assert _names(select_runs(files, max_runs=4, max_gb=1000.0)) == [
        "run0.raw",
        "run2.raw",
        "run4.raw",
        "run6.raw",
    ]


def test_select_runs_default_max_gb_is_100() -> None:
    files = [_pf(f"run{i}.raw", 40 * GB) for i in range(7)]
    assert _names(select_runs(files, max_runs=4)) == ["run0.raw", "run6.raw"]


def test_select_runs_unknown_sizes_count_zero() -> None:
    files = [_pf(f"run{i}.raw", None) for i in range(7)]
    got = select_runs(files, max_runs=4, max_gb=0.001)
    assert _names(got) == ["run0.raw", "run2.raw", "run4.raw", "run6.raw"]


def test_select_runs_mixed_unknown_sizes() -> None:
    sizes = [60 * GB, None, None, None, None, None, 60 * GB]
    files = [_pf(f"run{i}.raw", s) for i, s in enumerate(sizes)]
    # k = 4 -> 0,2,4,6 = 120 GB > 100; k = 3 -> 0,3,6 = 120 GB; k = 2 -> 0,6 = 120 GB; k = 1 -> 60
    assert _names(select_runs(files, max_runs=4, max_gb=100.0)) == ["run0.raw"]


def test_select_runs_always_at_least_one() -> None:
    assert _names(select_runs([_pf("huge.raw", 500 * GB)], max_gb=100.0)) == ["huge.raw"]
    files = [_pf(f"run{i}.raw", 500 * GB) for i in range(3)]
    assert _names(select_runs(files, max_gb=1.0)) == ["run0.raw"]


def test_select_runs_is_deterministic() -> None:
    files = [_pf(f"run{i:02d}.raw", (i % 3) * 20 * GB) for i in range(25)]
    a = select_runs(files, max_runs=9, max_gb=150.0)
    b = select_runs(list(reversed(files)), max_runs=9, max_gb=150.0)
    assert a == b


def _reference_select(files: list[PrideFile], max_runs: int, max_gb: float) -> list[str]:
    """SPEC §7 literally (Python's built-in ``round``)."""
    srt = sorted(files, key=lambda f: f.name)
    n = len(srt)
    k = min(max_runs, n)
    while True:
        idx = [0] if k == 1 else [round(i * (n - 1) / (k - 1)) for i in range(k)]
        chosen = [srt[i] for i in idx]
        if k == 1 or sum(f.size_bytes or 0 for f in chosen) / GB <= max_gb:
            return [f.name for f in chosen]
        k -= 1


_NAMES = st.lists(
    st.text(alphabet="abcdefghij0123456789_", min_size=1, max_size=8),
    min_size=1,
    max_size=40,
    unique=True,
)
_SIZE = st.one_of(st.none(), st.integers(min_value=0, max_value=60 * GB))


@settings(max_examples=200, deadline=None)
@given(
    data=st.data(),
    names=_NAMES,
    max_runs=st.integers(min_value=1, max_value=15),
    max_gb=st.floats(min_value=0.0, max_value=300.0, allow_nan=False),
)
def test_select_runs_properties(
    data: st.DataObject, names: list[str], max_runs: int, max_gb: float
) -> None:
    sizes = data.draw(st.lists(_SIZE, min_size=len(names), max_size=len(names)))
    files = [_pf(n + ".raw", s) for n, s in zip(names, sizes, strict=True)]
    srt = sorted(files, key=lambda f: f.name)
    got = select_runs(files, max_runs=max_runs, max_gb=max_gb)

    # a subsequence of the name-sorted input
    positions = [srt.index(f) for f in got]
    assert positions == sorted(set(positions))
    assert 1 <= len(got) <= min(max_runs, len(files))
    # the first sorted file is always taken; the last whenever two or more are
    assert got[0] == srt[0]
    if len(got) >= 2:
        assert got[-1] == srt[-1]
    total = sum(f.size_bytes or 0 for f in got) / GB
    assert len(got) == 1 or total <= max_gb
    # deterministic and independent of the input order
    shuffled = data.draw(st.permutations(files))
    assert select_runs(list(shuffled), max_runs=max_runs, max_gb=max_gb) == got
    assert [f.name for f in got] == _reference_select(files, max_runs, max_gb)


# ---------------------------------------------------------------- tdf_acquisition


def _tdf_file(tmp_path: Path, types: list[int]) -> Path:
    d = tmp_path / "run.d"
    d.mkdir()
    path = d / "analysis.tdf"
    con = sqlite3.connect(path)
    con.execute(
        "CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL NOT NULL, Polarity CHAR(1), "
        "ScanMode INTEGER, MsMsType INTEGER NOT NULL, NumScans INTEGER, NumPeaks INTEGER)"
    )
    con.execute("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
    con.execute("INSERT INTO GlobalMetadata VALUES ('AcquisitionSoftware', 'timsControl')")
    con.executemany(
        "INSERT INTO Frames VALUES (?, ?, '+', 9, ?, 927, 1000)",
        [(i + 1, 0.1 * i, t) for i, t in enumerate(types)],
    )
    con.commit()
    con.close()
    return path


@pytest.mark.parametrize(
    ("types", "expected"),
    [
        ([0, 9, 9, 9, 0, 9, 9, 9], "DIA"),
        ([9, 9], "DIA"),
        ([0, 8, 8, 8, 0, 8], "DDA"),
        ([8], "DDA"),
        ([0, 8, 9, 0], "unknown"),
        ([8, 9], "unknown"),
        ([0, 0, 0], "unknown"),
        ([], "unknown"),
    ],
)
def test_tdf_acquisition(tmp_path: Path, types: list[int], expected: str) -> None:
    path = _tdf_file(tmp_path, types)
    assert tdf_acquisition(path) == expected
    assert tdf_acquisition(str(path)) == expected


def test_tdf_acquisition_is_read_only(tmp_path: Path) -> None:
    path = _tdf_file(tmp_path, [0, 9, 9])
    before = path.read_bytes()
    listing = sorted(p.name for p in path.parent.iterdir())
    assert tdf_acquisition(path) == "DIA"
    assert path.read_bytes() == before
    assert sorted(p.name for p in path.parent.iterdir()) == listing


# ---------------------------------------------------------------- mzml_acquisition

_MZML_HEAD = """<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0" id="t">
  <cvList count="1">
    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" \
version="4.1.79" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>
  </cvList>
  <run id="t" defaultInstrumentConfigurationRef="IC1">
    <spectrumList count="{n}" defaultDataProcessingRef="dp">
"""
_MZML_TAIL = """    </spectrumList>
  </run>
</mzML>
"""


def _spectrum(i: int, level: int, lower: float | None = None, upper: float | None = None) -> str:
    out = [
        f'      <spectrum index="{i}" id="controllerType=0 controllerNumber=1 scan={i + 1}" '
        'defaultArrayLength="0">',
        f'        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{level}"/>',
        '        <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>',
        '        <scanList count="1"><scan><cvParam cvRef="MS" accession="MS:1000016" '
        f'name="scan start time" value="{0.01 * i:.3f}" unitCvRef="UO" '
        'unitAccession="UO:0000031" unitName="minute"/></scan></scanList>',
    ]
    if level == 2:
        out.append('        <precursorList count="1"><precursor><isolationWindow>')
        out.append(
            '          <cvParam cvRef="MS" accession="MS:1000827" '
            'name="isolation window target m/z" value="500.0"/>'
        )
        if lower is not None:
            out.append(
                '          <cvParam cvRef="MS" accession="MS:1000828" '
                f'name="isolation window lower offset" value="{lower}"/>'
            )
        if upper is not None:
            out.append(
                '          <cvParam cvRef="MS" accession="MS:1000829" '
                f'name="isolation window upper offset" value="{upper}"/>'
            )
        out.append(
            '        </isolationWindow><selectedIonList count="1"><selectedIon>'
            '<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="500.0"/>'
            "</selectedIon></selectedIonList></precursor></precursorList>"
        )
    out.append('        <binaryDataArrayList count="0"/>')
    out.append("      </spectrum>")
    return "\n".join(out) + "\n"


def _mzml(spectra: list[tuple[int, float | None, float | None]], *, tail: bool = True) -> str:
    body = "".join(_spectrum(i, lv, lo, hi) for i, (lv, lo, hi) in enumerate(spectra))
    return _MZML_HEAD.replace("{n}", str(len(spectra))) + body + (_MZML_TAIL if tail else "")


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    if name.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


MS1 = (1, None, None)


def _ms2(lower: float | None, upper: float | None) -> tuple[int, float | None, float | None]:
    return (2, lower, upper)


def test_mzml_acquisition_dia_wide_windows(tmp_path: Path) -> None:
    spectra = [MS1, *[_ms2(12.5, 12.5)] * 6, MS1, *[_ms2(12.5, 12.5)] * 6]
    assert mzml_acquisition(_write(tmp_path, "dia.mzML", _mzml(spectra))) == "DIA"


def test_mzml_acquisition_dda_narrow_windows(tmp_path: Path) -> None:
    spectra = [MS1, *[_ms2(0.8, 0.8)] * 10, MS1, *[_ms2(0.8, 0.8)] * 10]
    assert mzml_acquisition(_write(tmp_path, "dda.mzML", _mzml(spectra))) == "DDA"


def test_mzml_acquisition_width_is_lower_plus_upper(tmp_path: Path) -> None:
    # lower 1.0 + upper 4.5 = 5.5 Th (DIA); twice either offset alone would say otherwise
    spectra = [MS1, *[_ms2(1.0, 4.5)] * 5]
    assert mzml_acquisition(_write(tmp_path, "asym.mzML", _mzml(spectra))) == "DIA"
    spectra = [MS1, *[_ms2(4.0, 0.5)] * 5]  # 4.5 Th: DDA
    assert mzml_acquisition(_write(tmp_path, "asym2.mzML", _mzml(spectra))) == "DDA"


def test_mzml_acquisition_threshold_5_is_dia(tmp_path: Path) -> None:
    spectra = [MS1, *[_ms2(2.5, 2.5)] * 3]
    assert mzml_acquisition(_write(tmp_path, "five.mzML", _mzml(spectra))) == "DIA"


def test_mzml_acquisition_uses_the_median(tmp_path: Path) -> None:
    spectra = [MS1, _ms2(0.8, 0.8), _ms2(0.8, 0.8), _ms2(0.8, 0.8), _ms2(200, 200), _ms2(200, 200)]
    assert mzml_acquisition(_write(tmp_path, "med.mzML", _mzml(spectra))) == "DDA"
    spectra = [MS1, _ms2(0.8, 0.8), _ms2(0.8, 0.8), _ms2(10, 10), _ms2(10, 10), _ms2(10, 10)]
    assert mzml_acquisition(_write(tmp_path, "med2.mzML", _mzml(spectra))) == "DIA"


def test_mzml_acquisition_no_ms2_is_unknown(tmp_path: Path) -> None:
    spectra = [MS1] * 5
    assert mzml_acquisition(_write(tmp_path, "ms1.mzML", _mzml(spectra))) == "unknown"


def test_mzml_acquisition_no_offsets_is_unknown(tmp_path: Path) -> None:
    spectra = [MS1, *[_ms2(None, None)] * 5]
    assert mzml_acquisition(_write(tmp_path, "nooff.mzML", _mzml(spectra))) == "unknown"


def test_mzml_acquisition_empty_spectrum_list_is_unknown(tmp_path: Path) -> None:
    assert mzml_acquisition(_write(tmp_path, "empty.mzML", _mzml([]))) == "unknown"


def test_mzml_acquisition_gzip(tmp_path: Path) -> None:
    dia = [MS1, *[_ms2(12.5, 12.5)] * 4]
    dda = [MS1, *[_ms2(0.7, 0.7)] * 4]
    assert mzml_acquisition(_write(tmp_path, "dia.mzML.gz", _mzml(dia))) == "DIA"
    assert mzml_acquisition(_write(tmp_path, "dda.mzML.gz", _mzml(dda))) == "DDA"


def test_mzml_acquisition_indexed_thermorawfileparser_fixture() -> None:
    """An indexedmzML wrapper as ThermoRawFileParser writes it (gzip, 25 Th windows)."""
    assert mzml_acquisition(DATA / "dia_small.mzML.gz") == "DIA"
    assert mzml_acquisition(str(DATA / "dia_small.mzML.gz")) == "DIA"


def test_mzml_acquisition_reads_only_the_first_max_spectra_ms2(tmp_path: Path) -> None:
    # 3 narrow MS2 first, then 10 wide ones
    spectra = [MS1, *[_ms2(0.8, 0.8)] * 3, MS1, *[_ms2(12.5, 12.5)] * 10]
    path = _write(tmp_path, "mix.mzML", _mzml(spectra))
    assert mzml_acquisition(path, max_spectra=3) == "DDA"
    assert mzml_acquisition(path) == "DIA"
    assert mzml_acquisition(path, max_spectra=200) == "DIA"


def test_mzml_acquisition_ms1_spectra_do_not_count_towards_max(tmp_path: Path) -> None:
    spectra = [*[MS1] * 5, *[_ms2(12.5, 12.5)] * 3]
    assert mzml_acquisition(_write(tmp_path, "ms1first.mzML", _mzml(spectra)), max_spectra=3) == (
        "DIA"
    )


def test_mzml_acquisition_streams_and_stops_early(tmp_path: Path) -> None:
    """The file is cut off after the third MS2 spectrum: a streaming reader never reaches it."""
    text = _mzml([MS1, *[_ms2(0.8, 0.8)] * 3], tail=False)
    text += '      <spectrum index="4" id="controllerType=0 controllerNumber=1 scan=5" defaultArr'
    path = _write(tmp_path, "truncated.mzML", text)
    assert mzml_acquisition(path, max_spectra=3) == "DDA"


# ---------------------------------------------------------------- resolve_acquisition


@pytest.mark.parametrize(
    ("data_modes", "metadata", "mode", "note"),
    [
        (["DIA", "DIA", "DIA"], "DIA", "DIA", None),
        (["DIA"], "unknown", "DIA", None),
        (["DDA", "DDA"], "DDA", "DDA", None),
        (["DDA", "DDA"], "unknown", "DDA", None),
        (["DIA", "DIA"], "DDA", "DIA", "metadata_disagrees"),
        (["DDA"], "DIA", "DDA", "metadata_disagrees"),
        (["DIA", "DDA"], "DIA", "mixed", "mixed_acquisition"),
        (["DDA", "DIA", "DDA"], "unknown", "mixed", "mixed_acquisition"),
        (["unknown", "unknown"], "DIA", "DIA", "from_metadata"),
        (["unknown"], "DDA", "DDA", "from_metadata"),
        (["unknown", "unknown"], "unknown", "unknown", "undetermined"),
        (["unknown"], "unknown", "unknown", "undetermined"),
    ],
)
def test_resolve_acquisition(
    data_modes: list[str], metadata: str, mode: str, note: str | None
) -> None:
    got_mode, got_note = resolve_acquisition(data_modes, metadata)
    assert got_mode == mode
    if note is None:
        assert not got_note  # no note: "" (or None)
    else:
        assert got_note == note


def test_resolve_acquisition_accepts_any_sequence() -> None:
    assert resolve_acquisition(("DIA", "DIA"), "DIA")[0] == "DIA"
