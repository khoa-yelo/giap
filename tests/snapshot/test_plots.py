"""M9-25 work item 2: one pytest-mpl snapshot per plot instance on both frozen fixtures
(docs/specs/M9-25-snapshot-baselines.md; GIAP_SPEC.md §10.1 "tolerance RMS <= 2"; index
decision 22).

``pytest tests/snapshot --mpl`` compares each figure with
``tests/snapshot/baseline/<fixture>-<stem>.png`` (RMS <= 2, 100 dpi); without ``--mpl`` the
figures are only rendered (examples 1 and 4). Cases, fixtures and naming:
``tests/helpers/snapshot_fixture.py``. A fixture that is not committed (or fails to render)
gives one failing ``<fixture>-setup`` case, never an empty parametrisation.

The other tests here need no ``--mpl``: the fixtures are rich (every plot non-empty on
``sim_plots``) and the empty-state tree renders every plot; every case has a baseline and every
baseline a case; size limits (spec work items 1-2, R11).
"""

from __future__ import annotations

import gzip
import importlib.util
import struct
from pathlib import Path

import pytest
from helpers import snapshot_fixture as sf

CASES = sf.cases()


def _param(case: sf.Case) -> object:
    if case.error is not None:
        return pytest.param(case, id=case.id)
    mark = pytest.mark.mpl_image_compare(
        filename=case.filename,
        tolerance=sf.TOLERANCE,
        savefig_kwargs={"dpi": sf.SNAPSHOT_DPI},
        style="default",
    )
    return pytest.param(case, id=case.id, marks=mark)


@pytest.mark.parametrize("case", [_param(c) for c in CASES])
def test_plot_snapshot(case: sf.Case) -> object:
    if case.error is not None:
        pytest.fail(case.error)
    from giap.viz import style

    style.apply()
    return sf.figure(case)


# ---------------------------------------------------------------- the fixtures


def _require(fixture: str) -> None:
    why = sf.missing(fixture)
    if why is not None:
        pytest.fail(why)


@pytest.mark.parametrize("fixture", sf.FIXTURES)
def test_fixture_is_committed(fixture: str) -> None:
    _require(fixture)
    runs = [p for p in (sf.fixture_cmp(fixture) / "runs").iterdir() if p.is_dir()]
    assert len(runs) >= 2, "CMP/runs needs two or more runs"
    for name in sf.COMPARE_TABLES:
        assert (sf.fixture_cmp(fixture) / "tables" / f"{name}.tsv").is_file(), name


def test_make_snapshot_fixture_script_exists() -> None:
    path = sf.REPO / "scripts" / "dev" / "make_snapshot_fixture.py"
    assert path.is_file(), "scripts/dev/make_snapshot_fixture.py (spec work item 1) is missing"
    spec = importlib.util.spec_from_file_location("make_snapshot_fixture", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert callable(getattr(mod, "main", None))
    assert mod.__doc__ and "giap run" in mod.__doc__


def test_rich_fixture_every_plot_draws_something() -> None:
    from giap.viz import base, checklist

    _require(sf.RICH)
    results = sf.styled(lambda: sf.raw_results(sf.RICH))
    assert set(results) == set(base.registered())
    assert len(results) == 16
    skipped = {p: r.reason for p, r in results.items() if isinstance(r, base.Skip)}
    assert skipped == {}
    empty = sorted(p for p, r in results.items() if any(checklist.is_empty_state(i) for i in r))
    assert empty == []
    for pid in ("V03", "V13"):
        loci = [i for i in results[pid] if sf.is_locus(pid, i.suffix)]
        assert len(loci) >= sf.LOCUS_KEEP, f"{pid}: {len(loci)} locus instances"


def _table(fixture: str, key: str) -> object:
    from giap.io.run_view import RunView

    frame = RunView.from_outdir(sf.fixture_out(fixture)).table(key)
    assert frame is not None, f"{fixture}: {key} absent"
    return frame


def test_rich_fixture_plants_the_spec_features() -> None:
    """Spec work item 1: planted island, prophage and BGC with low detection, flagged operon
    pairs, ALT sites, novel ORFs, function and determinants tables."""
    _require(sf.RICH)
    assert len(_table(sf.RICH, "islands")) >= 1
    regions = _table(sf.RICH, "regions_stats")
    for kind in ("prophage", "bgc"):
        rows = regions[regions["region_type"] == kind]
        assert len(rows) >= 1, kind
        assert (rows["o_minus_e"].astype(float) < 0).any(), f"no {kind} with low detection"
    assert len(_table(sf.RICH, "operon_flags")) >= 1
    variants = _table(sf.RICH, "variants")
    assert (variants["call"] == "ALT").any()
    novel = _table(sf.RICH, "novel_orfs")
    assert int(novel["called"].fillna(False).astype(bool).sum()) >= sf.LOCUS_KEEP
    n_enrich = len(_table(sf.RICH, "enrichment_gsea")) + len(
        _table(sf.RICH, "enrichment_detection")
    )
    assert n_enrich >= 1
    assert len(_table(sf.RICH, "determinants_model")) >= 1
    assert len(_table(sf.RICH, "genes_tsv")) >= 100


#: plots whose primary table is header-only on ``sim_empty`` and so draw ``base.empty_state``
EMPTY_STATE_PIDS = ("V02", "V03", "V05", "V07", "V08", "V10", "V11", "V12", "V13", "V14", "V15",
                    "V16")  # fmt: skip


def test_empty_fixture_is_header_only() -> None:
    _require(sf.EMPTY)
    roots = [sf.fixture_out(sf.EMPTY), sf.fixture_cmp(sf.EMPTY)]
    tsvs = [p for r in roots for p in sorted(r.rglob("*.tsv*")) if "support" not in p.parts]
    assert tsvs
    for path in tsvs:
        raw = path.read_bytes()
        text = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode("utf-8")
        lines = [ln for ln in text.splitlines() if ln.strip()]
        assert len(lines) == 1, f"{path.relative_to(sf.REPO)}: {len(lines) - 1} data rows"


def test_empty_fixture_renders_every_plot_in_its_empty_state() -> None:
    from giap.viz import base, checklist

    _require(sf.EMPTY)
    results = sf.styled(lambda: sf.raw_results(sf.EMPTY))
    skipped = {p: r.reason for p, r in results.items() if isinstance(r, base.Skip)}
    assert skipped == {}
    not_empty = [p for p in EMPTY_STATE_PIDS if not all(map(checklist.is_empty_state, results[p]))]
    assert not_empty == []


# ---------------------------------------------------------------- baselines


def test_every_case_has_a_baseline() -> None:
    errors = [c.error for c in CASES if c.error is not None]
    assert errors == []
    absent = [c.filename for c in CASES if not (sf.BASELINE_DIR / c.filename).is_file()]
    assert absent == [], f"{len(absent)} baselines missing under tests/snapshot/baseline"


def test_every_baseline_has_a_case() -> None:
    assert sf.BASELINE_DIR.is_dir(), "tests/snapshot/baseline is not committed"
    known = {c.filename for c in CASES}
    orphans = sorted(p.name for p in sf.BASELINE_DIR.iterdir() if p.name not in known)
    assert orphans == []


def test_baselines_cover_all_16_plots_on_both_fixtures() -> None:
    from giap.steps.plots import load_plot_modules
    from giap.viz import base

    load_plot_modules()
    pids = set(base.registered())
    for fixture in sf.FIXTURES:
        have = {c.pid for c in CASES if c.fixture == fixture and c.error is None}
        assert have == pids, f"{fixture}: missing {sorted(pids - have)}"
        for pid, prefix in (("V03", "V03_locus_"), ("V13", "V13_novel_orfs_locus_")):
            stems = [c.stem for c in CASES if c.fixture == fixture and c.pid == pid]
            loci = [s for s in stems if s.startswith(prefix)]
            # "V03 and V13 locus plots: the first two instances only" (spec work item 2)
            assert len(loci) == (sf.LOCUS_KEEP if fixture == sf.RICH else 0), stems


def _png_dpi(path: Path) -> float | None:
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n", path.name
    pos = 8
    while pos < len(data):
        (length,) = struct.unpack(">I", data[pos : pos + 4])
        kind = data[pos + 4 : pos + 8]
        if kind == b"pHYs":
            x, _, unit = struct.unpack(">IIB", data[pos + 8 : pos + 17])
            return x * 0.0254 if unit == 1 else None
        pos += 12 + length
    return None


def test_baselines_are_small_100_dpi_pngs() -> None:
    assert sf.BASELINE_DIR.is_dir(), "tests/snapshot/baseline is not committed"
    pngs = sorted(sf.BASELINE_DIR.iterdir())
    assert pngs
    for path in pngs:
        assert path.suffix == ".png", path.name
        assert path.stat().st_size < sf.MAX_BASELINE_BYTES, path.name
        dpi = _png_dpi(path)
        assert dpi is not None and abs(dpi - sf.SNAPSHOT_DPI) < 0.5, (path.name, dpi)


def test_snapshot_fixture_within_r11_limits() -> None:
    assert sf.SNAPSHOT.is_dir(), "tests/data/snapshot is not committed"
    files = [p for p in sf.SNAPSHOT.rglob("*") if p.is_file()]
    total = sum(p.stat().st_size for p in files)
    assert total < sf.MAX_FIXTURE_BYTES, f"tests/data/snapshot is {total} bytes"
    assert [p.name for p in files if p.stat().st_size > sf.MAX_FILE_BYTES] == []
    # R11 "compress text fixtures with gzip": no plain TSV above 16 kB
    big_plain = [
        str(p.relative_to(sf.SNAPSHOT))
        for p in files
        if ".tsv" in p.name and p.stat().st_size > 16_000 and p.read_bytes()[:2] != b"\x1f\x8b"
    ]
    assert big_plain == []


def test_snapshot_fixture_has_provenance() -> None:
    prov = sf.SNAPSHOT / "PROVENANCE.md"
    assert prov.is_file(), "tests/data/snapshot/PROVENANCE.md is missing"
    text = prov.read_text(encoding="utf-8")
    assert "make_snapshot_fixture.py" in text
    assert "giapsim" in text
