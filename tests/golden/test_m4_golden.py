"""The committed M4 tables of the 50-kb fixture (GIAP_SPEC.md §11.4; M4-20 example 2).

Fixture: ``build_m4_fixture(work, seed=GOLDEN_SEED, size="sim50k")`` (the M0 50-kb genome of
``simulate(SimConfig(), GOLDEN_SEED)`` as an M4-07 layout, plus the synthetic prophage block and
AMR set), then ``run_gene_table`` -> ``run_detectability`` -> ``run_context_stats`` ->
``run_sensitivity`` (``genes.tsv`` is written by ``context_stats``, M5-09 decision 3) with the
golden config. The four ``tables/*.tsv`` must equal ``tests/data/golden/sim50k/m4/`` byte for
byte, with one exception.

**Lead decision 9: the four float columns of ``regions_stats.tsv`` named in ``NUMERIC_COLUMNS``
are compared numerically, everything else in that file byte for byte.**
``giap.context.regions.residualise`` (ADR-0018) solves its least squares with ``np.linalg.lstsq``,
a LAPACK SVD that is not bit-portable: one fixed 51x4 design gives ``beta[0]`` =
``0x1.85e19f7a40364p-6`` on iron and ``...377p-6`` on wheat, the three CI legs disagree with each
other, and thread pinning does not help (probe 1863133), so no recording machine can make the
matrix green. Everything the claim is really about -- every count, p-value, q-value, flag, the
other three tables and TE7 (two runs on one machine) -- stays byte-identical. Writer-level
rounding is deferred to the M6 writers (decision 10).

``$GIAP_DB_DIR`` is a test database with the fixture ``cog-20.def.tab`` (tests/helpers/
gene_table_fixture.py), never the installed one: otherwise ``cog_category`` -- and so the golden
genes table -- would depend on whether ``giap setup --cog`` has run on the machine.

Detection on 50 kb: the genome has 51 genes, 30 of them with peptide evidence, and the smallest
attainable picked gene q-value is 1/29 = 0.0345. At the default ``stats.alpha.identification`` =
0.01 *nothing* is ever detected and every downstream number is a zero
(``test_default_alpha_detects_nothing_on_50_kb`` pins that known limitation, the one the MVP-01
run 1861972 hit on sim100k); the golden therefore runs at ``alpha.identification`` = 0.05, where
29 of the 51 genes are detected (lead decision 1).

**M4-21: a second case, ``sim100k_real``,** runs the same four steps on the outputs of the real
M2/M3 steps instead of the giapsim contract fixture: giapsim SIM100K with its precomputed layers
and proteomics -> ``giap run --executor local`` (``parse_bakta``, TUs, the regions layer,
``context_layers``, the M3 mapping steps; tests/helpers/m4_real.py) -> the M4 chain on that
``OUT/work``. Its files live in ``tests/data/golden/sim100k_real/m4/``; ``sim50k/`` stays the
M4-only regression. The spec calls the directory ``sim50k_real``, but the 50-kb giapsim genome
plants no region at all (its precomputed ``regions.tsv`` is header-only), so its
``regions_stats.tsv`` would be empty; SIM100K plants a prophage, a BGC, a CRISPR array and a
plasmid and is the genome of every sim100k recording and of the MVP-01 run. Same golden config
(alpha 0.05: SIM100K's picked-q floor is 1/85 = 0.0118 > 0.01), same decision 9 tolerance.

Regenerating (on purpose, reviewed in the PR diff as data, with the statistician's note):
``GIAP_UPDATE_GOLDEN=1 pytest tests/golden/test_m4_golden.py`` writes the files and skips.
"""

from __future__ import annotations

import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc

GOLDEN_ROOT = Path(__file__).resolve().parents[1] / "data" / "golden"
GOLDEN = GOLDEN_ROOT / "sim50k" / "m4"
#: case -> committed directory; ``sim100k_real`` goes through the real M2/M3 path (M4-21)
CASES = {"sim50k": GOLDEN, "sim100k_real": GOLDEN_ROOT / "sim100k_real" / "m4"}
#: region types each case must test (the real SIM100K run has no AMR set)
REGION_TYPES = {"sim50k": {"prophage", "amr"}, "sim100k_real": {"prophage"}}
#: the 50-kb genome cannot reach the default 0.01 (see the module docstring)
GOLDEN_ALPHA = 0.05
#: regions_stats.tsv columns downstream of the LAPACK least squares of ``residualise``
#: (lead decision 9): compared numerically, not byte for byte
NUMERIC_COLUMNS = ("hl_shift_log2", "effect", "ci_low", "ci_high")
#: the SVD difference is 1-2 ulp (~1e-16 relative); this bound is still ~7 orders tighter than
#: any difference a changed statistic could hide
REL_TOL, ABS_TOL = 1e-9, 1e-12
NA = ("NA", "")


def _close(want: str, got: str) -> bool:
    if want in NA or got in NA:
        return want in NA and got in NA
    a, b = float(want), float(got)
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL)


def golden_config(seed: int) -> Any:
    """The shipped defaults (``stats.permutations`` and ``stats.bootstrap`` = 1,000) with the
    golden seed and ``alpha.identification`` = 0.05."""
    return mc.config(seed, stats={"alpha": {"identification": GOLDEN_ALPHA}})


@pytest.fixture(scope="module", params=sorted(CASES))
def case(request: pytest.FixtureRequest) -> str:
    return str(request.param)


@pytest.fixture(scope="module")
def golden(case: str) -> Path:
    return CASES[case]


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory, golden_seed: int, case: str) -> dict[str, Path]:
    root = tmp_path_factory.mktemp(f"m4_golden_{case}")
    db = mc.test_db(root / "db")
    if case == "sim100k_real":
        import m4_real as mr

        work = mr.giap_run(mr.sim100k(root / "sim"), root / "OUT", db, golden_seed)
    else:
        work = mc.build(root / "work", seed=golden_seed, size="sim50k")
    out = mc.run_m4(work, golden_config(golden_seed), db)
    if os.environ.get("GIAP_UPDATE_GOLDEN") == "1":
        dest = CASES[case]
        dest.mkdir(parents=True, exist_ok=True)
        for name, path in out.items():
            shutil.copy(path, dest / name)
        pytest.skip(f"golden files written to {dest}")
    return out


def test_m4_tables_equal_the_committed_golden_files(run: dict[str, Path], golden: Path) -> None:
    """Byte identity for every table but ``regions_stats.tsv`` (decision 9, next test)."""
    missing = [n for n in mc.OUTPUT_TABLES if not (golden / n).is_file()]
    assert not missing, f"golden files missing in {golden}: {missing}"
    for name, path in sorted(run.items()):
        if name == "regions_stats.tsv":
            continue
        want = golden / name
        if path.read_bytes() != want.read_bytes():
            pytest.fail(
                f"{name} differs from the golden file; first lines (n, want, got):"
                f" {mc.first_diff(want, path)}"
            )


def test_regions_stats_equals_the_golden_file_up_to_the_lstsq_floats(
    run: dict[str, Path], golden: Path
) -> None:
    """Lead decision 9: ``NUMERIC_COLUMNS`` are compared with ``REL_TOL`` / ``ABS_TOL``, every
    other cell -- counts, p-values, q-values, methods, labels, flags -- byte for byte."""
    want, got = mc.read(golden / "regions_stats.tsv"), mc.read(run["regions_stats.tsv"])
    assert list(got.columns) == list(want.columns), "regions_stats header changed"
    assert len(got) == len(want), f"{len(got)} rows, golden has {len(want)}"
    exact = [c for c in want.columns if c not in NUMERIC_COLUMNS]
    for col in exact:
        bad = [
            (i, a, b) for i, (a, b) in enumerate(zip(want[col], got[col], strict=True)) if a != b
        ][:5]
        assert not bad, f"regions_stats column {col} differs (row, want, got): {bad}"
    for col in NUMERIC_COLUMNS:
        assert col in want.columns, col
        bad = [
            (i, a, b)
            for i, (a, b) in enumerate(zip(want[col], got[col], strict=True))
            if not _close(a, b)
        ][:5]
        assert not bad, (
            f"regions_stats column {col} differs by more than rel {REL_TOL} / abs {ABS_TOL}"
            f" (row, want, got): {bad}; only 1-2 ulp of LAPACK noise is allowed here"
        )


def test_the_golden_regions_stats_floats_are_not_all_na(run: dict[str, Path], golden: Path) -> None:
    """The tolerance above must not be hiding an empty comparison."""
    want = mc.read(golden / "regions_stats.tsv")
    for col in NUMERIC_COLUMNS:
        assert [v for v in want[col] if v not in NA], f"{col} is NA in every golden row"


def test_the_float_tolerance_accepts_ulps_and_nothing_larger() -> None:
    """``REL_TOL`` must cover the LAPACK difference (1-2 ulp, the lead's two observed values) and
    still reject anything a changed statistic would produce."""
    iron, wheat = float.fromhex("0x1.85e19f7a40364p-6"), float.fromhex("0x1.85e19f7a40377p-6")
    assert iron != wheat and _close(repr(iron), repr(wheat))
    assert _close("-1.549019607843137", "-1.5490196078431372")
    assert not _close("0.5", "0.5000001"), "a 2e-7 relative change must fail"
    assert not _close("0.0", "1e-9"), "ABS_TOL must not swallow a 1e-9 absolute change"
    assert _close("NA", "NA") and not _close("NA", "0.5") and not _close("0.5", "")


def test_the_golden_tables_are_not_vacuous(run: dict[str, Path], golden: Path, case: str) -> None:
    """A golden of zeros would pass for ever without testing anything: the committed run must
    detect genes, test both planted regions and rerun both sensitivity variants."""
    genes = mc.read(golden / "genes.tsv") if (golden / "genes.tsv").is_file() else None
    assert genes is not None, f"{golden / 'genes.tsv'} is not committed"
    detected = [v for v in genes["detection_status"] if v == "detected"]
    assert 0 < len(detected) < len(genes), f"{len(detected)} detected of {len(genes)}"
    rs = mc.read(golden / "regions_stats.tsv")
    assert set(rs["region_type"]) >= REGION_TYPES[case]
    assert any(v not in ("NA", "") for v in rs["p_value"]), "no region test ran"
    ds = mc.read(golden / "detection_summary.tsv")
    assert set(ds["section"]) >= {"tier", "basic", "detectability", "qc"}
    sv = mc.read(golden / "sensitivity.tsv")
    assert sorted(set(sv["variant"])) == ["missing=exclude", "missing=min"]


def test_default_alpha_detects_nothing_on_50_kb(tmp_path: Path, golden_seed: int) -> None:
    """Known limitation, pinned so a bigger fixture or a changed floor shows up here: with the
    default ``alpha.identification`` = 0.01 the 51-gene fixture detects no gene at all, because
    the smallest attainable picked q is 1/29 = 0.0345."""
    work = mc.build(tmp_path / "work", seed=golden_seed, size="sim50k")
    mc.run_m4(work, mc.config(golden_seed), mc.test_db(tmp_path / "db"))
    g = mc.genes(work)
    assert len(g) < 100, "fixture grew: revisit the alpha override of the golden run"
    assert int(g["detected"].astype(bool).sum()) == 0
    assert float(g["gene_q_value"].min()) > 0.01
