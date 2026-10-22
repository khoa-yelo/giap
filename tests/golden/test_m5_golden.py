"""The committed M5 tables of the two golden fixtures (GIAP_SPEC.md §11.4; M5-11 examples 1-2).

Cases (tests/helpers/m5_chain.py):

* ``sim50k`` -- ``build_m4_fixture(work, seed=GOLDEN_SEED, size="sim50k")`` and the M4 golden
  config (``alpha.identification`` = 0.05, M4-20 decision 1). ``islands.tsv`` is header-only:
  a 50-kb genome has about 25 TUs, too few for an island (lead decision 39 keeps it, it tests the
  empty path).
* the island case ``m5_chain.ISLAND`` -- an oriC at 0.8 and one planted +3 SD island of 20 TUs,
  the shipped defaults (B = 1,000, ``block_len: auto``). The spec names it ``sim1000_island``; the
  helper's docstring has the measurements behind ``sim4500_island`` (lead decision pending).

For each case ``tests/data/golden/<case>/m5/`` holds ``autocorrelation.tsv``, ``islands.tsv`` and
``genes_context.tsv`` (``locus_tag`` plus the three genes.tsv columns M5 fills: ``island_id``,
``replichore``, ``replichore_position``; gzipped for the island case, R11). The run must equal them
byte for byte (gzip compared decompressed). The M4 tables of ``sim50k`` stay compared by
tests/golden/test_m4_golden.py and, minus the three M5 columns, by
tests/integration/test_m5_pipeline.py.

**Decision 39: the island golden is checked against l' before it is frozen.** With
``GIAP_UPDATE_GOLDEN=1`` the files are written only when the run passes the l' check
(``m5_chain.l_prime_problems``: chr1 ``block_len: auto``, l' = min(20, 2 l-hat + 1) and shorter
than the planted island) and the example 2 check (one ``high`` island overlapping the planted TUs);
otherwise the regeneration fails and writes nothing. The same checks also run on every comparison.

Regenerating (on purpose, reviewed in the PR diff as data):
``GIAP_UPDATE_GOLDEN=1 pytest tests/golden/test_m5_golden.py`` writes the files and skips.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc
import m5_chain as m5

CASES = ("sim50k", m5.ISLAND.name)


@pytest.fixture(scope="module", params=CASES)
def case(request: pytest.FixtureRequest) -> str:
    return str(request.param)


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory, golden_seed: int, case: str) -> dict[str, Any]:
    root = tmp_path_factory.mktemp(f"m5_golden_{case}")
    db = mc.test_db(root / "db")
    if case == "sim50k":
        work = m5.build_sim50k(root / "work", golden_seed)
        tables = m5.run_m5(work, m5.sim50k_config(golden_seed), db)
    else:
        work = m5.build_island(root / "work", golden_seed)
        tables = m5.run_m5(work, m5.island_config(golden_seed), db)
    out = {"work": work, "tables": tables, "bytes": m5.produced(tables)}
    if m5.updating():
        if case != "sim50k":
            islands = mc.read(tables["islands.tsv"])
            bad = m5.l_prime_problems(work, islands) + m5.island_problems(work, islands)
            if bad:
                pytest.fail(f"refusing to freeze the {case} golden (decision 39): {bad}")
        dest = m5.golden_dir(case)
        for logical, fname in m5.golden_names(case).items():
            m5.write_golden(dest / fname, out["bytes"][logical])
        pytest.skip(f"golden files written to {dest}")
    return out


def test_m5_tables_equal_the_committed_golden_files(run: dict[str, Any], case: str) -> None:
    dest = m5.golden_dir(case)
    names = m5.golden_names(case)
    missing = [f for f in names.values() if not (dest / f).is_file()]
    assert not missing, f"golden files missing in {dest}: {missing}"
    for logical, fname in sorted(names.items()):
        want, got = m5.read_golden(dest / fname), run["bytes"][logical]
        if want != got:
            pytest.fail(
                f"{case}/{fname} differs from the golden file; first lines (n, want, got):"
                f" {m5.first_diff(want, got)}"
            )


def test_genes_context_golden_matches_the_m4_gene_list(run: dict[str, Any], case: str) -> None:
    """The golden genes_context rows are the run's genes.tsv rows in genes.tsv order."""
    dest = m5.golden_dir(case)
    path = dest / m5.golden_names(case)["genes_context.tsv"]
    if not path.is_file():
        pytest.fail(f"{path} is not committed")
    lines = m5.read_golden(path).decode().splitlines()
    assert lines[0].split("\t") == ["locus_tag", *m5.CONTEXT_COLUMNS]
    tags = [ln.split("\t")[0] for ln in lines[1:]]
    assert tags == list(mc.read(run["tables"]["genes.tsv"])["locus_tag"])


def test_sim50k_islands_golden_is_header_only() -> None:
    """Decision 39: the 50-kb islands golden tests the empty path (header, no rows)."""
    path = m5.golden_dir("sim50k") / "islands.tsv"
    assert path.is_file(), f"{path} is not committed"
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1, f"{len(lines) - 1} island rows on a 50-kb genome"
    stem = mc.descriptor("islands")
    assert lines[0].split("\t") == [f["name"] for f in stem["fields"]]


def test_island_golden_is_not_vacuous() -> None:
    """The committed island golden has exactly one ``high`` island, Moran's I and join-count rows
    under both nulls, the three ori-ter rows (``descriptive``), and filled gene columns."""
    dest, name = m5.golden_dir(m5.ISLAND.name), m5.golden_names(m5.ISLAND.name)
    for f in name.values():
        assert (dest / f).is_file(), f"{dest / f} is not committed"
    isl = mc.read(dest / "islands.tsv")
    assert len(isl) == 1 and list(isl["direction"]) == ["high"], isl.to_dict("records")
    ac = mc.read(dest / "autocorrelation.tsv")
    moran = ac[ac["statistic"] != "oriter_slope"]
    assert set(moran["null_model"]) >= {"gene_label", "tu_block"}, set(moran["null_model"])
    assert any(v not in ("NA", "") for v in moran["p_value"]), "no autocorrelation p-value"
    ori = ac[ac["statistic"] == "oriter_slope"]
    assert len(ori) == 3, ori.to_dict("records")
    assert all("descriptive" in f.split(";") for f in ori["flags"])
    assert not any("no_oric" in f.split(";") for f in ori["flags"])
    lines = m5.read_golden(dest / name["genes_context.tsv"]).decode().splitlines()
    rows = [ln.split("\t") for ln in lines[1:]]
    in_island = [r for r in rows if r[1] not in ("NA", "")]
    assert in_island and all(r[1] == isl["island_id"].iloc[0] for r in in_island)
    assert sorted(r[0] for r in in_island) == sorted(isl["genes"].iloc[0].split(";"))
    assert all(r[2] in ("right", "left") for r in rows), "a chr1 gene without a replichore"
    assert {r[2] for r in rows} == {"right", "left"}


def test_island_golden_passes_the_l_prime_check(run: dict[str, Any], case: str) -> None:
    """Decision 39 on every comparison, not only at regeneration."""
    if case == "sim50k":
        pytest.skip("no island on 50 kb")
    islands = mc.read(run["tables"]["islands.tsv"])
    bad = m5.l_prime_problems(run["work"], islands) + m5.island_problems(run["work"], islands)
    assert not bad, bad


def test_golden_data_stays_under_one_megabyte(case: str) -> None:
    """R11: the committed M5 golden data of a case is < 1 MB."""
    dest = m5.golden_dir(case)
    files = [dest / f for f in m5.golden_names(case).values()]
    assert all(p.is_file() for p in files), [str(p) for p in files if not p.is_file()]
    total = sum(p.stat().st_size for p in files)
    assert total < m5.MAX_GOLDEN_BYTES, f"{case}: {total} bytes"


def test_gzip_golden_round_trip_is_reproducible(tmp_path: Path) -> None:
    """``write_golden`` gives the same .gz bytes twice (mtime 0) and ``read_golden`` inverts it."""
    data = b"locus_tag\tisland_id\nG1\tISL_1\n"
    a, b = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    m5.write_golden(a, data)
    m5.write_golden(b, data)
    assert a.read_bytes() == b.read_bytes()
    assert m5.read_golden(a) == data
