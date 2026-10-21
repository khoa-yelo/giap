"""The committed M4 tables of the 50-kb fixture (GIAP_SPEC.md §11.4; M4-20 example 2).

Fixture: ``build_m4_fixture(work, seed=GOLDEN_SEED, size="sim50k")`` (the M0 50-kb genome of
``simulate(SimConfig(), GOLDEN_SEED)`` as an M4-07 layout, plus the synthetic prophage block and
AMR set), then ``run_gene_table`` -> ``run_detectability`` -> ``run_context_stats`` ->
``run_sensitivity`` with ``GOLDEN_CONFIG``. The four ``tables/*.tsv`` must equal
``tests/data/golden/sim50k/m4/`` byte for byte.

``$GIAP_DB_DIR`` is a test database with the fixture ``cog-20.def.tab`` (tests/helpers/
gene_table_fixture.py), never the installed one: otherwise ``cog_category`` -- and so the golden
genes table -- would depend on whether ``giap setup --cog`` has run on the machine.

Detection on 50 kb: the genome has 47 genes and 29 of them carry peptide evidence, so the picked
gene q-value cannot fall below 1/29 = 0.0345. At the default ``stats.alpha.identification`` = 0.01
*nothing* is ever detected and every downstream number is a zero
(``test_default_alpha_detects_nothing_on_50_kb`` pins that known limitation, the one the MVP-01
run 1861972 hit on sim100k); the golden therefore runs at ``alpha.identification`` = 0.05, where
33 of the 47 genes are detected. Lead decision needed (M4-20 report, question 1).

Regenerating (on purpose, reviewed in the PR diff as data, with the statistician's note):
``GIAP_UPDATE_GOLDEN=1 pytest tests/golden/test_m4_golden.py`` writes the files and skips.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc

GOLDEN = Path(__file__).resolve().parents[1] / "data" / "golden" / "sim50k" / "m4"
#: the 50-kb genome cannot reach the default 0.01 (see the module docstring)
GOLDEN_ALPHA = 0.05


def golden_config(seed: int) -> Any:
    """The shipped defaults (``stats.permutations`` and ``stats.bootstrap`` = 1,000) with the
    golden seed and ``alpha.identification`` = 0.05."""
    return mc.config(seed, stats={"alpha": {"identification": GOLDEN_ALPHA}})


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory, golden_seed: int) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("m4_golden")
    work = mc.build(root / "work", seed=golden_seed, size="sim50k")
    return mc.run_m4(work, golden_config(golden_seed), mc.test_db(root / "db"))


def test_m4_tables_equal_the_committed_golden_files(run: dict[str, Path]) -> None:
    if os.environ.get("GIAP_UPDATE_GOLDEN") == "1":
        GOLDEN.mkdir(parents=True, exist_ok=True)
        for name, path in run.items():
            shutil.copy(path, GOLDEN / name)
        pytest.skip(f"golden files written to {GOLDEN}")
    missing = [n for n in mc.OUTPUT_TABLES if not (GOLDEN / n).is_file()]
    assert not missing, f"golden files missing in {GOLDEN}: {missing}"
    for name, path in sorted(run.items()):
        want = GOLDEN / name
        if path.read_bytes() != want.read_bytes():
            pytest.fail(
                f"{name} differs from the golden file; first lines (n, want, got):"
                f" {mc.first_diff(want, path)}"
            )


def test_the_golden_tables_are_not_vacuous(run: dict[str, Path]) -> None:
    """A golden of zeros would pass for ever without testing anything: the committed run must
    detect genes, test both planted regions and rerun both sensitivity variants."""
    genes = mc.read(GOLDEN / "genes.tsv") if (GOLDEN / "genes.tsv").is_file() else None
    assert genes is not None, f"{GOLDEN / 'genes.tsv'} is not committed"
    detected = [v for v in genes["detection_status"] if v == "detected"]
    assert 0 < len(detected) < len(genes), f"{len(detected)} detected of {len(genes)}"
    rs = mc.read(GOLDEN / "regions_stats.tsv")
    assert set(rs["region_type"]) >= {"prophage", "amr"}
    assert any(v not in ("NA", "") for v in rs["p_value"]), "no region test ran"
    ds = mc.read(GOLDEN / "detection_summary.tsv")
    assert set(ds["section"]) >= {"tier", "basic", "detectability", "qc"}
    sv = mc.read(GOLDEN / "sensitivity.tsv")
    assert sorted(set(sv["variant"])) == ["missing=exclude", "missing=min"]


def test_default_alpha_detects_nothing_on_50_kb(tmp_path: Path, golden_seed: int) -> None:
    """Known limitation, pinned so a bigger fixture or a changed floor shows up here: with the
    default ``alpha.identification`` = 0.01 the 47-gene fixture detects no gene at all, because
    the smallest attainable picked q is 1/29."""
    work = mc.build(tmp_path / "work", seed=golden_seed, size="sim50k")
    mc.run_m4(work, mc.config(golden_seed), mc.test_db(tmp_path / "db"))
    g = mc.genes(work)
    assert len(g) < 100, "fixture grew: revisit the alpha override of the golden run"
    assert int(g["detected"].astype(bool).sum()) == 0
    assert float(g["gene_q_value"].min()) > 0.01
