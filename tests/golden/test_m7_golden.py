"""The committed Module 5 and 6 tables of the 50-kb M7 fixture (docs/specs/M7-15-integration-
hardening.md item 1, example 1; GIAP_SPEC.md §11.2 golden layer, §11.4).

Lead decision 1 (2026-10-22): a second golden directory ``tests/data/golden/sim50k_m7/`` with the
M7 settings (tests/helpers/m7_chain.py: a 50-kb giapsim genome, SNV strain difference, de novo
peptides and their reference identifications, planted novel ORFs with the M7-12 knobs, every
layer precomputed); the sim50k goldens are not regenerated.

``giap run --executor local --modules core --seed GOLDEN_SEED`` on the reference genome must give
``pileup.tsv.gz``, ``variants.tsv``, ``variants.vcf``, ``novel_orfs.tsv`` and
``annotation_patch.gff3`` equal to the committed files byte for byte (gzip compared decompressed;
the VCF ``##source=giap-<version>`` line with the version masked; BLAS-dependent numeric columns,
none today, per the M5-11 decision 44 pattern: ``m7_chain.BLAS_NUMERIC``).

**Freezing checks.** With ``GIAP_UPDATE_GOLDEN=1`` the files are written only when giapsim has the
M7-12 knobs and the run is not vacuous (``m7_chain.golden_problems``: a ``variants.tsv`` row at a
planted SNV and one flagged ``strain_difference``, a called novel ORF on a planted one, a CDS in the
patch, a novel-ORF class whose M7-11b class q floor is <= alpha, pileup rows); otherwise the
regeneration fails and writes nothing. The same checks run on every comparison.

Regenerating (on purpose, reviewed in the PR diff as data):
``GIAP_UPDATE_GOLDEN=1 pytest tests/golden/test_m7_golden.py`` writes the files and skips.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m7_chain as m7


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory, golden_seed: int) -> dict[str, Any]:
    root = tmp_path_factory.mktemp("m7_golden")
    sim = m7.simulate(root / "sim", golden_seed)
    with m7.db_env(root):
        out = m7.run_local(sim, root / "RUN" / "OUT", golden_seed)
    got = {"sim": sim, "out": out, "bytes": m7.produced(out)}
    if m7.updating():
        if not m7.has_m7_12():
            pytest.fail("refusing to freeze sim50k_m7: giapsim lacks the M7-12 knobs")
        bad = m7.golden_problems(sim, out)
        if bad:
            pytest.fail(f"refusing to freeze sim50k_m7: {bad}")
        for rel, fname in m7.GOLDEN_NAMES.items():
            m7.write_golden(m7.GOLDEN_DIR / fname, got["bytes"][rel])
        pytest.skip(f"golden files written to {m7.GOLDEN_DIR}")
    return got


def test_example_1_tables_equal_the_committed_golden_files(run: dict[str, Any]) -> None:
    missing = [f for f in m7.GOLDEN_NAMES.values() if not (m7.GOLDEN_DIR / f).is_file()]
    assert not missing, f"golden files missing in {m7.GOLDEN_DIR}: {missing}"
    differ = []
    for rel, fname in sorted(m7.GOLDEN_NAMES.items()):
        want, got = m7.read_golden(m7.GOLDEN_DIR / fname), run["bytes"][rel]
        if not m7.golden_equal(rel, want, got):
            differ.append((fname, m7.first_diff(want, got)))
    assert not differ, f"differs from the golden (file, first lines n/want/got): {differ}"


def test_golden_fixture_uses_the_m7_12_settings() -> None:
    """Lead decision 1: "the M7 settings" -- the M7-12 novel-ORF length range and planted RBSs."""
    assert m7.has_m7_12(), "giapsim.SimConfig lacks the M7-12 knobs " + str(sorted(m7.M7_12_KNOBS))
    cfg = m7.sim_config()
    for key, value in m7.M7_12_KNOBS.items():
        assert getattr(cfg, key) == value, key


def test_run_is_not_a_vacuous_golden(run: dict[str, Any]) -> None:
    assert not m7.golden_problems(run["sim"], run["out"])


def test_golden_variants_and_novel_tables_are_not_vacuous() -> None:
    """The committed files themselves: rows in each, PASS-only VCF records."""
    for fname in m7.GOLDEN_NAMES.values():
        assert (m7.GOLDEN_DIR / fname).is_file(), f"{m7.GOLDEN_DIR / fname} is not committed"
    var = m7.read_golden(m7.GOLDEN_DIR / "variants.tsv").decode().splitlines()
    assert len(var) > 1
    nov = m7.read_golden(m7.GOLDEN_DIR / "novel_orfs.tsv").decode().splitlines()
    assert len(nov) > 1
    recs = [ln for ln in m7.read_golden(m7.GOLDEN_DIR / "variants.vcf").decode().splitlines()
            if ln and not ln.startswith("#")]  # fmt: skip
    assert recs and all(ln.split("\t")[6] == "PASS" for ln in recs)


def test_golden_data_stays_under_one_megabyte() -> None:
    """R11: the committed sim50k_m7 data is < 1 MB."""
    paths = [m7.GOLDEN_DIR / f for f in m7.GOLDEN_NAMES.values()]
    assert all(p.is_file() for p in paths), [str(p) for p in paths if not p.is_file()]
    total = sum(p.stat().st_size for p in paths)
    assert total < m7.MAX_GOLDEN_BYTES, total


def test_sim50k_goldens_are_not_touched() -> None:
    """Decision 1: the sim50k directory keeps its m4 / m5 / mapping goldens only."""
    sim50k = m7.GOLDEN_DIR.parent / "sim50k"
    assert sorted(p.name for p in sim50k.iterdir()) == ["m4", "m5", "mapping"]


# --------------------------------------------------------------------------- comparator


def test_mask_version_masks_only_the_source_line() -> None:
    data = b"##fileformat=VCFv4.3\n##source=giap-0.0.1\n#CHROM\tPOS\nchr1\t5\n"
    assert m7.mask_version(data) == (
        b"##fileformat=VCFv4.3\n##source=giap-VERSION\n#CHROM\tPOS\nchr1\t5\n"
    )
    assert m7.mask_version(data.replace(b"0.0.1", b"1.2.0")) == m7.mask_version(data)


def test_tsv_equal_loosens_only_the_listed_columns() -> None:
    want = b"a\tb\n0.12345\t0.5\n"
    near = b"a\tb\n0.12346\t0.5\n"
    assert not m7.tsv_equal(want, near)
    assert m7.tsv_equal(want, near, ("a",))
    assert not m7.tsv_equal(want, b"a\tb\n0.12345\t0.50001\n", ("a",))
    assert not m7.tsv_equal(want, b"a\tb\nNA\t0.5\n", ("a",))
    assert not m7.tsv_equal(want, b"a\tc\n0.12345\t0.5\n", ("a",))


def test_blas_numeric_names_only_golden_tsv_columns() -> None:
    for fname, cols in m7.BLAS_NUMERIC.items():
        assert fname in m7.GOLDEN_NAMES.values() and fname.endswith(".tsv"), fname
        path = m7.GOLDEN_DIR / fname
        if path.is_file():
            header = m7.read_golden(path).decode().splitlines()[0].split("\t")
            assert set(cols) <= set(header), (fname, cols)


def test_gzip_golden_is_compared_decompressed(tmp_path: Path) -> None:
    data = b"contig\tcodon_start\nchr1\t1\n"
    m7.write_golden(tmp_path / "p.tsv.gz", data)
    assert m7.read_golden(tmp_path / "p.tsv.gz") == data
    assert gzip.decompress((tmp_path / "p.tsv.gz").read_bytes()) == data
