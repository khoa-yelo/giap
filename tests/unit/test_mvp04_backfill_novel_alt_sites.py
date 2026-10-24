"""``giapsim.real_genome.backfill_novel_alt_sites`` (PR #324 review item 5).

Contract (module docstring and the function's docstring): on an input directory written before
``truth/novel_alt_sites.tsv`` existed, rerun ``simulate_on_genome`` with the same config, seed
and inputs in a scratch directory and, only when every other file (``provenance.json`` aside) is
byte-identical and neither side has a file the other lacks, add exactly the new table. A second
call is a no-op. Any mismatch (other seed, a changed byte, an extra or a missing file, an
existing table that differs) is a ``ValueError`` and leaves the directory byte-unchanged; no
scratch directory is left behind.
"""

from __future__ import annotations

import hashlib
import importlib
import shutil
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from helpers import mvp04_data as M

from giapsim.expression import IslandSpec

SEED = 20261023
TABLE = "truth/novel_alt_sites.tsv"


def _rg() -> Any:
    return importlib.import_module("giapsim.real_genome")


def _cfg() -> Any:
    base = importlib.import_module("giapsim.scenarios").k12_evidence()
    strain = replace(base.strain, snv_rate_per_kb=5.0)  # as test_mvp04_novel_alt_sites
    expr = replace(base.expression, islands=(IslandSpec(M.CONTIG, 2, 8, 3.0),))
    return replace(base, strain=strain, expression=expr, n_novel_per_class=20)


@pytest.fixture(scope="module")
def ex(tmp_path_factory: pytest.TempPathFactory) -> M.Excerpt:
    return M.make_excerpt(tmp_path_factory.mktemp("k12_excerpt_backfill"))


@pytest.fixture(scope="module")
def direct(ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A direct run (it writes the table): the reference bytes."""
    out = tmp_path_factory.mktemp("backfill_direct") / "inputs"
    _simulate(ex, out, SEED)
    assert (out / TABLE).stat().st_size > 0
    assert len(M.read_tsv(out / TABLE)) >= 1, "need a non-empty table to test the backfill"
    return out


def _simulate(ex: M.Excerpt, out: Path, seed: int) -> None:
    with M.capture_warnings():
        _rg().simulate_on_genome(_cfg(), seed, out, genome_fasta=ex.fasta,
                                 bakta_json=ex.bakta_json, tu_set=None,
                                 check_hashes=False)  # fmt: skip


def _old(direct: Path, tmp: Path) -> Path:
    """The direct run as written before the table existed."""
    d = tmp / "old"
    shutil.copytree(direct, d)
    (d / TABLE).unlink()
    return d


def _hashes(d: Path) -> dict[str, str]:
    return {p.relative_to(d).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(d.rglob("*")) if p.is_file()}  # fmt: skip


def _backfill(ex: M.Excerpt, d: Path, scratch: Path, seed: int = SEED) -> Path:
    scratch.mkdir(exist_ok=True)
    with M.capture_warnings():
        path: Path = _rg().backfill_novel_alt_sites(
            _cfg(), seed, d, genome_fasta=ex.fasta, bakta_json=ex.bakta_json, tu_set=None,
            check_hashes=False, scratch=scratch,
        )  # fmt: skip
    return path


def test_backfill_writes_exactly_the_new_table(ex: M.Excerpt, direct: Path, tmp_path: Path) -> None:
    d = _old(direct, tmp_path)
    before = _hashes(d)
    path = _backfill(ex, d, tmp_path / "scratch")
    after = _hashes(d)
    assert path == d / TABLE
    assert set(after) == set(before) | {TABLE}
    assert {k: v for k, v in after.items() if k != TABLE} == before
    assert (d / TABLE).read_bytes() == (direct / TABLE).read_bytes()
    assert not (d / (TABLE + ".part")).exists()
    assert list((tmp_path / "scratch").iterdir()) == []


def test_a_second_backfill_is_a_no_op(ex: M.Excerpt, direct: Path, tmp_path: Path) -> None:
    d = _old(direct, tmp_path)
    _backfill(ex, d, tmp_path / "scratch")
    before, mtime = _hashes(d), (d / TABLE).stat().st_mtime_ns
    assert _backfill(ex, d, tmp_path / "scratch") == d / TABLE
    assert _hashes(d) == before
    assert (d / TABLE).stat().st_mtime_ns == mtime, "an identical table is left as is"
    assert list((tmp_path / "scratch").iterdir()) == []


def test_backfill_of_a_direct_run_is_a_no_op(ex: M.Excerpt, direct: Path, tmp_path: Path) -> None:
    d = tmp_path / "full"
    shutil.copytree(direct, d)
    before = _hashes(d)
    _backfill(ex, d, tmp_path / "scratch")
    assert _hashes(d) == before


def _changed_peptides(d: Path) -> None:
    p = d / "peptides.tsv"
    b = bytearray(p.read_bytes())
    i = b.index(b"\n") + 1  # first byte of the first data row
    b[i] = ord("A") if b[i] != ord("A") else ord("C")
    p.write_bytes(bytes(b))


def _extra(d: Path) -> None:
    (d / "truth" / "extra.tsv").write_text("x\n")


def _missing(d: Path) -> None:
    (d / "truth" / "genes.tsv").unlink()


def _differing_table(d: Path) -> None:
    (d / TABLE).write_text("\t".join(_rg().NOVEL_ALT_COLUMNS) + "\n")  # header only


@pytest.mark.parametrize(
    "tamper", [_changed_peptides, _extra, _missing, _differing_table],
    ids=["changed_byte", "extra_file", "missing_file", "differing_table"],
)  # fmt: skip
def test_backfill_refuses_a_mismatch_and_changes_nothing(
    ex: M.Excerpt, direct: Path, tmp_path: Path, tamper: Any
) -> None:
    d = _old(direct, tmp_path)
    tamper(d)
    before = _hashes(d)
    with pytest.raises(ValueError):
        _backfill(ex, d, tmp_path / "scratch")
    assert _hashes(d) == before
    assert list((tmp_path / "scratch").iterdir()) == []


def test_backfill_refuses_another_seed_and_changes_nothing(
    ex: M.Excerpt, direct: Path, tmp_path: Path
) -> None:
    d = _old(direct, tmp_path)
    before = _hashes(d)
    with pytest.raises(ValueError):
        _backfill(ex, d, tmp_path / "scratch", seed=SEED + 1)
    assert _hashes(d) == before
    assert list((tmp_path / "scratch").iterdir()) == []


def test_backfill_refuses_a_non_input_directory(ex: M.Excerpt, tmp_path: Path) -> None:
    d = tmp_path / "empty"
    d.mkdir()
    with pytest.raises(ValueError):
        _backfill(ex, d, tmp_path / "scratch")
    assert _hashes(d) == {}
