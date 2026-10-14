"""Bakta command builder and DB selection (issue M2-03 example 10; GIAP_SPEC.md §7.1.1, §5.2)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

GENOME = Path("/w/genome.norm.fa")
OUT = Path("/w/bakta")
DB = Path("/db/bakta/db")
REPLICONS = Path("/w/replicons.tsv")


def _cmd(**kw: Any) -> list[str]:
    from giap.annotate.bakta_cmd import bakta_command

    args: dict[str, Any] = {
        "genome": GENOME,
        "outdir": OUT,
        "prefix": "P",
        "db": DB,
        "threads": 16,
        "all_circular": True,
        "complete": "auto",
        "replicons_tsv": REPLICONS,
        "regions_gff3": None,
        "extra_args": [],
    }
    args.update(kw)
    return [str(x) for x in bakta_command(**args)]


_HEAD = ["bakta", "--db", str(DB), "--output", str(OUT), "--prefix", "P", "--threads", "16",
         "--keep-contig-headers"]  # fmt: skip


def test_example10_all_circular_auto_complete() -> None:
    assert _cmd() == [*_HEAD, "--complete", str(GENOME)]


def test_example10_not_all_circular_uses_replicons_table() -> None:
    assert _cmd(all_circular=False) == [*_HEAD, "--replicons", str(REPLICONS), str(GENOME)]


def test_complete_true_forces_complete() -> None:
    assert _cmd(all_circular=False, complete=True) == [*_HEAD, "--complete", str(GENOME)]


def test_complete_false_uses_replicons_table() -> None:
    assert _cmd(complete=False) == [*_HEAD, "--replicons", str(REPLICONS), str(GENOME)]


def test_regions_appended() -> None:
    regions = Path("/w/regions.gff3")
    assert _cmd(regions_gff3=regions) == [
        *_HEAD, "--complete", "--regions", str(regions), str(GENOME)
    ]  # fmt: skip


def test_extra_args_last_before_genome() -> None:
    regions = Path("/w/regions.gff3")
    cmd = _cmd(regions_gff3=regions, extra_args=["--skip-plot", "--min-contig-length", "200"])
    assert cmd == [
        *_HEAD,
        "--complete",
        "--regions",
        str(regions),
        "--skip-plot",
        "--min-contig-length",
        "200",
        str(GENOME),
    ]


def test_returns_list_of_str() -> None:
    from giap.annotate.bakta_cmd import bakta_command

    cmd = bakta_command(
        genome=GENOME,
        outdir=OUT,
        prefix="P",
        db=DB,
        threads=4,
        all_circular=True,
        complete="auto",
        replicons_tsv=REPLICONS,
        regions_gff3=None,
        extra_args=(),
    )
    assert isinstance(cmd, list)
    assert all(isinstance(x, str) for x in cmd)


@pytest.mark.parametrize(
    ("present", "expected_dir", "expected_type"),
    [
        (["db", "db-light"], "db", "full"),
        (["db"], "db", "full"),
        (["db-light"], "db-light", "light"),
        ([], "db-light", "light"),
    ],
)
def test_db_auto(tmp_path: Path, present: list[str], expected_dir: str, expected_type: str) -> None:
    from giap.annotate.bakta_cmd import bakta_db_path

    for name in present:
        (tmp_path / "bakta" / name).mkdir(parents=True)
    path, db_type = bakta_db_path("auto", tmp_path)
    assert path == tmp_path / "bakta" / expected_dir
    assert db_type == expected_type
