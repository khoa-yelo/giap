"""Bakta command builder and DB selection (issue M2-03 example 10; GIPPA_SPEC.md §7.1.1, §5.2)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

GENOME = Path("/w/genome.norm.fa")
OUT = Path("/w/bakta")
DB = Path("/db/bakta/db")
REPLICONS = Path("/w/replicons.tsv")


def _cmd(**kw: Any) -> list[str]:
    from gippa.annotate.bakta_cmd import bakta_command

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


# M2-07 lead decision 11: --replicons <table> is always passed (the table carries the c/l topology
# and the replicon types); --complete is added only when every replicon is circular (or forced).
# The order of the mode flags among themselves is free; they sit after _HEAD and before
# --regions / extra args / the genome.


def _mode(cmd: list[str]) -> list[str]:
    """The flags between _HEAD and --regions/extra args/genome, as a sorted list of tokens."""
    assert cmd[: len(_HEAD)] == _HEAD, cmd
    rest = cmd[len(_HEAD) : -1]
    out: list[str] = []
    i = 0
    while i < len(rest) and rest[i] in {"--complete", "--replicons"}:
        if rest[i] == "--replicons":
            out += [rest[i], rest[i + 1]]
            i += 2
        else:
            out.append(rest[i])
            i += 1
    return out


def _pairs(mode: list[str]) -> set[tuple[str, ...]]:
    got: set[tuple[str, ...]] = set()
    i = 0
    while i < len(mode):
        if mode[i] == "--replicons":
            got.add((mode[i], mode[i + 1]))
            i += 2
        else:
            got.add((mode[i],))
            i += 1
    return got


def test_example10_all_circular_auto_complete_plus_replicons() -> None:
    cmd = _cmd()
    assert cmd[-1] == str(GENOME)
    assert _pairs(_mode(cmd)) == {("--complete",), ("--replicons", str(REPLICONS))}
    assert len(cmd) == len(_HEAD) + 4


def test_example10_not_all_circular_uses_replicons_table_only() -> None:
    cmd = _cmd(all_circular=False)
    assert cmd == [*_HEAD, "--replicons", str(REPLICONS), str(GENOME)]


def test_complete_true_forces_complete_and_keeps_replicons() -> None:
    cmd = _cmd(all_circular=False, complete=True)
    assert _pairs(_mode(cmd)) == {("--complete",), ("--replicons", str(REPLICONS))}
    assert len(cmd) == len(_HEAD) + 4


def test_complete_false_uses_replicons_table_only() -> None:
    assert _cmd(complete=False) == [*_HEAD, "--replicons", str(REPLICONS), str(GENOME)]


@pytest.mark.parametrize("all_circular", [True, False])
@pytest.mark.parametrize("complete", [True, False, "auto"])
def test_replicons_always_passed_decision_11(all_circular: bool, complete: Any) -> None:
    cmd = _cmd(all_circular=all_circular, complete=complete)
    assert cmd.count("--replicons") == 1, cmd
    assert cmd[cmd.index("--replicons") + 1] == str(REPLICONS)
    want_complete = all_circular if complete == "auto" else bool(complete)
    assert ("--complete" in cmd) is want_complete, cmd


def test_regions_appended() -> None:
    regions = Path("/w/regions.gff3")
    cmd = _cmd(regions_gff3=regions)
    mode = _mode(cmd)
    assert _pairs(mode) == {("--complete",), ("--replicons", str(REPLICONS))}
    assert cmd[len(_HEAD) + len(mode) :] == ["--regions", str(regions), str(GENOME)]


def test_extra_args_last_before_genome() -> None:
    regions = Path("/w/regions.gff3")
    cmd = _cmd(regions_gff3=regions, extra_args=["--skip-plot", "--min-contig-length", "200"])
    mode = _mode(cmd)
    assert _pairs(mode) == {("--complete",), ("--replicons", str(REPLICONS))}
    assert cmd[len(_HEAD) + len(mode) :] == [
        "--regions",
        str(regions),
        "--skip-plot",
        "--min-contig-length",
        "200",
        str(GENOME),
    ]


def test_returns_list_of_str() -> None:
    from gippa.annotate.bakta_cmd import bakta_command

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
    from gippa.annotate.bakta_cmd import bakta_db_path

    for name in present:
        (tmp_path / "bakta" / name).mkdir(parents=True)
    path, db_type = bakta_db_path("auto", tmp_path)
    assert path == tmp_path / "bakta" / expected_dir
    assert db_type == expected_type
