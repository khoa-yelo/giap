"""Bakta 1.12.2 command line and database selection (GIAP_SPEC.md §7.1.1, §5.2; M2-03).

Builds ``bakta --db DB --output OUT --prefix P --threads N --keep-contig-headers
[--complete] --replicons replicons.tsv [--regions regions.gff3] [extra args] genome.norm.fa``.
Bakta itself is never imported (R5); it runs as a separate process.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Literal

DbType = Literal["full", "light"]


def bakta_db_path(db_setting: str, db_dir: Path) -> tuple[Path, DbType]:
    """Resolve ``annotation.bakta.db``: ``auto`` picks ``bakta/db`` if present, else ``db-light``.

    ``full``/``light`` select the matching directory under ``db_dir/bakta``; anything else is a
    path to a Bakta DB, typed ``light`` when its name contains ``light``.
    """
    base = db_dir / "bakta"
    if db_setting == "auto":
        return (base / "db", "full") if (base / "db").is_dir() else (base / "db-light", "light")
    if db_setting == "full":
        return base / "db", "full"
    if db_setting == "light":
        return base / "db-light", "light"
    path = Path(db_setting)
    return path, "light" if "light" in path.name else "full"


def bakta_command(
    *,
    genome: Path,
    outdir: Path,
    prefix: str,
    db: Path,
    threads: int,
    all_circular: bool,
    complete: bool | Literal["auto"],
    replicons_tsv: Path,
    regions_gff3: Path | None,
    extra_args: Sequence[str],
) -> list[str]:
    """Bakta argv; ``--replicons`` always, ``--complete`` iff forced or (auto) all circular."""
    cmd = ["bakta", "--db", str(db), "--output", str(outdir), "--prefix", prefix]
    cmd += ["--threads", str(threads), "--keep-contig-headers"]
    use_complete = all_circular if complete == "auto" else bool(complete)
    cmd += [*(["--complete"] if use_complete else []), "--replicons", str(replicons_tsv)]
    if regions_gff3 is not None:
        cmd += ["--regions", str(regions_gff3)]
    cmd += [str(a) for a in extra_args]
    cmd.append(str(genome))
    return cmd
