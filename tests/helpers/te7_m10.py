"""Builders for the M10-06b tests (TE7: ``validation/pilot/te7.py``, the ``te_real.py te7``
subcommand, the TE7 half of the ``metrics.py`` hook; docs/specs/M10-06-te-real-and-te7.md).

Hand-made run trees, ``logs/trace.tsv`` files with real Nextflow process names,
``outputs.sha256`` lists, converted peptide tables and a temporary git repository for the
code-equivalence check. No network, nothing from ``$GIPPA_ROOT``, nothing derived from a DIA-NN
report.
"""

from __future__ import annotations

import csv
import gzip
import importlib
import json
import sys
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from helpers.te_real_m10 import REPO, edit_manifest, git, sha256_bytes

TE7_PY = REPO / "validation" / "pilot" / "te7.py"
KNOWN_COLUMNS = REPO / "validation" / "pilot" / "te7_known_columns.tsv"

__all__ = ["REPO", "edit_manifest", "git", "sha256_bytes"]


def load_te7() -> Any:
    """``validation.pilot.te7`` (M10-06b); fails the calling test (not an error) when absent."""
    import pytest

    if not TE7_PY.is_file():
        pytest.fail(
            f"{TE7_PY.relative_to(REPO)} does not exist (M10-06b: TE7 moves out of te_real.py)"
        )
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return importlib.import_module("validation.pilot.te7")


# ------------------------------------------------------------------ TE7 trees


def reader(b: bytes) -> Callable[[], bytes]:
    return lambda: b


def readers(files: Mapping[str, bytes]) -> dict[str, Callable[[], bytes]]:
    """A ``classify_tree`` file map: path -> zero-argument bytes reader."""
    return {p: reader(b) for p, b in files.items()}


def gz(data: bytes, mtime: int = 0) -> bytes:
    import io

    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=mtime, filename="") as fh:
        fh.write(data)
    return buf.getvalue()


GENES_TSV = (
    b"locus_tag\tgene\tabundance_log2\tn_peptides_unique\n"
    b"T_00010\tthrA\t10.0\t3\n"
    b"T_00020\tthrB\t11.5\t2\n"
    b"T_00030\tthrC\t12.25\t4\n"
    b"T_00040\tyaaA\tNA\t0\n"
)
#: the same table, three abundance_log2 cells moved by one ulp-scale step (<= 2e-15)
GENES_TSV_NUMERIC = (
    b"locus_tag\tgene\tabundance_log2\tn_peptides_unique\n"
    b"T_00010\tthrA\t10.000000000000002\t3\n"
    b"T_00020\tthrB\t11.500000000000002\t2\n"
    b"T_00030\tthrC\t12.250000000000002\t4\n"
    b"T_00040\tyaaA\tNA\t0\n"
)
ISLANDS_TSV = b"island_id\tcontig\tstart\tend\nI1\tchr1\t10\t500\nI2\tchr1\t900\t1500\n"
ISLANDS_TSV_REORDERED = b"island_id\tcontig\tstart\tend\nI2\tchr1\t900\t1500\nI1\tchr1\t10\t500\n"
REGIONS_TSV = (
    b"region_id\tmetric\tn_genes\teffect\tp_value\n"
    b"prophage_1\tdetected\t12\t0.25\t0.031\n"
    b"bgc_1\tdetected\t8\t-0.5\t0.2\n"
)
REGIONS_TSV_LAPACK = (
    b"region_id\tmetric\tn_genes\teffect\tp_value\n"
    b"prophage_1\tdetected\t12\t0.25000000000000006\t0.031000000000000003\n"
    b"bgc_1\tdetected\t8\t-0.5\t0.2\n"
)
#: known LAPACK-dependent columns for the pure tests (the committed file is checked separately)
KNOWN = {"regions_stats.tsv": {"effect", "p_value"}}
PILEUP = b"contig\tpos\tdepth\nchr1\t1\t3\nchr1\t2\t4\n"


def base_tree() -> dict[str, bytes]:
    """Example 7's reference tree (every class present)."""
    return {
        "tables/genes.tsv": GENES_TSV,
        "tables/islands.tsv": ISLANDS_TSV,
        "tables/regions_stats.tsv": REGIONS_TSV,
        "tables/pileup.tsv.gz": gz(PILEUP, 0),
        "mapping/peptide_hits.tsv": b"peptide\troute\nAAAK\texact\n",
        "annotation/genes.gff3": b"##gff-version 3\n",
        "annotation/tools/promotech/promotech.csv": b"chrom,start,score\nchr1,10,0.91\n",
        "plots/data/landscape.tsv": b"x\ty\n1\t2\n",
        "plots/landscape.svg": b"<svg></svg>\n",
        "schemas/datapackage.json": b'{"resources": []}\n',
        "work/ab/cdef12/.command.log": b"ref work\n",
        "logs/trace.tsv": b"task_id\thostname\n1\trye-02\n",
        ".nextflow.log": b"started ref\n",
        "run_manifest.json": b'{"outdir": "/ref"}\n',
        "report/report.html": b"<html>ref</html>\n",
        "preflight/preflight.json": b'{"x": 1}\n',
    }


def ex7_replicate() -> dict[str, bytes]:
    rep = dict(base_tree())
    rep["work/ab/cdef12/.command.log"] = b"rep work\n"
    rep["logs/trace.tsv"] = b"task_id\thostname\n1\twheat-03\n"
    rep[".nextflow.log"] = b"started rep\n"
    rep["run_manifest.json"] = b'{"outdir": "/rep1"}\n'
    rep["report/report.html"] = b"<html>rep</html>\n"
    rep["preflight/preflight.json"] = b'{"x": 2}\n'
    rep["annotation/tools/promotech/promotech.csv"] = b"chrom,start,score\nchr1,10,0.92\n"
    return rep


def ex8_replicate() -> dict[str, bytes]:
    rep = dict(base_tree())
    rep["tables/genes.tsv"] = GENES_TSV_NUMERIC
    rep["tables/pileup.tsv.gz"] = gz(PILEUP, 1_700_000_000)
    rep["tables/islands.tsv"] = ISLANDS_TSV_REORDERED
    rep["tables/regions_stats.tsv"] = REGIONS_TSV_LAPACK
    rep["mapping/extra.tsv"] = b"a\tb\n1\t2\n"
    return rep


def producers_for(path: str, step: str, ref: Sequence[str], rep: Sequence[str]) -> dict[str, Any]:
    return {path: {"producer": step, "cpu_classes": {"ref": list(ref), "rep1": list(rep)}}}


def write_tree(root: Path, files: Mapping[str, bytes]) -> Path:
    for rel, data in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
    return root


# ------------------------------------------------------------------ git repository (TE7 pair)


def make_repo(root: Path) -> dict[str, str]:
    """A temporary repository with three commits: ``base``, ``docs`` (a change under ``docs/``
    only) and ``src`` (a change under ``src/``). Returns name -> commit SHA."""
    root.mkdir(parents=True, exist_ok=True)
    git(root, "init", "-q")
    for rel, text in (
        ("src/gippa/a.py", "x = 1\n"),
        ("docs/notes.md", "a\n"),
        ("workflow/main.nf", "// nf\n"),
        ("pyproject.toml", "[project]\n"),
        ("uv.lock", "lock\n"),
        ("scripts/slurm/gippa_run.sbatch", "#!/bin/bash\n"),
    ):
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "base")
    shas = {"base": git(root, "rev-parse", "HEAD")}
    (root / "docs" / "notes.md").write_text("b\n", encoding="utf-8")
    git(root, "commit", "-q", "-am", "docs only")
    shas["docs"] = git(root, "rev-parse", "HEAD")
    (root / "src" / "gippa" / "a.py").write_text("x = 2\n", encoding="utf-8")
    git(root, "commit", "-q", "-am", "code")
    shas["src"] = git(root, "rev-parse", "HEAD")
    return shas


def write_run(
    out: Path,
    files: Mapping[str, bytes],
    *,
    commit: str,
    outdir_label: str,
    hosts: Sequence[str] = ("rye-02",),
    **manifest_over: Any,
) -> Path:
    """A run tree for ``te_real.py te7``: ``files`` plus a schema-valid ``run_manifest.json``
    (``git_commit`` = ``commit``, ``git_dirty`` false, ``resume`` false, no cached step) and a
    ``logs/trace.tsv`` with ``hosts``. ``manifest_over`` replaces top-level manifest keys."""
    from helpers.pilot_metrics_m10 import write_manifest, write_trace

    write_tree(
        out, {k: v for k, v in files.items() if k not in ("run_manifest.json", "logs/trace.tsv")}
    )
    path = write_manifest(
        out,
        genome_sha="a" * 64,
        peptides_sha="b" * 64,
        steps=[
            {
                "name": "gene_table",
                "status": "completed",
                "n_tasks": 1,
                "n_retries": 0,
                "wall_time_s": 30.0,
                "cpu_time_s": 25.0,
            },
            {
                "name": "context_stats",
                "status": "completed",
                "n_tasks": 1,
                "n_retries": 0,
                "wall_time_s": 20.0,
                "cpu_time_s": 18.0,
            },
        ],
    )
    man = json.loads(path.read_text(encoding="utf-8"))
    man.update(
        {
            "git_commit": commit,
            "git_dirty": False,
            "resume": False,
            "outdir": f"/runs/{outdir_label}",
            "command_line": ["gippa", "run", "--outdir", f"/runs/{outdir_label}"],
        }
    )
    man.update(manifest_over)
    path.write_text(json.dumps(man, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    write_trace(out, hosts)
    return out


CONV_COLS = (
    "run",
    "modified_peptide",
    "peptide",
    "charge",
    "is_decoy",
    "q_value",
    "intensity",
    "protein_ids",
)


def conv_frame(rows: Sequence[tuple[str, str, int, bool, float, float]]) -> pd.DataFrame:
    """A converted table: rows ``(run, modified_peptide, charge, is_decoy, q_value, intensity)``;
    the stripped ``peptide`` is ``modified_peptide`` without ``(…)`` modifications."""
    import re

    return pd.DataFrame(
        [
            {
                "run": r,
                "modified_peptide": m,
                "peptide": re.sub(r"\(.*?\)", "", m),
                "charge": c,
                "is_decoy": d,
                "q_value": q,
                "intensity": i,
                "protein_ids": "JMJJOG_00010",
            }
            for r, m, c, d, q, i in rows
        ],
        columns=list(CONV_COLS),
    )


def sums_text(sums: Mapping[str, str]) -> str:
    """An ``outputs.sha256`` list (``sha256sum`` format)."""
    return "".join(f"{h}  {name}\n" for name, h in sorted(sums.items()))


def write_conv(path: Path, frame: pd.DataFrame) -> Path:
    """The converted table as gzip TSV with a fixed mtime (M10-03 layout)."""
    from gippa.io.tsv import format_value

    lines = ["\t".join(frame.columns)]
    lines += ["\t".join(format_value(v) for v in row) for row in frame.itertuples(index=False)]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gz(("\n".join(lines) + "\n").encode("utf-8")))
    return path


# ------------------------------------------------------------------ traces and manifest inputs


def write_named_trace(out: Path, rows: Sequence[tuple[str, str]]) -> Path:
    """``OUT/logs/trace.tsv`` with one task per ``(Nextflow name, hostname)``, for example
    ``("GIPPA:VARIANTS:VARIANTS_STEP (chr1)", "rye-02")`` (the real trace's ``name`` column)."""
    from helpers.manifest_data import trace_row

    recs = [trace_row(n, hostname=h, task_id=str(i + 1)) for i, (n, h) in enumerate(rows)]
    path = out / "logs" / "trace.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]), delimiter="	")
        w.writeheader()
        w.writerows(recs)
    return path


def set_peptides_sha(out: Path, sha: str) -> None:
    """Set the ``peptides`` input sha256 of ``OUT/run_manifest.json``."""

    def fn(m: dict[str, Any]) -> None:
        for i in m["inputs"]:
            if i.get("role") == "peptides":
                i["sha256"] = sha

    edit_manifest(out, fn)
