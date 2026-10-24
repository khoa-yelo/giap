"""Builders for the M10-06 tests (docs/specs/M10-06-te-real-and-te7.md).

Everything is synthetic or hand-made: genes / variants / peptide frames, two hand-made run trees
for TE7, ``outputs.sha256`` lists, converted peptide tables, a temporary git repository for the
code-equivalence check, and the M10-05 synthetic pilot (``helpers.pilot_metrics_m10``) for the
CLI. No network, nothing from ``$GIAP_ROOT``, nothing derived from a DIA-NN report.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib
import json
import os
import subprocess
import sys
from collections.abc import Callable, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
TE_REAL_PY = REPO / "validation" / "pilot" / "te_real.py"
KNOWN_COLUMNS = REPO / "validation" / "pilot" / "te7_known_columns.tsv"
MINI_CONTAMINANTS = REPO / "tests" / "data" / "contaminants_mini.faa"
NISSLE_JSON = REPO / "tests/data/recorded/bakta/1.12.2-full/nissle" / ("nissle1917.json.gz")

CAUSE_CLASSES = ("pipeline_bug", "annotation", "search", "biology", "definition")


def load_te_real() -> Any:
    """``validation.pilot.te_real``; fails the calling test (not an error) when absent."""
    import pytest

    if not TE_REAL_PY.is_file():
        pytest.fail(f"{TE_REAL_PY.relative_to(REPO)} does not exist")
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return importlib.import_module("validation.pilot.te_real")


def run_te_real(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    """``python validation/pilot/te_real.py ARGS`` with ``PYTHONPATH=$REPO/src:$REPO``."""
    import pytest

    if not TE_REAL_PY.is_file():
        pytest.fail(f"{TE_REAL_PY.relative_to(REPO)} does not exist")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(REPO / "src"), str(REPO), str(REPO / "tests" / "sim")])
    return subprocess.run(
        [sys.executable, str(TE_REAL_PY), *args],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(cwd) if cwd else None,
        check=False,
    )


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def codes(warnings: Iterable[Any]) -> list[str]:
    """Warning codes: a warning is a string or a dict with ``code``."""
    return [w if isinstance(w, str) else str(w.get("code")) for w in warnings]


def tag_set(x: Any) -> set[str]:
    """Locus tags of a ``ribosomal_set`` result (a frame with ``locus_tag`` or an iterable)."""
    if isinstance(x, pd.DataFrame):
        return set(x["locus_tag"].astype(str))
    return {str(t) for t in x}


# ------------------------------------------------------------------ genes / variants frames

GENE_DEFAULTS: dict[str, Any] = {
    "gene": None,
    "product": "hypothetical protein",
    "feature_type": "cds",
    "contig": "chr1",
    "strand": "+",
    "crosses_origin": False,
    "replicon_type": "chromosome",
    "region_ids": None,
    "detection_status": "not_detected",
    "detected": False,
    "abundance_log2": None,
    "abundance_imputed": None,
    "abundance_type": "ibaq",
    "abundance_scale": "molar",
    "top3_log2": None,
    "n_peptides_unique": 0,
    "n_observable_peptides": 5,
}


def grow(locus_tag: str, **kw: Any) -> dict[str, Any]:
    """One ``genes.tsv``-like row; ``measured=<x>`` sets a detected, measured gene at x."""
    row = {"locus_tag": locus_tag, **GENE_DEFAULTS}
    if "measured" in kw:
        v = kw.pop("measured")
        row.update(
            detection_status="detected",
            detected=True,
            abundance_log2=float(v),
            abundance_imputed=False,
            n_peptides_unique=3,
        )
    row.update(kw)
    return row


def gframe(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    cols = ["locus_tag", *GENE_DEFAULTS]
    return pd.DataFrame([{c: r.get(c) for c in cols} for r in rows], columns=cols)


VARIANT_COLS = (
    "variant_id",
    "contig",
    "gene",
    "aa_pos",
    "ref_aa",
    "alt_aa",
    "alt_type",
    "call",
    "confidence",
    "n_ref",
    "n_alt",
    "class_error_rate",
    "source",
    "flags",
)


def vrow(
    vid: str,
    gene: str,
    aa_pos: int,
    ref: str,
    alt: str,
    conf: str,
    flags: str | None,
    *,
    n_alt: int = 1,
    alt_type: str = "substitution",
    call: str = "ALT",
) -> dict[str, Any]:
    return {
        "variant_id": vid,
        "contig": "chr1",
        "gene": gene,
        "aa_pos": aa_pos,
        "ref_aa": ref,
        "alt_aa": alt,
        "alt_type": alt_type,
        "call": call,
        "confidence": conf,
        "n_ref": 0,
        "n_alt": n_alt,
        "class_error_rate": "0.001",
        "source": "database",
        "flags": flags,
    }


def vframe(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(
        [{c: r.get(c) for c in VARIANT_COLS} for r in rows], columns=list(VARIANT_COLS)
    )


def ds_frame(id_mixups: float | None) -> pd.DataFrame:
    """A ``detection_summary.tsv``-like frame with the QC ``id_mixups`` row (NA when None)."""
    return pd.DataFrame(
        [
            {
                "section": "qc",
                "group_type": "run",
                "group": "all",
                "metric": "id_mixups",
                "value": float("nan") if id_mixups is None else float(id_mixups),
                "ci_low": float("nan"),
                "ci_high": float("nan"),
                "n": None if id_mixups is None else 10,
                "flags": "requires_variants" if id_mixups is None else None,
            }
        ]
    )


# Examples 1 and 2 (spec "Examples")
EX1_GENES = [
    grow("JMJJOG_00010", gene="thrA", contig="chr1", replicon_type="chromosome"),
    grow("JMJJOG_00020", gene="fdnG", contig="chr1", replicon_type="chromosome"),
    grow("JMJJOG_00030", gene="yaaA", contig="chr1", replicon_type="chromosome"),
    grow(
        "JMJJOG_04600",
        gene="pmA",
        contig="pMUT1",
        strand="-",
        replicon_type="plasmid",
        region_ids="prophage_2",
    ),
]


def ex1_variants() -> list[dict[str, Any]]:
    return [
        vrow("v1", "JMJJOG_00010", 1, "V", "M", "high", "start_codon;database_derived"),
        vrow("v2", "JMJJOG_00020", 196, "U", "*", "high", "stop_codon"),
        vrow("v3", "JMJJOG_00030", 41, "N", "D", "candidate", "database_derived"),
    ]


def ex2_extra() -> dict[str, Any]:
    return vrow(
        "v4", "JMJJOG_04600", 57, "N", "D", "high", "database_derived;mass_equivalent", n_alt=2
    )


# ------------------------------------------------------------------ peptide frames (TE5)


def prow(pep: str, ids: str, *, q: float | None = 0.001, decoy: bool = False) -> dict[str, Any]:
    return {
        "peptide": pep,
        "modified_peptide": pep,
        "protein_ids": ids,
        "q_value": q,
        "pep": 0.001,
        "is_decoy": decoy,
        "intensity": 1000.0,
        "run": "WT_1",
        "charge": 2,
        "source": "database",
        "search_db": "same",
    }


def pframe(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(list(rows))


def hits_frame(rows: Sequence[tuple[str, str]]) -> pd.DataFrame:
    """``peptide_hits.tsv``-like rows ``(peptide, route)`` (every row is an accepted placement)."""
    return pd.DataFrame(
        [
            {
                "peptide": p,
                "route": r,
                "locus_class": "cds_inframe",
                "gene": "JMJJOG_00010",
                "contig": "chr1",
                "start": 1,
                "end": 3 * len(p),
                "strand": "+",
            }
            for p, r in rows
        ],
        columns=["peptide", "route", "locus_class", "gene", "contig", "start", "end", "strand"],
    )


def unmapped_frame(rows: Sequence[tuple[str, str, str | None]]) -> pd.DataFrame:
    """``unmapped.tsv``-like peptide rows ``(peptide, label, contaminant_match)``."""
    return pd.DataFrame(
        [
            {
                "sequence_id": p,
                "sequence": p,
                "kind": "peptide",
                "label": lab,
                "reason": "no_hit",
                "contaminant_match": m,
            }
            for p, lab, m in rows
        ],
        columns=["sequence_id", "sequence", "kind", "label", "reason", "contaminant_match"],
    )


EX6_FAMILIES = {
    "trypsin": {"Cont_P00761"},
    "keratin": {"Cont_P04264", "Cont_P13645"},
    "keratin_associated": {"Cont_Q9BYQ7"},
}
ALBUMIN_PEPTIDE = "LVNELTEFAK"


def ex6_inputs(with_tolerant: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Example 6 (TE5, A2): P, peptide_hits, unmapped."""
    P = [
        prow("VATVSLPR", "Cont_P00761"),
        prow("SLNNQFASFIDK", "Cont_P04264"),
        prow("LSSPATLNSR", "Cont_P00761;JMJJOG_01234"),
        prow(ALBUMIN_PEPTIDE, "Cont_P02769"),
    ]
    hits = [("LSSPATLNSR", "exact")]
    unm: list[tuple[str, str, str | None]] = [
        ("VATVSLPR", "contaminant", "Cont_P00761"),
        ("SLNNQFASFIDK", "contaminant", "Cont_P13645"),
        (ALBUMIN_PEPTIDE, "contaminant", "Cont_P02769"),
    ]
    if with_tolerant:
        P.append(prow("LASYLDK", "Cont_P13645"))
        hits.append(("LASYLDK", "tolerant"))
    return pframe(P), hits_frame(hits), unmapped_frame(unm)


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


def git(repo: Path, *args: str) -> str:
    env = {
        **os.environ,
        "GIT_AUTHOR_NAME": "t",
        "GIT_AUTHOR_EMAIL": "t@example.org",
        "GIT_COMMITTER_NAME": "t",
        "GIT_COMMITTER_EMAIL": "t@example.org",
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_NOSYSTEM": "1",
    }
    res = subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, env=env, check=True
    )
    return res.stdout.strip()


def make_repo(root: Path) -> dict[str, str]:
    """A temporary repository with three commits: ``base``, ``docs`` (a change under ``docs/``
    only) and ``src`` (a change under ``src/``). Returns name -> commit SHA."""
    root.mkdir(parents=True, exist_ok=True)
    git(root, "init", "-q")
    for rel, text in (
        ("src/giap/a.py", "x = 1\n"),
        ("docs/notes.md", "a\n"),
        ("workflow/main.nf", "// nf\n"),
        ("pyproject.toml", "[project]\n"),
        ("uv.lock", "lock\n"),
        ("scripts/slurm/giap_run.sbatch", "#!/bin/bash\n"),
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
    (root / "src" / "giap" / "a.py").write_text("x = 2\n", encoding="utf-8")
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
            "command_line": ["giap", "run", "--outdir", f"/runs/{outdir_label}"],
        }
    )
    man.update(manifest_over)
    path.write_text(json.dumps(man, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    write_trace(out, hosts)
    return out


def edit_manifest(out: Path, fn: Callable[[dict[str, Any]], None]) -> None:
    path = out / "run_manifest.json"
    man = json.loads(path.read_text(encoding="utf-8"))
    fn(man)
    path.write_text(json.dumps(man, sort_keys=True, indent=2) + "\n", encoding="utf-8")


# ------------------------------------------------------------------ DIA-NN layer inputs

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
    from giap.io.tsv import format_value

    lines = ["\t".join(frame.columns)]
    lines += ["\t".join(format_value(v) for v in row) for row in frame.itertuples(index=False)]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gz(("\n".join(lines) + "\n").encode("utf-8")))
    return path


# ------------------------------------------------------------------ the synthetic pilot for the CLI

RIBO_N = 10
TE4_VALUES = {"rplL": 22.0, "rplJ": 20.1, "atpA": 19.0, "atpD": 19.6}


def write_contaminants(path: Path) -> Path:
    """``contaminants_mini.faa`` with the ``Cont_`` accessions of the real library (M10-03)."""
    text = MINI_CONTAMINANTS.read_text(encoding="utf-8").replace(">sp|P", ">sp|Cont_P")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def te_genes_rows(world: Any, *, ribo_missed: int = 0) -> list[dict[str, Any]]:
    """A2 ``genes.tsv`` rows: every coding gene detected and measured on a molar iBAQ scale; the
    first ``RIBO_N`` cds genes are 50S r-proteins at 25 (``ribo_missed`` of them not detected);
    the next four are rplL, rplJ, atpA, atpD at example 5's values."""
    from helpers.pilot_metrics_m10 import gene_row

    coding = sorted((g for g in world.genes if g.feature_type == "cds"), key=lambda g: g.locus_tag)
    special = {g.locus_tag: i for i, g in enumerate(coding[: RIBO_N + 4])}
    names = list(TE4_VALUES)
    rows = []
    for k, g in enumerate(sorted(world.genes, key=lambda x: x.locus_tag)):
        r = gene_row(g, status="detected", regions=world.regions)
        r.update(
            abundance_log2=float(5 + k % 15),
            abundance_imputed=False,
            abundance_type="ibaq",
            abundance_scale="molar",
            top3_log2=float(4 + k % 15),
            n_peptides_unique=2,
            n_observable_peptides=6,
            gene=None,
            product="hypothetical protein",
        )
        i = special.get(g.locus_tag)
        if i is not None and i < RIBO_N:
            r.update(product=f"50S ribosomal protein L{i + 1}", abundance_log2=25.0)
            if i < ribo_missed:
                r.update(
                    detection_status="not_detected",
                    detected=False,
                    abundance_log2=None,
                    abundance_imputed=None,
                )
        elif i is not None:
            name = names[i - RIBO_N]
            r.update(gene=name, abundance_log2=TE4_VALUES[name])
        rows.append(r)
    return rows


def te_variant_rows(world: Any, *, violation: bool = False) -> list[dict[str, Any]]:
    coding = sorted((g for g in world.genes if g.feature_type == "cds"), key=lambda g: g.locus_tag)
    rows = [
        vrow("v1", coding[20].locus_tag, 1, "V", "M", "high", "start_codon;database_derived"),
        vrow("v2", coding[21].locus_tag, 30, "U", "*", "high", "stop_codon"),
        vrow("v3", coding[22].locus_tag, 41, "N", "D", "candidate", "database_derived"),
    ]
    if violation:
        rows.append(
            vrow(
                "v4",
                coding[23].locus_tag,
                57,
                "N",
                "D",
                "high",
                "database_derived;mass_equivalent",
                n_alt=2,
            )
        )
    for r in rows:
        r.update(
            codon_start=1, strand="+", ref_codon="GTG", n_ref=0, nt_changes=1, class_error_n=10
        )
    return rows


def build_te_pilot(
    root: Path,
    world: Any,
    *,
    violation: bool = False,
    ribo_missed: int = 0,
    variants: bool = True,
    id_mixups: float | None = 2.0,
    a1_cont_unmapped: bool = False,
) -> tuple[Any, Path]:
    """The M10-05 synthetic pilot plus what TE2-TE5 read; returns (Pilot, contaminants FASTA)."""
    from helpers.pilot_metrics_m10 import (
        CONT_PEPTIDE,
        N_PEPTIDES,
        build_pilot,
        unmapped_row,
        write_table,
    )

    p = build_pilot(root, world)
    write_table(
        p.out_a2 / "tables" / "genes.tsv", "genes", te_genes_rows(world, ribo_missed=ribo_missed)
    )
    if variants:
        write_table(
            p.out_a2 / "tables" / "variants.tsv",
            "variants",
            te_variant_rows(world, violation=violation),
        )
    qc = {"section": "qc", "group_type": "run", "group": "all"}
    write_table(
        p.out_a2 / "tables" / "detection_summary.tsv",
        "detection_summary",
        [
            {**qc, "metric": "genome_inconsistent_rate", "value": None, "n": None},
            {
                **qc,
                "metric": "id_mixups",
                "value": id_mixups,
                "n": None if id_mixups is None else 4,
                "flags": "requires_variants" if id_mixups is None else None,
            },
        ],
    )
    if a1_cont_unmapped:
        unm = [unmapped_row(pep) for _, pep in p.peptides[2:3]]
        unm.append(unmapped_row(CONT_PEPTIDE, label="unmapped", reason="no_hit"))
        write_table(p.out_a1 / "tables" / "unmapped.tsv", "unmapped", unm)
    assert N_PEPTIDES >= 20
    return p, write_contaminants(root / "db" / "contaminants" / "contaminants.fasta")


def te_compute_args(p: Any, contam: Path, out: Path) -> list[str]:
    return [
        "compute",
        "--a2",
        str(p.out_a2),
        "--a1",
        str(p.out_a1),
        "--peptides-a1",
        str(p.pep_a1),
        "--peptides-a2",
        str(p.pep_a2),
        "--contaminants",
        str(contam),
        "--out",
        str(out),
    ]
