"""Build the M9-25 frozen snapshot fixtures (docs/specs/M9-25-snapshot-baselines.md item 1).

Usage: python scripts/dev/make_snapshot_fixture.py [--work DIR] [--dest DIR]

Simulates one giapsim genome (SIM100K's signals at 300 kb, a SNV-only strain, de novo peptides,
planted novel ORFs, synthetic GO/KEGG/COG xrefs), plants an abundance island, low prophage/BGC
detection and in-operon terminators, and runs ``giap run --executor local`` on it twice (run_b
on two thirds of the spectra) with every module and a synthetic test ``$GIAP_DB_DIR`` (never
the real one), then ``giap compare`` on the two runs. It writes, under ``DEST`` (default
``tests/data/snapshot``):

* ``sim_plots/OUT/``: the ``OUT/`` subset ``RunView`` reads, from the first run, with
  ``report.max_locus_plots: 2`` in ``resolved_config.yaml``; ``sim_plots/CMP/``: the compare
  tables and, per run, ``tables/genes.tsv`` and ``preflight/preflight.json``;
* ``sim_empty/``: the same tree with every TSV cut to its header and every Parquet to zero rows.

TSVs over 16 kB are gzipped in place (the readers detect gzip by magic bytes). Deterministic:
the same giapsim seed and run seeds give the same bytes (paths, time stamps and timings are
scrubbed). Details: tests/data/snapshot/PROVENANCE.md. Regenerate only in a PR that says why.
"""

from __future__ import annotations

import argparse
import dataclasses
import gzip
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / "src"), str(REPO / "tests"), str(REPO / "tests" / "sim")]

SIM_SEED = 20261014
RUN_SEEDS = {"run_a": 7, "run_b": 11}
GENOME_LENGTH = 300_000
NOVEL = (("intergenic", 6), ("cds_outframe", 6), ("antisense", 28))
SETS = (
    "report.max_locus_plots=2",
    "context.islands.detrend_span_frac=0.5",
    "determinants.permutations=19",
    "stats.bootstrap=199",
)
#: ``OUT/``-relative directories kept besides the ``RunView`` files
KEEP_DIRS = ("logs/steps", "plots/support")
#: ``RunView`` keys no plot or report page reads (left out to keep the fixture small)
DROP = ("variants_vcf", "pileup", "novel_orfs_faa", "novel_orfs_gff3", "annotation_patch",
        "gene_table", "genes_pre", "gene_context", "tu_folds", "detectability_summary",
        "tier_summary", "novel_class_fdr")  # fmt: skip
#: run_b keeps 2 of every 3 de novo peptide rows, so the two runs' detections differ (V16)
KEEP_EVERY = 3
#: the planted abundance island: ``ISLAND_GENES`` consecutive genes from ``ISLAND_FIRST`` (in
#: genome order) raised by ``ISLAND_SHIFT`` log units in the de novo peptides
ISLAND_FIRST, ISLAND_GENES, ISLAND_SHIFT = 120, 26, 12.0
#: planted in-operon terminators (confidence 1, above ``terminatornet_yes``) between the first
#: two genes of the first ``TERM_TUS`` multi-gene chr1 TUs after ``TERM_FROM`` (flagged pairs)
TERM_TUS, TERM_FROM = 4, 200_000
GZIP_ABOVE = 16_000
MODULES = "core,mobile,bgc,function"
#: ISO-8601 time stamps (events, manifest), fixed in the fixture so a rebuild gives the same bytes
STAMP = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d+)?(?:Z|[+-]\d\d:\d\d)?")
FIXED_STAMP = "2026-01-01T00:00:00+00:00"
#: machine-dependent numbers (step timings, free disk), set to 0 for the same reason
MEASURED = re.compile(r'"(cpu_time_s|wall_time_s|free_disk_bytes)": [-+.0-9eE]+')


def sim_config() -> Any:
    import giapsim
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.strain import StrainConfig

    return dataclasses.replace(
        giapsim.SIM100K,
        genome_length=GENOME_LENGTH,
        strain=StrainConfig(snv_rate_per_kb=1.0, indel_rate_per_kb=0.0),
        proteomics=ProteomicsConfig(),
        denovo=DenovoConfig(),
        n_novel_orfs=NOVEL,
        n_nterm_extensions=3,
        novel_orf_len_aa=(30, 80),
        plant_novel_rbs=True,
    )


def add_xrefs(src: Path, dest: Path) -> None:
    """``src`` (Bakta-like JSON) with synthetic COG/GO/KEGG ``db_xrefs`` on every CDS."""
    doc = json.loads(src.read_text(encoding="utf-8"))
    cds = [f for f in doc["features"] if f.get("type") == "cds"]
    for i, f in enumerate(cds):
        f["db_xrefs"] = [f"COG:COG{i % 40 + 1:04d}"]
        if i % 3:
            f["db_xrefs"] += [f"GO:{i % 30 + 1:07d}", f"KEGG:K{i % 25 + 1:05d}"]
    dest.write_text(json.dumps(doc), encoding="utf-8")


def write_db(root: Path, bakta: Path) -> Path:
    """The synthetic function DB for ``bakta``'s xrefs plus every preflight requirement."""
    import pandas as pd
    from helpers.function_db import ITEMS, PATHS, write_function_db
    from helpers.mapping_steps import write_test_db

    rows = [
        (f["locus"], *(x.split(":", 1) if not x.startswith("GO:") else ("GO", x)))
        for f in json.loads(bakta.read_text(encoding="utf-8"))["features"]
        for x in f.get("db_xrefs") or []
    ]
    write_function_db(root, pd.DataFrame(rows, columns=["locus_tag", "db", "id"]))
    write_test_db(root)
    doc = json.loads((root / "db_manifest.json").read_text(encoding="utf-8"))
    doc |= {n: {"version": "synthetic", "path": PATHS[n]} for n in ITEMS}
    (root / "db_manifest.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return root


def invoke(args: list[str]) -> None:
    from click.testing import CliRunner

    from giap.cli import cli

    res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    if res.exit_code != 0:
        raise SystemExit(f"giap {args[0]} exited {res.exit_code}:\n{res.output[-4000:]}")


def thin(src: Path, dest: Path) -> Path:
    """``src`` without every ``KEEP_EVERY``-th data row (the de novo and reference peptide
    files are paired row by row, so both lose the same spectra)."""
    lines = src.read_text(encoding="utf-8").splitlines(keepends=True)
    dest.write_text("".join(lines[:1] + [ln for i, ln in enumerate(lines[1:]) if i % KEEP_EVERY]),
                    encoding="utf-8")  # fmt: skip
    return dest


def _rows(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return lines[0].split("\t"), [ln.split("\t") for ln in lines[1:]]


def plant(sim: Any, work: Path) -> dict[str, Path]:
    """The de novo and the paired reference peptides with two planted signals, decided per
    spectrum from the reference peptide's source gene (``protein_ids``):

    * the island: intensities of the ``ISLAND_*`` chr1 CDS genes multiplied by
      ``2 ** ISLAND_SHIFT`` (the run's abundance is the peptides' iBAQ);
    * low detection: every spectrum of every second gene of each prophage and BGC region
      (``precomputed/regions.tsv``) dropped, so those regions' O - E is negative.
    """
    from helpers.m7_chain import denovo, reference_peptides

    feats = json.loads(Path(sim.reference_json).read_text(encoding="utf-8"))["features"]
    chr1 = [str(f["locus"]) for f in feats if f.get("type") == "cds" and f["sequence"] == "chr1"]
    island = set(chr1[ISLAND_FIRST : ISLAND_FIRST + ISLAND_GENES])
    _, regions = _rows(Path(sim.outdir) / "precomputed" / "regions.tsv")
    silent = {t for r in regions if r[1] in ("prophage", "bgc") for t in r[7].split(";")[::2]}
    head, ref = _rows(reference_peptides(sim))
    pid, run, spec = (head.index(c) for c in ("protein_ids", "run", "spectrum_id"))
    tags = {(r[run], r[spec]): set(r[pid].split(";")) for r in ref}
    out = {}
    for name, src in (("peptides", denovo(sim)), ("reference_peptides", reference_peptides(sim))):
        head, rows = _rows(src)
        run, spec, inten = (head.index(c) for c in ("run", "spectrum_id", "intensity"))
        kept = ["\t".join(head)]
        for r in rows:
            genes = tags.get((r[run], r[spec]), set())
            if genes & silent:
                continue
            if genes & island and r[inten] != "NA":
                r[inten] = f"{float(r[inten]) * 2**ISLAND_SHIFT:.6g}"
            kept.append("\t".join(r))
        out[name] = work / f"{name}_a.tsv"
        out[name].write_text("\n".join(kept) + "\n", encoding="utf-8")
    return out


def plant_terminators(sim: Any, dest: Path) -> Path:
    """``precomputed/terminators.tsv`` plus one confidence-1 call at the junction of the first two
    genes of each ``TERM_*`` TU (a concordant OP pair with a terminator between is flagged
    ``conflict_terminator_in_operon``)."""
    feats = json.loads(Path(sim.reference_json).read_text(encoding="utf-8"))["features"]
    ends = {str(f["locus"]): (int(f["start"]), int(f["stop"])) for f in feats if f["type"] == "cds"}
    _, tus = _rows(Path(sim.outdir) / "truth" / "tus.tsv")
    multi = [r for r in tus if r[1] == "chr1" and int(r[2]) > TERM_FROM and ";" in r[5]]
    src = Path(sim.outdir) / "precomputed" / "terminators.tsv"
    lines = src.read_text(encoding="utf-8").splitlines()
    for r in multi[:TERM_TUS]:
        up, down = sorted((ends[g] for g in r[5].split(";")[:2]), key=lambda se: se[0])
        mid = (up[1] + down[0]) // 2
        lines.append(f"chr1\t{mid - 10}\t{mid + 10}\t{r[4]}\t1.0\tTerminatorNet")
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return dest


def giap_run(sim: Any, work: Path, out: Path, seed: int, peptides: dict[str, Path]) -> Path:
    from helpers.m7_chain import run_config

    bakta = work / "reference_xrefs.json"
    cfg = run_config(sim, work / "run_config.yaml")
    terms = Path(sim.outdir).resolve() / "precomputed" / "terminators.tsv"
    text = cfg.read_text(encoding="utf-8").replace(str(Path(sim.reference_json).resolve()),
                                                   str(bakta.resolve()))  # fmt: skip
    text = text.replace(str(terms), str(plant_terminators(sim, work / "terminators.tsv")))
    cfg.write_text(text, encoding="utf-8")
    args = ["run", "--genome", str(sim.reference_fasta), "--peptides", str(peptides["peptides"]),
            "--reference-peptides", str(peptides["reference_peptides"]), "--abundance",
            str(sim.abundance), "--annotation", str(bakta), "--config", str(cfg),
            "--genome-relation", "related", "--modules", MODULES,
            "--outdir", str(out), "--seed", str(seed), "--executor", "local",
            "--threads", "2"]  # fmt: skip
    for s in SETS:
        args += ["--set", s]
    invoke(args)
    return out


def subset(out: Path, dest: Path) -> None:
    """Copy the files ``RunView`` reads (and ``KEEP_DIRS``) of ``out`` into ``dest``."""
    from giap.io.run_view import PATHS

    rels = {rel for k, rel in PATHS.items() if k not in DROP and (out / rel).is_file()}
    for d in KEEP_DIRS:
        rels |= {p.relative_to(out).as_posix() for p in (out / d).rglob("*") if p.is_file()}
    for rel in sorted(rels):
        (dest / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(out / rel, dest / rel)


def scrub(root: Path, work: Path) -> None:
    """Replace the scratch directory by ``WORK``, every time stamp by ``FIXED_STAMP`` and every
    ``MEASURED`` number by 0 in the text files under ``root``."""
    for p in sorted(root.rglob("*")):
        if p.is_file() and p.suffix in (".json", ".jsonl", ".yaml", ".tsv"):
            text = p.read_text(encoding="utf-8").replace(str(work.resolve()), "WORK")
            text = STAMP.sub(FIXED_STAMP, text.replace(str(work), "WORK"))
            p.write_text(MEASURED.sub(r'"\1": 0', text), encoding="utf-8")


def write_config(out: Path) -> None:
    """``out/resolved_config.yaml``: the defaults with the run's modules, relation and ``SETS``
    (no input paths; the plots read only these keys)."""
    import yaml

    from giap.config import resolve_config

    sets = [*SETS, f"modules=[{MODULES}]", "genome_relation=related"]
    cfg = resolve_config(overrides=sets).model_dump(mode="json")
    text = yaml.safe_dump(cfg, sort_keys=True)
    (out / "resolved_config.yaml").write_text(text, encoding="utf-8")


def compress(root: Path) -> None:
    """Gzip every TSV above ``GZIP_ABOVE`` bytes (``mtime=0``: reproducible): in place for the
    schema tables, as ``<name>.tsv.gz`` for the support files (read by suffix)."""
    for p in sorted(root.rglob("*.tsv")):
        raw = p.read_bytes()
        if len(raw) > GZIP_ABOVE and raw[:2] != b"\x1f\x8b":
            dest = p.with_name(p.name + ".gz") if "support" in p.parts else p
            p.unlink()
            dest.write_bytes(gzip.compress(raw, mtime=0))


def header_only(root: Path) -> None:
    """Cut every TSV under ``root`` (support files included, written as plain ``.tsv``) to its
    header and every Parquet to zero rows."""
    import pandas as pd

    for p in sorted(root.rglob("*")):
        if ".tsv" in p.name:
            raw = p.read_bytes()
            text = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode("utf-8")
            p.unlink()
            p.with_name(p.name.removesuffix(".gz")).write_text(
                text.split("\n", 1)[0] + "\n", encoding="utf-8"
            )
        elif p.suffix == ".parquet":
            pd.read_parquet(p).iloc[:0].to_parquet(p, index=False)


def build(work: Path, dest: Path) -> None:
    import giapsim

    sim = giapsim.simulate(sim_config(), SIM_SEED, work / "sim")
    add_xrefs(Path(sim.reference_json), work / "reference_xrefs.json")
    db = write_db(work / "db", work / "reference_xrefs.json")
    os.environ["GIAP_DB_DIR"] = str(db)
    planted = plant(sim, work)
    thinned = {k: thin(p, work / f"{k}_b.tsv") for k, p in planted.items()}
    peptides = {"run_a": planted, "run_b": thinned}
    runs = {lab: giap_run(sim, work, work / lab / "OUT", s, peptides[lab])
            for lab, s in RUN_SEEDS.items()}  # fmt: skip
    invoke(["compare", *map(str, runs.values()), "--labels", ",".join(runs),
            "--outdir", str(work / "CMP")])  # fmt: skip
    rich = dest / "sim_plots"
    shutil.rmtree(rich, ignore_errors=True)
    subset(runs["run_a"], rich / "OUT")
    for name in ("compare_detection", "compare_abundance", "compare_stability"):
        (rich / "CMP" / "tables").mkdir(parents=True, exist_ok=True)
        shutil.copyfile(work / "CMP" / "tables" / f"{name}.tsv",
                        rich / "CMP" / "tables" / f"{name}.tsv")  # fmt: skip
    for lab, out in runs.items():
        for rel in ("tables/genes.tsv", "preflight/preflight.json"):
            (rich / "CMP" / "runs" / lab / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(out / rel, rich / "CMP" / "runs" / lab / rel)
    write_config(rich / "OUT")
    scrub(rich, work)
    compress(rich)
    empty = dest / "sim_empty"
    shutil.rmtree(empty, ignore_errors=True)
    shutil.copytree(rich, empty)
    header_only(empty)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--work", type=Path, default=None, help="scratch directory (default: temp)")
    ap.add_argument("--dest", type=Path, default=REPO / "tests" / "data" / "snapshot")
    ns = ap.parse_args(argv)
    if ns.work is not None:
        ns.work.mkdir(parents=True, exist_ok=True)
        build(ns.work, ns.dest)
        return 0
    with tempfile.TemporaryDirectory() as tmp:
        build(Path(tmp), ns.dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
