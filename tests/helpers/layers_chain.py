"""M2-20 test inputs: every ``context_layers`` input produced by the merged steps (test author).

REAL-CONTRACT RULE: nothing here invents a table shape. ``build_chain`` runs the merged producer
steps through ``gippa.steps.registry.run_step`` on gippasim ``SIM100K`` (the genome of every
sim100k recording):

* ``parse_bakta`` (M2-04) on the simulator's Bakta JSON + genome;
* ``uniop_input`` + ``operon_pairs`` (M2-12) with the **recorded** ``uniop.pred``
  (``tests/data/recorded/uniop/3ed14dfc/sim100k``, recorded on exactly this ``uniop_input``);
* ``codon_usage`` (M2-17);
* ``terminator_layer`` (M2-15) with the **recorded** TerminatorNet TSV and TransTermHP ``.tt``
  of sim100k (default tool TerminatorNet + cross-check);
* ``region_layer`` (M2-19) with a geNOMAD directory and an antiSMASH JSON written in the recorded
  formats (headers / JSON skeleton copied from ``tests/data/recorded/{genomad,antismash}``): one
  provirus on ``plasmid_1`` (so genes lie in ``prophage_1`` and ``plasmid_1``) and one antiSMASH
  area on ``chr1``;
* ``promoter_windows`` + ``promoter_layer`` (M2-13) with a Promotech ``sequences_predictions.csv``
  in the recorded format (``CHROM<TAB>SEQ<TAB>PRED``, one row per ``w<index>``), scores a fixed
  full-precision function of the index;
* ``tx_layer`` (M2-14) with Promoter Calculator CSVs in the recorded 27-column format
  (``helpers.promoter_calculator_fixtures.write_pc_csv``), hits for every third TU, header-only
  CSVs (what the tool writes when it finds nothing) for the rest;
* ``ostir_input`` + ``rbs_layer`` (M2-16) with an ``ostir.csv`` in the recorded header, one row
  at the requested start (plus a second start, as OSTIR reports several) for every gene except
  every seventh one (those genes have no ``rbs`` row: edge case "missing from rbs -> NA").

``precomputed_operon_pairs`` reruns ``operon_pairs`` with an ``annotation.precomputed.operon_pairs``
TSV (M2-02 format) to pin the probabilities of spec examples 2 and 3.
"""

from __future__ import annotations

import csv
import importlib
import json
import shutil
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest

TESTS = Path(__file__).resolve().parents[1]
REC = TESTS / "data" / "recorded"
UNIOP_REC = REC / "uniop" / "3ed14dfc" / "sim100k"
TN_REC = REC / "terminators" / "terminatornet" / "sim100k" / "terminatornet.tsv"
TT_REC = REC / "terminators" / "transtermhp" / "sim100k" / "transtermhp.tt"
GENOMAD_REC = REC / "genomad" / "1.12.0" / "sim100k" / "genome_summary"
ANTISMASH_REC = REC / "antismash" / "8.0.4" / "sim100k" / "bakta.json"
OSTIR_HEADER = (REC / "headers" / "ostir.ostir.csv.txt").read_text().splitlines()[0]
PROMOTECH_HEADER = "CHROM\tSEQ\tPRED"
#: Promotech writes full-precision probabilities (e.g. 0.1773268671071446, recorded run); the
#: fixture scores carry more than the 3 decimals a BED score keeps (lead decision 12).
PROMOTECH_DENOM = 1_000_003
TAX = "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;"
PROVIRUS_PLASMID = (1_001, 6_000)  # 1-based closed on plasmid_1 (8 000 bp), as geNOMAD writes it
BGC_AREA = (20_000, 30_000)  # antiSMASH JSON: 0-based half-open on chr1
NO_RBS_EVERY = 7
PARTIAL_LOCUS = "GSIM_00120"
RIBOSOMAL_LOCI = tuple(
    f"GSIM_{n:05d}" for n in range(295, 345, 5)
)  # every 7th gene (genes.parquet order) has no OSTIR row

PARSE_OUT = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}
#: context_layers keys -> file in the chain directory (the work-directory names, M2-index).
LAYER_FILES = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "pairs_raw": "pairs_raw.parquet",
    "tus": "tus.parquet",
    "promoter_windows": "promoter_windows.parquet",
    "tx_rates": "tx_rates.parquet",
    "terminators": "terminators.parquet",
    "terminators_crosscheck": "terminators_crosscheck.tsv",
    "rbs": "rbs.parquet",
    "codon": "codon.parquet",
    "regions": "regions.parquet",
    "region_genes": "region_genes.parquet",
}
#: lead decision 4 (gc_local_dev): the genome FASTA is a required context_layers input.
GENOME_KEY = "genome"
#: §7.1.8 / §9.1 files in ``annotation/`` (``id_map.tsv`` only with a user GFF3/GBFF).
ANNOTATION_FILES = (
    "genes.gff3", "regions.tsv", "regions.bed", "promoters.bed", "terminators.bed",
    "terminators_crosscheck.tsv", "operon_pairs.tsv", "tus.tsv",
)  # fmt: skip


def fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    """A clean registry with every ``gippa.steps`` module (re)registered."""
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"gippa.steps.{n}"))
    return reg


@contextmanager
def registry() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield fresh_steps(mp)


@dataclass
class Chain:
    """Producer outputs of SIM100K; ``inputs`` = every context_layers key -> path."""

    root: Path
    sim: Any
    inputs: dict[str, Path] = field(default_factory=dict)

    def frame(self, key: str) -> pd.DataFrame:
        return pd.read_parquet(self.inputs[key])

    def subset(self, *drop: str) -> dict[str, Path]:
        return {k: v for k, v in self.inputs.items() if k not in drop}


def _run(
    reg: ModuleType, name: str, cfg: Any, inputs: Mapping[str, Any], outs: Mapping[str, Path]
) -> None:
    reg.run_step(name, cfg, dict(inputs), dict(outs))


def genomad_dir(root: Path) -> Path:
    """geNOMAD output (recorded layout, prefix ``genome``): plasmid_1 a plasmid with a provirus."""
    d = root / "genomad" / "genome_summary"
    d.mkdir(parents=True)

    def header(name: str) -> str:
        return (GENOMAD_REC / f"genome_{name}.tsv").read_text().splitlines()[0]

    a, b = PROVIRUS_PLASMID
    name = f"plasmid_1|provirus_{a}_{b}"
    (d / "genome_virus_summary.tsv").write_text(
        header("virus_summary") + "\n"
        f"{name}\t{b - a + 1}\tProvirus\t{a}-{b}\t4\t11\t0.9312\tNA\t2\t9.5\t{TAX}\n"
    )
    (d / "genome_plasmid_summary.tsv").write_text(
        header("plasmid_summary") + "\n"
        "plasmid_1\t8000\tNo terminal repeats\t6\t11\t0.9871\tNA\t1\t3.2\tNA\tNA\n"
    )
    (d / "genome_virus_genes.tsv").write_text(header("virus_genes") + "\n")
    return root / "genomad"


def antismash_dir(root: Path, genes: pd.DataFrame) -> Path:
    """antiSMASH output: the recorded sim100k JSON with one NRPS area on chr1."""
    a, b = BGC_AREA
    inside = genes[(genes["contig"] == "chr1") & (genes["start"] >= a) & (genes["end"] <= b)]
    assert len(inside), "SIM100K has genes in the antiSMASH test area"
    feats: list[dict[str, Any]] = [
        {"type": "region", "location": f"[{a}:{b}]", "qualifiers": {"product": ["NRPS"]}}
    ]
    for k, g in enumerate(inside.itertuples(index=False)):
        q: dict[str, Any] = {"locus_tag": [g.locus_tag]}
        if k == 0:
            q["gene_kind"] = ["biosynthetic"]
        feats.append(
            {"location": f"[{g.start}:{g.end}]({g.strand})", "type": "CDS", "qualifiers": q}
        )
    data = json.loads(ANTISMASH_REC.read_text())
    rec = next(r for r in data["records"] if r["id"] == "chr1")
    rec["areas"] = [{"start": a, "end": b, "products": ["NRPS"]}]
    rec["features"] = feats
    d = root / "antismash"
    d.mkdir(parents=True)
    (d / "bakta.json").write_text(json.dumps(data))
    (d / "chr1.region001.gbk").write_text("LOCUS chr1\n//\n")
    return d


def promotech_dir(root: Path, windows: Path) -> Path:
    """``promotech_000/sequences_predictions.csv`` for every ``w<index>`` of the window chunks."""
    out = root / "promotech_000"
    out.mkdir(parents=True)
    lines = [PROMOTECH_HEADER]
    for fa in sorted(windows.glob("windows_*.fa")):
        recs = fa.read_text().split(">")[1:]
        for rec in recs:
            head, seq = rec.split("\n", 1)
            idx = int(head.strip()[1:])
            score = ((idx * 7919) % PROMOTECH_DENOM) / PROMOTECH_DENOM
            lines.append(f"{head.strip()}\t{seq.replace(chr(10), '')}\t{score}")
    (out / "sequences_predictions.csv").write_text("\n".join(lines) + "\n")
    return out


def pc_csvs(root: Path, tus: pd.DataFrame) -> list[Path]:
    """One Promoter Calculator CSV per TU (stem = TU id); hits for every third TU."""
    from helpers.promoter_calculator_fixtures import header_only, write_pc_csv

    paths = []
    for k, t in enumerate(tus.itertuples(index=False)):
        n = int(t.upstream_end - t.upstream_start)
        path = root / "pc" / f"{t.tu_id}.csv"
        if k % 3 == 0 and n >= 20:
            hits = [(max(1, n - 60), "+", 1.5 + k), (max(1, n - 120), "+", 0.25 + k / 10)]
            write_pc_csv(path, hits)
        else:
            header_only(path)
        paths.append(path)
    return paths


def ostir_csv(root: Path, ostir_input: Path) -> Path:
    """``ostir.csv`` (recorded header) for every gene except every ``NO_RBS_EVERY``-th one."""
    from gippa.annotate.ostir import read_ostir_input

    rows = read_ostir_input(ostir_input)
    path = root / "ostir.csv"
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(OSTIR_HEADER.split(","))
        for i, r in enumerate(rows):
            if i % NO_RBS_EVERY == NO_RBS_EVERY - 1:
                continue
            codon = r.seq[r.start - 1 : r.start + 2]
            if r.start > 10:  # a second start reported by OSTIR (ignored by parse_ostir)
                w.writerow([r.id, "ATG", r.start - 9, "3.5", "4", "1.0", "-2.0", "-9.0", "0.1",
                            "0.0", "-1.194"])  # fmt: skip
            w.writerow([r.id, codon, r.start, f"{100.0 + 3 * i:.4f}", 5 + i % 7,
                        "0.5", f"{-12.0 + i % 5:.3f}", f"{-14.0 + i % 3:.1f}",
                        f"{0.1 * (i % 4):.4f}", f"{-0.2 * (i % 2):.1f}", "-1.194"])  # fmt: skip
    return path


def bakta_json_variant(src: Path, dst: Path) -> Path:
    """The simulator's Bakta JSON with Bakta's own fields set on a few CDS: ``truncated`` on
    ``PARTIAL_LOCUS`` (a partial gene) and ``50S ribosomal protein`` products on
    ``RIBOSOMAL_LOCI`` (the CAI reference set, gene flag ``cai_reference``)."""
    data = json.loads(src.read_text())
    seen = set()
    for f in data["features"]:
        if f.get("locus") == PARTIAL_LOCUS:
            f["truncated"] = "5'"
            seen.add(f["locus"])
        if f.get("locus") in RIBOSOMAL_LOCI:
            f["product"] = f"50S ribosomal protein L{RIBOSOMAL_LOCI.index(f['locus']) + 1}"
            seen.add(f["locus"])
    assert seen == {PARTIAL_LOCUS, *RIBOSOMAL_LOCI}, seen
    dst.write_text(json.dumps(data))
    return dst


def build_chain(root: Path) -> Chain:
    """Run every producer on SIM100K into ``root``; returns the context_layers inputs."""
    import gippasim
    from gippa.config import GippaConfig

    sim = gippasim.simulate(gippasim.SIM100K, gippasim.SIM100K_SEED, root / "sim")
    w = root / "work"
    w.mkdir(parents=True)
    bakta = bakta_json_variant(sim.bakta_json, root / "bakta.json")
    cfg = GippaConfig(seed=1)
    with registry() as reg:
        _run(reg, "parse_bakta", cfg, {"bakta_json": bakta, "genome": sim.genome_fasta},
             {k: w / v for k, v in PARSE_OUT.items()})  # fmt: skip
        core = {k: w / PARSE_OUT[k] for k in ("genes", "replicons")}
        _run(reg, "uniop_input", cfg, core,
             {"uniop_cds": w / "uniop_cds.gff", "uniop_map": w / "uniop_map.tsv"})  # fmt: skip
        shutil.copyfile(UNIOP_REC / "uniop.pred", w / "uniop.pred")
        _run(reg, "operon_pairs", cfg,
             {**core, "genome": sim.genome_fasta, "uniop_cds": w / "uniop_cds.gff",
              "uniop_map": w / "uniop_map.tsv", "uniop_pred": w / "uniop.pred"},
             {"pairs_raw": w / "pairs_raw.parquet", "tus": w / "tus.parquet",
              "tu_upstream": w / "tu_upstream", "operon_meta": w / "operon_meta.json"})  # fmt: skip
        _run(reg, "codon_usage", cfg,
             {"genes": w / "genes.parquet", "trna_counts": w / "trna_counts.tsv"},
             {"codon": w / "codon.parquet"})  # fmt: skip
        _run(reg, "terminator_layer", cfg,
             {"replicons": core["replicons"], "terminatornet": TN_REC, "transtermhp": TT_REC},
             {"terminators": w / "terminators.parquet",
              "terminators_crosscheck": w / "terminators_crosscheck.tsv"})  # fmt: skip
        genes = pd.read_parquet(w / "genes.parquet")
        tools = root / "tools"
        _run(reg, "region_layer", cfg,
             {**core, "features": w / "features.parquet", "genomad": genomad_dir(tools),
              "antismash": antismash_dir(tools, genes)},
             {"regions": w / "regions.parquet",
              "region_genes": w / "region_genes.parquet"})  # fmt: skip
        _run(reg, "promoter_windows", cfg, {**core, "genome": sim.genome_fasta},
             {"windows": w / "windows", "windows_map": w / "windows_map.parquet"})  # fmt: skip
        _run(reg, "promoter_layer", cfg,
             {"windows_map": w / "windows_map.parquet",
              "promotech": [promotech_dir(tools, w / "windows")]},
             {"promoter_windows": w / "promoter_windows.parquet"})  # fmt: skip
        pcs = pc_csvs(tools, pd.read_parquet(w / "tus.parquet"))
        _run(reg, "tx_layer", cfg,
             {"pc": pcs, "tus": w / "tus.parquet", "replicons": core["replicons"]},
             {"tx_rates": w / "tx_rates.parquet"})  # fmt: skip
        _run(reg, "ostir_input", cfg,
             {**core, "genome": sim.genome_fasta, "annotation_meta": w / "annotation_meta.json"},
             {"ostir_input": w / "ostir_input.csv"})  # fmt: skip
        _run(reg, "rbs_layer", cfg,
             {"ostir": ostir_csv(tools, w / "ostir_input.csv"),
              "ostir_input": w / "ostir_input.csv"},
             {"rbs": w / "rbs.parquet"})  # fmt: skip
    inputs = {k: w / v for k, v in LAYER_FILES.items()}
    inputs[GENOME_KEY] = Path(sim.genome_fasta)
    return Chain(root, sim, inputs)


def genome_seqs(chain: Chain) -> dict[str, str]:
    """contig -> sequence of the chain's genome FASTA (the ``genome`` input)."""
    from gippa.io.fasta import read_genome_fasta

    return {r.id: r.sequence for r in read_genome_fasta(chain.inputs[GENOME_KEY])}


def precomputed_operon_pairs(
    chain: Chain, out: Path, overrides: Mapping[tuple[str, str], float]
) -> tuple[Path, Path]:
    """Rerun ``operon_pairs`` with the chain's pairs as an ``annotation.precomputed.operon_pairs``
    TSV (M2-02 format), ``overrides`` replacing p of given (gene_a, gene_b) pairs. Returns the new
    (pairs_raw, tus) parquet paths."""
    from gippa.config import GippaConfig
    from gippa.io.tsv import write_tsv

    pairs = chain.frame("pairs_raw")
    rows = [
        (a, b, float(overrides.get((a, b), p)), s)
        for a, b, p, s in pairs[["gene_a", "gene_b", "uniop_p", "source"]].itertuples(index=False)
    ]
    missing = set(overrides) - {(a, b) for a, b, _, _ in rows}
    assert not missing, f"not adjacent same-strand pairs: {missing}"
    out.mkdir(parents=True, exist_ok=True)
    tsv = out / "operon_pairs.tsv"
    write_tsv(tsv, ["gene_a", "gene_b", "uniop_p", "source"], rows)
    cfg = GippaConfig.model_validate(
        {"seed": 1, "annotation": {"precomputed": {"operon_pairs": str(tsv)}}}
    )
    w = chain.inputs["genes"].parent
    with registry() as reg:
        _run(reg, "operon_pairs", cfg,
             {"genes": chain.inputs["genes"], "replicons": chain.inputs["replicons"],
              "genome": chain.sim.genome_fasta, "uniop_cds": w / "uniop_cds.gff",
              "uniop_map": w / "uniop_map.tsv"},
             {"pairs_raw": out / "pairs_raw.parquet", "tus": out / "tus.parquet",
              "tu_upstream": out / "tu_upstream",
              "operon_meta": out / "operon_meta.json"})  # fmt: skip
    return out / "pairs_raw.parquet", out / "tus.parquet"


def read_layers_frames(chain: Chain, **replace: Path | None) -> dict[str, pd.DataFrame | None]:
    """The ``gene_layers`` arguments as DataFrames read from the producer parquet files;
    ``replace`` swaps a key's file (``None`` = layer missing)."""
    paths: dict[str, Path | None] = {**chain.inputs, **replace}
    keys = ("genes", "tus", "pairs_raw", "promoter_windows", "tx_rates", "terminators", "rbs",
            "codon", "region_genes", "replicons")  # fmt: skip
    return {k: (None if paths[k] is None else pd.read_parquet(paths[k])) for k in keys}


def call_gene_layers(
    frames: Mapping[str, pd.DataFrame | None], genome: Mapping[str, str] | None = None
) -> pd.DataFrame:
    """``gene_layers`` with the spec's positional order (M2-20 public interface) and the
    keyword ``genome`` (contig -> sequence; lead decision 4, ``gc_local_dev``)."""
    from gippa.annotate.layers import gene_layers

    return gene_layers(
        frames["genes"], frames["tus"], frames["pairs_raw"], frames["promoter_windows"],
        frames["tx_rates"], frames["terminators"], frames["rbs"], frames["codon"],
        frames["region_genes"], frames["replicons"], genome=genome,
    )  # fmt: skip
