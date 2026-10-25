"""Hand-built MVP-04 run and giapsim input directories for ``scripts/dev/mvp04_summary.py``
(docs/specs/MVP-04-ecoli-evidence-run.md, acceptance gates 1-16).

``build(tmp)`` writes a minimal published GIAP output directory and a giapsim input directory
(plus an identical second giapsim call for gate 9) in which **every gate passes and every value
lies in its expected range**. Each file carries only the columns the summary reads. Tests then
edit one file to move one gate. Nothing here imports giap or giapsim.

The numbers (``N`` below) are chosen so that every gate's value can be checked by hand:

* 4,000 protein genes ``K12_00000`` ... (step 5), the first 10 ``partial_excluded``: 3,990 eligible;
* genes 10..2529 have a true peptide (2,520, share 0.6316, gate 8);
* detected: genes 10..2389 (true, 2,380) and 2530..2549 (false, a protein null each): 2,400
  detected, realised FDP 20 / 2,400 = 0.00833 (gates 6, 7);
* genes 2550..2589 are decoy wins (gene q NA, decoy ``rev_<gene>`` at q 0.001), 2590..2594 have a
  decoy above ``q_value_max`` (gene q NA; not a win, ADR-0047) (gate 5);
* tier ``decoy_picked``, f = 1, q_floor = 1 / 2,600 (gates 2-4);
* 50 extra true peptides on genes 10..59 map by the tolerant route, 30 null peptides are
  unmapped: exact share 2,540 / 2,620 (gate 10);
* 60 truth ALT sites (50 with 2 true covering peptides), 60 matching ``high`` rows (gate 11);
* 120 truth ORFs per novel class (110 with a true peptide), 40 correct calls per class (gate 12);
* 2,000 two-gene TUs; 4 planted islands (lead decision 20): one high of 20 TUs (``HIGH``, index 0)
  and three lows of 8 TUs, all four called in their direction (gate 13; the lows are exempt, so
  calling them is allowed but not required). Every island TU is scored (no partial genes), so
  w = 20 / 8 and q = 1,590 / 3,990 undetected eligible genes;
* 3,999 adjacent pairs, 2,000 of them OP; 2,000 operons (gate 14);
* a report with content in the five module sections (gate 15); 5 h wall, 8 CPU-h (gate 16).
"""

from __future__ import annotations

import csv
import json
import shutil
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

CONTIG = "U00096.3"
N_GENES = 4000
N_PARTIAL = 10
TRUE_GENES = range(10, 2530)
DETECTED_TRUE = range(10, 2390)
FALSE_DETECTED = range(2530, 2550)
DECOY_WINS = range(2550, 2590)
DECOY_ABOVE_QMAX = range(2590, 2595)
TOLERANT_GENES = range(10, 60)
N_UNMAPPED = 30
N_ALT = 60
N_ALT_MULTI = 50
NOVEL_CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
N_NOVEL_TRUTH = 120
N_NOVEL_WITH_PEPTIDE = 110
N_NOVEL_CALLS = 40
N_TUS = N_GENES // 2
#: Planted islands in truth order: the high one first (``HIGH``), then the three lows.
ISLAND_FIRST_TUS = (700, 100, 400, 1300)
ISLAND_EFFECT_SD = (3.0, -3.0, -3.0, -3.0)
ISLAND_N_TUS = (20, 8, 8, 8)
HIGH = 0
LOWS = (1, 2, 3)
Q_FLOOR = 1 / 2600
REPORT_SECTIONS = ("detection", "genome-context", "operons", "variants", "novel-orfs")
ALL_SECTIONS = ("summary", "inputs-qc", "mapping", *REPORT_SECTIONS, "function",
                "determinants", "sensitivity", "methods", "provenance")  # fmt: skip
TRUTH_TABLES = ("peptides", "genes", "islands", "alt_sites", "novel_alt_sites", "novel_orfs")
#: ``mapping.tolerant.min_len`` in ``resolved_config.yaml`` (gate 11 recall; peptides are 8 aa)
TOLERANT_MIN_LEN = 7

_AA = "ACDEFGHILMNPQSTVWY"  # no K / R: every peptide ends in its one K


def locus(i: int) -> str:
    return f"K12_{i * 5:05d}"


def peptide(i: int) -> str:
    """A distinct tryptic-looking peptide per integer ``i`` (7 residues + K)."""
    s = ""
    for _ in range(7):
        s += _AA[i % len(_AA)]
        i //= len(_AA)
    return s + "K"


def gene_coords(i: int) -> tuple[int, int, str]:
    """1-based closed (start, end) and strand of gene ``i``; TU ``k`` = genes 2k, 2k+1."""
    start = 1 + i * 1000
    return start, start + 899, "+" if (i // 2) % 2 == 0 else "-"


def tu_id(i: int) -> str:
    return f"TU_{i // 2 + 1}"


def tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows([["NA" if v is None else v for v in r] for r in rows])
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def edit_tsv(path: Path, fn: Callable[[list[dict[str, str]]], list[dict[str, str]]]) -> None:
    """Rewrite ``path`` with ``fn(rows)``; the header is kept."""
    with path.open(newline="", encoding="utf-8") as fh:
        header = next(csv.reader(fh, delimiter="\t"))
    rows = fn(read_tsv(path))
    tsv(path, header, [[r.get(h, "NA") for h in header] for r in rows])


def append_tsv(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    edit_tsv(path, lambda old: [*old, *({k: str(v) for k, v in r.items()} for r in rows)])


def edit_json(path: Path, fn: Callable[[dict[str, Any]], None]) -> None:
    data = json.loads(path.read_text())
    fn(data)
    path.write_text(json.dumps(data, indent=2) + "\n")


@dataclass(frozen=True)
class Dirs:
    out: Path  # published GIAP outputs
    inputs: Path  # the giapsim call
    again: Path  # the second giapsim call with the same seed (gate 9)


# ---------------------------------------------------------------- inputs (giapsim)


def _novel_truth(k: int, cls: str) -> dict[str, Any]:
    """Truth ORF ``k`` of ``cls``: 1-based closed [start, end] including the stop codon."""
    plus = k % 2 == 0
    start = 5_000_000 + NOVEL_CLASSES.index(cls) * 100_000 + k * 600  # beyond the genes
    end = start + 299
    stop = end - 2 if plus else start  # lowest coordinate of the stop codon (M7-12 decision 16)
    return {
        "orf_id": f"NOVEL_{cls}_{k:04d}",
        "class": cls,
        "contig": CONTIG,
        "start": start,
        "end": end,
        "strand": "+" if plus else "-",
        "aa": "M" + "A" * 98,
        "stop": stop,
        "true_start": start if plus else end,
        "n_true_peptides": 1 if k < N_NOVEL_WITH_PEPTIDE else 0,
        "rbs_planted": "false",
    }


def _write_inputs(inp: Path) -> None:
    peps: list[list[Any]] = []
    truth: list[list[Any]] = []
    n = 0

    def add(
        pep: str, entry: str, decoy: bool, q: float, wrong: bool, origin: str, tag: str
    ) -> None:
        peps.append([pep, entry, q, q * 10, round(5 - q * 100, 6), "true" if decoy else "false",
                     "NA" if decoy else 1000, "run1"])  # fmt: skip
        if not decoy:
            # a CDS peptide covers its gene's truth ALT codon (gene start + 30; gate 11 recall)
            a = gene_coords(int(tag[4:]) // 5)[0] + 24 if origin == "cds" else "NA"
            b = a + 23 if origin == "cds" else "NA"
            truth.append([pep, "true" if wrong else "false", origin, tag, CONTIG, a, b])

    for i in TRUE_GENES:
        add(peptide(n), locus(i), False, 0.0001, False, "cds", locus(i))
        n += 1
    for i in TOLERANT_GENES:  # an SNV peptide: true, maps by the tolerant route
        add(peptide(n), locus(i), False, 0.0001, False, "cds", locus(i))
        n += 1
    for i in FALSE_DETECTED:  # lead decision 9: a protein null of an undetected locus
        add(peptide(n), locus(i), False, 0.0001, True, "undetected_product", locus(i))
        n += 1
    for _ in range(N_UNMAPPED):
        add(peptide(n), "GIAPORF_antisense_U00096.3_9001_9150_m", False, 0.002, True,
            "sixframe_noncoding", "NA")  # fmt: skip
        n += 1
    for i in DECOY_WINS:
        add(peptide(n), "rev_" + locus(i), True, 0.001, True, "", "")
        n += 1
    for i in DECOY_ABOVE_QMAX:
        add(peptide(n), "rev_" + locus(i), True, 0.05, True, "", "")
        n += 1
    tsv(inp / "peptides.tsv",
        ["peptide", "protein_ids", "q_value", "pep", "score", "is_decoy", "intensity", "run"],
        peps)  # fmt: skip
    tsv(inp / "truth" / "peptides.tsv",
        ["peptide", "is_wrong", "origin", "locus_tag", "contig", "start", "end"],
        truth)  # fmt: skip
    tsv(inp / "truth" / "genes.tsv", ["locus_tag", "tu_id", "contig", "start", "end", "strand"],
        [[locus(i), tu_id(i), CONTIG, *gene_coords(i)] for i in range(N_GENES)])  # fmt: skip
    tsv(
        inp / "truth" / "islands.tsv",
        [
            "island_id",
            "replicon",
            "first_tu",
            "n_tus",
            "effect_sd",
            "effect_log2",
            "first_tu_id",
            "last_tu_id",
        ],
        [
            [
                f"ISL_{k + 1}",
                CONTIG,
                t,
                n,
                e,
                e * 1.6,
                f"TU_{t}",
                f"TU_{t + n - 1}",
            ]
            for k, (t, n, e) in enumerate(
                zip(ISLAND_FIRST_TUS, ISLAND_N_TUS, ISLAND_EFFECT_SD, strict=True)
            )
        ],
    )
    tsv(inp / "truth" / "alt_sites.tsv",
        ["contig", "codon_start", "strand", "gene", "ref_aa", "alt_aa", "n_true_peptides",
         "frameshift_upstream"],
        [[CONTIG, alt_codon(k), "+", locus(10 + k), "A", "V", 2 if k < N_ALT_MULTI else 1,
          "false"] for k in range(N_ALT)])  # fmt: skip
    tsv(inp / "truth" / "novel_alt_sites.tsv",
        ["contig", "codon_start", "strand", "orf_id", "ref_aa", "alt_aa", "n_true_peptides"],
        [])  # fmt: skip
    rows = [_novel_truth(k, c) for c in NOVEL_CLASSES for k in range(N_NOVEL_TRUTH)]
    tsv(inp / "truth" / "novel_orfs.tsv", list(rows[0]), [list(r.values()) for r in rows])
    (inp / "provenance.json").write_text(
        json.dumps({"seed": 20261023, "tu_source": "regulondb", "decoys": "scored_null"}, indent=2)
        + "\n"
    )
    (inp / "config.yaml").write_text("seed: 20261023\nsample_name: k12_evidence\n")
    (inp / "abundance.tsv").write_text("locus_tag\tlog_abundance\n")


def alt_codon(k: int) -> int:
    """1-based lowest coordinate of truth ALT site ``k`` (inside gene 10 + k, + frame)."""
    return gene_coords(10 + k)[0] + 30


# ---------------------------------------------------------------- published GIAP outputs


def _genes_tsv(out: Path) -> None:
    rows = []
    for i in range(N_GENES):
        if i < N_PARTIAL:
            status, det, q = "partial_excluded", "false", None
        elif i in DETECTED_TRUE or i in FALSE_DETECTED:
            status, det, q = "detected", "true", 0.001
        elif i in TRUE_GENES:
            status, det, q = "not_detected", "false", 0.05
        else:
            status, det, q = "not_detected", "false", None
        rows.append([locus(i), "cds", det, "decoy_picked", status, q, "NA"])
    tsv(out / "tables" / "genes.tsv",
        ["locus_tag", "feature_type", "detected", "detection_tier", "detection_status",
         "gene_q_value", "flags"], rows)  # fmt: skip


def _detection_summary(out: Path) -> None:
    n_det = len(DETECTED_TRUE) + len(FALSE_DETECTED)
    tier = [("n_detected", n_det), ("n_detected_2pep", n_det), ("n_decoys_paired", 45),
            ("n_decoys_unpaired", 0), ("decoy_paired_share", 1.0), ("decoy_win_weight", 1.0),
            ("q_floor", Q_FLOOR)]  # fmt: skip
    rows = [["tier", "tier", "decoy_picked", m, v, None, None, None, ""] for m, v in tier]
    rows += [["tier", "q_threshold", "0.01", "n_detected", n_det, None, None, None, ""]]
    tsv(out / "tables" / "detection_summary.tsv",
        ["section", "group_type", "group", "metric", "value", "ci_low", "ci_high", "n", "flags"],
        rows)  # fmt: skip


def _mapping(out: Path, inp: Path) -> None:
    hits, unmapped = [], []
    truth = {r["peptide"]: r for r in read_tsv(inp / "truth" / "peptides.tsv")}
    tolerant = set()
    seen: dict[str, int] = {}
    for r in read_tsv(inp / "peptides.tsv"):
        if r["is_decoy"] == "true":
            continue
        gene = r["protein_ids"]
        seen[gene] = seen.get(gene, 0) + 1
        if seen[gene] == 2:
            tolerant.add(r["peptide"])
    for r in read_tsv(inp / "peptides.tsv"):
        if r["is_decoy"] == "true":
            continue
        pep, t = r["peptide"], truth[r["peptide"]]
        if t["origin"] == "sixframe_noncoding":
            unmapped.append([pep, pep, "peptide", "unmapped", "no_hit", "NA"])
            continue
        route = "tolerant" if pep in tolerant else "exact"
        hits.append([pep, route, 1 if route == "tolerant" else 0, "cds_inframe", "unique",
                     r["protein_ids"]])  # fmt: skip
    tsv(out / "mapping" / "peptide_hits.tsv",
        ["peptide", "route", "n_diff", "locus_class", "uniqueness", "gene"], hits)  # fmt: skip
    tsv(
        out / "mapping" / "unmapped.tsv",
        ["sequence_id", "sequence", "kind", "label", "reason", "contaminant_match"],
        unmapped,
    )


VCF_HEADER = (
    "##fileformat=VCFv4.3\n"
    f"##contig=<ID={CONTIG},length=4641652>\n"
    '##INFO=<ID=GENE,Number=1,Type=String,Description="gene">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _variants(out: Path) -> None:
    rows = []
    for k in range(N_ALT):
        rows.append([f"VAR_{k + 1}", CONTIG, alt_codon(k), "+", locus(10 + k), 11, "GCC", "A",
                     "V", "substitution", "ALT", "high", 0, 2, 1, "GTC",
                     "snv;strain_difference;database_derived"])  # fmt: skip
    for k in range(5):  # candidates: never counted as high
        rows.append([f"VAR_{N_ALT + k + 1}", CONTIG, alt_codon(k) + 300, "+", locus(10 + k), 111,
                     "GCC", "A", "T", "substitution", "ALT", "candidate", 0, 1, 1, "ACC",
                     "snv;strain_difference;database_derived"])  # fmt: skip
    tsv(out / "tables" / "variants.tsv",
        ["variant_id", "contig", "codon_start", "strand", "gene", "aa_pos", "ref_codon", "ref_aa",
         "alt_aa", "alt_type", "call", "confidence", "n_ref", "n_alt", "nt_changes",
         "alt_codon_vcf", "flags"], rows)  # fmt: skip
    recs = "".join(
        f"{CONTIG}\t{alt_codon(k) + 1}\tVAR_{k + 1}\tC\tT\t.\tPASS\tGENE={locus(10 + k)}\n"
        for k in range(N_ALT)
    )
    (out / "tables" / "variants.vcf").write_text(VCF_HEADER + recs)


def novel_call(t: Mapping[str, Any], *, called: bool = True, fdr: str = "0.005") -> list[Any]:
    """A ``novel_orfs.tsv`` row for truth ORF ``t``; GIAP's ``end`` is the last base of the stop
    in transcription direction (the lowest coordinate on ``-``)."""
    plus = t["strand"] == "+"
    start_longest = t["start"] if plus else t["end"]
    end = t["end"] if plus else t["start"]
    return [f"GIAPNOVEL_{t['orf_id']}", t["class"], t["contig"], t["strand"], start_longest, end,
            fdr, "true" if called else "false", "NA"]  # fmt: skip


NOVEL_HEADER = ["orf_id", "class", "contig", "strand", "start_longest", "end", "class_fdr",
                "called", "flags"]  # fmt: skip


def _novel(out: Path) -> None:
    rows = [novel_call(_novel_truth(k, c)) for c in NOVEL_CLASSES for k in range(N_NOVEL_CALLS)]
    tsv(out / "tables" / "novel_orfs.tsv", NOVEL_HEADER, rows)


def island_genes(k: int) -> list[str]:
    first, n = ISLAND_FIRST_TUS[k], ISLAND_N_TUS[k]
    return [locus(i) for i in range(N_GENES) if first <= i // 2 + 1 < first + n]


ISLAND_HEADER = ["island_id", "direction", "contig", "start", "end", "n_tus", "n_genes", "genes",
                 "p_value", "q_value", "flags"]  # fmt: skip


def island_call(k: int, genes: Sequence[str], direction: str | None = None) -> list[Any]:
    d = direction or ("high" if ISLAND_EFFECT_SD[k] > 0 else "low")
    return [f"island_{k + 1}", d, CONTIG, 1, 2, ISLAND_N_TUS[k], len(genes), ";".join(genes),
            0.001, "NA", "NA"]  # fmt: skip


def _context(out: Path) -> None:
    tsv(out / "tables" / "islands.tsv", ISLAND_HEADER,
        [island_call(k, island_genes(k)) for k in range(len(ISLAND_FIRST_TUS))])  # fmt: skip
    ac = []
    for null, val, p in (("gene_label", 0.41, 0.001), ("tu_block", 0.22, 0.001)):
        for lag in (1, 2):
            ac.append([CONTIG, "morans_i", lag, null, val / lag, p, "moran_i_permutation", 999,
                       "NA"])  # fmt: skip
    tsv(out / "tables" / "autocorrelation.tsv",
        ["contig", "statistic", "lag", "null_model", "value", "p_value", "method", "n_perm",
         "flags"], ac)  # fmt: skip
    prophage = [locus(i) for i in range(3000, 3010)]
    tsv(out / "annotation" / "regions.tsv",
        ["region_id", "region_type", "contig", "start", "end", "source_tool", "score", "genes",
         "notes"],
        [["crispr_1", "crispr", CONTIG, 62388, 62799, "Bakta", "NA", "NA", ""],
         ["prophage_1", "prophage", CONTIG, 3_000_001, 3_010_000, "geNOMAD", "0.98",
          ";".join(prophage), ""]])  # fmt: skip
    tsv(out / "tables" / "regions_stats.tsv",
        ["region_id", "region_type", "subset", "contig", "n_genes", "n_eligible", "stat", "effect",
         "p_value", "q_value", "null_model", "flags"],
        [["crispr_1", "crispr", "all", CONTIG, 0, 0, "o_minus_e", "NA", "NA", "NA",
          "circular_shift", "small_region;few_shifts"],
         ["prophage_1", "prophage", "all", CONTIG, 10, 10, "o_minus_e", 1.5, 0.02, 0.04,
          "circular_shift", "NA"],
         ["prophage_1", "prophage", "all", CONTIG, 10, 10, "hl_shift", 0.8, 0.03, 0.04,
          "circular_shift", "NA"]])  # fmt: skip


def _operons(out: Path) -> None:
    pairs = []
    for i in range(N_GENES - 1):
        same_tu = i // 2 == (i + 1) // 2
        cls = "OP" if same_tu else "NO"
        test = "both_detected_t" if same_tu and i < 2000 else "none"
        pairs.append([f"PAIR_{i + 1}", cls, locus(i), locus(i + 1), test, "NA"])
    tsv(out / "tables" / "pairs.tsv", ["pair_id", "class", "gene_a", "gene_b", "test", "flags"],
        pairs)  # fmt: skip
    tsv(out / "tables" / "operons.tsv", ["tu_id", "contig", "strand", "genes", "n_genes"],
        [[f"TU_{k + 1}", CONTIG, "+", f"{locus(2 * k)};{locus(2 * k + 1)}", 2]
         for k in range(N_TUS)])  # fmt: skip
    tsv(out / "tables" / "operon_flags.tsv",
        ["pair_id", "gene_a", "gene_b", "uniop_class", "flag", "q_value"],
        [["PAIR_1", locus(0), locus(1), "OP", "unexplained_discordance", 0.02]])  # fmt: skip


def section_html(sid: str, content: str) -> str:
    return f'<section id="{sid}">\n<h2>{sid}</h2>\n{content}\n</section>\n'


FIGURE = '<div class="plot" data-plot="V07">\n<figure>\n<svg xmlns="http://www.w3.org/2000/svg"></svg>\n</figure>\n</div>'
TABLE = '<table data-table="genes_tsv">\n<thead><tr><th>locus_tag</th></tr></thead>\n<tbody><tr><td>K12_00050</td></tr></tbody>\n</table>'  # noqa: E501
PLACEHOLDER = '<div class="plot" data-plot="V12">\n<p class="note">V12: plot file absent</p>\n</div>\n<p class="note">variants.tsv: not produced in this run.</p>'  # noqa: E501


def report_html(sections: Mapping[str, str]) -> str:
    body = "".join(section_html(s, sections.get(s, "")) for s in ALL_SECTIONS)
    return (
        f'<!DOCTYPE html>\n<html lang="en">\n<body>\n<h1>GIAP report</h1>\n{body}</body>\n</html>\n'
    )


def _report(out: Path) -> None:
    content = {s: FIGURE + "\n" + TABLE for s in REPORT_SECTIONS}
    (out / "report").mkdir(parents=True, exist_ok=True)
    (out / "report" / "report.html").write_text(report_html(content))


TRACE_HEADER = ["task_id", "hash", "native_id", "name", "tag", "status", "exit", "attempt",
                "submit", "duration", "realtime", "%cpu", "peak_rss", "peak_vmem", "rchar",
                "wchar", "workdir", "hostname"]  # fmt: skip


def trace_row(
    k: int, name: str, realtime: str, cpu: str, rss: str, status: str = "COMPLETED"
) -> list[Any]:
    return [
        k,
        f"aa/{k:06d}",
        1870000 + k,
        name,
        "-",
        status,
        0,
        1,
        "2026-10-24 00:00:00.000",
        realtime,
        realtime,
        cpu,
        rss,
        "20 GB",
        "1 MB",
        "1 MB",
        f"/w/{k}",
        "barley-01",
    ]


def _manifest_and_logs(out: Path) -> None:
    (out / "resolved_config.yaml").write_text(
        f"mapping:\n  tolerant:\n    enabled: true\n    min_len: {TOLERANT_MIN_LEN}\n"
    )
    manifest = {
        "exit_code": 0,
        "status": "completed",
        "started_at": "2026-10-24T00:00:00+00:00",
        "finished_at": "2026-10-24T05:00:00+00:00",
        "detection_tier": "decoy_picked",
        "steps": [
            {"name": n, "status": "completed"}
            for n in (
                "parse_bakta",
                "map_peptides",
                "gene_table",
                "context_stats",
                "operons",
                "variants",
                "novel",
                "plots",
            )
        ],
        "flags": ["bakta_db_light_absent_ok"],
        "warnings": ["context_layers:x.y"],
    }
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "logs").mkdir(parents=True, exist_ok=True)
    (out / "logs" / "events.jsonl").write_text(json.dumps(
        {"data": {"flag": "layer_missing:antismash"}, "event": "context_layers.layer_missing",
         "level": "WARNING", "step": "context_layers"}) + "\n")  # fmt: skip
    tsv(
        out / "logs" / "trace.tsv",
        TRACE_HEADER,
        [
            trace_row(1, "ANNOTATE:TERMINATORNET (k12)", "1h", "200.0%", "10 GB"),
            trace_row(2, "ANNOTATE:PROMOTER_CALCULATOR (TU_1)", "30m", "400.0%", "2 GB"),
            trace_row(3, "ANNOTATE:PROMOTER_CALCULATOR (TU_2)", "30m", "400.0%", "3 GB"),
            trace_row(4, "GENES:GENE_TABLE (k12)", "2h", "100.0%", "12 GB"),
            trace_row(5, "GENES:GENE_TABLE (k12)", "9h", "100.0%", "90 GB", status="CACHED"),
        ],
    )


def build(tmp: Path) -> Dirs:
    """The all-pass run; see the module docstring."""
    out, inp, again = tmp / "run", tmp / "inputs", tmp / "inputs_te7"
    _write_inputs(inp)
    shutil.copytree(inp, again)
    _genes_tsv(out)
    _detection_summary(out)
    _mapping(out, inp)
    _variants(out)
    _novel(out)
    _context(out)
    _operons(out)
    _report(out)
    _manifest_and_logs(out)
    return Dirs(out, inp, again)


def set_tier_metric(out: Path, metric: str, value: Any) -> None:
    """Set the ``tier``/``tier`` row ``metric`` of ``detection_summary.tsv``."""

    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            if r["group_type"] == "tier" and r["metric"] == metric:
                r["value"] = str(value)
        return rows

    edit_tsv(out / "tables" / "detection_summary.tsv", fn)
