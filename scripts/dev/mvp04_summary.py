"""MVP-04 evidence-run summary (docs/specs/MVP-04-ecoli-evidence-run.md "The run", acceptance
gates 1-16; lead decisions 11, 15, 19, 20 and 22).

Reads only published GIPPA outputs and the gippasim inputs (plus a second gippasim call with the same
seed for gate 9) and writes ``<outdir>/evidence_summary.json`` (strict JSON: NaN as null):
``{"gates": {"1": entry, ..., "16": entry}, "all_gates_ok", "outside_expected"}``. Each entry has
``value``, ``expected``, ``gate``, ``gate_ok``, ``in_expected`` and ``status`` (``pass``,
``outside_expected`` = gate met but value outside the expected range, reported; ``fail``) plus
per-gate keys; the gate definitions are pinned in tests/unit/test_mvp04_summary_gates.py. Each
``gate<N>`` only reads files and never raises on a missing output (the gate then fails).
Gate 13 (decision 20) gates the one planted high island; the three planted lows are exempt
(structural cap) and reported with their ceiling from ``gippa.context.islands``.

Gate 7 (part 1, ``gate7()``):

* a gene is **true** when it is the ``locus_tag`` of a truth row with ``is_wrong = false``
  (decision 11; a protein-entry null's truth row carries its locus tag but ``is_wrong = true``);
* ``realised_gene_fdp`` = detected genes that are not true / detected genes;
* ``false_detected_by_channel`` counts the false-detected genes with at least one
  ``peptide_hits`` row (``gene`` = the gene) of a null peptide of the channel; a gene may count
  in several (decision 15): ``protein_null`` (truth ``origin = undetected_product``),
  ``nterm_extension_sixframe_null`` (``origin = sixframe_noncoding`` whose target rows in the
  input carry a ``GIPPAORF_nterm_extension_*`` entry), ``tolerant_sixframe_null`` (any other
  six-frame null whose hit on the gene has ``route = tolerant``), and ``other`` (none of these).

    python scripts/dev/mvp04_summary.py OUTDIR INPUTS [--inputs-again DIR]
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import importlib.util
import json
import math
import re
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

from gippa.context.islands import low_direction_ceiling

PROTEIN_NULL_ORIGIN = "undetected_product"  # MVP-04 lead decision 9
SIXFRAME_NULL_ORIGIN = "sixframe_noncoding"
NTERM_PREFIX = "GIPPAORF_nterm_extension_"
CHANNELS = ("protein_null", "nterm_extension_sixframe_null", "tolerant_sixframe_null", "other")
SUMMARY_NAME = "evidence_summary.json"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _ids(cell: str | None) -> list[str]:
    return [x for x in (cell or "").split(";") if x and x != "NA"]


def _true(cell: str | None) -> bool:
    return (cell or "").strip().lower() == "true"


def gate7_from_rows(
    genes: Iterable[Mapping[str, str]],
    hits: Iterable[Mapping[str, str]],
    peptides: Iterable[Mapping[str, str]],
    truth: Iterable[Mapping[str, str]],
) -> dict[str, Any]:
    """Gate 7 from the four tables' rows (module doc)."""
    detected = {g["locus_tag"] for g in genes if _true(g.get("detected"))}
    true_genes: set[str] = set()
    origin: dict[str, str] = {}
    for t in truth:
        if _true(t.get("is_wrong")):
            if t.get("origin") in (PROTEIN_NULL_ORIGIN, SIXFRAME_NULL_ORIGIN):
                origin[t["peptide"]] = str(t["origin"])
        else:
            true_genes.update(_ids(t.get("locus_tag")))
    nterm = {
        p["peptide"]
        for p in peptides
        if not _true(p.get("is_decoy"))
        and any(e.startswith(NTERM_PREFIX) for e in _ids(p.get("protein_ids")))
    }
    false_genes = detected - true_genes
    by: dict[str, set[str]] = {k: set() for k in CHANNELS}
    for h in hits:
        gene, pep = h.get("gene", ""), h["peptide"]
        if gene not in false_genes or pep not in origin:
            continue
        if origin[pep] == PROTEIN_NULL_ORIGIN:
            by["protein_null"].add(gene)
        elif pep in nterm:
            by["nterm_extension_sixframe_null"].add(gene)
        elif h.get("route") == "tolerant":
            by["tolerant_sixframe_null"].add(gene)
    by["other"] = false_genes - set().union(*(by[k] for k in CHANNELS[:-1]))
    n_det, n_false = len(detected), len(false_genes)
    return {
        "n_detected": n_det,
        "n_false_detected": n_false,
        "realised_gene_fdp": n_false / n_det if n_det else float("nan"),
        "false_detected_by_channel": {k: len(by[k]) for k in CHANNELS},
    }


def gate7(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Gate 7 of the run in ``outdir`` against the gippasim inputs in ``inputs``."""
    outdir, inputs = Path(outdir), Path(inputs)
    return gate7_from_rows(
        _rows(outdir / "tables" / "genes.tsv"),
        _rows(outdir / "mapping" / "peptide_hits.tsv"),
        _rows(inputs / "peptides.tsv"),
        _rows(inputs / "truth" / "peptides.tsv"),
    )


# ---------------------------------------------------------------- part 2: gates 1-16

ALPHA = 0.01  # stats.alpha.identification
Q_VALUE_MAX = 0.01  # input.q_value_max
FORBIDDEN_WARNINGS = ("decoy_unscored", "decoy_pairing_low", "q_floor_above_alpha",
                      "no_eligible_genes")  # fmt: skip
REPORT_SECTIONS = ("detection", "genome-context", "operons", "variants", "novel-orfs")
RUNTIME_KEY = "mvp04_k12_evidence"
NOVEL_CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
TRUTH_TABLES = ("peptides", "genes", "islands", "alt_sites", "novel_alt_sites", "novel_orfs")
BUDGET = {"wall_h": 8.0, "cpu_h": 150.0, "rss_gb": 64.0}  # spec "Runtime budget"
Range = list[float | None]


def status_of(gate_ok: bool, in_expected: bool) -> str:
    """``pass``, ``outside_expected`` (gate met, value outside its range: reported) or ``fail``."""
    if not gate_ok:
        return "fail"
    return "pass" if in_expected else "outside_expected"


def _entry(
    value: Any, expected: Any, gate: str, ok: bool, inside: bool, **extra: Any
) -> dict[str, Any]:
    ok, inside = bool(ok), bool(inside)
    return {"value": value, "expected": expected, "gate": gate, "gate_ok": ok,
            "in_expected": inside, "status": status_of(ok, inside), **extra}  # fmt: skip


def _within(v: float | None, lo: float | None, hi: float | None) -> bool:
    return v is not None and (lo is None or v >= lo) and (hi is None or v <= hi)


def _read(path: Path) -> list[dict[str, str]]:
    """Rows of a TSV, ``[]`` when it is absent (the gate then fails)."""
    try:
        return _rows(path)
    except OSError:
        return []


def _json(path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _num(cell: str | None) -> float | None:
    try:
        v = float(cell)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _int(cell: str | None) -> int | None:
    v = _num(cell)
    return int(v) if v is not None else None


def _eligible(outdir: Path) -> list[dict[str, str]]:
    genes = _read(Path(outdir) / "tables" / "genes.tsv")
    return [g for g in genes if g.get("detection_status") != "partial_excluded"]


def _true_genes(inputs: Path) -> set[str]:
    out: set[str] = set()
    for t in _read(Path(inputs) / "truth" / "peptides.tsv"):
        if not _true(t.get("is_wrong")):
            out.update(_ids(t.get("locus_tag")))
    return out


def _tier(outdir: Path) -> str | None:
    tiers = {g.get("detection_tier", "") for g in _read(Path(outdir) / "tables" / "genes.tsv")}
    tiers.discard("")
    if not tiers:
        return None
    return tiers.pop() if len(tiers) == 1 else "mixed"


def _tier_metric(outdir: Path, metric: str) -> float | None:
    rows = [r for r in _read(Path(outdir) / "tables" / "detection_summary.tsv")
            if r.get("group_type") == "tier" and r.get("metric") == metric]  # fmt: skip
    tier = _tier(outdir)
    rows.sort(key=lambda r: r.get("group") != tier)  # the run's tier first
    return _num(rows[0].get("value")) if rows else None


def gate1(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Run completes: exit 0, ``status = completed``, no failed step."""
    m = _json(Path(outdir) / "run_manifest.json") or {}
    failed = [str(s.get("name")) for s in m.get("steps", []) if s.get("status") == "failed"]
    code, status = m.get("exit_code"), m.get("status")
    ok = code == 0 and status == "completed" and not failed
    return _entry(code, 0, "exit 0; no step failed", ok, ok, run_status=status,
                  failed_steps=failed)  # fmt: skip


def gate2(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Detection tier ``decoy_picked`` and none of ``FORBIDDEN_WARNINGS``."""
    outdir = Path(outdir)
    m = _json(outdir / "run_manifest.json") or {}
    texts = [str(x) for k in ("flags", "warnings") for x in m.get(k, []) or []]
    with contextlib.suppress(OSError):
        texts += (outdir / "logs" / "events.jsonl").read_text().splitlines()
    # exact tokens only: novel_orfs.class_q_floor_above_alpha is a different (per-class) warning
    pats = {t: re.compile(r"(?<![A-Za-z0-9_])" + re.escape(t) + r"(?![A-Za-z0-9_])")
            for t in FORBIDDEN_WARNINGS}  # fmt: skip
    found = {t for t, pat in pats.items() if any(pat.search(x) for x in texts)}
    q_floor = _tier_metric(outdir, "q_floor")
    if q_floor is not None and q_floor > ALPHA:
        found.add("q_floor_above_alpha")
    tier = _tier(outdir)
    ok = tier == "decoy_picked" and not found
    return _entry(tier, "decoy_picked", "decoy_picked; no " + ", ".join(FORBIDDEN_WARNINGS), ok,
                  tier == "decoy_picked", warnings_found=sorted(found))  # fmt: skip


def gate3(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Paired share f >= 0.8."""
    f = _tier_metric(outdir, "decoy_paired_share")
    return _entry(f, [0.98, 1.0], "f >= 0.8", f is not None and f >= 0.8, _within(f, 0.98, 1.0))


def gate4(outdir: Path, inputs: Path) -> dict[str, Any]:
    """T = 1 / q_floor >= 110."""
    q = _tier_metric(outdir, "q_floor")
    t = round(1 / q) if q else None
    ok = q is not None and q > 0 and q <= 1 / 110 + 1e-15
    return _entry(t, [2300, 3100], "T >= 110 (q_floor <= 1/110)", ok, _within(t, 2300, 3100),
                  q_floor=q)  # fmt: skip


def gate5(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Decoy wins: eligible genes with gene q NA and a paired decoy at q <= Q_VALUE_MAX."""
    paired: set[str] = set()
    for p in _read(Path(inputs) / "peptides.tsv"):
        q = _num(p.get("q_value"))
        if _true(p.get("is_decoy")) and q is not None and q <= Q_VALUE_MAX:
            paired.update(e[4:] for e in _ids(p.get("protein_ids")) if e.startswith("rev_"))
    genes = [g for g in _eligible(outdir) if g["locus_tag"] in paired]
    wins = sum(1 for g in genes if _num(g.get("gene_q_value")) is None)
    return _entry(wins, [20, 300], ">= 1 decoy win", wins >= 1, _within(wins, 20, 300),
                  n_paired_genes=len(genes))  # fmt: skip


def gate6(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Detected at alpha: >= 1000 and <= 0.90 x eligible."""
    n = _tier_metric(outdir, "n_detected")
    elig = _eligible(outdir)
    true = _true_genes(inputs)
    n_true = sum(1 for g in elig if g["locus_tag"] in true)
    value = int(n) if n is not None else None
    ok = n is not None and n >= 1000 and n <= 0.90 * len(elig)
    lo, hi = 0.85 * n_true, float(n_true)
    return _entry(value, [lo, hi], ">= 1000 and <= 0.90 x eligible", ok, _within(n, lo, hi),
                  n_eligible=len(elig), n_truth_detected=n_true)  # fmt: skip


def gate7_entry(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Gate 7: ``gate7()`` plus the entry keys; FDP <= 0.02."""
    try:
        res = gate7(outdir, inputs)
    except (OSError, KeyError):
        res = {}
    fdp = res.get("realised_gene_fdp")
    fdp = fdp if isinstance(fdp, float) and math.isfinite(fdp) else None
    return {**res, **_entry(fdp, [0.0, 0.01], "realised gene FDP <= 0.02",
                            fdp is not None and fdp <= 0.02, _within(fdp, 0.0, 0.01))}  # fmt: skip


def gate8(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Truth-detected share of eligible genes in [0.55, 0.70]."""
    elig = _eligible(outdir)
    true = _true_genes(inputs)
    n_true = sum(1 for g in elig if g["locus_tag"] in true)
    share = n_true / len(elig) if elig else None
    inside = _within(share, 0.55, 0.70)
    return _entry(share, [0.55, 0.70], "share in [0.55, 0.70]", inside, inside,
                  n_truth_detected=n_true, n_eligible=len(elig))  # fmt: skip


def _compared(d: Path) -> set[str]:
    rel = {f"truth/{p.name}" for p in (d / "truth").glob("*.tsv")}
    return rel | ({"peptides.tsv"} if (d / "peptides.tsv").is_file() else set())


def gate9(inputs: Path, again: Path | None) -> dict[str, Any]:
    """TE7 on the inputs: ``peptides.tsv`` and ``truth/*.tsv`` byte-identical across two calls."""
    text = "byte-identical peptides.tsv and truth/*.tsv"
    required = {"peptides.tsv", *(f"truth/{t}.tsv" for t in TRUTH_TABLES)}
    if again is None:
        return _entry(None, [], text, False, False)
    a, b = Path(inputs), Path(again)
    fa, fb = _compared(a), _compared(b)
    diff = sorted((fa ^ fb) | {r for r in fa & fb if (a / r).read_bytes() != (b / r).read_bytes()})
    ok = not diff and required <= fa and required <= fb
    return _entry(diff, [], text, ok, not diff)


def gate10(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Mapping shares over distinct non-decoy peptides; tolerant share > 0."""
    outdir = Path(outdir)
    target, decoy = set(), set()
    for p in _read(Path(inputs) / "peptides.tsv"):
        (decoy if _true(p.get("is_decoy")) else target).add(p["peptide"])
    routes: dict[str, set[str]] = {}
    for h in _read(outdir / "mapping" / "peptide_hits.tsv"):
        routes.setdefault(h["peptide"], set()).add(h.get("route", ""))
    for u in _read(outdir / "mapping" / "unmapped.tsv"):
        if u.get("kind") == "peptide":
            routes.setdefault(u.get("sequence", ""), set())
    peps = [p for p in routes if p not in decoy - target]
    n = len(peps)
    exact = sum(1 for p in peps if "exact" in routes[p])
    tol = sum(1 for p in peps if "tolerant" in routes[p] and "exact" not in routes[p])
    counts = {"exact": exact, "tolerant": tol, "unmapped": n - exact - tol}
    shares = {k: (v / n if n else None) for k, v in counts.items()}
    ok = bool(n) and tol > 0
    return _entry(shares["exact"], [0.95, 1.0], "tolerant share > 0", ok,
                  _within(shares["exact"], 0.95, 1.0), shares=shares, n_peptides=n)  # fmt: skip


def _site(contig: str | None, codon: str | None, aa: str | None) -> tuple[str, int | None, str]:
    return (contig or "", _int(codon), aa or "")


def _vcf_valid(path: Path) -> bool:
    try:
        import pysam

        with pysam.VariantFile(str(path)) as vf:
            for _ in vf:
                pass
    except Exception:  # any parse error means invalid
        return False
    return True


_IL = str.maketrans({"I": "J", "L": "J"})
NOVEL_STARTS = frozenset({"GTG", "TTG", "CTG"})  # non-ATG starts reading M as initiator


def _fold(s: str | None) -> str:
    """I, L and J as one token (GIPPA's folded I/L ALT residue is ``J``, M7-02 rule 5)."""
    return (s or "").translate(_IL)


def _tolerant_min_len(outdir: Path) -> int | None:
    """``mapping.tolerant.min_len`` of ``OUTDIR/resolved_config.yaml``; ``None`` when unreadable."""
    import yaml

    try:
        cfg = yaml.safe_load((Path(outdir) / "resolved_config.yaml").read_text())
        v = cfg["mapping"]["tolerant"]["min_len"]
    except (OSError, yaml.YAMLError, KeyError, TypeError):
        return None
    return int(v) if isinstance(v, int | float) and not isinstance(v, bool) else None


def _novel_starts(inputs: Path) -> set[tuple[str, str, int | None]]:
    """(contig, strand, 1-based lowest base) of every truth novel ORF's start codon."""
    out = set()
    for o in _read(Path(inputs) / "truth" / "novel_orfs.tsv"):
        ts = _int(o.get("true_start"))
        low = ts if o.get("strand") == "+" or ts is None else ts - 2
        out.add((o.get("contig") or "", o.get("strand") or "", low))
    return out


def _initiator_row(r: Mapping[str, str], starts: set[tuple[str, str, int | None]]) -> bool:
    """An ALT M over an initiator: a ``start_codon`` flag, or a GTG/TTG/CTG truth novel start."""
    if "start_codon" in _ids(r.get("flags")):
        return True
    return (r.get("alt_aa") == "M" and (r.get("ref_codon") or "").upper() in NOVEL_STARTS
            and (r.get("contig") or "", r.get("strand") or "", _int(r.get("codon_start")))
            in starts)  # fmt: skip


def _recall_sites(
    truth: Sequence[Mapping[str, str]], peptides: Sequence[Mapping[str, str]], min_len: int
) -> set[tuple[str, int | None, str]]:
    """Folded truth sites with >= 2 distinct I/L-folded correct peptides of length >= ``min_len``
    of the site's locus (``gene`` or ``orf_id``) whose 1-based closed interval holds the codon."""
    by: dict[tuple[str, str], list[tuple[int, int, str]]] = {}
    for p in peptides:
        if p.get("origin") == "decoy" or _true(p.get("is_wrong")):
            continue
        seq, a, b = p.get("peptide") or "", _int(p.get("start")), _int(p.get("end"))
        if len(seq) < min_len or a is None or b is None:
            continue
        for locus in _ids(p.get("locus_tag")):
            by.setdefault((locus, p.get("contig") or ""), []).append((a, b, _fold(seq)))
    out = set()
    for t in truth:
        cs, contig = _int(t.get("codon_start")), t.get("contig") or ""
        if cs is None:
            continue
        locus = t.get("gene") or t.get("orf_id") or ""
        seen = {s for a, b, s in by.get((locus, contig), ()) if a <= cs and cs + 2 <= b}
        if len(seen) >= 2:
            out.add(_site(contig, t.get("codon_start"), _fold(t.get("alt_aa"))))
    return out


def gate11(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Variants: ``high`` ALT rows against ``truth/alt_sites.tsv`` + ``novel_alt_sites.tsv``.

    Residues compare with I/L folded (ALT ``J`` = truth ``I`` or ``L``); initiator rows (ALT M
    at a ``start_codon``-flagged or GTG/TTG/CTG novel-ORF start: gippasim leaves initiators out
    of truth) leave the precision denominator (``n_initiator_m_excluded``); recall sites have
    >= 2 true peptides of ``truth/peptides.tsv`` of at least ``mapping.tolerant.min_len``
    (``resolved_config.yaml``; no config, no recall). ``value`` is every high row."""
    outdir, inputs = Path(outdir), Path(inputs)
    rows = _read(outdir / "tables" / "variants.tsv")
    truth = [*_read(inputs / "truth" / "alt_sites.tsv"),
             *_read(inputs / "truth" / "novel_alt_sites.tsv")]  # fmt: skip
    sites = {_site(t.get("contig"), t.get("codon_start"), _fold(t.get("alt_aa"))) for t in truth}
    starts = _novel_starts(inputs)
    high_rows = [r for r in rows if r.get("confidence") == "high"]
    init = [r for r in high_rows if _initiator_row(r, starts)]
    called = {_site(r.get("contig"), r.get("codon_start"), _fold(r.get("alt_aa")))
              for r in high_rows}  # fmt: skip
    scored = [_site(r.get("contig"), r.get("codon_start"), _fold(r.get("alt_aa")))
              for r in high_rows if not _initiator_row(r, starts)]  # fmt: skip
    n_high, n_true = len(high_rows), sum(1 for h in scored if h in sites)
    precision = n_true / len(scored) if scored else None
    min_len = _tolerant_min_len(outdir)
    multi = (_recall_sites(truth, _read(inputs / "truth" / "peptides.tsv"), min_len)
             if min_len is not None else set())  # fmt: skip
    recall = len(multi & called) / len(multi) if multi else None
    need = {"strain_difference", "database_derived"}
    missing = sum(1 for r in rows if not need <= set(_ids(r.get("flags"))))
    vcf = _vcf_valid(outdir / "tables" / "variants.vcf")
    ok = (n_high >= 30 and precision is not None and precision >= 0.95 and recall is not None
          and recall >= 0.90 and missing == 0 and vcf)  # fmt: skip
    inside = _within(n_high, 60, 300) and precision is not None and precision >= 0.97
    return _entry(n_high, [60, 300], "n_high >= 30; precision >= 0.95; recall >= 0.90; flags; VCF",
                  ok, inside, precision=precision, recall=recall, n_recall_sites=len(multi),
                  n_initiator_m_excluded=len(init), recall_min_len=min_len,
                  n_rows_missing_flags=missing, vcf_valid=vcf)  # fmt: skip


def gate12(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Novel ORF calls at class q <= ALPHA against the truth ORFs (contig, strand, stop)."""
    truth = _read(Path(inputs) / "truth" / "novel_orfs.tsv")
    keys = {(t.get("contig"), t.get("strand"), _int(t.get("stop"))) for t in truth}
    classes: dict[str, dict[str, Any]] = {}
    calls = _read(Path(outdir) / "tables" / "novel_orfs.tsv")
    for c in NOVEL_CLASSES:
        tc = [t for t in truth if t.get("class") == c]
        with_pep = sum(1 for t in tc if (_num(t.get("n_true_peptides")) or 0) >= 1)
        n_calls = n_true = 0
        for r in calls:
            fdr = _num(r.get("class_fdr"))
            if r.get("class") != c or not _true(r.get("called")) or fdr is None or fdr > ALPHA:
                continue
            end = _int(r.get("end"))
            stop = None if end is None else end - 2 if r.get("strand") == "+" else end
            n_calls += 1
            n_true += (r.get("contig"), r.get("strand"), stop) in keys
        classes[c] = {"n_truth": len(tc), "n_truth_with_peptide": with_pep, "n_calls": n_calls,
                      "n_true_calls": n_true, "precision": n_true / n_calls if n_calls else None,
                      "below_1_over_alpha": with_pep < 1 / ALPHA}  # fmt: skip
    live = [k for k in classes.values() if k["n_calls"] >= 1]
    ok = len(live) >= 2 and all(k["precision"] >= 0.95 for k in live)
    inside = all(k["n_calls"] >= 30 for k in classes.values() if not k["below_1_over_alpha"])
    return _entry(len(live), ">= 30 calls per class with >= 100 truth ORFs with a peptide",
                  ">= 2 non-vacuous classes; precision >= 0.95 each", ok, inside,
                  classes=classes)  # fmt: skip


def _tu_num(tu: str | None) -> int | None:
    return _int((tu or "").rsplit("_", 1)[-1])


def gate13(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Context: lag-1 tu_block Moran's I, the planted high island recovered, region stats present.

    Lead decisions 20 and 22: every planted high island must be recovered by a ``high`` call
    covering >= 50 % of its genes; a spurious call shares no gene with any planted island (lows
    included). The lows are exempt (structural cap) and reported in ``low_ceilings`` with w (the
    island's TUs holding >= 1 eligible ``tables/genes.tsv`` gene), q (the undetected share of the
    eligible genes) and the raw ``low_direction_ceiling(q, w)``.
    """
    outdir, inputs = Path(outdir), Path(inputs)
    tu_genes: dict[int, set[str]] = {}
    for g in _read(inputs / "truth" / "genes.tsv"):
        n = _tu_num(g.get("tu_id"))
        if n is not None:
            tu_genes.setdefault(n, set()).add(g["locus_tag"])
    planted: list[tuple[str, str, dict[int, set[str]]]] = []
    for isl in _read(inputs / "truth" / "islands.tsv"):
        lo, hi = _tu_num(isl.get("first_tu_id")), _tu_num(isl.get("last_tu_id"))
        tus = {k: v for k, v in tu_genes.items() if lo is not None and hi is not None
               and lo <= k <= hi}  # fmt: skip
        d = "high" if (_num(isl.get("effect_sd")) or 0) > 0 else "low"
        planted.append((isl.get("island_id", ""), d, tus))
    calls = [(c.get("direction"), set(_ids(c.get("genes"))))
             for c in _read(outdir / "tables" / "islands.tsv")]  # fmt: skip
    highs = [set().union(*tus.values()) for _, d, tus in planted if d == "high"]
    n_rec = sum(1 for gs in highs if gs and any(
        cd == "high" and len(cg & gs) >= 0.5 * len(gs) for cd, cg in calls))  # fmt: skip
    all_planted = set().union(*(v for _, _, tus in planted for v in tus.values()))
    n_spur = sum(1 for _, cg in calls if not cg & all_planted)
    eligible = _eligible(outdir)
    scored = {g["locus_tag"] for g in eligible}
    n_und = sum(1 for g in eligible if not _true(g.get("detected")))
    q = n_und / len(eligible) if eligible else math.nan
    low_ceilings: dict[str, dict[str, Any]] = {}
    for iid, d, tus in planted:
        if d == "low":
            w = sum(1 for v in tus.values() if v & scored)
            low_ceilings[iid] = {"w": w, "q": q, "ceiling": low_direction_ceiling(q, w)}
    ac = next((r for r in _read(outdir / "tables" / "autocorrelation.tsv")
               if r.get("statistic") == "morans_i" and _int(r.get("lag")) == 1
               and r.get("null_model") == "tu_block"), {})  # fmt: skip
    mi, p = _num(ac.get("value")), _num(ac.get("p_value"))
    stats = _read(outdir / "tables" / "regions_stats.tsv")
    types = {r.get("region_type", "") for r in _read(outdir / "annotation" / "regions.tsv")}
    by_type = {t: [s for s in stats if s.get("region_type") == t] for t in types}
    na_only = [t for t, ss in by_type.items()
               if any((_num(s.get("n_eligible")) or 0) >= 3 for s in ss)
               and not any(_num(s.get("effect")) is not None for s in ss)]  # fmt: skip
    missing = sorted(t for t, ss in by_type.items() if not ss)
    regions_ok = bool(types) and not missing and not na_only
    ok = (mi is not None and mi > 0 and p is not None and p <= 0.01 and bool(highs)
          and n_rec == len(highs) and n_spur <= 2 and regions_ok)  # fmt: skip
    return _entry(n_rec, [1, 1], "Moran's I > 0, p <= 0.01; every planted high island; <= 2 "
                  "spurious; regions_stats rows; lows exempt", ok,
                  _within(n_rec, 1, 1) and n_spur == 0, n_recovered_high=n_rec,
                  n_spurious=n_spur, low_exempt=True, low_ceilings=low_ceilings,
                  autocorrelation={"morans_i": mi, "p_value": p}, regions_ok=regions_ok,
                  regions_missing_types=missing)  # fmt: skip


def gate14(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Operons: tests run, no ``too_few_op_pairs``, >= 1000 TUs; RegulonDB recovery reported."""
    outdir, inputs = Path(outdir), Path(inputs)
    pairs = _read(outdir / "tables" / "pairs.tsv")
    operons = _read(outdir / "tables" / "operons.tsv")
    n_op = sum(1 for r in pairs if r.get("class") == "OP")
    tests_run = any(r.get("test", "none") != "none" for r in pairs)
    too_few = any("too_few_op_pairs" in _ids(r.get("flags")) for r in pairs)
    recovery: float | str | None = "not available"
    if (_json(inputs / "provenance.json") or {}).get("tu_source") == "regulondb":
        tus: dict[str, set[str]] = {}
        for g in _read(inputs / "truth" / "genes.tsv"):
            tus.setdefault(g.get("tu_id", ""), set()).add(g["locus_tag"])
        multi = [frozenset(v) for v in tus.values() if len(v) >= 2]
        called = {frozenset(_ids(o.get("genes"))) for o in operons}
        recovery = sum(1 for t in multi if t in called) / len(multi) if multi else None
    ok = tests_run and not too_few and len(operons) >= 1000
    return _entry(n_op, [1500, None], "tests not all none; no too_few_op_pairs; >= 1000 TUs", ok,
                  _within(n_op, 1500, None), n_pairs=len(pairs), n_tus=len(operons),
                  tests_run=tests_run, too_few_op_pairs=too_few,
                  regulondb_recovery=recovery)  # fmt: skip


def gate15(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Report: the five module sections hold a figure or a table; no q-floor caveat."""
    try:
        page = (Path(outdir) / "report" / "report.html").read_text()
    except OSError:
        page = ""
    body = dict(re.findall(r'<section id="([^"]+)"[^>]*>(.*?)</section>', page, re.S))
    full = [s for s in REPORT_SECTIONS
            if "<svg" in body.get(s, "") or "<table data-table=" in body.get(s, "")]  # fmt: skip
    caveat = "q_floor_above_alpha" in page
    ok = bool(page) and len(full) == len(REPORT_SECTIONS) and not caveat
    return _entry(len(full), [5, 5], "five module cards with content; no q-floor caveat", ok,
                  len(full) == 5, empty_sections=[s for s in REPORT_SECTIONS if s not in full],
                  q_floor_caveat=caveat)  # fmt: skip


def _measure_runtimes() -> Any:
    path = Path(__file__).resolve().parents[1] / "measure_runtimes.py"
    spec = importlib.util.spec_from_file_location("_mvp04_measure_runtimes", path)
    assert spec is not None and spec.loader is not None, path
    mod = sys.modules.get(spec.name) or importlib.util.module_from_spec(spec)
    if spec.name not in sys.modules:  # dataclasses need the module registered
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
    return mod


def gate16(outdir: Path, inputs: Path) -> dict[str, Any]:
    """Runtime (recorded, never failed): wall h, CPU-h, peak RSS against the budget."""
    outdir = Path(outdir)
    mr = _measure_runtimes()
    trace = _read(outdir / "logs" / "trace.tsv")
    tasks = mr.parse_trace(trace)
    cpu_h = 0.0
    for r in trace:
        if r.get("status", "").strip().upper() in mr._SKIP:
            continue
        secs = mr.seconds(r.get("realtime", "")) or mr.seconds(r.get("duration", "")) or 0.0
        cpu_h += secs / 3600 * (_num((r.get("%cpu") or "").rstrip("%")) or 0.0) / 100
    rss = max((t.rss_bytes for t in tasks), default=0.0) / 1024**3
    m = _json(outdir / "run_manifest.json") or {}
    try:
        span = dt.datetime.fromisoformat(m["finished_at"]) - dt.datetime.fromisoformat(
            m["started_at"]
        )
        wall: float | None = span.total_seconds() / 3600
    except (KeyError, TypeError, ValueError):
        wall = None
    over = ((wall or 0.0) > BUDGET["wall_h"] or cpu_h > BUDGET["cpu_h"]
            or rss > BUDGET["rss_gb"])  # fmt: skip
    rows = mr.summarise(tasks, RUNTIME_KEY, "?") if tasks else []
    return _entry(wall, [3, 8], "recorded; overrun -> backlog, never a failure", True,
                  _within(wall, 3, 8), cpu_h=cpu_h, peak_rss_gb=rss, over_budget=over,
                  rows=rows)  # fmt: skip


def summarise(outdir: Path, inputs: Path, inputs_again: Path | None = None) -> dict[str, Any]:
    """Every gate's entry, ``all_gates_ok`` and the gates outside their expected range."""
    gates: dict[str, dict[str, Any]] = {}
    for n in range(1, 17):
        if n == 7:
            gates["7"] = gate7_entry(outdir, inputs)
        elif n == 9:
            gates["9"] = gate9(inputs, inputs_again)
        else:
            gates[str(n)] = globals()[f"gate{n}"](outdir, inputs)
    return {"gates": gates, "all_gates_ok": all(e["gate_ok"] for e in gates.values()),
            "outside_expected": [n for n, e in gates.items()
                                 if e["status"] == "outside_expected"]}  # fmt: skip


def _strict(x: Any) -> Any:
    """NaN and infinities become null (strict JSON)."""
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, dict):
        return {k: _strict(v) for k, v in x.items()}
    if isinstance(x, list | tuple):
        return [_strict(v) for v in x]
    return x


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("outdir", type=Path)
    ap.add_argument("inputs", type=Path)
    ap.add_argument("--inputs-again", type=Path, default=None,
                    help="the second gippasim call with the same seed (gate 9)")  # fmt: skip
    a = ap.parse_args(argv)
    res = _strict(summarise(a.outdir, a.inputs, a.inputs_again))
    path = a.outdir / SUMMARY_NAME
    path.write_text(json.dumps(res, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
