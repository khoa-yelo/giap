"""TE7: a second A2 run, byte for byte (GIAP layer), and the DIA-NN replicate (an upstream
finding, never the verdict); GIAP_SPEC.md §11.7; spec docs/specs/M10-06-te-real-and-te7.md.

Pure functions plus ``run``, which ``te_real.py te7`` (the CLI) calls::

    python validation/pilot/te_real.py te7 --ref OUT_A --rep OUT_B --repo GIAP_REPO \\
        [--diann-a DIANN_A --diann-b DIANN_B --peptides-a PEP_A2 --peptides-b CONV_B] \\
        --out METRICS_DIR

``te7_files.tsv`` holds row-level data (paths and hashes) and is never committed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
import re
import subprocess
from collections import Counter
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from itertools import zip_longest
from pathlib import Path, PurePosixPath
from typing import Any, TypeAlias

import numpy as np
import pandas as pd
from validation.pilot import metrics as pm
from validation.pilot.te_real import _out_dir, _s, _spec_sha, _true, _write

CODE_PATHS = ("src/", "workflow/", "pyproject.toml", "uv.lock", "scripts/slurm/giap_run.sbatch")
EXCLUDED_TOP = ("work", "logs", "report", "preflight")
L3_KEY = ["run", "modified_peptide", "charge", "is_decoy"]
#: Nextflow process aliases whose registry step is not the lowercased alias minus ``_step``
STEP_ALIASES = {"report_html": "report"}
FILES_COLS = ["path", "class", "sha256_a", "sha256_b", "equal"]
Reader = Callable[[], bytes]
Frame: TypeAlias = pd.DataFrame | Callable[[], pd.DataFrame] | None


# ---------------------------------------------------------------- GIAP layer


def path_class(path: str) -> str:
    top = path.split("/", 1)[0]
    if top in EXCLUDED_TOP or top.startswith(".nextflow") or path == "run_manifest.json":
        return "excluded"
    if path.startswith("annotation/tools/"):
        return "external"
    if path.startswith("plots/") and not path.startswith("plots/data/"):
        return "figures"
    return "gate"


def _numeric(la: list[str], lb: list[str]) -> dict[str, Any] | None:
    """Same header and row count, and every differing cell parses as a number on both sides."""
    if len(la) != len(lb) or not la or la[0] != lb[0]:
        return None
    header = la[0].split("\t")
    cols: set[str] = set()
    n, max_abs, max_rel = 0, 0.0, 0.0
    for x, y in zip(la[1:], lb[1:], strict=True):
        cx, cy = x.split("\t"), y.split("\t")
        if x == y:
            continue
        if not len(cx) == len(cy) == len(header):
            return None
        for h, u, v in zip(header, cx, cy, strict=True):
            if u == v:
                continue
            try:
                fu, fv = float(u), float(v)
            except ValueError:
                return None
            n, d, scale = n + 1, abs(fu - fv), max(abs(fu), abs(fv))
            cols.add(h)
            if math.isfinite(d):
                max_abs, max_rel = max(max_abs, d), max(max_rel, d / scale if scale else 0.0)
    order = [h for h in dict.fromkeys(header) if h in cols]
    return {"columns": order, "n_cells": n, "max_abs": max_abs, "max_rel": max_rel} if n else None


def _diagnose(
    path: str,
    ba: bytes | None,
    bb: bytes | None,
    labels: Sequence[str],
    known: Mapping[str, Collection[str]],
    producer: Mapping[str, Any],
) -> dict[str, Any]:
    d: dict[str, Any] = dict.fromkeys(
        ("sub_class", "columns", "n_cells", "max_abs", "max_rel", "first_differing_line")
    )
    d |= {"path": path, "producer": producer.get("producer", "unknown"),
          "cpu_classes": producer.get("cpu_classes")}  # fmt: skip
    if ba is None or bb is None:
        return d | {"class": f"only_in:{labels[0] if bb is None else labels[1]}"}
    if path.endswith(".gz"):
        try:
            ba, bb = gzip.decompress(ba), gzip.decompress(bb)
        except (OSError, EOFError):
            return d | {"class": "content"}
        if ba == bb:
            return d | {"class": "gzip_header_only"}
    try:
        la, lb = ba.decode("utf-8").splitlines(), bb.decode("utf-8").splitlines()
    except UnicodeDecodeError:
        return d | {"class": "content"}
    if la != lb and la[:1] == lb[:1] and sorted(la[1:]) == sorted(lb[1:]):
        return d | {"class": "row_order_only"}
    num = _numeric(la, lb)
    if num is not None:
        name = PurePosixPath(path).name
        listed = set(known.get(name, known.get(name.removesuffix(".gz"), ())))
        cc = producer.get("cpu_classes") or {}
        ca, cb = (set(cc.get(x) or ()) for x in labels[:2])
        lapack = bool(ca) and bool(cb) and ca != cb and set(num["columns"]) <= listed  # S1
        return d | num | {"class": "numeric_only", "sub_class": "known_lapack" if lapack else None}
    line = next((i + 1 for i, (x, y) in enumerate(zip_longest(la, lb)) if x != y), None)
    return d | {"class": "content", "first_differing_line": line}


def classify_tree(
    files_a: Mapping[str, Reader],
    files_b: Mapping[str, Reader],
    known_columns: Mapping[str, Collection[str]],
    producers: Mapping[str, Any],
    *,
    labels: Sequence[str] = ("ref", "rep1"),
) -> dict[str, Any]:
    """te7.json's ``giap`` object: the byte comparison of the gate files decides TE7. The
    te7_files.tsv rows ride along under ``_rows`` (popped before writing)."""
    counts: Counter[str] = Counter()
    found: dict[str, list[Any]] = {"gate": [], "external": [], "figures": []}
    rows: list[list[Any]] = []
    for p in sorted(set(files_a) | set(files_b)):
        cls = path_class(p)
        counts[cls] += 1
        if cls == "excluded":
            continue
        ba = files_a[p]() if p in files_a else None
        bb = files_b[p]() if p in files_b else None
        ha, hb = (None if x is None else hashlib.sha256(x).hexdigest() for x in (ba, bb))
        rows.append([p, cls, ha, hb, ha is not None and ha == hb])
        if ha is not None and ha == hb:
            continue
        found[cls].append(
            _diagnose(p, ba, bb, labels, known_columns, producers.get(p) or {})
            if cls == "gate"
            else {"path": p}
        )
    return {
        "status": "fail" if found["gate"] else "pass",
        "differing": found["gate"],
        "external_differing": found["external"],
        "figures_differing": found["figures"],
        **{f"n_{c}": counts[c] for c in ("gate", "external", "figures", "excluded")},
        "_rows": rows,
    }


# ---------------------------------------------------------------- DIA-NN layer


def _keyed(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    for c in ("run", "modified_peptide", "charge"):
        d[c] = d[c].map(_s)
    d["is_decoy"] = [_true(v) for v in d["is_decoy"]]
    for c in ("q_value", "intensity"):
        d[c] = pd.to_numeric(d[c], errors="coerce").astype(float)
    d["acc"] = ~d["is_decoy"] & (d["q_value"] <= pm.Q_ACCEPT)
    return d.drop_duplicates(L3_KEY).set_index(L3_KEY)


def _jaccard(a: set[Any], b: set[Any]) -> float | None:
    return len(a & b) / len(a | b) if a | b else None


def _differs(x: pd.Series, y: pd.Series) -> int:
    u, v = x.to_numpy(dtype=float), y.to_numpy(dtype=float)
    return int((~((u == v) | (np.isnan(u) & np.isnan(v)))).sum())


def _l3(conv_a: pd.DataFrame, conv_b: pd.DataFrame) -> dict[str, Any]:
    ia, ib = _keyed(conv_a), _keyed(conv_b)
    both = ia.index.intersection(ib.index)
    A, B = ia.loc[both], ib.loc[both]
    xa, xb = A["intensity"].to_numpy(), B["intensity"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.log2(xb / xa)
    ok = (xa > 0) & (xb > 0) & np.isfinite(ratio)
    return {
        "only_a": len(ia) - len(both),
        "only_b": len(ib) - len(both),
        "both": len(both),
        "jaccard_precursors": _jaccard(set(ia.index[ia["acc"]]), set(ib.index[ib["acc"]])),
        "jaccard_peptides": _jaccard(
            set(ia.loc[ia["acc"], "peptide"]), set(ib.loc[ib["acc"], "peptide"])
        ),
        "q_value_differs": _differs(A["q_value"], B["q_value"]),
        "acceptance_flips": int((A["acc"].to_numpy() != B["acc"].to_numpy()).sum()),
        "intensity_differs": _differs(A["intensity"], B["intensity"]),
        "max_abs_log2_intensity_ratio": float(np.abs(ratio[ok]).max()) if ok.any() else None,
    }


def diann_layer(
    sums_a: Mapping[str, str],
    sums_b: Mapping[str, str],
    conv_a: Frame,
    conv_b: Frame,
    *,
    conv_sha: tuple[str, str],
) -> dict[str, Any]:
    """L0 (file list), L1 (report sha256), L2 (converted table sha256), L3 (aggregates). The
    converted tables may be frames or zero-argument readers; they are read only for L3."""
    names = sorted(set(sums_a) | set(sums_b))
    equal = [n for n in names if n in sums_a and sums_a[n] == sums_b.get(n)]
    l0 = {"files_equal": equal, "files_differ": [n for n in names if n not in equal]}
    report = [n for n in names if n.endswith("report.parquet")]
    if report and set(report) <= set(equal):
        return {"status": "identical_report", "l0": l0, "l3": None}
    if conv_sha[0] == conv_sha[1]:
        return {"status": "identical_converted", "l0": l0, "l3": None}
    a, b = (c() if callable(c) else c for c in (conv_a, conv_b))
    if a is None or b is None:
        raise ValueError("diann_layer: converted tables are needed when L2 differs")
    return {"status": "nondeterministic", "l0": l0, "l3": _l3(a, b)}


# ---------------------------------------------------------------- pair checks


def pair_problems(ref: Mapping[str, Any], rep: Mapping[str, Any]) -> list[str]:
    """The manifest fields that stop two runs from forming a TE7 pair (empty: a pair). Both runs
    must be fresh; outdir, command line and timestamps may differ (top level, not compared)."""
    fresh: dict[str, Callable[[Mapping[str, Any]], bool]] = {
        "exit_code": lambda m: m.get("exit_code") == 0,
        "resume": lambda m: m.get("resume") is False,
        "git_dirty": lambda m: m.get("git_dirty") is False,
        "steps.cached": lambda m: all(s.get("status") != "cached" for s in m.get("steps", [])),
    }
    same: dict[str, Callable[[Mapping[str, Any]], Any]] = {
        "inputs.sha256": lambda m: sorted((i.get("role"), i.get("sha256"))
                                          for i in m.get("inputs", [])),
        "seeds": lambda m: m.get("seeds"),
        "databases": lambda m: sorted((d.get("name"), d.get("version"), d.get("checksum"))
                                      for d in m.get("databases") or []),
        "config": lambda m: m.get("config"),
    }  # fmt: skip
    bad = [f"{lab}.{k}" for lab, m in (("ref", ref), ("rep", rep))
           for k, ok in fresh.items() if not ok(m)]  # fmt: skip
    return bad + [k for k, f in same.items() if f(ref) != f(rep)]


# ---------------------------------------------------------------- file I/O


def _tree(root: Path) -> dict[str, Reader]:
    out: dict[str, Reader] = {}
    for d, _dirs, names in os.walk(root):
        for n in names:
            p = Path(d) / n
            if p.is_file():
                out[p.relative_to(root).as_posix()] = p.read_bytes
    return out


def _step_classes(out: Path) -> tuple[list[str], dict[str, list[str]]]:
    """Trace hosts and, per registry step, its sorted CPU classes (``VARIANTS_STEP`` -> variants,
    ``REPORT_HTML`` -> report; S3)."""
    from scripts.measure_runtimes import node_class, process_name

    hosts: set[str] = set()
    per: dict[str, set[str]] = {}
    for r in pm._rows(out / "logs" / "trace.tsv"):
        h = (r.get("hostname") or "").strip()
        if h and h != "-":
            name = process_name(r.get("name", ""))
            step = STEP_ALIASES.get(name, name.removesuffix("_step"))
            hosts.add(h)
            per.setdefault(step, set()).add(node_class(h))
    return sorted(hosts), {k: sorted(v) for k, v in sorted(per.items())}


def _producers(
    paths: Iterable[str], classes: Mapping[str, Mapping[str, list[str]]]
) -> dict[str, Any]:
    """OUT path -> its producing registered step (output key = file name with ``_`` for ``.``, or
    its stem) and that step's CPU classes per run."""
    import giap.steps
    from giap.steps.registry import all_steps

    giap.steps.load_step_modules()
    owner = {k: s.name for s in all_steps() for k in s.outputs}
    out: dict[str, Any] = {}
    for p in paths:
        name = PurePosixPath(p).name
        step = owner.get(name.replace(".", "_")) or owner.get(name.split(".", 1)[0])
        if step is not None:
            cc = {k: list(v.get(step, [])) for k, v in classes.items()}
            out[p] = {"producer": step, "cpu_classes": cc}
    return out


def _known_columns() -> dict[str, set[str]]:
    text = Path(__file__).with_name("te7_known_columns.tsv").read_text(encoding="utf-8")
    out: dict[str, set[str]] = {}
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        out.setdefault(PurePosixPath(r["table"]).name, set()).add(r["column"])
    return out


def _git_equivalent(repo: Path, a: str, b: str) -> bool:
    res = subprocess.run(["git", "-C", str(repo), "diff", "--quiet", a, b, "--", *CODE_PATHS],
                         capture_output=True, text=True, check=False)  # fmt: skip
    if res.returncode not in (0, 1):
        raise pm.Precondition(
            f"git_commit: cannot compare {a} and {b} in {repo}: {res.stderr.strip()}"
        )
    return res.returncode == 0


def _diann(a: argparse.Namespace, pair_peptides: Sequence[str | None]) -> dict[str, Any]:
    given = [a.diann_a, a.diann_b, a.peptides_a, a.peptides_b]
    if not any(given):
        return {"status": "not_run", "reason": "not_given", "l0": None, "l3": None, "nodes": []}
    if not all(given):
        raise pm.Precondition("--diann-a, --diann-b, --peptides-a and --peptides-b go together")
    shas = (pm.sha256(pm._need(a.peptides_a)), pm.sha256(pm._need(a.peptides_b)))
    if any(s != shas[0] for s in pair_peptides):  # S2: end_to_end rests on this
        raise pm.Precondition(f"{a.peptides_a}: sha256 differs from the pair's peptides input")
    sums, nodes = [], []
    for d in (a.diann_a, a.diann_b):
        text = pm._need(d / "outputs.sha256").read_text(encoding="utf-8")
        sums.append({ln.split(None, 1)[1].strip().lstrip("*"): ln.split(None, 1)[0]
                     for ln in text.splitlines() if ln.strip()})  # fmt: skip
        prov = d / "provenance.json"
        nodes.append(pm._json(prov).get("node") if prov.is_file() else None)

    def read(p: Path) -> Callable[[], pd.DataFrame]:
        return lambda: pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)

    conv = (read(a.peptides_a), read(a.peptides_b))
    return diann_layer(sums[0], sums[1], *conv, conv_sha=shas) | {"nodes": nodes}


def te7(a: argparse.Namespace) -> tuple[dict[str, Any], list[list[Any]]]:
    """te7.json and the te7_files.tsv rows (reads only; raises on a precondition)."""
    _out_dir(a.out)
    (ref, ref_job), (rep, rep_job) = pm.resolve_attempt(a.ref), pm.resolve_attempt(a.rep)
    mans = [pm._json(pm._need(d) / "run_manifest.json") for d in (ref, rep)]
    if problems := pair_problems(*mans):
        raise pm.Precondition("not a TE7 pair: " + ", ".join(problems))
    commits = [str(m.get("git_commit") or "") for m in mans]
    if not all(commits) or not _git_equivalent(a.repo, *commits):
        raise pm.Precondition(f"git_commit: {commits[0]} and {commits[1]} are not code-equivalent")
    diann = _diann(a, [pm._input_sha(m, "peptides") for m in mans])
    labels = [m.group(0) if (m := re.search(r"rep\d+", d.name)) else default
              for d, default in ((ref, "ref"), (rep, "rep1"))]  # fmt: skip
    labels = labels if labels[0] != labels[1] else ["ref", "rep1"]
    pair, classes = {}, {}
    for lab, out, man, job in zip(labels, (ref, rep), mans, (ref_job, rep_job), strict=True):
        hosts, classes[lab] = _step_classes(out)
        pair[lab] = {"outdir": str(out.resolve()), "git_commit": man.get("git_commit"), "job": job,
                     "hosts": hosts, "cpu_classes_by_step": classes[lab]}  # fmt: skip
    fa, fb = _tree(ref), _tree(rep)
    gate = [p for p in sorted(set(fa) | set(fb)) if path_class(p) == "gate"]
    giap = classify_tree(fa, fb, _known_columns(), _producers(gate, classes), labels=labels)
    rows = giap.pop("_rows")
    e2e = giap["status"] == "pass" and diann["status"] in (
        "identical_report",
        "identical_converted",
    )
    doc = {"schema": "giap-pilot-te7/1", "spec": _spec_sha(), "pair": pair, "giap": giap,
           "diann": diann, "end_to_end": e2e}  # fmt: skip
    return doc, rows


def run(a: argparse.Namespace) -> dict[str, str]:
    """``te_real.py te7``: write te7.json and te7_files.tsv; the statuses for the summary line."""
    doc, rows = te7(a)
    _write(a.out, "te7.json", doc, "te7_files.tsv", FILES_COLS, rows)
    return {"te7": doc["giap"]["status"], "diann": doc["diann"]["status"]}
