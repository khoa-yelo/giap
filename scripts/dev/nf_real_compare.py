"""Compare one real Nextflow process run with its recording (scripts/slurm/nf_process_real.sbatch).

    nf_real_compare.py PROCESS OUT REC     # prints "ok" or "mismatch (<reason>)"

OUT holds the run's published outputs; REC the recorded outputs in the same layout. Every recorded
file must exist in OUT with byte-identical content, except these documented normalisations:
  BAKTA        bakta/<sample>.json: features without "id" (sorted canonical JSON) and the replicons
               (id, type, topology, length); the other Bakta files embed run dates and paths
  GENOMAD      genomad/genome_summary/ only; JSON files without "start_time"
  ANTISMASH    antismash/bakta.json without "timings"; every "enabled_types" and "gene_functions"
               list and every list in a "definition_domains" map sorted (antiSMASH builds them from
               Python sets: hash-seed order; seen on nissle1917); the GenBank/zip embed run dates
  TRANSTERMHP  without the "--pval-conf=<path>" header line (the module resolves expterm.dat itself)
Stdlib only.
"""

from __future__ import annotations

import json
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

Norm = Callable[[str, bytes], bytes]


def _canon(x: Any) -> str:
    return json.dumps(x, sort_keys=True)


def _files(root: Path) -> dict[str, Path]:
    return {str(p.relative_to(root)): p for p in sorted(root.rglob("*")) if p.is_file()}


def compare_trees(got: Path, want: Path, norm: Norm | None = None) -> str:
    """'' if every file under `want` is under `got` with equal (normalised) bytes, else why not."""
    g, w = _files(got), _files(want)
    if not w:
        return f"no recorded files under {want}"
    missing = sorted(set(w) - set(g))
    if missing:
        return f"{len(missing)} of {len(w)} recorded files missing, e.g. {missing[0]}"
    for name, p in w.items():
        a, b = g[name].read_bytes(), p.read_bytes()
        if norm is not None:
            a, b = norm(name, a), norm(name, b)
        if a != b:
            return f"{name} differs"
    return ""


def drop_json_key(key: str) -> Norm:
    def norm(name: str, data: bytes) -> bytes:
        if not name.endswith(".json"):
            return data
        d = json.loads(data)
        d.pop(key, None)
        return _canon(d).encode()

    return norm


def drop_pval_conf(name: str, data: bytes) -> bytes:
    return b"".join(x for x in data.splitlines(True) if not x.startswith(b"--pval-conf="))


def _replicons(d: dict[str, Any]) -> list[tuple[Any, ...]]:
    return [(s["id"], s.get("type"), s.get("topology"), s["length"]) for s in d["sequences"]]


def _features(d: dict[str, Any]) -> list[str]:
    return sorted(_canon({k: v for k, v in f.items() if k != "id"}) for f in d["features"])


_SET_LISTS = ("enabled_types", "gene_functions")  # written from Python sets: order varies by run


def _sort_set_value(k: str, v: Any) -> Any:
    if k in _SET_LISTS:
        return sorted(v)
    if k == "definition_domains" and isinstance(v, dict):  # {rule: [domain, ...]} from sets
        return {r: sorted(d) for r, d in v.items()}
    return _sort_enabled_types(v)


def _sort_enabled_types(x: Any) -> Any:
    if isinstance(x, dict):
        return {k: _sort_set_value(k, v) for k, v in x.items()}
    if isinstance(x, list):
        return [_sort_enabled_types(v) for v in x]
    return x


def antismash_json(data: bytes) -> str:
    d = json.loads(data)
    d.pop("timings", None)
    return _canon(_sort_enabled_types(d))


def compare_bakta(out: Path, rec: Path) -> str:
    runs = sorted(out.glob("bakta/*.json"))
    if len(runs) != 1:
        return f"expected one bakta/*.json, found {len(runs)}"
    got, want = json.loads(runs[0].read_text()), json.loads((rec / runs[0].name).read_text())
    if _replicons(got) != _replicons(want):
        return f"replicons {_replicons(got)} vs recording {_replicons(want)}"
    g, w = _features(got), _features(want)
    if g == w:
        return ""
    return f"features: run {len(g)}, recording {len(w)}, {len(set(g) ^ set(w))} differ"


def compare(proc: str, out: Path, rec: Path) -> str:
    if proc == "BAKTA":
        return compare_bakta(out, rec)
    if proc == "GENOMAD":
        sub = "genomad/genome_summary"
        return compare_trees(out / sub, rec / sub, drop_json_key("start_time"))
    if proc == "ANTISMASH":
        name = "antismash/bakta.json"
        a, b = (antismash_json((d / name).read_bytes()) for d in (out, rec))
        return "" if a == b else f"{name} differs"
    if proc == "TRANSTERMHP":
        return compare_trees(out, rec, drop_pval_conf)
    return compare_trees(out, rec)


def main(argv: list[str]) -> None:
    proc, out, rec = argv[1], Path(argv[2]), Path(argv[3])
    try:
        reason = compare(proc, out, rec)
    except (OSError, ValueError, KeyError) as e:
        reason = f"{type(e).__name__}: {e}"
    print("ok" if not reason else f"mismatch ({reason})")


if __name__ == "__main__":
    main(sys.argv)
