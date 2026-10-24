"""Shared contract and readers for the M9-17 tests (docs/specs/M9-17-plot-support-analysis.md;
ADR-0043, the "ADR-M9-A" of the spec).

Contract choices the test author made where the spec leaves them open (each one is reported to
the lead):

* **Step wiring.** Each producing step declares one new **optional output** that is a
  directory: ``function_support`` (step ``function``), ``determinants_support`` (step
  ``determinants``) and ``context_stats_support`` (step ``context_stats``). When the key is
  wired the step always creates the directory (``run_step`` checks every declared output
  exists), and writes its support files somewhere under it as ``<name>.tsv`` or
  ``<name>.tsv.gz`` (ADR-0043 rule 2: exactly one of the two; gzip above 2^20 bytes). The tests
  find a file anywhere under the directory, so ``<dir>/<name>.tsv`` and
  ``<dir>/support/<name>.tsv`` (``write_support(name, df, dir)``) both pass. Without the key the
  step writes the same ``tables/`` bytes as before (acceptance "No change to any tables/ file
  bytes").
* **Publishing.** ``giap.runner.local.publish`` puts each file at
  ``OUT/plots/support/<name>.tsv[.gz]``, where ``RunView.from_outdir(OUT).support(name)`` reads
  it.
* **Cells.** §9.2 through ``giap.io.tsv`` (``%.6g`` floats, ``NA``, ``true``/``false``), so
  identities the spec states to 1e-9 in memory are checked on the files to the ``%.6g``
  rounding.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

#: step -> its new optional output key (a directory)
SUPPORT_KEYS: Mapping[str, str] = {
    "function": "function_support",
    "determinants": "determinants_support",
    "context_stats": "context_stats_support",
}

#: support file -> columns, verbatim from the spec's "Support files" table
COLUMNS: Mapping[str, tuple[str, ...]] = {
    "gsea_running": (
        "collection", "term", "rank_position", "tu_id", "score", "in_set", "running_es",
    ),
    "determinants_shap": ("locus_tag", "feature", "feature_value", "shap"),
    "determinants_pd": ("feature", "grid_value", "pd_mean", "pd_low", "pd_high"),
    "null_regions": ("region_id", "contig", "stat", "shift_k", "null_value", "observed"),
    "null_islands": ("contig", "direction", "replicate", "max_scan_s", "observed_max"),
    "null_determinants": ("model", "permutation", "cv_r2", "observed_cv_r2_mean"),
}  # fmt: skip
#: Q2 (lead decision 2): the percentile-grid file names its third column ``percentile``
GRID_COLUMNS: tuple[str, ...] = (
    "region_id", "contig", "stat", "percentile", "null_value", "observed",
)  # fmt: skip
#: Q2: the exhaustive file is written while it is under 5 MB (1 MB = 2^20 bytes, ADR-0043)
NULL_REGIONS_MAX_BYTES = 5 * 2**20
#: the PD features (V15(b))
PD_FEATURES: tuple[str, ...] = ("cai", "tai", "tir_ostir", "promotech_max")
PD_GRID_POINTS = 21  # 5 %, 10 %, ..., 95 % quantiles plus min and max
SHAP_TOP = 20
GSEA_TOP = 6
#: relative rounding of a ``%.6g`` cell (half a unit in the 6th significant digit, doubled)
REL = 1e-5


def find(root: Path, name: str) -> Path | None:
    """The one ``<name>.tsv`` / ``<name>.tsv.gz`` under ``root`` (None when there is none)."""
    if not root.exists():
        return None
    hits = sorted(
        p for p in root.rglob(f"{name}.tsv*") if p.name in (f"{name}.tsv", f"{name}.tsv.gz")
    )
    assert len(hits) <= 1, f"more than one {name} support file: {hits}"
    return hits[0] if hits else None


def need(root: Path, name: str) -> Path:
    path = find(root, name)
    assert path is not None, (
        f"no {name}.tsv[.gz] under {root}: {sorted(p.name for p in root.rglob('*'))}"
        if root.exists()
        else f"support directory {root} was not created"
    )
    return path


def names(root: Path) -> list[str]:
    """Support file stems under ``root``."""
    if not root.exists():
        return []
    return sorted(p.name.split(".")[0] for p in root.rglob("*.tsv*") if p.is_file())


def read(path: Path) -> pd.DataFrame:
    """A support TSV (plain or gz) as strings, ``NA`` verbatim."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, compression="infer")


def header(path: Path) -> list[str]:
    return list(read(path).columns)


def raw_text(path: Path) -> str:
    import gzip

    data = path.read_bytes()
    if path.name.endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode("utf-8")


def check_written_form(path: Path) -> None:
    """ADR-0043 rule 2: ``.tsv.gz`` exactly when the TSV exceeds 2^20 bytes; gzip ``mtime=0``;
    one header row, ``\\n`` endings."""
    text = raw_text(path)
    size = len(text.encode("utf-8"))
    if path.name.endswith(".gz"):
        assert size > 2**20, f"{path.name}: gzipped at {size} bytes"
        assert path.read_bytes()[4:8] == b"\x00\x00\x00\x00", "gzip MTIME must be 0"
    else:
        assert size <= 2**20, f"{path.name}: {size} bytes should be gzipped"
    assert text.endswith("\n") and "\r" not in text
    assert not text.startswith("#"), "§9.2: one header row, no comment lines"


def num(cell: str) -> float:
    if cell in ("", "NA"):
        return math.nan
    if cell == "INF":
        return math.inf
    if cell == "-INF":
        return -math.inf
    return float(cell)


def nums(col: pd.Series) -> np.ndarray:
    return np.array([num(c) for c in col], dtype=float)


def close(a: float, b: float, rel: float = REL, abs_: float = 1e-9) -> bool:
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def step_spec(name: str) -> Any:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    return reg.get(name)


def run_with_spec_outputs(step: str, cfg: Any, inputs: Mapping[str, Path], workdir: Path) -> None:
    """``run_step`` with every declared output at ``workdir / key``, as ``run_local`` does."""
    from giap.steps import registry as reg

    spec = step_spec(step)
    workdir.mkdir(parents=True, exist_ok=True)
    reg.run_step(step, cfg, dict(inputs), {k: workdir / k for k in spec.outputs})


def publish_and_view(workdir: Path, out: Path) -> Any:
    """``publish`` the work directory to ``out`` and open it as ``giap report`` does."""
    from giap.io.run_view import RunView
    from giap.runner.local import publish

    (out / "tables").mkdir(parents=True, exist_ok=True)
    publish(workdir, out)
    return RunView.from_outdir(out)


def tables_bytes(paths: Mapping[str, Path]) -> dict[str, bytes]:
    return {k: p.read_bytes() for k, p in sorted(paths.items()) if p.is_file()}


def tree_bytes(root: Path) -> dict[str, bytes]:
    """Relative path -> bytes of every file under ``root``."""
    return {p.relative_to(root).as_posix(): p.read_bytes()
            for p in sorted(root.rglob("*")) if p.is_file()}  # fmt: skip
