"""Image and data-TSV writers for plot instances (M9-01 behaviour 7; M9-index decision 6).

``write_instance`` writes ``Vnn_<name>[_<suffix>].<fmt>`` for each format and
``data/Vnn_<name>[_<suffix>].tsv`` through ``giap.io.tsv`` (§9.2, ADR-0006). Image metadata that
carries a date or a software version is dropped, so two renders of one instance are
byte-identical on one platform (``svg.hashsalt`` is fixed here, whether or not
``giap.viz.style.apply`` ran). The figures are ``matplotlib.figure.Figure`` objects that pyplot
never manages, so nothing needs closing. The ``Instance.coords`` columns are written ``+ 1``
(1-based closed data TSVs, §6.5).
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib

from giap.io.tsv import write_tsv
from giap.viz.base import Instance

#: per-format ``savefig`` metadata that removes dates and software strings
METADATA: dict[str, dict[str, Any]] = {
    "svg": {"Date": None},
    "pdf": {"CreationDate": None, "ModDate": None},
    "png": {"Software": None},
}
#: fixed SVG id salt (``giap.viz.style`` sets the same value)
HASHSALT = "giap"
_UNSAFE = re.compile(r"[^A-Za-z0-9._-]")


def stem(pid: str, name: str, suffix: str | None) -> str:
    """``Vnn_<name>`` or ``Vnn_<name>_<suffix>`` with the suffix sanitised to ``[A-Za-z0-9._-]``."""
    return f"{pid}_{name}" if suffix is None else f"{pid}_{name}_{_UNSAFE.sub('_', suffix)}"


def write_instance(
    inst: Instance, pid: str, name: str, plots_dir: Path, formats: Sequence[str], dpi: int
) -> list[Path]:
    """Write the instance's images and data TSV; returns every path written."""
    unknown = [f for f in formats if f not in METADATA]
    if unknown:
        raise ValueError(f"unknown plot formats {unknown}; known: {sorted(METADATA)}")
    base = stem(pid, name, inst.suffix)
    plots_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    with matplotlib.rc_context({"svg.hashsalt": HASHSALT}):
        for fmt in formats:
            path = plots_dir / f"{base}.{fmt}"
            inst.figure.savefig(path, format=fmt, dpi=dpi, metadata=METADATA[fmt])
            written.append(path)
    tsv = plots_dir / "data" / f"{base}.tsv"
    tsv.parent.mkdir(parents=True, exist_ok=True)
    frame = inst.data.copy()
    for col in inst.coords:
        frame[col] = frame[col] + 1
    write_tsv(tsv, [str(c) for c in frame.columns], frame.itertuples(index=False, name=None))
    written.append(tsv)
    return written


def remove_plot_files(plots_dir: Path, pid: str) -> None:
    """Remove every image and data TSV of ``pid`` (all instances); nothing else is touched."""
    targets: list[tuple[Path, frozenset[str]]] = [
        (plots_dir, frozenset(METADATA)),
        (plots_dir / "data", frozenset({"tsv"})),
    ]
    for folder, suffixes in targets:
        if not folder.is_dir():
            continue
        for path in folder.glob(f"{pid}_*"):
            if path.is_file() and path.suffix.lstrip(".") in suffixes:
                path.unlink()
