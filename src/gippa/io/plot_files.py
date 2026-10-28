"""Image and data-TSV writers for plot instances (M9-01 behaviour 7; M9-index decision 6).

``write_instance`` writes ``Vnn_<name>[_<suffix>].<fmt>`` for each format and
``data/Vnn_<name>[_<suffix>].tsv`` through ``gippa.io.tsv`` (§9.2, ADR-0006). Image metadata that
carries a date or a software version is dropped, so two renders of one instance are
byte-identical on one platform (``svg.hashsalt`` is fixed here, whether or not
``gippa.viz.style.apply`` ran). matplotlib hashes full-precision floats into SVG clip-path, marker
and hatch ids, so ``stable_svg_ids`` renumbers them and 1-ulp layout noise cannot change the
bytes. An SVG carries the finding title as its ``<title>`` (M9-15). The figures are
``matplotlib.figure.Figure`` objects that pyplot never manages, so nothing needs closing. The
``Instance.coords`` columns are written ``+ 1`` (1-based closed data TSVs, §6.5).
``read_plot_files`` reads them back for the report (M9-02).
"""

from __future__ import annotations

import html
import io
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

import matplotlib
import matplotlib.image
import numpy as np

from gippa.io.tsv import write_tsv
from gippa.viz.base import Instance

#: per-format ``savefig`` metadata that removes dates and software strings
METADATA: dict[str, dict[str, Any]] = {
    "svg": {"Date": None},
    "pdf": {"CreationDate": None, "ModDate": None},
    "png": {"Software": None},
}
#: fixed SVG id salt (``gippa.viz.style`` sets the same value)
#: the pre-rename salt stays: a new salt changes every SVG id (ADR-0064 item 7); the hashed ids
#: are renumbered by ``stable_svg_ids`` anyway
HASHSALT = "giap"
_UNSAFE = re.compile(r"[^A-Za-z0-9._-]")
#: matplotlib's hashed SVG ids (clip paths, markers, hatches) and their renamed stems
_HASHED_ID = re.compile(rb'id="([pmh][0-9a-f]{10})"')
_HASHED_USE = re.compile(rb'(id="|url\(#|href="#)([pmh][0-9a-f]{10})(?=[")])')
_ANY_ID = re.compile(rb'\bid="([^"]*)"')
_ID_STEM = {b"p": b"clip", b"m": b"marker", b"h": b"hatch"}


def stem(pid: str, name: str, suffix: str | None) -> str:
    """``Vnn_<name>`` or ``Vnn_<name>_<suffix>`` with the suffix sanitised to ``[A-Za-z0-9._-]``."""
    return f"{pid}_{name}" if suffix is None else f"{pid}_{name}_{_UNSAFE.sub('_', suffix)}"


def stable_svg_ids(svg: bytes) -> bytes:
    """``svg`` with every hashed id renamed ``clip<n>``/``marker<n>``/``hatch<n>`` (per prefix, in
    order of first definition) and every ``id="…"``, ``url(#…)`` (also inside ``style=``) and
    ``href="#…"`` rewritten; text content is left alone.

    The new names skip any id already in the file, so they cannot collide.
    """
    hashed = list(dict.fromkeys(_HASHED_ID.findall(svg)))
    taken = set(_ANY_ID.findall(svg)) - set(hashed)
    counters: dict[bytes, int] = {}
    mapping: dict[bytes, bytes] = {}
    for old in hashed:
        prefix = old[:1]
        n = counters.get(prefix, 0)
        while (new := _ID_STEM[prefix] + str(n).encode()) in taken:
            n += 1
        counters[prefix] = n + 1
        taken.add(new)
        mapping[old] = new

    def _sub(m: re.Match[bytes]) -> bytes:
        old = m.group(2)
        return m.group(1) + mapping[old] if old in mapping else m.group(0)

    return _HASHED_USE.sub(_sub, svg)


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
            meta = METADATA[fmt] | ({"Title": inst.title} if fmt == "svg" else {})
            if fmt == "svg":
                buf = io.BytesIO()
                inst.figure.savefig(buf, format=fmt, dpi=dpi, metadata=meta)
                path.write_bytes(stable_svg_ids(buf.getvalue()))
            else:
                inst.figure.savefig(path, format=fmt, dpi=dpi, metadata=meta)
            written.append(path)
    tsv = plots_dir / "data" / f"{base}.tsv"
    tsv.parent.mkdir(parents=True, exist_ok=True)
    frame = inst.data.copy()
    for col in inst.coords:
        frame[col] = frame[col] + 1
    write_tsv(tsv, [str(c) for c in frame.columns], frame.itertuples(index=False, name=None))
    written.append(tsv)
    return written


def read_png(src: Path | BinaryIO) -> np.ndarray:
    """A PNG as float RGBA in [0, 1], shape ``(h, w, 4)`` (the M9-24 colour-blind check)."""
    img = np.asarray(matplotlib.image.imread(src, format="png"), dtype=float)
    if img.shape[-1] == 3:
        img = np.concatenate([img, np.ones((*img.shape[:2], 1))], axis=-1)
    return img


def write_png(path: Path, rgba: np.ndarray) -> None:
    """Write float RGBA in [0, 1] as a PNG without a software stamp."""
    path.parent.mkdir(parents=True, exist_ok=True)
    matplotlib.image.imsave(path, rgba, format="png", metadata=METADATA["png"])


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


_SVG_TITLE = re.compile(r"<title>(.*?)</title>", re.S)


def svg_title(svg: str) -> str | None:
    """The finding title ``write_instance`` stores as the SVG's ``<title>`` (M9-15 decision 7)."""
    m = _SVG_TITLE.search(svg)
    return None if m is None else html.unescape(m.group(1))


#: image formats in the order the report links them
IMAGE_FORMATS = ("svg", "png", "pdf")
_PLOT_FILE = re.compile(r"^(V\d{2})_.+$")


@dataclass(frozen=True)
class InstanceFiles:
    """One instance's files under ``plots/``: ``images`` and ``data`` are ``plots/``-relative."""

    stem: str
    svg: str | None  # the SVG text, when the SVG exists
    images: tuple[str, ...]
    data: str | None


@dataclass(frozen=True)
class PlotFiles:
    pid: str
    instances: tuple[InstanceFiles, ...]


def read_plot_files(plots_dir: Path) -> dict[str, PlotFiles]:
    """Every ``Vnn_*`` image and data TSV under ``plots_dir``, by plot ID then stem (M9-02)."""
    stems: dict[str, dict[str, list[str]]] = {}
    found: list[tuple[str, str]] = []
    if plots_dir.is_dir():
        found += [(p.name, p.suffix[1:]) for p in plots_dir.iterdir() if p.is_file()]
    if (plots_dir / "data").is_dir():
        found += [(f"data/{p.name}", "tsv") for p in (plots_dir / "data").glob("*.tsv")]
    for rel, fmt in found:
        name = Path(rel).stem
        if _PLOT_FILE.match(name) and fmt in (*IMAGE_FORMATS, "tsv"):
            stems.setdefault(name, {}).setdefault(fmt, []).append(rel)
    out: dict[str, list[InstanceFiles]] = {}
    for name in sorted(stems):
        by = stems[name]
        svg = by.get("svg")
        inst = InstanceFiles(
            stem=name,
            svg=(plots_dir / svg[0]).read_text(encoding="utf-8") if svg else None,
            images=tuple(by[f][0] for f in IMAGE_FORMATS if f in by),
            data=by["tsv"][0] if "tsv" in by else None,
        )
        out.setdefault(name[:3], []).append(inst)
    return {pid: PlotFiles(pid, tuple(insts)) for pid, insts in sorted(out.items())}
