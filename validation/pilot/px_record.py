"""The ProteomeXchange record of PXD078459 -> the three wild-type runs to download (M10-01).

Stdlib only, so ``python -I validation/pilot/px_record.py RECORD.xml --out MANIFEST.tsv`` runs on
any node. ``parse_record`` and ``select_wild_type`` are pure; only ``main`` reads and writes files.

Exit codes of the CLI: 0 manifest written; 2 fewer than three wild-type runs, or a record that
names no timsTOF instrument (no manifest is written then).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

MANIFEST_COLUMNS = (
    "run",
    "file_name",
    "url",
    "url_record",
    "bytes",
    "checksum_record",
    "sha256",
    "sample",
    "selected_because",
)
_CHECKSUMS = {"md5": "md5", "sha-1": "sha1", "sha1": "sha1", "sha-256": "sha256"}
_WT = re.compile(r"(?<![a-z])(wt|wildtype|wild[ _-]type)(?![a-z])", re.I)
_MOD_TOKENS = {"delta", "ko", "mut", "cured"}
_RUN_DIR = re.compile(r"^(.*?\.d)/", re.I)


class SelectionError(ValueError):
    """The record does not hold three wild-type runs."""


@dataclass(frozen=True)
class RecordFile:
    file_name: str
    uri: str
    kind: str  # "raw" | "not_raw"
    bytes: int | None
    checksum: str | None  # "<algorithm>:<lowercase hex>"
    sample: str | None
    annotations: tuple[tuple[str, str], ...] = ()  # (characteristic, value) pairs


@dataclass(frozen=True)
class Selected:
    run: str
    file_name: str
    selected_because: str
    file: RecordFile


@dataclass(frozen=True)
class Rejected:
    file_name: str
    reason: str  # "not_wild_type" | "beyond_three" | "not_raw"


@dataclass(frozen=True)
class Selection:
    selected: tuple[Selected, ...]
    rejected: tuple[Rejected, ...]


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _params(elem: ET.Element) -> list[tuple[str, str]]:
    return [
        (p.get("name", ""), p.get("value", ""))
        for p in elem
        if _local(p.tag) in ("cvParam", "userParam")
    ]


def parse_record(xml: bytes) -> list[RecordFile]:
    """Every ``DatasetFile`` of a ProteomeXchange XML record, in record order."""
    root = ET.fromstring(xml)
    files = []
    for elem in root.iter():
        if _local(elem.tag) != "DatasetFile":
            continue
        uri, kind, size, checksum = "", "not_raw", None, None
        notes: list[tuple[str, str]] = []
        for name, value in _params(elem):
            key = name.strip().lower()
            if key.endswith("uri"):
                if not uri or "raw file" in key:
                    uri = value
                if "raw file" in key:
                    kind = "raw"
            elif key == "file size" and value.strip().isdigit():
                size = int(value)
            elif key in _CHECKSUMS and value.strip():
                checksum = f"{_CHECKSUMS[key]}:{value.strip().lower()}"
            elif m := re.fullmatch(r"characteristics\[(.+)\]", key):
                notes.append((m.group(1).strip(), value))
        file_name = elem.get("name") or uri.rstrip("/").rsplit("/", 1)[-1]
        sample = "; ".join(f"characteristics[{k}]={v}" for k, v in notes) or None
        files.append(RecordFile(file_name, uri, kind, size, checksum, sample, tuple(notes)))
    return files


def instruments(xml: bytes) -> list[str]:
    """The names of the record's ``Instrument`` cvParams."""
    root = ET.fromstring(xml)
    return [
        name
        for elem in root.iter()
        if _local(elem.tag) == "Instrument"
        for name, _ in _params(elem)
    ]


def natural_key(name: str) -> list[int | str]:
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", name)]


def _is_modified(text: str) -> bool:
    if "Δ" in text or "::" in text:
        return True
    return any(t.lower() in _MOD_TOKENS for t in re.split(r"[^0-9A-Za-z]+", text))


def _is_wild_type(text: str) -> bool:
    return bool(_WT.search(text)) and not _is_modified(text)


def _unit(file_name: str) -> str:
    """The run a file belongs to: its ``<name>.d/`` folder, else the file itself (rule 4)."""
    m = _RUN_DIR.match(file_name)
    return m.group(1) if m else file_name


def _wild_type_reason(unit: str, members: Sequence[RecordFile]) -> str | None:
    genotypes = [(k, v) for f in members for k, v in f.annotations if k.lower() == "genotype"]
    if genotypes:  # the annotation wins over the file name
        ok = all(_is_wild_type(v) for _, v in genotypes)
        return f"sdrf:{genotypes[0][0]}={genotypes[0][1]}" if ok else None
    m = _WT.search(unit)
    return f"token:{m.group(1)}" if m and not _is_modified(unit) else None


def select_wild_type(files: Sequence[RecordFile]) -> Selection:
    """The first three wild-type runs by natural sort; every other file with its reason."""
    rejected = [Rejected(f.file_name, "not_raw") for f in files if f.kind != "raw"]
    units: dict[str, list[RecordFile]] = {}
    for f in files:
        if f.kind == "raw":
            units.setdefault(_unit(f.file_name), []).append(f)
    candidates: list[tuple[str, str]] = []
    for unit in sorted(units, key=natural_key):
        reason = _wild_type_reason(unit, units[unit])
        if reason is None:
            rejected += [Rejected(f.file_name, "not_wild_type") for f in units[unit]]
        else:
            candidates.append((unit, reason))
    if len(candidates) < 3:
        names = ", ".join(u for u, _ in candidates)
        raise SelectionError(f"fewer than three wild-type runs: {names}")
    selected = []
    for i, (unit, reason) in enumerate(candidates[:3], start=1):
        selected += [Selected(f"WT_{i}", f.file_name, reason, f) for f in units[unit]]
    for unit, _ in candidates[3:]:
        rejected += [Rejected(f.file_name, "beyond_three") for f in units[unit]]
    return Selection(tuple(selected), tuple(rejected))


def _cell(value: object) -> str:
    text = "NA" if value is None or value == "" else str(value)
    return re.sub(r"[\t\r\n]+", " ", text)


def manifest_rows(selection: Selection) -> list[dict[str, str]]:
    return [
        {
            "run": s.run,
            "file_name": _cell(s.file_name),
            "url": _cell(s.file.uri),
            "url_record": _cell(s.file.uri),
            "bytes": _cell(s.file.bytes),
            "checksum_record": _cell(s.file.checksum),
            "sha256": "NA",
            "sample": _cell(s.file.sample),
            "selected_because": _cell(s.selected_because),
        }
        for s in selection.selected
    ]


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("record", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    xml = args.record.read_bytes()
    names = instruments(xml)
    print(f"instrument: {', '.join(names) or 'NA'}")
    if not any("timstof" in n.lower() for n in names):
        print("the record does not describe timsTOF data", file=sys.stderr)
        return 2
    try:
        selection = select_wild_type(parse_record(xml))
    except SelectionError as err:
        print(err, file=sys.stderr)
        return 2
    for s in selection.selected:
        print(f"selected {s.run}\t{s.file_name}\t{s.selected_because}")
    for r in selection.rejected:
        print(f"not selected\t{r.file_name}\t{r.reason}")
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE
        )
        writer.writeheader()
        writer.writerows(manifest_rows(selection))
    return 0


if __name__ == "__main__":
    sys.exit(main())
