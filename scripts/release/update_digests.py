"""Pin GHCR image digests in ``workflow/conf/containers.config`` (M11-06 item 3; ADR-0008).

Usage:

    python scripts/release/update_digests.py --digests digests.json
        [--config workflow/conf/containers.config] [--nextflow-config workflow/nextflow.config]
        --date YYYY-MM-DD --source-commit <git sha>

``digests.json`` maps process names to ``ghcr.io/khoa-yelo/giap-<tool>@sha256:<64 hex>``. For each
entry the ``withName: '<PROCESS>' { container = 'docker://…' }`` line inside ``process { }`` is
replaced in place, or added before the block's closing brace; the tool leaves the "Conda-only for
now" comment list; one ``//   ghcr.io/khoa-yelo/giap-<tool>:<sha7> (built <date> from <commit>)``
line per image is added to (or replaced in) the "Tags the digests were resolved from" comment.
Every other line stays byte-identical, and a second run with the same input changes nothing.

With ``--nextflow-config`` the newly pinned processes are also appended to the apptainer profile's
``withName: '…' { conda = '' }`` selector, so the pinned set equals the conda-off set (item 4).

The whole batch is validated before anything is written: a malformed digest, a tag instead of a
digest, an unknown process, another tool's or another registry's image, or any UniOP entry (R5)
exits 1 with a one-line reason on stderr and leaves every file untouched.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO / "workflow" / "conf" / "containers.config"
REGISTRY = "ghcr.io/khoa-yelo"

# process -> image tool; the _NOVEL aliases are the second-round runs in
# workflow/subworkflows/local/novel.nf (same module, same image).
PROCESS_TOOL: dict[str, str] = {
    "PROMOTECH": "promotech",
    "PROMOTECH_NOVEL": "promotech",
    "PROMOTER_CALCULATOR": "promoter_calculator",
    "TERMINATORNET": "terminatornet",
    "BACTERMFINDER": "bactermfinder",
    "OSTIR": "ostir",
    "OSTIR_NOVEL": "ostir",
}
TOOL_ORDER: tuple[str, ...] = tuple(dict.fromkeys(PROCESS_TOOL.values()))

_DIGEST = re.compile(r".+@sha256:[0-9a-f]{64}")
_CONDA_HEADER = "// Conda-only for now"
_CONDA_LIST = "//   uniop"
_TAGS_HEADER = "// Tags the digests were resolved from"
_CONT = "//   "
_CONDA_OFF = re.compile(
    r"^(?P<pre>\s*withName: ')(?P<names>[A-Z0-9_|]+)(?P<post>' \{ conda = '' \})$"
)


class DigestError(ValueError):
    """One bad entry: the whole batch is rejected."""


def validate(digests: Mapping[str, object]) -> dict[str, str]:
    """Return the pins in canonical process order, or raise DigestError for the first bad entry."""
    for process, ref in digests.items():
        if process.upper() == "UNIOP":
            raise DigestError("UNIOP must not be pinned to a public image (unlicensed; R5)")
        if process not in PROCESS_TOOL:
            known = ", ".join(PROCESS_TOOL)
            raise DigestError(f"unknown process {process} (pinnable: {known})")
        if not isinstance(ref, str) or not _DIGEST.fullmatch(ref):
            raise DigestError(f"bad digest for {process}: {ref!r} (need <image>@sha256:<64 hex>)")
        expected = f"{REGISTRY}/giap-{PROCESS_TOOL[process]}@sha256:"
        if not ref.startswith(expected):
            raise DigestError(f"wrong image for {process}: {ref} (need {expected}<64 hex>)")
    return {p: str(digests[p]) for p in PROCESS_TOOL if p in digests}


def pin_line(process: str, ref: str) -> str:
    return f"  withName: '{process}' {{ container = 'docker://{ref}' }}"


def _process_block(lines: list[str]) -> tuple[int, int]:
    try:
        start = lines.index("process {")
    except ValueError:
        raise DigestError("no `process {` block in the config") from None
    end = next((i for i in range(start + 1, len(lines)) if lines[i] == "}"), None)
    if end is None:
        raise DigestError("unterminated `process {` block in the config")
    return start, end


def _set_pins(lines: list[str], pins: Mapping[str, str]) -> list[str]:
    out = list(lines)
    for process, ref in pins.items():
        start, end = _process_block(out)
        pat = re.compile(rf"^\s*withName:\s*'{process}'\s*\{{")
        hits = [i for i in range(start + 1, end) if pat.match(out[i])]
        if len(hits) > 1:
            raise DigestError(f"{process} is pinned more than once in the config")
        if hits:
            out[hits[0]] = pin_line(process, ref)
        else:
            out.insert(end, pin_line(process, ref))
    return out


def _drop_conda_only(lines: list[str], tools: set[str]) -> list[str]:
    out = list(lines)
    header = next((i for i, ln in enumerate(out) if ln.startswith(_CONDA_HEADER)), None)
    if header is None:
        raise DigestError("no `Conda-only for now` comment in the config")
    at = next((i for i in range(header, len(out)) if out[i].startswith(_CONDA_LIST)), None)
    if at is None:
        raise DigestError("no conda-only tool list (`//   uniop, …`) in the config")
    names = [n.strip() for n in out[at][len(_CONT) :].split(",") if n.strip()]
    out[at] = _CONT + ", ".join(n for n in names if n not in tools)
    return out


def _record_tags(lines: list[str], tools: list[str], date: str, commit: str) -> list[str]:
    out = list(lines)
    header = next((i for i, ln in enumerate(out) if ln.startswith(_TAGS_HEADER)), None)
    if header is None:
        raise DigestError("no `Tags the digests were resolved from` comment in the config")
    for tool in tools:
        prefix = f"{_CONT}{REGISTRY}/giap-{tool}:"
        line = f"{prefix}{commit[:7]} (built {date} from {commit})"
        end = header + 1
        while end < len(out) and out[end].startswith(_CONT):
            end += 1
        existing = [i for i in range(header + 1, end) if out[i].startswith(prefix)]
        if existing:
            out[existing[0]] = line
        else:
            out.insert(end, line)
    return out


def update_containers(text: str, pins: Mapping[str, str], date: str, commit: str) -> str:
    """Return the new containers.config text for already validated pins."""
    lines = text.splitlines()
    tools = [t for t in TOOL_ORDER if any(PROCESS_TOOL[p] == t for p in pins)]
    lines = _set_pins(lines, pins)
    lines = _drop_conda_only(lines, set(tools))
    lines = _record_tags(lines, tools, date, commit)
    return "\n".join(lines) + "\n"


def update_conda_off(text: str, processes: list[str]) -> str:
    """Append processes to the single `withName: '…' { conda = '' }` selector in nextflow.config."""
    lines = text.splitlines()
    hits = [i for i, ln in enumerate(lines) if _CONDA_OFF.match(ln)]
    if len(hits) != 1:
        raise DigestError(
            f"expected one `{{ conda = '' }}` selector in nextflow.config, found {len(hits)}"
        )
    m = _CONDA_OFF.match(lines[hits[0]])
    assert m is not None
    names = m.group("names").split("|")
    names += [p for p in processes if p not in names]
    lines[hits[0]] = m.group("pre") + "|".join(names) + m.group("post")
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def _write(path: Path, text: str) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    os.chmod(tmp, path.stat().st_mode & 0o777)
    os.replace(tmp, path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--digests", type=Path, required=True)
    ap.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    ap.add_argument("--nextflow-config", type=Path, default=None)
    ap.add_argument("--date", required=True, help="build date, YYYY-MM-DD (UTC)")
    ap.add_argument("--source-commit", required=True, help="the commit the images were built from")
    args = ap.parse_args(argv)
    try:
        try:
            dt.date.fromisoformat(args.date)
        except ValueError:
            raise DigestError(f"bad --date {args.date!r} (need YYYY-MM-DD)") from None
        if not re.fullmatch(r"[0-9a-f]{7,40}", args.source_commit):
            raise DigestError(f"bad --source-commit {args.source_commit!r} (need a hex git sha)")
        raw = json.loads(args.digests.read_text())
        if not isinstance(raw, dict) or not raw:
            raise DigestError("digests.json must be a non-empty object {process: image@sha256:…}")
        pins = validate(raw)
        new_cfg = update_containers(args.config.read_text(), pins, args.date, args.source_commit)
        new_nf = None
        if args.nextflow_config is not None:
            new_nf = update_conda_off(args.nextflow_config.read_text(), list(pins))
    except (DigestError, OSError, json.JSONDecodeError) as exc:
        print(f"update_digests: {exc}", file=sys.stderr)
        return 1
    if new_cfg != args.config.read_text():
        _write(args.config, new_cfg)
    nf = args.nextflow_config
    if new_nf is not None and nf is not None and new_nf != nf.read_text():
        _write(nf, new_nf)
    for process, ref in pins.items():
        print(f"pinned {process} {ref}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
