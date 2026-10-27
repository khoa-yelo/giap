"""Pin GHCR image digests in ``workflow/conf/containers.config`` (M11-06 item 3; ADR-0008).

Usage:

    python scripts/release/update_digests.py --digests digests.json
        [--config workflow/conf/containers.config] [--nextflow-config workflow/nextflow.config]
        --date YYYY-MM-DD --source-commit <git sha> [--provenance prov.json]

``digests.json`` maps process names to ``ghcr.io/khoa-yelo/gippa-<tool>@sha256:<64 hex>``. For each
entry the ``withName: '<PROCESS>' { container = '<ref>' }`` line inside ``process { }`` (a bare
``<ref>``, no ``docker://`` prefix; ADR-0059) is replaced in place, or added before the block's
closing brace; the tool leaves the "Conda-only for now" comment list; one
``//   ghcr.io/khoa-yelo/gippa-<tool>:<sha7> (built <date> from <commit>)`` line per image is added
to (or replaced in) the "Tags the digests were resolved from" comment.
Every other line stays byte-identical, and a second run with the same input changes nothing.

``workflow/conf/images.json`` (``--images``, default next to ``--config``) gets one entry per
pinned tool (M11-14 item 2): ``key``/``tag`` from ``gippa.io.image_key`` over ``--root`` (default
the repository), the bare ``ref``, the ``processes`` whose ``withName`` line carries that ref,
``built_from`` (``--source-commit``) and the registry ``licence``; other entries are kept.
``FGSEA`` has no Nextflow process: its ref goes to ``images.json`` only (``"processes": []``).

``--provenance prov.json`` (``{tool: {"commit": <sha>, "date": "YYYY-MM-DD"}}``) names the build
of an image whose tag already existed (M11-14b): its ``built_from`` and its tag-record line carry
that commit and date instead of ``--source-commit`` and ``--date``.

With ``--nextflow-config`` the newly pinned processes are also appended to the one
``withName: '…' { conda = '' }`` selector, so the pinned set equals the conda-off set (item 4).
The selector is looked up in the given file and in every file it pulls in through
``includeConfig`` (recursively, paths relative to the including file, as Nextflow resolves them):
since M11-14 item 3b it lives in ``workflow/conf/container_common.config``, shared by the
``apptainer`` and ``docker`` profiles. Exactly one selector must exist across those files.

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
sys.path.insert(0, str(REPO / "src"))
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
# images with no Nextflow process: pinned in images.json only (M11-14, fgsea)
IMAGE_ONLY: dict[str, str] = {"FGSEA": "fgsea"}
TOOL_ORDER: tuple[str, ...] = tuple(dict.fromkeys(PROCESS_TOOL.values()))

_DIGEST = re.compile(r".+@sha256:[0-9a-f]{64}")
_CONDA_HEADER = "// Conda-only for now"
_CONDA_LIST = "//   uniop"
_TAGS_HEADER = "// Tags the digests were resolved from"
_CONT = "//   "
_CONDA_OFF = re.compile(
    r"^(?P<pre>\s*withName: ')(?P<names>[A-Z0-9_|]+)(?P<post>' \{ conda = '' \})$"
)
_INCLUDE = re.compile(r"""^\s*includeConfig\s+['"](?P<path>[^'"]+)['"]\s*(?://.*)?$""")


class DigestError(ValueError):
    """One bad entry: the whole batch is rejected."""


def validate(digests: Mapping[str, object]) -> dict[str, str]:
    """Return the pins in canonical process order, or raise DigestError for the first bad entry."""
    for process, ref in digests.items():
        if process.upper() == "UNIOP":
            raise DigestError("UNIOP must not be pinned to a public image (unlicensed; R5)")
        tool = PROCESS_TOOL.get(process) or IMAGE_ONLY.get(process)
        if tool is None:
            known = ", ".join([*PROCESS_TOOL, *IMAGE_ONLY])
            raise DigestError(f"unknown process {process} (pinnable: {known})")
        if not isinstance(ref, str) or not _DIGEST.fullmatch(ref):
            raise DigestError(f"bad digest for {process}: {ref!r} (need <image>@sha256:<64 hex>)")
        expected = f"{REGISTRY}/gippa-{tool}@sha256:"
        if not ref.startswith(expected):
            raise DigestError(f"wrong image for {process}: {ref} (need {expected}<64 hex>)")
    refs: dict[str, set[str]] = {}
    for process, ref in digests.items():
        refs.setdefault(PROCESS_TOOL.get(process) or IMAGE_ONLY[process], set()).add(str(ref))
    if split := sorted(t for t, r in refs.items() if len(r) > 1):
        raise DigestError(f"one image, several digests: {', '.join(split)} (pin one ref per tool)")
    return {p: str(digests[p]) for p in [*PROCESS_TOOL, *IMAGE_ONLY] if p in digests}


def pin_line(process: str, ref: str) -> str:
    return f"  withName: '{process}' {{ container = '{ref}' }}"


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


def _record_tags(
    lines: list[str],
    tools: list[str],
    date: str,
    commit: str,
    provenance: Mapping[str, tuple[str, str]] | None = None,
) -> list[str]:
    out = list(lines)
    prov = provenance or {}
    header = next((i for i, ln in enumerate(out) if ln.startswith(_TAGS_HEADER)), None)
    if header is None:
        raise DigestError("no `Tags the digests were resolved from` comment in the config")
    for tool in tools:
        prefix = f"{_CONT}{REGISTRY}/gippa-{tool}:"
        c, d = prov.get(tool, (commit, date))
        line = f"{prefix}{c[:7]} (built {d} from {c})"
        end = header + 1
        while end < len(out) and out[end].startswith(_CONT):
            end += 1
        existing = [i for i in range(header + 1, end) if out[i].startswith(prefix)]
        if existing:
            out[existing[0]] = line
        else:
            out.insert(end, line)
    return out


def update_containers(
    text: str,
    pins: Mapping[str, str],
    date: str,
    commit: str,
    provenance: Mapping[str, tuple[str, str]] | None = None,
) -> str:
    """Return the new containers.config text for already validated pins.

    ``provenance`` maps a tool to ``(commit, date)`` of an image built in an earlier run (its tag
    already existed); every other tool is recorded with ``date`` and ``commit``.
    """
    lines = text.splitlines()
    pins = {p: r for p, r in pins.items() if p in PROCESS_TOOL}
    if not pins:
        return text
    tools = [t for t in TOOL_ORDER if any(PROCESS_TOOL[p] == t for p in pins)]
    lines = _set_pins(lines, pins)
    lines = _drop_conda_only(lines, set(tools))
    lines = _record_tags(lines, tools, date, commit, provenance)
    return "\n".join(lines) + "\n"


def update_conda_off(text: str, processes: list[str]) -> str:
    """Append processes to the single `withName: '…' { conda = '' }` selector of a config text."""
    lines = text.splitlines()
    hits = [i for i, ln in enumerate(lines) if _CONDA_OFF.match(ln)]
    if len(hits) != 1:
        raise DigestError(
            f"expected one `{{ conda = '' }}` selector in the config, found {len(hits)}"
        )
    m = _CONDA_OFF.match(lines[hits[0]])
    assert m is not None
    names = m.group("names").split("|")
    names += [p for p in processes if p not in names]
    lines[hits[0]] = m.group("pre") + "|".join(names) + m.group("post")
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def update_images(
    text: str | None,
    pins: Mapping[str, str],
    config: str,
    root: Path,
    commit: str,
    provenance: Mapping[str, tuple[str, str]] | None = None,
) -> str:
    """Return the new images.json text: one entry per pinned tool, other entries kept.

    ``built_from`` is ``commit``, or the commit in ``provenance`` for an image built earlier.
    """
    prov = provenance or {}
    from gippa.annotate.tools import TOOLS
    from gippa.io.image_key import key, tag

    images: dict[str, object] = json.loads(text) if text else {}
    pinned = re.compile(r"^\s*withName: '([A-Z0-9_]+)' \{ container = '([^']+)' \}$")
    lines = [m for ln in config.splitlines() if (m := pinned.match(ln))]
    refs: dict[str, str] = {}  # tool -> ref of its first pinned process (canonical order)
    for process, ref in pins.items():
        refs.setdefault(PROCESS_TOOL.get(process) or IMAGE_ONLY[process], ref)
    for tool, ref in refs.items():
        k = key(tool, root)
        images[tool] = {
            "key": k,
            "tag": tag(TOOLS[tool].version, k),
            "ref": ref,
            "processes": sorted(m.group(1) for m in lines if m.group(2) == ref),
            "built_from": prov.get(tool, (commit, ""))[0],
            "licence": TOOLS[tool].licence,
        }
    return json.dumps(dict(sorted(images.items())), indent=2, ensure_ascii=False) + "\n"


def parse_provenance(raw: object) -> dict[str, tuple[str, str]]:
    """``{tool: {"commit": <hex sha>, "date": YYYY-MM-DD}}`` -> ``{tool: (commit, date)}``."""
    if not isinstance(raw, dict):
        raise DigestError("provenance must be an object {tool: {commit, date}}")
    out: dict[str, tuple[str, str]] = {}
    for tool, entry in raw.items():
        if tool not in TOOL_ORDER and tool not in IMAGE_ONLY.values():
            raise DigestError(f"provenance for unknown tool {tool!r}")
        if not isinstance(entry, dict):
            raise DigestError(f"provenance of {tool} must be {{commit, date}}")
        commit, date = str(entry.get("commit", "")), str(entry.get("date", ""))
        if not re.fullmatch(r"[0-9a-f]{7,40}", commit):
            raise DigestError(f"bad provenance commit for {tool}: {commit!r}")
        try:
            dt.date.fromisoformat(date)
        except ValueError:
            raise DigestError(f"bad provenance date for {tool}: {date!r}") from None
        out[str(tool)] = (commit, date)
    return out


def read_provenance(path: Path) -> dict[str, tuple[str, str]]:
    """Read a ``--provenance`` file (see :func:`parse_provenance`)."""
    return parse_provenance(json.loads(path.read_text()))


def count_conda_off(text: str) -> int:
    """Number of `withName: '…' { conda = '' }` selector lines in a config text."""
    return sum(1 for ln in text.splitlines() if _CONDA_OFF.match(ln))


def included_configs(text: str) -> list[str]:
    """The ``includeConfig '<path>'`` targets of a config text, in order (``//`` lines skipped)."""
    return [m.group("path") for ln in text.splitlines() if (m := _INCLUDE.match(ln))]


def config_closure(root: Path) -> list[Path]:
    """``root`` and every config it includes, recursively, each once, in first-seen order."""
    seen: list[Path] = []
    todo = [root.resolve()]
    while todo:
        path = todo.pop(0)
        if path in seen:
            continue
        if not path.is_file():
            raise DigestError(f"includeConfig target not found: {path}")
        seen.append(path)
        todo.extend((path.parent / inc).resolve() for inc in included_configs(path.read_text()))
    return seen


def selector_file(nextflow_config: Path) -> Path:
    """The one file among ``nextflow_config`` and its includes that holds the conda-off selector."""
    files = config_closure(nextflow_config)
    counts = {p: count_conda_off(p.read_text()) for p in files}
    total = sum(counts.values())
    if total != 1:
        raise DigestError(
            f"expected one `{{ conda = '' }}` selector in {nextflow_config} and its "
            f"includeConfig files, found {total}"
        )
    return next(p for p, n in counts.items() if n == 1)


def _write(path: Path, text: str) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    os.chmod(tmp, path.stat().st_mode & 0o777 if path.exists() else 0o644)
    os.replace(tmp, path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--digests", type=Path, required=True)
    ap.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    ap.add_argument("--nextflow-config", type=Path, default=None)
    ap.add_argument("--images", type=Path, default=None, help="default: next to --config")
    ap.add_argument("--root", type=Path, default=REPO, help="the tree the image keys come from")
    ap.add_argument("--date", required=True, help="build date, YYYY-MM-DD (UTC)")
    ap.add_argument("--source-commit", required=True, help="the commit the images were built from")
    ap.add_argument(
        "--provenance",
        type=Path,
        default=None,
        help="JSON {tool: {commit, date}} for images built in an earlier run (existing tags)",
    )
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
        prov = read_provenance(args.provenance) if args.provenance is not None else {}
        new_cfg = update_containers(
            args.config.read_text(), pins, args.date, args.source_commit, prov
        )
        new_nf = None
        nf: Path | None = None
        if args.nextflow_config is not None:
            nf = selector_file(args.nextflow_config)
            new_nf = update_conda_off(nf.read_text(), [p for p in pins if p in PROCESS_TOOL])
        images: Path = args.images or args.config.parent / "images.json"
        old_images = images.read_text() if images.is_file() else None
        new_images = update_images(old_images, pins, new_cfg, args.root, args.source_commit, prov)
    except (DigestError, OSError, ValueError) as exc:
        print(f"update_digests: {exc}", file=sys.stderr)
        return 1
    if new_cfg != args.config.read_text():
        _write(args.config, new_cfg)
    if new_nf is not None and nf is not None and new_nf != nf.read_text():
        _write(nf, new_nf)
    if new_images != old_images:
        _write(images, new_images)
    for process, ref in pins.items():
        print(f"pinned {process} {ref}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
