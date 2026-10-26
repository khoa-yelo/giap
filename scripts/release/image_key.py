"""Image keys and tags for the GHCR tools (M11-14 item 2; ADR-0060 §1).

Usage:

    python scripts/release/image_key.py [--root DIR] [--images PATH] TOOL ...
    python scripts/release/image_key.py [--root DIR] [--images PATH] --plan
    python scripts/release/image_key.py [--root DIR] [--images PATH] --check

``TOOL ...`` prints ``<tool>\\t<key>\\t<tag>`` per tool. ``--plan`` prints the GHCR tools whose key
differs from ``images.json`` (``will rebuild: …``). ``--check`` exits 1 with
``stale image: <tool> (inputs changed since <key>)`` for every ``images.json`` entry whose inputs
changed. The key itself lives in ``giap.io.image_key``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from giap.io.image_key import key, tag  # noqa: E402

UNPUBLISHED = {"uniop"}  # unlicensed, never published (R5)


def ghcr_tools(root: Path) -> list[str]:
    base = root / "workflow" / "containers"
    return sorted(p.name for p in base.iterdir() if p.is_dir() and p.name not in UNPUBLISHED)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--root", type=Path, default=REPO)
    ap.add_argument("--images", type=Path, default=None)
    ap.add_argument("tools", nargs="*", default=[])
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)
    if sum([bool(args.tools), args.plan, args.check]) != 1:
        ap.error("give TOOL ..., --plan or --check (exactly one)")
    root: Path = args.root
    images_path: Path = args.images or root / "workflow" / "conf" / "images.json"
    if args.tools:
        from giap.annotate.tools import TOOLS

        for tool in args.tools:
            if tool not in TOOLS:
                print(f"image_key: unknown tool {tool}", file=sys.stderr)
                return 2
            k = key(tool, root)
            print(f"{tool}\t{k}\t{tag(TOOLS[tool].version, k)}")
        return 0
    images = json.loads(images_path.read_text()) if images_path.is_file() else {}
    if args.plan:
        lines = []
        for tool in ghcr_tools(root):
            old = images.get(tool, {}).get("key")
            new = key(tool, root)
            if old != new:
                lines.append(f"will rebuild: {tool} ({old or 'unpinned'} → {new})")
        print("\n".join(lines) if lines else "will rebuild: none")
        return 0
    stale = [t for t, e in sorted(images.items()) if key(t, root) != e["key"]]
    for tool in stale:
        print(f"stale image: {tool} (inputs changed since {images[tool]['key']})")
    return 1 if stale else 0


if __name__ == "__main__":
    sys.exit(main())
