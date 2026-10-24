"""Generate the report's methods snippets from ``docs/methods`` (M9-16 behaviour 2).

Every ``<topic>.md`` holding a ``<!-- report:begin -->`` ... ``<!-- report:end -->`` block becomes
``<topic>.html``: one ``<p>`` per blank-line-separated paragraph (lines joined by a space), one
``<ul>`` of ``<li>`` per run of ``- `` lines (indented lines continue the item), text
HTML-escaped. ``--check`` writes nothing and exits 1 naming each topic whose snippet differs from
its source, is missing, or has no source (an orphan). Standard library only.

    python scripts/dev/sync_report_snippets.py [--check] [--docs DIR] [--out DIR]
"""

from __future__ import annotations

import argparse
import sys
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BEGIN, END = "<!-- report:begin -->", "<!-- report:end -->"


def block(text: str) -> str | None:
    """The marked block of a Markdown document, or ``None``."""
    start = text.find(BEGIN)
    if start < 0:
        return None
    stop = text.find(END, start)
    return text[start + len(BEGIN) : stop if stop >= 0 else len(text)]


def to_html(src: str) -> str:
    out: list[str] = []
    para: list[str] = []
    items: list[list[str]] = []

    def flush() -> None:
        if para:
            out.append(f"<p>{escape(' '.join(para), quote=False)}</p>")
            para.clear()
        if items:
            lis = "".join(f"<li>{escape(' '.join(i), quote=False)}</li>" for i in items)
            out.append(f"<ul>{lis}</ul>")
            items.clear()

    for raw in src.splitlines():
        line = raw.strip()
        if not line:
            flush()
        elif line.startswith("- "):
            if para:
                flush()
            items.append([line[2:].strip()])
        elif items and raw[:1].isspace():
            items[-1].append(line)
        else:
            if items:
                flush()
            para.append(line)
    flush()
    return "\n".join(out) + "\n"


def expected(docs: Path) -> dict[str, str]:
    found = {}
    for md in sorted(docs.glob("*.md")):
        body = block(md.read_text(encoding="utf-8"))
        if body is not None:
            found[md.stem] = to_html(body)
    return found


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--docs", type=Path, default=ROOT / "docs" / "methods")
    ap.add_argument("--out", type=Path, default=ROOT / "src" / "giap" / "report" / "snippets")
    args = ap.parse_args(argv)
    want = expected(args.docs)
    have = {p.stem: p for p in sorted(args.out.glob("*.html"))} if args.out.is_dir() else {}
    if args.check:
        bad = [t for t, h in want.items() if t not in have or have[t].read_text("utf-8") != h]
        bad += [f"{t} (orphan: no report block in {args.docs.name}/{t}.md)"
                for t in have if t not in want]  # fmt: skip
        for topic in bad:
            print(f"snippet out of sync: {topic}")
        return 1 if bad else 0
    args.out.mkdir(parents=True, exist_ok=True)
    for topic, text in want.items():
        (args.out / f"{topic}.html").write_text(text, encoding="utf-8")
    for topic, path in have.items():
        if topic not in want:
            path.unlink()
    return 0


if __name__ == "__main__":
    sys.exit(main())
