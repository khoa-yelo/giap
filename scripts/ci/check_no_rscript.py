#!/usr/bin/env python3
"""Fail if any file under src/ mentions ``Rscript`` (M8-12 licence boundary; GIAP_SPEC.md R5).

goseq and BiasedUrn live only in the conda env giap-golden-r; R fgsea runs in giap-fgsea, and
src/ starts R only through workflow/bin/giap_fgsea.R (ADR-0038). Every file is scanned, not only
``*.py``. One exception (M11-14 14c): in giap/io/fgsea_r.py, the image argv list that starts
``["Rscript", BAKED_SCRIPT`` (optionally with ``"--vanilla"`` between) runs that same script in a
separate process; any other mention there, ``-e`` included, still fails.
Usage: check_no_rscript.py [ROOT] (default src).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

NEEDLE = b"Rscript"
IMAGE_LAUNCHER = ("giap", "io", "fgsea_r.py")
IMAGE_ARGV = re.compile(rb'\["Rscript", (?:"--vanilla", )?BAKED_SCRIPT[,\]]')


def violations(root: Path) -> list[str]:
    found = []
    files = (p for p in root.rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    for path in sorted(files):  # bytecode compiled from a scanned .py repeats its strings
        for lineno, line in enumerate(path.read_bytes().split(b"\n"), 1):
            allowed = (
                path.parts[-3:] == IMAGE_LAUNCHER
                and line.count(NEEDLE) == 1
                and IMAGE_ARGV.search(line) is not None
                and b'"-e"' not in line
            )
            if NEEDLE in line and not allowed:
                found.append(f"{path}:{lineno}: mentions Rscript")
    return found


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "src")
    bad = violations(root)
    for line in bad:
        print(line)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
