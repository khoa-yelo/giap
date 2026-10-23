#!/usr/bin/env python3
"""Fail if any file under src/ mentions ``Rscript`` (M8-12 licence boundary; GIAP_SPEC.md R5).

goseq and BiasedUrn live only in the conda env giap-golden-r; R fgsea runs in giap-fgsea, and
src/ starts R only through workflow/bin/giap_fgsea.R (ADR-0038). Every file is scanned, not only
``*.py``.  Usage: check_no_rscript.py [ROOT] (default src).
"""

from __future__ import annotations

import sys
from pathlib import Path

NEEDLE = b"Rscript"


def violations(root: Path) -> list[str]:
    found = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        for lineno, line in enumerate(path.read_bytes().split(b"\n"), 1):
            if NEEDLE in line:
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
