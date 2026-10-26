#!/usr/bin/env python3
"""Fail if any file under src/ mentions ``Rscript`` (M8-12 licence boundary; GIAP_SPEC.md R5).

goseq and BiasedUrn live only in the conda env giap-golden-r; R fgsea runs in giap-fgsea, and
src/ starts R only through workflow/bin/giap_fgsea.R (ADR-0038). Every file is scanned, not only
``*.py``. One exception (M11-14 14c): giap/io/fgsea_r.py may name ``Rscript`` on a line that also
names the baked-in giap_fgsea.R (``BAKED_SCRIPT`` or "baked-in script"), the image argv that runs
that same script in a separate process.  Usage: check_no_rscript.py [ROOT] (default src).
"""

from __future__ import annotations

import sys
from pathlib import Path

NEEDLE = b"Rscript"
IMAGE_LAUNCHER = ("giap", "io", "fgsea_r.py")
BAKED = (b"BAKED_SCRIPT", b"baked-in script")


def violations(root: Path) -> list[str]:
    found = []
    files = (p for p in root.rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    for path in sorted(files):  # bytecode compiled from a scanned .py repeats its strings
        for lineno, line in enumerate(path.read_bytes().split(b"\n"), 1):
            if NEEDLE in line and not (
                path.parts[-3:] == IMAGE_LAUNCHER and any(b in line for b in BAKED)
            ):
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
