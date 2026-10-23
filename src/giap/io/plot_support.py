"""Writer of the plot support tables under ``plots/support/`` (M9-03; ADR-0043).

``write_support`` writes the frame as given (already 1-based, ``giap.viz.support``) through
``giap.io.tsv`` (§9.2 cells). A TSV larger than ``GZIP_THRESHOLD`` bytes is replaced by
``<name>.tsv.gz`` with a fixed gzip header (``mtime=0``); exactly one of the two exists after a
write. ``giap report`` reads these files and never calls this writer.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from giap.io._open import open_text_write
from giap.io.tsv import write_tsv

__all__ = ["GZIP_THRESHOLD", "SUPPORT_DIR", "write_support"]

#: ``plots/`` subdirectory of the support tables
SUPPORT_DIR = "support"
#: ADR-0043 decision 2 / lead decision 6: 1 MB = 2^20 bytes of uncompressed TSV
GZIP_THRESHOLD = 2**20


def write_support(name: str, df: pd.DataFrame, plots_dir: Path) -> Path:
    """Write ``plots_dir/support/<name>.tsv`` (or ``.tsv.gz`` above the threshold); returns it."""
    base = plots_dir / SUPPORT_DIR
    base.mkdir(parents=True, exist_ok=True)
    plain, packed = base / f"{name}.tsv", base / f"{name}.tsv.gz"
    write_tsv(plain, [str(c) for c in df.columns], df.itertuples(index=False, name=None))
    if plain.stat().st_size <= GZIP_THRESHOLD:
        packed.unlink(missing_ok=True)
        return plain
    with open_text_write(packed) as fh:
        fh.write(plain.read_text(encoding="utf-8"))
    plain.unlink()
    return packed
