"""Promoter Calculator command line (GIAP_SPEC.md §7.1.1; M2-14).

Promoter Calculator is GPL: it is never imported, only run as a separate process from its own
environment, once per TU upstream region with ``-j 1`` (regions run in parallel, G11). Kept in sync
with ``workflow/modules/local/promoter_calculator/main.nf``.
"""

from __future__ import annotations

from pathlib import Path


def promoter_calculator_command(region_fa: Path, out_csv: Path) -> list[str]:
    """``promoter-calculator -i <region.fa> -o <region.csv> -j 1``."""
    return ["promoter-calculator", "-i", str(region_fa), "-o", str(out_csv), "-j", "1"]
