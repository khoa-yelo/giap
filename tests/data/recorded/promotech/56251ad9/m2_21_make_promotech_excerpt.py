"""Regenerate the K-12 Promotech excerpt (M2-21; run by hand with ``$GIAP_ROOT`` set).

Source: ``$GIAP_ROOT/recordings/promotech/k12_mg1655/`` (GIAP-written ``w<index>`` windows, targeted
``-s -m RF-HOT`` runs). Keeps the first ``N`` records of ``windows_000.fa`` and the last ``N`` of
``windows_009.fa`` with their prediction rows copied byte-for-byte (header line kept), so the
excerpt covers the lowest and the highest window indices of the run.
"""

from __future__ import annotations

import os
from pathlib import Path

N = 100
SRC = Path(os.environ["GIAP_ROOT"]) / "recordings" / "promotech" / "k12_mg1655"
OUT = Path(__file__).resolve().parent / "k12_mg1655_excerpt"


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for chunk, head in (("000", True), ("009", False)):
        fa = (SRC / "windows" / f"windows_{chunk}.fa").read_text().splitlines(keepends=True)
        pred = (SRC / f"promotech_{chunk}" / "sequences_predictions.csv").read_text()
        plines = pred.splitlines(keepends=True)
        assert len(fa) == 2 * (len(plines) - 1)
        fa_keep = fa[: 2 * N] if head else fa[-2 * N :]
        p_keep = plines[1 : N + 1] if head else plines[-N:]
        (OUT / f"windows_{chunk}.fa").write_text("".join(fa_keep))
        (OUT / f"windows_{chunk}.sequences_predictions.csv").write_text(plines[0] + "".join(p_keep))


if __name__ == "__main__":
    main()
