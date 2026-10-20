"""Regenerate the K-12 UniOP excerpt (M2-21; run by hand with ``$GIAP_ROOT`` set, not by pytest).

Source: ``$GIAP_ROOT/recordings/uniop/k12_mg1655/`` (one contig, U00096.3, genes in start order).
Keeps the first ``N`` genes of ``uniop_cds.gff`` / ``uniop_map.tsv`` and ``uniop.pred`` rows
``1 .. N-1`` (title and header kept), all lines byte-for-byte: row ``k`` pairs gene ``k`` with gene
``k+1``, so every kept row has both partners in the excerpt.
"""

from __future__ import annotations

import os
from pathlib import Path

N = 300
SRC = Path(os.environ["GIAP_ROOT"]) / "recordings" / "uniop" / "k12_mg1655"
OUT = Path(__file__).resolve().parent / "k12"


def main() -> None:
    OUT.mkdir(exist_ok=True)
    gff = (SRC / "uniop_cds.gff").read_text().splitlines(keepends=True)
    mp = (SRC / "uniop_map.tsv").read_text().splitlines(keepends=True)
    pred = (SRC / "uniop.pred").read_text().splitlines(keepends=True)
    assert gff[0].startswith("##gff-version")
    starts = [int(ln.split("\t")[3]) for ln in gff[1:]]
    assert starts == sorted(starts) and len({ln.split("\t")[0] for ln in gff[1:]}) == 1
    (OUT / "uniop_cds.gff").write_text("".join(gff[: N + 1]))
    (OUT / "uniop_map.tsv").write_text("".join(mp[: N + 1]))
    (OUT / "uniop.pred").write_text("".join(pred[: 2 + N - 1]))


if __name__ == "__main__":
    main()
