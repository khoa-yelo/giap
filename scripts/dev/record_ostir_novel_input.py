"""Write the rule-9 OSTIR input of M7-08: three sim100k novel candidates, one multi-start.

Usage: python scripts/dev/record_ostir_novel_input.py <outdir>
(PYTHONPATH=src:tests/sim:tests). Simulates gippasim SIM100K with the M7-08 test configuration
(tests/helpers/novel_steps_fixture.py: M3-19 proteomics, planted novel ORFs and N-terminal
extensions), runs ``gippa run --executor local`` on it and keeps three rows of the
``ostir_input_novel.csv`` that ``novel_candidates`` wrote: the first multi-start row (``end >
start``) and the first two other rows, in file order. Writes ``<outdir>/ostir_input.csv``.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path


def main(argv: list[str]) -> int:
    from helpers import novel_steps_fixture as nf

    out = Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        sim = nf.simulate(root / "sim")
        with nf.db_env(root):
            code, log = nf.run_gippa(sim, root / "OUT")
        if code != 0:
            print(log, file=sys.stderr)
            return code
        lines = (root / "OUT" / "work" / "ostir_input_novel").read_text().splitlines()
    header, rows = lines[0], lines[1:]
    multi = next(r for r in rows if int(r.split(",")[3]) > int(r.split(",")[2]))
    keep = [multi, *[r for r in rows if r != multi][:2]]
    keep.sort(key=rows.index)
    (out / "ostir_input.csv").write_text("\n".join([header, *keep]) + "\n", encoding="utf-8")
    print(f"{len(rows)} candidate rows; kept {[r.split(',')[0] for r in keep]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
