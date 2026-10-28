# tests/data/m11_16c/table11_digests.json

Golden sha256 digests of gippasim outputs at genetic code 11, used by
`tests/unit/test_m11_16c_identity.py` (M11-16 item 8, AC 14: code 11 stays byte-identical).

- Source: the simulator in `tests/sim/gippasim/` and `tests/helpers/genetic_code_fixture.py` at
  main 1f7203d3 (branch point of `spec/translation-tables`; unchanged on main through 8c619d80),
  before any M11-16c change. Generated 2026-10-28 by the 16c test-author.
- Contents: `simulate` (every file `gippasim.simulate` writes, per case and seed 1-3, output
  directory replaced by `<OUT>`), `layout` (each `Layout` frame as CSV), `m11_14w_fixture`
  (`synthetic_genome(seed, gc, code)` at codes 4 and 11). Checked stable under different
  `PYTHONHASHSEED` values.
- Never regenerate in 16c: a mismatch is a regression. Regenerating later needs a lead decision
  that the change in bytes is intended.

Generator (run with `PYTHONPATH=src:tests/sim:tests` from the repository root, pre-16c code):

```python
import json, sys, tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0, "tests/unit")
import test_m11_16c_identity as m
from gippasim.simulate import simulate
from gippasim.layout import LayoutConfig, simulate_layout
out = {"simulate": {}, "layout": {}, "m11_14w_fixture": {}}
for name, mk in m._cases().items():
    out["simulate"][name] = {}
    for s in m.SEEDS:
        d = Path(tempfile.mkdtemp())
        simulate(mk(), s, d)
        out["simulate"][name][str(s)] = m._digests(d)
lc = {"default": LayoutConfig(n_genes_chromosome=300),
      "draft_origin": LayoutConfig(n_genes_chromosome=300, plasmids=(20,), linear_contigs=(10, 10), origin_crossing=True, oric_frac=0.8)}
for name, cfg in lc.items():
    out["layout"][name] = {str(s): m._layout_digests(simulate_layout(cfg, np.random.default_rng(s))) for s in m.SEEDS}
for s in (1, 7):
    for gc in (0.28, 0.5):
        for c in (4, 11):
            out["m11_14w_fixture"][f"{s}:{gc}:{c}"] = m._fixture_digest(s, gc, c)
out["_provenance"] = "gippasim at main 1f7203d3 (pre-M11-16c); generated 2026-10-28 by the 16c test-author; never regenerate in 16c"
p = Path("tests/data/m11_16c/table11_digests.json")
p.parent.mkdir(parents=True, exist_ok=True)
p.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
print(p, p.stat().st_size)
```

The generator ran with `mk()` (no `genetic_code` key, which did not exist yet); the test passes
`genetic_code=11` explicitly.
