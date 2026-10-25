"""The mapping steps with annotated stops (M7-16a rule 8; spec Q10 index format).

- ``sixframe_index`` gains the input ``genes`` (M2-04 ``genes.parquet``), writes index format 2
  with the target and per-decoy-type cells, and emits one INFO event ``mapping.stop_sites`` with
  ``data = {"terminal", "recoded", "not_a_stop", "cells_target", "cells_decoy": {kind: n}}``;
- an empty gene model gives zero sites, the same frames and the event with zeros;
- ``map_tolerant`` loads the cells and passes them to both routes (a read-through peptide gets its
  cost-1 target row); a format-1 index makes it raise ("rebuild the six-frame index").

World: ``helpers.readthrough_index.world()`` (6 terminal sites, 7 target cells: the stop at
``b:[3, 6) -`` has an origin copy; ``A7`` ends in GCT, so ``not_a_stop`` = 1).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from helpers.readthrough_index import FIXTURE, SEED, inputs, world
from helpers.readthrough_toys import model, write_annotation

from giap.config import resolve_config
from giap.mapping.sixframe import build_sixframe, decode
from giap.steps import load_step_modules
from giap.steps import registry as reg

EVENT = "mapping.stop_sites"


@pytest.fixture(autouse=True)
def _steps() -> None:
    load_step_modules()


def _cfg():  # type: ignore[no-untyped-def]
    return resolve_config(overrides=[f"seed={SEED}"])


def _six(tmp: Path, ins: dict[str, Path]) -> Path:
    out = tmp / "sixframe.npz"
    spec = reg.get("sixframe_index")
    reg.run_step("sixframe_index", _cfg(), {k: ins[k] for k in spec.inputs}, {"sixframe": out})
    return out


def _events(caplog: pytest.LogCaptureFixture) -> list[logging.LogRecord]:
    return [r for r in caplog.records if getattr(r, "event", None) == EVENT]


def test_sixframe_index_declares_the_genes_input() -> None:
    spec = reg.get("sixframe_index")
    assert {"genome", "replicons", "peptides", "genes"} <= set(spec.inputs)
    assert spec.outputs == ("sixframe",)


def test_sixframe_index_writes_format_2_cells_and_the_event(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    from giap.mapping.sixframe import load_stop_cells

    caplog.set_level(logging.INFO)
    out = _six(tmp_path, inputs(tmp_path))
    with np.load(out) as npz:
        meta = json.loads(str(npz["meta"]))
    assert meta["format"] == 2 and len(meta["stop_sites"]) == 6
    cells = load_stop_cells(out)
    assert len(cells.sites) == 6
    assert sum(len(c) for c in cells.target) == 7
    assert set(cells.decoys) == {"reverse", "shuffle"}
    ev = _events(caplog)
    assert len(ev) == 1 and ev[0].levelno == logging.INFO
    assert ev[0].data == {  # type: ignore[attr-defined]
        "terminal": 6, "recoded": 0, "not_a_stop": 1, "cells_target": 7,
        "cells_decoy": {"reverse": 7, "shuffle": 7},
    }  # fmt: skip


def test_empty_gene_model_gives_zero_sites_and_the_same_frames(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    from giap.mapping.sixframe import load_stop_cells

    caplog.set_level(logging.INFO)
    ins = inputs(tmp_path)
    _genome, reps, _gm = world()
    ins["genes"] = write_annotation(model(reps, []), tmp_path / "empty")["genes"]
    out = _six(tmp_path, ins)
    cells = load_stop_cells(out)
    assert cells.sites == () and all(len(c) == 0 for c in cells.target)
    assert all(len(c) == 0 for d in cells.decoys.values() for c in d)
    assert _events(caplog)[0].data == {  # type: ignore[attr-defined]
        "terminal": 0, "recoded": 0, "not_a_stop": 0, "cells_target": 0,
        "cells_decoy": {"reverse": 0, "shuffle": 0},
    }  # fmt: skip
    with np.load(out) as new, np.load(FIXTURE) as old:
        for key in ("target_aa", "reverse_aa", "shuffle_aa"):
            assert new[key].tobytes() == old[key].tobytes(), key


def _rt_peptides(tmp: Path) -> tuple[Path, str, tuple[int, int]]:
    """A database peptide over A1's TAA (``a:[297, 300) +``, frame +0, p = 99) with the
    stop read as W, plus one exact interior peptide; returns (path, peptide, genomic window)."""
    genome, reps, _ = world()
    index = build_sixframe(genome, reps)
    f0 = decode(index.frames[0].aa)
    assert f0[99] == "*"
    a = 92
    assert f0[a : a + 15].count("*") == 1 and "X" not in f0[a : a + 15]
    pep = f0[a:99] + "W" + f0[100 : a + 15]
    inner = next(f0[s : s + 12] for s in range(0, 80) if "*" not in f0[s : s + 12])
    df = pd.DataFrame({
        "peptide": [pep, inner], "source": ["database"] * 2, "q_value": [0.001] * 2,
        "pep": [0.01] * 2, "score": [np.nan] * 2, "is_decoy": [False] * 2,
        "search_db": ["other"] * 2,
    })  # fmt: skip
    path = tmp / "peptides_rt.parquet"
    df.to_parquet(path, index=False)
    return path, pep, (3 * a, 3 * (a + 15))


def _map(tmp: Path, six: Path, peptides: Path) -> pd.DataFrame:
    cfg = _cfg()
    hx, ht = tmp / "hits_exact.parquet", tmp / "hits_tolerant.parquet"
    reg.run_step("map_exact", cfg, {"peptides": peptides, "sixframe": six}, {"hits_exact": hx})
    reg.run_step("map_tolerant", cfg, {"peptides": peptides, "sixframe": six, "hits_exact": hx},
                 {"hits_tolerant": ht})  # fmt: skip
    return pd.read_parquet(ht)


def test_map_tolerant_passes_the_cells(tmp_path: Path) -> None:
    ins = inputs(tmp_path)
    ins["peptides"], pep, (s, e) = _rt_peptides(tmp_path)
    six = _six(tmp_path, ins)
    hits = _map(tmp_path, six, ins["peptides"])
    t = hits[(hits["peptide"] == pep) & (hits["search"] == "target")]
    assert len(t) == 1
    r = t.iloc[0]
    assert (r["contig"], r["start"], r["end"], r["strand"]) == ("a", s, e, "+")
    assert (r["n_diff"], r["diff_types"], r["route"]) == (1, "substitution", "tolerant")


def test_map_tolerant_rejects_a_format_1_index(tmp_path: Path) -> None:
    ins = inputs(tmp_path)
    ins["peptides"], _pep, _ = _rt_peptides(tmp_path)
    with pytest.raises(reg.StepError, match="rebuild the six-frame index"):
        _map(tmp_path, FIXTURE, ins["peptides"])
