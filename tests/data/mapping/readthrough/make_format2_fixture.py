"""Regenerate sixframe_format2.npz (run on a checkout whose sixframe_index writes format 2).

This is the M7-16 work file of ``helpers.readthrough_index.world()`` (with genes, ``seed=7``)
written on main at 1f7203d3, before M11-16b. It pins table-11 byte identity (every array and the
meta apart from ``format`` and ``translation_table``) and the D3 format-2 load rule.

PYTHONPATH=src:tests/sim:tests python tests/data/mapping/readthrough/make_format2_fixture.py
"""

import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
from helpers.readthrough_index import run_sixframe_index

FIXTURE2 = Path(__file__).resolve().parent / "sixframe_format2.npz"

with tempfile.TemporaryDirectory() as d:
    out = run_sixframe_index(Path(d), with_genes=True)
    with np.load(out) as npz:
        assert json.loads(str(npz["meta"]))["format"] == 2, "not a format-2 writer"
    shutil.copyfile(out, FIXTURE2)
print(FIXTURE2)
