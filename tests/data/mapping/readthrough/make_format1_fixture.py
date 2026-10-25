"""Regenerate sixframe_format1.npz (run on a checkout whose sixframe_index writes format 1).

PYTHONPATH=src:tests/sim:tests python tests/data/mapping/readthrough/make_format1_fixture.py
"""

import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
from helpers.readthrough_index import FIXTURE, run_sixframe_index

with tempfile.TemporaryDirectory() as d:
    out = run_sixframe_index(Path(d), with_genes=False)
    with np.load(out) as npz:
        assert json.loads(str(npz["meta"]))["format"] == 1, "not a format-1 writer"
    shutil.copyfile(out, FIXTURE)
print(FIXTURE)
