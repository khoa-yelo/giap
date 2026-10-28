"""M2-15b: a candidate at window 0 of a strand (the contig edge) is kept by ``chunk``.

#455 review N3: no fake sequence in ``test_m2_15b_tnet_shard_cli.py`` has a candidate at window 0,
so a ``j >= 0`` -> ``j > 0`` mutant of the candidate filter in ``gippa_terminatornet_shard``
(upstream ``scan_oneprocess_1`` returns ``j = -1`` for "no candidate"; 0 is a real window) survives.
The fake scorer makes a candidate wherever the window starts with ``G``, and a call (pred >= 0.5)
unless the next base is ``T``: a record that starts with ``GA`` and ends with ``TC`` has a called
candidate at window 0 of both strands.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from unit.test_m2_15b_tnet_shard_cli import (  # noqa: F401  (fixture home)
    _chunk,
    _fasta,
    _plan,
    _rc,
    _reference,
    _run,
    _seq,
    _workdir,
    home,
)


@pytest.mark.parametrize(("length", "core"), [(900, 300), (400, 250_000)])
def test_window_0_candidate_is_called_on_both_strands(
    home: Path,  # noqa: F811
    tmp_path: Path,
    length: int,
    core: int,
) -> None:
    s = "GA" + _seq(length - 4, length) + "TC"
    wd = _workdir(home, tmp_path, _fasta([("edge", s)]))
    assert _run(home, wd, "plan", "--core", str(core)).returncode == 0
    for strand, seq in (("+", s), ("-", _rc(s))):
        assert seq.startswith("GA")
        ref = _reference(seq)
        assert ref[0][0] == 0 + 3, "the reference must call window 0 (hairpin_start 3)"
        rows = [x for x in _plan(wd) if x["strand"] == strand]
        assert rows[0]["core_start"] == "0"
        concat: list[tuple] = []
        for x in rows:
            payload, _js = _chunk(home, wd, x["chunk_id"])
            concat += [tuple(c) for c in payload["calls"]]
        first = [c for c in concat if c[0] == 3]
        assert first, f"{strand}: the window-0 candidate was dropped"
        assert type(first[0][7]) is np.float64
        assert concat == ref, f"{strand}: sharded calls differ from the single-chunk scan"
