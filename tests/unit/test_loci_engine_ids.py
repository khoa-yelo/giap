"""M3-17 lead decision 10 (docs/specs/M3-17-mapping-steps.md, "Lead decisions after code review"):
``annotate_peptide_hits`` takes ``protein_ids`` as contract C9 gives it — a list (Parquet hands back
a numpy array) — as well as the ``;``-joined string, like ``contaminants._as_ids``. The
``engine_protein_mismatch`` rules of M3-14 (tests/unit/test_loci.py) are unchanged whatever the
container. Fixtures come from tests/unit/test_loci.py (hit at [130, 160) + lies in G1).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
import pandas as pd
import pytest
from test_loci import annotate, flags_of, hit_row, peps  # type: ignore[import-not-found]

KEY = ("MISM", 130, 160, "+")

CONTAINERS: dict[str, Callable[[list[str]], Any]] = {
    "list": list,
    "tuple": tuple,
    "ndarray": lambda ids: np.asarray(ids, dtype=object),
    "str_ndarray": lambda ids: np.asarray(ids, dtype=str),
}


def _peptides(ids: Any, search_db: str = "same") -> pd.DataFrame:
    rows = [hit_row("MISM", 130, 160, "+")]
    df = peps(rows, search_db=search_db)
    df["protein_ids"] = pd.Series([ids], index=df.index, dtype=object)
    return df


@pytest.mark.parametrize("container", sorted(CONTAINERS))
@pytest.mark.parametrize(
    ("ids", "flagged"),
    [
        (["GSIM_00005"], False),
        (["GSIM_00010", "GSIM_00005"], False),
        (["rev_GSIM_00010", "GSIM_00005"], False),
        (["rev_GSIM_00005"], True),  # a decoy entry never vouches for the gene
        (["CON__GSIM_00005", "GSIM_00010"], True),
        (["GSIM_00010"], True),
        ([], True),
    ],
)
def test_engine_ids_accept_list_like_protein_ids(
    container: str, ids: list[str], flagged: bool
) -> None:
    rows = [hit_row("MISM", 130, 160, "+")]
    out = annotate(rows, peptides=_peptides(CONTAINERS[container](ids)))
    assert ("engine_protein_mismatch" in flags_of(out[KEY]["flags"])) == flagged


@pytest.mark.parametrize("container", sorted(CONTAINERS))
def test_list_and_joined_string_agree(container: str) -> None:
    rows = [hit_row("MISM", 130, 160, "+")]
    for ids in (["GSIM_00010", "GSIM_00005"], ["GSIM_00010"]):
        a = annotate(rows, peptides=_peptides(CONTAINERS[container](ids)))
        b = annotate(rows, peptides=_peptides(";".join(ids)))
        assert a[KEY]["flags"] == b[KEY]["flags"], ids


def test_list_ids_without_same_strain_never_flag() -> None:
    rows = [hit_row("MISM", 130, 160, "+")]
    out = annotate(rows, peptides=_peptides(np.asarray(["GSIM_00010"]), search_db="other"))
    assert "engine_protein_mismatch" not in flags_of(out[KEY]["flags"])
