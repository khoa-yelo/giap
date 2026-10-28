"""``gippa run --executor local`` stages the five function resources (M8-06 item 7, 06b).

docs/specs/M8-06-function-step-and-run.md, item 7 "Resource staging": ``cli._local_roots`` stages
each of ``gippa.setupdb.installed.FUNCTION_RESOURCES`` that ``installed.resource`` resolves, under
its setup item name, and only when ``"function" in cfg.modules`` (``--modules core`` copies no
32 MB OBO). The ``resource()`` rule (M4-15 decision 5): the setup manifest lists the item *and*
the file exists. The COG 2020 definition table (``cog2020_def``, M4-15) is staged as before,
whatever the modules.

The database is the synthetic one of ``tests/helpers/function_db.py`` (never the real
``$GIPPA_DB_DIR``), with the manifest rewritten in the per-name form of
``helpers.mapping_steps.write_test_db``, so preflight's run-mode environment check passes (every
requirement listed) while the function items are listed only as each test says. Preflight runs
for real on a 6-kb random genome with three proteins.
"""

from __future__ import annotations

import json
import random
import tempfile
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import pytest
from helpers.function_db import FUNCTION_RESOURCES, ITEMS, PATHS, FunctionDb, write_function_db
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

XREFS = pd.DataFrame({"locus_tag": ["g1", "g1", "g1", "g2"], "db": ["COG", "GO", "KEGG", "COG"],
                      "id": ["COG0001", "GO:0000001", "K00001", "COG0002"]})  # fmt: skip
FN = ["core", "function"]


def _db(
    root: Path, items: Iterable[str] = ITEMS, listed: Iterable[str] | None = None
) -> FunctionDb:
    """The synthetic function DB (``items`` written, ``listed`` in the manifest, default
    ``items``) plus every preflight requirement, in one per-name manifest."""
    from helpers.mapping_steps import write_test_db

    items = list(items)
    names = items if listed is None else list(listed)
    db = write_function_db(root, XREFS, items=items, listed=names)
    write_test_db(root)  # rewrites db_manifest.json: every requirement of every module
    doc = json.loads((root / "db_manifest.json").read_text(encoding="utf-8"))
    for name in ITEMS:
        doc.pop(name, None)
    doc |= {n: {"version": "synthetic", "path": PATHS[n]} for n in names}
    (root / "db_manifest.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return db


def _evidence(root: Path) -> tuple[Path, Path]:
    rng = random.Random(3)
    root.mkdir(parents=True, exist_ok=True)
    genome = root / "genome.fa"
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    genome.write_text(">chr1 [topology=circular]\n" + seq + "\n", encoding="utf-8")
    aa = "ACDEFGHIKLMNPQRSTVWY"
    prot = root / "proteins.faa"
    prot.write_text("".join(f">p{i}\nM" + "".join(rng.choice(aa) for _ in range(150)) + "\n"
                            for i in range(3)), encoding="utf-8")  # fmt: skip
    return genome, prot


@pytest.fixture
def evidence(tmp_path: Path) -> tuple[Path, Path]:
    return _evidence(tmp_path / "in")


@pytest.fixture(autouse=True)
def _no_real_db(monkeypatch: pytest.MonkeyPatch) -> None:
    """``$GIPPA_DB_DIR`` never points at the real databases."""
    monkeypatch.delenv("GIPPA_DB_DIR", raising=False)


def _roots(
    out: Path,
    genome: Path,
    prot: Path,
    db: Path | None,
    modules: list[str],
    *,
    env_db: Path | None = None,
    db_option: bool = True,
) -> dict[str, Path]:
    """``_local_roots`` with ``--db-dir db`` (or none) and ``$GIPPA_DB_DIR`` = ``env_db or db``."""
    from gippa.cli import _local_roots
    from gippa.config import GippaConfig

    inputs: dict[str, Path | None] = {"proteins": prot, "peptides": None, "abundance": None,
                                      "annotation": None, "reference_peptides": None,
                                      "db_dir": db if db_option else None}  # fmt: skip
    with pytest.MonkeyPatch.context() as mp:
        env = env_db or db
        if env is not None:
            mp.setenv("GIPPA_DB_DIR", str(env))
        return _local_roots(GippaConfig(modules=modules), out, genome, inputs)


def _same(a: Path, b: Path) -> bool:
    return Path(a).resolve() == Path(b).resolve()


def test_function_on_stages_the_five_resources_under_their_item_names(
    tmp_path: Path, evidence: tuple[Path, Path]
) -> None:
    db = _db(tmp_path / "db")
    roots = _roots(tmp_path / "OUT", *evidence, db.root, FN)
    for name in FUNCTION_RESOURCES:
        assert name in roots, f"{name} not staged: {sorted(roots)}"
        assert _same(roots[name], db.root / PATHS[name]), name
    assert _same(roots["cog2020_def"], db.root / PATHS["cog2020_def"])


def test_function_on_reads_gippa_db_dir_when_no_db_dir_option(
    tmp_path: Path, evidence: tuple[Path, Path]
) -> None:
    db = _db(tmp_path / "db")
    roots = _roots(tmp_path / "OUT", *evidence, db.root, FN, db_option=False)
    assert set(FUNCTION_RESOURCES) <= set(roots), sorted(roots)


def test_modules_core_stages_no_function_resource(
    tmp_path: Path, evidence: tuple[Path, Path]
) -> None:
    """``--modules core`` copies no OBO; the COG definition table is still staged (M4-15)."""
    db = _db(tmp_path / "db")
    roots = _roots(tmp_path / "OUT", *evidence, db.root, ["core"])
    assert not set(FUNCTION_RESOURCES) & set(roots), sorted(roots)
    assert "cog2020_def" in roots


def test_other_optional_modules_without_function_stage_nothing(
    tmp_path: Path, evidence: tuple[Path, Path]
) -> None:
    db = _db(tmp_path / "db")
    roots = _roots(tmp_path / "OUT", *evidence, db.root, ["core", "mobile", "bgc"])
    assert not set(FUNCTION_RESOURCES) & set(roots), sorted(roots)


def test_only_the_installed_subset_is_staged(tmp_path: Path, evidence: tuple[Path, Path]) -> None:
    """Spec "only some resources (e.g. COG only)": only ``cog2020_fun`` (and ``cog2020_def``)."""
    db = _db(tmp_path / "db", items=("cog2020_fun", "cog2020_def"))
    roots = _roots(tmp_path / "OUT", *evidence, db.root, FN)
    assert [k for k in FUNCTION_RESOURCES if k in roots] == ["cog2020_fun"], sorted(roots)


def test_listed_but_deleted_file_is_not_staged(tmp_path: Path, evidence: tuple[Path, Path]) -> None:
    """Spec "a resource listed in the manifest but its file deleted": not installed."""
    db = _db(tmp_path / "db", items=[n for n in ITEMS if n != "go_basic"], listed=ITEMS)
    roots = _roots(tmp_path / "OUT", *evidence, db.root, FN)
    assert "go_basic" not in roots
    assert {"goslim_generic", "kegg_ko_pathway", "kegg_pathway_names", "cog2020_fun"} <= set(roots)


def test_file_present_but_not_listed_is_not_staged(
    tmp_path: Path, evidence: tuple[Path, Path]
) -> None:
    """The ``resource()`` rule: a file that setup never recorded is not installed."""
    db = _db(tmp_path / "db", listed=[n for n in ITEMS if n != "kegg_ko_pathway"])
    roots = _roots(tmp_path / "OUT", *evidence, db.root, FN)
    assert "kegg_ko_pathway" not in roots
    assert "go_basic" in roots


def test_no_database_at_all_stages_nothing_and_does_not_fail(
    tmp_path: Path, evidence: tuple[Path, Path]
) -> None:
    roots = _roots(tmp_path / "OUT", *evidence, None, FN)
    assert not set(FUNCTION_RESOURCES) & set(roots), sorted(roots)
    assert {"genome", "preflight"} <= set(roots)


@settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    written=st.sets(st.sampled_from(ITEMS)),
    listed=st.sets(st.sampled_from(ITEMS)),
    function_on=st.booleans(),
)
def test_staged_resources_are_written_and_listed_iff_function_on(
    written: set[str], listed: set[str], function_on: bool
) -> None:
    """Invariant (item 7): the staged function keys are written ∩ listed when ``function`` is on,
    else none; each staged path is the item's ``db_sources.yaml`` path under the database."""
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        genome, prot = _evidence(root / "in")
        db = _db(root / "db", items=written, listed=listed)
        roots = _roots(root / "OUT", genome, prot, db.root, FN if function_on else ["core"])
        got = {k for k in FUNCTION_RESOURCES if k in roots}
        want = (written & listed & set(FUNCTION_RESOURCES)) if function_on else set()
        assert got == want
        for k in got:
            assert _same(roots[k], db.root / PATHS[k])
