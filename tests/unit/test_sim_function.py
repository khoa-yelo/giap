"""gippasim function and determinants scenarios (issue M8-14; GIPPA_SPEC.md §7.7.5, §7.8.3, §8.4,
§8.5 attacks 2 and 3, §11.3, §11.4; docs/specs/M8-14-sim-function-determinants.md).

The spec's ``tests/sim/tests/test_function_sim.py`` lives in tests/unit like the other gippasim
tests (the simulator tree is the implementer's). Examples 1, 3 and 6's "over seeds" claims (and the
"detection depends only on p_true" half of example 2) are gated in
tests/calibration/test_sim_function_calibration.py; the ontology invariants are Hypothesis
properties in tests/property/test_sim_function_props.py.

Interface assumed (the spec gives only the scenario signatures ``-> Pair``; reported to the lead):

- ``gippasim.scenarios``: ``ora_null(n)``, ``ora_attack_confounded(n, member_shift)``,
  ``gsea_null(n)``, ``gsea_attack_operon_packed(n, n_sets)``,
  ``determinants_null(n, *, noise_sd=1.0, features="planted")`` and
  ``determinants_planted(n, coefs, noise_sd, *, features="planted")``, each returning
  ``Pair = (LayoutConfig, ExpressionConfig)``; ``coefs`` keys are ``promoter``, ``rbs``, ``cai``
  (blocks ``transcription``, ``translation``, ``codon``); ``features="gaussian"`` draws iid N(0, 1)
  features (unit tests only, lead decision 2).
- ``gippasim.function_sim``:

  * ``RESERVED_PREFIXES``: the reserved ID prefixes (``GO:9``, ``K9``, ``map9``; decision 1);
  * ``CONFOUND_COVARIATE``: the M4-13 covariate carrying ``ora_attack_confounded``'s shift,
    ``"n_obs"`` (column ``n_obs_true``), realised by K/R edits of the members' ``aa`` (lead
    decisions 4 and 6, B1);
  * ``synthetic_ontology(n_terms, depth, rng, *, genes=())`` -> ``(obo_text, table)``, ``table``
    a C3-shaped (``locus_tag``, ``db``, ``id``) gene -> GO term table over ``genes``;
  * ``draw(pair, rng)`` -> a world with ``layout`` (gippasim Layout; the synthetic GO/KEGG/COG
    xrefs in its C3 ``xrefs``), ``expression`` (gippasim Expression), ``obo`` (str),
    ``kegg_links`` (``ko``, ``pathway``: ``ko:K9xxxx``, ``path:map9xxxx``), ``sets`` (one row per
    member: ``set_id``, ``collection``, ``locus_tag``, ``confounded``), ``features`` (determinants
    scenarios: ``locus_tag``, ``promoter``, ``rbs``, ``cai``) and ``truth`` (a mapping of frames);
  * ``write_truth(world, outdir)`` -> ``{name: path}``, writing ``truth/determinants.tsv`` (one row
    per feature: ``feature``, ``block``, ``beta``, ``true_r2``, ``residual_var``) and
    ``truth/determinants_genes.tsv`` (``locus_tag``, ``dominant_block``) for the determinants
    scenarios; the same frames are ``world.truth["determinants"]`` / ``["determinants_genes"]``.

Readings: "true R² from the β and feature variances" (example 5) is, for independent unit-variance
z features, ``Σβ² / (Σβ² + residual_var)``; the per-gene dominant block is the block of
``argmax_j |β_j z_ij|`` with ``z`` the population (ddof 0) z-score over all genes; "TU sizes
identical to the base layout" (example 6) is read as the null keeping the planted scenario's
layout and detection pattern for the same seed.
"""

from __future__ import annotations

import math
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
GIPPASIM = ROOT / "tests" / "sim" / "gippasim"
sys.path.insert(0, str(ROOT / "tests" / "helpers"))
SCENARIOS = ("ora_null", "ora_attack_confounded", "gsea_null", "gsea_attack_operon_packed",
             "determinants_null", "determinants_planted")  # fmt: skip
EX5 = {"promoter": 0.5, "rbs": 0.4, "cai": 0.3}
BLOCKS = {"promoter": "transcription", "rbs": "translation", "cai": "codon"}
FEATURES = ("promoter", "rbs", "cai")
BAD_IDS = ["ora_n0", "gsea_neg", "shift_nan", "shift_inf", "n_sets0", "coef_key", "noise_neg",
           "coef_nan", "features_mode", "det_n0"]  # fmt: skip
AA_DERIVED = ("aa", "nt")  # layout.genes columns a K/R edit of a member may change
COG_LETTERS = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
GO_COLLECTIONS = ("go_bp", "go_mf", "go_cc")

RE_GO = re.compile(r"^GO:9\d{6}$")
RE_KO = re.compile(r"^K9\d{4}$")
RE_MAP = re.compile(r"^map9\d{4}$")
# any GO / KO / KEGG-map shaped token, reserved or not
ANY_GO = re.compile(r"GO:\d{7}")
ANY_KO = re.compile(r"(?<![A-Za-z0-9])K\d{5}(?!\d)")
ANY_MAP = re.compile(r"(?<![A-Za-z0-9])map\d{5}(?!\d)")


def _fs() -> Any:
    import gippasim.function_sim

    return gippasim.function_sim


def _sc() -> Any:
    import gippasim.scenarios

    return gippasim.scenarios


def draw(pair: Any, seed: int) -> Any:
    return _fs().draw(pair, np.random.default_rng(seed))


def _genes(world: Any) -> pd.DataFrame:
    return world.expression.genes.reset_index(drop=True)


def _tu_of(layout: Any) -> dict[str, str]:
    return {str(g): str(t) for t, gs in zip(layout.tus["tu_id"], layout.tus["genes"], strict=True)
            for g in gs}  # fmt: skip


def _members(world: Any) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for s, g in zip(world.sets["set_id"], world.sets["locus_tag"], strict=True):
        out.setdefault(str(s), set()).add(str(g))
    return out


def _confounded(world: Any) -> set[str]:
    s = world.sets
    return set(s.loc[s["confounded"].astype(bool), "locus_tag"].astype(str))


def _obo_ids(obo: str) -> list[str]:
    return re.findall(r"^id:\s*(\S+)", obo, flags=re.M)


def _read_go(obo: str, tmp: Path) -> Any:
    from gippa.io.function_resources import read_go

    path = tmp / "synthetic.obo"
    path.write_text(obo)
    return read_go(path)


def _z(x: Any) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    return (a - a.mean()) / a.std()


def _frames_equal(a: Any, b: Any) -> None:
    for name in ("genes", "tus", "xrefs", "truth"):
        pd.testing.assert_frame_equal(getattr(a.layout, name), getattr(b.layout, name))
    pd.testing.assert_frame_equal(a.expression.genes, b.expression.genes)


@pytest.fixture(scope="module")
def ora4000() -> Any:
    return draw(_sc().ora_null(4000), 1)


@pytest.fixture(scope="module")
def attack4000() -> Any:
    return draw(_sc().ora_attack_confounded(4000, -1.5), 1)


@pytest.fixture(scope="module")
def gsea4000() -> Any:
    return draw(_sc().gsea_null(4000), 2)


@pytest.fixture(scope="module")
def packed4000() -> Any:
    return draw(_sc().gsea_attack_operon_packed(4000, 20), 3)


@pytest.fixture(scope="module")
def planted3000() -> Any:
    return draw(_sc().determinants_planted(3000, EX5, 1.0, features="gaussian"), 4)


# ----------------------------------------------------------------------------- interface


def test_every_scenario_exists_and_returns_a_pair() -> None:
    from gippasim.expression import ExpressionConfig
    from gippasim.layout import LayoutConfig

    sc = _sc()
    pairs = [sc.ora_null(600), sc.ora_attack_confounded(600, -1.0), sc.gsea_null(600),
             sc.gsea_attack_operon_packed(600, 5), sc.determinants_null(600),
             sc.determinants_planted(600, EX5, 1.0)]  # fmt: skip
    for name in SCENARIOS:
        assert callable(getattr(sc, name)), name
        assert name in sc.__all__, name
    for p in pairs:
        assert isinstance(p, tuple) and len(p) == 2
        assert isinstance(p[0], LayoutConfig) and isinstance(p[1], ExpressionConfig)
        assert p[0].n_genes_chromosome == 600


def test_function_sim_exports() -> None:
    fs = _fs()
    for name in ("RESERVED_PREFIXES", "synthetic_ontology", "draw", "write_truth"):
        assert hasattr(fs, name), name
    prefixes = fs.RESERVED_PREFIXES
    values = set(prefixes.values()) if isinstance(prefixes, dict) else set(prefixes)
    assert values == {"GO:9", "K9", "map9"}


def test_function_sim_imports_nothing_from_gippa() -> None:
    """Sequence-free gippasim modules import nothing from ``gippa`` (as layout and expression)."""
    code = (
        "import sys, gippasim.function_sim, gippasim.scenarios; "
        "print(sorted(m for m in sys.modules if m == 'gippa' or m.startswith('gippa.')))"
    )
    env_path = f"{GIPPASIM.parent}"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         env={"PYTHONPATH": env_path, "PATH": "/usr/bin:/bin"},
                         check=True)  # fmt: skip
    assert out.stdout.strip() == "[]", out.stdout


@pytest.mark.parametrize(
    "call",
    [
        lambda sc: sc.ora_null(0),
        lambda sc: sc.gsea_null(-5),
        lambda sc: sc.ora_attack_confounded(1000, math.nan),
        lambda sc: sc.ora_attack_confounded(1000, math.inf),
        lambda sc: sc.gsea_attack_operon_packed(1000, 0),
        lambda sc: sc.determinants_planted(1000, {"promoter": 0.5, "bogus": 1.0}, 1.0),
        lambda sc: sc.determinants_planted(1000, EX5, -1.0),
        lambda sc: sc.determinants_planted(1000, {"promoter": math.nan}, 1.0),
        lambda sc: sc.determinants_planted(1000, EX5, 1.0, features="uniform"),
        lambda sc: sc.determinants_null(0),
    ],
    ids=BAD_IDS,
)
def test_bad_arguments_are_refused(call: Any) -> None:
    with pytest.raises(ValueError):
        call(_sc())


# ----------------------------------------------------------------------------- ontology (§1)


def test_synthetic_ontology_parses_with_goatools_and_m803_reader(tmp_path: Path) -> None:
    from goatools.obo_parser import GODag

    genes = [f"G{i:04d}" for i in range(1000)]
    obo, table = _fs().synthetic_ontology(120, 5, np.random.default_rng(0), genes=genes)
    path = tmp_path / "s.obo"
    path.write_text(obo)
    dag = GODag(str(path), optional_attrs={"relationship"}, load_obsolete=False, prt=None)
    live = {k for k, t in dag.items() if k == t.id}
    assert len(live) == 120
    assert all(RE_GO.match(t) for t in live), sorted(live)[:5]
    m803 = _read_go(obo, tmp_path)
    assert set(m803.names) == live
    assert {m803.namespace[t] for t in live} == {"biological_process", "molecular_function",
                                                 "cellular_component"}  # fmt: skip
    # is_a and part_of both used
    assert any(t.parents for k, t in dag.items() if k == t.id)
    assert any(getattr(t, "relationship", {}).get("part_of") for k, t in dag.items()
               if k == t.id)  # fmt: skip
    assert list(table.columns) == ["locus_tag", "db", "id"]
    assert set(table["db"]) == {"GO"}
    assert set(table["id"]) <= live
    assert set(table["locus_tag"]) <= set(genes)


def test_synthetic_ontology_is_a_dag_of_the_requested_depth(tmp_path: Path) -> None:
    obo, _ = _fs().synthetic_ontology(200, 6, np.random.default_rng(1), genes=[])
    dag = _read_go(obo, tmp_path)
    memo: dict[str, int] = {}

    def height(t: str, stack: frozenset[str] = frozenset()) -> int:
        assert t not in stack, f"cycle through {t}"
        if t not in memo:
            ps = dag.parents[t]
            memo[t] = 0 if not ps else 1 + max(height(p, stack | {t}) for p in ps)
        return memo[t]

    longest = max(height(t) for t in dag.names)
    assert longest in (5, 6), longest  # depth in edges or in levels
    for t, ps in dag.parents.items():
        assert ps <= set(dag.names), t  # no dangling parent
        assert {dag.namespace[p] for p in ps} <= {dag.namespace[t]} or ps == frozenset()


def test_synthetic_ontology_annotates_30_to_60_percent_of_genes() -> None:
    genes = [f"G{i:05d}" for i in range(4000)]
    for seed in range(5):
        _, table = _fs().synthetic_ontology(150, 5, np.random.default_rng(seed), genes=genes)
        share = table["locus_tag"].nunique() / len(genes)
        assert 0.3 - 0.02 <= share <= 0.6 + 0.02, (seed, share)


def test_deep_ontology_keeps_terms_with_zero_genes(tmp_path: Path) -> None:
    """Degenerate row: deeper than the gene count supports -> 0-gene terms kept in the OBO."""
    genes = [f"G{i}" for i in range(10)]
    obo, table = _fs().synthetic_ontology(300, 12, np.random.default_rng(2), genes=genes)
    dag = _read_go(obo, tmp_path)
    assert len(dag.names) == 300
    assert table["id"].nunique() < 300  # most terms have no gene, still in the OBO


def test_synthetic_ontology_same_seed_same_output() -> None:
    genes = [f"G{i}" for i in range(500)]
    a = _fs().synthetic_ontology(80, 4, np.random.default_rng(7), genes=genes)
    b = _fs().synthetic_ontology(80, 4, np.random.default_rng(7), genes=genes)
    c = _fs().synthetic_ontology(80, 4, np.random.default_rng(8), genes=genes)
    assert a[0] == b[0]
    pd.testing.assert_frame_equal(a[1], b[1])
    assert (a[0], a[1].to_csv()) != (c[0], c[1].to_csv())


# ----------------------------------------------------------------------------- reserved IDs


def _world_ids(world: Any) -> dict[str, set[str]]:
    x = world.layout.xrefs
    kl = world.kegg_links
    return {
        "go": set(x.loc[x["db"] == "GO", "id"].astype(str)) | set(_obo_ids(world.obo)),
        "ko": set(x.loc[x["db"] == "KEGG", "id"].astype(str))
        | {str(k).removeprefix("ko:") for k in kl["ko"]},
        "map": {str(p).removeprefix("path:") for p in kl["pathway"]},
    }


@pytest.mark.parametrize("which", ["ora", "attack", "gsea", "packed"])
def test_every_synthetic_id_is_in_a_reserved_range(
    which: str, ora4000: Any, attack4000: Any, gsea4000: Any, packed4000: Any
) -> None:
    world = {"ora": ora4000, "attack": attack4000, "gsea": gsea4000, "packed": packed4000}[which]
    ids = _world_ids(world)
    assert ids["go"] and ids["ko"] and ids["map"]
    assert not [t for t in ids["go"] if not RE_GO.match(t)]
    assert not [t for t in ids["ko"] if not RE_KO.match(t)]
    assert all(RE_MAP.match(t) for t in ids["map"]), sorted(ids["map"])
    assert set(world.layout.xrefs["db"]) <= {"COG", "GO", "KEGG"}
    # the KEGG link/name formats: every link row is ``ko:K9xxxx`` -> ``path:map9xxxx``
    assert all(str(k).startswith("ko:K9") for k in world.kegg_links["ko"])
    assert all(str(p).startswith("path:map9") for p in world.kegg_links["pathway"])
    assert 2 <= len(ids["map"]) <= 20
    # every annotated ID is defined: GO in the OBO, KO in the link table
    x = world.layout.xrefs
    assert set(x.loc[x["db"] == "GO", "id"]) <= set(_obo_ids(world.obo))
    assert set(x.loc[x["db"] == "KEGG", "id"]) <= {str(k).removeprefix("ko:")
                                                   for k in world.kegg_links["ko"]}  # fmt: skip


def _fixture_ids() -> dict[str, set[str]]:
    """GO / KO / map IDs that appear in the test GO/KEGG fixtures."""
    files = [*sorted((ROOT / "tests" / "helpers").glob("*.py")),
             *sorted((ROOT / "tests" / "unit" / "function").glob("*.py")),
             ROOT / "tests" / "unit" / "test_setup_resources.py",
             ROOT / "tests" / "property" / "test_collections_props.py",
             *sorted((ROOT / "tests" / "data" / "recorded" / "fgsea").rglob("*.tsv")),
             *sorted((ROOT / "tests" / "data" / "expected").rglob("*.tsv"))]  # fmt: skip
    out: dict[str, set[str]] = {"go": set(), "ko": set(), "map": set()}
    for f in files:
        if not f.is_file():
            continue
        text = f.read_text(errors="replace")
        out["go"] |= set(ANY_GO.findall(text))
        out["ko"] |= set(ANY_KO.findall(text))
        out["map"] |= set(ANY_MAP.findall(text))
    return out


def test_fixture_scan_finds_the_known_fixture_ids() -> None:
    """The collision scan below is not vacuous: it sees the M8-03 fixture's real-range IDs."""
    seen = _fixture_ids()
    assert {"GO:0008150", "GO:0000101", "GO:9999999"} <= seen["go"]
    assert {"K00001", "K99999"} <= seen["ko"]
    assert {"map00010", "map01100"} <= seen["map"]


def test_no_synthetic_id_collides_with_a_real_fixture_id(
    ora4000: Any, attack4000: Any, gsea4000: Any, packed4000: Any
) -> None:
    """Lead decision 1: no synthetic ID equals a real (non-reserved) ID of the test fixtures;
    nor the fixtures' "absent" sentinel ``GO:9999999``."""
    from function_resources_data import ABSENT  # tests/helpers

    fixture = _fixture_ids()
    real = {
        "go": {t for t in fixture["go"] if not RE_GO.match(t)} | {ABSENT},
        "ko": {t for t in fixture["ko"] if not RE_KO.match(t)},
        "map": {t for t in fixture["map"] if not RE_MAP.match(t)},
    }
    assert real["go"] and real["ko"] and real["map"]
    for world in (ora4000, attack4000, gsea4000, packed4000):
        ids = _world_ids(world)
        for kind in ("go", "ko", "map"):
            assert not ids[kind] & real[kind], (kind, sorted(ids[kind] & real[kind])[:5])


def test_simulator_files_carry_no_real_go_kegg_literal() -> None:
    """Acceptance: no real GO/KEGG content in the simulator's function files (a grep for IDs
    outside the reserved prefixes)."""
    for name in ("function_sim.py", "scenarios.py"):
        path = GIPPASIM / name
        assert path.is_file(), name
        text = path.read_text()
        bad = [t for t in ANY_GO.findall(text) if not RE_GO.match(t)]
        bad += [t for t in ANY_KO.findall(text) if not RE_KO.match(t)]
        bad += [t for t in ANY_MAP.findall(text) if not RE_MAP.match(t)]
        assert not bad, (name, bad[:5])
        assert "rest.kegg.jp" not in text and "geneontology.org" not in text, name


# ----------------------------------------------------------------------------- ORA (examples 1, 2)


def test_ora_null_detection_is_covariates_only(ora4000: Any) -> None:
    from gippasim.expression import DEFAULT_BETA

    _, ecfg = _sc().ora_null(4000)
    assert ecfg.detection == "covariates_only"
    assert dict(ecfg.beta) == dict(DEFAULT_BETA)
    assert not ecfg.region_detect_shift


def test_ora_null_categories_cover_the_collections(ora4000: Any) -> None:
    s = ora4000.sets
    assert list(s.columns[:4]) == ["set_id", "collection", "locus_tag", "confounded"]
    assert {"cog", "kegg", *GO_COLLECTIONS} <= set(s["collection"])
    assert not s["confounded"].astype(bool).any()
    assert set(s["locus_tag"]) <= set(ora4000.layout.genes["locus_tag"])
    cog = s.loc[s["collection"] == "cog", "set_id"].astype(str)
    assert set(cog) <= COG_LETTERS
    go = s.loc[s["collection"].isin(GO_COLLECTIONS), "set_id"].astype(str)
    assert all(RE_GO.match(t) for t in go)
    kegg = s.loc[s["collection"] == "kegg", "set_id"].astype(str)
    assert all(RE_MAP.match(t) for t in kegg)


def test_ora_null_members_have_the_same_p_true_as_non_members(ora4000: Any) -> None:
    """One-seed sanity of example 1 (the level claim is in the calibration file): the largest
    category's mean p_true is within 0.05 of the rest's."""
    g = _genes(ora4000)
    biggest = max(_members(ora4000).values(), key=len)
    inside = g["locus_tag"].isin(biggest)
    assert inside.sum() >= 100
    d = g.loc[inside, "p_detect_true"].mean() - g.loc[~inside, "p_detect_true"].mean()
    assert abs(d) < 0.05, d


def _logit(p: Any) -> np.ndarray:
    x = np.asarray(p, dtype=float)
    return np.log(x) - np.log1p(-x)


def _covariate_logit(world: Any) -> tuple[np.ndarray, np.ndarray]:
    """(M4-09 covariate logit ``b0 + b1 log len + b2 log1p n_obs + b3 GRAVY`` from the world's
    *recorded* covariates — ``aa`` length, ``n_obs_true``, ``gravy_true`` — and the mask of genes
    with a protein sequence)."""
    from gippasim.expression import DEFAULT_BETA as B

    g = _genes(world)
    aa = world.layout.genes.set_index("locus_tag").loc[g["locus_tag"], "aa"]
    length = aa.map(lambda a: len(a) if isinstance(a, str) else 0).to_numpy(float)
    ok = length > 0
    n_obs = g["n_obs_true"].to_numpy(float)
    gravy = g["gravy_true"].to_numpy(float)
    logit = np.full(len(g), B["intercept"])
    logit[ok] += (B["log_length"] * np.log(length[ok]) + B["log1p_n_obs"] * np.log1p(n_obs[ok])
                  + B["gravy"] * gravy[ok])  # fmt: skip
    return logit, ok


def _recomputed_covariates(world: Any) -> pd.DataFrame:
    """GIPPA's gene-table covariates (``gippa.genetable.covariates.covariates`` with the default
    ``abundance.enzyme`` / ``observable_len``) recomputed from each gene's ``aa``, indexed like
    ``_genes(world)``; ``has_aa`` marks genes with an amino-acid sequence (pseudogenes without one
    are translated from ``nt`` by GIPPA but carry no covariates in gippasim, so they are
    excluded)."""
    from gippa.config import GippaConfig
    from gippa.genetable.covariates import covariates

    ab = GippaConfig().abundance
    g = _genes(world)
    aa = world.layout.genes.set_index("locus_tag").loc[g["locus_tag"], "aa"].tolist()
    rows = []
    for a in aa:
        if isinstance(a, str) and a:
            c = covariates(a, ab.enzyme, tuple(ab.observable_len))
            rows.append((True, c.length_aa, c.n_observable_peptides, c.gravy))
        else:
            rows.append((False, 0, 0, math.nan))
    return pd.DataFrame(rows, columns=["has_aa", "length_aa", "n_obs", "gravy"])


def _recomputed_logit(world: Any) -> tuple[np.ndarray, np.ndarray]:
    """(M4-09 covariate logit from the covariates *recomputed from ``aa``*, has-aa mask)."""
    from gippasim.expression import DEFAULT_BETA as B

    c = _recomputed_covariates(world)
    ok = c["has_aa"].to_numpy(bool)
    logit = np.full(len(c), B["intercept"])
    logit[ok] += (B["log_length"] * np.log(c["length_aa"].to_numpy(float)[ok])
                  + B["log1p_n_obs"] * np.log1p(c["n_obs"].to_numpy(float)[ok])
                  + B["gravy"] * c["gravy"].to_numpy(float)[ok])  # fmt: skip
    return logit, ok


def _shift_pair() -> tuple[Any, Any, np.ndarray]:
    """(world at member_shift -1.5, world at 0, member mask over ``_genes``), same seed."""
    a = draw(_sc().ora_attack_confounded(1500, -1.5), 21)
    b = draw(_sc().ora_attack_confounded(1500, 0.0), 21)
    conf = _confounded(a)
    assert conf and conf == _confounded(b)
    pd.testing.assert_series_equal(_genes(a)["locus_tag"], _genes(b)["locus_tag"])
    return a, b, _genes(a)["locus_tag"].isin(conf).to_numpy()


def test_example2_population_mean_p_true_at_most_half(ora4000: Any) -> None:
    """Example 2 as a population ratio (lead decision 4): the -1.5 shift applied to *every* gene
    halves the mean p_true at least (statistician: 0.4655, SD 0.002 over seeds). Independent of
    which genes the confounded category draws."""
    p = _genes(ora4000)["p_detect_true"].to_numpy(float)
    q = 1.0 / (1.0 + np.exp(-(_logit(p) - 1.5)))
    assert q.mean() <= 0.5 * p.mean(), (q.mean(), p.mean(), q.mean() / p.mean())


def test_example2_confounded_members_are_less_detectable(attack4000: Any) -> None:
    g = _genes(attack4000)
    conf = _confounded(attack4000)
    assert conf, "no confounded category"
    inside = g["locus_tag"].isin(conf)
    assert 20 <= inside.sum() < len(g) / 2
    pin = g.loc[inside, "p_detect_true"].mean()
    pout = g.loc[~inside, "p_detect_true"].mean()
    assert pin < pout, (pin, pout)
    _, ecfg = _sc().ora_attack_confounded(4000, -1.5)
    assert ecfg.detection == "covariates_only"  # detection through p_true only


def test_confounding_covariate_is_recorded() -> None:
    """B1 (lead decision 6): ``n_obs`` carries the shift, and every member's recorded
    ``n_obs_true`` is GIPPA's observable-peptide count recomputed from its ``aa`` (M8-15 refits
    p_detect from the gene table, which sees only the sequences); likewise every gene with an
    ``aa``."""
    assert _fs().CONFOUND_COVARIATE == "n_obs"
    a, _, inside = _shift_pair()
    rec = _genes(a)["n_obs_true"].to_numpy(np.int64)
    c = _recomputed_covariates(a)
    ok = c["has_aa"].to_numpy(bool)
    new = c["n_obs"].to_numpy(np.int64)
    assert (inside & ok).sum() >= 20
    np.testing.assert_array_equal(rec[inside & ok], new[inside & ok])
    np.testing.assert_array_equal(rec[ok], new[ok])


def test_attack_shift_is_carried_by_the_members_covariate() -> None:
    """B1 (lead decision 6): between shift -1.5 and 0 (same seed) only the confounded members'
    ``aa`` (and ``nt``, derived from it) differ in the layout, and only their ``n_obs_true`` /
    ``p_detect_true`` / ``detected_true`` in the expression; non-members are identical; K/R edits
    leave every member's GRAVY unchanged; movable members (``n_obs`` > 0 at shift 0) lose
    observable peptides, immovable ones stay members, unchanged."""
    from gippa.genetable.covariates import protein_sequence

    a, b, inside = _shift_pair()
    conf = _confounded(a)
    la, lb = a.layout.genes.set_index("locus_tag"), b.layout.genes.set_index("locus_tag")
    pd.testing.assert_index_equal(la.index, lb.index)
    lm = la.index.isin(conf)
    pd.testing.assert_frame_equal(la[~lm], lb[~lm])
    keep = [c for c in la.columns if c not in AA_DERIVED]
    pd.testing.assert_frame_equal(la.loc[lm, keep], lb.loc[lm, keep])
    for name in ("replicons", "tus", "operon_pairs", "regions", "region_genes", "truth", "xrefs",
                 "features"):  # fmt: skip
        pd.testing.assert_frame_equal(getattr(a.layout, name), getattr(b.layout, name))
    # nt stays the coding sequence of the edited aa wherever it was at shift 0
    for lt in la.index[lm]:
        aa_b, nt_b, ft = lb.at[lt, "aa"], lb.at[lt, "nt"], lb.at[lt, "feature_type"]
        if isinstance(aa_b, str) and protein_sequence(None, nt_b, ft) == aa_b.strip("*"):
            assert protein_sequence(None, la.at[lt, "nt"], ft) == str(la.at[lt, "aa"]).strip("*"), (
                lt
            )

    ga, gb = _genes(a), _genes(b)
    pd.testing.assert_frame_equal(ga[~inside].reset_index(drop=True),
                                  gb[~inside].reset_index(drop=True))  # fmt: skip
    keep = [c for c in ga.columns if c not in ("n_obs_true", "p_detect_true", "detected_true")]
    pd.testing.assert_frame_equal(ga.loc[inside, keep], gb.loc[inside, keep])  # gravy_true too

    ca, cb = _recomputed_covariates(a), _recomputed_covariates(b)
    ok = cb["has_aa"].to_numpy(bool)
    m = inside & ok
    np.testing.assert_allclose(ca["gravy"].to_numpy(float)[m], cb["gravy"].to_numpy(float)[m],
                               rtol=0, atol=1e-12)  # fmt: skip
    na, nb = ca["n_obs"].to_numpy(np.int64), cb["n_obs"].to_numpy(np.int64)
    movable = inside & ok & (nb > 0)
    assert movable.sum() >= 20
    # M8-14 decision 8: a few short sORFs with too few K/R cannot lose a peptide; they stay put.
    assert np.all(na[movable] <= nb[movable]), "a member gained observable peptides"
    assert np.mean(na[movable] < nb[movable]) >= 0.95, "too few movable members lost peptides"
    still = la.index[lm & ~la.index.isin(ga.loc[movable, "locus_tag"])]
    pd.testing.assert_frame_equal(la.loc[still], lb.loc[still])  # immovable: untouched
    assert set(still) <= conf  # ... and still counted as members


def test_attack_kr_edits_permute_codons_only() -> None:
    """M8-14 review B1: per edited member the K/R move permutes whole codons, so ``nt`` keeps its
    length, its bases past ``3 * len(aa)`` and its codon multiset."""
    a, b, _ = _shift_pair()
    conf = _confounded(a)
    la, lb = a.layout.genes.set_index("locus_tag"), b.layout.genes.set_index("locus_tag")
    edited = [lt for lt in conf if la.at[lt, "aa"] != lb.at[lt, "aa"]]
    assert len(edited) >= 20
    for lt in edited:
        na, nb, n = la.at[lt, "nt"], lb.at[lt, "nt"], len(str(lb.at[lt, "aa"]))
        if not (isinstance(nb, str) and len(nb) >= 3 * n):
            continue  # no usable nt at shift 0: only aa moves
        assert len(na) == len(nb), lt
        assert na[3 * n :] == nb[3 * n :], lt
        assert sorted(na[3 * i : 3 * i + 3] for i in range(n)) == sorted(
            nb[3 * i : 3 * i + 3] for i in range(n)
        ), lt


@pytest.mark.parametrize("shift", [-4.0, 2.0])
def test_unreachable_member_shift_raises(shift: float) -> None:
    """M8-14 review B2: a shift K/R moves cannot deliver (outside about [-2, +0.5], or more than
    ``SHIFT_TOL`` off on average) raises rather than silently under-delivering."""
    assert _fs().SHIFT_TOL == 0.15
    with pytest.raises(ValueError, match="not reachable"):
        draw(_sc().ora_attack_confounded(1500, shift), 21)


@pytest.mark.parametrize("shift", [-1.5, 0.0])
def test_attack_logit_is_computed_from_covariates_only(shift: float) -> None:
    """B1: p_true is the M4-09 logistic of the recorded covariates for every gene, members
    included; no direct logit term remains, so an M4-13 fit on the covariates can see the shift."""
    w = draw(_sc().ora_attack_confounded(1500, shift), 21)
    logit, ok = _covariate_logit(w)
    assert ok.sum() > 1000
    p = _genes(w)["p_detect_true"].to_numpy(float)
    np.testing.assert_allclose(p[ok], 1.0 / (1.0 + np.exp(-logit[ok])), rtol=0, atol=1e-9)
    # lead decision 6: the recorded covariates are the ones GIPPA recomputes from the sequences
    c = _recomputed_covariates(w)
    assert np.array_equal(c["has_aa"].to_numpy(bool), ok)
    g = _genes(w)
    np.testing.assert_array_equal(g["n_obs_true"].to_numpy(np.int64)[ok],
                                  c["n_obs"].to_numpy(np.int64)[ok])  # fmt: skip
    np.testing.assert_allclose(g["gravy_true"].to_numpy(float)[ok],
                               c["gravy"].to_numpy(float)[ok], rtol=0, atol=1e-12)  # fmt: skip
    rl, _ = _recomputed_logit(w)
    np.testing.assert_allclose(_logit(p)[ok], rl[ok], rtol=0, atol=1e-8)


def test_example2_confounded_categories_are_named_categories(attack4000: Any) -> None:
    s = attack4000.sets
    conf_sets = set(s.loc[s["confounded"].astype(bool), "set_id"])
    assert conf_sets
    # a confounded set is confounded on all its rows (a category property, not a gene one)
    per = s.groupby("set_id")["confounded"].agg(lambda v: set(v.astype(bool)))
    assert all(len(v) == 1 for v in per), per[per.map(len) > 1]


def test_member_shift_zero_reduces_to_ora_null() -> None:
    """Degenerate row: ``member_shift = 0`` gives ``ora_null``'s detection truth."""
    for seed in (11, 12):
        a = draw(_sc().ora_attack_confounded(1500, 0.0), seed)
        b = draw(_sc().ora_null(1500), seed)
        pd.testing.assert_frame_equal(a.layout.genes, b.layout.genes)
        for col in ("p_detect_true", "detected_true"):
            pd.testing.assert_series_equal(_genes(a)[col], _genes(b)[col])


def test_member_shift_moves_the_members_logit_by_the_shift() -> None:
    """B1 (lead decision 6): every member's logit is the covariate model's logit from the
    covariates recomputed from its ``aa`` (1e-8); non-members' and immovable members' logits do
    not move; the mean shift over movable members (``n_obs`` > 0 at shift 0) is -1.5 +- 0.15 (an
    integer ``n_obs`` cannot hit -1.5 exactly)."""
    a, b, inside = _shift_pair()
    ga, gb = _genes(a), _genes(b)
    rl, _ = _recomputed_logit(a)
    np.testing.assert_allclose(_logit(ga["p_detect_true"])[inside], rl[inside], rtol=0, atol=1e-8)
    diff = _logit(ga["p_detect_true"]) - _logit(gb["p_detect_true"])
    np.testing.assert_allclose(diff[~inside], 0.0, atol=1e-8)
    cb = _recomputed_covariates(b)
    movable = inside & cb["has_aa"].to_numpy(bool) & (cb["n_obs"].to_numpy() > 0)
    assert movable.sum() >= 20
    np.testing.assert_allclose(diff[inside & ~movable], 0.0, atol=1e-8)
    assert np.all(diff[movable] <= 1e-8), "a movable member moved towards higher detectability"
    assert np.mean(diff[movable] < 0) >= 0.95, "too few movable members moved"
    mean = float(diff[movable].mean())
    assert abs(mean + 1.5) <= 0.15, (mean, int(movable.sum()), int(inside.sum()))


# ----------------------------------------------------------------------------- GSEA (examples 3, 4)


def test_gsea_null_set_sizes_spread_over_10_to_500_tus(gsea4000: Any) -> None:
    tu = _tu_of(gsea4000.layout)
    sizes = [len({tu[g] for g in m}) for m in _members(gsea4000).values()]
    assert len(sizes) >= 100  # M8-15: >= 100 sets per replicate
    assert min(sizes) >= 10 and max(sizes) <= 500, (min(sizes), max(sizes))
    assert min(sizes) < 50 and max(sizes) > 200, (min(sizes), max(sizes))  # spread
    assert not gsea4000.sets["confounded"].astype(bool).any()


def test_gsea_null_abundance_has_tu_structure(gsea4000: Any) -> None:
    """Abundance from the M4-08 model (TU level + within-TU noise): within-TU ICC well above 0."""
    g = _genes(gsea4000)
    tu = _tu_of(gsea4000.layout)
    g = g.assign(tu=g["locus_tag"].map(tu))
    multi = g.groupby("tu")["abundance_log2_true"].transform("size") >= 2
    within = g[multi].groupby("tu")["abundance_log2_true"].var().mean()
    total = g["abundance_log2_true"].var()
    assert within < 0.6 * total, (within, total)


def test_gsea_null_membership_uncorrelated_with_abundance_one_seed(gsea4000: Any) -> None:
    """One-seed sanity of example 3: the mean over sets of Spearman(membership, abundance) is
    near 0 (|mean| < 0.03; the level claim is in the calibration file)."""
    from scipy import stats

    g = _genes(gsea4000)
    rhos = []
    for m in _members(gsea4000).values():
        ind = g["locus_tag"].isin(m).to_numpy(float)
        rhos.append(stats.spearmanr(ind, g["abundance_log2_true"]).statistic)
    assert abs(float(np.mean(rhos))) < 0.03, np.mean(rhos)


def test_example4_every_packed_set_is_testable_after_tu_collapse(packed4000: Any) -> None:
    """Example 4 as amended (lead decision 4, B2): each packed set holds 10-500 TUs (inside
    ``gsea_size``): 1-3 packed operons of 10-30 genes plus single-gene TUs, whole TUs only."""
    tu = _tu_of(packed4000.layout)
    tus_frame = packed4000.layout.tus
    size = {str(t): len(gs) for t, gs in zip(tus_frame["tu_id"], tus_frame["genes"], strict=True)}
    members = _members(packed4000)
    assert len(members) == 20
    for sid, m in members.items():
        tus = {tu[g] for g in m}
        assert 10 <= len(tus) <= 500, (sid, len(tus))
        assert m == {g for g, t in tu.items() if t in tus}, (sid, "a TU is split")
        packed = sorted(size[t] for t in tus if size[t] > 1)
        assert 1 <= len(packed) <= 3, (sid, packed)
        assert all(10 <= k <= 30 for k in packed), (sid, packed)  # packed operons, 10-30 genes
        assert len(tus) - len(packed) >= 7, (sid, "too few single-gene TUs")


def test_example4_packed_tus_share_their_tu_level(packed4000: Any) -> None:
    """Gene-level ranking sees each packed TU as many correlated genes: within-TU spread is
    small against the between-TU spread of the packed TUs (the shared TU level is the attack)."""
    g = _genes(packed4000).set_index("locus_tag")
    tu = _tu_of(packed4000.layout)
    size = packed4000.layout.tus.set_index("tu_id")["genes"].map(len)
    packed = {tu[x] for m in _members(packed4000).values() for x in m if size[tu[x]] >= 10}
    rows = g.loc[[x for x, t in tu.items() if t in packed]]
    by = rows.groupby(rows.index.map(tu))["abundance_log2_true"]
    assert by.mean().std() > 1.5 * by.std().mean(), (by.mean().std(), by.std().mean())


def test_example4_packed_tus_have_no_extra_effect(packed4000: Any) -> None:
    """B2: packed TUs draw their TU level like every other chr1 TU (no extra N(0, 3^2) effect).
    One-seed sanity (the exchangeability claim is gated in the attack file): the SD of the
    packed TUs' ``tu_log2`` is within 1.4x of the other chr1 TUs' (old design: about 1.65x)."""
    g = _genes(packed4000)
    tu = _tu_of(packed4000.layout)
    tus = packed4000.layout.tus
    chr1 = set(tus.loc[tus["contig"] == "chr1", "tu_id"].astype(str))
    size = tus.set_index("tu_id")["genes"].map(len)
    packed = {tu[x] for m in _members(packed4000).values() for x in m if size[tu[x]] >= 10}
    lvl = g.assign(tu=g["locus_tag"].map(tu)).groupby("tu")["tu_log2"].first()
    inside = lvl.index.isin(packed)
    rest = lvl.index.isin(chr1) & ~inside
    assert inside.sum() >= 20
    ratio = lvl[inside].std() / lvl[rest].std()
    assert ratio < 1.4, ratio


# ----------------------------------------------------------------------------- determinants


def test_example5_recorded_true_r2_is_analytic(planted3000: Any, tmp_path: Path) -> None:
    paths = _fs().write_truth(planted3000, tmp_path)
    assert Path(paths["determinants"]) == tmp_path / "truth" / "determinants.tsv"
    t = pd.read_csv(paths["determinants"], sep="\t")
    assert list(t.columns[:5]) == ["feature", "block", "beta", "true_r2", "residual_var"]
    assert list(t["feature"]) == list(FEATURES)
    assert dict(zip(t["feature"], t["block"], strict=True)) == BLOCKS
    assert dict(zip(t["feature"], t["beta"], strict=True)) == pytest.approx(EX5, abs=1e-12)
    assert t["true_r2"].nunique() == 1 and t["residual_var"].nunique() == 1
    s = sum(b * b for b in EX5.values())
    resid = float(t["residual_var"].iloc[0])
    assert resid >= 1.0 - 1e-12  # at least the noise variance noise_sd^2
    assert float(t["true_r2"].iloc[0]) == pytest.approx(s / (s + resid), abs=1e-6)
    pd.testing.assert_frame_equal(t, planted3000.truth["determinants"], check_dtype=False)


def test_example5_response_follows_the_planted_model(planted3000: Any) -> None:
    g = _genes(planted3000).set_index("locus_tag")
    f = planted3000.features.set_index("locus_tag").loc[g.index]
    assert list(planted3000.features.columns[:4]) == ["locus_tag", *FEATURES]
    z = np.column_stack([_z(f[c]) for c in FEATURES])
    beta = np.array([EX5[c] for c in FEATURES])
    y = g["abundance_log2_true"].to_numpy(float)
    resid = y - y.mean() - z @ beta
    t = planted3000.truth["determinants"]
    want = float(t["residual_var"].iloc[0])
    assert resid.var() == pytest.approx(want, rel=0.15), (resid.var(), want)
    for j in range(3):
        assert abs(np.corrcoef(resid, z[:, j])[0, 1]) < 0.06, FEATURES[j]
    # realised R^2 of the planted signal within MC error of the recorded truth
    realised = (z @ beta).var() / y.var()
    assert realised == pytest.approx(float(t["true_r2"].iloc[0]), abs=0.06)


def test_example5_dominant_block_per_gene(planted3000: Any, tmp_path: Path) -> None:
    paths = _fs().write_truth(planted3000, tmp_path)
    d = pd.read_csv(paths["determinants_genes"], sep="\t")
    assert list(d.columns[:2]) == ["locus_tag", "dominant_block"]
    f = planted3000.features.set_index("locus_tag")
    assert set(d["locus_tag"]) == set(f.index)
    contrib = np.column_stack([np.abs(EX5[c] * _z(f[c])) for c in FEATURES])
    want = dict(zip(f.index, np.array([BLOCKS[c] for c in FEATURES])[contrib.argmax(1)],
                    strict=True))  # fmt: skip
    got = dict(zip(d["locus_tag"], d["dominant_block"], strict=True))
    assert got == want


def test_coefs_all_zero_true_r2_is_zero(tmp_path: Path) -> None:
    w = draw(_sc().determinants_planted(800, dict.fromkeys(FEATURES, 0.0), 1.0,
                                        features="gaussian"), 5)  # fmt: skip
    t = w.truth["determinants"]
    assert (t["beta"] == 0).all()
    assert (t["true_r2"] == 0).all()


def test_coefs_all_zero_equals_null_in_distribution() -> None:
    """Degenerate row: zero coefs -> the null's response law (pooled mean and SD, 4 seeds)."""
    zero = dict.fromkeys(FEATURES, 0.0)
    a = [_genes(draw(_sc().determinants_planted(1500, zero, 1.0), s))["abundance_log2_true"]
         for s in range(4)]  # fmt: skip
    b = [_genes(draw(_sc().determinants_null(1500, noise_sd=1.0), 100 + s))["abundance_log2_true"]
         for s in range(4)]  # fmt: skip
    ya, yb = pd.concat(a), pd.concat(b)
    assert abs(ya.mean() - yb.mean()) < 0.5 * yb.std()
    assert 0.85 < ya.std() / yb.std() < 1.15, (ya.std(), yb.std())


def test_example6_null_keeps_the_layout_and_detection_pattern() -> None:
    for seed in (31, 32):
        n = draw(_sc().determinants_null(3000), seed)
        p = draw(_sc().determinants_planted(3000, EX5, 1.0), seed)
        pd.testing.assert_frame_equal(n.layout.tus, p.layout.tus)
        pd.testing.assert_series_equal(_genes(n)["detected_true"], _genes(p)["detected_true"])
        assert list(n.features.columns[:4]) == ["locus_tag", *FEATURES]


def test_example6_null_response_has_tu_structure() -> None:
    w = draw(_sc().determinants_null(3000), 33)
    g = _genes(w).assign(tu=lambda d: d["locus_tag"].map(_tu_of(w.layout)))
    multi = g.groupby("tu")["abundance_log2_true"].transform("size") >= 2
    within = g[multi].groupby("tu")["abundance_log2_true"].var().mean()
    assert within < 0.8 * g["abundance_log2_true"].var()


def test_planted_features_come_from_the_simulator_layers() -> None:
    """Lead decision 2: the default (recovery-run) features are the planted signals, not
    Gaussians: CAI lies in (0, 1]; promoter and RBS scores are not N(0, 1) draws."""
    from scipy import stats

    w = draw(_sc().determinants_planted(1500, EX5, 1.0), 6)
    f = w.features
    assert f["cai"].between(0.0, 1.0, inclusive="right").all()
    for c in FEATURES:
        assert f[c].notna().all() and f[c].nunique() > 10, c
    gauss = all(stats.kstest(f[c], "norm").pvalue > 1e-3 for c in FEATURES)
    assert not gauss
    t = w.truth["determinants"]
    assert 0.0 < float(t["true_r2"].iloc[0]) < 1.0


def test_small_gene_counts_are_allowed() -> None:
    """Degenerate row: n_genes < 500 allowed (M8-15 fast tier uses 400)."""
    w = draw(_sc().determinants_null(400), 7)
    assert len(w.layout.genes) == 400
    assert len(w.features) == 400


# ----------------------------------------------------------------------------- example 7, TE7


@pytest.mark.parametrize("name", SCENARIOS)
def test_example7_same_seed_identical_outputs(name: str, tmp_path: Path) -> None:
    sc = _sc()
    pair = {
        "ora_null": lambda: sc.ora_null(1200),
        "ora_attack_confounded": lambda: sc.ora_attack_confounded(1200, -1.5),
        "gsea_null": lambda: sc.gsea_null(1200),
        "gsea_attack_operon_packed": lambda: sc.gsea_attack_operon_packed(1200, 5),
        "determinants_null": lambda: sc.determinants_null(1200),
        "determinants_planted": lambda: sc.determinants_planted(1200, EX5, 1.0),
    }[name]
    a, b, c = draw(pair(), 42), draw(pair(), 42), draw(pair(), 43)
    _frames_equal(a, b)
    assert a.obo == b.obo
    pd.testing.assert_frame_equal(a.sets, b.sets)
    pd.testing.assert_frame_equal(a.kegg_links, b.kegg_links)
    assert set(a.truth) == set(b.truth)
    for k in a.truth:
        pd.testing.assert_frame_equal(a.truth[k], b.truth[k])
    if a.features is not None:
        pd.testing.assert_frame_equal(a.features, b.features)
    pa = _fs().write_truth(a, tmp_path / "a")
    pb = _fs().write_truth(b, tmp_path / "b")
    assert {k: Path(v).read_bytes() for k, v in pa.items()} == {
        k: Path(v).read_bytes() for k, v in pb.items()
    }
    differs = not a.expression.genes.equals(c.expression.genes) or a.obo != c.obo
    assert differs, "seed 43 reproduced seed 42"


def _coverage_active() -> bool:
    cov = sys.modules.get("coverage")
    return cov is not None and cov.Coverage.current() is not None


# Example 7's 1 s budget is pinned on Linux without coverage, in CPU time: the macOS CI runner,
# shared by the xdist workers, took 1.1-1.8 s wall clock (#228's run 40023186211), and a Linux
# coverage shard 1.06 s (#273); Slurm's fast suite (no --cov) and Linux CI pass.
@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="timing budget pinned on Linux")
@pytest.mark.skipif(_coverage_active(), reason="timing budget not measured under coverage")
@pytest.mark.parametrize("name", SCENARIOS)
def test_example7_under_one_second_for_4000_genes(name: str) -> None:
    sc = _sc()
    pair = {
        "ora_null": lambda: sc.ora_null(4000),
        "ora_attack_confounded": lambda: sc.ora_attack_confounded(4000, -1.5),
        "gsea_null": lambda: sc.gsea_null(4000),
        "gsea_attack_operon_packed": lambda: sc.gsea_attack_operon_packed(4000, 20),
        "determinants_null": lambda: sc.determinants_null(4000),
        "determinants_planted": lambda: sc.determinants_planted(4000, EX5, 1.0),
    }[name]()
    best = math.inf
    for s in range(3):
        t0 = time.process_time()  # CPU time: robust to waiting for a core under -n auto
        draw(pair, s)
        best = min(best, time.process_time() - t0)
    assert best < 1.0, best
