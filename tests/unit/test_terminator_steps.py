"""Steps ``transterm_input`` and ``terminator_layer`` (M2-15; M2-index step contract rows
``transterm_input`` and ``terminator_layer``; GIAP_SPEC.md §6.2, §6.5, §7.1.5; Example 9).

Run through ``giap.steps.registry.run_step`` after ``parse_bakta`` on giapsim's SIM100K (the genome
the recorded sim100k TerminatorNet/TransTermHP runs used; same contig names and lengths).

Conventions pinned here where the contract table is silent:

* Step keys are not fixed by the spec; inputs are resolved by name (a key containing
  ``terminatornet`` gets ``terminatornet.tsv``, ``transterm``/``tt`` the ``.tt`` file, ``coords``
  ``genes.coords``; ``genes``/``replicons``/``genome`` the parse_bakta outputs). Outputs: the key
  containing ``crosscheck`` is ``terminators_crosscheck.tsv``; the other ``terminator_layer`` output
  is ``terminators.parquet``.
* ``terminators.parquet`` (work table, 0-based half-open, §6.5) has at least the M2-02 precomputed
  ``terminators`` columns ``contig, start, end, strand, confidence, source`` and holds the primary
  tool's calls (+ BacTermFinder when enabled), never the cross-check tool's: ``terminator_between``
  counts every non-BacTermFinder call, so TransTermHP calls in the default layer would leak into
  Module 4 flags.
* ``terminators_crosscheck.tsv`` header is exactly the Q2 columns; ``all`` is the last row.
* Example 9 (``--terminators transtermhp``): the layer holds TransTermHP calls with CONF/100; since
  ``run_step`` requires every declared output, "cross-check table not written" is accepted as either
  an absent file or a header-only file.
"""

from __future__ import annotations

import importlib
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow.parquet as pq
import pytest

import giapsim
from giap.annotate.bakta import parse_bakta_json
from giap.config import GiapConfig, resolve_config
from giap.model.genes import Gene

DATA = Path(__file__).resolve().parents[1] / "data" / "recorded" / "terminators"
TN_SIM = DATA / "terminatornet" / "sim100k" / "terminatornet.tsv"
TT_SIM = DATA / "transtermhp" / "sim100k" / "transtermhp.tt"
PARSE_OUT = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}
LAYER_COLS = {"contig", "start", "end", "strand", "confidence", "source"}
XC_HEADER = ["contig", "n_terminatornet", "n_transtermhp", "n_matched", "share_matched"]


def _fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"giap.steps.{n}"))
    return reg


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield _fresh_steps(mp)


class Sim:
    def __init__(self, res: Any, work: Path) -> None:
        self.res = res
        self.work = work
        self.genes: list[Gene] = list(parse_bakta_json(res.bakta_json).model.genes)
        self.avail: dict[str, Path] = {k: work / v for k, v in PARSE_OUT.items()}
        self.avail["genome"] = res.genome_fasta


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Sim:
    out = tmp_path_factory.mktemp("sim100k")
    res = giapsim.simulate(giapsim.SIM100K, giapsim.SIM100K_SEED, out / "sim")
    work = out / "work"
    with pytest.MonkeyPatch.context() as mp:
        reg = _fresh_steps(mp)
        reg.run_step(
            "parse_bakta",
            GiapConfig(seed=1),
            {"bakta_json": res.bakta_json, "genome": res.genome_fasta},
            {k: work / v for k, v in PARSE_OUT.items()},
        )
    return Sim(res, work)


def _resolve(key: str, avail: dict[str, Path], extra: dict[str, Path]) -> Path | None:
    if key in avail:
        return avail[key]
    k = key.lower()
    if "terminatornet" in k:
        return extra.get("terminatornet")
    if "transterm" in k or k in {"tt", "transtermhp"}:
        return extra.get("transtermhp")
    if "coords" in k:
        return extra.get("coords")
    if "genome" in k:
        return avail["genome"]
    return None


def _inputs(reg: ModuleType, name: str, sim: Sim, extra: dict[str, Path]) -> dict[str, Path]:
    spec = reg.get(name)
    ins: dict[str, Path] = {}
    for key in spec.inputs:
        p = _resolve(key, sim.avail, extra)
        assert p is not None, f"{name}: cannot provide required input {key!r}"
        ins[key] = p
    for key in spec.optional_inputs:
        p = _resolve(key, sim.avail, extra)
        if p is not None:
            ins[key] = p
    return ins


def _layer_outputs(reg: ModuleType, out: Path) -> tuple[dict[str, Path], Path, Path]:
    keys = reg.get("terminator_layer").outputs
    xc = [k for k in keys if "crosscheck" in k]
    other = [k for k in keys if "crosscheck" not in k]
    assert len(xc) == 1 and len(other) == 1, f"terminator_layer outputs: {keys}"
    outs = {xc[0]: out / "terminators_crosscheck.tsv", other[0]: out / "terminators.parquet"}
    return outs, outs[other[0]], outs[xc[0]]


def _run_layer(
    reg: ModuleType, sim: Sim, out: Path, cfg: GiapConfig | None = None
) -> tuple[Path, Path]:
    extra = {"terminatornet": TN_SIM, "transtermhp": TT_SIM}
    outs, parquet, xc = _layer_outputs(reg, out)
    reg.run_step(
        "terminator_layer",
        cfg or GiapConfig(seed=1),
        _inputs(reg, "terminator_layer", sim, extra),
        outs,
    )
    return parquet, xc


def _tn_rows() -> list[tuple[str, int, int, str, float]]:
    out = []
    for ln in TN_SIM.read_text().splitlines():
        f = ln.split("\t")
        if len(f) == 11:
            out.append((f[0], int(f[8]) - 1, int(f[9]), f[10], float(f[7])))
    return sorted(out)


def _tt_rows() -> list[tuple[str, int, int, str, float]]:
    out, contig = [], ""
    for ln in TT_SIM.read_text().splitlines():
        if ln.startswith("SEQUENCE "):
            contig = ln.split()[1]
        elif ln.startswith("  TERM "):
            f = ln.split()
            a, b = int(f[2]), int(f[4])
            out.append((contig, min(a, b) - 1, max(a, b), f[5], int(f[7]) / 100))
    return sorted(out)


def _layer(path: Path) -> list[dict[str, Any]]:
    rows = pq.read_table(path).to_pylist()
    assert rows, "terminators.parquet is empty"
    assert set(rows[0]) >= LAYER_COLS, f"missing columns {LAYER_COLS - set(rows[0])}"
    return rows


# ------------------------------------------------------------------- registration


def test_steps_are_registered_per_the_contract_table(steps: ModuleType) -> None:
    ti = steps.get("transterm_input")
    assert ti.kind == "py" and ti.precomputed == "terminators"
    assert len(ti.outputs) == 1
    tl = steps.get("terminator_layer")
    assert tl.kind == "py" and tl.precomputed == "terminators"
    assert any("crosscheck" in k for k in tl.outputs)
    assert any("bactermfinder" in k.lower() for k in tl.optional_inputs)
    for ext in ("terminatornet", "transtermhp", "bactermfinder"):
        spec = steps.get(ext)
        assert spec.kind == "ext" and spec.precomputed == "terminators"
    assert steps.get("bactermfinder").module == "opt_in"


# ------------------------------------------------------------------- transterm_input


def test_transterm_input_writes_the_coords_of_every_gene(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    (key,) = steps.get("transterm_input").outputs
    out = tmp_path / "genes.coords"
    steps.run_step(
        "transterm_input",
        GiapConfig(seed=1),
        _inputs(steps, "transterm_input", sim, {}),
        {key: out},
    )
    want = []
    for g in sim.genes:
        iv = g.interval
        a, b = (iv.start + 1, iv.end) if iv.strand == "+" else (iv.end, iv.start + 1)
        want.append(f"{g.locus_tag} {a} {b} {iv.contig}")
    assert sorted(out.read_text().splitlines()) == sorted(want)
    assert out.read_text().endswith("\n")


def test_transterm_input_is_deterministic(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    (key,) = steps.get("transterm_input").outputs
    ins = _inputs(steps, "transterm_input", sim, {})
    steps.run_step("transterm_input", GiapConfig(seed=1), ins, {key: tmp_path / "a.coords"})
    steps.run_step("transterm_input", GiapConfig(seed=99), ins, {key: tmp_path / "b.coords"})
    assert (tmp_path / "a.coords").read_bytes() == (tmp_path / "b.coords").read_bytes()


# ------------------------------------------------------------------- terminator_layer


def test_layer_holds_terminatornet_calls_zero_based(
    steps: ModuleType, sim: Sim, tmp_path: Path
) -> None:
    parquet, _ = _run_layer(steps, sim, tmp_path)
    rows = _layer(parquet)
    assert {r["source"] for r in rows} == {"TerminatorNet"}  # TransTermHP is cross-check only
    got = sorted((r["contig"], r["start"], r["end"], r["strand"], r["confidence"]) for r in rows)
    assert got == _tn_rows()
    assert ("chr1", 1887, 1907, "+", 0.8812014060630281) in got  # "chr1 TERM_1 … 1888 1907 +"


def test_layer_crosscheck_table(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    _, xc = _run_layer(steps, sim, tmp_path)
    lines = xc.read_text().splitlines()
    assert lines[0].split("\t") == XC_HEADER
    rows = {f[0]: f[1:] for f in (ln.split("\t") for ln in lines[1:])}
    assert lines[-1].split("\t")[0] == "all"
    tn, tt = _tn_rows(), _tt_rows()
    matched = sum(
        any(
            o[0] == p[0] and o[3] == p[3] and max(0, max(o[1], p[1]) - min(o[2], p[2])) <= 20
            for o in tt
        )
        for p in tn
    )
    n_chr1 = sum(p[0] == "chr1" for p in tn)
    assert [int(v) for v in rows["chr1"][:3]] == [n_chr1, 30, matched]
    assert float(rows["chr1"][3]) == pytest.approx(matched / n_chr1, abs=1e-4)
    assert [int(v) for v in rows["plasmid_1"][:3]] == [3, 0, 0]
    assert [int(v) for v in rows["all"][:3]] == [133, 30, matched]
    assert float(rows["all"][3]) == pytest.approx(matched / 133, abs=1e-4)


def test_layer_is_byte_deterministic(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    p1, x1 = _run_layer(steps, sim, tmp_path / "a", GiapConfig(seed=1))
    p2, x2 = _run_layer(steps, sim, tmp_path / "b", GiapConfig(seed=7))
    assert p1.read_bytes() == p2.read_bytes()
    assert x1.read_bytes() == x2.read_bytes()


def test_example9_transtermhp_as_primary(steps: ModuleType, sim: Sim, tmp_path: Path) -> None:
    cfg = resolve_config(None, ["annotation.terminators.tool=transtermhp", "seed=1"])
    parquet, xc = _run_layer(steps, sim, tmp_path, cfg)
    rows = _layer(parquet)
    assert {r["source"] for r in rows} == {"TransTermHP"}
    got = sorted((r["contig"], r["start"], r["end"], r["strand"]) for r in rows)
    want = _tt_rows()
    assert got == [w[:4] for w in want]
    assert sorted(r["confidence"] for r in rows) == pytest.approx(sorted(w[4] for w in want))
    assert not xc.exists() or len(xc.read_text().splitlines()) <= 1
