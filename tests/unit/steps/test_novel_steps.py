"""The novel-ORF steps and the second round under the local executor (issue M7-08;
docs/specs/M7-08-novel-steps-and-second-round.md, ADR-0024; GIAP_SPEC.md §6.2, §6.4, §7.6.2).

Covered: the registration of ``novel_candidates``, ``promotech_novel`` (ext),
``novel_promoter_layer``, ``ostir_novel`` (ext), ``novel_rbs_layer`` and ``novel_orfs``; the local
plan (rule 5, example 6, "no peptides"); rules 1-4 on the work files of one real local run, with
hand-written Promotech / OSTIR outputs and hand-built layers; examples 1 and 6 through ``giap run``.

Readings pinned by the test author (see tests/helpers/novel_steps_fixture.py for the position
convention of ``rbs_novel``):

* ``novel_rbs_layer``: the spec declares inputs ``ostir_novel`` and ``ostir_input_novel`` only, but
  the OSTIR input ids (``NOVEL_<contig>_<lowest+1>_<p|m>``) do not give the genome position of a
  minus-strand start; the step may therefore declare further ``novel_candidates`` outputs as
  inputs (the registration test allows exactly that superset).
* rule 4: ``promoter_score`` of an ORF with one start is the maximum of the Promotech scores of
  that start's windows over ``novel_promoter_windows`` and the first-round ``promoter_windows``,
  joined on ``(contig, strand, end3p)`` -- the window map's ``end3p`` and the layer's ``end3p``
  share the M2-13 convention, and a TSV layer carries the same numbers 1-based (M2-02).
* rule 5 / example 6: the exit code of a missing layer is the one ``giap run --executor local``
  gives every missing layer today (``ExitCode.INVALID_INPUT`` = 1, pinned by
  tests/unit/test_local_executor.py); the spec's "exit 2" contradicts its own "exactly as a missing
  ``rbs`` is today" and is read as the latter.
* the run flag ``layer_missing:<key>`` is a WARNING record of a ``giap.steps`` logger, as
  ``gene_table`` and ``context_layers`` emit theirs.
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import novel_steps_fixture as nf

from giap.exitcodes import ExitCode

NOVEL_STEPS = ("novel_candidates", "promotech_novel", "novel_promoter_layer", "ostir_novel",
               "novel_rbs_layer", "novel_orfs")  # fmt: skip
CANDIDATE_OUTPUTS = {"novel_candidates", "novel_windows", "novel_windows_map", "ostir_input_novel"}
# name -> (kind, required inputs, optional inputs, outputs, precomputed)
SPEC: dict[str, tuple[str, set[str], set[str], set[str], str | None]] = {
    "novel_candidates": (
        "py",
        {"mapping", "peptides", "genes", "features", "replicons", "genome", "annotation_meta"},
        {"contaminants"}, CANDIDATE_OUTPUTS, None,
    ),
    "promotech_novel": ("ext", {"novel_windows"}, set(), {"promotech_novel"}, "promoters"),
    "novel_promoter_layer": ("py", {"novel_windows_map", "promotech_novel"}, set(),
                             {"novel_promoter_windows"}, "promoters"),
    "ostir_novel": ("ext", {"ostir_input_novel"}, set(), {"ostir_novel"}, "rbs_novel"),
    "novel_rbs_layer": ("py", {"ostir_novel", "ostir_input_novel"}, set(), {"rbs_novel"},
                        "rbs_novel"),
    "novel_orfs": (
        "py", {"novel_candidates", "gene_table", "genes", "replicons", "genome"},
        {"novel_promoter_windows", "rbs_novel", "promoter_windows", "terminators",
         "peptides", "features", "hits_exact", "hits_tolerant", "placement_costs"},  # M7-11b
        {"novel_orfs", "novel_orfs_gff3", "novel_orfs_faa", "annotation_patch",
         "novel_class_fdr"}, None,
    ),
}  # fmt: skip
NOVEL_LOCUS_CLASSES = {"nterm_extension", "cterm_extension", "intergenic", "cds_outframe",
                       "antisense", "pseudogene", "sorf"}  # fmt: skip


def _reg() -> Any:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    return reg


# ---------------------------------------------------------------- registration


@pytest.mark.parametrize("name", NOVEL_STEPS)
def test_step_is_registered_as_specified(name: str) -> None:
    spec = _reg().get(name)
    kind, req, opt, outs, pre = SPEC[name]
    assert spec.kind == kind
    assert spec.module == "core"
    assert set(spec.outputs) == outs
    assert spec.precomputed == pre
    assert set(spec.optional_inputs) == opt
    if name == "novel_rbs_layer":  # may also read novel_candidates outputs (module docstring)
        assert req <= set(spec.inputs) <= req | CANDIDATE_OUTPUTS, spec.inputs
    else:
        assert set(spec.inputs) == req
    if kind == "py":
        assert spec.func is not None and spec.func.__module__ == "giap.steps.novel"
    else:
        assert spec.func is None


def test_novel_orfs_needs_the_gene_table() -> None:
    """Lead decision 3: NOVEL runs only with a gene table (M7-09's RBS quantile)."""
    assert "gene_table" in _reg().get("novel_orfs").inputs


# ---------------------------------------------------------------- local plan (rule 5)


def _layers_cfg(tmp_path: Path, *, rbs_novel: bool = True, promoters: bool = True) -> Any:
    from giap.config import resolve_config

    keys = ["bakta_json", "operon_pairs", "tx_rates", "terminators", "rbs", "regions"]
    keys += ["promoters"] if promoters else []
    keys += ["rbs_novel"] if rbs_novel else []
    over = []
    for k in keys:
        p = tmp_path / f"{k}.tsv"
        p.write_text("x\n")
        over.append(f"annotation.precomputed.{k}={p}")
    if not rbs_novel:
        over.append("annotation.precomputed.rbs_novel=null")
    return resolve_config(None, [*over, "modules=[core]"])


ROOTS = {"genome", "peptides", "bakta_json"}


def _plan(cfg: Any, roots: set[str]) -> list[str]:
    from giap.runner import local

    reg = _reg()
    return local.plan_local(cfg, [s.name for s in reg.all_steps()], roots)


def test_local_plan_replaces_the_second_round_with_layers(tmp_path: Path) -> None:
    plan = _plan(_layers_cfg(tmp_path), ROOTS)
    assert "novel_candidates" in plan and "novel_orfs" in plan
    for name in ("promotech_novel", "novel_promoter_layer", "ostir_novel", "novel_rbs_layer"):
        assert name not in plan, name
    assert plan.index("detectability") < plan.index("novel_orfs"), plan
    assert plan.index("novel_candidates") < plan.index("novel_orfs"), plan


def test_example_6_plan_without_rbs_novel_names_the_key(tmp_path: Path) -> None:
    from giap.runner import local

    with pytest.raises(local.LayerMissing) as info:
        _plan(_layers_cfg(tmp_path, rbs_novel=False), ROOTS)
    assert info.value.key == "annotation.precomputed.rbs_novel"


def test_without_promoters_layer_the_local_plan_names_promoters(tmp_path: Path) -> None:
    from giap.runner import local

    with pytest.raises(local.LayerMissing) as info:
        _plan(_layers_cfg(tmp_path, promoters=False), ROOTS)
    assert info.value.key == "annotation.precomputed.promoters"


def test_no_peptides_no_novel_step_and_no_rbs_novel_needed(tmp_path: Path) -> None:
    """Degenerate case "no peptides": the NOVEL block is not entered, so a run without peptides
    needs no ``rbs_novel`` layer."""
    _reg().get("novel_candidates")
    plan = _plan(_layers_cfg(tmp_path, rbs_novel=False), {"genome", "bakta_json"})
    assert not set(plan) & set(NOVEL_STEPS), plan


# ---------------------------------------------------------------- one real local run


@pytest.fixture(scope="module")
def run1(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, Path, Path]:
    root = tmp_path_factory.mktemp("m7_08_sim100k")
    sim = nf.simulate(root / "sim")
    out = root / "OUT"
    with nf.db_env(root):
        code, log = nf.run_giap(sim, out)
    assert code == 0, log
    return sim, out, out / "work"


@pytest.fixture(scope="module")
def cfg1(run1: tuple[Any, Path, Path]) -> Any:
    from giap.config import resolve_config

    return resolve_config(Path(run1[0].config_yaml), [f"seed={nf.SEED}"])


@pytest.fixture(scope="module")
def orfs1(run1: tuple[Any, Path, Path]) -> list[Any]:
    sim, _, work = run1
    orfs = nf.library_orfs(work, sim)
    assert any(o.strand == "+" and len(o.starts) >= 2 for o in orfs)
    assert any(o.strand == "-" and len(o.starts) >= 2 for o in orfs)
    assert sum(len(o.starts) == 1 for o in orfs) >= 2
    return orfs


@pytest.fixture(scope="module")
def cands(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path_factory: pytest.TempPathFactory
) -> dict[str, Path]:
    """``novel_candidates`` run on the work files of ``run1`` into a fresh directory."""
    _, _, work = run1
    out = tmp_path_factory.mktemp("cands")
    return nf.run("novel_candidates", cfg1, nf.step_inputs("novel_candidates", work), out)


# ---------------------------------------------------------------- example 1


def test_example_1_four_files_published(run1: tuple[Any, Path, Path]) -> None:
    _, out, _ = run1
    for name in nf.TABLES:
        assert (out / "tables" / name).is_file(), name


def test_example_1_one_row_per_reconstructed_orf(
    run1: tuple[Any, Path, Path], orfs1: list[Any]
) -> None:
    from giap.io.novel_tables import read_novel_orfs

    _, out, work = run1
    table = read_novel_orfs(out / "tables" / "novel_orfs.tsv")
    assert sorted(table["orf_id"]) == sorted(o.orf_id for o in orfs1)
    cand = pd.read_parquet(work / "novel_candidates")
    assert set(cand["orf_id"]) == set(table["orf_id"])
    nf.assert_valid_novel_orfs(out / "tables" / "novel_orfs.tsv")


def test_example_1_best_ostir_start_set_where_a_start_exists(
    run1: tuple[Any, Path, Path],
) -> None:
    """The simulator's ``rbs_novel`` layer scores every start codon, so every ORF with a start
    gets ``start_best_ostir`` and ``ostir_tir``: the layer value at that start, at least the
    value at the longest start."""
    sim, out, _ = run1
    rows = nf.tsv_dicts(out / "tables" / "novel_orfs.tsv")
    layer = {(r["contig"], r["strand"], r["start"]): float(r["tir_ostir"] or "nan")
             for r in nf.tsv_dicts(nf.rbs_novel_layer(sim))}  # fmt: skip
    with_start = [r for r in rows if r["start_longest"] is not None]
    assert with_start
    for r in rows:
        if r["start_longest"] is None:
            assert r["start_best_ostir"] is None and r["ostir_tir"] is None, r["orf_id"]
            continue
        assert r["start_best_ostir"] is not None, r["orf_id"]
        assert r["ostir_tir"] is not None, r["orf_id"]
        best = layer[(r["contig"], r["strand"], r["start_best_ostir"])]
        assert float(r["ostir_tir"]) == pytest.approx(best), r["orf_id"]
        assert best >= layer[(r["contig"], r["strand"], r["start_longest"])], r["orf_id"]
        assert "ostir_missing" not in (r["flags"] or "").split(";"), r["orf_id"]


def test_example_1_gff3_files_are_valid(run1: tuple[Any, Path, Path]) -> None:
    from helpers.novel_outputs_fixture import check_gff3

    sim, out, _ = run1
    reps = dict(nf.annotation_of(sim).model.replicons)
    feats = check_gff3(out / "tables" / "novel_orfs.gff3", reps)
    assert len(feats) == len(nf.tsv_dicts(out / "tables" / "novel_orfs.tsv"))
    check_gff3(out / "tables" / "annotation_patch.gff3", reps)
    faa = (out / "tables" / "novel_orfs.faa").read_text(encoding="utf-8")
    assert faa.count(">") >= sum(
        r["start_longest"] is not None for r in nf.tsv_dicts(out / "tables" / "novel_orfs.tsv")
    )


# ---------------------------------------------------------------- example 6


def test_example_6_local_run_without_rbs_novel_exits_naming_the_key(
    run1: tuple[Any, Path, Path], tmp_path: Path
) -> None:
    sim = run1[0]
    with nf.db_env(tmp_path):
        code, log = nf.run_giap(sim, tmp_path / "OUT", "--set",
                                "annotation.precomputed.rbs_novel=null")  # fmt: skip
    assert code == ExitCode.INVALID_INPUT, log
    assert "annotation.precomputed.rbs_novel is required with --executor local" in log, log
    assert not (tmp_path / "OUT" / "tables" / "novel_orfs.tsv").exists()


# ---------------------------------------------------------------- rule 1: novel_candidates


def test_rule_1_candidates_one_row_per_orf_and_peptide(
    cands: dict[str, Path], orfs1: list[Any]
) -> None:
    df = pd.read_parquet(cands["novel_candidates"])
    assert {"orf_id", "peptide"} <= set(df.columns)
    got = {(o, p.upper().replace("I", "L")) for o, p in zip(df["orf_id"], df["peptide"],
                                                            strict=True)}  # fmt: skip
    assert got == {(o.orf_id, p) for o in orfs1 for p in o.peptides}
    assert not df.duplicated(["orf_id", "peptide"]).any()


def test_rule_1_windows_and_map_are_the_m7_06_second_round(
    run1: tuple[Any, Path, Path], cfg1: Any, cands: dict[str, Path], orfs1: list[Any]
) -> None:
    from giap.novel.candidates import promoter_windows, window_chunks

    sim = run1[0]
    ann = nf.annotation_of(sim)
    chunks, want = window_chunks(promoter_windows(orfs1, ann.model, dict(ann.sequences), cfg1))
    d = cands["novel_windows"]
    assert d.is_dir()
    files = sorted(d.glob("*.fa"))
    assert [f.name for f in files] == [f"windows_{i:03d}.fa" for i in range(len(chunks))]
    headers = [ln[1:].strip() for f in files for ln in f.read_text().splitlines()
               if ln.startswith(">")]  # fmt: skip
    assert headers == [h for c in chunks for h, _ in c]
    got = pd.read_parquet(cands["novel_windows_map"])
    cols = list(want.columns)
    assert set(cols) <= set(got.columns), got.columns
    assert list(got[cols].itertuples(index=False, name=None)) == list(
        want.itertuples(index=False, name=None)
    )
    assert set(got["orf_id"]) == {o.orf_id for o in orfs1 if o.starts}


def test_rule_1_ostir_input_is_one_multi_start_row_per_orf_with_a_start(
    run1: tuple[Any, Path, Path], cfg1: Any, cands: dict[str, Path], orfs1: list[Any]
) -> None:
    from giap.annotate.ostir import read_ostir_input
    from giap.novel.candidates import ostir_rows

    sim, _, work = run1
    ann = nf.annotation_of(sim)
    anti_sd = json.loads((work / "annotation_meta").read_text())["anti_sd"]
    want = ostir_rows(orfs1, ann.model, dict(ann.sequences), anti_sd, cfg1)
    path = cands["ostir_input_novel"]
    assert path.read_text().splitlines()[0] == nf.OSTIR_INPUT_HEADER
    assert read_ostir_input(path) == want
    assert any(r.end > r.start for r in want)


# ---------------------------------------------------------------- rule 2: the two layer steps


def _scores(windows: Path) -> dict[str, float]:
    out = {}
    for fa in sorted(windows.glob("windows_*.fa")):
        for ln in fa.read_text().splitlines():
            if ln.startswith(">"):
                h = ln[1:].strip()
                out[h] = round(((int(h[1:]) * 7919) % 1_000_003) / 1_000_003, 6)
    return out


def _promotech_dirs(root: Path, windows: Path) -> list[Path]:
    """One ``promotech_NNN/sequences_predictions.csv`` per window chunk (the M2-13 layout)."""
    scores = _scores(windows)
    dirs = []
    for fa in sorted(windows.glob("windows_*.fa")):
        d = root / fa.stem.replace("windows", "promotech")
        d.mkdir(parents=True)
        lines = [nf.PROMOTECH_HEADER]
        recs = fa.read_text().split(">")[1:]
        for rec in recs:
            head, seq = rec.split("\n", 1)
            lines.append(f"{head.strip()}\t{seq.replace(chr(10), '')}\t{scores[head.strip()]}")
        (d / "sequences_predictions.csv").write_text("\n".join(lines) + "\n")
        dirs.append(d)
    return dirs


def test_rule_2_novel_promoter_layer_is_promoter_layer_on_the_novel_map(
    run1: tuple[Any, Path, Path], cfg1: Any, cands: dict[str, Path], tmp_path: Path
) -> None:
    from giap.annotate.promoters import promoter_layer

    work = run1[2]
    dirs = _promotech_dirs(tmp_path / "pt", cands["novel_windows"])
    ins = nf.step_inputs("novel_promoter_layer", work,
                         {**cands, "promotech_novel": dirs})  # fmt: skip
    outs = nf.run("novel_promoter_layer", cfg1, ins, tmp_path / "o")
    got = pd.read_parquet(outs["novel_promoter_windows"])
    want = promoter_layer(pd.read_parquet(cands["novel_windows_map"]),
                          _scores(cands["novel_windows"]))  # fmt: skip
    assert list(got.columns) == ["contig", "strand", "end3p", "score"]
    assert list(got.itertuples(index=False, name=None)) == list(
        want.itertuples(index=False, name=None)
    )


START_CODONS = ("ATG", "GTG", "TTG")


def _window_starts(r: Any) -> list[int]:
    """1-based window positions of the in-frame start codons in ``[start, end]``."""
    return [p for p in range(r.start, r.end + 1, 3) if r.seq[p - 1 : p + 2] in START_CODONS]


def test_rule_2_novel_rbs_layer_keys_every_scored_start_by_position(
    run1: tuple[Any, Path, Path], cfg1: Any, cands: dict[str, Path], orfs1: list[Any],
    tmp_path: Path,
) -> None:  # fmt: skip
    from giap.annotate.ostir import read_ostir_input

    sim, _, work = run1
    genome = nf.genome_of(sim)
    by_id = {o.orf_id: o for o in orfs1}
    rows, want, dropped = [], {}, None
    for i, r in enumerate(read_ostir_input(cands["ostir_input_novel"])):
        o = by_id[r.id]
        ps = _window_starts(r)
        assert len(ps) == len(o.starts), r.id  # every window start codon is an ORF start
        for j, (p, x) in enumerate(zip(ps, o.starts, strict=True)):
            tir = 1000.0 * (i + 1) + p
            if dropped is None and len(ps) >= 2 and j == len(ps) - 1:
                dropped = (o.contig, o.strand, nf.tpos(o.strand, x))  # OSTIR wrote no row here
                continue
            rows.append(nf.ostir_out_row(r.id, r.seq[p - 1 : p + 2], p, tir))
            want[(o.contig, o.strand, nf.tpos(o.strand, x))] = tir
        if i == 0:  # a row at an unrequested (out-of-frame) position is ignored
            rows.append(nf.ostir_out_row(r.id, "ATG", ps[0] + 1, 1.0))
    rows.append(nf.ostir_out_row("NOVEL_chr1_1_p", "ATG", 51, 2.0))  # an id nobody asked for
    assert dropped is not None
    csv = nf.write_ostir_csv(tmp_path / "ostir_novel.csv", rows)
    ins = nf.step_inputs("novel_rbs_layer", work, {**cands, "ostir_novel": csv})
    outs = nf.run("novel_rbs_layer", cfg1, ins, tmp_path / "o")
    got = pd.read_parquet(outs["rbs_novel"])
    assert tuple(got.columns) == nf.RBS_NOVEL_COLUMNS
    keys = list(zip(got["contig"], got["strand"], got["start"].astype(int), strict=True))
    assert len(keys) == len(set(keys))
    assert dict(zip(keys, got["tir_ostir"], strict=True)) == pytest.approx(want)
    assert dropped not in keys
    for (contig, strand, start), codon, sd in zip(keys, got["start_codon"], got["sd_dg"],
                                                  strict=True):  # fmt: skip
        assert nf.codon_at(genome[contig], strand, start) == codon
        assert sd == pytest.approx(-want[(contig, strand, start)] / 1000, abs=1e-4)


def test_rule_2_empty_ostir_output_gives_an_empty_layer(
    run1: tuple[Any, Path, Path], cfg1: Any, cands: dict[str, Path], tmp_path: Path
) -> None:
    """Degenerate case: OSTIR wrote no row at any requested start."""
    csv = nf.write_ostir_csv(tmp_path / "ostir_novel.csv", [])
    ins = nf.step_inputs("novel_rbs_layer", run1[2], {**cands, "ostir_novel": csv})
    got = pd.read_parquet(nf.run("novel_rbs_layer", cfg1, ins, tmp_path / "o")["rbs_novel"])
    assert len(got) == 0 and tuple(got.columns) == nf.RBS_NOVEL_COLUMNS


# ---------------------------------------------------------------- rule 3: empty second round


def _no_candidate_mapping(work: Path, dst: Path) -> Path:
    from giap.io.tsv import read_tsv, write_tsv

    shutil.copytree(work / "mapping", dst)
    header, rows = read_tsv(dst / "peptide_hits.tsv")
    k = header.index("locus_class")
    kept = [r for r in rows if r[k] not in NOVEL_LOCUS_CLASSES]
    assert len(kept) < len(rows)
    write_tsv(dst / "peptide_hits.tsv", header, kept)
    return dst


def test_rule_3_no_candidate_empty_second_round_and_empty_tables(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    from helpers.novel_outputs_fixture import NOVEL_ORFS_COLUMNS

    work = run1[2]
    mapping = _no_candidate_mapping(work, tmp_path / "mapping")
    c = nf.run("novel_candidates", cfg1,
               nf.step_inputs("novel_candidates", work, {"mapping": mapping}),
               tmp_path / "c")  # fmt: skip
    assert len(pd.read_parquet(c["novel_candidates"])) == 0
    assert c["novel_windows"].is_dir() and not any(c["novel_windows"].iterdir())
    assert len(pd.read_parquet(c["novel_windows_map"])) == 0
    assert c["ostir_input_novel"].read_text() == nf.OSTIR_INPUT_HEADER + "\n"

    skip = ("novel_promoter_windows", "rbs_novel", "promoter_windows", "terminators")
    ins = nf.step_inputs("novel_orfs", work, {"novel_candidates": c["novel_candidates"]}, skip)
    o = nf.run("novel_orfs", cfg1, ins, tmp_path / "o")
    assert o["novel_orfs"].read_text() == "\t".join(NOVEL_ORFS_COLUMNS) + "\n"
    for key in ("novel_orfs_gff3", "annotation_patch"):
        lines = o[key].read_text().splitlines()
        assert lines[0] == "##gff-version 3", key
        assert all(ln.startswith("#") for ln in lines), key
    assert ">" not in o["novel_orfs_faa"].read_text()


# ---------------------------------------------------------------- rule 4: novel_orfs joins


def _rule4_layers(
    cands: dict[str, Path], orfs: list[Any]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Hand-built second-round layers and the expected joins.

    * rbs: the ORF ``tie`` (multi-start) has its two most upstream starts tied at the top
      (-> the more upstream); every other multi-start ORF peaks at its most downstream start;
      the single-start ORF ``none`` has no row (-> NA, ``ostir_missing``).
    * promoters: windows owned by single-start ORF ``p`` score 0.6 in the second round and one
      of them 0.9 in the first round (-> 0.9); single-start ORF ``q``'s own windows 0.4 (-> 0.4);
      every other window 0.2; one first-round row lies at no candidate window.
    """
    multi = [o for o in orfs if len(o.starts) >= 2]
    single = [o for o in orfs if len(o.starts) == 1]
    assert multi and len(single) >= 3
    tie, none_, p, q = multi[0], single[0], single[1], single[2]
    rbs_rows, best = [], {}
    for o in orfs:
        if not o.starts or o is none_:
            continue
        n = len(o.starts)
        tirs = (
            [5000.0, 5000.0] + [100.0] * (n - 2) if o is tie else [100.0 + 10 * j for j in range(n)]
        )
        for x, tir in zip(o.starts, tirs, strict=True):  # start_codon: _fill_codons
            rbs_rows.append([o.contig, o.strand, nf.tpos(o.strand, x), tir, -tir / 1000, -10.5,
                             0.25, -0.5, 7, None])  # fmt: skip
        k = 0 if o is tie else n - 1
        best[o.orf_id] = (nf.tpos(o.strand, o.starts[k]), tirs[k])
    rbs = pd.DataFrame(rbs_rows, columns=list(nf.RBS_NOVEL_COLUMNS))
    wmap = pd.read_parquet(cands["novel_windows_map"])
    uniq = wmap.drop_duplicates(["contig", "strand", "end3p"])
    score = {p.orf_id: 0.6, q.orf_id: 0.4}
    novel = pd.DataFrame({
        "contig": uniq["contig"], "strand": uniq["strand"],
        "end3p": uniq["end3p"].astype("int64"),
        "score": [score.get(o, 0.2) for o in uniq["orf_id"]],
    }).reset_index(drop=True)  # fmt: skip
    pw = uniq[uniq["orf_id"] == p.orf_id].iloc[len(uniq[uniq["orf_id"] == p.orf_id]) // 2]
    taken = set(zip(uniq["contig"], uniq["strand"], uniq["end3p"], strict=True))
    far = next(e for e in range(1, 10_000) if (p.contig, "+", e) not in taken)
    first = pd.DataFrame({"contig": [pw["contig"], p.contig], "strand": [pw["strand"], "+"],
                          "end3p": [int(pw["end3p"]), far], "score": [0.9, 0.99]})  # fmt: skip
    want = {"best": best, "none": none_.orf_id, "promoter": {p.orf_id: 0.9, q.orf_id: 0.4}}
    return rbs, novel, first, want


def _fill_codons(rbs: pd.DataFrame, genome: dict[str, str]) -> pd.DataFrame:
    rbs = rbs.copy()
    rbs["start_codon"] = [nf.codon_at(genome[c], s, int(p)) for c, s, p in
                          zip(rbs["contig"], rbs["strand"], rbs["start"], strict=True)]  # fmt: skip
    return rbs.astype({"start": "int64", "spacer_bp": "int64"})


def _check_rule4(table: pd.DataFrame, want: dict[str, Any]) -> None:
    by = {r["orf_id"]: r for r in table.to_dict("records")}
    for oid, (pos, tir) in want["best"].items():
        assert by[oid]["start_best_ostir"] == pos, oid
        assert by[oid]["ostir_tir"] == pytest.approx(tir), oid
    r = by[want["none"]]
    assert pd.isna(r["start_best_ostir"]) and pd.isna(r["ostir_tir"])
    assert "ostir_missing" in str(r["flags"]).split(";")
    for oid, s in want["promoter"].items():
        assert by[oid]["promoter_score"] == pytest.approx(s), oid


@pytest.fixture(scope="module")
def rule4(
    run1: tuple[Any, Path, Path], cands: dict[str, Path], orfs1: list[Any]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    rbs, novel, first, want = _rule4_layers(cands, orfs1)
    return _fill_codons(rbs, nf.genome_of(run1[0])), novel, first, want


def _parquet_cfg(sim: Any) -> Any:
    from giap.config import resolve_config

    return resolve_config(Path(sim.config_yaml), [
        f"seed={nf.SEED}", "annotation.precomputed.promoters=null",
        "annotation.precomputed.rbs_novel=null", "annotation.precomputed.terminators=null",
    ])  # fmt: skip


def test_rule_4_joins_scores_by_position_from_parquet_layers(
    run1: tuple[Any, Path, Path], cands: dict[str, Path], rule4: Any, tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:  # fmt: skip
    from giap.io.novel_tables import read_novel_orfs

    sim, _, work = run1
    rbs, novel, first, want = rule4
    paths = {}
    for key, df in (("rbs_novel", rbs), ("novel_promoter_windows", novel),
                    ("promoter_windows", first)):  # fmt: skip
        paths[key] = tmp_path / f"{key}.parquet"
        df.to_parquet(paths[key], index=False)
    ins = nf.step_inputs("novel_orfs", work, {"novel_candidates": cands["novel_candidates"],
                                              **paths}, skip=("terminators",))  # fmt: skip
    with caplog.at_level(logging.WARNING, logger="giap.steps"):
        o = nf.run("novel_orfs", _parquet_cfg(sim), ins, tmp_path / "o")
    _check_rule4(read_novel_orfs(o["novel_orfs"]), want)
    msgs = [r.getMessage() for r in caplog.records]
    assert any("layer_missing:terminators" in m for m in msgs), msgs
    assert not any("layer_missing:rbs_novel" in m for m in msgs), msgs


def _tsv_layers(tmp_path: Path, rule4: Any) -> tuple[Path, Path]:
    """The same layers as M2-02 TSVs (1-based positions): one promoters file serves both rounds
    (the first-round 0.9 row replaces the second-round score at its position)."""
    from giap.io.tsv import write_tsv

    rbs, novel, first, _ = rule4
    rbs_tsv = tmp_path / "rbs_novel.tsv"
    write_tsv(rbs_tsv, list(nf.RBS_NOVEL_COLUMNS),
              [[c, s, int(p) + 1, *rest] for c, s, p, *rest in
               rbs.itertuples(index=False, name=None)])  # fmt: skip
    merged = pd.concat([first, novel]).drop_duplicates(["contig", "strand", "end3p"])
    merged = merged.sort_values(["contig", "end3p", "strand"])
    prom_tsv = tmp_path / "promoters.tsv"
    write_tsv(prom_tsv, ["contig", "strand", "end3p", "score"],
              [list(r) for r in merged.itertuples(index=False, name=None)])  # fmt: skip
    return rbs_tsv, prom_tsv


def test_rule_5_novel_orfs_accepts_the_precomputed_tsv_layers(
    run1: tuple[Any, Path, Path], cands: dict[str, Path], rule4: Any, tmp_path: Path
) -> None:
    """Under the local executor ``novel_orfs`` receives the layers' TSVs (``load_precomputed``)
    instead of the Parquet files; the join is the same."""
    from giap.config import resolve_config
    from giap.io.novel_tables import read_novel_orfs

    sim, _, work = run1
    rbs_tsv, prom_tsv = _tsv_layers(tmp_path, rule4)
    cfg = resolve_config(Path(sim.config_yaml), [
        f"seed={nf.SEED}", f"annotation.precomputed.promoters={prom_tsv}",
        f"annotation.precomputed.rbs_novel={rbs_tsv}", "annotation.precomputed.terminators=null",
    ])  # fmt: skip
    ins = nf.step_inputs("novel_orfs", work, {
        "novel_candidates": cands["novel_candidates"], "rbs_novel": rbs_tsv,
        "novel_promoter_windows": prom_tsv, "promoter_windows": prom_tsv,
    }, skip=("terminators",))  # fmt: skip
    o = nf.run("novel_orfs", cfg, ins, tmp_path / "o")
    _check_rule4(read_novel_orfs(o["novel_orfs"]), rule4[3])


def test_rule_4_missing_optional_layers_give_na_and_run_flags(
    run1: tuple[Any, Path, Path], cands: dict[str, Path], orfs1: list[Any], tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:  # fmt: skip
    from giap.io.novel_tables import read_novel_orfs

    sim, _, work = run1
    skip = ("novel_promoter_windows", "rbs_novel", "promoter_windows", "terminators")
    ins = nf.step_inputs("novel_orfs", work, {"novel_candidates": cands["novel_candidates"]}, skip)
    with caplog.at_level(logging.WARNING, logger="giap.steps"):
        o = nf.run("novel_orfs", _parquet_cfg(sim), ins, tmp_path / "o")
    t = read_novel_orfs(o["novel_orfs"])
    assert sorted(t["orf_id"]) == sorted(x.orf_id for x in orfs1)
    for col in ("start_best_ostir", "ostir_tir", "promoter_score", "terminator_downstream",
                "terminator_conf", "terminator_distance_nt"):  # fmt: skip
        assert t[col].isna().all(), col
    with_start = {x.orf_id for x in orfs1 if x.starts}
    for r in t.to_dict("records"):
        if r["orf_id"] in with_start:
            assert "ostir_missing" in str(r["flags"]).split(";"), r["orf_id"]
    msgs = " | ".join(r.getMessage() for r in caplog.records)
    assert "layer_missing:rbs_novel" in msgs and "layer_missing:terminators" in msgs, msgs
    assert any(f"layer_missing:{k}" in msgs for k in
               ("novel_promoter_windows", "promoters", "promoter_windows")), msgs  # fmt: skip


def test_rule_4_outputs_byte_identical_on_rerun(
    run1: tuple[Any, Path, Path], cands: dict[str, Path], rule4: Any, tmp_path: Path
) -> None:
    """TE7 at the step level: the same inputs give byte-identical files."""
    sim, _, work = run1
    rbs_tsv, prom_tsv = _tsv_layers(tmp_path, rule4)
    from giap.config import resolve_config

    cfg = resolve_config(Path(sim.config_yaml), [
        f"seed={nf.SEED}", f"annotation.precomputed.promoters={prom_tsv}",
        f"annotation.precomputed.rbs_novel={rbs_tsv}",
    ])  # fmt: skip
    ins = nf.step_inputs("novel_orfs", work, {
        "novel_candidates": cands["novel_candidates"], "rbs_novel": rbs_tsv,
        "novel_promoter_windows": prom_tsv, "promoter_windows": prom_tsv,
    })  # fmt: skip
    a = nf.run("novel_orfs", cfg, ins, tmp_path / "a")
    b = nf.run("novel_orfs", cfg, ins, tmp_path / "b")
    for key in a:
        assert a[key].read_bytes() == b[key].read_bytes(), key
