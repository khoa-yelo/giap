"""M10-05 acceptance criterion 10: the pilot metrics on two real local GIPPA runs of gippasim.

Spec: docs/specs/M10-05-pilot-metrics.md, acceptance criterion 10 (``slow``; Slurm via
``scripts/slurm/test_branch.sbatch``, never the session allocation).

Two ``gippa run --executor local`` runs on gippasim SIM100K (``tests/helpers/m4_real.py``,
``pi_wrong = 0`` so every simulated identification is a true one):

* "A2": the simulator's same-strain peptides plus two planted ones, an out-of-frame translation
  inside a ``+`` strand CDS and the N-terminal product (M-initial) of a strict non-ATG gene that
  the simulator did not already emit;
* "A1": the simulator's table with ``search_db = other``, UniProt-like ``protein_ids`` and a few
  residues substituted so that tolerant and unmapped peptides exist.

Then ``gippa compare``, ``metrics.py compute`` and ``render``. Oracle: TE1's only miss is the
planted out-of-frame peptide; the non-ATG N-terminal products in D2 (the planted one and the
simulator's own, structural misses before M11-20) map by the initiator credit (ADR-0068), so
every other peptide is exact in frame. Metric 2's A1 ``tolerant + unmapped`` equals
``genome_inconsistent_rate * n``; metric 3's prophage, BGC and plasmid subsets equal gippasim's
planted region gene sets; metric 7 has one entry per manifest step. A test ``$GIPPA_DB_DIR``
only; no network.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc
import m4_real as mr
from pilot_metrics_m10 import load_te_proxies, run_cli, write_json

pytestmark = pytest.mark.slow

SEED = 20261023


def _translate_out_of_frame(genome: dict[str, str], g: Any) -> str | None:
    """A 12-aa stop-free stretch of frame +1 inside ``+`` strand gene ``g``, or None."""
    from gippa.mapping.sixframe import translate11_text

    iv = g.interval
    seq = genome[iv.contig][iv.start + 1 : iv.end - 2]
    prot = translate11_text(seq)
    for i in range(10, len(prot) - 12):
        cand = prot[i : i + 12]
        if "*" not in cand and "X" not in cand:
            return cand
    return None


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    from click.testing import CliRunner

    from gippa.cli import cli
    from gippa.io.fasta import read_genome_fasta

    root = tmp_path_factory.mktemp("m10_05_sim_pilot")
    db = mc.test_db(root / "db")
    sim = mr.sim100k(root / "sim", pi_wrong=0.0)
    te = load_te_proxies()
    genes, _ = te.read_coding_genes(Path(sim.bakta_json), Path(sim.genome_fasta))
    genome = {r.id: r.sequence for r in read_genome_fasta(Path(sim.genome_fasta))}
    base = pd.read_csv(sim.peptides, sep="\t", dtype=str, keep_default_na=False)
    base["modified_peptide"] = base["peptide"]
    present = set(base["peptide"])

    non_atg = te.non_atg_nterm(genes)
    planted_nonatg = next(p for p in sorted(non_atg) if p not in present and 7 <= len(p) <= 30)
    plus = [g for g in sorted(genes, key=lambda x: x.locus_tag)
            if g.feature_type == "cds" and g.interval.strand == "+" and len(g.aa) > 60]  # fmt: skip
    oof_gene, oof = next((g, s) for g in plus if (s := _translate_out_of_frame(genome, g)))
    template = base.iloc[0].to_dict()
    planted = [
        {**template, "peptide": oof, "modified_peptide": oof, "protein_ids": oof_gene.locus_tag,
         "q_value": "0.0001", "spectrum_id": "run1.planted1"},
        {**template, "peptide": planted_nonatg, "modified_peptide": planted_nonatg,
         "protein_ids": non_atg[planted_nonatg], "q_value": "0.0001",
         "spectrum_id": "run1.planted2"},
    ]  # fmt: skip
    a2 = pd.concat([base, pd.DataFrame(planted)], ignore_index=True)
    pep_a2 = root / "inputs" / "a2.tsv"
    pep_a2.parent.mkdir(parents=True)
    a2.to_csv(pep_a2, sep="\t", index=False)

    a1 = base.copy()
    a1["search_db"] = "other"
    accessions = {t: f"P{i:05d}" for i, t in enumerate(sorted(set(a1["protein_ids"])))}
    a1["protein_ids"] = a1["protein_ids"].map(accessions)
    distinct = sorted(set(a1["peptide"]))
    swap = {}
    for k, p in enumerate(distinct[:13]):
        if k < 10:  # one substitution: tolerant (or unmapped)
            mid = len(p) // 2
            new = "W" if p[mid] != "W" else "Y"
            swap[p] = p[:mid] + new + p[mid + 1 :]
        else:  # reversed: unmapped
            swap[p] = p[::-1]
    a1["peptide"] = a1["peptide"].map(lambda s: swap.get(s, s))
    a1["modified_peptide"] = a1["peptide"]
    pep_a1 = root / "inputs" / "a1.tsv"
    a1.to_csv(pep_a1, sep="\t", index=False)

    out_a1, out_a2, cmp = root / "out_a1", root / "out_a2", root / "cmp_a1_a2"
    mr.gippa_run(dataclasses.replace(sim, peptides=pep_a2), out_a2, db, SEED)
    mr.gippa_run(dataclasses.replace(sim, peptides=pep_a1), out_a1, db, SEED)
    with mc.db_env(db):
        r = CliRunner().invoke(cli, ["compare", str(out_a1), str(out_a2), "--labels", "A1,A2",
                                     "--outdir", str(cmp)])  # fmt: skip
    assert r.exit_code == 0, r.output

    convs, sums = {}, {}
    for arm, table in (("a1", a1), ("a2", a2)):
        q = pd.to_numeric(table["q_value"], errors="coerce")
        target = ~table["is_decoy"].str.lower().isin(["true", "1"])
        n = int(table.loc[target & (q <= 0.01), "peptide"].nunique())
        counts = {"distinct_peptides_target_q01": n, "q_crossings": 0}
        convs[arm] = write_json(root / "inputs" / f"{arm}.convert.json", {"counts": counts})
        sums[arm] = write_json(root / "diann" / arm / "summary.json", {"distinct_stripped_q01": n})

    k12 = root / "k12.fasta"
    coding = sorted((g for g in genes if g.feature_type in ("cds", "sorf")),
                    key=lambda x: x.locus_tag)  # fmt: skip
    k12.write_text("".join(f">sp|Q{i:05d}|K12\n{g.aa}\n" for i, g in enumerate(coding[::2])),
                   encoding="utf-8")  # fmt: skip

    metrics = root / "metrics"
    res = run_cli(
        "compute", "--a1", str(out_a1), "--a2", str(out_a2), "--cmp", str(cmp),
        "--peptides-a1", str(pep_a1), "--peptides-a2", str(pep_a2),
        "--convert-a1", str(convs["a1"]), "--convert-a2", str(convs["a2"]),
        "--summary-a1", str(sums["a1"]), "--summary-a2", str(sums["a2"]),
        "--bakta-json", str(sim.bakta_json), "--genome", str(sim.genome_fasta),
        "--k12-fasta", str(k12), "--out", str(metrics), cwd=root,
    )  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr
    rend = run_cli("render", "--metrics", str(metrics / "metrics.json"),
                   "--out", str(metrics / "REPORT.md"), cwd=root)  # fmt: skip
    return {"root": root, "sim": sim, "metrics": metrics, "render": rend, "a2": a2,
            "oof": oof, "planted_nonatg": planted_nonatg, "non_atg": non_atg,
            "out_a1": out_a1, "out_a2": out_a2}  # fmt: skip


def _doc(w: dict[str, Any]) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads((w["metrics"] / "metrics.json").read_text(encoding="utf-8"))
    return doc


def _d2(w: dict[str, Any]) -> set[str]:
    from gippa.config import GippaConfig
    from gippa.io.peptides import read_peptide_table, validate_peptides
    from gippa.steps.mapping import select_queries

    man = json.loads((w["out_a2"] / "run_manifest.json").read_text(encoding="utf-8"))
    cfg = GippaConfig.model_validate(man["config"])  # the run's own config (review S7)
    raw, _ = read_peptide_table(w["root"] / "inputs" / "a2.tsv")
    table, _ = validate_peptides(raw, cfg.input)  # GIPPA's mapping front end (review S2)
    return set(select_queries(table, cfg)["peptide"].astype(str))  # no Cont_ rows


def test_te1_misses_are_the_planted_truth(world: dict[str, Any]) -> None:
    te1 = _doc(world)["te1"]
    d2 = _d2(world)
    structural = d2 & set(world["non_atg"])
    assert world["planted_nonatg"] in structural
    assert te1["n_peptides"] == len(d2)
    assert te1["non_atg_nterm"]["peptides"] == len(structural)
    assert te1["non_atg_nterm"]["misses"] == 0  # M11-20: credited
    assert te1["n_peptides"] - te1["exact_inframe"] == 1
    exc = pd.read_csv(world["metrics"] / "te1_exceptions.tsv", sep="\t", dtype=str,
                      keep_default_na=False)  # fmt: skip
    assert set(exc["peptide"]) == {world["oof"]}


def test_metric_2_a1_equals_genome_inconsistent_rate(world: dict[str, Any]) -> None:
    m2 = _doc(world)["m2"]["A1"]
    gi = m2["genome_inconsistent"]
    assert m2["counts"]["tolerant"] + m2["counts"]["unmapped"] > 0  # the substitutions bite
    assert gi["count"] == round(gi["value"] * gi["n"])
    assert m2["counts"]["tolerant"] + m2["counts"]["unmapped"] == gi["count"]
    assert gi["consistent"] is True


def _subset(cov: Any, name: str) -> dict[str, Any]:
    if isinstance(cov, dict) and "subsets" in cov:
        cov = cov["subsets"]
    if isinstance(cov, dict):
        return cov[name]  # type: ignore[no-any-return]
    return next(e for e in cov if e.get("subset") == name)  # type: ignore[no-any-return]


def test_metric_3_subsets_equal_the_planted_regions(world: dict[str, Any]) -> None:
    genes = mc.read(world["out_a2"] / "tables" / "genes.tsv")
    g = set(genes.loc[genes["feature_type"].isin(["cds", "sorf"]), "locus_tag"])
    regions = pd.read_csv(world["sim"].truth["regions"], sep="\t", dtype=str,
                          keep_default_na=False)  # fmt: skip
    planted: dict[str, set[str]] = {"prophage": set(), "bgc": set(), "plasmid": set()}
    for rtype, tags in zip(regions["region_type"], regions["genes"], strict=True):
        if rtype in planted and tags not in ("", "NA"):
            planted[rtype] |= set(tags.split(";"))
    m3 = _doc(world)["m3"]["A2"]
    for name, tags in planted.items():
        assert tags & g, name
        assert _subset(m3, name)["genes"] == len(tags & g), name
    total = sum(_subset(m3, k)["genes"] for k in ("strain_specific", "k12_shared", "no_peptides"))
    assert total == len(g) == _subset(m3, "all")["genes"]


def test_metric_7_one_entry_per_manifest_step(world: dict[str, Any]) -> None:
    doc = _doc(world)
    for arm, out in (("A1", world["out_a1"]), ("A2", world["out_a2"])):
        man = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
        assert [s["name"] for s in doc["m7"][arm]["steps"]] == [s["name"] for s in man["steps"]]
        assert doc["m7"][arm]["hosts"] == []  # local executor: no trace.tsv


def test_hypotheses_and_render(world: dict[str, Any]) -> None:
    doc = _doc(world)
    assert doc["h1"]["verdict"] in {"supported", "not_supported", "not_evaluable"}
    assert doc["h3"]["verdict"] in {"supported", "not_supported", "not_evaluable"}
    rend = world["render"]
    assert rend.returncode == 0, rend.stdout + rend.stderr
    text = (world["metrics"] / "REPORT.md").read_text(encoding="utf-8")
    assert world["oof"] not in text and world["planted_nonatg"] not in text
    assert "**missing:" not in text  # every table is present: no marker (review B2)
    assert (
        doc["gippa_commit"]
        == json.loads((world["out_a2"] / "run_manifest.json").read_text(encoding="utf-8"))[
            "git_commit"
        ]
    )  # review S4
