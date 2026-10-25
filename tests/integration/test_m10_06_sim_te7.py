"""M10-06 acceptance criterion 12: TE7 and TE5 on two real local GIAP runs of giapsim.

Spec: docs/specs/M10-06-te-real-and-te7.md, acceptance criterion 12 (``slow``; Slurm via
``scripts/slurm/test_branch.sbatch``, never the session allocation).

Two ``giap run --executor local`` runs of giapsim SIM100K (``pi_wrong = 0``) with the same seed
in two directories, on a peptide table with planted trypsin and keratin digests from
``tests/data/contaminants_mini.faa`` (accessions renamed to the real library's ``Cont_`` form;
the same file is the test ``$GIAP_DB_DIR`` contaminant library). Then ``te_real.py te7``
(``pass``), ``te_real.py compute`` (TE5 ``pass``; ``|S| - genome_exact`` equals the planted family
peptides that the independent oracle ``tests/helpers/contaminant_oracle.py`` finds nowhere in
the genome; TE2 runs on the run's ``variants.tsv``), and a copy of one run with one
``tables/genes.tsv`` cell changed (``fail``, ``numeric_only`` on that column). No network.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import gene_table_fixture as gt
import m4_real as mr
from contaminant_oracle import hits_genome, il, six_frames, tryptic_digest
from helpers.te_real_m10 import REPO, edit_manifest, git, run_te_real, write_contaminants
from mapping_steps import write_test_db

pytestmark = pytest.mark.slow

SEED = 20261024
FAMILY_SOURCES = ("Cont_P00761", "Cont_P04264", "Cont_P13645")  # trypsin, KRT1, KRT10


def _read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif line.strip():
            seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    from giap.io.protein_fasta import resolve_protein_id

    root = tmp_path_factory.mktemp("m10_06_sim_te7")
    contam = write_contaminants(root / "lib" / "contaminants.fasta")
    db = write_test_db(root / "db", contaminants=contam)
    gt.write_cog_def(db / "cog" / "cog-20.def.tab")  # the COG fixture, as m4_chain.test_db
    man = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    man["cog2020_def"] = {"version": "2020", "path": "cog/cog-20.def.tab"}
    (db / "db_manifest.json").write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")
    sim = mr.sim100k(root / "sim", pi_wrong=0.0)
    genome = _read_fasta(Path(sim.genome_fasta))
    frames = [il(f) for seq in genome.values() for f in six_frames(seq)]

    lib = {resolve_protein_id(k)[0]: v for k, v in _read_fasta(contam).items()}
    planted: dict[str, str] = {}
    for acc in FAMILY_SOURCES:
        for p in tryptic_digest(lib[acc], 10, 30):  # >= 10 aa: no chance tolerant placement
            planted.setdefault(p, acc)
    assert len(planted) >= 20
    base = pd.read_csv(sim.peptides, sep="\t", dtype=str, keep_default_na=False)
    targets = (
        base[~base["is_decoy"].str.lower().isin(["true", "1"])] if "is_decoy" in base else base
    )
    template = targets.iloc[0].to_dict()
    rows = [
        {
            **template,
            "peptide": p,
            "protein_ids": acc,
            "q_value": "0.0001",
            **({"modified_peptide": p} if "modified_peptide" in base else {}),
            **({"spectrum_id": f"run1.cont{i}"} if "spectrum_id" in base else {}),
        }
        for i, (p, acc) in enumerate(sorted(planted.items()))
    ]
    table = pd.concat([base, pd.DataFrame(rows)], ignore_index=True)
    pep = root / "inputs" / "a2.tsv"
    pep.parent.mkdir(parents=True)
    table.to_csv(pep, sep="\t", index=False)

    s = dataclasses.replace(sim, peptides=pep)
    rep1, rep2 = root / "runs" / "rep1", root / "runs" / "rep2"
    mr.giap_run(s, rep1, db, SEED)
    mr.giap_run(s, rep2, db, SEED)
    head = git(REPO, "rev-parse", "HEAD")
    for out in (rep1, rep2):
        # the test checkout may be dirty or the manifest commit unset: the pair under test is
        # the byte comparison, so pin both manifests to the checkout's HEAD, clean
        edit_manifest(out, lambda m: m.update(git_commit=head, git_dirty=False))
    not_in_genome = {p for p in planted if not hits_genome(p, frames)}
    return {
        "root": root,
        "rep1": rep1,
        "rep2": rep2,
        "pep": pep,
        "contam": contam,
        "planted": planted,
        "not_in_genome": not_in_genome,
    }


def _te7(w: dict[str, Any], ref: Path, rep: Path, out: Path) -> dict[str, Any]:
    res = run_te_real(
        "te7", "--ref", str(ref), "--rep", str(rep), "--repo", str(REPO), "--out", str(out)
    )
    assert res.returncode == 0, res.stdout + res.stderr
    doc: dict[str, Any] = json.loads((out / "te7.json").read_text(encoding="utf-8"))
    return doc


def test_two_runs_same_seed_pass_te7(world: dict[str, Any]) -> None:
    doc = _te7(world, world["rep1"], world["rep2"], world["root"] / "m_te7")
    g = doc["giap"]
    assert g["status"] == "pass", g["differing"]
    assert g["n_gate"] > 0 and g["n_excluded"] > 0


def test_te5_pass_and_s_is_the_planted_truth(world: dict[str, Any]) -> None:
    out = world["root"] / "m_compute"
    res = run_te_real(
        "compute",
        "--a2",
        str(world["rep1"]),
        "--a1",
        str(world["rep2"]),
        "--peptides-a1",
        str(world["pep"]),
        "--peptides-a2",
        str(world["pep"]),
        "--contaminants",
        str(world["contam"]),
        "--out",
        str(out),
    )
    assert res.returncode == 0, res.stdout + res.stderr
    doc = json.loads((out / "te_real.json").read_text(encoding="utf-8"))
    for arm in ("A1", "A2"):
        a = doc["te5"][arm]
        assert a["n_s"] == len(world["planted"]), arm
        assert a["n_s"] - a["counts"]["genome_exact"] == len(world["not_in_genome"]), arm
        assert a["status"] == "pass", a
    assert doc["te5"]["status"] == "pass"
    te2 = doc["te2"]
    if (world["rep1"] / "tables" / "variants.tsv").is_file():
        with (world["rep1"] / "tables" / "variants.tsv").open(newline="", encoding="utf-8") as fh:
            n = sum(1 for _ in csv.DictReader(fh, delimiter="\t"))
        assert te2["status"] in ("pass", "fail")
        assert te2["counts"]["alt_sites"] == n
    else:
        assert te2["status"] == "not_evaluable"
    raw = (out / "te_real.json").read_text(encoding="utf-8")
    assert not any(p in raw for p in world["planted"])


def test_one_changed_genes_cell_is_numeric_only(world: dict[str, Any]) -> None:
    rep3 = world["root"] / "runs" / "rep3"
    shutil.copytree(world["rep2"], rep3)
    path = rep3 / "tables" / "genes.tsv"
    lines = path.read_text(encoding="utf-8").split("\n")
    header = lines[0].split("\t")
    changed = None
    for i in range(1, len(lines)):
        cells = lines[i].split("\t")
        for col in ("abundance_log2", "length_aa", "mw_da", "gravy"):
            j = header.index(col) if col in header else -1
            if j >= 0 and cells[j] not in ("", "NA"):
                cells[j] = repr(float(cells[j]) + 0.5)
                lines[i] = "\t".join(cells)
                changed = col
                break
        if changed:
            break
    assert changed, "no numeric cell to change in genes.tsv"
    path.write_text("\n".join(lines), encoding="utf-8")
    doc = _te7(world, world["rep1"], rep3, world["root"] / "m_te7_changed")
    g = doc["giap"]
    assert g["status"] == "fail"
    d = {x["path"]: x for x in g["differing"]}
    assert list(d) == ["tables/genes.tsv"]
    assert d["tables/genes.tsv"]["class"] == "numeric_only"
    assert d["tables/genes.tsv"]["columns"] == [changed]
    assert d["tables/genes.tsv"]["n_cells"] == 1
