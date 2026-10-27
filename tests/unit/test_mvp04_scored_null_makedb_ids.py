"""gippasim ``decoys = "scored_null"``: target rows name makedb entries whatever ``wrong_source``
(MVP-04 gate 12; statistician diagnosis local/stat/mvp04-gates-2-12.md; M7-11b decision 9 in the
module docstring of tests/sim/gippasim/proteomics.py; MVP-04 step 5: the simulated search is against
``gippa makedb --six-frame-orfs``).

``scenarios.k12_evidence`` uses ``decoys = "scored_null"`` with the default ``wrong_source =
"shuffled"``. Run 20261023_a1 showed the true novel-ORF target rows carrying gippasim's truth IDs
``NOVEL_xxxx``: GIPPA then finds no makedb entry for them (``no_makedb_pair``) and the class FDR of
almost every true ORF is unestimated. In scored_null mode the search database is the makedb FASTA,
so (as in ``sixframe_noncoding``):

* a true novel ORF's peptide names the makedb entry it lies in (``GIPPAORF_…``), never ``NOVEL_``;
* a true protein's peptide reaching upstream of its annotated start names the gene's
  ``nterm_extension`` entry;
* the scored_null wrong-peptide draw does not depend on ``wrong_source`` (it already draws its own
  protein-entry and six-frame nulls), so with the fix the scored_null outputs are the same for
  ``shuffled`` and ``sixframe_noncoding``;
* ``sixframe_noncoding`` outputs (unscored and scored_null) stay byte for byte as at main 98f89786.

Small and fast: ``novel_orf_power(30, "100k")`` (about 0.7 s a run), checked against the makedb
database GIPPA itself builds (``m7_round3_oracle.makedb_entries``). The pairing check reproduces
GIPPA's own picked class-FDR key (``gippa.novel.fdr``: contig, strand, stop codon) on the IDs the
search reports, which is all ``makedb_id_map`` builds entries from.
"""

from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from helpers import m7_round3_oracle as Q
from helpers import m7_sim_oracle as O

import gippasim
import gippasim.scenarios as S
from gippa.novel import fdr as F

SEED = 30
FILES = ("peptides.tsv", "truth/peptides.tsv", "truth/proteomics.json")
#: sha256 at main 98f89786 (before the fix): novel_orf_power(30, "100k"), seed 30.
SIXFRAME_SHA = {
    "peptides.tsv": "547680c738e0978829319627b8c1148e99acef6e7993b97e67b26bce12209e94",
    "truth/peptides.tsv": "3585a5247fd095a8588aa9807f198c5cc2b2f2676c15c3dd31e1eb5da6c34fdf",
    "truth/proteomics.json": "473dbecb56a4cd26f8772b256e7f1a673f87d1e270d98d0d838512107235e360",
}
#: the same with ``decoys = "scored_null"``.
SIXFRAME_SCORED_SHA = {
    "peptides.tsv": "fa7842e8436565f26a1a69652a71f8301e8e23ceb9e828d42685e7133210efa4",
    "truth/peptides.tsv": "75d534320e244ae6416760eb80922a52cf62491c674c0d76c360e7712fbc24be",
    "truth/proteomics.json": "3ffbe3d247a246d5c7fada56220b67308dffe16d66d31d923d6041ebee9c0428",
}
NULL_ORIGINS = {"undetected_product", "sixframe_noncoding"}


def config(wrong_source: str, decoys: str = "scored_null") -> Any:
    c = S.novel_orf_power(SEED, "100k")
    return replace(c, proteomics=replace(c.proteomics, decoys=decoys, wrong_source=wrong_source))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def outs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """scored_null runs with the k12_evidence default (``shuffled``) and with six-frame."""
    res = {}
    for ws in ("shuffled", "sixframe_noncoding"):
        out = tmp_path_factory.mktemp(f"scored_{ws}")
        gippasim.simulate(config(ws), SEED, out)
        res[ws] = out
    return res


@pytest.fixture(scope="module")
def entries(outs: dict[str, Path], tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """``gippa makedb --six-frame-orfs`` on the shuffled run (same genome and annotation)."""
    return Q.makedb_entries(outs["shuffled"], tmp_path_factory.mktemp("makedb"))


def _targets(out: Path) -> dict[str, set[str]]:
    """Target peptide -> every protein ID its rows report."""
    ids: dict[str, set[str]] = {}
    for r in O.rows(out / "peptides.tsv"):
        if r["is_decoy"] != "true":
            ids.setdefault(r["peptide"], set()).update(Q.ids(r["protein_ids"]))
    return ids


def _novel(out: Path) -> dict[str, dict[str, str]]:
    return {r["orf_id"]: r for r in O.rows(out / "truth" / "novel_orfs.tsv")}


def _novel_rows(out: Path) -> list[dict[str, str]]:
    """Truth rows of true peptides of planted novel ORFs (``locus_tag = NOVEL_…``)."""
    novel = _novel(out)
    return [t for t in O.rows(out / "truth" / "peptides.tsv")
            if t["is_wrong"] != "true" and t["locus_tag"] in novel]  # fmt: skip


# ---------------------------------------------------------------- novel ORF target IDs


@pytest.mark.parametrize("ws", ["shuffled", "sixframe_noncoding"])
def test_novel_orf_targets_name_no_truth_id(outs: dict[str, Path], ws: str) -> None:
    out = outs[ws]
    targets = _targets(out)
    bad = sorted({i for ids in targets.values() for i in ids if i.startswith("NOVEL_")})
    assert not bad, f"{ws}: {len(bad)} target protein IDs are gippasim truth IDs: {bad[:5]}"


@pytest.mark.parametrize("ws", ["shuffled", "sixframe_noncoding"])
def test_novel_orf_targets_name_the_makedb_entry_they_lie_in(
    outs: dict[str, Path], entries: dict[str, Any], ws: str
) -> None:
    out = outs[ws]
    targets = _targets(out)
    rows = [t for t in _novel_rows(out) if t["peptide"] in targets]
    assert len(rows) >= 50, f"{ws}: only {len(rows)} detected novel-ORF truth rows"
    idx = Q.EntryIndex(entries, Q.lengths(out))
    unknown = sorted({i for t in rows for i in targets[t["peptide"]]
                      if i.startswith("GIPPAORF_") and i not in entries})  # fmt: skip
    assert not unknown, f"{ws}: IDs that name no makedb entry: {unknown[:5]}"
    missing = []
    for t in rows:
        home = {e.db_id for e in idx.entry_of(t["contig"], int(t["start"]), int(t["end"]),
                                              t["strand"])}  # fmt: skip
        if not home & targets[t["peptide"]]:
            missing.append((t["peptide"], t["locus_tag"], sorted(targets[t["peptide"]])))
    assert not missing, (
        f"{ws}: {len(missing)} / {len(rows)} novel-ORF target peptides do not name the makedb "
        f"entry they lie in: {missing[:3]}"
    )


@pytest.mark.parametrize("ws", ["shuffled", "sixframe_noncoding"])
def test_true_novel_orfs_pair_with_a_reported_makedb_entry(outs: dict[str, Path], ws: str) -> None:
    """Integration-light: GIPPA's picked class FDR pairs an ORF with an entry by (contig, strand,
    stop) among the class-tagged IDs the search reported (``makedb_id_map``). Every planted
    non-``nterm_extension`` ORF with a detected peptide must find one (else ``no_makedb_pair``)."""
    out = outs[ws]
    targets = _targets(out)
    lengths = Q.lengths(out)
    reported = {F._strip(i, F.DEFAULT_DECOY_PREFIXES) for ids in targets.values() for i in ids}
    ents = F._entries(None, reported)
    keys = {(e.contig, e.strand, F._stop(e, lengths)) for e in ents.values()
            if e.cls != "nterm_extension"}  # fmt: skip
    found: dict[str, bool] = {}
    novel = _novel(out)
    for t in _novel_rows(out):
        o = novel[t["locus_tag"]]
        if o["class"] == "nterm_extension" or t["peptide"] not in targets:
            continue
        key = (o["contig"], o["strand"], (int(o["stop"]) - 1) % lengths[o["contig"]])
        found[o["orf_id"]] = key in keys
    assert len(found) >= 20, f"{ws}: only {len(found)} planted ORFs with a detected peptide"
    unpaired = sorted(k for k, v in found.items() if not v)
    assert not unpaired, (
        f"{ws}: {len(unpaired)} / {len(found)} would be no_makedb_pair: {unpaired[:5]}"
    )


# ---------------------------------------------------------------- N-terminal extensions


@pytest.mark.parametrize("ws", ["shuffled", "sixframe_noncoding"])
def test_upstream_nterm_peptide_names_its_nterm_extension_entry(
    outs: dict[str, Path], entries: dict[str, Any], ws: str
) -> None:
    out = outs[ws]
    targets = _targets(out)
    feats = {r["feature_id"]: r for r in O.rows(out / "truth" / "features.tsv")}
    nterm = {}  # gene locus tag -> its nterm_extension entry IDs
    for e in entries.values():
        if e.cls == "nterm_extension" and e.locus_tag:
            nterm.setdefault(e.locus_tag, set()).add(e.db_id)
    checked, bad = 0, []
    for t in O.rows(out / "truth" / "peptides.tsv"):
        f = feats.get(t["locus_tag"])
        if t["origin"] != "gene" or f is None or not f["true_start"] or t["peptide"] not in targets:
            continue
        up = (int(t["start"]) < int(f["start"]) if f["strand"] == "+"
              else int(t["end"]) > int(f["end"]))  # fmt: skip
        if not up or t["locus_tag"] not in nterm:
            continue
        checked += 1
        if not nterm[t["locus_tag"]] & targets[t["peptide"]]:
            bad.append((t["peptide"], t["locus_tag"], sorted(targets[t["peptide"]])))
    assert checked >= 5, f"{ws}: only {checked} detected upstream N-terminal peptides"
    assert not bad, (
        f"{ws}: {len(bad)} / {checked} upstream peptides do not name their gene's "
        f"nterm_extension entry: {bad[:3]}"
    )


# ---------------------------------------------------------------- the wrong-peptide draw


def test_scored_null_outputs_do_not_depend_on_wrong_source(outs: dict[str, Path]) -> None:
    """The scored_null nulls are drawn from both sources whatever ``wrong_source``; once the
    target IDs follow makedb in both, the two runs are the same file for file."""
    diff = [f for f in FILES
            if sha(outs["shuffled"] / f) != sha(outs["sixframe_noncoding"] / f)]  # fmt: skip
    assert not diff, f"scored_null outputs differ between shuffled and sixframe_noncoding: {diff}"


@pytest.mark.parametrize("ws", ["shuffled", "sixframe_noncoding"])
def test_scored_null_wrong_draw_is_protein_and_sixframe_nulls_only(
    outs: dict[str, Path], ws: str
) -> None:
    """No shuffled (or sixframe_noncoding-mode) wrong source in scored_null: every wrong truth row
    is a protein-entry null (locus tag ID) or a six-frame null (``GIPPAORF_`` ID)."""
    out = outs[ws]
    targets = _targets(out)
    wrong = [t for t in O.rows(out / "truth" / "peptides.tsv")
             if t["is_wrong"] == "true" and t["origin"] != "decoy"]  # fmt: skip
    assert wrong
    assert {t["origin"] for t in wrong} <= NULL_ORIGINS, {t["origin"] for t in wrong}
    for t in wrong:
        ids = targets.get(t["peptide"], set())
        if not ids:  # every PSM of this null went to its decoy (rule 6 coin)
            continue
        if t["origin"] == "sixframe_noncoding":
            assert all(i.startswith("GIPPAORF_") for i in ids), (t, ids)
        else:
            assert not any(i.startswith(("GIPPAORF_", "NOVEL_")) for i in ids), (t, ids)


# ---------------------------------------------------------------- sixframe_noncoding unchanged


@pytest.mark.parametrize(
    ("decoys", "expected"),
    [("unscored", SIXFRAME_SHA), ("scored_null", SIXFRAME_SCORED_SHA)],
)
def test_sixframe_noncoding_outputs_unchanged(
    tmp_path: Path, decoys: str, expected: dict[str, str]
) -> None:
    gippasim.simulate(config("sixframe_noncoding", decoys), SEED, tmp_path)
    assert {f: sha(tmp_path / f) for f in FILES} == expected
